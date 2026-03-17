# =============================================================================
# openCOSMO-RS_cpp Dockerfile
# Multi-stage build: Python module (.so) + CLI binary
# SIMD level auto-detected from host CPU at build time
# =============================================================================

# ---------------------------------------------------------------------------
# Stage 1: Builder
# ---------------------------------------------------------------------------
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        python3-dev \
        python3-pip \
        libgomp1 \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir "pybind11[global]" numpy

# ---------- detect best SIMD level the *build* host supports ---------------
# Precedence: FMA > AVX > SSE3 (fallback)
# Produces /tmp/simd_define  (e.g. -D__FMA__)
#      and /tmp/simd_march   (e.g. -mfma -mavx)
RUN SIMD="__SSE3__"; MARCH="-msse3"; \
    if grep -qw avx /proc/cpuinfo 2>/dev/null; then \
        SIMD="__AVX__"; MARCH="-mavx"; \
    fi; \
    if grep -qw fma /proc/cpuinfo 2>/dev/null; then \
        SIMD="__FMA__"; MARCH="-mfma -mavx"; \
    fi; \
    echo "$SIMD"  > /tmp/simd_define && \
    echo "$MARCH" > /tmp/simd_march  && \
    echo "Detected SIMD: $SIMD ($MARCH)"

WORKDIR /src
COPY . .

# -- populate empty eigen submodule (need ≥3.4 for Eigen::indexing::all) ------
RUN rm -rf eigen && \
    git clone --depth 1 \
        https://gitlab.com/libeigen/eigen.git eigen

# -- populate empty pybind11 submodule via shallow clone --------------------
RUN rm -rf pybind11 && \
    git clone --depth 1 https://github.com/pybind/pybind11.git pybind11

# -- patch CMakeLists.txt ---------------------------------------------------
# 1) Remove the hard-coded SIMD COMPILE_FLAGS line; we pass flags via CMAKE_CXX_FLAGS
# 2) Add project root to include dirs so #include "nlohmann/json.hpp" resolves
RUN sed -i '/COMPILE_FLAGS -D__AVX__/d' CMakeLists.txt && \
    sed -i 's|set(MY_INCLUDES|set(MY_INCLUDES\n    "${CMAKE_CURRENT_SOURCE_DIR}"|' CMakeLists.txt

# -- build Python module (.so) ---------------------------------------------
RUN SIMD_DEF=$(cat /tmp/simd_define) && SIMD_MARCH=$(cat /tmp/simd_march) && \
    echo "Building Python module with: -D${SIMD_DEF} ${SIMD_MARCH}" && \
    mkdir build_py && cd build_py && \
    cmake .. -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_CXX_FLAGS="-D${SIMD_DEF} ${SIMD_MARCH}" && \
    cmake --build . --config Release -j"$(nproc)" && \
    echo "Python module built:" && ls -lh openCOSMORS*.so

# -- build CLI binary -------------------------------------------------------
RUN SIMD_DEF=$(cat /tmp/simd_define) && SIMD_MARCH=$(cat /tmp/simd_march) && \
    echo "Building CLI binary with: -D${SIMD_DEF} ${SIMD_MARCH}" && \
    mkdir build_cli && cd build_cli && \
    cmake .. -DBINARY= -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_CXX_FLAGS="-D${SIMD_DEF} ${SIMD_MARCH}" && \
    cmake --build . --config Release -j"$(nproc)" && \
    echo "CLI binary built:" && ls -lh openCOSMORS

# ---------------------------------------------------------------------------
# Stage 2: Slim runtime image
# ---------------------------------------------------------------------------
FROM ubuntu:22.04 AS runtime

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        libgomp1 \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir numpy

WORKDIR /app

# -- copy CLI binary -------------------------------------------------------
COPY --from=builder /src/build_cli/openCOSMORS /usr/local/bin/openCOSMORS

# -- copy Python module -----------------------------------------------------
COPY --from=builder /src/build_py/openCOSMORS*.so /app/bindings/

# -- copy example / data files -----------------------------------------------
COPY --from=builder /src/bindings/*.orcacosmo  /app/bindings/
COPY --from=builder /src/bindings/run_example.py /app/bindings/

# -- record which SIMD was used ----------------------------------------------
COPY --from=builder /tmp/simd_define /app/.simd_level

ENV PYTHONPATH=/app/bindings

WORKDIR /app/bindings

CMD ["python3", "run_example.py"]
