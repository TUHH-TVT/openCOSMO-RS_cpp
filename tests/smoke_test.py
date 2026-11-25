"""Smoke test for built wheels: import the module and run the bundled example."""
import contextlib
import io
import os
import runpy

import openCOSMORS  # noqa: F401

example = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bindings", "run_example.py")
buffer = io.StringIO()
with contextlib.redirect_stdout(buffer):
    runpy.run_path(example, run_name="__main__")
output = buffer.getvalue()

assert "ERROR:" not in output, output
assert "ln_gamma" in output, output
assert "average_surface_energies" in output, output
print(output)
print("smoke test passed")
