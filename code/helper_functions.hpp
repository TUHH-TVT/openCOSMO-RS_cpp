/*
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @author: Simon Mueller, 2022
*/


#pragma once
#include "types.hpp"
#include <fstream>

std::function<void(std::string)> display;
std::function<void(std::string, unsigned long)> displayTime;

static inline std::string leftTrim(std::string val) {
	val.erase(val.begin(), std::find_if(val.begin(), val.end(), [](unsigned char c) {
		return !std::isspace(c);
		}));

	return val;
}

static inline std::string rightTrim(std::string val) {
	val.erase(std::find_if(val.rbegin(), val.rend(), [](unsigned char c) {
		return !std::isspace(c);
		}).base(), val.end());
	return val;
}

static inline std::string trim(std::string val) {
	val = rightTrim(val);
	val = leftTrim(val);
	return val;
}

static inline std::string replace(std::string haystack, const std::string& from, const std::string& to) {
	size_t start_pos = haystack.find(from);
	if (start_pos == std::string::npos)
		return haystack;
	haystack.replace(start_pos, from.length(), to);
	return haystack;
}

static inline bool startsWith(std::string& haystack, const std::string& needle) {
	return haystack.rfind(needle, 0) == 0;
}

static inline bool endsWith(std::string& haystack, const std::string& needle) {

	if (needle.size() > haystack.size()) return false;

	return std::equal(needle.rbegin(), needle.rend(), haystack.rbegin());
}

static inline std::vector<std::string> split(const std::string& haystack, char delimiter)
{
	std::vector<std::string> retVal;
	std::stringstream _stringstream(haystack);
	std::string token;
	while (std::getline(_stringstream, token, delimiter)) {
		retVal.push_back(token);
	}
	return retVal;
}
static inline double round_and_truncate(double number_val, int n) {
	number_val = round(number_val / pow(10, n)) * pow(10, n);
	return number_val;
}

// Returns next greater multiple of 8
// if x is divisible by 8 it returns x
static inline int RoundUpToNextMultipleOfEight(int x) {
	return ((x + 7) & (-8));
}

// Returns previous multiple of 8
// if x is divisible by 8 it returns x
static inline int RoundDownToNextMultipleOfEight(int x) {
	if (x % 8 == 0) {
		return x;
	}

	return RoundUpToNextMultipleOfEight(x) - 8;
}

template<typename T>
std::string convertPointerToString(T* pointer) {
	std::ostringstream oss;
	oss << pointer;
	return oss.str();
}

template <typename T>
void Write1DArraytoFile(std::string Path, T* Array_ptr, int N_rows, int N_cols, bool transpose = false, int print_to_row_N = -1, int print_to_col_N = -1, int n_repeat = 1) {

	if (print_to_row_N == -1) {
		print_to_row_N = N_rows;
	}
	if (print_to_col_N == -1) {
		print_to_col_N = N_cols;
	}

	FILE* fp = NULL;
	fp = fopen(Path.c_str(), "w");
	if (fp == NULL) {
		throw std::runtime_error("Could not open file: " + Path);
	}

	if (transpose) {
		for (int h = 0; h < n_repeat; h++) {
			for (int i = 0; i < print_to_col_N; i++) {
				for (int j = 0; j < print_to_row_N; j++) {
					fprintf(fp, "%14.6e  ", Array_ptr[j * N_cols + i]);
				}
				fprintf(fp, "\n");
			}
		}
	}
	else {
		for (int h = 0; h < n_repeat; h++) {
			for (int i = 0; i < print_to_row_N; i++) {
				for (int j = 0; j < print_to_col_N; j++) {
					fprintf(fp, "%14.6e  ", Array_ptr[i * N_cols + j]);
				}
				fprintf(fp, "\n");
			}
		}
	}

	fclose(fp);
}

void WriteEigenMatrixtoFile(std::string Path, Eigen::MatrixXd m) {
	FILE* fp = NULL;
	fp = fopen(Path.c_str(), "w");
	if (fp == NULL) {
		throw std::runtime_error("Could not open file: " + Path);
	}

	for (int i = 0; i < m.rows(); i++) {
		for (int j = 0; j < m.cols(); j++) {
			fprintf(fp, "%.6e ", m(i, j));
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

void WriteExtendedSigmaProfiletoFile(std::string Path, segmentTypeCollection& segments) {

	FILE* fp = NULL;
	fp = fopen(Path.c_str(), "w");
	if (fp == NULL) {
		throw std::runtime_error("Could not open file: " + Path);
	}

	for (int i = 0; i < segments.size(); i++) {
		fprintf(fp, "%4d  %14.6e  %14.6e  %4d  %3d  %3d  %14.6e\n", i, segments.SegmentTypeSigma[i], segments.SegmentTypeSigmaCorr[i], segments.SegmentTypeAtomicNumber[i], \
			segments.SegmentTypeHBtype[i], segments.SegmentTypeGroup[i], segments.SegmentTypeAreas[i][0]);
	}

	fclose(fp);
}