/*
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @author: Simon Mueller, 2022
*/


#pragma once
#include "general.hpp"
#include <fstream>
#include <iomanip>

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
void Write1DArrayToFile(std::string path, const T* array_ptr, int n_rows, int n_cols,
	bool transpose = false, int print_to_row_n = -1, int print_to_col_n = -1,
	int n_repeat = 1) {

	// Adjust print_to_row_n and print_to_col_n if they are not provided
	if (print_to_row_n == -1) {
		print_to_row_n = n_rows;
	}
	if (print_to_col_n == -1) {
		print_to_col_n = n_cols;
	}

	// Open file using RAII (std::ofstream)
	std::ofstream file(path.data());
	if (!file) {
		throw std::runtime_error("Could not open file: " + std::string(path));
	}

	// Set formatting options: fixed point, 6 decimal precision
	file << std::fixed << std::setprecision(6);

	// Loop through the data, handling transposition if needed
	if (transpose) {
		for (int h = 0; h < n_repeat; ++h) {
			for (int i = 0; i < print_to_col_n; ++i) {
				for (int j = 0; j < print_to_row_n; ++j) {
					file << std::setw(14) << array_ptr[j * n_cols + i] << "  ";
				}
				file << "\n";
			}
		}
	}
	else {
		for (int h = 0; h < n_repeat; ++h) {
			for (int i = 0; i < print_to_row_n; ++i) {
				for (int j = 0; j < print_to_col_n; ++j) {
					file << std::setw(14) << array_ptr[i * n_cols + j] << "  ";
				}
				file << "\n";
			}
		}
	}

	// Check if any errors occurred during writing
	if (!file) {
		throw std::runtime_error("Error occurred while writing to file: " + std::string(path));
	}
}

void WriteEigenMatrixToFile(std::string path, const Eigen::MatrixXd& m) {
	std::ofstream file(path.data());
	if (!file) {
		throw std::runtime_error("Could not open file: " + std::string(path));
	}

	file.precision(6);  // Set precision to 6 digits
	file << std::scientific;  // Use scientific notation

	for (int i = 0; i < m.rows(); ++i) {
		for (int j = 0; j < m.cols(); ++j) {
			file << m(i, j) << " ";
		}
		file << "\n";
	}

	if (!file) {
		throw std::runtime_error("Error occurred while writing to file: " + std::string(path));
	}
}

void WriteExtendedSigmaProfileToFile(std::string path, segmentTypeCollection& segments) {

	std::ofstream file(path.data());
	if (!file) {
		throw std::runtime_error("Could not open file: " + std::string(path));
	}

	file << std::fixed << std::setprecision(6);

	// Write segment data to file
	for (int i = 0; i < segments.size(); ++i) {
		file << std::setw(4) << i << "  "
			<< std::setw(14) << segments.SegmentTypeSigma[i] << "  "
			<< std::setw(14) << segments.SegmentTypeSigmaCorr[i] << "  "
			<< std::setw(4) << segments.SegmentTypeAtomicNumber[i] << "  "
			<< std::setw(3) << segments.SegmentTypeHBtype[i] << "  "
			<< std::setw(3) << segments.SegmentTypeGroup[i] << "  "
			<< std::setw(14) << segments.SegmentTypeAreas[i][0] << "\n";
	}

	if (!file) {
		throw std::runtime_error("Error occurred while writing to file: " + std::string(path));
	}
}