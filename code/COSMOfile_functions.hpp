/*
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @author: Simon Mueller, 2022
*/


#pragma once
#include <stdarg.h>
#include <string>
#include <fstream>
#include "helper_functions.hpp"

#define NUM_ARGS(...) std::tuple_size<decltype(std::make_tuple(__VA_ARGS__))>::value
#define parse_line(haystack, format, ...)  (_parse_line(NUM_ARGS(__VA_ARGS__), haystack, format,  __VA_ARGS__))

void _parse_line(int numberOfArgumentsToParse, const std::string& haystack, const char* format, ...) {
    const char* str = haystack.c_str();
    va_list args;

    va_start(args, format);
    int numberOfArgumentsFound = vsscanf(str, format, args);
    va_end(args);

    if (numberOfArgumentsFound != numberOfArgumentsToParse) {
        throw std::runtime_error("The following line could not be parsed: " + haystack);
    }
}

std::string periodicTableElements[118] = { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", \
                                            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", \
                                            "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", \
                                            "Cs", "Ba", \
                                            "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", \
                                            "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", \
                                            "Fr", "Ra", \
                                            "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", \
                                            "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og" };

int findAtomicNumberByName(std::string& name) {
    for (int i = 0; i < 118; i++) {
        if (periodicTableElements[i] == name) {
            return i + 1;
        }
    }
    throw std::runtime_error("Could not find the atomic number for followng atom: " + name);
}

std::string scan_for(std::ifstream& haystack_file, std::string needle, std::string mode = "start", int numberOfLinesToSkip = 0, bool throwErrorIfNotFound = true) {
    
    std::string currentLine;
    std::string matchingLine = "";

    while (std::getline(haystack_file, currentLine))
    {
        currentLine = trim(currentLine);

        if (mode == "start")
            if (startsWith(currentLine, needle)) {
                matchingLine = currentLine;
                break;
            }

        if (mode == "end")
            if (endsWith(currentLine, needle)) {
                matchingLine = currentLine;
                break;
            }
    }

    for (int i = 0; i < numberOfLinesToSkip; i++)
        std::getline(haystack_file, currentLine);

    if (matchingLine != "")
        return matchingLine;

    if (throwErrorIfNotFound) {
        throw std::runtime_error("the following string was not found while reading haystack_file: " + needle);
    }
    else {
        return "";
    }
}

molecule getMoleculeFromTurbomoleCOSMOfile(std::string& path) {

    std::ifstream cosmoFile(path);

    if (cosmoFile.fail()) {
        throw std::runtime_error("The following COSMOfile could not be opened: " + path);
    }

    std::string currentLine;
    std::string lastLine;

    int section = 1;

    std::vector<double> atomPositions_X;
    std::vector<double> atomPositions_Y;
    std::vector<double> atomPositions_Z;
    std::vector<double> atomRadii;
    std::vector<int> atomAtomicNumbers;

    std::vector<double> segmentPositions_X;
    std::vector<double> segmentPositions_Y;
    std::vector<double> segmentPositions_Z;
    std::vector<int> segmentAtomIndices;
    std::vector<double> segmentAreas;
    std::vector<double> segmentSigmas;

    molecule newMolecule;

    std::getline(cosmoFile, currentLine);
    std::getline(cosmoFile, currentLine);
    newMolecule.qmMethod = replace(trim(currentLine), "prog.: ", "");

    while (std::getline(cosmoFile, currentLine))
    {
        currentLine = trim(currentLine);

        if (startsWith(lastLine, "#atom")) {
            section = 2;
        }

        if (startsWith(currentLine, "$coord_car")) {
            section = 3;
        }

        if (startsWith(lastLine, "#  n   atom ")) {
            section = 4;
            std::getline(cosmoFile, currentLine);
            std::getline(cosmoFile, currentLine);
        }

        if (section == 1) {

            if (startsWith(currentLine, "area")) {
                currentLine = replace(trim(currentLine), " ", "");
                parse_line(currentLine, " %*[^=]=%lf", &(newMolecule.Area));
                newMolecule.Area = pow(0.529177249, 2) * newMolecule.Area;
            }

            if (startsWith(currentLine, "volume")) {
                currentLine = replace(trim(currentLine), " ", "");
                parse_line(currentLine, " %*[^=]=%lf", &(newMolecule.Volume));
                newMolecule.Volume = pow(0.529177249, 3) * newMolecule.Volume;
            }

            if (startsWith(currentLine, "Total energy + OC corr.")) {
                currentLine = replace(trim(currentLine), " ", "");
                parse_line(currentLine, " %*[^=]=%lf", &(newMolecule.epsilonInfinityTotalEnergy));
            }
        }

        if (section == 2) {
            double atomPosition_X, atomPosition_Y, atomPosition_Z, atomRadius;
            char temp_atomName[3];
            parse_line(currentLine, "%*i %lf %lf %lf %3s %lf", &atomPosition_X, &atomPosition_Y, &atomPosition_Z, &temp_atomName, &atomRadius);
            atomPositions_X.push_back(0.529177249 * atomPosition_X);
            atomPositions_Y.push_back(0.529177249 * atomPosition_Y);
            atomPositions_Z.push_back(0.529177249 * atomPosition_Z);

            temp_atomName[0] = toupper(temp_atomName[0]);
            std::string atomName(temp_atomName);
            atomName = trim(atomName);

            atomAtomicNumbers.push_back(findAtomicNumberByName(atomName));

            atomRadii.push_back(atomRadius);
        }

        if (section == 4) {
            int segmentAtomIndex;
            double segmentPosition_X, segmentPosition_Y, segmentPosition_Z, segmentArea, segmentSigma;
            parse_line(currentLine, "%*i %i %lf %lf %lf %*lf %lf %lf", &segmentAtomIndex, &segmentPosition_X, &segmentPosition_Y, &segmentPosition_Z, &segmentArea, &segmentSigma);
            segmentAtomIndices.push_back(segmentAtomIndex - 1);

            segmentPositions_X.push_back(0.529177249 * segmentPosition_X);
            segmentPositions_Y.push_back(0.529177249 * segmentPosition_Y);
            segmentPositions_Z.push_back(0.529177249 * segmentPosition_Z);

            segmentAreas.push_back(segmentArea);

            segmentSigmas.push_back(segmentSigma);
        }

        lastLine = currentLine;
    }
    cosmoFile.close();

    newMolecule.atomPositions = Eigen::MatrixXd::Zero(atomPositions_X.size(), 3);

    for (int i = 0; i < atomPositions_X.size(); i++) {
        newMolecule.atomPositions(i, 0) = atomPositions_X[i];
        newMolecule.atomPositions(i, 1) = atomPositions_Y[i];
        newMolecule.atomPositions(i, 2) = atomPositions_Z[i];
    }
    newMolecule.atomAtomicNumbers = Eigen::Map<Eigen::VectorXi>(atomAtomicNumbers.data(), atomAtomicNumbers.size());
    newMolecule.atomRadii = Eigen::Map<Eigen::VectorXd>(atomRadii.data(), int(atomRadii.size()));

    newMolecule.segmentPositions = Eigen::MatrixXd::Zero(segmentPositions_X.size(), 3);

    for (int i = 0; i < segmentPositions_X.size(); i++) {
        newMolecule.segmentPositions(i, 0) = segmentPositions_X[i];
        newMolecule.segmentPositions(i, 1) = segmentPositions_Y[i];
        newMolecule.segmentPositions(i, 2) = segmentPositions_Z[i];
    }

    newMolecule.segmentAtomIndices = Eigen::Map<Eigen::VectorXi>(segmentAtomIndices.data(), segmentAtomIndices.size());
    newMolecule.segmentAreas = Eigen::Map<Eigen::VectorXd>(segmentAreas.data(), segmentAreas.size());
    newMolecule.segmentSigmas = Eigen::Map<Eigen::VectorXd>(segmentSigmas.data(), segmentSigmas.size());

    return newMolecule;
}

molecule getMoleculeFromORCACOSMOfile(std::string& path) {

    std::ifstream cosmoFile(path);

    if (cosmoFile.fail()) {
        throw std::runtime_error("The following COSMOfile could not be opened: " + path);
    }

    std::string currentLine;
    std::string lastLine;

    std::vector<double> atomPositions_X;
    std::vector<double> atomPositions_Y;
    std::vector<double> atomPositions_Z;
    std::vector<double> atomRadii;
    std::vector<int> atomAtomicNumbers;

    std::vector<double> segmentPositions_X;
    std::vector<double> segmentPositions_Y;
    std::vector<double> segmentPositions_Z;
    std::vector<int> segmentAtomIndices;
    std::vector<double> segmentAreas;
    std::vector<double> segmentSigmas;

    molecule newMolecule;

    std::getline(cosmoFile, currentLine);

    std::vector<std::string> parts = split(currentLine, ':');
    newMolecule.qmMethod = replace(trim(parts[1]), "COSMO", "CPCM");
    newMolecule.name = trim(parts[0]);

    scan_for(cosmoFile, "#ENERGY");
    std::getline(cosmoFile, currentLine);
    currentLine = trim(replace(currentLine, "FINAL SINGLE POINT ENERGY", ""));
    newMolecule.epsilonInfinityTotalEnergy = std::stod(currentLine);

    scan_for(cosmoFile, "#XYZ_FILE");

    int numberOfAtoms;
    std::getline(cosmoFile, currentLine);
    parse_line(currentLine, "%i", &numberOfAtoms);
    std::getline(cosmoFile, currentLine);

    for (int i = 0; i < numberOfAtoms; i++) {

        std::getline(cosmoFile, currentLine);

        double atomPosition_X, atomPosition_Y, atomPosition_Z;
        char temp_atomName[3];
        parse_line(currentLine, "%3s %lf %lf %lf", &temp_atomName , &atomPosition_X, &atomPosition_Y, &atomPosition_Z);
        atomPositions_X.push_back(atomPosition_X);
        atomPositions_Y.push_back(atomPosition_Y);
        atomPositions_Z.push_back(atomPosition_Z);

        temp_atomName[0] = toupper(temp_atomName[0]);
        std::string atomName(temp_atomName);
        atomName = trim(atomName);

        atomAtomicNumbers.push_back(findAtomicNumberByName(atomName));
    }
    newMolecule.atomAtomicNumbers = Eigen::Map<Eigen::VectorXi>(atomAtomicNumbers.data(), atomAtomicNumbers.size());

    currentLine = scan_for(cosmoFile, "# Volume", "end");
    parse_line(currentLine, " %lf", &(newMolecule.Volume));
    newMolecule.Volume = pow(0.529177249, 3) * newMolecule.Volume;

    currentLine = scan_for(cosmoFile, "# Area", "end");
    parse_line(currentLine, " %lf", &(newMolecule.Area));
    newMolecule.Area = pow(0.529177249, 2) * newMolecule.Area;

    currentLine = scan_for(cosmoFile, "# CPCM dielectric energy", "end");
    float uncorrectedDialectricEnergy = 0.0;
    parse_line(currentLine, "%f", &uncorrectedDialectricEnergy);

    scan_for(cosmoFile, "# CARTESIAN COORDINATES (A.U.) + RADII (A.U.)", "start", 1);
    for (int i = 0; i < numberOfAtoms; i++) {
        std::getline(cosmoFile, currentLine);
        double atomRadius;
        parse_line(currentLine, "%*lf %*lf %*lf %lf", &atomRadius);
        atomRadii.push_back(atomRadius * 0.529177249);
    }
    newMolecule.atomRadii = Eigen::Map<Eigen::VectorXd>(atomRadii.data(), int(atomRadii.size()));

    scan_for(cosmoFile, "# SURFACE POINTS (A.U.)", "start", 2);
    while (std::getline(cosmoFile, currentLine)) {

        currentLine = trim(currentLine);

        if (currentLine != "") {

            int segmentAtomIndex;
            double segmentPosition_X, segmentPosition_Y, segmentPosition_Z, segmentArea, segmentCharge;
            parse_line(currentLine, "%lf %lf %lf %lf %*lf %lf %*lf %*lf %*lf %i", &segmentPosition_X, &segmentPosition_Y, &segmentPosition_Z, &segmentArea, &segmentCharge, &segmentAtomIndex);
            segmentAtomIndices.push_back(segmentAtomIndex);

            segmentPositions_X.push_back(0.529177249 * segmentPosition_X);
            segmentPositions_Y.push_back(0.529177249 * segmentPosition_Y);
            segmentPositions_Z.push_back(0.529177249 * segmentPosition_Z);

            segmentArea = 0.529177249 * 0.529177249 * segmentArea;
            segmentAreas.push_back(segmentArea);

            segmentSigmas.push_back(segmentCharge / segmentArea);
        }
        else {
            break;
        }
    }

    std::string matchingLine = scan_for(cosmoFile, "#COSMO_corrected", "start", 0, false);
    if (matchingLine != "") {
        std::getline(cosmoFile, currentLine);
        if (!startsWith(currentLine, "Corrected dielectric energy   =")) {
            throw std::runtime_error("Could not find the corrected dielectric energy entry in the orcacosmo file.");
        }
        currentLine = trim(replace(currentLine, "Corrected dielectric energy   =", ""));
        float correctedDielectricEnergy = 0.0;
        parse_line(currentLine, "%f", &correctedDielectricEnergy);
        newMolecule.epsilonInfinityTotalEnergy = newMolecule.epsilonInfinityTotalEnergy - uncorrectedDialectricEnergy + correctedDielectricEnergy;

        std::getline(cosmoFile, currentLine);
        std::getline(cosmoFile, currentLine);

        int segmentIndex = 0;
        while (std::getline(cosmoFile, currentLine)) {

            currentLine = trim(currentLine);

            if (currentLine == "##################################################")
                break;

            if (currentLine != "") {

                double correctedSegmentCharge;
                parse_line(currentLine, "%lf", &correctedSegmentCharge);

                segmentSigmas[segmentIndex] = correctedSegmentCharge / segmentAreas[segmentIndex];
                segmentIndex += 1;
            }
        }

        if (segmentIndex != segmentSigmas.size()) {
            throw std::runtime_error("Not enough corrected charges where found parsing the following file: " + path);
        }
    }
    cosmoFile.close();

    newMolecule.atomPositions = Eigen::MatrixXd::Zero(atomPositions_X.size(), 3);

    for (int i = 0; i < atomPositions_X.size(); i++) {
        newMolecule.atomPositions(i, 0) = atomPositions_X[i];
        newMolecule.atomPositions(i, 1) = atomPositions_Y[i];
        newMolecule.atomPositions(i, 2) = atomPositions_Z[i];
    }
    newMolecule.atomAtomicNumbers = Eigen::Map<Eigen::VectorXi>(atomAtomicNumbers.data(), atomAtomicNumbers.size());
    newMolecule.atomRadii = Eigen::Map<Eigen::VectorXd>(atomRadii.data(), int(atomRadii.size()));

    newMolecule.segmentPositions = Eigen::MatrixXd::Zero(segmentPositions_X.size(), 3);
    for (int i = 0; i < segmentPositions_X.size(); i++) {
        newMolecule.segmentPositions(i, 0) = segmentPositions_X[i];
        newMolecule.segmentPositions(i, 1) = segmentPositions_Y[i];
        newMolecule.segmentPositions(i, 2) = segmentPositions_Z[i];
    }

    newMolecule.segmentAtomIndices = Eigen::Map<Eigen::VectorXi>(segmentAtomIndices.data(), segmentAtomIndices.size());
    newMolecule.segmentAreas = Eigen::Map<Eigen::VectorXd>(segmentAreas.data(), segmentAreas.size());
    newMolecule.segmentSigmas = Eigen::Map<Eigen::VectorXd>(segmentSigmas.data(), segmentSigmas.size());

    return newMolecule;
}
