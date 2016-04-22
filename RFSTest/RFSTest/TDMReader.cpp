#include "TDMReader.h"
#include <iostream>
#include <cstring>
#include <sstream>

using namespace std;

TDMReader::TDMReader(): headerRead(false), buffer(""), header("") {}

/**
 * <summary> Constructor of the TDMReader class. </summary>
 * <par> Initializes the input stream. </par>
 * <param name = "_filename"> Input filename. </param>
 */
TDMReader::TDMReader(string _filename): headerRead(false), buffer(""), header("")
{   
    inputFile.open(_filename.c_str());

    if (!inputFile.is_open())
        cout << "TDMReader: Failed to open input file!" << endl;
}

/**
 * <summary> Destructor of the TDMReader class instance. </summary>
 * <par> Closes the input stream. </par>
 */
TDMReader::~TDMReader()
{
    inputFile.close();
}

/**
 * <summary> Reads the header of the CAMRa TDM file. </summary>
 * <par> The header ends when the DATA_START delimeter is encountered. </par>
 */
void TDMReader::readHeader() {

    // TODO: add eof() check (?)
    if (inputFile.is_open()) {

        // Look for the data delimeters
        while (!headerRead && !inputFile.eof()) {

            getline(inputFile, buffer);

            if (strcmp(buffer.substr(0, DATA_DELIM_START_LEN).c_str(), DATA_DELIM_START) == 0) {
                headerRead = true;
                buffer.clear();
                break;
            }

            header += (buffer + "\n");
            buffer.clear();
        }
    }
}

/**
 * <summary> Reads a data entry from the TDM file. </summary>
 * <returns> A vector<double> containing range, azimuth, elevation,
 * signal strength, cross-polar signal strength and date. </returns>
 */
vector<double> TDMReader::readEntry() {

    if (!headerRead) {
        readHeader();
        headerRead = true;
    }

    vector<double> result;
    vector<string> elements;
    string date, time;

    // Read the first line
    buffer.clear();
    getline(inputFile, buffer);
    split(buffer,' ', elements);

    if (elements.size() == 1) {
        if (strcmp(elements.at(0).c_str(), DATA_DELIM_END) == 0) {
            cout << "TDMReader: End of file reached." << endl;
            result.clear();
            return result;
        }
    }

    result.push_back( stod( elements.back()));

	// Range, azimuth, elevation, ss, cpss
    for (int i = 0; i < 4; i++) {
        buffer.clear();
        getline(inputFile, buffer);
        split(buffer,' ', elements);
        result.push_back( stod( elements.back()));
    }

	// date
    date = elements.at(2);
    split (date,'T',elements);

    date = elements.at(0);		// Date YY-MM-DD
    time = elements.at(1);		// Time HH:MM:SS

    split (date,'-',elements);

    for (int i = 0; i< 3; i++)
        result.push_back(stod(elements.at(i)));

    split(time,':',elements);

    for (int i = 0; i < 3; i++)
        result.push_back(stod(elements.at(i)));

    return result;
}

/**
 * <summary> Reads a data entry to an Eigen::VectorXd. </summary>
 * <returns> A VectorXd containing range, azimuth, elevation,
 * signal strength, cross-polar signal strength and date. </returns>
 */
Eigen::VectorXd TDMReader::readEntryEigen()
{
	vector<double> v = readEntry();
	Eigen::VectorXd e(v.size());
	for (size_t i = 0; i < v.size(); i++)
		e(i) = v[i];

	return e;
}

/**
 * <summary> Accessor of the header of the current TDM file. </summary>
 * <returns> A reference to the string containing the header. </returns>
 */
string TDMReader::getHeader() {

    if (!headerRead) {
        readHeader();
        headerRead = true;
    }

    return this->header;
}

/**
 * <summary> Splits a string according to the delimeters specified. </summary>
 * Source: http://stackoverflow.com/questions/236129/split-a-string-in-c
 * <returns> A vector of strings containing the elements of the initial string </returns>
 */
vector<string>& TDMReader::split(const string &s, char delim, vector<string> &elements) {

    std::stringstream ss(s);
    std::string item;

    if (elements.size() != 0)
        elements.clear();

    while (std::getline(ss, item, delim)) 
        if (item.size() > 0)
            elements.push_back(item);

    return elements;
}
