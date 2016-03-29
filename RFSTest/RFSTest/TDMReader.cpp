#include "TDMReader.h"
#include <iostream>
#include <cstring>
#include <sstream>

using namespace std;

TDMReader::TDMReader(): headerRead(false), buffer(""), header("") {}

/**
 * Constructor of the TDMReader class instance.
 * Initializes the input.
 * @brief TDMReader::TDMReader
 * @param filename
 */
TDMReader::TDMReader(string filename): headerRead(false), buffer(""), header("")
{   
    inputFile.open(filename.c_str());

    if (!inputFile.is_open())
        cout << "Failed to open input file!" << endl;
}

/**
 * Destructor of the TDMReader class instance.
 * Closes the input file.
 * @brief TDMReader::~TDMReader
 */
TDMReader::~TDMReader()
{
    inputFile.close();
}

/**
 * Reads the header of the CAMRa TDM file. The header ends when the data starts.
 * @brief TDMReader::readHeader
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
 * Reads one data entry from the TDM file.
 * @brief TDMReader::readEntry
 * @return
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
            cout << "End of file reached." << endl;
            result.clear();
            return result;
        }
    }

    result.push_back( stod( elements.back()));

    for (int i = 0; i < 4; i++) {
        buffer.clear();
        getline(inputFile, buffer);
        split(buffer,' ', elements);
        result.push_back( stod( elements.back()));
    }

    date = elements.at(2);
    split (date,'T',elements);

    date = elements.at(0);
    time = elements.at(1);

    split (date,'-',elements);

    for (int i = 0; i< 3; i++)
        result.push_back(stod(elements.at(i)));

    split(time,':',elements);

    for (int i = 0; i < 3; i++)
        result.push_back(stod(elements.at(i)));

    return result;
}

Eigen::VectorXd TDMReader::readEntryEigen()
{
	vector<double> v = readEntry();
	Eigen::VectorXd e(v.size());
	for (size_t i = 0; i < v.size(); i++)
		e(i) = v[i];

	return e;
}

/**
 * Accessor of the header of the current TDM file.
 * Reads the header if it was not read previously.
 * @brief TDMReader::getHeader
 * @return
 */
string& TDMReader::getHeader() {

    if (!headerRead) {
        readHeader();
        headerRead = true;
    }

    return this->header;
}

/**
 * <summary> Splits a string according to the delimeters specified. </summary>
 * Source: http://stackoverflow.com/questions/236129/split-a-string-in-c
 * <return> A vector of strings containing the elements of the initial string </return>
 */
vector<string>& TDMReader::split(const string &s, char delim, vector<string> &elements) {\

    std::stringstream ss(s);
    std::string item;

    if (elements.size() != 0)
        elements.clear();

    while (std::getline(ss, item, delim)) {
        if (item.size() > 0)
            elements.push_back(item);
    }

    return elements;
}
