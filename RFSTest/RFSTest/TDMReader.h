#ifndef TDMREADER_H
#define TDMREADER_H

#include <fstream>
#include <vector>
#include <Eigen\Eigen>

using namespace std;

enum TDM_PARAMS {RANGE, AZI, ELE, SS, CPSS};

class TDMReader
{
public:
	TDMReader();
    TDMReader(string filename);
    ~TDMReader();

    vector<double> readEntry();
	Eigen::VectorXd readEntryEigen();

    string& getHeader();

private:
    ifstream inputFile;
    string buffer, header;
    bool headerRead;

    void readHeader();
    vector<string>& split(const string&, char, vector<string>&);

    static constexpr const char* DATA_DELIM_START   = "DATA_START";
    static constexpr const char* DATA_DELIM_END     = "DATA_STOP";
    static constexpr const char* META_DELIM_START   = "META_START";
    static constexpr const char* META_DELIM_END     = "META_STOP";

    const unsigned short DATA_DELIM_START_LEN = 10;
    const unsigned short DATA_DELIM_END_LEN = 9;
    const unsigned short META_DELIM_START_LEN = 10;
    const unsigned short META_DELIM_END_LEN = 9;
};

#endif // TDMREADER_H
