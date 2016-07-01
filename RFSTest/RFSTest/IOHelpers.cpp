#include "IOHelpers.h"
#include "KalmanFilter.h"

#include <ctime>
#include <iostream>
#include <sstream>
#include <dirent.h>
#include <yaml-cpp/yaml.h>

using namespace std;
using namespace Eigen;

/**
 * <summary> A namespace for input / output helpers. </summary>
 */
namespace IOHelpers 
{
	/*
	* <summary> Empty constructor of the TDMReader class. </summary>
	*/
	TDMReader::TDMReader() : headerRead(false), buffer(""), header("") {}

	/*
	* <summary> Default constructor of the TDMReader class. </summary>
	* <param name = "_filename"> String containing the filename of the .tdm source file. </param>
	*/
	TDMReader::TDMReader(string _filename) : headerRead(false), buffer(""), header("")
	{
		inputFile.open(_filename.c_str());

		if (!inputFile.is_open())
			std::cout << "TDMReader: Failed to open input file!" << endl;
	}

	/**
	 * <summary> Destructor of the TDMReader class. </summary>
	 */
	TDMReader::~TDMReader()
	{
		if (inputFile.is_open())
			inputFile.close();
	}

	/**
	 * <summary> Reads the header and stores it in the header variable. </summary>
	 */
	void TDMReader::readHeader()
	{
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
	 * <summary> Header string accessor. </summary>
	 * <returns> String containing the header of the current .tdm file. </returns>
	 */
	string TDMReader::getHeader()
	{
		if (inputFile.is_open())
		{
			if (!headerRead) {
				readHeader();
				headerRead = true;
			}
			return header;
		}

		return "";
	}

	/**
	 * <summary> Reads an entry from the CAMRa-originated .tdm file. </summary>
	 * <returns> A vector of values in the entry: range, azimuth, elevation,
	 * signal strength, cross-polar signal strength </returns>
	 */
	vector<double> TDMReader::readEntry() 
	{
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
		split(buffer, ' ', elements);

		if (elements.size() == 1) {
			if (strcmp(elements.at(0).c_str(), DATA_DELIM_END) == 0) {
				std::cout << "TDMReader: End of file reached." << endl;
				result.clear();
				return result;
			}
		}

		result.push_back(stod(elements.back()));

		// Range, azimuth, elevation, ss, cpss
		for (int i = 0; i < 4; i++) {
			buffer.clear();
			getline(inputFile, buffer);
			split(buffer, ' ', elements);
			result.push_back(stod(elements.back()));
		}

		// Date
		date = elements.at(2);
		split(date, 'T', elements);

		date = elements.at(0);		// Date YY-MM-DD
		time = elements.at(1);		// Time HH:MM:SS

		split(date, '-', elements);

		for (int i = 0; i< 3; i++)
			result.push_back(stod(elements.at(i)));

		split(time, ':', elements);

		for (int i = 0; i < 3; i++)
			result.push_back(stod(elements.at(i)));

		return result;
	}

	/*
	 * <summary> Returns the VectorXd instead of vector<double> after the call of readEntry. </summary>
	 * <returns> A VectorXd containing the entry values: range, azimith, elevation, 
	 * signal strength, cross-polar signal strength. </returns>
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
	 * <summary> Opens the input file stream of the .file if it was not invoked in the constructor. </summary>
	 * <param name = "_filename"> A string containing the filename of the specfied file. </param>
	 */
	void TDMReader::open(const string & _filename)
	{
		inputFile.open(_filename.c_str());

		if (!inputFile.is_open())
			std::cout << "TDMReader: Failed to open input file!" << endl;
	}

	/**
	 * <summary> Transforms all .tdm files stored in the specified directory in the .csv format. </summary>
	 * <param name = "_path"> Directory path. </param>
	 */
	void TDMReader::TDMDirToCSV(const string & _path)
	{
		DIR *dir = opendir(_path.c_str());
		if (!dir)
			cout << "Directory not found" << endl;
		else
		{
			dirent *entry;
			ofstream ofscsv;
			string filePath;
			vector<double> dataEntry;

			// Read all files
			while (entry = readdir(dir))
			{
				// Check for .tdm extension
				if (has_suffix(entry->d_name, ".tdm"))
				{
					cout << entry->d_name << endl;
					filePath = _path + "\\" + entry->d_name;
					TDMReader tdmr(filePath);

					filePath.replace(filePath.end() - 3, filePath.end(), "csv");
					ofscsv.open(filePath);

					tdmr.getHeader();

					while (true) {

						dataEntry = tdmr.readEntry();

						if (!dataEntry.size())
							break;

						for (size_t i = 0; i < dataEntry.size() - 1; i++)
							ofscsv << dataEntry[i] << ",";

						ofscsv << dataEntry[dataEntry.size() - 1] << endl;
					}

					ofscsv.close();
				}
			}
		}

		closedir(dir);
	}

	/**
	 * <summary> Converts single .tdm file to YAML format. </summary> 
	 */
	void TDMReader::TDMToYAML(const string & _filename)
	{

	}

	/**
	 * <summary> Checks if the specified string contains the selected suffix. </summary>
	 * <param name = "_s"> Source string. </param>
	 * <param name = "_suffix"> Suffix string to be found. </param>
	 * <returns> A boolean indicating if the source string contains the suffix substring. </returns>
	 */
	bool has_suffix(const string& _s, const string& _suffix)
	{
		return (_s.size() >= _suffix.size()) && equal(_suffix.rbegin(), _suffix.rend(), _s.rbegin());
	}

	/**
	 * <summary> Reads the yaml parameters for the RFS SSA application. </summary>
	 * <param name = "_filename"> A .yaml parameter file. </param>
	 * <returns> Parameters structure. </returns>
	 */
	parameters readParametersYAML(const string & _filename)
	{
		YAML::Node config = YAML::LoadFile(_filename);

		parameters params;

		// Filename
		params.filename = config["filename"].as<std::string>();
		params.inputType = config["input_type"].as<std::string>();

		// Sensor-related
		params.sensorType = config["sensor_configuration"]["type"].as<std::string>();
		params.observationType = config["sensor_configuration"]["observation_type"].as<std::string>();
		params.observationDim = config["sensor_configuration"]["observation_size"].as<size_t>();

		params.sensorPosition = VectorXd(3);
		params.sensorPosition << config["sensor_configuration"]["latitude"].as<double>(),
			config["sensor_configuration"]["longitude"].as<double>(),
			config["sensor_configuration"]["altitude"].as<double>();

		params.pD = config["sensor_configuration"]["pd"].as<double>();
		params.kappa = config["sensor_configuration"]["kappa"].as<double>();
		params.lambda = config["sensor_configuration"]["lambda"].as<double>();
		params.clutter = config["sensor_configuration"]["clutter"].as<double>();

		params.R = MatrixXd::Zero(params.observationDim, params.observationDim);
		for (size_t i = 0; i < params.observationDim; i++)
			params.R(i, i) = config["sensor_configuration"]["R"][i].as<double>();

		// Single-target filter related
		params.singleTargetFilterType = config["filter_parameters_single_target"]["type"].as<std::string>();
		params.stateDim = config["filter_parameters_single_target"]["state_dim"].as<size_t>();
		params.motionModel = config["filter_parameters_single_target"]["motion_model"].as<std::string>();
		params.dt = config["filter_parameters_single_target"]["dt"].as<double>();

		if (!strcmp(params.motionModel.c_str(), "cv")) {
			params.Q = KalmanFilter::getCVQ(params.stateDim, params.dt);
			params.F = KalmanFilter::getCVF(params.stateDim, params.dt);
		}

		// Multiple-target filter related
		params.multipleTargetFilterType = config["filter_parameters_multiple_target"]["type"].as<std::string>();
		params.gaussianMixtureMaxSize = config["filter_parameters_multiple_target"]["gm_max_size"].as<size_t>();
		params.birthType = config["filter_parameters_multiple_target"]["birth_type"].as<std::string>();
		params.birthSize = config["filter_parameters_multiple_target"]["birth_num"].as<size_t>();

		params.birthCovariance = MatrixXd::Zero(params.stateDim, params.stateDim);
		for (size_t i = 0; i < params.stateDim; i++)
			params.birthCovariance(i, i) = config["filter_parameters_multiple_target"]["birth_cov"][i].as<double>();

		params.mergeThreshold = config["filter_parameters_multiple_target"]["th_merge"].as<double>();
		params.pruneThreshold = config["filter_parameters_multiple_target"]["th_prune"].as<double>();
		params.estimateThreshold = config["filter_parameters_multiple_target"]["th_estimate"].as<double>();

		params.birthIntensity = config["filter_parameters_multiple_target"]["birth_intensity"].as<double>();
		params.pB = config["filter_parameters_multiple_target"]["pb"].as<double>();
		params.pS = config["filter_parameters_multiple_target"]["ps"].as<double>();
		params.qInit = config["filter_parameters_multiple_target"]["q_initial"].as<double>();

		params.epsilon = config["filter_parameters_multiple_target"]["beta_epsilon"].as<double>();

		params.H = MatrixXd::Zero(params.observationDim, params.stateDim);
		for (size_t i = 0; i < params.observationDim; i++)
			params.H(i, i) = 1;

		return params;
	}

	/**
	 * <summary> Output stream operator for parameters structure. </summary>
	 * <params = "_os"> A reference to the output stream. </params>
	 * <params = "_params"> A constant reference to the parameters structure. </params>
	 */
	ostream & operator<<(std::ostream & _os, const parameters & _params)
	{
		_os << "Filename: " << _params.filename << endl;
		_os << endl;

		_os << "Sensor type: " << _params.sensorType << endl;
		_os << "\tObservation type: " << _params.observationType << endl;
		_os << "\tObservation size: " << _params.observationDim << endl;
		_os << "\tSensor position: [" << _params.sensorPosition.transpose() << "]" << endl;
		_os << "\tpD: " << _params.pD << endl;
		_os << "\tlambda: " << _params.lambda << endl;
		_os << "\tkappa: " << _params.kappa << endl;
		_os << "\tclutter: " << _params.clutter << endl;
		_os << "\tR:" << endl << _params.R << endl;
		_os << endl;

		_os << "Single target filter: " << endl;
		_os << "\tFilter type: " << _params.singleTargetFilterType << endl;
		_os << "\tState dimension: " << _params.stateDim << endl;
		_os << "\tMotion model: " << _params.motionModel << endl;
		_os << "\tdt: " << _params.dt << endl;
		_os << "\tF: " << endl << _params.F << endl;
		_os << "\tQ: " << endl << _params.Q << endl;
		_os << endl;

		_os << "Multi-target filter: " << endl;
		_os << "\tFilter type: " << _params.multipleTargetFilterType << endl;
		_os << "\tGaussian Mixture Max. Size: " << _params.gaussianMixtureMaxSize << endl;
		_os << "\tBirth type: " << _params.birthType << endl;
		_os << "\tBirth size: " << _params.birthSize << endl;
		_os << "\tMerge threshold: " << _params.mergeThreshold << endl;
		_os << "\tPrune threshold: " << _params.pruneThreshold << endl;
		_os << "\tEstimate threshold: " << _params.estimateThreshold << endl;
		_os << "\tBirth intensity: " << _params.birthIntensity << endl;
		_os << "\tProbability of birth: " << _params.pB << endl;
		_os << "\tProbability of survival: " << _params.pS << endl;
		_os << "\tInitial q: " << _params.qInit << endl;
		_os << "\tBirth covariance: " << endl << _params.birthCovariance << endl;

		_os << "\tBeta epsilon: " << _params.epsilon;

		return _os;
	}

	/*
	* <summary> A function for printing Eigen vector with single spacing between components. </summary>
	* <param name = "_os"> A reference to the output stream. </param>
	* <param name = "_v"> A constant reference to the vector to print. </param>
	*/
	void printVector(ofstream& _os, const VectorXd & _v)
	{
		for (size_t i = 0; i < (size_t)_v.size() - 1; i++)
			_os << _v(i) << ",";
		_os << _v(_v.size() - 1);
	}

	/**
	* <summary> Splits a string according to the delimeters specified. </summary>
	* Source: http://stackoverflow.com/questions/236129/split-a-string-in-c
	* <returns> A vector of strings containing the elements of the initial string </returns>
	*/
	vector<string>& split(const string &s, char delim, vector<string> &elements) {

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

	/*
	 * <summary> Reads .netcdf processed entry from the .csv file . </summary>
	 * <param name = "_input"> Input .csv file. </param>
	 * <returns> A VectorXd containing entry values. </returns>
 	 */
	VectorXd readNetCDFProcessedEntry(ifstream& _input)
	{
		string buffer;
		getline(_input, buffer);
		vector<string> elements;
		split(buffer, ',', elements);
		VectorXd result(elements.size());

		if (elements.size() != 0)
		{
			for (size_t i = 0; i < elements.size(); i++)
				result[i] = stod(elements[i]);
		}

		return result;
	}

	/*
	 * <summary> Calculates the elapsed time in seconds between two clocks. </summary>
	 * <param name = "_start"> Start clock. </param>
	 * <param name = "_end"> End clock. </param>
	 * <return> Time diffrence between end and start in seconds. </return>
	 */
	double elapsedSeconds(const double& _start, const double& _end) 
	{ 
		return (_end - _start) / CLOCKS_PER_SEC; 
	}
}
