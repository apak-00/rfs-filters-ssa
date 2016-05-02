#include <iostream>
#include <Eigen/Core>

#include "Astro.h"
#include "GMPHDFilter.h"
#include "GMJoTTFilter.hpp"
#include "bGMJoTTFilter.hpp"
#include "TDMReader.h"

#include <fstream>
#include <iomanip>

#include <dirent.h>


using namespace std;
using namespace Eigen;

// Function for JoTT testing
template <typename T>
void testGMJoTT(const T& _filter, const size_t& _sDim, const size_t& _zDim, const double& _dt, const double& _pS, const double& _pB,
	const double& _q, const unsigned int& _nBirth, const double& _bIntensity, const string& _outputFileName);

template <typename T>
void testbGMJoTT(const T& _filter, const size_t& _sDim, const size_t& _zDim, const double& _dt, const double& _pS, const double& _pB,
	const double& _q, const unsigned int& _nBirth, const double& _bIntensity, const string& _outputFileName, const double& _epsilon);

void testTransformations();
void printVector(ofstream& _os, const VectorXd& v);

template<typename T>
void printMixtureAstroRAZEL(const T& _mixture, const Sensor& _sensor, const bool& _beta);
void prepareVariables(const size_t& _sDim, const size_t& _zDim, const double& _dt, const size_t& _gmmSize);
void testSTM();
void tdmToCSV(const string& _s);

vector<string>& split(const string &s, char delim, vector<string> &elements);
MatrixXd getInitCovCV(const size_t& _dim, const double& _pCov, const double& _vCov);
VectorXd readSimData();
VectorXd readAdditionalClutterMeasurements();
bool has_suffix(const string& s, const string& suffix);

Sensor sensor;
gaussian_mixture gmm;
beta_gaussian_mixture bgmm;
string filename, filenameSim, filenameClutter;
ifstream inputSim, inputClutter;
VectorXd lBound, uBound;
bool sim;
string buffer;

// 7000 4643 7088 6000+(Sim)
// Test
int main(int arcg, char** argv) 
{
	sim = false;
	
	filename = "dataset_2.tdm";
	filenameSim = "iss_sim.txt";
	
	double dt = dt = 0.0704;
	size_t sDim = 6, zDim = 3;

	MatrixXd F = KalmanFilter::getCVF(sDim, dt), Q = KalmanFilter::getCVQ(sDim, dt) * 0.01; // 1000 2500 good
	KalmanFilter kf(F, Q, dt);
	ExtendedKalmanFilter ekf(F, Q, dt);

	// sDim, zDim, dt, pS, pB, q, nB, bI, file
	testGMJoTT(ekf, sDim, zDim, dt, 1, 0.5, 0.01, 1, 1, "output_gmjott_6d.txt");

	//testbGMJoTT(ekf, sDim, zDim, dt, 1, 1e-5, 0, 1, 1, "output_bgmjott_6d.txt", 0.2);

	return 0; 
}

/**
 * <summary> Prepares variables </summary>
 */
void prepareVariables(const size_t& _sDim, const size_t& _zDim, const double& _dt, const size_t& _gmmSize)
{
	// Parameter intialization

	// Sensor, observation matrices
	MatrixXd R = MatrixXd::Identity(_zDim, _zDim), H = MatrixXd::Zero(_zDim, _sDim);

	if (_sDim == 2 && _zDim == 1) 
	{		// Range-only linear filtering
		 
		R << 1;										
		H << 1, 0;									
	}
	else if (_sDim == 6 && _zDim == 3) 
	{ // ECI state with RAZEL measurements
	
		R(0, 0) = 0.5;			// Range covariance 0.75
		R(1, 1) = 0.075;		// Azimuth covariance	0.075
		R(2, 2) = 0.075;		// Elevation covariance	0.075
		R = R * 1;		// 12
	}
		
	double pD = 0.8, lambda = 1, V = 1e-6;		// Probability of detection and clutter
	sensor =  Sensor(_zDim, _sDim, pD, lambda, V, R, H);

	VectorXd pos(3); 
	pos << 51.1483578, -1.4384458, 0.081;
	sensor.setPosition(pos);

	// GMM
	gmm = gaussian_mixture(_sDim, _gmmSize);
	bgmm = beta_gaussian_mixture(_sDim, _gmmSize);
}

/**
* <summary> GMJoTT filter test. </summary>
*/
template<typename T>
void testGMJoTT(const T& filter, const size_t& _sDim, const size_t& _zDim, const double& _dt, const double& _pS, const double& _pB, 
	const double& _q, const unsigned int& _nBirth, const double& _bIntensity, 
	const string& _outputFileName)
{
	VectorXd info, bearing(2), dummy = VectorXd::Zero(_sDim);
	vector<VectorXd> measurements;
	vector<gaussian_component> estimates;
	ofstream outputFile(_outputFileName); 				// Output file 
	MatrixXd iCov = getInitCovCV(_sDim, 3600, 10);		// Initial covariance matrix
	TDMReader reader(filename);							// TDM Reader for CAMRa

	if (sim)
		inputSim.open(filenameSim);

	// Misc
	bool dtCalculated = false;
	Astro::date datePrev, dateCurr;
	double dt = 0, deg2rad = M_PI / 180.0;

	// Additional measurements - added 30/03/2016
	string filenameClutter = "rand_meas.txt";
	inputClutter.open(filenameClutter);
	bool additionalClutter = false;
	VectorXd clutterMeasurement(3);

	// State, observation's dimensions and intial timestep
	prepareVariables(_sDim, _zDim, _dt, 100);

	GMJoTTFilter<ExtendedKalmanFilter> gmjottfilter(filter, _nBirth, _bIntensity, _pS, iCov, lBound, uBound, _q, _pB);

	// Main loop ------------------------------------------------------------------------
	for (size_t i = 0; i < 10000; i++)
	{
		if (sim) 
			info = readSimData();
		else 
			info = reader.readEntryEigen();

		if (!info.size())
			break;
		
		if (info(1) > 360.0 )
			info(1) -= 360.0;
		info(1) = 180.0 - info(1);		// Azimuth correction (?)

		dateCurr = Astro::date((int)info(5), (int)info(6), (int)info(7),
			(int)info(8), (int)info(9), info(10));

		if (!dtCalculated)
		{
			datePrev = dateCurr;
			dtCalculated = true;
		}
		else
		{
			dt = dateCurr - datePrev;
			gmjottfilter.updateKFTimestep(dt);
			datePrev = dateCurr;
		}

		measurements.clear();
		measurements.push_back(info.segment(0, _zDim));
		bearing << info(1), info(2);

		// Addtioinal clutter test
		if (additionalClutter) 
		{
			VectorXd ranges = readAdditionalClutterMeasurements();
			
			for (size_t j = 0; j < (size_t)ranges.size(); j++)
			{
				clutterMeasurement << ranges(j), bearing;
				measurements.push_back(clutterMeasurement);
			}
		}
		
		sensor.setZ(measurements);
		sensor.setBearing(bearing);

		// JoTT Filter
		//cout << endl << "\tStep: " << i << "\t " << r.transpose() << " " << "--------" << endl;
		//cout << "\t\t" << Astro::razelToTEME(r, sensor.getPosition(), sensor.getDateJD(), sensor.getLOD(), 
			//sensor.getXp(), sensor.getYp()).transpose() << endl;

		bool c = false;

		// JoTT Filter
		if (c)
			cout << "\tSTEP " << i << " ----- " << info(0) << " " << info(1) << " " << info(2) << " -----" << endl;

		// Predict
		gmjottfilter.predict(gmm, sensor);
		if (c)
		{
			cout << "\tPREDICT: q = " << gmjottfilter.getQ() << "\tEw = " << gmm.weightSum() << endl;
			printMixtureAstroRAZEL(gmm, sensor, false);
		}

		// Update
		gmjottfilter.update(gmm, sensor);
		if (c)
		{
			cout << "\tUPDATE: q = " << gmjottfilter.getQ() << "\tEw = " << gmm.weightSum() << endl;
			printMixtureAstroRAZEL(gmm, sensor, false);
		}

		if (c)
			cout << "\tMERGE: " << gmm.size() << "->";
		// Merge
		gmm.merge(0.8);
		if (c)
		{
			cout << gmm.size() << "\tEw = " << gmm.weightSum() << endl;
			printMixtureAstroRAZEL(bgmm, sensor, true);
		}

		if (c)
			cout << "\tPRUNE: " << gmm.size() << "->";
		// Prune
		gmm.prune(1e-8);
		if (c)
		{
			cout << gmm.size() << "\tEw = " << gmm.weightSum() << endl;
			printMixtureAstroRAZEL(gmm, sensor, true);
		}

		estimates = gmm.getEstimates(0.5);

		// File stream, Matlab output ----------
		// Measurement
		outputFile << i << " ";
		printVector(outputFile, info.segment(0, _zDim));
		
		// Estimates - ECEF
		outputFile << estimates.size() << " " << gmjottfilter.getQ() << " ";
		if (!estimates.empty()) 
			printVector(outputFile, Astro::temeToECEF(estimates[0].m, sensor.getDateJD(), 0, 0, 0));
		else 
			printVector(outputFile, dummy);
			
		// Date
		printVector(outputFile, info.segment(5, 6));

		// Estimates - RAZEL
		if (!estimates.empty()) 
			printVector(outputFile, Astro::temeToRAZEL(estimates[0].m, sensor.getPosition(), sensor.getDateJD(), 0, 0, 0));
		else
			printVector(outputFile, dummy);

		// Covariance
		if (!estimates.empty()) 
		{
			printVector(outputFile, estimates[0].P.diagonal());

			outputFile << estimates[0].P.determinant() << " "
				<< estimates[0].P.block<3, 3>(0, 0).determinant() << " "
				<< estimates[0].P.block<3, 3>(3, 3).determinant() << " ";
		}
		else 
		{
			printVector(outputFile, dummy);
			outputFile << "0 0 0 ";
		}
		
		if (!estimates.empty()) 
			outputFile << estimates[0].tag[1];
		else
			outputFile << "-1";

		outputFile << endl;

		// Console output
		cout << i << " " << estimates.size() << "/" << gmm.size() << " "; 

		if (!estimates.empty()) {
			cout << setprecision(2) << estimates[0].w << " [" << estimates[0].tag[1] << "][" << estimates[0].tag[2] << "] ";
			if (estimates[0].kindaConverged)
				cout << "S ";
			else
				cout << "E ";

			if (estimates[0].kindaConverged)
				cout << "Converged" << endl;
		}

		cout << endl;
	}

	outputFile.close();
	inputClutter.close();
}

template<typename T>
void testbGMJoTT(const T & _filter, const size_t & _sDim, const size_t & _zDim, const double & _dt, const double & _pS, const double & _pB, 
	const double & _q, const unsigned int & _nBirth, const double & _bIntensity, const string & _outputFileName, const double & _epsilon)
{
	VectorXd info, bearing(2), dummy = VectorXd::Zero(_sDim);
	vector<VectorXd> measurements;
	vector<beta_gaussian_component> estimates;
	ofstream outputFile(_outputFileName); 				// Output file 
	MatrixXd iCov = getInitCovCV(_sDim, 3600, 10);		// Initial covariance matrix
	TDMReader reader(filename);							// TDM Reader for CAMRa

	if (sim)
		inputSim.open(filenameSim);

	// Misc
	bool dtCalculated = false;
	Astro::date datePrev, dateCurr;
	double dt = 0, deg2rad = M_PI / 180.0;

	// Additional measurements - added 30/03/2016
	string filenameClutter = "rand_meas.txt";
	inputClutter.open(filenameClutter);
	bool additionalClutter = false;
	VectorXd clutterMeasurement(3);

	// State, observation's dimensions and intial timestep
	prepareVariables(_sDim, _zDim, _dt, 20);

	bGMJoTTFilter<ExtendedKalmanFilter> bgmjottfilter(_filter, _nBirth, _bIntensity, _pS, iCov, lBound, uBound, _q, _pB, _epsilon);

	// Main loop ------------------------------------------------------------------------
	for (size_t i = 0; i < 10000; i++)
	{
		if (sim)
			info = readSimData();
		else
			info = reader.readEntryEigen();

		if (!info.size())
			break;

		if (info(1) > 360.0)
			info(1) -= 360.0;
		info(1) = 180.0 - info(1);		// Azimuth correction (?)

		dateCurr = Astro::date((int)info(5), (int)info(6), (int)info(7),
			(int)info(8), (int)info(9), info(10));

		if (!dtCalculated)
		{
			datePrev = dateCurr;
			dtCalculated = true;
		}
		else
		{
			dt = dateCurr - datePrev;
			bgmjottfilter.updateKFTimestep(dt);
			datePrev = dateCurr;
		}

		measurements.clear();
		measurements.push_back(info.segment(0, _zDim));
		bearing << info(1), info(2);

		// Addtioinal clutter test
		if (additionalClutter)
		{
			VectorXd ranges = readAdditionalClutterMeasurements();

			for (size_t j = 0; j < (size_t)ranges.size(); j++)
			{
				clutterMeasurement << ranges(j), bearing;
				measurements.push_back(clutterMeasurement);
			}
		}

		sensor.setZ(measurements);
		sensor.setBearing(bearing);
		 
		bool c = false;

		if (i >= 3707)
		{
			c = true;
		}
			
		// JoTT Filter
		if (c)
			cout  << "\tSTEP " << i << " ----- " << info(0) << " " << info(1) << " " << info(2) << " -----" << endl;

		// Predict
		bgmjottfilter.predict(bgmm, sensor);
		if (c)
		{
			cout << "\tPREDICT: q = " << bgmjottfilter.getQ() << "\tEw = " << bgmm.weightSum() << endl;
			printMixtureAstroRAZEL(bgmm, sensor, true);
		}

		// Update
		bgmjottfilter.update(bgmm, sensor);
		if (c)
		{
			cout << "\tUPDATE: q = " << bgmjottfilter.getQ() << "\tEw = " << bgmm.weightSum() << endl;
			printMixtureAstroRAZEL(bgmm, sensor, true);
		}
		
		if (c)
			cout << "\tMERGE: " << bgmm.size() << "->";
		// Merge
		bgmm.merge(0.8);
		if (c)
		{
			cout << bgmm.size() << "\tEw = " << bgmm.weightSum() << endl;
			printMixtureAstroRAZEL(bgmm, sensor, true);
		}

		if (c)
			cout << "\tPRUNE: " << bgmm.size() << "->";
		// Prune
		bgmm.prune(1e-4);
		if (c) 
		{
			cout << bgmm.size() << "\tEw = " << bgmm.weightSum() << endl;
			printMixtureAstroRAZEL(bgmm, sensor, true);
		}
		
		estimates = bgmm.getEstimates(0.5);

		// File stream, Matlab output ----------
		// Measurement
		outputFile << i << " ";
		printVector(outputFile, info.segment(0, _zDim));

		// Estimates - ECEF
		outputFile << estimates.size() << " " << bgmjottfilter.getQ() << " ";
		if (!estimates.empty())
			printVector(outputFile, Astro::temeToECEF(estimates[0].m, sensor.getDateJD(), 0, 0, 0));
		else
			printVector(outputFile, dummy);

		// Date
		printVector(outputFile, info.segment(5, 6));

		// Estimates - RAZEL
		if (!estimates.empty())
			printVector(outputFile, Astro::temeToRAZEL(estimates[0].m, sensor.getPosition(), sensor.getDateJD(), 0, 0, 0));
		else
			printVector(outputFile, dummy);

		// Covariance
		if (!estimates.empty())
		{
			printVector(outputFile, estimates[0].P.diagonal());

			outputFile << estimates[0].P.determinant() << " "
				<< estimates[0].P.block<3, 3>(0, 0).determinant() << " "
				<< estimates[0].P.block<3, 3>(3, 3).determinant() << " ";
		}
		else
		{
			printVector(outputFile, dummy);
			outputFile << "0 0 0 ";
		}

		if (!estimates.empty())
			outputFile << estimates[0].tag[1] << " ";
		else 
			outputFile << "-1 ";

		if (!estimates.empty())
			outputFile << estimates[0].u / (estimates[0].u + estimates[0].w);
		else
			outputFile << "-1 ";

		outputFile << endl;

		// Console output
		cout << i << " [" << estimates.size() << "/" << bgmm.size() << "] " << setprecision(4) << info(0) << " q = " << setprecision(4) << bgmjottfilter.getQ() << " ";

		if (c)
			cout << " Estimates: " << estimates.size() << "/" << bgmm.size() << endl;

		if (!estimates.empty()) { 
			
			auto meanRAZEL = Astro::temeToRAZEL(estimates[0].m, sensor.getPosition(), sensor.getDateJD(), sensor.getLOD(), sensor.getXp(), sensor.getYp());
			
			cout << setprecision(2)  << "[" << estimates[0].tag[0] << "][" << estimates[0].tag[1] << "][" << estimates[0].tag[2] << "] " 
				<< estimates[0].w << setprecision(4) << " ["<< meanRAZEL(0) << "]";
			cout << " B: " << estimates[0].u / (estimates[0].u + estimates[0].v);
			
			if (estimates[0].kindaConverged)
				cout << " S";
			else
				cout << " E";

			if (estimates[0].kindaConverged)
				cout << "Converged" << endl;
		}

		cout << endl;
	}

	outputFile.close();
	inputClutter.close();
}

template<typename T>
void printMixtureAstroRAZEL(const T & _mixture, const Sensor & _sensor, const bool& _beta)
{
	cout << "\tMixture: [" << _mixture.size() << "/" << _mixture.nMax << "]" << endl;

	for (auto c : _mixture.components) 
	{
		auto meanRAZEL = Astro::temeToRAZEL(c.m, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
		cout << "[" << c.tag[0] << "][" << c.tag[1] << "][" << c.tag[2] << "] w: " << c.w
			<< "\t m: " << meanRAZEL(0) << " " << meanRAZEL(1) << " " << meanRAZEL(2) << " ";

		// For BGM
		//if (_beta)
		//	cout << "Beta: [" << c.u / (c.u + c.v) << "]";

		cout << endl;
	} 
}

/**
* ///////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

/**
 * <summary> Checks if the specified string contains the desired character sequence (suffix). </summary>
 * <param name = "_s"> The main string to be searched. </param>
 * <param name = "_suffix"> The required substring. </param>
 */
bool has_suffix(const string& _s, const string& _suffix)
{
	return (_s.size() >= _suffix.size()) && equal(_suffix.rbegin(), _suffix.rend(), _s.rbegin());
}

/**
 * <summary> A function for testing Shepperd state transition matrix. </summary>
 */
void testSTM() {

	ofstream stm;
	stm.open("test_stm.txt");

	VectorXd state(6), result(6), resultF(6); 
	MatrixXd shepperd;
	state << -1288.75878841, 4184.25775350, 5168.32743761, -7.52063480, -0.61162224, -1.37936734;

	shepperd = Astro::getShepperdMatrix(state, 10.0, result, Astro::MU_E);
	resultF = shepperd * state;

	cout << setprecision(8) << state.transpose() << endl;
	cout << result.transpose() << endl;
	cout << resultF.transpose() << endl;

	cout << endl;

	stm.close();
}

/*
 * <summary> A function for testing coordinate transformations. </summary>
 */
void testTransformations()
{
	double deg2rad = M_PI / 180.0;

	VectorXd razel(6), sez(6), razelRad(6), teme(6);
	Astro::date d(2014, 2, 28, 10, 22, 3.65);
	double jd = Astro::getJulianDay(d);

	razel << 638.615, 180 - 168.435, 16.1771, 1, 0.1, 0.2;
	teme << 3811.41, 413.827, 5739.91, 0.257854, -0.743775, -1.24866;

	cout << "TEME:\t" << teme.transpose() << endl;
	VectorXd ecef = Astro::temeToECEF(teme, jd, 0, 0, 0);
	cout << "ECEF:\t" << ecef.transpose() << endl;

}

/*
 * <summary> A function for printing Eigen vector with single spacing between components. </summary>
 * <param name = "_os"> A reference to the output stream. </param>
 * <param name = "_v"> A constant reference to the vector to print. </param>
 */
void printVector(ofstream& _os, const VectorXd & _v)
{
	for (size_t i = 0; i < (size_t)_v.size(); i++)
		_os << _v(i) << " ";
}

/**
 * <summary> Uniform birth process for range filtering </summary>
 * <param name = "_sensor"> A refenrce to the sensor. </param>
 * <param name = "_sensor"> Measurement number </param>
 * <returns> A 2-dimensional VectorXd containing range = 1000 and zero velocity. </returns>
 */
VectorXd uniformBirthRange2D(Sensor& _sensor, const size_t& _zNum)
{
	VectorXd birth(2);
	birth << 1000, 0;
	return birth;
}

/**
 * <summary> Uniform birth process for 6D filtering in ECI. </summary>
 * <param name = "_sensor"> A reference to the sensor. </param>
 * <param name = "_zNum"> Mesurement number. </param>
 * <returns> A VectorXd containing a ECI state vector with range = 1000 km along the current sensor bearing. </returns>
 */
VectorXd uniformBirthRange6DRAZEL(Sensor & _sensor, const size_t& _zNum)
{
	VectorXd birth = _sensor.getZ(_zNum);
	birth(0) = 1000;
	return Astro::razelToTEME(birth, _sensor.getPosition(), Astro::getJulianDay(_sensor.getDate()), _sensor.getXp(), sensor.getYp(), 2);
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

/**
 * <summary> Reads simulated data. </summary>
 * <returns> A VectorXd with the simulated rada measurement and corresponding timestamp. </returns>
 */
VectorXd readSimData()
{
	getline(inputSim, buffer);
	vector<string> elements;
	split(buffer, ' ', elements);
	VectorXd result(elements.size());

	if (elements.size() != 0)
	{
		for (size_t i = 0; i < elements.size(); i++)
			result[i] = stod(elements[i]);
	}

	return result;
}

/**
 * <summary> A function for reading additional clutter measurements generated in Matlab. </summary>
 * <retuns> A VectorXd containing clutter measurements. </returns>
 */
VectorXd readAdditionalClutterMeasurements()
{
	getline(inputClutter, buffer);
	vector<string> elements;
	split(buffer, ' ', elements);
	VectorXd result(elements.size());

	if (elements.size() != 0)
	{
		for (size_t i = 0; i < elements.size(); i++)
			result[i] = stod(elements[i]);
	}

	return result;
}

/**
 * <summary> A function for converting all files with .tdm extension to .csv file with raw data. </summary>
 * <param name = "_s"> A path to the folder containing .tdm files. </param>
 */
void tdmToCSV(const string & _s)
{
	DIR *dir = opendir(_s.c_str());
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
				filePath = _s + "\\" + entry->d_name;
				TDMReader tdmr(filePath);

				filePath.replace(filePath.end() - 3, filePath.end(), "csv");
				ofscsv.open(filePath);

				tdmr.getHeader();

				while (true) {

					dataEntry = tdmr.readEntry();

					if (!dataEntry.size())
						break;

					for (size_t i = 0; i < dataEntry.size()-1; i++) 
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
 * <summary> Gets the inital covariance for the constant velocity model. </summary>
 * <param name = "_dim"> Dimension of the state vector. </param>
 * <param name = "_pCov"> Covariance of the position. </param>
 * <param name = "_vCov"> Covariance of the velocity </param>
 * <returns> MatrixXd containing the initial covariance matrix (constant velocity model). </returns>
 */
MatrixXd getInitCovCV(const size_t & _dim, const double & _pCov, const double & _vCov)
{
	assert(_dim % 2 == 0 && "Works only for even-sized matrices");
	MatrixXd initCov(_dim, _dim), i = MatrixXd::Identity(_dim / 2, _dim / 2);
	initCov << i * _pCov, i * 0, i * 0, i * _vCov;
	return initCov;
}