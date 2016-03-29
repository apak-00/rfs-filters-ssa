#include <iostream>
#include <Eigen/Core>

#include "Astro.h"
#include "GMPHDFilter.h"
#include "GMJoTTFilter.hpp"
#include "TDMReader.h"

#include <fstream>
#include <iomanip>

using namespace std;
using namespace Eigen;

template <typename T>
void testGMJoTT(const T& _filter, const size_t& _sDim, const size_t& _zDim, const double& _dt, const double& _pS, const double& _pB,
	const double& _q, const unsigned int& _nBirth, const double& _bIntensity, const string& _outputFileName);

void testTransformations();
void printVector(ofstream& _os, const VectorXd& v);
void prepareVariables(const size_t& _sDim, const size_t& _zDim, const double& _dt, const size_t& _gmmSize);
void testSTM();
vector<string>& split(const string &s, char delim, vector<string> &elements);

MatrixXd getInitCovCV(const size_t& _dim, const double& _pCov, const double& _vCov);
VectorXd readSimData();

Sensor sensor;
gaussian_mixture gmm;
string filename, filenameSim;
ifstream inputSim;
VectorXd lBound, uBound;
bool sim;
string buffer;

// 7000 4643 7088 6000+(Sim)

int main(int arcg, char** argv) 
{
	sim = false;
	
	filename = "dataset_3.tdm";
	filenameSim = "iss_sim.txt";
	
	double dt = dt = 0.0704;
	size_t sDim = 6, zDim = 3;

	MatrixXd F = KalmanFilter::getCVF(sDim, dt), Q = KalmanFilter::getCVQ(sDim, dt) * 0.01;
	KalmanFilter kf(F, Q, dt);
	ExtendedKalmanFilter ekf(F, Q, dt);

	// sDim, zDim, dt, pS, pB, q, nB, bI, file
	testGMJoTT(ekf, sDim, zDim, dt, 1, 0.5, 0.01, 1, 1, "output_gmjott_6d.txt");

	//testSTM();

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
	
		R(0, 0) = 0.5;		// Range covariance
		R(1, 1) = 0.075;		// Azimuth covariance
		R(2, 2) = 0.075;		// Elevation covariance
		R *= 1;
	}
		
	double pD = 0.8, lambda = 1, V = 1e-6;		// Probability of detection and clutter
	sensor =  Sensor(_zDim, _sDim, pD, lambda, V, R, H);

	VectorXd pos(3); 
	pos << 51.1483578, -1.4384458, 0.081;
	sensor.setPosition(pos);

	// GMM
	gmm = gaussian_mixture(_sDim, _gmmSize);
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
		measurements.push_back(info.segment(0,_zDim));
		bearing << info(1), info(2);
		sensor.setZ(measurements);
		sensor.setBearing(bearing);

		// JoTT Filter
		//cout << endl << "\tStep: " << i << "\t " << r.transpose() << " " << "--------" << endl;
		//cout << "\t\t" << Astro::razelToTEME(r, sensor.getPosition(), sensor.getDateJD(), sensor.getLOD(), 
			//sensor.getXp(), sensor.getYp()).transpose() << endl;

		gmjottfilter.predict(gmm, sensor);
		//cout << "\tPredict: q = " << gmjottfilter.getQ() << endl << gmm;
		gmjottfilter.update(gmm, sensor);
		//cout << "\tUpdate: q = " << gmjottfilter.getQ() << endl << gmm;	
		gmm.merge(0.8);
		//cout << "\tMerge:" << endl << gmm;
		gmm.prune(1e-8);
		//cout << "\tPrune:" << endl << gmm;

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
		printVector(outputFile, info.segment(4, 6));

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
}

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

	/*
	for (size_t i = 10; i <= 110; i += 10) {
		shepperd =  Astro::getShepperdMatrix(state, (double)i, result, Astro::MU_E);
		cout << result.transpose() << endl;
		printVector(stm, result);
		stm << endl;
	}
	*/

	stm.close();
}

/**
* ///////////////////////////////////////////////////////////////////////////////////////////////////////////
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

void printVector(ofstream& _os, const VectorXd & v)
{
	for (size_t i = 0; i < v.size(); i++)
		_os << v(i) << " ";
}

/**
* <summary> Uniform birth process for range filtering </summary>
*/
VectorXd uniformBirthRange2D(Sensor& _sensor, const size_t& _zNum)
{
	VectorXd birth(2);
	birth << 1000, 0;

	return birth;
}

/**
* <summary> Uniform birth process for 6D filtering in ECI. </summary>
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
* <return> A vector of strings containing the elements of the initial string </return>
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
* <summary> Gets the inital covariance for the constant velocity model. </summary>
*/
MatrixXd getInitCovCV(const size_t & _dim, const double & _pCov, const double & _vCov)
{
	assert(_dim % 2 == 0 && "Works only for even-sized matrices");
	MatrixXd initCov(_dim, _dim), i = MatrixXd::Identity(_dim / 2, _dim / 2);
	initCov << i * _pCov, i * 0, i * 0, i * _vCov;
	return initCov;
}