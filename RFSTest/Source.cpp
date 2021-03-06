#include <iostream>
#include <fstream>
#include <iomanip>		
#include <memory>
#include <Eigen/Core>
#include <yaml-cpp/yaml.h>

#include <ctime> 
#include <numeric>

#include "Debug.h"
#include "Astro.h"
#include "GMPHDFilter.h"
#include "GMJoTTFilter.hpp"
#include "BGMJoTTFilter.hpp"
#include "SMCJoTTFilter.hpp"

#include "ExtendedKalmanFilter.h"
#include "UnscentedKalmanFilter.h"

#include "IOHelpers.h"

using namespace std;
using namespace Eigen;
using namespace IOHelpers;

const double DEG2RAD = M_PI / 180.0;

// JoTT testing
void testFilter(parameters& _p);

// GM JoTT Filter
template <typename MultiTargetFilter, typename Mixture>
void runFilter(MultiTargetFilter& _filter, Sensor& _sensor, Mixture& _mixture, parameters& _p);

// SMC JoTT Filter
template <typename MultiTargetFilter, typename Swarm>
void runSMCJoTTFilter(MultiTargetFilter& _filter, Sensor& _sensor, Swarm& _mixture, parameters& _p);

// PHD testing
template <typename MultiTargetFilter, typename Mixture>
void runPHDFilter(MultiTargetFilter& _filter, Sensor& _sensor, Mixture& _mixture, parameters& _p);

// Temporary testing functions
void testTransformations();
void testSingleTargetFilter(parameters & _p);

// Temporary observation function
VectorXd observeTEMEToRAZEL(const VectorXd& _teme, const Sensor& _sensor)
{
	return Astro::temeToRAZEL(_teme, _sensor.getPosition(), _sensor.getDateJD(),
		_sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
}

// Main
int main(int arcg, char** argv) 
{	
	Astro::date d(2014, 02, 28, 10, 20, 0.899995);
	auto jd = Astro::getJulianDay(d);
	auto gmst = Astro::getGMST(jd);
	cout.precision(16);
	cout << "JD: " << jd << endl;
	cout << "GMST: " << gmst << endl;

	VectorXd eci(6), gd(3);
	eci << 2190.6626281, -2131.7118564, 6188.48993845, 6.22542053, -2.93149252, -3.20660876;
	gd << 51.1483578  * M_PI / 180.0, -1.4384458 * M_PI / 180.0, 0.081;
	auto ecef = Astro::temeToECEF(eci, jd, 0, 0, 0);
	auto gdecef = Astro::geodeticToECEF(gd);
	auto sez = Astro::ecefToSEZ(ecef, gd);
	auto razel = Astro::sezToRAZEL(sez);

	cout << "ECI:" << eci.transpose() << endl;
	cout << "ECEF:" << ecef.transpose() << endl;
	cout << "GD: " << gdecef.transpose() << endl;
	cout << "SEZ: "  << sez.transpose() << endl;
	cout << "RAZEL: " << razel.transpose() << endl;




	//string filename_params = "config/config-smcjott-2018.yaml";
	//parameters p = readParametersYAML(filename_params);
	//testFilter(p);
	//testSingleTargetFilter(p);
 	return 0;
}

/*
 * <summary> Tests single-target filters (KF, EKF, UKF) </summary>
 * <param name = "_p"> YAML parameters. </param>
 */
void testSingleTargetFilter(parameters & _p) {

	VectorXd entry, measurement(3), temp(6), bearing(2), sensorPos;
	Astro::date dateCurr, datePrev;
	double dt, azRad, elRad;
	shared_ptr<KalmanFilter> kf;
	ifstream input;
	ofstream output;
	bool debug_ = true;

	// Filter choice
	if (!strcmp(_p.singleTargetFilterType.c_str(), "kf"))
		kf = make_shared<KalmanFilter>(_p.F, _p.Q * 0.1, _p.dt);
	else if (!strcmp(_p.singleTargetFilterType.c_str(), "ekf"))
		kf = make_shared<ExtendedKalmanFilter>(_p.F, _p.Q * 0.1, _p.dt);
	else if (!strcmp(_p.singleTargetFilterType.c_str(), "ukf"))
	{
		kf = make_shared<UnscentedKalmanFilter>(_p.Q * 0.1, _p.ukfSigmaSamplerW, 0);
		kf->setT(_p.dt);
	} 

	// Sensor 
	Sensor sensor = Sensor(_p.observationDim, _p.stateDim, _p.pD, _p.lambda, 1e-6, _p.R, _p.H);
	sensor.setPosition(_p.sensorPosition);
	// Additional reference system-related parameters
	sensor.setXp(0);
	sensor.setYp(0);
	sensor.setLOD(0);

	// Temporary
	_p.R(1, 1) = _p.R(1, 1) * DEG2RAD;
	_p.R(2, 2) = _p.R(2, 2) * DEG2RAD;

	sensorPos = sensor.getPosition();
	auto measurements = sensor.getZ();

	input.open("Data/sim_input.csv");

	if (!input.is_open())
		cout << "Can't open input file!" << endl;

	output.open("Results/sim_output.csv");

	// Read first entry
	entry = readSimulatedEntry(input);
	datePrev = Astro::date((int)entry(0), (int)entry(1), (int)entry(2), (int)entry(3), (int)entry(4), entry(5));
	sensor.setDate(datePrev);
	azRad = entry(7);
	elRad = entry(8);
	measurement << entry(6), azRad, elRad;		// Simulated angles are already in radians

	cout << "First entry: " << entry.transpose() << endl << endl;
	cout << "Sensor position: " << "Lat: " << sensorPos(0) / DEG2RAD	// Latitude and longitude are in radians too
		<< " Lon: " << sensorPos(1) / DEG2RAD
		<< " Alt: " << sensorPos(2) << endl;
	cout << "First measurement (rad): " << measurement.transpose() << endl;

	temp << Astro::razelToTEME(measurement, sensor.getPosition(), sensor.getDateJD(), sensor.getLOD(), sensor.getXp(), sensor.getYp()), 0, 0, 0;

	cout << "Init: " << temp.transpose() << endl;
	cout << endl;

	// Initiate the gaussian component
	gaussian_component gc(temp, _p.birthCovariance, 0, 0);

	// Main loop
	for (size_t i = 0; i < 1e4 ; i++) 
	{
		// Read entry
		entry = readSimulatedEntry(input);

		if (!entry.size())
			break;

		// Date
		dateCurr = Astro::date((int)entry(0), (int)entry(1), (int)entry(2), (int)entry(3), (int)entry(4), entry(5));
		dt = dateCurr - datePrev;
		kf->setT(dt);
		datePrev = dateCurr;

		// Measurements
		bearing << entry(7), entry(8);
		measurement = entry.segment(6, _p.observationDim);
		measurements.clear();
		measurements.push_back(measurement);

		sensor.setDate(dateCurr);
		sensor.setBearing(bearing);
		sensor.setZ(measurements);

		if (i >= 0) 
		{
			debug_ = false;
			kf->debug = true;
		}

		// Filter
		kf->predict(gc);
		if (debug_)
			cout << "Predict: " << gc.m.transpose() << endl;
		kf->update(gc, sensor, 0);
		if (debug_)
			cout << "Update: " << gc.m.transpose() << endl;

		temp = Astro::temeToRAZEL(gc.m, sensor.getPosition(), sensor.getDateJD(), sensor.getLOD(), sensor.getXp(), sensor.getYp());

		cout << i << " " << temp.head(3).transpose() << " | " << measurements[0].transpose() << endl;
		output << dateCurr << ",";
		printVector(output, temp);
		output << ",";
		printVector(output, gc.m);
		output << "," << gc.P.determinant() << "," << gc.P.block<3,3>(0,0).determinant() << "," << gc.P.block<3,3>(3,3).determinant() << endl;
	}

	input.close();
	output.close();
}

/*
* <summary> Tests JoTT filters (GM-JoTT, BGM-JoTT, r-JoTT) </summary>
* <param name = "_p"> YAML parameters. </param>
*/
void testFilter(parameters & _p)
{
#ifdef MY_DEBUG
	clock_t begin = clock();
#endif
	shared_ptr<KalmanFilter> kf;

	// Choose the single-target filter
	if (!strcmp(_p.singleTargetFilterType.c_str(), "kf"))
		kf = make_shared<KalmanFilter>(_p.F, _p.Q * 0.1, _p.dt);
	else if (!strcmp(_p.singleTargetFilterType.c_str(), "ekf"))
		kf = make_shared<ExtendedKalmanFilter>(_p.F, _p.Q * 0.1, _p.dt);
	else if (!strcmp(_p.singleTargetFilterType.c_str(), "ukf"))
	{
		kf = make_shared<UnscentedKalmanFilter>(_p.Q * 0.1,  _p.ukfSigmaSamplerW, 0);
		kf->setT(_p.dt);
	}
		
	// Sensor 
	Sensor sensor = Sensor(_p.observationDim, _p.stateDim, _p.pD, _p.lambda, 1e-6, _p.R, _p.H);
	sensor.setPosition(_p.sensorPosition);
	// Additional reference system-related parameters
	sensor.setXp(0);
	sensor.setYp(0);
	sensor.setLOD(0);

	// Temporary
	_p.R(1, 1) = _p.R(1, 1) * DEG2RAD * DEG2RAD;
	_p.R(2, 2) = _p.R(2, 2) * DEG2RAD * DEG2RAD;

#ifdef MY_DEBUG
	clock_t end = clock();
	double elapsed = double(end - begin) / CLOCKS_PER_SEC;
	//std::cout << "Load time:\t" << elapsed << " sec" << std::endl;
#endif

	// Choose multi-target filter
	if (!strcmp(_p.multipleTargetFilterType.c_str(), "gmjott")) 
	{
		GMJoTTFilter gmjott(kf, _p.birthSize, _p.birthIntensity, _p.pS, _p.birthCovariance,
			VectorXd::Zero(_p.stateDim), VectorXd::Zero(_p.stateDim), _p.qInit, _p.pB);
		gaussian_mixture gm(_p.stateDim, _p.gaussianMixtureMaxSize);
		runFilter(gmjott, sensor, gm, _p);
	}
	else if (!strcmp(_p.multipleTargetFilterType.c_str(), "bgmjott")) 
	{
		BGMJoTTFilter bgmjott(kf, _p.birthSize, _p.birthIntensity, _p.pS, _p.birthCovariance,
			VectorXd::Zero(_p.stateDim), VectorXd::Zero(_p.stateDim), _p.qInit, _p.pB, _p.epsilon);
		beta_gaussian_mixture bgm(_p.stateDim, _p.gaussianMixtureMaxSize);
		runFilter(bgmjott, sensor, bgm, _p);
	} 
	else if (!strcmp(_p.multipleTargetFilterType.c_str(), "smcjott"))
	{
		// Propagation function
		std::function<VectorXd(const VectorXd&, const double&, const VectorXd&)> propagate = Astro::integrationPrediction;
		std::function<VectorXd(const VectorXd&, const Sensor&)> observe = observeTEMEToRAZEL;

		SMCJoTTFilter smcjott(propagate, observe, _p.birthSize, _p.birthIntensity, _p.pS, _p.pB, _p.qInit, _p.noiseAcceleration);
		particle_mixture ps(_p.partcileSwarmSize, _p.stateDim);

		// New: 22/11/2016 Added temporary parameters for birth testing
		smcjott.setLowerBirthBoundRange(_p.lowerBirthBoundRange);
		smcjott.setUpperBirthBoundRange(_p.upperBirthBoundRange);
		smcjott.setBirthSigmaRange(_p.birthSigmaRange);

		runSMCJoTTFilter(smcjott, sensor, ps, _p);

	}
	else if (!strcmp(_p.multipleTargetFilterType.c_str(), "gmphd")) 
	{
		GMPHDFilter gmphd(kf, _p.birthSize, _p.birthIntensity, _p.pS, _p.birthCovariance,
			VectorXd::Zero(_p.stateDim), VectorXd::Zero((_p.stateDim)));
		gaussian_mixture gm(_p.stateDim, _p.gaussianMixtureMaxSize);
		runPHDFilter(gmphd, sensor, gm, _p);
	}
}

/*
* <summary> Temporary transformation testing function (RAZEL, ECEF, ECI) </summary>
*/
void testTransformations()
{
	VectorXd geo(3);
	geo << 51.1483578 * M_PI / 180.0, -1.4384458 * M_PI / 180.0, 0.081;

	VectorXd teme(6);
	double jd = Astro::getJulianDay(Astro::date(2016, 7, 20, 21, 23, 59.99999464));
	teme << -561.05679138, -5463.10658244, 4157.25287649, 1.701447528, 4.385651194, 5.979623125;
	VectorXd ecef = Astro::temeToECEF(teme, jd, 0, 0, 0);
	VectorXd sez = Astro::temeToSEZ(teme, geo, jd, 0, 0, 0);
	VectorXd razel = Astro::temeToRAZEL(teme, geo, jd, 0, 0, 0);

	cout << "TEME: " << teme.transpose() << endl << endl;

	cout << "ECEF : " << ecef.transpose() << endl;
	cout << "SEZ  : " << sez.transpose() << endl;
	cout << "RAZEL: " << razel.transpose() << endl;
	cout << endl;

	VectorXd sez1 = Astro::razelToSEZ(razel);
	VectorXd ecef1 = Astro::sezToECEF(sez1, geo);
	VectorXd teme1 = Astro::razelToTEME(razel, geo, jd, 0, 0, 0);

	cout << "-> SEZ : " << sez1.transpose() << endl;
	cout << "-> ECEF: " << ecef1.transpose() << endl;
	cout << "-> TEME: " << teme1.transpose() << endl;
}
 
/**
* <summary> JoTT filter test. GM </summary>
*/
template <typename MultiTargetFilter, typename Mixture>
void runFilter(MultiTargetFilter& _filter, Sensor& _sensor, Mixture& _mixture, parameters& _p)
{
	VectorXd info(11), infoTemp, bearing(2), temp(3);
	auto measurements = _sensor.getZ();
	auto estimates = _mixture.getEstimates(1);
	bool dtCalculated = false;
	Astro::date datePrev, dateCurr;

	double deg2rad = M_PI / 180.0;

	// IO
	TDMReader tdmReader;
	bool tdm = false, sim = false;
	//YAML::Emitter result, resultSmoothing;
	//result << YAML::BeginSeq;
	//resultSmoothing << YAML::BeginSeq;
	ofstream outputCSV("Results/result_gmjott.csv");
	//outputYAML("Results/result_gmjott.yaml"),
	//outputYAMLSmoothing("Results/result_gmjott_smoothing.yaml");		// Output file 
	ifstream input;

#ifdef MY_DEBUG
	// Debug
	vector<double> elapsed(6);
	clock_t start, end;
	double elapsedTotal;
#endif

	if (!strcmp(_p.inputType.c_str(), "tdm")) {
		tdmReader.open(_p.filename);
		if (!tdmReader.isOpen())
		{
			cout << ".tdm file is not open." << endl;
			return;
		}

		tdm = true;
	}
	else if (!strcmp(_p.inputType.c_str(), "csv_netcdf_cfar"))
	{
		input.open(_p.filename);
		if (!input.is_open())
		{
			cout << ".netcdf processed file is not open." << endl;
			return;
		}
	}
	else if (!strcmp(_p.inputType.c_str(), "sim"))
	{
		input.open(_p.filename);
		if (!input.is_open())
		{
			cout << ".csv simulation data file is not open." << endl;
			return;
		}
		sim = true;
	}

	// Main loop 
	for (size_t i = 0; i < 2e5; i++)
	{
#ifdef MY_DEBUG
		start = clock();
#endif
		// Read input
		if (tdm) {
			info = tdmReader.readEntryEigen();

			if (!info.size())
				break;
		}
		// This option also reads simulated files, 
		else
		{
			infoTemp = readNetCDFProcessedEntry(input);
			
			if (!infoTemp.size())
				break;

			//cout << infoTemp.segment(0, 6).transpose() << endl;

			// Y, M, D, H, M, S, R(gt/0), Az, Ele, n, z1 , .. zN
			info << infoTemp(6), infoTemp(7), infoTemp(8), 0, 0, infoTemp.segment(0, 6);
		}

		if (!info.size())
			break;
		
		if (info(1) > 360.0 )
			info(1) -= 360.0;
		info(1) = 180.0 - info(1);		// Azimuth correction (?)

		info(1) *= deg2rad;
		info(2) *= deg2rad;

		dateCurr = Astro::date((int)info(5), (int)info(6), (int)info(7),
			(int)info(8), (int)info(9), info(10));

		if (!dtCalculated)
		{
			datePrev = dateCurr;
			dtCalculated = true;
			continue;
		}

		_filter.setT(dateCurr - datePrev);
		datePrev = dateCurr;

		measurements.clear();

		// Check input type
		if (tdm)
			measurements.push_back(info.segment(0, _p.observationDim));
		else
			for (size_t j = 0; j < infoTemp(9); j++) 
			{
				temp << infoTemp(10 + j), info(1), info(2);
				measurements.push_back(temp);
			}

		bearing << info(1), info(2);

		_sensor.setDate(dateCurr);
		_sensor.setZ(measurements);
		_sensor.setBearing(bearing);

#ifdef MY_DEBUG
		bool output_ = false;
			
		end = clock();
		elapsed[0] = elapsedSeconds(start, end);

		if (output_)
			cout << "[" << i << "] " << setprecision(4) << _filter.getQ() << " Meas: " << measurements[0].transpose() << endl;

		// Predict
		start = clock();

		_filter.predict(_mixture, _sensor);
		
		if (output_)
		{
			cout << "Predict:" << endl;
			printMixtureRAZEL(_mixture, _sensor);
		}

		end = clock(); 
		elapsed[1] = elapsedSeconds(start, end);
		
		// Update
		start = clock();

		if (output_)
		{
			//cout << "Update:" << endl;
			//printMixtureRAZEL(_mixture, _sensor);
		}

		_filter.update(_mixture, _sensor);
		end = clock();
		elapsed[2] = elapsedSeconds(start, end);
		
		// Merge
		start = clock();
		_mixture.merge(_p.mergeThreshold);
		end = clock();
		elapsed[3] = elapsedSeconds(start, end);

		// Prune
		start = clock();
		_mixture.prune(_p.pruneThreshold);
		end = clock();
		elapsed[4] = elapsedSeconds(start, end);
		 
		if (output_)
		{
			cout << "Update (m/p):" << endl;
			printMixtureRAZEL(_mixture, _sensor);
		}

		start = clock();
		estimates = _mixture.getEstimates(_p.estimateThreshold);
		end = clock();
		elapsed[5] = elapsedSeconds(start, end);

		elapsedTotal = accumulate(elapsed.begin(), elapsed.end(), 0.0);			// Total execution time

		//cout << setprecision(5) << "Measurements:\t" << elapsed[0] << " sec" << endl;
		//cout << "Predict:\t" << elapsed[1] << " sec" << endl;
		//cout << "Update:\t\t" << elapsed[2] << " sec" << endl;
		//cout << "Merge:\t\t" << elapsed[3] << " sec" << endl;
		//cout << "Prune:\t\t" << elapsed[4] << " sec" << endl;
		//cout << "Estimates:\t" << elapsed[5] << " sec" << endl;
		//cout << "Total: \t\t" << elapsedTotal << " sec" << endl;

		//system("pause");
		
#else
 		_filter.predict(_mixture, _sensor);
		_filter.update(_mixture, _sensor);
		_mixture.merge(_p.mergeThreshold);
		_mixture.prune(_p.pruneThreshold);
		estimates = _mixture.getEstimates(_p.estimateThreshold);
#endif

		printEstimatesToCSVFull(_filter, outputCSV, estimates, _sensor, i, 0);
		//printEstimatesToYAMLFull(_filter, result, estimates, _sensor, i, 0);
		//printMixtureToYAML(_filter, resultSmoothing, _mixture, _sensor, 10);

		// Console output
		double z;
		if (!_sensor.getZ().size())
			z = 0;
		else
			z = _sensor.getZ(0)(0);

		cout << i << " " << z << " [" << estimates.size() << "/" << _mixture.size() << "][" << _filter.getQ() << "]";
		
		if (!estimates.empty()) 
		{
			VectorXd temp = Astro::temeToRAZEL(estimates[0].m, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
			cout << estimates[0].w << " [" << estimates[0].tag[1] << "][" << estimates[0].tag[2] << "] " << temp(0);
		} 
		
		cout << endl;
#ifdef MY_DEBUG
		//cout << endl; 
#endif
	}

	//result << YAML::EndSeq;
	//outputYAML << result.c_str();

	//resultSmoothing << YAML::EndSeq;
	//outputYAMLSmoothing << resultSmoothing.c_str();

	// Close the output files
	outputCSV.close();
	//outputYAML.close();
	//outputYAMLSmoothing.close();
}

/**
* <summary> GMPHD filter test. </summary>
*/
template <typename MultiTargetFilter, typename Mixture>
void runPHDFilter(MultiTargetFilter& _filter, Sensor& _sensor, Mixture& _mixture, parameters& _p)
{
	VectorXd info(11), infoTemp, bearing(2), temp(3);
	auto measurements = _sensor.getZ();
	auto estimates = _mixture.getEstimates(1);
	bool dtCalculated = false;
	Astro::date datePrev, dateCurr;

	double deg2rad = M_PI / 180.0;

	// IO
	TDMReader tdmReader;
	bool tdm = false;
	//YAML::Emitter result;
	ofstream outputCSV("Results/result_gmjott.csv");
	//outputYAML("Results/result_gmjott.yaml");		// Output file 
	ifstream input;

#ifdef MY_DEBUG
	// Debug
	vector<double> elapsed(6);
	clock_t start, end;
	double elapsedTotal;
#endif

	if (!strcmp(_p.inputType.c_str(), "tdm")) {
		tdmReader.open(_p.filename);
		if (!tdmReader.isOpen())
		{
			cout << ".tdm file is not open." << endl;
			return;
		}

		tdm = true;
	}
	else if (!strcmp(_p.inputType.c_str(), "csv_netcdf_cfar"))
	{
		input.open(_p.filename);
		if (!input.is_open())
		{
			cout << ".netcdf processed file is not open." << endl;
			return;
		}

	}

	// Main loop 
	for (size_t i = 0; i < 2e5; i++)
	{
#ifdef MY_DEBUG
		start = clock();
#endif
		// Read input
		if (tdm) {
			info = tdmReader.readEntryEigen();
			if (!info.size())
				break;
		}
		else
		{
			infoTemp = readNetCDFProcessedEntry(input);

			if (!infoTemp.size())
				break;
			info << 0, infoTemp(7), infoTemp(8), 0, 0, infoTemp.segment(1, 6);
			cout << infoTemp.transpose() << endl;
			//cout << info.transpose() << endl;
		}

		if (!info.size())
			break;

		if (info(1) > 360.0)
			info(1) -= 360.0;
		info(1) = 180.0 - info(1);		// Azimuth correction (?)

		info(1) *= deg2rad;
		info(2) *= deg2rad;

		dateCurr = Astro::date((int)info(5), (int)info(6), (int)info(7),
			(int)info(8), (int)info(9), info(10));

		if (!dtCalculated)
		{
			datePrev = dateCurr;
			dtCalculated = true;
			continue;
		}

		_filter.setT(dateCurr - datePrev);
		datePrev = dateCurr;

		measurements.clear();
		// Check input type
		if (tdm)
			measurements.push_back(info.segment(0, _p.observationDim));
		else
			for (size_t j = 0; j < infoTemp(9); j++) {
				temp << infoTemp(10 + j), info(1), info(2);
				measurements.push_back(temp);
			}

		cout << "Measurements" << endl;
		for (auto mm : measurements)
			cout << mm.transpose() << endl;


		bearing << info(1), info(2);

		_sensor.setDate(dateCurr);
		_sensor.setZ(measurements);
		_sensor.setBearing(bearing);

#ifdef MY_DEBUG
		bool output_ = false;


		end = clock();
		elapsed[0] = elapsedSeconds(start, end);

		if (output_)
			cout << "[" << i << "] " << setprecision(4) << " Meas: " << measurements[0].transpose() << endl;

		// Predict
		start = clock();

		_filter.predict(_mixture, _sensor);

		if (output_)
		{
			cout << "Predict:" << endl;
			printMixtureRAZEL(_mixture, _sensor);
		}

		end = clock();
		elapsed[1] = elapsedSeconds(start, end);

		// Update
		start = clock();

		if (output_)
		{
			//cout << "Update:" << endl;
			//printMixtureRAZEL(_mixture, _sensor);
		}

		_filter.update(_mixture, _sensor);
		end = clock();
		elapsed[2] = elapsedSeconds(start, end);

		// Merge
		start = clock();
		_mixture.merge(_p.mergeThreshold);
		end = clock();
		elapsed[3] = elapsedSeconds(start, end);

		// Prune
		start = clock();
		_mixture.prune(_p.pruneThreshold);
		end = clock();
		elapsed[4] = elapsedSeconds(start, end);

		if (output_)
		{
			cout << "Update (m/p):" << endl;
			printMixtureRAZEL(_mixture, _sensor);
		}

		start = clock();
		estimates = _mixture.getEstimates(_p.estimateThreshold);
		end = clock();
		elapsed[5] = elapsedSeconds(start, end);

		elapsedTotal = accumulate(elapsed.begin(), elapsed.end(), 0.0);			// Total execution time

		//cout << setprecision(5) << "Measurements:\t" << elapsed[0] << " sec" << endl;
		//cout << "Predict:\t" << elapsed[1] << " sec" << endl;
		//cout << "Update:\t\t" << elapsed[2] << " sec" << endl;
		//cout << "Merge:\t\t" << elapsed[3] << " sec" << endl;
		//cout << "Prune:\t\t" << elapsed[4] << " sec" << endl;
		//cout << "Estimates:\t" << elapsed[5] << " sec" << endl;
		//cout << "Total: \t\t" << elapsedTotal << " sec" << endl;

		//system("pause");

#else
		_filter.predict(_mixture, _sensor);
		_filter.update(_mixture, _sensor);
		_mixture.merge(_p.mergeThreshold);
		_mixture.prune(_p.pruneThreshold);
		estimates = _mixture.getEstimates(_p.estimateThreshold);
#endif

		printEstimatesToCSVFull(_filter, outputCSV, estimates, _sensor, i, 0);
		//printEstimatesToYAMLFull(_filter, result, estimates, _sensor, i, 0);

		// Console output
		double z;
		if (!_sensor.getZ().size())
			z = 0;
		else
			z = _sensor.getZ(0)(0);

		cout << i << " " << z << " [" << estimates.size() << "/" << _mixture.size() << "]";

		if (!estimates.empty())
		{
			VectorXd temp = Astro::temeToRAZEL(estimates[0].m, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
			cout << estimates[0].w << " [" << estimates[0].tag[1] << "][" << estimates[0].tag[2] << "] " << temp(0);
		}

		cout << endl;
	}

	//outputYAML << result.c_str();

	outputCSV.close();
	//outputYAML.close();
}

/**
 * SMC JoTT Testing
 */
template<typename MultiTargetFilter, typename Swarm>
void runSMCJoTTFilter(MultiTargetFilter & _filter, Sensor & _sensor, Swarm & _mixture, parameters & _p)
{
	VectorXd info(11), infoTemp, bearing(4), bearingPrev(4), temp(3), avgEstimate, avgEstimateRAZEL, gt(3), gtRate(3), gtPrev(3);
	auto measurements = _sensor.getZ();
	bool dtCalculated = false;
	Astro::date datePrev, dateCurr;

	double deg2rad = M_PI / 180.0;

	// IO
	bool tdm = false, sim = false;
	TDMReader tdmReader;
	ifstream input;
	ofstream outputCSV("Results/result_smcjott.csv");

#ifdef MY_DEBUG
	// Debug
	vector<double> elapsed(6);
	//clock_t start, end;
	//double elapsedTotal;
#endif

	if (!strcmp(_p.inputType.c_str(), "tdm")) {
		tdmReader.open(_p.filename);
		if (!tdmReader.isOpen())
		{
			cout << ".tdm file is not open." << endl;
			return;
		}

		tdm = true;
	}
	else if (!strcmp(_p.inputType.c_str(), "csv_netcdf_cfar"))
	{
		input.open(_p.filename);
		if (!input.is_open())
		{
			cout << ".netcdf processed file is not open." << endl;
			return;
		}
	}
	else if (!strcmp(_p.inputType.c_str(), "sim"))
	{
		input.open(_p.filename);
		if (!input.is_open())
		{
			cout << ".csv simulation data file is not open." << endl;
			return;
		}
		sim = true;
	}

	for (size_t i = 0; i < 2e5; i++)
	{
		// Read input
		if (tdm) {
			info = tdmReader.readEntryEigen();
			if (!info.size())
				break;
		}
		// This option also reads simulated files, 
		else
		{
			infoTemp = readNetCDFProcessedEntry(input);
			if (!infoTemp.size())
				break;

			// Y, M, D, H, M, S, R(gt/0), Az, Ele, n, z1 , .. zN
			info << infoTemp(6), infoTemp(7), infoTemp(8), 0, 0, infoTemp.segment(0, 6);
			gt << infoTemp(6), infoTemp(7), infoTemp(8);
			
		}
		
		if (!info.size())
			break;

		if (info(1) > 360.0)
			info(1) -= 360.0;
		info(1) = 180.0 - info(1);		// Azimuth correction (?)

		info(1) *= deg2rad;
		info(2) *= deg2rad;

		dateCurr = Astro::date((int)info(5), (int)info(6), (int)info(7),
			(int)info(8), (int)info(9), info(10));

		if (!dtCalculated)
		{
			datePrev = dateCurr;
			dtCalculated = true;

			bearing << info(1), info(2), 0, 0;
			bearingPrev = bearing;
			_sensor.setBearing(bearing);
			_sensor.setGT(gt);

			continue;
		}

		_filter.setT(dateCurr - datePrev);
		datePrev = dateCurr;

		measurements.clear();

		// Check input type
		if (tdm)
			measurements.push_back(info.segment(0, _p.observationDim));
		else
			for (size_t j = 0; j < infoTemp(9); j++)
			{
				temp << infoTemp(10 + j), info(1), info(2);
				measurements.push_back(temp);
			}
		
		// Calculate bearing rate (temporary)
		bearingPrev = _sensor.getBearing();
		gtPrev = _sensor.getGT();
		bearing << info(1), info(2), (info(1) - bearingPrev(0)) / _filter.getT(), (info(2) - bearingPrev(1)) / _filter.getT();
		gtRate << (info(0) - gtPrev(0)) / _filter.getT(), bearing(2), bearing(3);

		_sensor.setGT(gt);
		_sensor.setGTRate(gtRate);
		_sensor.setDate(dateCurr);
		_sensor.setZPrev(_sensor.getZ());
		_sensor.setZ(measurements);
		_sensor.setBearing(bearing);

		// TEMPORARY
		if (i == 1)
		{
			//VectorXd lBound(3), uBound(3);
			//double t = 0.1;
			//lBound << 500, info(1) - t, info(2) - t;
			//uBound << 1500, info(1) + t, info(2) + t;

			VectorXd bounds(6);
			bounds << 50, 0.1, 0.1, 0.1, 0.00001, 0.00001;

			_mixture.populateRandomRAZEL2(_sensor, bounds, true);

			cout << "Mixture size: " << _mixture.size() << endl;
		}

		_filter.predict(_mixture, _sensor);
		_filter.update(_mixture, _sensor);

		// Console output
		double z;
		if (!_sensor.getZ().size())
			z = 0;
		else
			z = _sensor.getZ(0)(0);

		avgEstimate = _mixture.getWeightedAverage();
		avgEstimateRAZEL = Astro::temeToRAZEL(avgEstimate, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
		cout << i << " " << _filter.getQ() << " " << info(0) << " [" << _mixture.size() << "]" << avgEstimateRAZEL.transpose().head(3);
		cout << endl;

		outputCSV << i << ", " << _filter.getQ() << ", " << z << ", ";

		for (size_t j = 0; j < 3; j++)
		{
			outputCSV << avgEstimateRAZEL(j) << ", ";
		}

		outputCSV << 0 << std::endl;

	}

	// Close the output files
	outputCSV.close();

}