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

#include "ExtendedKalmanFilter.h"
#include "UnscentedKalmanFilter.h"

#include "IOHelpers.h"

using namespace std;
using namespace Eigen;
using namespace IOHelpers;

// Function for JoTT testing
void testFilter(parameters& _p);
template <typename MultiTargetFilter, typename Mixture>
void runFilter(MultiTargetFilter& _filter, Sensor& _sensor, Mixture& _mixture, parameters& _p);

// Test
int main(int arcg, char** argv) 
{	
	string filename_params = "config.yaml";
	parameters p = readParametersYAML(filename_params);
	testFilter(p);

	return 0; 
}

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
	else if (!strcmp(_p.singleTargetFilterType.c_str(), "ekf"))
		kf = make_shared<UnscentedKalmanFilter>(_p.Q * 0.1);

	// Sensor 
	Sensor sensor = Sensor(_p.observationDim, _p.stateDim, _p.pD, _p.lambda, 1e-6, _p.R, _p.H);
	sensor.setPosition(_p.sensorPosition);
	// Additional reference system-related parameters
	sensor.setXp(0);
	sensor.setYp(0);
	sensor.setLOD(0);

#ifdef MY_DEBUG
	clock_t end = clock();
	double elapsed = double(end - begin) / CLOCKS_PER_SEC;
	cout << "Load time:\t" << elapsed << " sec" << endl;
#endif

	// Choose multi-target filter
	if (!strcmp(_p.multipleTargetFilterType.c_str(), "gmjott")) {
		GMJoTTFilter gmjott(kf, _p.birthSize, _p.birthIntensity, _p.pS, _p.birthCovariance,
			VectorXd::Zero(_p.stateDim), VectorXd::Zero(_p.stateDim), _p.qInit, _p.pB);
		gaussian_mixture gm(_p.stateDim, _p.gaussianMixtureMaxSize);
		runFilter(gmjott, sensor, gm, _p);
	}
	else if (!strcmp(_p.multipleTargetFilterType.c_str(), "bgmjott")) {
		BGMJoTTFilter bgmjott(kf, _p.birthSize, _p.birthIntensity, _p.pS, _p.birthCovariance,
			VectorXd::Zero(_p.stateDim), VectorXd::Zero(_p.stateDim), _p.qInit, _p.pB, _p.epsilon);
		beta_gaussian_mixture bgm(_p.stateDim, _p.gaussianMixtureMaxSize);
		runFilter(bgmjott, sensor, bgm, _p);
	}
}
 
/**
* <summary> GMJoTT filter test. </summary>
*/
template <typename MultiTargetFilter, typename Mixture>
void runFilter(MultiTargetFilter& _filter, Sensor& _sensor, Mixture& _mixture, parameters& _p)
{
	VectorXd info(11), infoTemp, bearing(2), temp(3);
	auto measurements = _sensor.getZ();
	auto estimates = _mixture.getEstimates(1);
	bool dtCalculated = false;
	Astro::date datePrev, dateCurr;

	// IO
	TDMReader tdmReader;
	bool tdm;
	YAML::Emitter result;
	ofstream outputCSV("Results/result_gmjott.csv") , outputYAML("Results/result_gmjott.yaml");		// Output file 
	ifstream input;

#ifdef MY_DEBUG
	// Debug
	vector<double> elapsed(6);
	clock_t start, end;
	double elapsedTotal;
#endif

	if (!strcmp(_p.inputType.c_str(), "tdm")) {
		tdmReader.open(_p.filename);
		tdm = true;
	}
	else if (!strcmp(_p.inputType.c_str(), "csv_netcdf_cfar"))
		input.open(_p.filename);

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
		}

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

		bearing << info(1), info(2);

		_sensor.setDate(dateCurr);
		_sensor.setZ(measurements);
		_sensor.setBearing(bearing);

#ifdef MY_DEBUG
		end = clock();
		elapsed[0] = elapsedSeconds(start, end); 


		// Predict
		start = clock();
		_filter.predict(_mixture, _sensor);
		end = clock();
		elapsed[1] = elapsedSeconds(start, end);

		// Update
		start = clock();
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

		start = clock();
		estimates = _mixture.getEstimates(_p.estimateThreshold);
		end = clock();
		elapsed[5] = elapsedSeconds(start, end);

		elapsedTotal = accumulate(elapsed.begin(), elapsed.end(), 0.0);			// Total execution time


		cout << setprecision(5) << "Measurements:\t" << elapsed[0] << " sec" << endl;
		cout << "Predict:\t" << elapsed[1] << " sec" << endl;
		cout << "Update:\t\t" << elapsed[2] << " sec" << endl;
		cout << "Merge:\t\t" << elapsed[3] << " sec" << endl;
		cout << "Prune:\t\t" << elapsed[4] << " sec" << endl;
		cout << "Estimates:\t" << elapsed[5] << " sec" << endl;
		cout << "Total: \t\t" << elapsedTotal << " sec" << endl;

		system("pause");

#else
		_filter.predict(_mixture, _sensor);
		_filter.update(_mixture, _sensor);
		_mixture.merge(_p.mergeThreshold);
		_mixture.prune(_p.pruneThreshold);
		estimates = _mixture.getEstimates(_p.estimateThreshold);
#endif

		printEstimatesToCSVFull(_filter, outputCSV, estimates, _sensor, i, 0);
		printEstimatesToYAMLFull(_filter, result, estimates, _sensor, i, 0);

		// Console output
		cout << i << " [" << estimates.size() << "/" << _mixture.size() << "][" << _filter.getQ() << "]\t"
			<< dateCurr.hour << ":" << dateCurr.minute << ":" << dateCurr.sec << " ";

		if (!estimates.empty()) 
			cout << estimates[0].w << " [" << estimates[0].tag[1] << "][" << estimates[0].tag[2] << "] ";
		
		cout << endl;
	}

	outputYAML << result.c_str();

	outputCSV.close();
	outputYAML.close();
}