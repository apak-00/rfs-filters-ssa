#pragma once
#include <vector>
#include <fstream>
#include <Eigen/Core>
#include <yaml-cpp/yaml.h>

#include "Sensor.h"

using namespace std;
using namespace Eigen;

/*
 * <summary> Input/output helpers namespace </summary>
 */
namespace IOHelpers
{
	
	/**
	 * <summary> Parameters structure for storing testing configuration for the Gaussian Mixture 
	 * Random Finite Sets Filters. </summary>
	 */
	struct parameters {

		std::string filename;				// Input filename
		std::string inputType;				// Input type (e.g. .tdm, .csv, .ncdf)

		// Sensor-related
		std::string sensorType;				// Sensor type (only radar atm)
		std::string observationType;		// Observation type (only range, azimuth, elevation atm)
		VectorXd sensorPosition;			// Sensor position (WGS-84 latitude, longitude, altitude (m))
		size_t observationDim;				// Dimensionality of the observation space
		double pD;							// Probability of target detection for constant-pD filters
		double kappa;						// Clutter intensity
		double lambda;						// Expected number of false alarms
		double clutter;						// Clutter parameter (to be specified)
		MatrixXd R;							// Observation noise matrix
		MatrixXd H;							// Observation matrix (for standard Kalman Filter)

		// Single-target filter related
		std::string singleTargetFilterType;	// Type of the single-target filter (e.g. KF, EKF, UKF)
		std::string motionModel;			// Motion model (only constant velocity model is used atm)
		size_t stateDim;					// Dimensionality of the state space
		MatrixXd Q;							// Process noise matrix
		MatrixXd F;							// State transition matrix (for KF and EKF)
		double dt;							// Timestep

		// Multiple-target filter related
		std::string multipleTargetFilterType;	// Type of the multi-target filter (e.g. GMPHD, GMJoTT, etc.)
		size_t gaussianMixtureMaxSize;			// Maximum size of th Gaussian Mixture
		std::string birthType;					// Birth type (uniform birth is used atm, can be changed to random)
		size_t birthSize;						// Number of birth components per timestep
		MatrixXd birthCovariance;				// Initial covariance matrix of the birth components
		double mergeThreshold;					// Merge threshold of the Gaussian Mixture
		double pruneThreshold;					// Prune threshold of the Gaussian Mixture
		double estimateThreshold;				// Estimation threshold of the Gaussian Mixture
		double birthIntensity;					// Overall intensity of the birth components
		double pB;								// Probability of spontaneous target birth
		double pS;								// Probability of target survival
		double qInit;							// Initial pribability of target existence (for JoTT filter)

		double epsilon;							// Epsilon increment for the Beta Gaussian Components

		// Temporary parameters
		double clutterMultiplierTemp;
		double ukfSigmaSamplerW;

		// Other
		double pmx;
		double pmy;
		double lod;
	};

	/*
	 * <summary> Enum used in the TDMReader class for signle entry field identification. </summary>
	 */
	enum TDM_PARAMS { RANGE, AZI, ELE, SS, CPSS };

	/*
	 * <summary> TDMReader class. Reads and processes tracking data message (.tdm) files. </summary>
	 */
	class TDMReader
	{
	public:

		TDMReader();
		TDMReader(string _filename);
		~TDMReader();

		void open(const string& _filename);
		string getHeader();
		vector<double> readEntry();
		Eigen::VectorXd readEntryEigen();

		bool isOpen();

		static void TDMDirToCSV(const string & _s);
		static void TDMToYAML(const string & _filename);

	private:
		ifstream inputFile;
		string buffer, header;
		bool headerRead;

		void readHeader();

		// .tdm file delimeters
		static constexpr const char* DATA_DELIM_START = "DATA_START";
		static constexpr const char* DATA_DELIM_END = "DATA_STOP";
		static constexpr const char* META_DELIM_START = "META_START";
		static constexpr const char* META_DELIM_END = "META_STOP";

		const unsigned short DATA_DELIM_START_LEN = 10;
		const unsigned short DATA_DELIM_END_LEN = 9;
		const unsigned short META_DELIM_START_LEN = 10;
		const unsigned short META_DELIM_END_LEN = 9;
	};

	bool has_suffix(const string& _s, const string& _suffix);
	vector<string>& split(const string &s, char delim, vector<string> &elements);
	double elapsedSeconds(const double& _start, const double& _end);

	ostream & operator<<(std::ostream & _os, const parameters & _params);
	parameters readParametersYAML(const string & _filename);
	VectorXd readNetCDFProcessedEntry(ifstream& _input);

	void printVector(ofstream& _os, const VectorXd & _v);

	/*
	 * <summary> Estimate output to .csv file. </summary>
	 * <param name = "_filter"> A reference to the GMRFSFilter. </param>
	 * <param name = "_os"> A reference to the output stream. </param>
	 * <param name = "_estimates"> A reference to the vector of components (templated). </param>
	 * <param name = "_sensor"> A reference to the sensor. </param>
	 * <param name = "_id"> Iteration count </param>
	 * <param name = "_elapsed"> Iteration execution time. Defaults to zero. </param>
 	 */
	template<typename T, typename F>
	void printEstimatesToCSVFull(F& _filter, ofstream& _os, const vector<T>& _estimates, Sensor& _sensor, const size_t& _id, const double& _elapsed = 0)
	{
		using fType = F;
		auto d = _sensor.getDate();
		VectorXd z;

		if (_sensor.getZ().size())
		{
			z = _sensor.getZ(0);
		}
		else
		{
			z = VectorXd::Zero(_sensor.getZDim());
		}

		size_t sDim = _sensor.getSDim();

		VectorXd sensorPosTEME = Astro::ecefToTEME(Astro::geodeticToECEF(_sensor.getPosition()), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
		VectorXd dummy = VectorXd::Zero(sDim);

		if (!_estimates.size())
		{
			// Date and first measurement
			_os << d << "," << z(0) << "," << z(1) << "," << z(2) << ",";

			// Sensor position TEME
			printVector(_os, sensorPosTEME);
			_os << ",";

			// Estimates' size and q
			_os << _estimates.size() << "," << _filter.getQ() << ",";

			// Estimate RAZEL
			printVector(_os, dummy);
			_os << ",";

			// Estimate TEME
			printVector(_os, T::getEmptyInfo());
		}
		else
		{
			for (auto i : _estimates)
			{
				// Date and first measurement
				_os << d << "," << z(0) << "," << z(1) << "," << z(2) << ",";

				// Sensor position TEME
				printVector(_os, sensorPosTEME);
				_os << ",";

				// Estimates' size and q
				_os << _estimates.size() << "," << _filter.getQ() << ",";

				// Estimate RAZEL
				printVector(_os, Astro::temeToRAZEL(i.m, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp()));
				_os << ",";

				// Estimate TEME
				//printVector(_os, i.m);
				_os << i;
			}
		}

		if (false)
		{
			_os << ",";
			VectorXd tempZ(6);
			tempZ << _sensor.getZ(0), 0, 0, 0;
			VectorXd tempZTEME = Astro::razelToTEME(tempZ, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
			printVector(_os, tempZTEME);
			//cout << tempZTEME.transpose() << endl;
		}

		_os << endl;
		 
	}

	/*
	* <summary> Estimate output to .yaml file. </summary>
	* <param name = "_filter"> A reference to the GMRFSFilter. </param>
	* <param name = "_yamle"> A reference to the associated YAML emitter. </param>
	* <param name = "_estimates"> A reference to the vector of components (templated). </param>
	* <param name = "_sensor"> A reference to the sensor. </param>
	* <param name = "_id"> Iteration count </param>
	* <param name = "_elapsed"> Iteration execution time. Defaults to zero. </param>
	*/
	template<typename T, typename F>
	inline void printEstimatesToYAMLFull(F& _filter, YAML::Emitter& _yamle, const vector<T>& _estimates, Sensor& _sensor, const size_t& _id, const double& _elapsed = 0) 
	{
		auto d = _sensor.getDate();
		VectorXd z = _sensor.getZ(0);
		size_t sDim = _sensor.getSDim();
		VectorXd sensorPosTEME = Astro::ecefToTEME(Astro::geodeticToECEF(_sensor.getPosition()), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
		
		_yamle << YAML::Flow;

		for (auto i : _estimates)
		{
			// Date and first measurement
			_yamle << YAML::BeginSeq << d.year << d.month << d.day << d.hour << d.minute << d.sec
				<< z(0) << z(1) << z(2);

			// Sensor position TEME
			for (size_t j = 0; j < sDim; j++)
				_yamle << sensorPosTEME(j);

			_yamle << _estimates.size() << _filter.getQ();

			// Estimate RAZEL
			VectorXd eRAZEL = Astro::temeToRAZEL(i.m, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
			for (size_t j = 0; j < sDim; j++)
				_yamle << eRAZEL(j);

			// Estimate TEME
			for (size_t j = 0; j < sDim; j++)
				_yamle << i.m(j);

			// Miscellaneous
			_yamle << i.w << i.tag[1] << i.P.determinant() << i.P.block<3, 3>(0, 0).determinant() << i.P.block<3, 3>(3, 3).determinant();

			_yamle << YAML::EndSeq;
		}
	}


	template<typename M, typename F>
	inline void printMixtureToYAML(F& _filter, YAML::Emitter& _yamle, const M& _mixture, Sensor& _sensor, const size_t& _printSize)
	{
		auto d = _sensor.getDate();
		VectorXd z = _sensor.getZ(0);
		size_t sDim = _sensor.getSDim();
		VectorXd sensorPosTEME = Astro::ecefToTEME(Astro::geodeticToECEF(_sensor.getPosition()), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
		
		size_t n = _mixture.size();
		if (n > _printSize)
			n = _printSize;

		for (size_t i = 0; i < n; i++)
		{
			auto gc = _mixture[i];

			_yamle << YAML::Flow << YAML::BeginSeq << d.year << d.month << d.day << d.hour << d.minute << d.sec
				<< z(0) << z(1) << z(2);
			
			// Sensor position TEME
			for (size_t j = 0; j < sDim; j++)
				_yamle << sensorPosTEME(j);

			// Estimate RAZEL
			VectorXd eRAZEL = Astro::temeToRAZEL(gc.m, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
			for (size_t j = 0; j < sDim; j++)
				_yamle << eRAZEL(j);

			// Estimate TEME
			for (size_t j = 0; j < sDim; j++)
				_yamle << gc.m(j);

			// Miscellaneous
			_yamle << gc.w << gc.tag[1] << gc.P.determinant() << gc.P.block<3, 3>(0, 0).determinant() << gc.P.block<3, 3>(3, 3).determinant();

			_yamle << YAML::EndSeq;
		}
	}

	template<typename Mixture>
	inline void printMixtureRAZEL(Mixture& _m, Sensor& _sensor)
	{
		VectorXd temp;
		for (size_t i = 0; i < _m.size(); i++)
		{
			temp = Astro::temeToRAZEL(_m[i].m, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
			cout << "\t" << i << " " << setprecision(4)  <<  _m[i].w << " " << temp.head(3).transpose() << endl;
		}
	}

	VectorXd readSimulatedEntry(ifstream& _input);

	
}