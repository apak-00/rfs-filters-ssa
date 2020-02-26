#pragma once
#include <Eigen/Dense>
#include "Sensor.h"
#include "MixtureModels.h"
#include "MathHelpers.h"

#include <functional>

class SMCJoTTFilter
{
protected:

	std::random_device rDev;
	std::mt19937 generator;

	double q, qPred;

	size_t nBirthComponents;
	double birthWeight;
	double pS, pB;
	double dt;

	// Temporary
	double lowerBirthBoundRange;
	double upperBirthBoundRange;
	double birthSigmaRange;

	// Propagation noise (temporary)
	VectorXd noise;

	// Birth
	particle_mixture b;

	std::function<VectorXd(const VectorXd&, const double&, const VectorXd&)> propagate;
	std::function<VectorXd(const VectorXd&, const Sensor&)> observe;

public:
	SMCJoTTFilter(std::function<VectorXd(const VectorXd&, const double&, const VectorXd&)> _propagate, 
		std::function<VectorXd(const VectorXd&, const Sensor&)> _observe, 
		const size_t& _nBirthComponents, const double & _birthWeight, const double & _pm, const double& _pB, const double & _q, 
		VectorXd _noise) 
		: nBirthComponents(_nBirthComponents), pS(_pm),  q(_q), pB(_pB), noise(_noise), b(_nBirthComponents, 6), birthWeight(_birthWeight)
	{
		propagate = _propagate;
		observe = _observe;

		generator = std::mt19937(rDev());
	}

	auto getQ() { return q; }
	auto getQPred() { return qPred; }
	auto getT() { return dt; }
	void setT(const double& _dt) { dt = _dt; }

	// Temporary
	// TODO: For the future: custom birth patterns
	const auto getLowerBirthBoundRange() { return lowerBirthBoundRange; }
	const auto getUpperBirthBoundRange() { return upperBirthBoundRange; }
	const auto getBirthSigma() { return birthSigmaRange; }

	void setLowerBirthBoundRange(const decltype(lowerBirthBoundRange)& _lowerBirthBoundRange) { lowerBirthBoundRange = _lowerBirthBoundRange; }
	void setUpperBirthBoundRange(const decltype(upperBirthBoundRange)& _lowerBirthBoundRange) { upperBirthBoundRange = _lowerBirthBoundRange; }
	void setBirthSigmaRange(const decltype(birthSigmaRange)& _birthSigmaRange) { birthSigmaRange = _birthSigmaRange; }

	/*
	 * <summary> Prediction for the Sequential Monte-Carlo JoTT filter. </summary>
	 * <param name = "_pm"> A particle mixture to predict. </param>
	 * <param name = "_sensor"> A sensor to take the bering from in order to perform the birth. </param>
	 */
	void predict(particle_mixture& _pm, Sensor& _sensor)
	{		
		// Debug
		//std::cout << "w = " << _pm.components[0].w << " | "<< _pm.components[0].m.transpose() << std::endl;

		double t, wb;

		// Probability of target existence prediction (28)
		qPred = pB * (1 - q) + pS * q;	

		t = pS * q / qPred;

		// Mean prediction
		for (auto &p : _pm.components)
		{
			propagate(p.m, dt, noise);
			p.w *= t;
		}
						
		// Birth
		wb = birthWeight * pB * (1 - q) / qPred;
		auto zPrev = _sensor.getZPrev();
		VectorXd tempBirth(6);

		// If there are previous measurements, use their location to place birth components
		if (zPrev.size()) 
		{
			// Number of birth components for each measurement
			size_t birthSizeZ = zPrev.size();

			// Number of mesurements @ previous timestep times number of particles to be born
			b = particle_mixture(nBirthComponents * birthSizeZ, _sensor.getSDim(), wb);

			VectorXd bounds(6);
			bounds << 100, 0.1, 0.1, 0.1, 0.00001, 0.00001;
			VectorXd posRAZEL(3), velRAZEL(3);
			posRAZEL = _sensor.getGT();
			velRAZEL = _sensor.getGTRate();

			for (size_t i = 0; i < birthSizeZ; i++) 
			{
				// The particles are born within a normal distribution centered at the position of the previous measurement
				// At the moment, only range uncertainty is taken into account
				// TODO: Change the constant standard deviation
				//std::normal_distribution<double> distribution(0, 10);
				//std::normal_distribution<double> distributionV(0, 10);
				

				std::uniform_real_distribution<double> dr(posRAZEL(0) - bounds(0), posRAZEL(0) + bounds(0));		// Distribution for Range
				std::uniform_real_distribution<double> da(posRAZEL(1) - bounds(1), posRAZEL(1) + bounds(1));		// Distribution for Azimuth
				std::uniform_real_distribution<double> de(posRAZEL(2) - bounds(2), posRAZEL(2) + bounds(2));		// Distribution for Elevation

				std::uniform_real_distribution<double> ddr(velRAZEL(0) - bounds(3), velRAZEL(0) + bounds(3));		// Distribution for Range rate
				std::uniform_real_distribution<double> dda(velRAZEL(1) - bounds(4), velRAZEL(1) + bounds(4));		// Distribution for Azimuth rate
				std::uniform_real_distribution<double> dde(velRAZEL(2) - bounds(5), velRAZEL(2) + bounds(5));		// Distribution for Elevation rate

				for (size_t j = 0; j < nBirthComponents; j++)
				{
					tempBirth <<  dr(generator), da(generator), de(generator), ddr(generator), dda(generator), dde(generator);
					b.components[i * nBirthComponents + j].m = Astro::razelToTEME(tempBirth, _sensor.getPosition(), _sensor.getDateJD(),
						_sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
				}
			}
		}
		// If no measurements are received in the previous timestep, perform random birth
		else
		{
			// Just a fixed number of particles in the 750 - 1250km range
			// TODO: Change the fixed numbers. Again.
			b = particle_mixture(nBirthComponents, _sensor.getSDim(), wb);
			std::uniform_real_distribution<double> distribution(lowerBirthBoundRange, upperBirthBoundRange);

			for (size_t i = 0; i < nBirthComponents; i++)
			{
				auto bearing = _sensor.getBearing();
				tempBirth << distribution(generator), bearing(0), bearing(1), 0, 0, 0;
				b.components[i].m = Astro::razelToTEME(tempBirth, _sensor.getPosition(), _sensor.getDateJD(),
					_sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
			}
		}

		// Concatenate predicted and birth components
		//std::cout << _pm.size() << " + ";
		//std::cout << b.size() << " = ";
		_pm = _pm + b;
		//std::cout << _pm.size() << std::endl;

	}

	/*
	 * <summary> Update for the Sequential Monte-Carlo JoTT filter. </summary>
	 * <param name = "_pm"> A particle mixture to update. </param>
	 * <param name = "_sensor"> A sensor to take the measurement and other information from. </param>
	 */
	void update(particle_mixture& _pm, Sensor& _sensor)
	{
		//// Debug before
		//std::cout << "Pre-update weigths (" << _pm.size() << "): " << std::endl;
		//for (size_t i = 0; i < 10; i++) {
		//	std::cout << "U- w = " << _pm.components[i].w << " | " << _pm.components[i].m.transpose() << std::endl;
		//}

		_pm.normalizeWeights();

		double cz = 1.0 / 57903 * 200;
		auto pD = _sensor.getPD();
		double I2 = 0, wSum = 0, delta_k = 0, g = 0;
		std::vector<double> ll(_pm.size(), 0);		// Likelihoods
		size_t n = _pm.size();

		// Current measurement and predicted measurement
		VectorXd z = VectorXd::Zero(3);
		std::vector<VectorXd> zPred(_pm.size());

		// Pre-calculate predicted measurements
		for (size_t i = 0; i < _pm.size(); i++)
			zPred[i] = observe(_pm.components[i].m, _sensor);

		// Loop for each measurement ...
		for (size_t i = 0; i < _sensor.getZ().size(); i++) 
		{
			z = _sensor.getZ(i);
			n = _pm.size();

			// Loop for each particle ...
			for (size_t j = 0; j < _pm.components.size(); j++)
			{
				//g = MathHelpers::gaussianLikelihood(z, zPred[j], _sensor.getS());	// Compute likelihood
				g = MathHelpers::gaussianLikelihood(z, zPred[j].head(3), _sensor.getR());
				ll[j] += g;		// Likelihood sum
				I2 += g * _pm.components[j].w;		// Approximate I2 integral (85)
			}
		}

		delta_k = pD * (1 - I2 / cz);		// delta_k approximation (86)
		q = (1 - delta_k) / (1 - q * delta_k) * q;

		// Update the weights (87)
		for (size_t i = 0; i < _pm.components.size(); i++)
			_pm.components[i].w = (1 - pD + pD * ll[i] / cz) * _pm.components[i].w;

		std::sort(_pm.components.begin(), _pm.components.end(), [](particle a, particle b) {return b.w < a.w; });

		// Truncate
		if (_pm.size() >= _pm.nMax) {
			std::vector<particle> truncated(_pm.components.begin(), _pm.components.begin() + _pm.nMax);
			_pm.components = truncated;
		}
		
		// Normalize the weights
		_pm.normalizeWeights();
		
		// Sampling importance re-sampling
		// If the effective number of particles is less than a threshold, resample:
		auto effectiveNumberOfParticles = _pm.getEffectiveN();
		if (effectiveNumberOfParticles < 50) //nThresold
			if (_pm.size() <= _pm.nMax)
				_pm.resampleITS(_pm.size());
			else 
				_pm.resampleITS(_pm.nMax);

		//std::cout << "Post-update weigths (" << _pm.size() << "): " << std::endl;
		//// Debug after
		//for (size_t i = 0; i < 10; i++) {
		//	std::cout << "U- w = " << _pm.components[i].w << " | " << _pm.components[i].m.transpose() << std::endl;
		//}
	}
};

