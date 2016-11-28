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
		const size_t& _nBirthComponents, const double & _birthWeight, const double & _pS, const double& _pB, const double & _q, 
		VectorXd _noise) 
		: nBirthComponents(_nBirthComponents), pS(_pS),  q(_q), pB(_pB), noise(_noise), b(_nBirthComponents, 6)
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

	void setLowerBirthBoundRange(const decltype(lowerBirthBoundRange)& _lowerBirthBoundRange) 
	{
		lowerBirthBoundRange = _lowerBirthBoundRange;
	}

	void setUpperBirthBoundRange(const decltype(upperBirthBoundRange)& _lowerBirthBoundRange)
	{
		upperBirthBoundRange = _lowerBirthBoundRange;
	}

	void setBirthSigmaRange(const decltype(birthSigmaRange)& _birthSigmaRange)
	{
		birthSigmaRange = _birthSigmaRange;
	}

	/*
	 * <summary> Prediction </summary>
	 */
	void predict(particle_mixture& _ps, Sensor& _sensor)
	{		
		double t, wb;

		// Probability of target existence prediction (28)
		qPred = pB * (1 - q) + pS * q;	

		t = pS * q / qPred;

		// Mean prediction
		for (auto &p : _ps.components)
		{
			propagate(p.m, dt, noise);
			p.w *= t;
		}
						
		// Birth
		wb = birthWeight * pB * (1 - q) / qPred;
		auto zPrev = _sensor.getZPrev();
		VectorXd tempBirth(6);

		// If there are previous measurements, us their location to place birth components
		if (!zPrev.size()) 
		{
			// Number of birth components for each measurement
			size_t birthSizeZ = nBirthComponents * zPrev.size();

			// Number of mesurements @ previous timestep times number of particles to be born
			b = particle_mixture(birthSizeZ, _sensor.getSDim(), wb);

			for (size_t i = 0; i < birthSizeZ; i++) 
			{
				// The particles are born within a normal distribution centered at the position of the previous measurement
				// At the moment, only range uncertainty is taken into account
				// TODO: Change the constant standard deviation
				std::normal_distribution<double> distribution(0, 10);

				for (size_t j = 0; j < nBirthComponents; j++)
				{
					tempBirth << distribution(generator) + zPrev[i](0), zPrev[i](1), zPrev[i](2),  0, 0, 0;
					b.components[i * nBirthComponents + j].m = Astro::razelToTEME(tempBirth, _sensor.getPosition(), _sensor.getDateJD(),
						_sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
				}
			}

			// Concatenate predicted and birth components
			_ps = _ps + b;
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
				tempBirth << distribution(generator) + zPrev[i](0), _sensor.getBearing(), 0, 0, 0;
				b.components[i].m = Astro::razelToTEME(tempBirth, _sensor.getPosition(), _sensor.getDateJD(),
					_sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
			}

		}
	}

	/*
	 * <summary> Update </summary>
	 */
	void update(particle_mixture& _ps, Sensor& _sensor)
	{
		double nThreshold = 200;
		double cz = 1.0 / 57903 * 200;
		auto pD = _sensor.getPD();
		double I1 = _ps.weightSum() * pD, I2, wSum = 0, delta_k = 0, g = 0;
		std::vector<double> ll(_ps.size(), 0);		// Likelihoods
		size_t n = _ps.size();

		VectorXd z = VectorXd::Zero(3);
		std::vector<VectorXd> zPred(_ps.size());

		// Pre-calculate predicted measurements
		for (size_t i = 0; i < _ps.size(); i++)
			zPred[i] = observe(_ps.components[i].m, _sensor);

		// For each measurement
		for (size_t i = 0; i < _sensor.getZ().size(); i++) 
		{
			z = _sensor.getZ(i);
			n = _ps.size();
			I2 = 0;

			// For each particle
			for (size_t j = 0; j < _ps.components.size(); j++)
			{
				g = MathHelpers::gaussianLikelihood(z, zPred[j], _sensor.getS());	// Compute likelihood
				ll[j] += g;		// Likelihood sum
				I2 += pD * g * _ps.components[j].w;		// Approximate I2 integral (85)
			}

			delta_k += I2;
		}

		delta_k = I1 - delta_k / cz;		// delta_k approximation (86)
		q = (1 - delta_k) / (1 - q * delta_k) * q;

		// Update the weights (87)
		for (size_t i = 0; i < _ps.components.size(); i++)
			_ps.components[i].w = (1 - pD + pD * ll[i] / cz) * _ps.components[i].w;
		
		// Normalize the weights
		_ps.normalizeWeights();
		
		// Sampling importance re-sampling
		// If the effective number of particles is less than a threshold, resample:
		if (_ps.getEffectiveN() < nThreshold)
			_ps.resampleITS(1000);
	}
};

