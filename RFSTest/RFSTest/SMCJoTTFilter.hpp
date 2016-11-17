#pragma once
#include <Eigen/Dense>
#include "Sensor.h"
#include "gmm.h"
#include "MathHelpers.h"

#include <functional>

class SMCJoTTFilter
{
protected:
	double q, qPred;

	size_t nBirthComponents;
	double birthWeight;
	double pS, pB;
	double dt;

	// Birth
	particle_swarm<particle> b;

	std::function<VectorXd(const VectorXd&, const double&)> propagate;
	std::function<VectorXd(const VectorXd&, const Sensor&)> observe;

public:
	auto getQ() { return q; }
	auto getQPred() { return qPred; }
	auto getT() { return dt; }

	SMCJoTTFilter(std::function<VectorXd(const VectorXd&, const double&)> _propagate, 
		std::function<VectorXd(const VectorXd&, const Sensor&)> _observe, 
		const size_t& _nBirthComponents, const double & _birthWeight, const double & _pS, const double& _pB, const double & _q) 
		: nBirthComponents(_nBirthComponents), pS(_pS),  q(_q), pB(_pB)
	{
		propagate = _propagate;
		observe = _observe;
	}

	/*
	 * <summary> Prediction</summary>
	 */
	void predict(particle_swarm<particle>& _ps, Sensor& _sensor)
	{		
		double t, wb;

		// Probability of target existence prediction (28)
		qPred = pB * (1 - q) + pS * q;	

		t = pS * q / qPred;

		// Mean prediction
		for (auto &p : _ps.particles)
		{
			propagate(p.m, dt);
			p.w *= t;
		}
						
		// Birth
		wb = birthWeight * pB * (1 - q) / qPred;
		b = particle_swarm<particle>(nBirthComponents, _sensor.getSDim(), wb);

		for (size_t i = 0; i < nBirthComponents; i++)
		{
			b.particles[i].m << 1000 + rand() % 500 - 250, _sensor.getBearing(), 0, 0, 0;
			b.particles[i].m = Astro::razelToTEME(b.particles[i].m, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
		}
			
		// Concatenate predicted and birth components
		_ps = _ps + b;

	}

	void update(particle_swarm<particle>& _ps, Sensor& _sensor)
	{
		double nThreshold = 200;
		double cz = 1.0 / 57903 * 200;
		auto pD = _sensor.getPD();
		double I1 = _ps.weightSum() * pD, I2, wSum = 0, delta_k = 0, g = 0;
		std::vector<double> ll(_ps.size(), 0);		// Likelihoods
		size_t n = _ps.size();

		VectorXd z = VectorXd::Zero(3), zPred = VectorXd::Zero(3);
		
		// For each measurement
		for (size_t i = 0; i < _sensor.getZ().size(); i++) 
		{
			z = _sensor.getZ(i);
			n = _ps.size();
			I2 = 0;

			// For each particle
			for (size_t j = 0; j < _ps.particles.size(); j++)
			{
				zPred = observe(_ps.particles[j].m, _sensor);
				g = MathHelpers::gaussianLikelihood(z, zPred, _sensor.getS());	// Compute likelihood
				ll[j] += g;		// Likelihood sum

				I2 += pD * g  * _ps.particles[j].w;		// Approximate I2 integral (85)
			}

			delta_k += I2;
		}

		delta_k = I1 - delta_k / cz;		// delta_k approximation (86)
		q = (1 - delta_k) / (1 - q * delta_k) * q;

		// Update the weights (87)
		for (size_t i = 0; i < _ps.particles.size(); i++)
			_ps.particles[i].w = (1 - pD + pD * ll[i] / cz) * _ps.particles[i].w;
		
		// Normalize the weights
		wSum = _ps.weightSum;
		for (auto &p : _ps.particles)
			p.w /= wSum;

		// TODO: MCMC Move
		
		// Sampling importance re-sampling
		double nEff = 0;

		for (auto p : _ps.particles)
			nEff += p.w * p.w;

		nEff = 1 / nEff;

		// If the effective number of particles is less than a threshold, resample:
		if (nEff < nThreshold)
		{

		}

	}
};

