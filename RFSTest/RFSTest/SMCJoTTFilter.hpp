#pragma once
#include "RFSFilter.h"
#include <Eigen/Dense>
#include "Sensor.h"
#include "gmm.h"

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
	std::vector<double> birthRanges;

	std::function<VectorXd(VectorXd&, double)> propagate;
	std::function<VectorXd(VectorXd, Sensor)> observe;

public:
	auto getQ() { return q; }
	auto getQPred() { return qPred; }
	auto getT() { return dt; }

	SMCJoTTFilter(std::function<VectorXd(VectorXd, double)> _propagate, const size_t& _nBirthComponents, const double & _birthWeight,
		const double & _pS, const double& _pB, const double & _q) : nBirthComponents(_nBirthComponents), pS(_pS),  q(_q), pB(_pB),
		b(_nBirthComponents)
	{
		propagate = _propagate;
		birthRanges.resize(nBirthComponents);
	}

	void predict(particle_swarm<particle>& _ps, Sensor& _sensor)
	{
		VectorXd randomBirth = VectorXd::Zero(6);
		
		qPred = pB * (1 - q) + pS * q;		// Probability of target existence prediction

		for (auto &p : _ps.particles)		// Mean prediction
			propagate(p.m, dt);				

		for (size_t i = 0; i < nBirthComponents; i++)		// Birth
		{
			birthRanges[i];
		}
	}
};

