// TODO: Complete

#pragma once
#include <Eigen/Dense>
#include "Sensor.h"
#include "MixtureModels.h"
#include "MathHelpers.h"

class BernoulliSmoother
{
protected:
	double pS;						// Probability of target survival
	double pB;						// Probabilirt of target birth
	double rCond;

public:
	BernoulliSmoother();

	void smooth(gaussian_bernoulli_model<gaussian_mixture> _gbmPrior, gaussian_bernoulli_model<gaussian_mixture> _gbmPosterior);
};

