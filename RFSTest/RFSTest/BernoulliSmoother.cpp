#include "BernoulliSmoother.h"

BernoulliSmoother::BernoulliSmoother() {}

void BernoulliSmoother::smooth(gaussian_bernoulli_model<gaussian_mixture> _gbmPrior, gaussian_bernoulli_model<gaussian_mixture> _gbmPosterior)
{
	double alphaR, betaR, alphaS, betaS, tempOne, tempTwo, smoothedR;

	tempOne = (1 - rCond) / (1 - _gbmPosterior.rPredicted);
	tempTwo = rCond / _gbmPosterior.rPredicted;

	alphaR = (1 - pB) * tempOne;
	betaR = pB * tempTwo;
	alphaS = (1 - pS) * tempOne;
	betaS = pS * tempTwo;

	// TODO:
	// Implement 23, 43, 45;
}
