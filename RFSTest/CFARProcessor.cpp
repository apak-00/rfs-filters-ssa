#include "CFARProcessor.h"

// TODO: To be implemented

/**
 * <summary> Default constructror. </summary>
 */
CFARProcessor::CFARProcessor() : falsePositiveRate(1e-4), innerWindowRadius(0.5), outerWindowRadius(5.0), 
minRange(1.0), maxRange(200.0) {}

/*
 * <summary> Processes 1D radar scan. </summary>
 * <param name = "_scan"> Vector containing the source data. </param>
 * <returns>  </returns>
 */
std::vector<double> CFARProcessor::processCACFAR1D(const std::vector<double>& _scan)
{
	// Auxilary variables
	double rangeResolution = (maxRange - minRange) / _scan.size();
	const int innerWindow = (int)(innerWindowRadius / rangeResolution);
	const int outerWindow = (int)(outerWindowRadius / rangeResolution);
	const int windowFactor = (outerWindow - outerWindow) * 2;
	const double tau = windowFactor * (pow(falsePositiveRate, -1.0 / windowFactor) - 1);

	// The output pointcloud
	std::vector<double> processedScan;
	std::vector<double> accumulatedIntensities(_scan.size());

	// Pre-accumulate intensities. This is similar to the integral images approach
	accumulatedIntensities[0] = pow(10, _scan[0] / 10);
	for (size_t j = 1; j < _scan.size(); j++)
		accumulatedIntensities[j] = accumulatedIntensities[j - 1] + pow(10, _scan[j] / 10);
	
	// Iterate for each possible point
	for (int j = outerWindow; j < _scan.size() - outerWindow; j++)
	{
		double point_real_power = pow(10, _scan[j] / 10);

		// With pre-accumulated intensities, an average becomes three sums
		double window_average = accumulatedIntensities[j + outerWindow]
			- accumulatedIntensities[j - outerWindow]
			- accumulatedIntensities[j + innerWindow]
			+ accumulatedIntensities[j - innerWindow];
		window_average = window_average / (2 * (outerWindow - innerWindow));

		// Check if the current point is above the predefined threshold
		if (point_real_power > tau * window_average)
		{
			const double base = falsePositiveRate / (windowFactor * (1 + point_real_power / window_average)) + 1;
			double intensity = pow(base, -windowFactor);
			processedScan.push_back(intensity);
		}
	}
	
	return processedScan;
}
