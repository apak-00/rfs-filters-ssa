#include "CFARProcessor.h"



CFARProcessor::CFARProcessor()
{
	fpr_ = 1E-4;
	inner_window_radius_ = 0.5;
	outer_window_radius_ = 5.0;
	cloud_min_range_ = 1.0;
	cloud_max_range_ = 200.0;
	angle_granularity_ = 0.0290077572003;
	range_granularity_ = 0.471542968784;
	operation_mode_ = -1;
}


CFARProcessor::~CFARProcessor()
{
}

std::vector<double> CFARProcessor::process_cloud_1D(const std::vector<double>& spectrum_cloud)
{
	double range_resolution = (cloud_max_range_ - cloud_min_range_)
		/ spectrum_cloud.size();
	const int inner_window = (int)(inner_window_radius_ / range_resolution);
	const int outer_window = (int)(outer_window_radius_ / range_resolution);
	const int window_factor = (outer_window - inner_window) * 2;
	const double tau = window_factor * (pow(fpr_, -1.0 / window_factor) - 1);

	// The output pointcloud
	std::vector<double> processed_cloud;
	std::vector<double> accum_intensities(spectrum_cloud.size());

	// Pre-accumulate intensities. This is similar to the integral images approach
	accum_intensities[0] = pow(10, spectrum_cloud[0] / 10);
	for (size_t j = 1; j < spectrum_cloud.size(); j++)
		accum_intensities[j] = accum_intensities[j - 1] + pow(10, spectrum_cloud[j] / 10);
	
	// Iterate for each possible point
	for (int j = outer_window; j < spectrum_cloud.size() - outer_window; j++)
	{
		double point_real_power = pow(10, spectrum_cloud[j] / 10);

		// With pre-accumulated intensities, an average becomes three sums
		double window_average = accum_intensities[j + outer_window]
			- accum_intensities[j - outer_window]
			- accum_intensities[j + inner_window]
			+ accum_intensities[j - inner_window];
		window_average = window_average / (2 * (outer_window - inner_window));

		// Check if the current point is above the predefined threshold
		if (point_real_power > tau * window_average)
		{
			const double base = fpr_ / (window_factor*(1 + point_real_power / window_average)) + 1;
			double intensity = pow(base, -window_factor);
			processed_cloud.push_back(intensity);
		}
	}
	
	return processed_cloud;
}
