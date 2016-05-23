#pragma once
#include <vector>
class CFARProcessor
{
private:

	int operation_mode_;			// Mode of operation
	double fpr_;					// Expected false positive rate
	double cloud_min_range_;		// Closest distance
	double cloud_max_range_;		// Far distance
	double inner_window_radius_;	// Radius of the inner area around the interest point
	double outer_window_radius_;	// Radius of the outer area around the interest point
	double angle_granularity_;		// Granularity of the angle
	double range_granularity_;		// Granularity of the range

public:
	CFARProcessor();
	~CFARProcessor();

	std::vector<double> process_cloud_1D(const std::vector<double>& _scan);
};

