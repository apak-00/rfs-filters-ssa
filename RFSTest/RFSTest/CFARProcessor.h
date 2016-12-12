#pragma once
#include <vector>

// TODO: Yet to be implemented

/**
 * <summary> Cell-averaging constant false alarm rate processor. </summary>
 */
class CFARProcessor
{
public:
	CFARProcessor();

	std::vector<double> processCACFAR1D(const std::vector<double>& _scan);

private:
	double falsePositiveRate;	// Expected false positive rate
	double minRange;			// Closest distance
	double maxRange;			// Far distance
	double innerWindowRadius;	// Radius of the inner area around the interest point
	double outerWindowRadius;	// Radius of the outer area around the interest point
};

