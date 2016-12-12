#include <algorithm>
#include "MixtureModels.h"
#include "MathHelpers.h"

// Temp
#include <iostream>

/**
 * <summary> An empty constructor fo the Gaussian Mixture. </summary>
 */
gaussian_mixture::gaussian_mixture() : g_mixture() {}

/**
 * <summary> Constructor that initializes an empty Gaussian Mixture of speified dimension. </summary>
 * <param name = "_dim"> Dimensionality of the stored Gaussian components. </param>
 * <param name = "_nMax"> Maximum number of the Gaussian componentes. </param>
 */
gaussian_mixture::gaussian_mixture(const size_t & _dim, const size_t & _nMax)
	: g_mixture(_dim, _nMax) {}

/**
 * <summary> Copy constructor of the Gaussian Mixture. </summary>
 * <param name = "_gm"> An instance of Gaussian Mixture to copy from. </param>
 */
gaussian_mixture::gaussian_mixture(const gaussian_mixture & _gm)
	: g_mixture(_gm) {}

/**
 * <summary> Merge of the Gaussian Components within threshold range. </summary>
 * <param name = "_mergeThreshold"> Merging distance threshold. </param>
 */
void gaussian_mixture::merge(const double & _mergeThreshold)
{
	std::vector<gaussian_component> temp;
	
	while (components.size() > 0) {

		// Find the component with the maximum weight
		auto max = std::max_element(components.begin(), components.end(),
			[](const gaussian_component& a, const gaussian_component& b) { return a.w < b.w; });

		temp.push_back(*max);
		components.erase(max);

		for (auto i = components.begin(); i != components.end();) {

			double distance = MathHelpers::hellinger(temp.back(), *i);

			if (distance < _mergeThreshold) {
				temp.back() = temp.back() + *i;
				// The component with the highest weight keeps the track (???)
				if (temp.back().tag[2] < i->tag[2]) {
					//temp.back().tag[1] = i->tag[1];
					temp.back().tag[2] = i->tag[2];
				}

				i = components.erase(i);
			}
			else
				i++;
		}
	}

	components = temp;
}

 VectorXd gaussian_mixture::getWeightsVector()
 {
	 VectorXd result(components.size());

	 for (size_t i = 0; i < components.size(); i++)
		 result(i) = components[i].w;
	 return result;
 }

/* ----------------------------------------------------------------------- */
/* --------------------- Gaussian Component ------------------------------ */
/* ----------------------------------------------------------------------- */

/**
 * <summary> Default constructor. </summary>
 * <par> Initializes an empty component. </par>
 */
gaussian_component::gaussian_component()
{
	m = VectorXd();
	P = MatrixXd();
	w = 0;
	initTag();
	kindaConverged = false;
}

/**
 * <summary> Constructor that takes all of the parameters. </summary>
 * <param name = "_m"> The mean. </param>
 * <param name = "_P"> The covariance. </param>
 * <param name = "_w"> The weight. </param>
 * <param name = "_id"> The id. </param>
 */
gaussian_component::gaussian_component(const decltype(m)& _m, const decltype(P)& _P, const decltype(w)& _w, const int& _id)
	: m(_m), P(_P), w(_w), kindaConverged(false)
{
	initTag(_id);
}

/**
 * <summary> Copy constructor. </summary>
 * <param> An instance of Gaussian component to copy from. </param>
 */
gaussian_component::gaussian_component(const gaussian_component & _gc)
	: m(_gc.m), P(_gc.P), w(_gc.w), kindaConverged(_gc.kindaConverged)
{
	for (size_t i = 0; i < 3; i++)
		tag[i] = _gc.tag[i];
}

/**
 * <summary> Addition operator overload. </summary>
 * <param name = "_gc"> Added component. </param> 
 */
gaussian_component gaussian_component::operator+(const gaussian_component & _gc) const
{
	assert(m.size() == _gc.m.size() && "Terms are of different dimensions.");

	gaussian_component result(*this);
	VectorXd d = m - _gc.m;

	result.w = w + _gc.w;
	result.m = (m * w + _gc.m * _gc.w) / result.w;
	result.P = (P * w + _gc.P * _gc.w) / result.w + d * d.transpose() * w * _gc.w;

	result.kindaConverged = kindaConverged || _gc.kindaConverged;

	return result;
}

/**
 * <summary> Initializes an empty tag comp onent. </summary>
 */
void gaussian_component::initTag()
{
	tag[0] = 0;
	tag[1] = 0;
	tag[2] = 0;
}

/**
 * <summary> Initializes a tag component according to the specified id. </summary>
 * <param name = "_id"> New of the component. </param>
 */
void gaussian_component::initTag(const int& _id) {
	tag[0] = _id;
	tag[1] = 0;
	tag[2] = 0;
}

/**
 * <summary> Stream operator overloading for Gaussian Component. </summary>
 * <param name = "_os"> A reference to the output stream. </param>
 * <param name = "_gc"> Gaussian component for the output. </param> 
 */
std::ostream & operator<<(std::ostream & _os, const gaussian_component & _gc)
{
	for (size_t i = 0, iSize = _gc.m.size(); i <= iSize; i++)
		_os << _gc.m(i) << ",";
	
	_os << _gc.w << "," << _gc.tag[1] << "," 
		<< _gc.P.determinant() << "," << _gc.P.block<3, 3>(0, 0).determinant() << "," << _gc.P.block<3, 3>(3, 3).determinant();

	return _os ;
}

std::ostream & operator<<(std::ostream & _os, const beta_gaussian_component & _gc)
{
	for (size_t i = 0, iSize = _gc.m.size(); i <= iSize; i++)
		_os << _gc.m(i) << ",";

	_os << _gc.w << "," << _gc.tag[1] << ","
		<< _gc.P.determinant() << "," << _gc.P.block<3, 3>(0, 0).determinant() << "," << _gc.P.block<3, 3>(3, 3).determinant() << ","
		<< _gc.u << "," << _gc.v << "," << _gc. u / (_gc.u + _gc.v);

	return _os;
}

/*
 * <summary> An empty constructor for the beta-Gaussian component structure. </summary>
 * <par> Initializes everything to zero, and u and v parameters of beta distribution to one. </par>
 */
beta_gaussian_component::beta_gaussian_component() : gaussian_component()
{
	u = 1;
	v = 1;
}

/**
 * <summary> Default constructor for the beta-Gaussian component structure. </summary>
 * <param name = "_m"> The mean. </param>
 * <param name = "_P"> The covariance. </param>
 * <param name = "_w"> The weight. </param>
 * <param name = "_id"> The id. </param>
 * <param name = "_u"> First parameter of the beta distribution. </param>
 * <param name = "_v"> Second parameter of the beta distribution. </param>
 */
beta_gaussian_component::beta_gaussian_component(const decltype(m)& _m, const decltype(P)& _P, const decltype(w)& _w, const int & _id, 
	const double & _u, const double & _v) : gaussian_component (_m, _P, _w, _id)
{
	u = _u;
	v = _v;
}

/**
 * <summary> Copy constructor for the beta-Gaussian component structure. </summary>
 * <param name = "_bgc"> A beta-Gaussian component to copy from. </param>
 */
beta_gaussian_component::beta_gaussian_component(const beta_gaussian_component & _bgc) : gaussian_component(_bgc)
{
	u = _bgc.u;
	v = _bgc.v;
}

/**
 * <summary> Addition operator overload for beta-Gaussian component. </summary>
 * <param name = "_bgc"> A constant reference to the beta-Gaussian component to add. </param>
 */
beta_gaussian_component beta_gaussian_component::operator+(const beta_gaussian_component & _bgc) const
{
	// I was confused how to call the parent operator properly and still preserve access to inital variables.
	// So the first part is simply copied from the gaussian_component + operator.

	assert(m.size() == _bgc.m.size() && "Terms are of different dimensions.");

	beta_gaussian_component result(*this);
	VectorXd d = m - _bgc.m;

	result.w = w + _bgc.w;
	result.m = (m * w + _bgc.m * _bgc.w) / result.w;
	result.P = (P * w + _bgc.P * _bgc.w) / result.w + d * d.transpose() * w * _bgc.w;

	result.kindaConverged = kindaConverged || _bgc.kindaConverged;

	// Beta component part

	// Do it according to the Mahler's book or in the simple manner
	if (mahler) 
	{
		double theta0, v0, v1, v2, mu0, mu1, mu2;

		mu1 = getBetaMean(u, v);
		mu2 = getBetaMean(_bgc.u, _bgc.v);
		v1 = getBetaVariance(u, v, mu1);
		v2 = getBetaVariance(_bgc.u, _bgc.v, mu2);

		mu0 = (w * mu1 + _bgc.w * mu2) / result.w;
		// TODO: Check the typo (?)
		v0 = -mu0 * mu0 + (w * (mu1 * mu1 + v1) + _bgc.w * (mu2 * mu2 + v2)) / result.w;

		theta0 = mu0 * (1 - mu0) / v0 - 1;

		result.u = theta0 * mu0;
		result.w = theta0 * (1 - mu0);
	}
	// Assign the beta parameters of the component with the highest weight
	else
	{
		result.u = w * u + _bgc.w * _bgc.u;
		result.v = w * v + _bgc.w * _bgc.v;

		/*
		if (w > _bgc.w)
		{
			result.u = u;
			result.v = v;
		}
		else 
		{
			result.u = _bgc.u;
			result.v = _bgc.v;
		}
		*/
	}

	return result;
}

/**
 * <summary> Obtain the mean of the beta component with the specified parameters. </summary>
 * <param name = "_u"> a (alpha) parameter of the beta distribution. </param>
 * <param name = "_v"> b (beta) parameter of the beta distribution. </param>
 * <returns> The mean of the beta distribution. </returns>
 */
double beta_gaussian_component::getBetaMean(const double & _u, const double & _v)
{
	return _u / (_u + _v);
}

/**
* <summary> Obtain the variance of the beta component with the specified parameters. </summary>
* <param name = "_u"> a (alpha) parameter of the beta distribution. </param>
* <param name = "_v"> b (beta) parameter of the beta distribution. </param>
* <param name = "_bMu"> Pre-computed mean of the beta distribution. </param>
* <returns> The variance of the beta distribution. </returns>
*/
double beta_gaussian_component::getBetaVariance(const double & _u, const double & _v, const double & _bMu)
{
	return _bMu * (1 - _bMu) / (_u + _v + 1);
}

/**
 * <summary> An empty constructor of the beta-Gaussian Mixture. </summary>
 */
beta_gaussian_mixture::beta_gaussian_mixture() : g_mixture() {}

/**
 * <summary >A constructor of the Beta Gaussian Mixture that takes the dimensionality of the stored components 
 * and the maximumum number of componentes as the parameters of the filter. </summary>
 * <param name = "_dim"> Dimensionality of the stored components. </param>
 * <param name = "_nMax"> Maximum number of components in the mixture. </param>
 */
beta_gaussian_mixture::beta_gaussian_mixture(const size_t & _dim, const size_t & _nMax) : g_mixture(_dim, _nMax) {}

/**
* <summary> Copy constructor of the beta-Gaussian Mixture. </summary>
* <param name = "_gm"> An instance of beta-Gaussian Mixture to copy from. </param>
*/
beta_gaussian_mixture::beta_gaussian_mixture(const beta_gaussian_mixture & _bgm) : g_mixture(_bgm) {}

/**
 * <summary> Merge procedure of the Beta Gaussian Mixture. </summary>
 * <par> Hellinger distance is used as a merging criteria. </par>
 * <param name = "_mergeThreshold"> Merge threshold</param>
 */
void beta_gaussian_mixture::merge(const double & _mergeThreshold)
{
	std::vector<beta_gaussian_component> temp;

	while (components.size() > 0) {

		// Find the component with the maximum weight
		auto max = std::max_element(components.begin(), components.end(),
			[](const gaussian_component& a, const gaussian_component& b) { return a.w < b.w; });

		temp.push_back(*max);
		components.erase(max);

		for (auto i = components.begin(); i != components.end();) {

			double distance = MathHelpers::hellinger(temp.back(), *i);

			if (distance < _mergeThreshold) {
				temp.back() = temp.back() + *i;
				// The component with the highest weight keeps the track (???)
				if (temp.back().tag[2] < i->tag[2]) {
					//temp.back().tag[1] = i->tag[1];
					temp.back().tag[2] = i->tag[2];
				}

				i = components.erase(i);
			}
			else
				i++;
		}
	}

	components = temp;
}

/*
 * TODO: To be implemented
 * <summary> Calculates Hellinger distance between Beta-Gaussian components. </summary>
 */
double beta_gaussian_mixture::betaHellinger(const beta_gaussian_component & _bgc1, const beta_gaussian_component & _bgc2)
{
	return 0.0;
}

// ---------- Particles ----------
/*
 * <summary> Constructors for particle. </summary>
 */
particle::particle() : w(0) {}		// Empty constructor
particle::particle(const size_t& _dim, const double& _w) : m(VectorXd::Zero(_dim)), w(_w) {}	// Zero constructor with dimenstion and weight
particle::particle(const VectorXd& _m, const double& _w) : m(_m), w(_w) {}						// Constructor with mean and weight

/*
 * <summary> Empty constructor. </summary>
 */
particle_mixture::particle_mixture() : mixture(){}

/*
 * <summary> Copy constructor. </summary>
 */
particle_mixture::particle_mixture(const particle_mixture & _pm) : mixture(_pm) {}

/*
* <summary> Default constructor. (?) </summary>
*/
particle_mixture::particle_mixture(const size_t & _nMax, const size_t & _dim, const double & _w) : mixture(_dim, _nMax)
{
	particle p(_dim, _w);
	components = std::vector<particle>(_nMax, p);
}

/*
 * <summary> Addition operator overload for particle swarm. <summary>
 */
particle_mixture particle_mixture::operator+(const particle_mixture& _pc) const
{
	// TODO: Optimize?
	particle_mixture result(*this);
	result.components.insert(std::end(components), std::begin(_pc.components), std::end(_pc.components));
	return result;
}
// Define the static variables
std::random_device particle_mixture::rDev;
std::mt19937 particle_mixture::generator = std::mt19937(rDev());
std::uniform_real_distribution<double> particle_mixture::itsd = std::uniform_real_distribution<double>(0.0, 1.0);

/*
* <summary> Returns effective number of particles. </summary>
*/
double particle_mixture::getEffectiveN()
{
	double result = 0.0;
	for (auto p : components)
		result += p.w * p.w;

	return 1.0 / result;
}

/**
* <summary> Inverse transform sampling </summary>
*/
void particle_mixture::resampleITS(const size_t& _size)
{
	double rIdx;
	std::vector<double> cdf(components.size());		// Cumulative distribution function
	std::vector<particle> resampled(_size);

	if (weightSum() != 1)
		normalizeWeights();

	// Sort elements in the descending order
	std::sort(components.begin(), components.end(), [](particle a, particle b) {return b.w < a.w; });

	// Calculate the CDF
	cdf[0] = components[0].w;
	for (size_t i = 1; i < components.size(); i++)
		cdf[i] += cdf[i - 1] + components[i].w;

	// Inverse Transform Sampling
	for (size_t i = 0; i < _size; i++)
	{
		rIdx = itsd(generator);				// Generate a random number between 0 and 1
											// TODO: Optimize
											// Get the closest particle index from the cdf
		auto closest = std::min_element(cdf.begin(), cdf.end(),
			[rIdx](double x, double y) {return abs(x - rIdx) < abs(y - rIdx); });
		resampled[i] = *closest;
	}

	components = resampled;
}

void particle_mixture::populateRandomRAZEL(const Sensor& _sensor, const VectorXd& _lBound, const VectorXd& _uBound)
{
	if (components.size() == 0)
		return;

	VectorXd randRAZEL = VectorXd::Zero(6), randTEME = VectorXd::Zero(6);
	double weight = 1.0 / components.size();

	std::uniform_real_distribution<double> dr(_lBound(0), _uBound(0));		// Distribution for Range
	std::uniform_real_distribution<double> da(_lBound(1), _uBound(1));		// Distribution for Azimuth
	std::uniform_real_distribution<double> de(_lBound(2), _lBound(2));		// Distribution for Elevation

	for (size_t i = 0; i < components.size(); i++)
	{
		randRAZEL << dr(generator), da(generator), de(generator), 0, 0, 0;
		randTEME = Astro::razelToTEME(randRAZEL, _sensor.getPosition(), _sensor.getDateJD(),
			_sensor.getLOD(), _sensor.getXp(), _sensor.getYp());

		components[i].m = randTEME;
		components[i].w = weight;
	}
}

VectorXd particle_mixture::getWeightedAverage()
{
	normalizeWeights();

	VectorXd result = VectorXd::Zero(6);

	for (auto c : components)
		result += c.m * c.w;

	return result;
}

