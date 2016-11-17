#include <algorithm>
#include "gmm.h"

// Temp
#include <iostream>

/**
 * <summary> An empty constructor fo the Gaussian Mixture. </summary>
 */
gaussian_mixture::gaussian_mixture() : mixture() {}

/**
 * <summary> Constructor that initializes an empty Gaussian Mixture of speified dimension. </summary>
 * <param name = "_dim"> Dimensionality of the stored Gaussian components. </param>
 * <param name = "_nMax"> Maximum number of the Gaussian componentes. </param>
 */
gaussian_mixture::gaussian_mixture(const size_t & _dim, const size_t & _nMax)
	: mixture(_dim, _nMax) {}

/**
 * <summary> Constructor that initializes a Gaussian Mixture with n random elements. </summary>
 * <param name = "_dim"> Dimensionality of the stored Gaussian components. </param>
 * <param name = "_nMax"> Maximum number of the Gaussian componentes. </param>
 * <param name = "_lBound"> Lower state bound (for random birth). </param>
 * <param name = "_uBound"> Upper state bound (for random birth). </param>
 * <param name = "_iCov"> Initial covariance matrix for the new components. </param>
 * <param name = "_iWeight"> Initial weight for the new components. </param>
 */
gaussian_mixture::gaussian_mixture(const size_t & _dim, const size_t & _n, const size_t & _nMax,
	const VectorXd& _lBound, const VectorXd& _uBound, const MatrixXd& _iCov, const double& _iWeight)
	: mixture(_dim, _n, _nMax)
{
	for (size_t i = 0; i < components.size(); i++) 
		components[i] = gaussian_component(randVec(_lBound, _uBound), _iCov, _iWeight, idCounter++);
}

/**
 * <summary> Copy constructor of the Gaussian Mixture. </summary>
 * <param name = "_gm"> An instance of Gaussian Mixture to copy from. </param>
 */
gaussian_mixture::gaussian_mixture(const gaussian_mixture & _gm)
	: mixture(_gm) {}

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

			double distance = hellinger(temp.back(), *i);

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
	for (size_t i = 0; i < _gc.m.size(); i++)
		_os << _gc.m(i) << ",";
	
	_os << _gc.w << "," << _gc.tag[1] << "," 
		<< _gc.P.determinant() << "," << _gc.P.block<3, 3>(0, 0).determinant() << "," << _gc.P.block<3, 3>(3, 3).determinant();

	return _os ;
}

std::ostream & operator<<(std::ostream & _os, const beta_gaussian_component & _gc)
{
	for (size_t i = 0; i < _gc.m.size(); i++)
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
beta_gaussian_mixture::beta_gaussian_mixture() : mixture() {}

/**
 * <summary >A constructor of the Beta Gaussian Mixture that takes the dimensionality of the stored components 
 * and the maximumum number of componentes as the parameters of the filter. </summary>
 * <param name = "_dim"> Dimensionality of the stored components. </param>
 * <param name = "_nMax"> Maximum number of components in the mixture. </param>
 */
beta_gaussian_mixture::beta_gaussian_mixture(const size_t & _dim, const size_t & _nMax) : mixture(_dim, _nMax) {}

/*
 * <summary> Main constructor of the beta-Gaussian Mixture. </summary>
 * <par> Initializes a beta-Gaussian mixture with a specified number of random beta-Gaussian components. 
 * The alpha and beta parameters of the beta distribution are equal to one. </par>
 *
 * <param name = "_dim"> Dimensionality of the stored beta-Gaussian components. </param>
 * <param name = "_nMax"> Maximum number of the beta-Gaussian componentes. </param>
 * <param name = "_lBound"> Lower state bound (for random birth). </param>
 * <param name = "_uBound"> Upper state bound (for random birth). </param>
 * <param name = "_iCov"> Initial covariance matrix for the new components. </param>
 * <param name = "_iWeight"> Initial weight for the new components. </param>
 */
beta_gaussian_mixture::beta_gaussian_mixture(const size_t & _dim, const size_t & _n, const size_t & _nMax, 
	const VectorXd & _lBound, const VectorXd & _uBound, const MatrixXd & _iCov, const double & _iWeight)
	: mixture(_dim, _n, _nMax)
{
	for (size_t i = 0; i < components.size(); i++)
		components[i] = beta_gaussian_component(randVec(_lBound, _uBound), _iCov, _iWeight, idCounter++, 1, 1);
}

/**
* <summary> Copy constructor of the beta-Gaussian Mixture. </summary>
* <param name = "_gm"> An instance of beta-Gaussian Mixture to copy from. </param>
*/
beta_gaussian_mixture::beta_gaussian_mixture(const beta_gaussian_mixture & _bgm) : mixture (_bgm)
{
}

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

			double distance = hellinger(temp.back(), *i);

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
 * <summary> Constructors for particle swarm. </summary>
 */
template<typename T>
particle_swarm<T>::particle_swarm() {}

template<typename T>
particle_swarm<T>::particle_swarm(const size_t & _n, const size_t & _dim, const double & _w)
{
	particle p(_dim, _w);
	particles = std::vector(_n, p);
}

template<typename T>
particle_swarm<T>::particle_swarm(const particle_swarm<T>& _pc) : particles(_pc.particles) {}

/*
 * <summary> Addition operator overload for particle swarm. <summary>
 */
template<typename T>
particle_swarm<T> particle_swarm<T>::operator+(const particle_swarm<T>& _pc) const
{
	// TODO: Optimize?
	particle_swarm<T> result(this);
	result.particles.insert(std::end(particles), std::begin(_pc.particles), std::end(_pc.particles));
	return result;
}

/*
 * <summary> Returns the size of the particles vector. </summary> 
 */
template<typename T>
size_t particle_swarm<T>::size()
{
	return particles.size();
}

/*
 * <summary> Returns the sum of the weights of the particles. </summary>
 */
template<typename T>
double particle_swarm<T>::weightSum()
{
	double result = std::accumulate(particles.begin(), particles.end(), 0.0, [](double _sum, const T& _p) { return _sum + _p.w; });
	//for (auto p : particles)
	//	result += p.w;
	return result;
}

/*
 * <summary> Normalizes the weights of the particles. </summary>
 */
template<typename T>
void particle_swarm<T>::normalize()
{
	double wSum = weightSum();
	if (wSum != 1)
		for (auto &p : particles)
			p.w /= wSum;
}

/**
 * <summary> Inverse transform sampling </summary>
 */
template<typename T>
void particle_swarm<T>::resampleITS(const size_t& _size)
{
	double rIdx;
	std::vector<double> cdf(particles.size());		// Cumulative distribution function
	vector<T> resampled(_size);

	if (weightSum() != 1)
		normalize();

	// Sort elements in the descending order
	std::sort(particles.begin(), particles.end(), [](T a, T b) {return b.w < a.w; });

	// Calculate the CDF
	cdf[0] = particles[0].w;
	for (size_t i = 1; i < particles.size(); i++)
		cdf[i] += cdf[i - 1] + particles[i].w;

	// Inverse Transform Sampling
	// Random number generator;
	std::random_device rDev;
	std::mt19937 generator(rDev());
	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	for (size_t i = 0; i < _size; i++)
	{
		rIdx = distribution(generator);		// Generate a random number between 0 and 1
		// TODO: Optimize
		// Get the closest particle index from the cdf
		auto closest = std::min_element(particles.begin(), particles.end(),
			[rIdx](T x, T y) {return abs(x.w - rIdx) < abs(y.w - tIdx); });
		resampled[i] = *closest;
	}

	particles = resampled;
}
