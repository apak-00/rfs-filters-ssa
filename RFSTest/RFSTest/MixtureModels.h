#pragma once
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <Eigen\Core>
#include <Eigen/LU>
#include "Sensor.h"

#include <iostream>

using namespace Eigen;

/*
* <summary> Mixture interface. </summary>
* <par> General interface for different mixtures (particle, Gaussian, etc.). </par>
*/
template <typename T>
struct mixture {

	mixture() : nMax(0), dimension(0) {};		// Empty constructor
	mixture(const size_t& _dim, const size_t& _nMax) : nMax(_nMax), dimension(_dim) {};
	mixture(const mixture& _m) : components(_m.components), nMax(_m.nMax), dimension(_m.dimension) {};		// Copy constructor

	std::vector<T> components;		// Vector with mixture components

	size_t nMax;					// Max. number of components, non-strict (?)
	size_t dimension;				// Dimensionality of the components

	/* Operator overloading */
	T operator [] (size_t i) const { return components[i]; }
	T& operator[] (size_t i) { return components[i]; }

	virtual void addComponent(const T& _c);

	/* Miscellaneous */
	auto size() const { return components.size(); };
	auto dim() const { return dimension; };
	double weightSum() const;
	void normalizeWeights();
};

/*
 * <summary> An interface for Gaussian Mixtures. </summary>
 */
template <typename T>
struct g_mixture : mixture<T>{

	unsigned int idCounter;			// A counter to track Gaussian components
	unsigned int trackCounter;		// A counter to track the tracks of the Gaussian components
	
	g_mixture();					// Empty constructor
	g_mixture(const size_t& _dim, const size_t& _nMax);
	g_mixture(const g_mixture& _gm);		// Copy constructor

	virtual void addComponent(const T& _c) override;
	virtual void merge(const double& _mergeThreshold) = 0;
	virtual void prune(const double& _pruneThreshold);
	virtual auto getEstimates(const double& _estimateThreshold) -> decltype(components);
};

template <typename T>
std::ostream& operator << (std::ostream& _os, const g_mixture<T>& _gm);

/**
 *	<summary> Gaussian component. </summary>
 */
struct gaussian_component {

	VectorXd m;		// Mean
	MatrixXd P;		// Covariance
	double w;		// Weight

	int tag[3];		// id, track id, track length

	bool kindaConverged;

	/* Constructors */
	gaussian_component();
	gaussian_component(const decltype(m)& _m, const decltype(P)& _P, const decltype(w)& _w, const int& _id);
	gaussian_component(const gaussian_component& _gc);

	/* Operator overloading */
	bool operator > (const gaussian_component& _gc) const { return w > _gc.w; }
	bool operator < (const gaussian_component& _gc) const { return w < _gc.w; }
	bool operator == (const gaussian_component _gc) const { return w == _gc.w; }
	bool operator >= (const gaussian_component& _gc) const { return w >= _gc.w; }
	bool operator <= (const gaussian_component& _gc) const { return w <= _gc.w; }

	gaussian_component operator+ (const gaussian_component& _gc) const;

	/* Miscellaneous */
	void initTag();
	void initTag(const int& id);

	static VectorXd getEmptyInfo() { return VectorXd::Zero(11); }
};

std::ostream& operator << (std::ostream& _os, const gaussian_component& _gc);

/**
 * <summary> Gaussian Mixture </summary>
 */
struct gaussian_mixture : g_mixture<gaussian_component> {
	 
	/* Constructors */
	gaussian_mixture();
	gaussian_mixture(const size_t& _dim, const size_t& _nMax);
	gaussian_mixture(const gaussian_mixture& _gm);

	/* GMM-related functions */
	virtual void merge(const double& _mergeThreshold);

	/* Other additional functions*/
	VectorXd getWeightsVector();
};

// Beta Gaussian Mixture

struct beta_gaussian_component : gaussian_component {

	double u, v;	// Beta distribution components
	bool mahler = false;

	/* Constructors */
	beta_gaussian_component();
	beta_gaussian_component(const decltype(m)& _m, const decltype(P)& _P, const decltype(w)& _w, const int& _id, const double& _u, const double& _v);
	beta_gaussian_component(const beta_gaussian_component& _bgc);

	beta_gaussian_component operator+ (const beta_gaussian_component& _bgc) const;

	static double getBetaMean(const double& _u, const double & _v);
	static double getBetaVariance(const double& _u, const double& _v, const double& _bMu);

	static VectorXd getEmptyInfo() { return VectorXd::Zero(14); }
};

std::ostream& operator << (std::ostream& _os, const beta_gaussian_component& _gc);

/**
 * <summary> Beta-Gaussian Mixture. </summary>
 */
struct beta_gaussian_mixture : g_mixture<beta_gaussian_component> {

	beta_gaussian_mixture();
	beta_gaussian_mixture(const size_t& _dim, const size_t& _nMax);
	beta_gaussian_mixture(const beta_gaussian_mixture& _gm);

	/* GMM-related functions */
	virtual void merge(const double& _mergeThreshold);

	static double betaHellinger(const beta_gaussian_component & _bgc1, const beta_gaussian_component & _bgc2);
};

/**
 * <summary> An empty constructor for the mixture interface. </summary>
 */
template<typename T>
inline g_mixture<T>::g_mixture() : mixture(), idCounter(0), trackCounter(0) {}

/**
 * <summary> A constructor that initializes an empty beta-Gaussian Mixture. </summary>
 * <param name = "_dim"> Dimensinality of the stored components. </param>
 * <param name = "_nMax"> Maximum number of components. </param>
 */
template<typename T>
inline g_mixture<T>::g_mixture(const size_t & _dim, const size_t & _nMax)
	: mixture(_dim, _nMax), idCounter(0), trackCounter(0) {}

/**
 * <summary> A copy constructor of the Mixture interface. </summary>
 * <param name = "_m"> A mixture to copy from. </param>
 */
template<typename T>
inline g_mixture<T>::g_mixture(const g_mixture & _mixture) : mixture(_mixture),
	idCounter(_mixture.idCounter), trackCounter(_mixture.trackCounter) {}

/**
* <summary> Adds a component to the mixture. </summary>
* <param name = "_c"> A component to be added. </param>
*/
template<typename T>
inline void mixture<T>::addComponent(const T & _component)
{
	assert(_component.m.size() == dimension && "Attempt to add a mixture component of different dimension ");
	components.push_back(_component);
}

/**
 * <summary> Adds a component to the mixture. </summary>
 * <param name = "_c"> A component to be added. </param>
 */
template<typename T>
inline void g_mixture<T>::addComponent(const T & _component)
{
	assert(_component.m.size() == dimension && "Attempt to add a mixture component of different dimension ");

	components.push_back(_component);

	if (!_component.tag[0])
		components.back().tag[0] = idCounter++;
}

/**
* <summary> Pruning of the Gaussian Components with weights under the threshold. </summary>
* <param name = _weightThreshold> Pruning weight threshold. </param>
*/
template<typename T>
inline void g_mixture<T>::prune(const double & _pruneThreshold)
{
	auto pruned = std::remove_if(components.begin(), components.end(),
		[&_pruneThreshold](const T& c) { return c.w < _pruneThreshold; });
	components.erase(pruned, components.end());

	if (components.size() > nMax) 
		components.resize(nMax);

	double wSum = 0;

	for (auto &c : components)	 
		wSum += c.w;

	for (auto &c : components)
		c.w /= wSum;
}

/**
 * <summary> Obtain weight-based estimates from the mixture. </summary>
 * <param name = "_estimateThreshold"> Weight threshold. </param>
 * <returns> A vector of components with the weight above threshold. </returns>
 */
template<typename T>
inline auto g_mixture<T>::getEstimates(const double & _estimateThreshold) -> decltype(components)
{
	auto below = std::remove_if(components.begin(), components.end(),
		[&_estimateThreshold](const T& gc) { return gc.w < _estimateThreshold; });

	for (auto i = components.begin(); i != below; i++) {
		if (i->tag[1] == 0)
			i->tag[1] = trackCounter++;
		i->tag[2]++;
	}

	return std::vector<T>(components.begin(), below);
}

/**
 * <summary> Get sum of the weights of all components of the mixture. </summary>
 * <returns> Sum of the weights of the components. </returns>
 */
template<typename T>
inline double mixture<T>::weightSum() const
{
	double sum = 0;
	for (auto &c : components)
		sum += c.w;
	//double sum = std::accumulate(components.begin(), components.end(), 0.0, [](double _sum, const T& _p) { return _sum + _p.w; });
	return sum;
}

template<typename T>
inline void mixture<T>::normalizeWeights()
{
	double wSum = weightSum();
	if (wSum != 1)
		for (auto &c : components)
			c.w /= wSum;
}

/**
* <summary> Stream operator overloading for Gaussian Component. </summary>
* <param name = "_os"> A reference to the output stream. </param>
* <param name = "_gm"> Gaussian Mixture for the output. </param>
*/
template<typename T>
inline std::ostream & operator<<(std::ostream & _os, const g_mixture<T> & _m)
{
	// TODO: insert return statement here
	_os << _m.size() << "/" << _m.nMax << std::endl;

	for (auto c : _m.components)
		_os << c << std::endl;

	return _os;
}

/*
 * <summary> Gaussian Bernoulli model for smoother. <summary>
 * <date> October 14, 2016 </date>
 */
template <typename T>
struct gaussian_bernoulli_model
{
	T gmixture;				// Gaussian Mixture
	double r;				// Probability of target existence (q)
	double rPredicted;		// Predicted probability of target existence (qPred)
};

/*
* <summary> Particle for SMC. </summary>
*/
struct particle
{
	VectorXd m;			// Mean
	double w;			// Weight

	particle();
	particle(const size_t& dim, const double& _w = 0);
	particle(const VectorXd& _m, const double& _w = 0);
};

/*
* <summary> Particle swarm for SMC. </summary>
*/
struct particle_mixture : mixture <particle>
{
	static std::random_device rDev;
	static std::mt19937 generator;
	static std::uniform_real_distribution<double> itsd;

	particle_mixture();
	particle_mixture(const particle_mixture& _pm);
	particle_mixture(const size_t& _n, const size_t& _dim, const double& _w = 0);

	particle_mixture operator+ (const particle_mixture& _gc) const;

	double getEffectiveN();
	void resampleITS(const size_t& _size);
	void populateRandomRAZEL(const Sensor& _sensor, const VectorXd& _lBound, const VectorXd& _uBound);
};