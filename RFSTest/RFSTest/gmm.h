#pragma once
#include <vector>
#include <Eigen\Core>
#include "MathHelpers.h"

using namespace Eigen;

/*
 * <summary> Mixture interface. </summary>
 */
template <typename T>
struct mixture {

	size_t nMax;
	size_t dimension;
	unsigned int idCounter;							// Id Counter
	unsigned int trackCounter;

	std::vector<T> components;		// Vector with components

	/* Operator overloading */
	T operator [] (size_t i) const { return components[i]; }
	T& operator[] (size_t i) { return components[i]; }

	virtual void merge(const double& _mergeThreshold) = 0;
	virtual void prune(const double& _pruneThreshold) = 0;
};

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

};

std::ostream& operator << (std::ostream& _os, const gaussian_component& _gc);

/**
 * <summary> Gaussian Mixture </summary>
 */
struct gaussian_mixture {

	size_t nMax;									// Maximum number of components
	size_t dimension;								// Dimensionality of the components
	unsigned int idCounter;							// Id Counter
	unsigned int trackCounter;
	std::vector<gaussian_component> components;		// Vector with components
	 
	/* Constructors */
	gaussian_mixture();
	gaussian_mixture(const size_t& _dim, const size_t& _nMax);
	gaussian_mixture(const size_t& _dim, const size_t& _n, const size_t& _nMax, 
		const VectorXd& _lBound, const VectorXd& _uBound, const MatrixXd& _iCov, const double& _iWeight);
	gaussian_mixture(const gaussian_mixture& _gm);

	/* Operator overloading */
	gaussian_component operator [] (size_t i) const { return components[i]; }
	gaussian_component& operator[] (size_t i)  { return components[i]; }

	/* GMM-related functions */
	void merge(const double& _mergeThreshold);
	void prune(const double& _pruneThreshold);
	void addComponent(const VectorXd& _m, const MatrixXd&, const double& _w);
	void addComponent(const gaussian_component& _gc);
	auto getEstimates(const double& _estimateThreshold) -> decltype(components);

	std::vector<double> getWeights();

	/* Miscellaneous */
	size_t size() const;
	size_t dim() const;

	/* Static functions (temporary?) */
	static VectorXd randVec(const VectorXd & _lowerBound, const VectorXd & upperBound);
	static double hellinger(const VectorXd & v1, const MatrixXd & p1, const VectorXd & v2, const MatrixXd & p2);
};

std::ostream& operator << (std::ostream& _os, const gaussian_mixture& _gm);

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

};


struct beta_gaussian_mixture {

};