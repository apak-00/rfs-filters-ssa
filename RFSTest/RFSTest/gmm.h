#pragma once
#include <vector>
#include <Eigen\Core>
#include <Eigen/LU>


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

	mixture();
	mixture(const size_t& _dim, const size_t& _nMax);
	mixture(const size_t& _dim, const size_t& _n, const size_t& _nMax);
	mixture(const mixture& _gm);

	/* Operator overloading */
	T operator [] (size_t i) const { return components[i]; }
	T& operator[] (size_t i) { return components[i]; }

	virtual void addComponent(const T& _c);
	virtual void merge(const double& _mergeThreshold) = 0;
	virtual void prune(const double& _pruneThreshold);
	virtual auto getEstimates(const double& _estimateThreshold) -> decltype(components);

	/* Miscellaneous */
	auto size() const { return components.size(); };
	auto dim() const { return dimension; };
	double weightSum() const;

	/* Static functions (temporary?) */
	static double hellinger(const T& _c1, const T & _c2);
	static VectorXd randVec(const VectorXd & _lowerBound, const VectorXd & upperBound);
};

template <typename T>
std::ostream& operator << (std::ostream& _os, const mixture<T>& _gm);

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
struct gaussian_mixture : mixture<gaussian_component> {
	 
	/* Constructors */
	gaussian_mixture();
	gaussian_mixture(const size_t& _dim, const size_t& _nMax);
	gaussian_mixture(const size_t& _dim, const size_t& _n, const size_t& _nMax, 
		const VectorXd& _lBound, const VectorXd& _uBound, const MatrixXd& _iCov, const double& _iWeight);
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
struct beta_gaussian_mixture : mixture<beta_gaussian_component> {

	beta_gaussian_mixture();
	beta_gaussian_mixture(const size_t& _dim, const size_t& _nMax);
	beta_gaussian_mixture(const size_t& _dim, const size_t& _n, const size_t& _nMax,
		const VectorXd& _lBound, const VectorXd& _uBound, const MatrixXd& _iCov, const double& _iWeight);
	beta_gaussian_mixture(const beta_gaussian_mixture& _gm);

	/* GMM-related functions */
	virtual void merge(const double& _mergeThreshold);

	static double betaHellinger(const beta_gaussian_component & _bgc1, const beta_gaussian_component & _bgc2);
};

/**
 * <summary> An empty constructor for the mixture interface. </summary>
 */
template<typename T>
inline mixture<T>::mixture() : nMax(0), dimension(0), idCounter(0), trackCounter(0) {}

/**
 * <summary> A constructor that initializes an empty beta-Gaussian Mixture. </summary>
 * <param name = "_dim"> Dimensinality of the stored components. </param>
 * <param name = "_nMax"> Maximum number of components. </param>
 */
template<typename T>
inline mixture<T>::mixture(const size_t & _dim, const size_t & _nMax) 
	: dimension(_dim), nMax(_nMax), idCounter(1), trackCounter(1) {}

/**
* <summary> A constructor that initializes an empty beta-Gaussian Mixture of the specified size. </summary>
* <param name = "_dim"> Dimensinality of the stored components. </param>
* <param name = "_n"> Initial number of components. </param>
* <param name = "_nMax"> Maximum number of components. </param>
*/
template<typename T>
inline mixture<T>::mixture(const size_t & _dim, const size_t & _n, const size_t & _nMax) 
	: dimension(_dim), nMax(_nMax), idCounter(1), trackCounter(1)
{
	components.resize(_n);
}

/**
 * <summary> A copy constructor of the Mixture interface. </summary>
 * <param name = "_m"> A mixture to copy from. </param>
 */
template<typename T>
inline mixture<T>::mixture(const mixture & _mixture) 
	: dimension(_mixture.dimension), nMax(_mixture.nMax), idCounter(_mixture.idCounter), 
	components(_mixture.components), trackCounter(_mixture.trackCounter) {}

/**
 * <summary> Adds a component to the mixture. </summary>
 * <param name = "_c"> A component to be added. </param>
 */
template<typename T>
inline void mixture<T>::addComponent(const T & _component)
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
inline void mixture<T>::prune(const double & _pruneThreshold)
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
inline auto mixture<T>::getEstimates(const double & _estimateThreshold) -> decltype(components)
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

	return sum;
}
/**
* <summary> Hellinger distance. </summary>
* <param name = "_v1"> First mean. </param>
* <param name = "_v1"> First mean's covariance. </param>
* <param name = "_v1"> Second mean. </param>
* <param name = "_v1"> Second mean's covariance. </param>
* <returns> Hellinger distance. </returns>
*/
template<typename T>
inline double mixture<T>::hellinger(const T& _c1, const T & _c2)
{
	VectorXd vDiff = _c1.m - _c2.m;
	MatrixXd pSum = _c1.P + _c2.P;
	double epsilon = (-0.25 * vDiff.transpose() * pSum.inverse() * vDiff)(0, 0);
	return 1 - sqrt(sqrt((_c1.P * _c2.P).determinant()) / (0.5 * pSum).determinant()) * exp(epsilon);
}

/**
 * <summary> [Old, Temporary] Returns a VectorXd containing random values within the specified range. </summary>
 * <param name = "_lowerBound"> </param>
 * <param name = "_lowerBound"> </param>
 * <returns> A VectorXd with random values. </returns>
 */
template<typename T>
inline VectorXd mixture<T>::randVec(const VectorXd & _lowerBound, const VectorXd & upperBound)
{
	size_t l = _lowerBound.size();
	VectorXd result = VectorXd::Zero(l);

	for (size_t i = 0; i < l; i++)
		result(i) = (double)(upperBound(i) - _lowerBound(i)) * (double)rand() / (double)RAND_MAX - (upperBound(i) - _lowerBound(i)) / 2;

	return result;
}

/**
* <summary> Stream operator overloading for Gaussian Component. </summary>
* <param name = "_os"> A reference to the output stream. </param>
* <param name = "_gm"> Gaussian Mixture for the output. </param>
*/
template<typename T>
inline std::ostream & operator<<(std::ostream & _os, const mixture<T> & _m)
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

	particle(const VectorXd& _m, const double& _w);
};

/*
 * <summary> Particle swarm for SMC </summary>
 */
template<typename T>
struct particle_swarm
{
	std::vector<T> particles;

	particle_swarm(const size_t& _n, const VectorXd& _mean = VectorXd(), const double& _weight = 0);
	particle_swarm(const particle_swarm<T>& _pc);

	particle_swarm<T> operator+ (const particle_swarm<T>& _gc) const;

};