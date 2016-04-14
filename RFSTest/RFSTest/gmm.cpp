#include <algorithm>
#include "gmm.h"

// Temp
#include <iostream>

/**
 * <summary> An empty constructor fo the Gaussian Mixture. </summary>
 */
gaussian_mixture::gaussian_mixture() : nMax(0), dimension(0), idCounter(0), trackCounter(0){}

/**
 * <summary> Constructor that initializes an empty Gaussian Mixture of speified dimension. </summary>
 * <param name = "_dim"> Dimensionality of the stored Gaussian components. </param>
 * <param name = "_nMax"> Maximum number of the Gaussian componentes. </param>
 */
gaussian_mixture::gaussian_mixture(const size_t & _dim, const size_t & _nMax)
	: dimension (_dim), nMax(_nMax), idCounter(1), trackCounter(1) {}

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
	: dimension(_dim), nMax(_nMax), idCounter(1), trackCounter(1)
{
	components.resize(_n);
	for (size_t i = 0; i < components.size(); i++) 
		components[i] = gaussian_component(randVec(_lBound, _uBound), _iCov, _iWeight, idCounter++);
}

/**
 * <summary> Copy constructor of the Gaussian Mixture. </summary>
 * <param name = "_gm"> An instance of Gaussian Mixture to copy from. </param>
 */
gaussian_mixture::gaussian_mixture(const gaussian_mixture & _gm)
	: dimension (_gm.dimension), nMax(_gm.nMax), idCounter( _gm.idCounter), components(_gm.components),
	trackCounter(_gm.trackCounter) {}

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

			double distance = hellinger(temp.back().m, temp.back().P, i->m, i->P);

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

/**
 * <summary> Pruning of the Gaussian Components with weights under the threshold. </summary>
 * <param name = _weightThreshold> Pruning weight threshold. </param>
 */
void gaussian_mixture::prune(const double & _pruneThreshold)
{
	auto pruned = std::remove_if(components.begin(), components.end(),
		[&_pruneThreshold](const gaussian_component& gc) { return gc.w < _pruneThreshold; });
	components.erase(pruned, components.end());
}

/**
 * <summary> Accessor of the components with weight greater than threshold </summary>
 * <param name = "_estimateThreshold"> Estimation weight threshold. </param>
 * <returns> Returns the Gaussian components with weight above threshold. </returns>
 */
auto gaussian_mixture::getEstimates(const double& _estimateThreshold) -> decltype(components)
{
	auto below = std::remove_if(components.begin(), components.end(),
		[&_estimateThreshold](const gaussian_component& gc) { return gc.w < _estimateThreshold; });

	for (auto i = components.begin(); i != below; i++) {
		if (i->tag[1] == 0)
			i->tag[1] = trackCounter++;
		i->tag[2]++;
	}

	return std::vector<gaussian_component>(components.begin(), below);
}

/**
 * <summary> Accessor of the weights variable. </summary>
 * <returns> And std::vector<double> containing the weights of the components. </returns>
 */
std::vector<double> gaussian_mixture::getWeights()
{
	std::vector<double> result;
	result.resize(components.size());

	for (size_t i = 0; i < components.size(); i++)
		result[i] = components[i].w;
	
	return result;
}

/**
 * <summary> Returns the number of Gaussian Components in the Gaussian Mixture. </summary>
 * <returns> The number of elements in the Gaussian Mixture. </returns>
 */
size_t gaussian_mixture::size() const
{
	return components.size();
}

/**
 * <summary> Returns the dimensionality of the Gaussian Components'. </summary>
 * <returns> The dimensionality of the components. </returns>
 */
size_t gaussian_mixture::dim() const
{
	return dimension;
}

/**
 * <summary> Adds a component. </summary>
 * <param = "_m"> Mean of the added component. </param>
 * <param = "_P"> Covariance of the added component. </param>
 * <param = "_w"> Weight of the added component. </param>
 */
void gaussian_mixture::addComponent(const VectorXd & _m, const MatrixXd & _P, const double & _w)
{
	assert(_m.size() == dimension && "Attempt to add a GMM component of different dimension ");

	if (components.size() < nMax)
		components.push_back(gaussian_component(_m, _P, _w, idCounter++));
}

/**
 * <summary> Adds a component. </summary>
 * <param = "_gc"> An instance of the component to be added. </param>
 */
void gaussian_mixture::addComponent(const gaussian_component & _gc)
{
	assert(_gc.m.size() == dimension && "Attempt to add a GMM component of different dimension ");

	if (components.size() < nMax) {

		components.push_back(_gc);

		if (!_gc.tag[0])
			components.back().tag[0] = idCounter++;
	}
}

/**
 * <summary> Random state vector for component initialization. </summary>
 * <param name = "_lowerBound"> The lower bound for the random vector. </param>
 * <param name = "_lowerBound"> The upper bound for the random vector. </param>
 * <returns> The random VectorXd with the specified boundraies. </returns>
 */
VectorXd gaussian_mixture::randVec(const VectorXd & _lowerBound, const VectorXd & upperBound)
{
	{
		size_t l = _lowerBound.size();
		VectorXd result = VectorXd::Zero(l);

		for (size_t i = 0; i < l; i++)
			result(i) = (double)(upperBound(i) - _lowerBound(i)) * (double)rand() / (double)RAND_MAX - (upperBound(i) - _lowerBound(i)) / 2;

		return result;
	}
}

/**
 * <summary> Hellinger distance. </summary>
 * <param name = "_v1"> First mean. </param>
 * <param name = "_v1"> First mean's covariance. </param>
 * <param name = "_v1"> Second mean. </param>
 * <param name = "_v1"> Second mean's covariance. </param>
 * <returns> Hellinger distance. </returns>
 */
inline double gaussian_mixture::hellinger(const VectorXd & _v1, const MatrixXd & _p1, const VectorXd & _v2, const MatrixXd & _p2)
{
	VectorXd vDiff = _v1 - _v2;
	MatrixXd pSum = _p1 + _p2;
	double epsilon = (-0.25 * vDiff.transpose() * pSum.inverse() * vDiff)(0, 0);
	return 1 - sqrt(sqrt((_p1 * _p2).determinant()) / (0.5 * pSum).determinant()) * exp(epsilon);
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
 * <summary> Initializes an empty tag component. </summary>
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
	return _os << _gc.tag[0] << " [" << _gc.tag[1] << "]" 
		<< "[" << _gc.tag[2] << "] " << _gc.w << " " << "("<< _gc.m.transpose() << ")";
}

/**
* <summary> Stream operator overloading for Gaussian Component. </summary>
* <param name = "_os"> A reference to the output stream. </param>
* <param name = "_gm"> Gaussian Mixture for the output. </param>
*/
std::ostream & operator<<(std::ostream & _os, const gaussian_mixture & _gm)
{
	_os << _gm.size() << "/" << _gm.nMax << std::endl;

	for (auto gc : _gm.components) 
		_os << gc << std::endl;
	
	return _os;
}
