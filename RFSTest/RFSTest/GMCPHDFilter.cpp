// TODO: Complete

#include "GMCPHDFilter.h"
#include "MathHelpers.h"

GMCPHDFilter::GMCPHDFilter(std::shared_ptr<KalmanFilter> _kf, const unsigned int & _nBirthComponents, const double & _birthIntensity, 
	const double & _pS, const MatrixXd & _iCov, const VectorXd & _lBound, const VectorXd & _uBound, const size_t& _NMax) : 
	nBirthComponents(_nBirthComponents), birthIntensity(_birthIntensity), pS(_pS), initialCovariance(_iCov), lowerBound(_lBound), upperBound(_uBound),
	NMax(_NMax)
{
	filter = _kf;

	// CPHD
	cardinality = VectorXd::Zero(_NMax);
}

void GMCPHDFilter::predict(gaussian_mixture & _gmm, Sensor & _sensor)
{
	double range;
	std::vector<double> birthRanges;
	VectorXd birth(_gmm.dim());

	// Prediction for all of the components
	for (auto &gc : _gmm.components) {
		filter->predict(gc);
		gc.w *= pS;
	}

	// New target birth
	double initialWeight = birthIntensity / nBirthComponents;

	if (_sensor.getZ().size() != 0)
	{
		for (size_t i = 0; i < nBirthComponents; i++)
			birthRanges.push_back((double)(i + 2) * 200);

		if (nBirthComponents == 1)
		{
			birthRanges[0] = 1000; // + rand() % 500 - 250;
		}

		for (size_t i = 0; (i < nBirthComponents) && (_gmm.size() < _gmm.nMax); i++)
		{
			// Uniform birth test
			range = birthRanges[i];

			if (_gmm.dim() == 2)
				birth << range, 0;
			else if (_gmm.dim() == 6)
			{
				VectorXd m(_gmm.dim());
				m << range, _sensor.getBearing(), 0, 0, 0;
				birth = Astro::razelToTEME(m, _sensor.getPosition(), _sensor.getDateJD(), _sensor.getLOD(), _sensor.getXp(), _sensor.getYp());
			}
			_gmm.addComponent(gaussian_component(birth, initialCovariance, initialWeight, _gmm.idCounter++));
		}

		// Cardinality prediction
		// Survivning components
		// TODO: Optimize

		VectorXd survivedCardinalityPredicted = VectorXd::Zero(NMax), cardinalityPredicted = VectorXd::Zero(NMax + 1), terms;

		size_t idxj, idxl, idxn;
		
		// Log sums ----------
		// TODO: Optimize 
		VectorXd tempLogSum = VectorXd::LinSpaced(NMax, 1.0, (double) NMax);		// 1, 2, ... NMax		
		for (size_t t = 0; t < NMax; t++)
			tempLogSum(t) = log(tempLogSum(t));

		// Integral sum
		// ln(1) = 0
		for (size_t t = 1; t < NMax; t++)
			tempLogSum(t) = tempLogSum(t) + tempLogSum(t - 1);
		// Log sum end ----------

		for (size_t j = 0; j < NMax; j++)
		{
			idxj = j + 1;
			terms = VectorXd::Zero(NMax + 1);

			for (size_t l = 0; l < NMax; l++)
			{
				idxl = l + 1;
				terms(idxl) = exp(tempLogSum(l) - tempLogSum(j) - tempLogSum(l-j) + j * log(pS) + (l - j) * log(1 - pS)) * cardinality(idxl);
			}

			cardinality(idxj) = terms.sum();
		}

		for (size_t n = 0; n < NMax; n++)
		{
			idxn = n + 1;
			terms = VectorXd::Zero(NMax + 1);
			for (size_t j = 0; j < n; j++)
			{
				idxj = j + 1;
				terms(idxj) = exp(-birthIntensity + (n-j) * log(birthIntensity) - tempLogSum(n-j)) * survivedCardinalityPredicted(idxj);
			}
			cardinalityPredicted(idxn) = terms.sum();
		}

		// Predicted cardinality normalization
		cardinalityPredicted = cardinalityPredicted / cardinalityPredicted.sum();
	}
}

void GMCPHDFilter::update(gaussian_mixture & _gmm, Sensor & _sensor)
{
	double cz = 1.0 / 57903 * 200;

	auto pD = _sensor.getPD();
	size_t n0 = _gmm.size();
	size_t zSize = _sensor.getZ().size();

	MatrixXd qz = MatrixXd::Zero(zSize, _gmm.size());
	VectorXd predWeights = _gmm.getWeightsVector();

	for (size_t i = 0; i < zSize; i++)
	{

		for (size_t j = 0; j < n0; j++) 
		{
			gaussian_component gct(_gmm[j]);
			filter->update(gct, _sensor, i);

			auto q = (1.0 / sqrt(pow(2.0 * M_PI, _sensor.getZDim()) * _sensor.getS().determinant()))
				* exp(-0.5 * MathHelpers::mahalanobis(_sensor.getZ(i), _sensor.getPredictedZ(), _sensor.getS()));

			qz(i, j) = q;
		}
	}

	// Pre-calculation for the elementary symmetric functions
	VectorXd xiVals = VectorXd::Zero(zSize);
	double tempProd;

	for (size_t l = 0; l < zSize; l++)
	{
		tempProd = (predWeights.transpose() * qz.row(l))(0);
		xiVals(l) = pD * tempProd / cz;
	}
		
	VectorXd esfValsE = MathHelpers::esf(xiVals), esfValsD = MatrixXd::Zero(zSize, zSize);

	for (size_t l = 0; l < zSize; l++)
	{
		//esfValsD.row(l) = MathHelpers::esf();
	}

	
}
