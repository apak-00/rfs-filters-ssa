#ifndef UNSCENTED_SAMPLER
#define UNSCENTED_SAMPLER
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <cmath>

//#define DEBUG
#ifdef DEBUG
#include <iostream>
#endif
namespace icl
{
	using namespace Eigen;
	using namespace std;

	template<int dim = Eigen::Dynamic, typename scalar = double>
	struct unscented_sampler
	{
		typedef Eigen::Matrix<scalar, dim, 1> mean_type;
		typedef Eigen::Matrix<scalar, dim, dim> cov_type;

		virtual void get_points(const mean_type &, const cov_type &, const scalar &, vector<Eigen::VectorXd> &, vector<scalar> &) = 0;
	};

	template<int dim = Eigen::Dynamic, typename scalar = double>
	struct standard_unscented_sampler : unscented_sampler<dim, scalar>
	{
		using typename unscented_sampler<dim, scalar>::mean_type;
		using typename unscented_sampler<dim, scalar>::cov_type;

		void get_points(const mean_type &mean_, const cov_type &cov_, const scalar &w0_, vector<Eigen::VectorXd> &points, vector<scalar> &weights)
		{
			size_t dim_point = mean_.rows();
			size_t num_sigma = 2 * dim_point + 1;
			points.resize(num_sigma);
			weights.resize(num_sigma);
			// Fill the sigma weights
			scalar w1 = (1.0 - w0_) / (2.0 * (scalar)dim_point);
			weights[0] = w0_;
			points[0] = mean_;
			fill(weights.begin() + 1, weights.end(), w1);
			cov_type sqS = (dim_point / (1.0 - w0_) * cov_).llt().matrixL();
			for (size_t i = 0; i < dim_point; i++)
			{
				points[1 + i] = mean_ + sqS.col(i);
			}
			for (size_t i = 0; i < dim_point; i++)
			{
				points[1 + i + dim_point] = mean_ - sqS.col(i);
			}
		}
	};

	template<int dim = Eigen::Dynamic, typename scalar = double>
	struct simplex_unscented_sampler : unscented_sampler<dim, scalar>
	{
		// Spherical simplex
		using typename unscented_sampler<dim, scalar>::mean_type;
		using typename unscented_sampler<dim, scalar>::cov_type;
		typedef pair<size_t, size_t> sigma_hashmap_key;
		typedef std::map<sigma_hashmap_key, Eigen::VectorXd> sigma_hashmap;
		typedef typename sigma_hashmap::iterator sigma_hashmap_iterator;
		typedef typename sigma_hashmap::const_iterator sigma_hashmap_const_iterator;
		sigma_hashmap sigma_points_simplex;

		void get_points(const mean_type &mean_, const cov_type &cov_, const scalar &w0_, vector<Eigen::VectorXd> &points, vector<scalar> &weights)
		{
			size_t dim_point = mean_.rows();
			size_t num_sigma = dim_point + 2;
			points.resize(num_sigma);
			weights.resize(num_sigma);
			// Fill the sigma weights
			scalar w1 = (1.0 - w0_) / ((scalar)dim_point + 1);
			weights[0] = w0_;
			fill(weights.begin() + 1, weights.end(), w1);
			cov_type sqS = cov_.llt().matrixL();
			for (size_t i = 0; i < num_sigma; i++)
			{
				mean_type sigma;
				calc_sigma_point_simplex(i, dim_point, weights, sigma);
				points[i] = mean_ + sqS*sigma;
			}
		}

		template<typename D>
		void calc_sigma_point_simplex(const size_t &i, const size_t &j, const vector<scalar> &weights, Eigen::MatrixBase<D> &output)
		{
			sigma_hashmap_const_iterator elem = sigma_points_simplex.find(make_pair(i, j));
			if (elem != sigma_points_simplex.end())
			{
				output = elem->second;
				return;
			}
			output = Eigen::VectorXd::Zero(j);
			if (j == 1)
			{
				if (i == 0)
				{
					output(0) = 0.0;
				}
				else if (i == 1)
				{
					output(0) = -1.0 / sqrt(2.0 * weights[1]);
				}
				else if (i == 2)
				{
					output(0) = 1.0 / sqrt(2.0 * weights[1]);
				}
			}
			else
			{
				Eigen::VectorXd part = Eigen::VectorXd::Zero(j - 1);
				if (i == 0)
				{
					calc_sigma_point_simplex(0, j - 1, weights, part);
					output.block(0, 0, j - 1, 1) = part;
					output[j - 1] = 0.0;
				}
				else if (i <= j)
				{
					calc_sigma_point_simplex(i, j - 1, weights, part);
					output.block(0, 0, j - 1, 1) = part;
					output[j - 1] = -1.0 / sqrt(j * (j + 1) * weights[j]);
				}
				else if (i == j + 1)
				{
					output.block(0, 0, j - 1, 1) = Eigen::VectorXd::Zero(j - 1);
					output[j - 1] = j / sqrt(j * (j + 1) * weights[j]);
				}
				else
				{
					throw std::runtime_error("Bad parameters for the simplex sigma point calculation");
				}
			}
			sigma_points_simplex[make_pair(i, j)] = output;
		}
	};

	template<int dim = Eigen::Dynamic, typename scalar = double>
	struct fourth_order_unscented_sampler : unscented_sampler<dim, scalar>
	{
		using typename unscented_sampler<dim, scalar>::mean_type;
		using typename unscented_sampler<dim, scalar>::cov_type;

		std::vector<Eigen::VectorXd> points_01;
		std::vector<scalar> weights_01;

		void get_points(const mean_type &mean_, const cov_type &cov_, const scalar &w0_, vector<Eigen::VectorXd> &points, vector<scalar> &weights)
		{
			size_t dim_point = mean_.rows();
			size_t num_sigma = 2 * dim_point * dim_point + 1;
			points.resize(num_sigma);
			weights.resize(num_sigma);
			if (points_01.size() == 0)
			{
				calculate_points_01(dim_point, w0_);
			}
			// Create the sigma points and then project them using the SQRT(cov)
			cov_type sqS = cov_.llt().matrixL();
			for (size_t i = 0; i < num_sigma; i++)
			{
				points[i] = mean_ + sqS * points_01[i];
				weights[i] = weights_01[i];
			}
		}

		void calculate_points_01(const size_t &dim_point, double w2)
		{
			size_t num_sigma = 2 * dim_point * dim_point + 1;
			points_01.reserve(num_sigma);
			weights_01.reserve(num_sigma);

			/*
			const double w0 = (dim_point * dim_point - 3 * dim_point + 6) / 6;
			const double w1 = -((scalar)dim_point - 2) / 6;
			const double w2 = 1.0 / 12.0;
			const double s1 = sqrt(3.0);
			const double s2 = sqrt(3.0);
			*/
			const double n = dim_point;
			// This is required to avoid ill-defined sqrt in w1
			if (w2 > (n - 1) / 12 || w2 < 0)
			{
				w2 = 1 / 12;
			}
			const double w0 = -((24 * n*n - 24 * n) * w2*w2 + (2 * n*n*n - 12 * n*n + 14 * n - 12) * w2 + n - 1) / (12 * w2 - n + 1);
			const double w1 = ((2 * n*n - 8 * n + 8) * w2) / (12 * w2 - n + 1);
			const double s1 = -sqrt(-(12 * w2 - n + 1) / ((n - 2) * w2)) / 2;
			const double s2 = 1 / (2 * sqrt(w2));

			points_01.push_back(Eigen::VectorXd::Zero(dim_point));
			weights_01.push_back(w0);
			for (size_t i = 0; i < dim_point; i++)
			{
				Eigen::VectorXd point = Eigen::VectorXd::Zero(dim_point);
				point[i] = s1;
				points_01.push_back(point);
				points_01.push_back(-point);
				weights_01.push_back(w1);
				weights_01.push_back(w1);
			}
			for (size_t i = 0; i < dim_point; i++)
			{
				for (size_t j = i + 1; j < dim_point; j++)
				{
					for (int k = 1; k > -2; k -= 2)
					{
						for (int l = 1; l > -2; l -= 2)
						{
							Eigen::VectorXd point = Eigen::VectorXd::Zero(dim_point);
							point[i] = k*s2;
							point[j] = l*s2;
							points_01.push_back(point);
							weights_01.push_back(w2);
						}
					}
				}
			}
		}
	};
}
#endif