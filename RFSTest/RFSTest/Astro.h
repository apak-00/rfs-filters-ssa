#pragma once
#include <Eigen/Eigen>

using namespace Eigen;

namespace Astro {

	/** Temporary structure for storing the date and time with millisecond precision **/
	struct date {

		int year, month, day, hour, min;
		double sec;

		date() : year(0), month(0), day(0), hour(0), min(0), sec(0) {}
		date(const int& _yy, const int& _mm, const int& _dd, const int& _hr, const int& _nMin, const double& _nSec) :
			year(_yy), month(_mm), day(_dd), hour(_hr), min(_nMin), sec(_nSec) {}
		date(const date& _d) : year(_d.year), month(_d.month), day(_d.day), hour(_d.hour), min(_d.min), sec(_d.sec) {}

		// TODO: Fix
		double toMS() const {
			return sec + min * 60.0 + hour * 3600.0;
		}

		double operator- (const date& _d) const {
			return toMS() - _d.toMS();
		}
	};

	const double R_EQ = 6378.137;				// Major semiaxis (equatorial radius)
	const double E = 8.1819190842622e-2;		// First eccentricity
	const double E2 = 6.69437999014e-3;			// First eccentricity squared
	const double T0 = 2451545.0;				// J2000 reference epoch (days)

	const double MU_E = 3.986004418e5;				// Earth Gravitational Constant (km^3 / s^2)
	const double MU_C = 1.327124400e11;				// Sun Gravitational Constant (km^3 / s^2)

	/* Main Transformations */
	VectorXd razelToSEZ(const VectorXd& _razel);
	VectorXd sezToRAZEL(const VectorXd& _sez);

	VectorXd razelToSEZVallado(const VectorXd& _razel);
	VectorXd sezToRAZELVallado(const VectorXd& _sez);

	VectorXd geodeticToECEF(const VectorXd& _geo);

	VectorXd sezToECEF(const VectorXd& _sez, const VectorXd& _geo);
	VectorXd ecefToSEZ(const VectorXd& _ecef, const VectorXd& _geo);

	VectorXd ecefToTEME(const VectorXd& _ecef, const double& _jd, const double& _lod, const double& _xp, const double& _yp);
	VectorXd temeToECEF(const VectorXd& _teme, const double& _jd, const double& _lod, const double& _xp, const double& _yp);

	/* Other Transformations*/
	MatrixXd getSEZToRAZELJacobian(const VectorXd& _sez, const size_t& _zDim);
	MatrixXd getSEZToRAZELJacobianFADBAD(const VectorXd& _sez, const size_t& _zDim);
	MatrixXd getSEZToTEMECovTfMat(const VectorXd& _geo, const double& _jd, const double& _xp, const double& _yp, const size_t& _dim);

	/* Derived */
	VectorXd sezToTEME(const VectorXd& _sez, const VectorXd& _geo, const double& _jd, const double& _lod, const double& _xp, const double& _yp);
	VectorXd temeToSEZ(const VectorXd& _teme, const VectorXd& _geo, const double& _jd, const double& _lod, const double& _xp, const double& _yp);

	VectorXd razelToTEME(const VectorXd& _razel, const VectorXd& _geo, const double& _jd, const double& _lod, const double& _xp, const double& _yp);
	VectorXd temeToRAZEL(const VectorXd& _teme, const VectorXd& _geo, const double& _jd, const double& _lod, const double& _xp, const double& _yp);

	/* Auxilary */
	MatrixXd getPolarMotionMatrix(const double& _xp, const double& _yp);

	double getJ2000Ref(const double _jdn);
	double getGMST(const double& _jd);
	double getJulianDay(const date& _date);
	double getJulianDay(const int& _year, const int& _month, const int& _day, const int& _hour, const int& _minute, const double& _sec);

	Matrix3d rotX(const double& _angle);
	Matrix3d rotY(const double& _angle);
	Matrix3d rotZ(const double& _angle);

	/* Other */
	VectorXd integrationPrediction(const VectorXd& _state, const double& _dt);
	MatrixXd getShepperdMatrix(const VectorXd& _m0, const double& _t, VectorXd& _m, const double& _mu);
}