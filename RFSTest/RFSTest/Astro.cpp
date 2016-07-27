#include "math.h"
#include "Astro.h"
#include "fadiff.h"
#include <boost\numeric\odeint\integrate\integrate.hpp>

#include <cmath>
#include <iostream>
#include <iomanip>

typedef std::vector<double> state_type;

/* SEZ to RAZEL transformation functions for automatic differentiation. */

/**
 * <summary> Cartesian -> spherical range. </summary>
 * <returns> Range in spherical coordinate system. </returns>
 */
fadbad::F<double> xSEZtoRAZEL(const fadbad::F<double> _x, const fadbad::F<double> _y, const fadbad::F<double> _z,
	const fadbad::F<double> _xdot, const fadbad::F<double> _ydot, const fadbad::F<double> _zdot) 
{
	return fadbad::sqrt(_x * _x + _y* _y + _z * _z);
}

/**
* <summary> Cartesian -> spherical azimuth. </summary>
* <returns> Azimuth in spherical coordinate system. </returns>
*/
fadbad::F<double> ySEZtoRAZEL(const fadbad::F<double> _x, const fadbad::F<double> _y, const fadbad::F<double> _z,
	const fadbad::F<double> _xdot, const fadbad::F<double> _ydot, const fadbad::F<double> _zdot) 
{
	return fadbad::atan2(_y, _x);
}

/**
* <summary> Cartesian -> spherical elevation. </summary>
* <returns> Elevation in spherical coordinate system. </returns>
*/
fadbad::F<double> zSEZtoRAZEL(const fadbad::F<double> _x, const fadbad::F<double> _y, const fadbad::F<double> _z,
	const fadbad::F<double> _xdot, const fadbad::F<double> _ydot, const fadbad::F<double> _zdot) 
{
	return fadbad::atan2(_z, sqrt(_x * _x + _y * _y));
}

/**
* <summary> Cartesian -> to spherical range rate. </summary>
* <returns> Range rate in spherical coordinate system. </returns>
*/
fadbad::F<double> xdotSEZtoRAZEL(const fadbad::F<double> _x, const fadbad::F<double> _y, const fadbad::F<double> _z,
	const fadbad::F<double> _xdot, const fadbad::F<double> _ydot, const fadbad::F<double> _zdot) 
{
	return (_x * _xdot + _y * _ydot + _z * _zdot) /
		fadbad::sqrt(_x * _x + _y* _y + _z * _z);
}

/**
* <summary> Cartesian -> to spherical azimuth rate. </summary>
* <returns> Azimuth rate in spherical coordinate system. </returns>
*/
fadbad::F<double> ydotSEZtoRAZEL(const fadbad::F<double> _x, const fadbad::F<double> _y, const fadbad::F<double> _z,
	const fadbad::F<double> _xdot, const fadbad::F<double> _ydot, const fadbad::F<double> _zdot) 
{
	return (_x * _ydot - _y * _xdot) / (_x * _x + _y * _y);
}

/**
* <summary> Cartesian -> spherical elevation rate. </summary>
* <returns> Elevation rate in spherical coordinate system. </returns>
*/
fadbad::F<double> zdotSEZtoRAZEL(const fadbad::F<double> _x, const fadbad::F<double> _y, const fadbad::F<double> _z,
	const fadbad::F<double> _xdot, const fadbad::F<double> _ydot, const fadbad::F<double> _zdot) 
{
	return  (_zdot * (_x * _x + _y * _y) - _z * (_x * _xdot + _y * _ydot)) / ((_x * _x + _y * _y + _z * _z) * (fadbad::sqrt(_x * _x + _y * _y)));
}

/**
* <summary> Function defining the derivative of the 6-dimensional orbital state vector with respect to time. </summary>
* <par> Simplest gravitational model is assumed with no gravitational acceleration. </par>
* <par> Function is defined according to Boost's odeint library conventions. State type is defined to be 
* std::vector<double> of dimesion 6. </par>
*
* <param name = "x"> A constant reference to initial state vector </param>
* <param name = "dxdt"> A reference to the variable that should contain the derivative of the initial state vector with respect to time </param>
*/
void integrateOrbit(const state_type &x, state_type &dxdt, const double /* t */)
{
	double r3 = pow(sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]), 3.0);

	dxdt.resize(6);

	dxdt[0] = x[3];
	dxdt[1] = x[4];
	dxdt[2] = x[5];
	dxdt[3] = -Astro::MU_E * x[0] / r3;
	dxdt[4] = -Astro::MU_E * x[1] / r3;
	dxdt[5] = -Astro::MU_E * x[2] / r3;
}

/*
 * <summary> Astro namespace. Contains basic operations for coordinate transformations. </summary>
 */
namespace Astro 
{
	/*
	 * <summary> Output stream operator overload for the date struct. </summary>
	 * <par> Output format depends on the specified enum. </par>
	 * <param name = "_os"> A reference to the output stream. </param>
	 */
	std::ostream & operator<<(std::ostream & _os, const date & _d)
	{
		switch (_d.outputType) {
		case 0: _os << _d.year << "," << _d.month << "," << _d.day << "," << _d.hour << "," << _d.minute << "," << _d.sec;
			break;
		}
		return _os;
	}

	/**
	* <summary> Range, Azimuth, Elevation to Cartesian South-East-Zenith (SEZ) system. </summary>
	* <par> Similar to MATLAB sph2cart. </par>
	* <param name = "_razel"> A VectorXd containing the position information in Range-Azimuth-Elevation frame.
	* Physical azimuth measured from positive X-axis (S) is assumed). </param>
	* <returns> VectorXd in the South-East-Zenith frame. </returns>
	*/
	VectorXd razelToSEZ(const VectorXd& _razel) 
	{
		VectorXd sez(_razel.size());

		double sinAz = sin(_razel(1)), cosAz = cos(_razel(1)),
			sinEl = sin(_razel(2)), cosEl = cos(_razel(2));

			// Slightly different from Vallado's book
			sez(0) = _razel(0) * cosEl * cosAz;
			sez(1) = _razel(0) * cosEl * sinAz;
			sez(2) = _razel(0) * sinEl;

			if (_razel.size() == 6)
			{
				sez(3) = _razel(3) * cosEl * cosAz - sez(2) * cosAz * _razel(5) - sez(1) * _razel(4);
				sez(4) = _razel(3) * cosEl * sinAz - sez(2) * sinAz * _razel(5) + sez(0) * _razel(4);
				sez(5) = _razel(3) * sinEl + _razel(0) * _razel(5) * cosEl;
			}

		return sez;
	}

	/**
	* <summary> Cartesian South-East-Zenith to Range, Azimuth, Elevation. </summary>
	* <par> Similar to MATLAB cart2sph. </par>
	* <param name = "_sez"> VectorXd containing the position in South-East-Zenith frame. </param>
	* <returns> VectorXd with the position in the Range-Azimuth-Elevation frame. </returns>
	*/
	VectorXd sezToRAZEL(const VectorXd& _sez) 
	{
		VectorXd razel(_sez.size());

		double s2 = _sez(0) * _sez(0) + _sez(1) * _sez(1);
		double r2 = s2 + _sez(2) * _sez(2);
		double s = sqrt(s2);

		razel(0) = sqrt(r2);
		razel(1) = atan2(_sez(1), _sez(0));
		razel(2) = atan2(_sez(2), s);

		if (_sez.size() == 6)
		{
			razel(3) = (_sez(0) * _sez(3) + _sez(1) * _sez(4) + _sez(2) * _sez(5)) / razel(0);
			razel(4) = (_sez(0) * _sez(4) - _sez(1) * _sez(3)) / s2;
			razel(5) = (_sez(5) * s2 - _sez(2) * (_sez(0) * _sez(3) + _sez(1) * _sez(4))) / (r2 * s);
		}

		return razel;
	}

	/**
	* <summary> Range, Azimuth, Elevation to Cartesian South-East-Zenith (SEZ) system. D. A. Vallado's book. </summary>
	* <par> Similar to MATLAB sph2cart. </par>
	* <param name = "_razel"> A VectorXd containing the position information in Range-Azimuth-Elevation frame.
	* Physical azimuth measured from positive X-axis (S) is assumed). </param>
	* <returns> VectorXd in the South-East-Zenith frame. </returns>
	*/
	VectorXd razelToSEZVallado(const VectorXd & _razel)
	{
		VectorXd sez(_razel.size());

		double sinAz = sin(_razel(1)), cosAz = cos(_razel(1)),
			sinEl = sin(_razel(2)), cosEl = cos(_razel(2));

		sez(0) = -_razel(0) * cosEl * cosAz;		// Minus here
		sez(1) = _razel(0) * cosEl * sinAz;
		sez(2) = _razel(0) * sinEl;

		if (_razel.size() == 6)
		{
			sez(3) = -_razel(3) * cosEl * cosAz + sez(2) * cosAz *  _razel(5) + sez(1) * _razel(4);		// Minus here
			sez(4) = _razel(3) * cosEl * sinAz - sez(2) * sinAz *  _razel(5) + sez(0) * _razel(4);
			sez(5) = _razel(3) * sinEl + _razel(0) *  _razel(5) * cosEl;
		}

		return sez;
	}

	/**
	* <summary> Cartesian South-East-Zenith to Range, Azimuth, Elevation. D. A. Vallado's book. </summary>
	* <par> Similar to MATLAB cart2sph. </par>
	* <param name = "_sez"> VectorXd containing the position in South-East-Zenith frame. </param>
	* <returns> VectorXd with the position in the Range-Azimuth-Elevation frame. </returns>
	*/
	VectorXd sezToRAZELVallado(const VectorXd & _sez)
	{
		VectorXd razel(_sez.size());

		double s2 = _sez(0) * _sez(0) + _sez(1) * _sez(1);
		double r2 = s2 + _sez(2) * _sez(2);
		double s = sqrt(s2);

		razel(0) = sqrt(r2);
		razel(1) = atan2(-_sez(1) / s, _sez(0) / s);
		razel(2) = asin(_sez(2) / razel(0));		// Asin here

		if (_sez.size() == 6)
		{
			razel(3) = (_sez(0) * _sez(3) + _sez(1) * _sez(4) + _sez(2) * _sez(5)) / razel(0);
			razel(4) = (_sez(1) * _sez(3) - _sez(0) * _sez(4)) / s2;
			razel(5) = (_sez(5) - razel(3) * sin(razel(2))) / s;
		}

		return razel;
	}

	/**
	* <summary> Calculates the Jacobian of the transformation between the South-East-Zenith and
	* Range-Azimuth-Elevation frames. </summary>
	* <param name = "_sez"> VectorXd in the South-East-Zenith frame. </param>
	* <param name = "_zDim"> The dimensionality of the observation. </param>
	* <returns> MatrixXd containing the partial derivatives of transformation. </returns>
	*/
	MatrixXd getSEZToRAZELJacobian(const VectorXd& _sez, const size_t& _zDim) 
	{
		MatrixXd j = MatrixXd::Zero(_zDim, _sez.size());

		double r2 = _sez(0) * _sez(0) + _sez(1) * _sez(1) + _sez(2) * _sez(2);
		double r = sqrt(r2);
		double rxy2 = _sez(0) * _sez(0) + _sez(1) * _sez(1);
		double rxy = sqrt(rxy2);

		// Range
		j(0, 0) = _sez(0) / r;
		j(0, 1) = _sez(1) / r;
		j(0, 2) = _sez(2) / r;

		// Azimuth
		j(1, 0) = -_sez(1) / rxy2;		// - y / (x^2 + y^2)
		j(1, 1) = _sez(0) / rxy2;		//	 x / (x^2 + y^2)

		// Elevation
		j(2, 0) = -_sez(2) * _sez(0) / (r2 * rxy);
		j(2, 1) = -_sez(2) * _sez(1) / (r2 * rxy);
		j(2, 2) = rxy / r2;

		// TODO: 6-dim
		if (_zDim == 6) 
		{
			// Range rate partial derivatives
			double pdotv = _sez(0) * _sez(3) + _sez(1) * _sez(4) + _sez(2) * _sez(5);

			j(3, 0) = (_sez(3) * r - pdotv * _sez(0) / r) / r2;
			j(3, 1) = (_sez(4) * r - pdotv * _sez(1) / r) / r2;
			j(3, 2) = (_sez(5) * r - pdotv * _sez(2) / r) / r2;
			j(3, 3) = _sez(0) / r;
			j(3, 4) = _sez(1) / r;
			j(3, 5) = _sez(2) / r;

			// Azimuth rate partial derivatives
			double rxy4 = rxy2 * rxy2;

			j(4, 0) = (_sez(4) * rxy2 - 2 * _sez(0) * (_sez(0) * _sez(4) - _sez(1) * _sez(3))) / rxy4;
			j(4, 1) = (-_sez(3) * rxy2 - 2 * _sez(1) * (_sez(0) * _sez(4) - _sez(1) * _sez(3))) / rxy4;
			//j(4, 2) = 0;
			j(4, 3) = -_sez(1) / rxy2;
			j(4, 4) = _sez(0) / rxy2;
			//j(4, 5) = 0;

			// Elevation rate partial derivatives
			double n = _sez(5) * rxy2 - _sez(2) * (_sez(0) * _sez(3) + _sez(1) * _sez(4)),
				d = r2 * rxy;
			double d2 = d * d;

			j(5, 0) = ((2 * _sez(0) * _sez(5) - _sez(2) * _sez(3)) * d - n * (2 * _sez(0) * rxy + _sez(0) / rxy)) / d2;
			j(5, 1) = ((2 * _sez(1) * _sez(5) - _sez(2) * _sez(4)) * d - n * (2 * _sez(1) * rxy + _sez(1) / rxy)) / d2;;
			j(5, 2) = ((_sez(0) * _sez(3) + _sez(1) * _sez(4)) * d - n * 2 * _sez(2) * rxy) / d2;
			j(5, 3) = -_sez(2) * _sez(0) / d;
			j(5, 4) = -_sez(2) * _sez(1) / d;
			j(5, 5) = rxy2 / d;
		}

		return j;
	}

	/**
	* <summary> Calculates the covariance transformation between SEZ and TEME Cartesian frames.
	* Linear transformations are not taken into account. </summary>
	* <param name = geo> VectorXd containing geodetic coordinates (latitude, longitude (degrees)
	* and altitude (km)) </param>
	* <param name = _jd> A date_ms structure contatining the date of the transofrmation. </param>
	* <param name = _xp> Polar motion coefficient for x axis (radians). </param>
	* <param name = _yp> Polar motion coefficient for x axis (radians). </param>
	* <param name = _dim> The dimensionality of the required transformation. </param>
	* <returns> Transformation MatrixXd </returns>
	*/
	MatrixXd getSEZToTEMECovTfMat(const VectorXd & _geo, const double & _jd, const double & _xp, const double & _yp, const size_t & _dim)
	{
		double gmst = getGMST(_jd);

		MatrixXd rotGd3 = rotZ(_geo(1)) * rotY(M_PI / 2.0 - _geo(0));
		MatrixXd rotPM3 = getPolarMotionMatrix(_xp, _yp);
		MatrixXd rotGMST3 = rotZ(gmst);

		MatrixXd tf = rotGMST3 * rotPM3 * rotGd3;

		if (_dim == 3)
			return tf;
		else if (_dim == 6) 
		{
			MatrixXd tf6 = MatrixXd::Zero(_dim, _dim);
			tf6.block<3, 3>(0, 0) = tf;
			tf6.block<3, 3>(3, 3) = tf;
			return tf6;
		}

		return MatrixXd();
	}

	/**
	 * <summary> Geodetic latitude, longitude, altitude to Earth-Centered, Earth-Fixed frame. </summary>
	 * <param name = "_geo"> Position in geodetic latitude, longitude, altitude. </param>
	 * <returns> VectorXd in Earth-Centered, Earth-Fixed frame. </returns>
	 */
	VectorXd geodeticToECEF(const VectorXd& _geo) 
	{
		VectorXd ecef(6);

		double s = sin(_geo(0));
		double N = R_EQ / sqrt(1.0 - E2 * pow(s, 2));			// Radius of curvature in prime meridian
		double t = (N + _geo(2)) * cos(_geo(0));

		ecef << t * cos(_geo(1)), t * sin(_geo(1)), ((1 - E2) * N + _geo(2)) * s, 0, 0, 0;     // Zero velocity

		return ecef;
	}

	/**
	 * <summary> South-East-Zenith Cartesian to Earth-Centered, Earth-Fixed Cartesian. </summary>
	 * <param name = "_sez"> VectorXd in South-East-Zenith frame. <param>
	 * <param name = "_geo"> Position of the site in geodetic latitude, longitude, altitude. </param>
	 * <returns> VectorXd in Earth-Centered, Earth-Fixed frame. </returns>
	 */
	VectorXd sezToECEF(const VectorXd & _sez, const VectorXd & _geo)
	{
		VectorXd site = geodeticToECEF(_geo), ecef(_sez.size());
		MatrixXd rot = rotZ(_geo(1)) * rotY(M_PI / 2.0 - _geo(0));

		if (_sez.size() == 3) 
		{
			ecef << rot * _sez;
			ecef += site.head(3);
		}
		else if (_sez.size() == 6) 
		{
			ecef << rot * _sez.head(3), rot * _sez.tail(3);
			ecef += site;
		}

		return ecef;
	}

	/**
	 * <summary> Eart-Centered, Earth-Fixed to South-East-Zenith Cartesian Frame. </summary>
	 * <param name = "_ecef"> VectorXd in Earth-Centered, Earth-Fixed frame. </param>
	 * <param name = "_geo"> Position of the site in geodetic latitude, longitude, altitude. </param>
	 * <returns> VectorXd in South, East, Zenith frame. </returns>
	 */
	VectorXd ecefToSEZ(const VectorXd & _ecef, const VectorXd & _geo)
	{
		MatrixXd rot = rotY(-(M_PI / 2.0 - _geo(0))) * rotZ(-_geo(1));
		VectorXd sez;

		if (_ecef.size() == 3) 
		{
			sez = _ecef - geodeticToECEF(_geo).head(3);
			sez << rot * sez;
		}
		else if (_ecef.size() == 6) 
		{
			sez = _ecef - geodeticToECEF(_geo);
			sez << rot * sez.head(3), rot * sez.tail(3);
		}

		return sez;
	}

	/**
	* <summary> Converts ECEF (WGS 84) coordinates to Earth-Centered Inertial True Equator, Mean Equinox. </summary>
	* <param name = "_ecef"> VectorXd in Earth-Centered, Earth-Fixed frame. </param>
	* <param name = "_jd"> Julian Day number. </param>
	* <param name = "_lod"> Excessive length of day (LOD). </param>
	* <param name = "_xp"> X polar motion coefficient (rad). </param>
	* <param name = "_yp"> Y polar motion coefficient (rad). </param>
	* <returns> A VectorXd with Earth-Centered Inertial (X axis towards mean equinox) Cartesian coordinates. </returns>
	*/
	VectorXd ecefToTEME(const VectorXd & _ecef, const double & _jd, const double & _lod, const double & _xp, const double & _yp)
	{
		double gmst = getGMST(_jd), thetaSa = 7.29211514670698e-05 * (1.0 - _lod / 86400.0);

		VectorXd teme(_ecef.size());
		MatrixXd pm = getPolarMotionMatrix(_xp, _yp), rotGMST = rotZ(gmst);

		Vector3d omegaEarth,
			posPEF = pm * _ecef.head(3); // Position in Pseudo Earth-Fixed Frame
		omegaEarth << 0, 0, thetaSa;

		if (_ecef.size() == 3)
			teme << rotGMST * posPEF;
		else if (_ecef.size() == 6)
			teme << rotGMST * posPEF, rotGMST * (pm * _ecef.tail(3) + omegaEarth.cross(posPEF));

		return teme;
	}

	/**
	* <summary> Converts TEME (True Equator, Mean Equinox) coordinates to ECEF (WGS 84). </summary>
	* <param name = "_teme"> VectorXd in True Equator, Mean Equinox frame. </param>
	* <param name = "_jd"> Julian Day number. </param>
	* <param name = "_lod"> Excessive length of day (LOD). </param>
	* <param name = "_xp"> X polar motion coefficient (rad). </param>
	* <param name = "_yp"> Y polar motion coefficient (rad). </param>
	* <returns> A VectorXd with ECEF Cartesian coordinates. </returns>
	*/
	VectorXd temeToECEF(const VectorXd & _teme, const double & _jd, const double & _lod, const double & _xp, const double & _yp)
	{
		double gmst = getGMST(_jd);

		VectorXd ecef(_teme.size());
		Vector3d tempThetaSa;
		MatrixXd pm = getPolarMotionMatrix(_xp, _yp).inverse();

		MatrixXd rotGMST(3, 3);
		rotGMST = rotZ(-gmst);		// st from Vallado's code (Matlab) // WRONG???

		tempThetaSa << 0, 0, 7.29211514670698e-05 * (1.0 - _lod / 86400.0);
		 
		if (_teme.size() == 3)
			ecef = pm * rotGMST * _teme.head(3);

		else if (_teme.size() == 6) 
		{
			Vector3d temp = rotGMST * _teme.head(3);
			ecef << pm * temp, pm * (rotGMST * _teme.tail(3) - tempThetaSa.cross(temp));
		}

		return ecef;
	}

	/**
	 * <summary> South, East, Zenith to True Equator, Mean Equinox. </summary>
	 * <param name = "_sez"> VectorXd in True Equator, Mean Equinox frame. </param>
	 * <param name = "_geo"> Position of the site in geodetic latitude, longitude, altitude. </param>
	 * <param name = "_jd"> Julian Day number. </param>
	 * <param name = "_lod"> Excessive length of day (LOD). </param>
	 * <param name = "_xp"> X polar motion coefficient (rad). </param>
	 * <param name = "_yp"> Y polar motion coefficient (rad). </param>
	 * <returns> VectorXd in True Equator, Mean Equinox frame. </returns>
	 */
	VectorXd sezToTEME(const VectorXd & _sez, const VectorXd& _geo, const double & _jd, const double & _lod, const double & _xp, const double & _yp)
	{
		return ecefToTEME(sezToECEF(_sez, _geo), _jd, _lod, _xp, _yp);
	}

	/**
	* <summary> True Equator, Mean Equinox to South, East, Zenith. </summary>
	* <param name = "_teme"> VectorXd in South, East, Zenith frame. </param>
	* <param name = "_geo"> Position of the site in geodetic latitude, longitude, altitude. </param>
	* <param name = "_jd"> Julian Day number. </param>
	* <param name = "_lod"> Excessive length of day (LOD). </param>
	* <param name = "_xp"> X polar motion coefficient (rad). </param>
	* <param name = "_yp"> Y polar motion coefficient (rad). </param>
	* <returns> VectorXd in South, East, Zenith frame. </returns>
	*/
	VectorXd temeToSEZ(const VectorXd & _teme, const VectorXd & _geo, const double & _jd, const double & _lod, const double & _xp, const double & _yp)
	{
		return ecefToSEZ(temeToECEF(_teme, _jd, _lod, _xp, _yp), _geo);
	}

	/**
	* <summary> Range, Azimuth, Elevation to True Equator, Mean Equinox. </summary>
	* <param name = "_razel"> VectorXd in Range, Azimuth, Elevation frame. </param>
	* <param name = "_geo"> Position of the site in geodetic latitude, longitude, altitude. </param>
	* <param name = "_jd"> Julian Day number. </param>
	* <param name = "_lod"> Excessive length of day (LOD). </param>
	* <param name = "_xp"> X polar motion coefficient (rad). </param>
	* <param name = "_yp"> Y polar motion coefficient (rad). </param>
	* <returns> VectorXd in True Equator, Mean Equinox frame. </returns>
	*/
	VectorXd razelToTEME(const VectorXd & _razel, const VectorXd & _geo, const double & _jd, const double & _lod, const double & _xp, const double & _yp)
	{
		return sezToTEME(razelToSEZ(_razel), _geo, _jd, _lod, _xp, _yp);
	}

	/**
	* <summary> True Equator, Mean Equinox to Range, Azimuth, Elevation. </summary>
	* <param name = "_teme"> VectorXd in South, East, Zenith frame. </param>
	* <param name = "_geo"> Position of the site in geodetic latitude, longitude, altitude. </param>
	* <param name = "_jd"> Julian Day number. </param>
	* <param name = "_lod"> Excessive length of day (LOD). </param>
	* <param name = "_xp"> X polar motion coefficient (rad). </param>
	* <param name = "_yp"> Y polar motion coefficient (rad). </param>
	* <returns> VectorXd in Range, Azimuth, Elevation frame. </returns>
	*/
	VectorXd temeToRAZEL(const VectorXd & _teme, const VectorXd & _geo, const double & _jd, const double & _lod, const double & _xp, const double & _yp)
	{
		return sezToRAZEL(temeToSEZ(_teme, _geo, _jd, _lod, _xp, _yp));
	}

	/**
	* <summary> Calculates the polar motion matrix according to the specified coefficients. </summary>
	* <param name = "xp"> x-axis polar motion coefficient (rad). </param>
	* <param name = "yp"> y-axis polar motion coefficient (rad). </param>
	* <returns> Polar motion rotation matrix. </returns>
	*/
	MatrixXd getPolarMotionMatrix(const double & _xp, const double & _yp)
	{
		MatrixXd pm(3, 3);

		double sinXp = sin(_xp), cosXp = cos(_xp), sinYp = sin(_yp), cosYp = cos(_yp);

		pm << cosXp, 0, -sinXp,       // Polar motion matrix
			sinXp * sinYp, cosYp, cosXp * sinYp,
			sinXp * cosYp, -sinYp, cosXp * cosYp;

		return pm;
	}

	/**
	* <summary> Get J2000 reference give Julian Date (JDN). </summary>
	* <param name = "_jdn"> Julian Day number. </para>m
	* <returns> The number of elapsed centuries since J2000 epoch. </returns>
	*/
	double getJ2000Ref(const double _jdn)
	{
		return (_jdn - T0) / 36525.0;
	}

	/**
	* <summary> Calculate Greenwich Mean Sidereal Time for a given date. </summary>
	* <param name = "_jd"> Julian Day number. </param>
	* <returns> Greenwich Mean Sidereal Time (radians). </returns>
	*/
	double getGMST(const double & _jd)
	{
		double T_UT1 = getJ2000Ref(_jd), thetaGMST, twoPi = 2 * M_PI;

		thetaGMST = -6.2e-6 * T_UT1 * T_UT1 * T_UT1 + 0.093104 * T_UT1 * T_UT1 + (876600.0 * 3600.0 + 8640184.812866) * T_UT1 + 67310.54841;
		thetaGMST = fmod(thetaGMST / 240.0 * M_PI / 180.0, twoPi);

		if (thetaGMST < 0.0)
			thetaGMST += twoPi;

		return thetaGMST;
	}

	/**
	* <summary> Calculate Julian Day Number for a given Gregorian calendar date. </summary>
	* <param name = "date"> Gregorian date. Astro::date structure. </param>
	* <returns> The Julian Day Number. </returns>
	*/
	double getJulianDay(const date & _date)
	{
		double B, C, JD;
		int month = _date.month;
		int year = _date.year;

		if (month == 1 || month == 2) {
			month = month + 12;
			year = year - 1;
		}

		int t = (int)floor(year * 0.01);		// Warning here

		B = 2 - t + floor(t * 0.25);

		C = (((_date.sec) / 60.0
			+ (double)_date.minute) / 60.0
			+ (double)_date.hour) / 24.0;

		JD = floor(365.25 * (year + 4716.0)) + floor(30.6001 * (month + 1))
			+ _date.day + B - 1524.5 + C;

		return JD;

	}

	/**
	* <summary> Calculate Julian Day Number for a given Gregorian calendar date. </summary>
	* <param name = "_year"> Year. </param>
	* <param name = "_month"> Month. </param>
	* <param name = "_day"> Day. </param>
	* <param name = "_hour"> Hour. </param>
	* <param name = "_minute"> Minute. </param>
	* <param name = "_sec"> Seconds (double precision). </param>
	* <returns> The Julian Day Number. </returns>
	*/
	double getJulianDay(const int & _year, const int & _month, const int & _day, const int & _hour, const int & _minute, const double & _sec)
	{
		return getJulianDay(date(_year, _month, _day, _hour, _minute, _sec));
	}

	/**
	* <summary> X-axis rotation matrix. </summary>
	* <param name = "_angle"> Rotation angle (radians). </angles>
	* <returns> X-axis rotation MatrixXd. </returns>
	*/
	Matrix3d rotX(const double& _angle) {
		double s = sin(_angle), c = cos(_angle);
		Matrix3d r;
		r << 1, 0, 0, 0, c, -s, 0, s, c;
		return r;
	}

	/**
	* <summary> Y-axis rotation matrix. </summary>
	* <param name = "_angle"> Rotation angle (radians). </angles>
	* <returns> Y-axis rotation matrix. </returns>
	*/
	Matrix3d rotY(const double& _angle) {
		double s = sin(_angle), c = cos(_angle);
		Matrix3d r;
		r << c, 0, s, 0, 1, 0, -s, 0, c;
		return r;
	}

	/**
	* <summary> Z-axis rotation matrix. </summary>
	* <param name = "_angle"> Rotation angle (radians). </angles>
	* <returns> Z-axis rotation matrix. </returns>
	*/
	Matrix3d rotZ(const double& _angle) {
		double s = sin(_angle), c = cos(_angle);
		Matrix3d r;
		r << c, -s, 0, s, c, 0, 0, 0, 1;
		return r;
	}

	/**
	 * <summary> Calculates the observation Jacobian using the FADBAD++ automatic differentiation library. </summary>
	 * <param name = "_sez"> VectorXd in the South-East-Zenith frame. </param>
	 * <param name = "_zDim"> The dimensionality of the observation. </param>
	 * <returns> MatrixXd containing the partial derivatives of transformation. </returns>
	 */
	MatrixXd getSEZToRAZELJacobianFADBAD(const VectorXd& _sez, const size_t& _zDim) {

		MatrixXd j = MatrixXd::Zero(_zDim, _sez.size());
		fadbad::F<double> x = _sez(0), y = _sez(1), z = _sez(2),
			xdot = _sez(3), ydot = _sez(4), zdot = _sez(5),
			rho, az, ele, rhodot, azdot, eledot;

		x.diff(0, 6);
		y.diff(1, 6);
		z.diff(2, 6);
		xdot.diff(3, 6);
		ydot.diff(4, 6);
		zdot.diff(5, 6);

		// Range
		rho = xSEZtoRAZEL(x, y, z, xdot, ydot, zdot);
		j(0, 0) = rho.d(0);
		j(0, 1) = rho.d(1);
		j(0, 2) = rho.d(2);

		// Azimuth
		az = ySEZtoRAZEL(x, y, z, xdot, ydot, zdot);
		j(1, 0) = az.d(0);
		j(1, 1) = az.d(1);
		j(1, 2) = az.d(2);

		// Elevation
		ele = zSEZtoRAZEL(x, y, z, xdot, ydot, zdot);
		j(2, 0) = ele.d(0);
		j(2, 1) = ele.d(1);
		j(2, 2) = ele.d(2);

		if (_zDim == 6) {
			
			// Range rate
			rhodot = xdotSEZtoRAZEL(x, y, z, xdot, ydot, zdot);
			j(3, 0) = rhodot.d(0);
			j(3, 1) = rhodot.d(1);
			j(3, 2) = rhodot.d(2);
			j(3, 3) = rhodot.d(3);
			j(3, 4) = rhodot.d(4);
			j(3, 5) = rhodot.d(5);

			// Azimuth rate
			azdot = ydotSEZtoRAZEL(x, y, z, xdot, ydot, zdot);
			j(4, 0) = azdot.d(0);
			j(5, 1) = azdot.d(1);
			j(6, 2) = azdot.d(2);
			j(7, 3) = azdot.d(3);
			j(8, 4) = azdot.d(4);
			j(9, 5) = azdot.d(5);

			// Elevation rate
			eledot = zdotSEZtoRAZEL(x, y, z, xdot, ydot, zdot);
			j(3, 0) = eledot.d(0);
			j(3, 1) = eledot.d(1);
			j(3, 2) = eledot.d(2);
			j(3, 3) = eledot.d(3);
			j(3, 4) = eledot.d(4);
			j(3, 5) = eledot.d(5);
		}

		return j;
	}

	/**
	* <summary> Propagates the step using numerical integration. </summary>
	* <param name = "_state"> The orbital state in the Earth-Centered Inertial frame. </param>
	* <param name = "_dt"> Desired timestep. </param>
	* <returns> VectorXd with the predicted state. </returns>
	*/
	VectorXd integrationPrediction(const VectorXd& _state, const double& _dt) {
		
		VectorXd result(6);

		state_type x0(6);

		x0[0] = _state(0);
		x0[1] = _state(1);
		x0[2] = _state(2);
		x0[3] = _state(3);
		x0[4] = _state(4);
		x0[5] = _state(5);

		boost::numeric::odeint::integrate(integrateOrbit, x0, 0.0, _dt, _dt / 100);

		result << x0[0], x0[1], x0[2], x0[3], x0[4], x0[5];

		return result;
	}

	/**
	* TODO: There's something wrong with the matrix
	 * <summary> Shepherd Matrix. </summary>
	 * <param name = "_m0"> Initial state vector (ECI frame). </param>
	 * <param name = "_t"> Desired timestep. </param>
	 * <param name = "_m"> A reference to store the predicted state vector. </param>
	 * <param name = "_mu"> Gravitational parameter. </param>
	 */
	MatrixXd getShepperdMatrix(const VectorXd & _m0, const double& _t, VectorXd & _m, const double& _mu)
	{
		MatrixXd theta = MatrixXd::Zero(_m0.size(), _m0.size());

		// Readability
		VectorXd mp0 = _m0.head(3), mv0 = _m0.tail(3), mp, mv;

		// Constants
		double r0 = mp0.norm(), nu0 = mp0.dot(mv0);
		double beta = 2 * _mu / r0 - mv0.dot(mv0);

		// Initialization
		double u = 0, deltaU = 0, P = 2 * M_PI * _mu * pow(beta, -1.5);
		int n = (int) floor(1.0 / P * (_t + P / 2.0 - 2.0 * nu0 / beta));
		deltaU = 2.0 * n * M_PI * pow(beta, -2.5);

		double t = 0, dt = t - _t, U, U0, U1, U2, U3, r;;
		size_t i = 0;
		bool qisbad = false;
		
		// Kepler iteration loop, beta > 0 (elliptic orbit)
		while (fabs(t - _t) > 0.0001) 
		{	
			i++;

			double q = beta * u * u / (1.0 + beta * u * u);

			if ((i > 25) || (q >= 1)) 
			{ 
				qisbad = true; 
				break; 
			}

			double G = 1.0, GPrev = 2.0, 
				A = 1.0, B = 1.0, nn = 0, k = -9, d = 15, l = 3;

			// Continued fraction
			while (fabs(G - GPrev) > 1e-14)
			{
				k = -k;
				l += 2;
				d += 4 * l;
				nn += (1 + k) * l;
				A = d / (d - nn * A * q);
				B *= (A - 1.0);
				GPrev = G;
				G += B;
			}

			double U0w2 = 1.0 - 2.0 * q,
				U1w2 = 2.0 * (1 - q) * u;

			U = 16.0 / 15.0 * pow(U1w2, 5) * G + deltaU;		// A.12
			U0 = 2.0 * U0w2 * U0w2 - 1;
			U1 = 2.0 * U0w2 * U1w2;
			U2 = 2.0 * U1w2 * U1w2;
			U3 = beta * U + U1 * U2 / 3.0;
			r = r0 * U0 + nu0 * U1 + _mu * U2;
			t = r0 * U1 + nu0 * U2 + _mu * U3;
			dt = t - _t;
			//u -= dt / 4.0 / (1-q) / r; 
			u -= dt / (1.0 - q) / (4.0*r + dt * beta * u);
		}
			
		if (qisbad) 
		{
			std::cerr << "Q not converged." << std::endl;
		} 
		else
		{
			double f = 1 - (_mu / r0) * U2;
			double g = r0 * U1 + nu0 * U2;
			double F = -_mu * U1 / (r * r0);
			double G = 1 - (_mu / r) * U2;
			
			mp = f * mp0 + g * mv0;
			mv = F * mp0 + G * mv0;
			_m << mp, mv;

			double W = g * U2 + 3.0 * _mu * U;

			MatrixXd M = MatrixXd::Zero(3,3);

			double r_2 = r * r, r0_2 = r0 * r0;
			double r_3 = r_2 * r, r0_3 = r0_2 * r0;

			M(0, 0) = (U0 / (r * r0) + 1.0 / r0_2 + 1.0 / r_2) * F - _mu * W / (r_3 * r0_3) ;
			M(0, 1) = F * U1 / r + (G - 1) / r_2;
			M(0, 2) = (G - 1) * U1 / r - _mu * W / r_3;

			M(1, 0) = -F * U1 / r0 - (f - 1) / r0_2;
			M(1, 1) = -F * U2;
			M(1, 2) = -(G - 1) * U2;

			M(2, 0) = (f - 1) * U1 / r0 - _mu * W / r0_3;
			M(2, 1) = (f - 1) * U2;
			M(2, 2) = g * U2 - W;

			MatrixXd I = MatrixXd::Zero(3,3), Phi_11, Phi_12, Phi_21, Phi_22;

			MatrixXd rv(2, 3), rvt, rv0(2, 3);
			rv0 << mp0.transpose(), mv0.transpose();
			rv << mp.transpose(), mv.transpose();
			rvt = rv.transpose();

			Phi_11 = f * I + rvt * M.block<2, 2>(1, 0) * rv0;
			Phi_12 = g * I + rvt * M.block<2, 2>(1, 1) * rv0;
			Phi_21 = F * I - rvt * M.block<2, 2>(0, 0) * rv0;
			Phi_22 = G * I - rvt * M.block<2, 2>(0, 1) * rv0;

			theta << Phi_11, Phi_12, Phi_21, Phi_22;
		}

		return theta;
	}
}