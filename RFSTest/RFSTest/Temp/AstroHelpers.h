#pragma once
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>

#define undefined 99999.1

using namespace Eigen;

/**
 * <summary> A set of basic astronomical and astrodynamical fucntions for the transformation of the celestial object's position. </summary>
 * <para> Topocentric Horizon, Earth-Centered, Earth-Fixed (ECEF) using World Geodetic System 1984 (WGS 84), Earth-Centered Inertial (ECI)
 * reference frames are present. Eigen library is used. </para>
 * <para> References: David. A. Vallado - Fundamentals of Astrodynamics and Applications, 4th Edition, 2013, Department of Defence World Geodetic
 * System 1984 Technical Report, 2000. </para>
 * <date> 2015-03-31 </date>
 * <author> Andrey Pak, VIBOT 8 </author>
 * 
 */
namespace AstroHelpers
{
    /** Temporary structure for storing time with millisecond precision **/
    struct date_ms {

        date_ms(): year(0), month(0), day(0), hour(0), min(0), sec(0)  {}
        date_ms(int yy, int mm, int dd, int hr, int nMin, double nSec):
            year(yy), month(mm), day(dd), hour(hr), min(nMin), sec(nSec) {}
		int year, month, day, hour, min;
		double sec;
	};

    /** Structure for storing measurement - Julian Day pair for Initial Orbit Determination **/
    struct iodmjd {

        iodmjd(): jd(0) {
            m = VectorXd::Zero(3);
        }
        iodmjd(double newJD, VectorXd newM): jd(newJD), m(newM) {}
        double jd;
        VectorXd m;
        bool operator < (const iodmjd &other) const { return jd < other.jd; }       // FOR DATE ONLY!!!
        bool operator == (const iodmjd &other) const { return  jd == other.jd && m == other.m; }

    };

    /** Constants **/
	const double T0			= 2451545.0;				// J2000 reference epoch (days)
    const double R_EQ		= 6378.137;					// Major semiaxis (equatorial radius)
	const double E			= 8.1819190842622e-2;		// First eccentricity
	const double E2			= 6.69437999014e-3;			// First eccentricity squared
    const double OMEGA_E	= 7.292115146706979e-5;	
    const double MU_E		= 398600.4418;				// Earth Gravitational Constant (km^3 / s^2)
	const double J2			= 0.00108263;				// J2 drag term

    /** Range, Azimuth, Elevation (Degrees) <-> South, East, Zenith (6) */
    VectorXd RAZELToSEZ(VectorXd razel);
	VectorXd SEZToRAZEL(VectorXd sez);

    /** South, East, Zenith <-> Earth-Centered, Earth-Fixed (6) */
	VectorXd RAZELToECEF(VectorXd sez, VectorXd geo, const int& type);
	VectorXd ECEFToRAZEL(VectorXd ecef, VectorXd geo, const int& type);

    /** Geodetic latitute, longitude and altitude to Earth-Centered, Earth-Fixed (6) */
    VectorXd GdToECEF(VectorXd geo, const int& type);                        

    /** Earth-Centered, Earth-Fixed <-> True Equator, Mean Equinox (ECI) */
	VectorXd ECEFToTEME(VectorXd ecef, const double& jd, const double& lod, const double& xp, const double& yp);
	VectorXd TEMEToECEF(VectorXd teme, const double& jd, const double& lod, const double& xp, const double& yp);

    /** Range, Azimuth, Elevation (Degrees) <-> True Equator, Mean Equinox */
    VectorXd RAZELToTEME(VectorXd radar, VectorXd geo, date_ms date, const double& xp, const double& yp,const int& type);
	VectorXd TEMEToRAZEL(VectorXd teme, VectorXd geo, date_ms date, const double& xp, const double& yp, const int& type);
\
    /** South, East, Zenith <-> True Equator, Mean Equinox (ECI) **/
    MatrixXd SEZToTEME(VectorXd sez, VectorXd geo, date_ms date, const double& xp, const double& yp, const int& type);
    MatrixXd TEMEToSEZ(VectorXd teme, VectorXd geo, date_ms date, const double& xp, const double& yp, const int& type);

    /** South, East, Zenith <-> Earth-Centered, Earth-fixed **/
    MatrixXd SEZToECEF (VectorXd sez, VectorXd geo, const int& type);
    MatrixXd ECEFToSEZ (VectorXd ecef, VectorXd geo, const int& type);

    /** Functions for EKF covariance transformations */
    MatrixXd getSEZToTEMECovTfMat(VectorXd geo, date_ms date, const double& xp, const double& yp, const double& dim);
    MatrixXd getTEMEToSEZCovTfMat(VectorXd geo, date_ms date, const double& xp, const double& yp, const double& dim);

    /** Functions directly adapted from D. Vallado's C++ code **/
    VectorXd rv2coe(VectorXd rv);
    VectorXd coe2rv(VectorXd coe);
    VectorXd newtonNu(const double& ecc, const double& nu);
	VectorXd newtonM(const double& ecc, const double& m);
	VectorXd pKepler(VectorXd rv, const double& dt, const double& nDot, const double& nDDot);
    Vector3d iodHGibbs(Vector3d, Vector3d, Vector3d, const double&, const double&, const double&);

    /** Rotation matrices **/
    Matrix3d rotX(const double& angleRad);
	Matrix3d rotY(const double& angleRad);									// Checked	30/03/2015
	Matrix3d rotZ(const double& angleRad);									// Checked	30/03/2015

    /** Angle between vectors **/
    double vAngle(VectorXd v1, VectorXd v2);
	
    /** Polar motion matrix */
    Matrix3d polarMotion(const double& xp, const double& yp);
	double getDayFraction(const int& hour, const int& minute, const double& sec);	// Checked	30/03/2015
    double getJulianDay(date_ms date);                                              // Checked	30/03/2015
	double getJulianDay(int year, int month, int day, int hour, int minute, double sec);			
    double getJ2000Ref(double jd);                                                  // Checked	30/03/2015
    double getGMST(double jd);                                                      // Checked (1)	30/03/2015

    /** Temporary **/
    double diffTimeHMS(date_ms start, date_ms end);

}

double getHellingerDistance(VectorXd v1, MatrixXd p1, VectorXd v2, MatrixXd p2);
