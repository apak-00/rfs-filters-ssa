#include "AstroHelpers.h"

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
namespace AstroHelpers {

    /**
    * <summary> X-axis rotation matrix. </summary>
    * Angles are specified in radians.
    */
    Matrix3d rotX(const double& angle) {
        double s = sin(angle), c = cos(angle);
        Matrix3d r;
        r << 1, 0, 0, 0, c, -s, 0, s, c;
        return r;
    }

    /**
    * <summary> Y-axis rotation matrix. </summary>
    * Angles are specified in radians.
    */
    Matrix3d rotY(const double& angle) {
        double s = sin(angle), c = cos(angle);
        Matrix3d r;
        r << c, 0, s, 0, 1, 0, -s, 0, c;
        return r;
    }

    /**
    * <summary> Z-axis rotation matrix. </summary>
    * Angles are specified in radians.
    */
    Matrix3d rotZ(const double& angle) {
        double s = sin(angle), c = cos(angle);
        Matrix3d r;
        r << c, -s, 0, s, c, 0, 0, 0, 1;
        return r;
    }

	/**
	 * <summary> Calculates the angle between two vectors. </summary>
	 * <returns> The angle between two specefied vectors of the same size. </returns>
	 */
    double vAngle(VectorXd v1, VectorXd v2) {
        assert(v1.size() == v2.size());
		double temp = v1.dot(v2) / (v1.norm() * v2.norm());
        return acos(temp);
    }
	
	/**
	 * <summary> Calculates the polar motion matrix according to the specified coefficients. </summary>
	 * <param name = "xp"> x-axis polar motion coefficient </param>
	 * <param name = "yp"> y-axis polar motion coefficient </param>
	 * <returns> Polar motion rotation matrix. </returns>
	 */
    Matrix3d polarMotion(const double& xp, const double& yp) {

        Matrix3d pm;

        double sinXp = sin(xp), cosXp = cos(xp), sinYp = sin(yp), cosYp = cos(yp);

        pm <<   cosXp,              0,        -sinXp,       // Polar motion matrix
                sinXp * sinYp,  cosYp, cosXp * sinYp,
                sinXp * cosYp, -sinYp, cosXp * cosYp;

        return pm;
    }

    /**
     * <summary> Converts spherical to Cartesian coordinates. </summary>
     * <param name = "sph"> VectorXd containing spherical coordinates (range, azimuth, elevation,
     * range rate, azimuth rate, elevation rate). Angles are specified in degrees,
     * rates are specified in degrees/second. </param>
     *
     * <returns> VectorXd containing Cartesian coordinates (x, y, z, vx, vy, vz). </returns>
     */
    VectorXd RAZELToSEZ(VectorXd sph) {

        VectorXd cart;

        if (sph.size() == 6) {

            cart = VectorXd::Zero(6);
            double sinEl, sinAz, cosEl, cosAz, degToRad = M_PI / 180.0;

            sph(1) *= degToRad; sph(2) *= degToRad;
            sph(4) *= degToRad; sph(5) *= degToRad;

            sinEl = sin(sph(2)); cosEl = cos(sph(2));
            sinAz = sin(sph(1)); cosAz = cos(sph(1));

            cart(0) = - sph(0) * cosEl * cosAz;
            cart(1) = sph(0) * cosEl * sinAz;
            cart(2) = sph(0) * sinEl;

            cart(3) = - cosEl * cosAz * sph(3)  +  cart(2) * cosAz * sph(5)  +  cart(1) * sph(4);
            cart(4) =  cosEl * sinAz * sph(3) - cart(2)*sph(5)*sinAz - cart(0)*sph(4);
            cart(5) =  sph(3)*sinEl + sph(0) * sph(5) * cosEl;

        } else if (sph.size() == 3) {

            cart = VectorXd::Zero(3);
            double sinEl, sinAz, cosEl, cosAz, degToRad = M_PI / 180.0;

            sph(1) *= degToRad; sph(2) *= degToRad;

            sinEl = sin(sph(2)); cosEl = cos(sph(2));
            sinAz = sin(sph(1)); cosAz = cos(sph(1));

            cart(0) = - sph(0) * cosEl * cosAz;
            cart(1) = sph(0) * cosEl * sinAz;
            cart(2) = sph(0) * sinEl;
        }

        return cart;
    }

	/**
	 * <summary> Converts Cartesian to spherical coordinates. </summary>
	 * <param name = "sez"> VectorXd contaiing Cartesian coordinates for position and velocity. </param>
	 * <returns> VectorXd containing spherical site coordinates (range, azimuth, elevation and their rates). </returns>
	 */
	VectorXd SEZToRAZEL(VectorXd sez) {

		VectorXd sph(6);
		double small = 0.0000001;

		/* ------------ calculate azimuth and elevation ------------- */
		double xy = sqrt(sez(0) * sez(0) + sez(1) * sez(1)), dxdy;

		Vector3d r, v;
		r << sez(0), sez(1), sez(2);
		v << sez(3), sez(4), sez(5);

		// Range
		sph(0) = r.norm();

		// Azimuth
		if (fabs(r(1)) < small) {
			if (xy < small) {
				dxdy = sqrt(v(0) * v(0) + v(1) * v(1));
				sph(1) = atan2(v(1) / dxdy, -v(0) / dxdy);
			}
			else 
			{
				if (r(0) > 0.0)
					sph(1) = M_PI;
				else
					sph(1) = 0.0;
			}
		}
		else 
			sph(1) = atan2(r(1) / xy, -r(0) / xy);

		if (sph(1) < 0)
			sph(1) += M_PI * 2;
		
		// Elevation
		if (xy < small)  // directly over the north pole
			sph(2) = copysign(1.0, r(2)) * M_PI / 2; // +- 90
		else
			sph(2) = asin(r(2) / sph(0));

		/* ----- Calculate range, azimuth and elevation rates ------- */
		// Range rate
		sph(3) = r.dot(v) / sph(0);

		// Azimuth rate
		if (fabs(xy * xy) > small)
			sph(4) = (v(0) * r(1) - v(1) * r(0)) / (xy * xy);
		else
			sph(4) = 0.0;

		// Elevation rate
		if (fabs(xy) > small)
			sph(5) = (v(2) - sph(3) * sin(sph(2))) / xy;
		else
			sph(5) = 0.0;

		double rad2deg = 180.0 / M_PI;

		// Temporary
		sph(1) *= rad2deg;
		sph(2) *= rad2deg;
		sph(4) *= rad2deg;
		sph(5) *= rad2deg;
 
		return sph;
	}

    /**
     * <summary> Transforms geodetic coordinates (latitude, longitude, height) to cartesian coordinates. </summary>
     * <param name = "geo"> VectorXd containing geodetic coordinates (latitude, longitude, height over ellisoid in meters). </param>
     * <param name = "type"> Types: 1 - spherical (not implemented), 2 - World Geodetic System 1984 (WGS 84) reference ellipsoid. </param>
     * <returns> VectorXd containing Cartesian Coordinates (x, y, z). </returns>
     */
    VectorXd GdToECEF(VectorXd geo, const int& type) {

        VectorXd cart(6);

        double latRad = geo(0) / 180.0 * M_PI, lonRad = geo(1) / 180.0 * M_PI;	// Degrees to radians
        double s = sin(latRad), c = cos(latRad);

        if (type == 1) {			// Simple spherical, not implemented
        } else if (type == 2) {		// WGS 84

            double N = R_EQ / sqrt(1.0 - E2 * pow(s,2) );			// Radius of curvature in prime meridian
            double t = (N + geo(2)) * c;
            cart << t * cos(lonRad), t * sin(lonRad), ((1 - E2) * N + geo(2)) * s, 0, 0, 0;     // Zero velocity
        }

        return cart;
    }

    /**
    * <summary> Transforms SEZ (Topocentric Horizon) to ECEF coordinate system of a given type. </summary>
    * <param name = "sez"> South-East-Zenith spherical coordinates (range, azimuth, elevation). </param>
    * <param name = "geo"> Geodetic coordinates (latitude, longitude). </param>
    * <param name = "type"> Type of geographical transformation used. Type 1 - geocentric (not implemented).
    * Type 2 - World Geodetic System 1984 (WGS 84) reference ellipsoid</param>
    * <returns> A VectorXd containing Cartesian coordinates. </returns>
    */
    VectorXd RAZELToECEF(VectorXd razel, VectorXd geo, const int& type) {

        VectorXd ecef(6), tempR(3), tempV(3), tempS(6);
        MatrixXd rot(3,3);

        double latRad = geo(0) / 180.0 * M_PI;
        double lonRad = geo(1) / 180.0 * M_PI;

        tempS = RAZELToSEZ(razel);

        tempR << tempS(0), tempS(1), tempS(2);
        tempV << tempS(3), tempS(4), tempS(5);
        rot = rotZ(lonRad) * rotY(M_PI / 2.0 - latRad);

        ecef << rot * tempR, rot * tempV;

        VectorXd observer = GdToECEF(geo, type);
        ecef += observer;

        return ecef;
    }

	/**
	* <summary> Transforms ECEF  to SEZ (Topocentric Horizon) coordinate system. </summary>
	* <param name = "ecef"> Earth-Centered, Earth-Fixed frame coordinates (position, velocity) </param>
	* <param name = "geo"> Geodetic coordinates (latitude, longitude). </param>
	* <param name = "type"> Type of geographical transformation used. Type 1 - geocentric (not implemented).
	* Type 2 - World Geodetic System 1984 (WGS 84) reference ellipsoid</param>
	* <returns> A VectorXd containing ECEF coordinates. </returns>
	*/
	VectorXd ECEFToRAZEL(VectorXd ecef, VectorXd geo, const int& type) {
		VectorXd tempR(3), tempV(3), temp(6);
		MatrixXd rot(3, 3);

		double latRad = geo(0) / 180.0 * M_PI;
		double lonRad = geo(1) / 180.0 * M_PI;

		temp = ecef - GdToECEF(geo, type);

		rot = rotY(-(M_PI / 2.0 - latRad)) * rotZ(-lonRad);
		
		tempR << temp(0), temp(1), temp(2);
		tempV << temp(3), temp(4), temp(5);

		temp << rot * tempR, rot * tempV;

		return SEZToRAZEL(temp);
	}

    /**
    * <summary> Converts ECEF (WGS 84) coordinates to ECI TEME for a certain date (rotation around Z axis). </summary>
    * <returns> A VectorXd with Earth-Centered Inertial (X axis towards mean equinox) Cartesian coordinates. </returns>
    */
    VectorXd ECEFToTEME(VectorXd ecef, const double& jd, const double& lod, const double& xp, const double& yp) {

        double gmst = getGMST(jd), thetaSa = 7.29211514670698e-05 * (1.0  - lod / 86400.0);
        Matrix3d pm = polarMotion(xp, yp);
        Matrix3d rotGMST = rotZ(gmst);
        VectorXd eciteme = VectorXd::Zero(6);
        Vector3d tempR, tempV, tempRTEME, tempVTEME, rPEF, tempThetaSa;

        tempThetaSa << 0, 0, thetaSa;

        tempR << ecef(0), ecef(1), ecef(2);
        tempV << ecef(3), ecef(4), ecef(5);

		rPEF = pm * tempR;		// Position in Pseudo Earth-Fixed Frame

        tempRTEME = rotGMST * rPEF;
        tempVTEME = rotGMST * (pm * tempV + tempThetaSa.cross(rPEF));

        eciteme << tempRTEME, tempVTEME;

        return eciteme;
    }

	/**
	* <summary> Converts TEME (True Equator, Mean Equinox) coordinates to ECEF (WGS 84) for a certain date (rotation around Z axis). </summary>
	* <returns> A VectorXd with ECEF Cartesian coordinates. </returns>
	*/
	VectorXd TEMEToECEF(VectorXd teme, const double& jd, const double& lod, const double& xp, const double& yp) {
		
		double gmst = getGMST(jd), thetaSa = 7.29211514670698e-05 * (1.0 - lod / 86400.0);
		VectorXd ecef(6); 
		Vector3d tempRTEME, tempVTEME, tempR, tempV, tempThetaSa, rPEF, vPEF;
		Matrix3d pm = polarMotion(xp, yp).inverse();
		Matrix3d rotGMST = rotZ(-gmst);

		tempThetaSa << 0, 0, thetaSa;

		tempRTEME << teme(0), teme(1), teme(2);
		tempVTEME << teme(3), teme(4), teme(5);

		tempR = rotGMST * tempRTEME;
		tempV = rotGMST * tempVTEME - tempThetaSa.cross(tempR);

		tempR = pm * tempR;
		tempV = pm * tempV;

		ecef << tempR, tempV;

		return ecef;
	}

    /**
     * <summary> Calculate elapsed day fraction. </summary>
     * <param name = "hour"> Elapsed hour. </param>
     * <param name = "min"> Elapsed minutes </param>
     * <param name = "sec"> Elapsed seconds (ms precision allowed). </param>
     * <returns>Fraction of the elapsed part of the day. </returns>
     */
    double getDayFraction(const int& hour, const int& min, const double& sec) {
        return ((double)hour * 3600.0 + (double)min * 60 + sec) / 86400.0;
    }

    /**
     * <summary> Calculate Julian Day Number for a given Gregorian calendar date. </summary>
     * <param name = "date"> Gregorian date. </param>
     * <returns> The Julian Day Number. </returns>
     */
    double getJulianDay(date_ms date) {

        double B, C, JD;
        int month = date.month;
        int year = date.year;

        if (month == 1 || month == 2) {
            month = month + 12;
            year = year - 1;
        }

        int t = (int)floor(year / 100.0);		// Warning here

        B = 2 - t + floor(t / 4);

        C = (((date.sec) / 60.0
            + (double)date.min) / 60.0
            + (double)date.hour) / 24.0;

        JD = floor(365.25 * (year + 4716.0)) + floor(30.6001 * (month + 1))
            + date.day + B - 1524.5 + C;

        return JD;
    }

	/**
	* <summary> Calculate Julian Day Number for a given Gregorian calendar date. </summary>
	* <returns> The Julian Day Number. </returns>
	*/
	double getJulianDay(int year, int month, int day, int hour, int minute, double sec) {

		double B, C, JD;

		if (month == 1 || month == 2) {
			month = month + 12;
			year = year - 1;
		}

		int t = (int)floor(year / 100.0);		// Warning here

		B = 2 - t + floor(t / 4);

		C = (((sec) / 60.0
			+ (double)minute) / 60.0
			+ (double)hour) / 24.0;

		JD = floor(365.25 * (year + 4716.0)) + floor(30.6001 * (month + 1))
			+ day + B - 1524.5 + C;

		return JD;
	}

    /**
     * <summary> Get J2000 reference give Julian Date (JDN). </summary>
     * <returns> The number of elapsed centuries since J2000 epoch. </returns>
     */
    double getJ2000Ref(double jdn) {
        return (jdn - T0) / 36525.0;
    }

    /**
     * <summary> Calculate Greenwich Mean Sidereal Time for a given date. </summary>
     * <returns> Greenwich Mean Sidereal Time (radians). </returns>
     */
    double getGMST(double jd) {

        double T_UT1 = getJ2000Ref(jd), thetaGMST, twoPi = 2 * M_PI;

        thetaGMST = -6.2e-6 * T_UT1 * T_UT1 * T_UT1 + 0.093104 * T_UT1 * T_UT1 + (876600.0 * 3600.0 + 8640184.812866) * T_UT1 + 67310.54841;
        thetaGMST = fmod(thetaGMST / 240.0 * M_PI / 180.0, twoPi);

        if (thetaGMST < 0.0)
            thetaGMST += twoPi;

        return thetaGMST;
    }

	/**
	 * <summary> Converts topocentric spherical radar coordinates (range, azimuth, elevation) to 
	 * True Equator, Mean Equinox Frame (Earth-Centered Inertial). </summary>
	 * <param name = "radar"> VectorXd containing range, azimuth, elevation. </param>
	 * <param name = "geo"> VectorXd containing geodetic latitude, longitude, altitude </param>
	 * <param name = "xp, yp"> Polar motion coefficient </param>
	 * <param name = "type"> Type for ECEF transformation (use 2 for WGS-84) </param>
	 * <returns> VectorXd containing teme coordinates for position and veloctiy. </returns>
	 */
    VectorXd RAZELToTEME(VectorXd radar, VectorXd geo, date_ms date, const double& xp, const double& yp,const int& type) {
        double jd = getJulianDay(date);
		VectorXd ecef = RAZELToECEF(radar, geo, type), teme;
        teme = ECEFToTEME(ecef, jd, 0, xp, yp);
		return teme;
    }

	/**
	* <summary> Converts True Equator, Mean Equinox Frame (Earth-Centered Inertial) coordinates to 
	* topocentric spherical radar coordinates (range, azimuth, elevation). </summary>
	* <param name = "teme"> TEME coordinates (postion, velocity). </param>
	* <param name = "geo"> VectorXd containing geodetic latitude, longitude, altitude </param>
	* <param name = "xp, yp"> Polar motion coefficient </param>
	* <param name = "type"> Type for ECEF transformation (use 2 for WGS-84) </param>
	* <returns> VectorXd containing topocentric horizon coordinates for (range, azimuth, elevation and their rates). </returns>
	*/
	VectorXd TEMEToRAZEL(VectorXd teme, VectorXd geo, date_ms date, const double& xp, const double& yp, const int& type) {
		double jd = getJulianDay(date);
		VectorXd ecef = TEMEToECEF(teme, jd, 0, xp, yp), razel;
		razel = ECEFToRAZEL(ecef, geo, 2);
		return razel;
	}

	/**
	 * <summary> Solves Kepler's equation with known true anomaly. </summary>
	 * Reference: David A. Vallado, Fundamentals of Astrodynamics and Applications, 4th ed., 2013
	 * <param name = ecc> Eccentircity </param>
	 * <param name = nu> True anomaly (radians) </param>
	 * <returns> VectorXd containing eccentric and mean anomalies. </returns>
	 */
    VectorXd newtonNu(const double& ecc, const double& nu) {
        double small = 0.00000001, sine, cose;
        double e0 = 999999.9, m = 999999.9;

        // --------------------------- circular ------------------------
        if ( fabs(ecc) < small  ) {
            m  = nu;
            e0 = nu;
        } else {
        // ---------------------- elliptical -----------------------
            if ( ecc < 1.0-small  ) {
                sine = (sqrt( 1.0 - ecc * ecc) * sin(nu) ) / ( 1.0 +ecc*cos(nu) );
                cose = (ecc + cos(nu)) / (1.0 + ecc * cos(nu));
                e0   = atan2(sine, cose);
                m    = e0 - ecc * sin(e0);
            } else
                 // -------------------- hyperbolic  --------------------
                if (ecc > 1.0 + small) {
                    if ((ecc > 1.0) && (fabs(nu) + 0.00001 < M_PI - acos(1.0/ecc))) {
                        sine = (sqrt(ecc * ecc - 1.0) * sin(nu)) / (1.0 + ecc*cos(nu));
                        e0   = asinh(sine);
                        m    = ecc * sinh(e0) - e0;
                    }
                } else
                    // ----------------- parabolic ---------------------
                    if (fabs(nu) < 168.0 * M_PI / 180.0) {
                        e0 = tan(nu * 0.5);
                        m  = e0 + (e0*e0*e0) / 3.0;
                    }
        }

        if (ecc < 1.0) {
            m = fmod(m, 2.0 *M_PI);
            if (m < 0.0)
                m = m + 2.0 *M_PI;
            e0 = fmod(e0, 2.0 *M_PI);
        }

        VectorXd result(2);
        result << e0, m;
        return result;
    }

	/**
	 * <summary> Finds the eccentric anomaly given the mean anomly using the Newton-Rhapson method. </summary>
	 * <param name = "ecc"> Eccentricity </param>
	 * <param name = "m"> Mean anomaly </param>
	 * <returns> VectorXd containing eccentric anomaly (e0), true anomaly (nu) </returns>
	 */
	VectorXd newtonM(const double& ecc, const double& m) {
		// Return values
		double e0, nu;

		const int numiter = 50;
		const double small = 0.00000001;       // small value for tolerances

		double e1, sinv, cosv, r1r = 0.0;
		int ktr;

		/* -------------------------- hyperbolic  ----------------------- */
		if ((ecc - 1.0) > small) {
			/* ------------  initial guess ------------- */
			if (ecc < 1.6)
			if (((m < 0.0) && (m > -M_PI)) || (m > M_PI))
				e0 = m - ecc;
			else
				e0 = m + ecc;
			else
			if ((ecc < 3.6) && (fabs(m) > M_PI)) // just edges)
				e0 = m - copysign(1.0, m) * ecc;
			else
				e0 = m / (ecc - 1.0); // best over 1.8 in middle

			ktr = 1;
			e1 = e0 + ((m - ecc * sinh(e0) + e0) / (ecc * cosh(e0) - 1.0));

			while ((fabs(e1 - e0) > small) && (ktr <= numiter)) {
				e0 = e1;
				e1 = e0 + ((m - ecc * sinh(e0) + e0) / (ecc * cosh(e0) - 1.0));
				ktr++;
			}

			/* ---------  find true anomaly  ----------- */
			sinv = -(sqrt(ecc * ecc - 1.0) * sinh(e1)) / (1.0 - ecc * cosh(e1));
			cosv = (cosh(e1) - ecc) / (1.0 - ecc * cosh(e1));
			nu = atan2(sinv, cosv);
		}
		else {
			/* ---------------------- parabolic ------------------------- */
			if (fabs(ecc - 1.0) < small) {
				e0 = r1r;
				ktr = 1;
				nu = 2.0 * atan(e0);
			}
			else{
				/* --------------------- elliptical --------------------- */
				if (ecc > small) {
					/* ------------  initial guess ------------- */
					if (((m < 0.0) && (m > -M_PI)) || (m > M_PI))
						e0 = m - ecc;
					else
						e0 = m + ecc;

					ktr = 1;
					e1 = e0 + (m - e0 + ecc * sin(e0)) / (1.0 - ecc * cos(e0));

					while ((fabs(e1 - e0) > small) && (ktr <= numiter))
					{
						ktr++;
						e0 = e1;
						e1 = e0 + (m - e0 + ecc * sin(e0)) / (1.0 - ecc * cos(e0));
					}

					/* ---------  find true anomaly  ----------- */
					sinv = (sqrt(1.0 - ecc * ecc) * sin(e1)) / (1.0 - ecc * cos(e1));
					cosv = (cos(e1) - ecc) / (1.0 - ecc * cos(e1));
					nu = atan2(sinv, cosv);
				}
				else {
					/* --------------------- circular --------------------- */
					ktr = 0;
					nu = m;
					e0 = m;
				}
			}
		}

		if (ktr > numiter)
			printf("newtonrhapson not converged in %3d iterations\n", numiter);

		VectorXd result(2);
		result << e0, nu;
		return result;
		
	}    // procedure newtonm


	/**
	 * <summary> Finds classical orbital elements given the geocentric equatorial position and velocity vectors. </summary>
	 * <param name = "rv"> VectorXd containing equatorial position and velocity values (km, km/s).  </param>
	 * <param name = "mu"> Gravitational constant (km^3 / (kg*s)). </param>
	 * <returns> VectorXd containing semilatus rectum (p), semimajor axis (a), eccentricity (e), inclination (incl), 
	 * longitude of the ascending node (omega), argument of perigee (argp), true anomaly(nu), mean anomaly(m), 
	 * argument of latitude (arglat), true longitude(truelon), longitude of periapsis (lonper). </returns>
	 */
    VectorXd rv2coe(VectorXd rv) {

        Vector3d r, v, hBar, nBar, eBar;
        double magR, magV, magH, magN, rDotV, c1, small, sme, hk, temp,
                twoPi = M_PI * 2.0, halfPi = M_PI / 2.0;
		char orbitType[3];
		strcpy(orbitType, "ei");

        // Return values
        double p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper;

        small  = 0.00000001;

        r << rv(0), rv(1), rv(2);
        v << rv(3), rv(4), rv(5);
        hBar = r.cross(v);

        magR = r.norm();
        magV = v.norm();
        magH = hBar.norm();

         // -----------  Find h n and e vectors   ----------------------
        if (magH > small) {
            nBar << -hBar(1), hBar(0), 0.0;
            magN = nBar.norm();
            c1 = magV * magV - MU_E / magR;
            rDotV = r.dot(v);
            eBar = (c1 * r - rDotV * v) / MU_E;
            ecc = eBar.norm();


            // ------------  Find a e and semi-latus rectum   --------------
            sme = magV * magV * 0.5 - MU_E / magR;
            if (fabs(sme) > small)
				a = -MU_E / (2.0 * sme);
            else {
                a = undefined;
                //std::cout << "Wrong a" << std::endl;
            }
            p = magH * magH / MU_E;

            // ------------  Find inclination   ----------------------------
            hk = hBar(2) / magH;
            incl = acos(hk);

            // ------------  Determine type of orbit for later use  --------
            // ------------  Elliptical, Parabolic, Hyperbolic inclined ----

            if (ecc < small) {
                // Circular equatorial
                if ((incl < small) || (fabs(incl - M_PI) < small))
					strcpy(orbitType, "ce");
                // Circular inclined
                else
					strcpy(orbitType, "ci"); 
            } else
                // Elliptical, parablic, hyperbolic equatorial
                if ((incl < small) || (fabs(incl - M_PI) < small))
					strcpy(orbitType, "ee");

            // ------------  Find longitude of ascending node --------------
            if (magN > small) {

                temp = nBar(0) / magN;
                if (fabs(temp) > 1.0)
                    temp = copysign(1.0, temp);
                omega = acos(temp);
                if (nBar(1) < 0.0)
                    omega = twoPi - omega;
            } else {
                omega = undefined;
                //std::cout << "Invalid omega" << std::endl;
            }

            // ------------  Find argument of perigee ----------------------
            if (strcmp(orbitType, "ei") == 0) {
                argp = vAngle(nBar, eBar);
                if (eBar(2) < 0.0)
                    argp = twoPi - argp;
            } else {
                argp = undefined;
                //std::cout << "Invalid argument of perigee" << std::endl;
            }

            // -------------  Find true anomaly at epoch    ----------------
            if (orbitType[0] == 'e' ) {

                nu =  vAngle( eBar,r);
                if ( rDotV < 0.0  )
                    nu = twoPi - nu;
            } else {
                nu = undefined;
                //std::cout << "Invalid true anomaly (nu)" << std::endl;
            }

            // -------------  Find argument of latitude - circular inclined -
            if ( strcmp(orbitType, "ci") == 0 ) {

                arglat = vAngle( nBar,r );
                if ( r(2) < 0.0  )
                    arglat= twoPi - arglat;
                m = arglat;
            } else {
                arglat = undefined;
                //std::cout << "Invalid argument of latitude ci" << std::endl;
            }

            // --------------  Find longitude of perigee - elliptical equatorial ----
            if  (( ecc > small ) && (strcmp(orbitType,"ee") == 0)) {

                temp = eBar(0) / ecc;
                if ( fabs(temp) > 1.0  )
                    temp = copysign(1.0, temp);
                lonper = acos( temp );

                if ( eBar(1) < 0.0  )
                    lonper = twoPi - lonper;
                if ( incl > halfPi )
                    lonper = twoPi - lonper;
            } else {
                lonper = undefined;
                //std::cout << "Invalid argument of longitude of perigee ee" << std::endl;
            }

            // -------- Find true longitude - circular equatorial ------
            if  (( magR > small ) && ( strcmp(orbitType, "ce") == 0 )) {

                temp = r(0) / magR;
                if (fabs(temp) > 1.0)
                    temp = copysign(1.0, temp);
                truelon= acos(temp);

                if ( r(1) < 0.0  )
                 truelon = twoPi - truelon;
                if ( incl > halfPi )
                 truelon= twoPi - truelon;
                 m = truelon;
            } else {
                truelon = undefined;
                //std::cout << "Invalid argument of longitude of perigee ee" << std::endl;
            }

            if (orbitType[0] == 'e') {
                VectorXd em = newtonNu(ecc, nu);
				m = em(1);
            }

        } else {

            p       = undefined;
            a       = undefined;
            ecc     = undefined;
            incl    = undefined;
            omega   = undefined;
            argp    = undefined;
            nu      = undefined;
            m       = undefined;
            arglat  = undefined;
            truelon = undefined;
            lonper  = undefined;

            std::cout << "Failed to transform classical orbital elements. " << std::endl;
        }

        VectorXd coe(11);
        coe << p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper;

        return coe;

    }

	/**
	 * <summary> Finds the position and velocity vectors in the equatorial frame given the set of
	 * classical orbital elements. </summary>
	 * <param name = "coe"> VectorXd containing semilatus rectum (p), eccentricity (ecc), inclination (incl),
	 * longitude of the ascending node (omega), argument of perigee (argp), true anomaly(nu), 
	 * argument of latitude(arglat), true longitude (truelon), longitude of periapsis (lonper) </param>
	 * <returns> VectorXd containing position and velocity values in the equatorial frame. </returns>
	 */
    VectorXd coe2rv(VectorXd coe) {
        double p = coe(0), ecc = coe(1), incl = coe(2), omega = coe(3), argp = coe(4), nu = coe(5),
                arglat = coe(6), truelon = coe(7), lonper = coe(8);

        double small = 0.0000001,  temp, sinnu, cosnu;
        Vector3d rPQW, vPQW, tempVec, r, v;

        // --------------------  implementation   ----------------------
        //       determine what type of orbit is involved and set up the
        //       set up angles for the special cases.
        // -------------------------------------------------------------
        if (ecc < small) {
           // ----------------  circular equatorial  ------------------
           if ((incl < small) | (fabs(incl - M_PI) < small)) {
               argp  = 0.0;
               omega = 0.0;
               nu    = truelon;
           } else {
               // --------------  circular inclined  ------------------
               argp = 0.0;
               nu   = arglat;
           }
        } else {
           // ---------------  elliptical equatorial  -----------------
           if (( incl < small) | (fabs(incl - M_PI) < small)) {
               argp  = lonper;
               omega = 0.0;
           }
        }

        // ----------  form pqw position and velocity vectors ----------
        cosnu = cos(nu);
        sinnu = sin(nu);
        temp = p / (1.0 + ecc*cosnu);
        rPQW << temp*cosnu, temp*sinnu, 0.0;
        if (fabs(p) < 0.00000001)
           p = 0.00000001;

        vPQW <<  -sinnu    * sqrt(MU_E/p), (ecc + cosnu) * sqrt(MU_E/p), 0.0;

        // ----------------  perform transformation to ijk  ------------
        tempVec = rotZ(argp) * rPQW;
        tempVec = rotX(incl) * tempVec;
        r = rotZ(omega) * tempVec;

        tempVec = rotZ(argp) * vPQW;
        tempVec = rotX(incl) * tempVec;
        v = rotZ(omega) * tempVec;

        /*
        // Original code
        rot3( rpqw   , -argp , tempvec );
        rot1( tempvec, -incl , tempvec );
        rot3( tempvec, -omega, r     );

        rot3( vpqw   , -argp , tempvec );
        rot1( tempvec, -incl , tempvec );
        rot3( tempvec, -omega, v     );
        */

        VectorXd result(6);
        result << r, v;
        return result;

    }

	/**
	 * <summary> Propagates satellite position and velocity vector over a give time period, 
	 * taking into account J2 pertubations. </summary>
	 * <param name = "rv"> VectorXd containing position and velocity values in the equatorial frame. </param>
	 * <param name = "dt"> Timestep in seconds. </param>
	 * <param name = "nDot"> Time rate of change of n. </param>
	 * <param name = "nDDot"> Time acceleration of change of n </param>
	 * <returns> VectorXd with the propagated step. </returns>
	 */
	VectorXd pKepler(VectorXd rv, const double& dt, const double& nDot, const double& nDDot) {
		
		VectorXd coe = rv2coe(rv);

        if (rv(0) == undefined) {
            //std::cout << "Failed rv2coe" << std::endl;
        }

		// coe << p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper;
		double p = coe(0), a = coe(1), ecc = coe(2), incl = coe(3), omega = coe(4), argp = coe(5),
			nu = coe(6), m = coe(7), arglat = coe(8), truelon = coe(9), lonper = coe(10);
		
		// Local variables
		double n, j2op2, omegaDot, argPDot, mDot, trueLonDot, argLatDot, lonPerDot,
			small = 0.00000001, twoPi = 2 * M_PI;

        double tt = pow(a, 3.0);
        n = sqrt(MU_E / tt);
		j2op2 = (n * 1.5 * R_EQ * R_EQ * J2) / (p*p);
		omegaDot = -j2op2 * cos(incl);
		argPDot = j2op2 * (2.0 - 2.5*sin(incl)*sin(incl));
		mDot = n;

		a = a - 2.0 * nDot * dt * a / (3.0 * n);
		ecc = ecc - 2.0*(1.0 - ecc)*nDot*dt / (3.0*n);
		p = a*(1.0 - ecc*ecc);
		
		// Circular
		if (ecc < small) {
			// ------------- Circular equatorial------------------
			if ((ecc < small) || (fabs(incl - M_PI) < small)) {
				trueLonDot = omegaDot + argPDot + mDot;
				truelon = truelon + trueLonDot * dt;
				truelon = fmod(truelon, twoPi);
			}
			// ------------- Circular inclined    ----------------
			else {
				omega = omega + omegaDot * dt;
				omega = fmod(omega, twoPi);
				argLatDot = argPDot + mDot;
				arglat = arglat + argLatDot * dt;
				arglat = fmod(arglat, twoPi);
			}
		}
		else {
			// -- Elliptical, Parabolic, Hyperbolic equatorial ---
			if ((incl < small) || (fabs(incl - M_PI) < small)) {
				lonPerDot = omegaDot + argPDot;
				lonper = lonper + lonPerDot * dt;
				lonper = fmod(lonper, twoPi);
				m = m + mDot*dt + nDot*dt*dt + nDDot*pow(dt, 3);
				m = fmod(m, twoPi);
				VectorXd e0nu = newtonM(ecc, m);
				nu = e0nu(1);

			}
			// -- Elliptical, Parabolic, Hyperbolic inclined ----
			else {
				omega = omega + omegaDot * dt;
				omega = fmod(omega, twoPi);
				argp = argp + argPDot  * dt;
				argp = fmod(argp, twoPi);
				m = m + mDot*dt + nDot*dt*dt + nDDot*pow(dt, 3);
				m = fmod(m, twoPi);
				VectorXd e0nu = newtonM(ecc, m);
				nu = e0nu(1);
			}
		}

		VectorXd temp(9), result;
		temp << p, ecc, incl, omega, argp, nu, arglat, truelon, lonper;
		
		result = coe2rv(temp);
		return result;
	}

	/**
	 * <summary> Determines the velocity component of the second position vector using Herrick-Gibbs initial orbit determination technique. </summary>
	 * <param name = r1> Position one </param>
	 * <param name = r2> Position two </param>
	 * <param name = r3> Position three </param>
	 * <param name = jd1> Julian Date of the postion one </param>
	 * <param name = jd2> Julian Date of the postion two </param>
	 * <param name = jd3> Julian Date of the postion three </param>
	 * <returns> VectorXd containing the velocity component of the second position. </returns>
	 */
	Vector3d iodHGibbs(Vector3d r1, Vector3d r2, Vector3d r3, const double& jd1, const double& jd2, const double& jd3) {

		// Return velocity vector & other returns
		Vector3d v2;
        double theta, theta1, copa;
		char error[12];

		Vector3d p, r1n;
		double   dt21, dt31, dt32, term1, term2, term3, tolangle;

		/* --------------------  initialize values   -------------------- */
		tolangle = 0.017452406;  // (1.0 deg in rad)

		strcpy(error, "Ok");

		theta = 0.0;
		theta1 = 0.0;

		dt21 = (jd2 - jd1) * 86400.0;
		dt31 = (jd3 - jd1) * 86400.0;   // differences in times
		dt32 = (jd3 - jd2) * 86400.0;

		/* ----------------------------------------------------------------
		*  determine if the vectors are coplanar.
		---------------------------------------------------------------- */
		p = r2.cross(r3); 
		p.normalize();
		r1n = r1;
		r1n.normalize();
		copa = asin(p.dot(r1n));

		if (fabs(p.dot(r1n)) > tolangle)
			strcpy(error, "Not coplanar");

		/* ----------------------------------------------------------------
		* check the size of the angles between the three position vectors.
		*   herrick gibbs only gives "reasonable" answers when the
		*   position vectors are reasonably close.  1.0 deg is only an estimate.
		---------------------------------------------------------------- */
		theta = vAngle(r1, r2);
		theta1 = vAngle(r2, r3);
		if ((theta > tolangle) || (theta1 > tolangle))
			strcpy(error, "Angle > tol");

		double magr1 = r1.norm(), magr2 = r2.norm(), magr3 = r3.norm();

		/* ------------ perform herrick-gibbs method to find v2 --------- */
		term1 = -dt32			* (1.0 / (dt21 * dt31) + MU_E / (12.0 * pow(magr1, 3)));
		term2 = (dt32 - dt21)	* (1.0 / (dt21 * dt32) + MU_E / (12.0 * pow(magr2, 3)));
		term3 = dt21			* (1.0 / (dt32 * dt31) + MU_E / (12.0 * pow(magr3, 3)));
		v2 = r1 * term1 + r2 * term2 + r3 * term3;

		std::cout << error << std::endl;

		return v2;
	}  

    /**
     * <summary> Calculates the covariance transformation between SEZ and TEME Cartesian frames.
     * Linear transformations are not taken into account. </summary>
     * <param name = geo> VectorXd containing geodetic coordinates (latitude, longitude (degrees)
     * and altitude (km)) </param>
     * <param name = date> A date_ms structure contatining the date of the transofrmation. </param>
     * <param name = xp> Polar motion coefficient for x axis (radians) </param>
     * <param name = yp> Polar motion coefficient for x axis (radians) </param>
     * @return MartixXd of the transformation
     */
    MatrixXd getSEZToTEMECovTfMat(VectorXd geo, date_ms date, const double& xp, const double& yp, const double& dim) {

        double jd = getJulianDay(date);
        double gmst = getGMST(jd);
        double latRad = geo(0) / 180.0 * M_PI;
        double lonRad = geo(1) / 180.0 * M_PI;

        MatrixXd tf;

        if (dim == 6) {

            tf = MatrixXd::Zero(6,6);

            MatrixXd rotGd3 = rotZ(lonRad) * rotY(M_PI / 2.0 - latRad);
            MatrixXd rotGd6 = MatrixXd::Zero(6,6);
            MatrixXd tempZero3 = MatrixXd::Zero(3,3);
            MatrixXd rotPM3 = polarMotion(xp, yp);
            MatrixXd rotPM6 = MatrixXd::Zero(6,6);
            MatrixXd rotGMST3 = rotZ(gmst);
            MatrixXd rotGMST6 = MatrixXd::Zero(6,6);

            // Rotation - geodetic 6x6
            rotGd6 << rotGd3, tempZero3, tempZero3, rotGd3;
            // Rotation - polar motion 6x6
            rotPM6 << rotPM3, tempZero3, tempZero3, rotPM3;
            // Rotation - GMST
            rotGMST6 << rotGMST3, tempZero3, tempZero3, rotGMST3;
            tf = rotGMST6 * rotPM6 * rotGd6;

        } else if (dim == 3){

            tf = MatrixXd::Zero(3,3);
            MatrixXd rotGd3 = rotZ(lonRad) * rotY(M_PI / 2.0 - latRad);
            MatrixXd rotPM3 = polarMotion(xp, yp);
            MatrixXd rotGMST3 = rotZ(gmst);

            tf = rotGMST3 * rotPM3 * rotGd3;

        } else {

            tf = MatrixXd::Zero(6,6);
        }


        return tf;
    }

    /**
     * <summary> Calculates the covariance transformation between TEME and SEZ Cartesian frames.
     * Linear transformations are not taken into account. </summary>
     * Can also be obtained by the inverting the matrix from getSEZToTEMECovTfMat function.
     * <param name = geo> VectorXd containing geodetic coordinates (latitude, longitude (degrees)
     * and altitude (km)) </param>
     * <param name = date> A date_ms structure contatining the date of the transofrmation. </param>
     * <param name = xp> Polar motion coefficient for x axis (radians) </param>
     * <param name = yp> Polar motion coefficient for x axis (radians) </param>
     * @return MartixXd of the transformation
     */
    MatrixXd getTEMEToSEZCovTfMat(VectorXd geo, date_ms date, const double& xp, const double& yp, const double& dim) {

        double jd = getJulianDay(date);
        double gmst = getGMST(jd);
        double latRad = geo(0) / 180.0 * M_PI;
        double lonRad = geo(1) / 180.0 * M_PI;

        MatrixXd tf;

        if (dim == 6) {

            tf = MatrixXd::Zero(6,6);

            MatrixXd rotGd3 = rotY( - M_PI / 2.0 + latRad) * rotZ(-lonRad);
            MatrixXd rotGd6 = MatrixXd::Zero(6,6);
            MatrixXd tempZero3 = MatrixXd::Zero(3,3);
            MatrixXd rotPM3 = polarMotion(xp, yp).inverse();
            MatrixXd rotPM6 = MatrixXd::Zero(6,6);
            MatrixXd rotGMST3 = rotZ(-gmst);
            MatrixXd rotGMST6 = MatrixXd::Zero(6,6);

            // Rotation - geodetic 6x6
            rotGd6 << rotGd3, tempZero3, tempZero3, rotGd3;
            // Rotation - polar motion 6x6
            rotPM6 << rotPM3, tempZero3, tempZero3, rotPM3;
            // Rotation - GMST
            rotGMST6 << rotGMST3, tempZero3, tempZero3, rotGMST3;
            tf = rotGd6 * rotPM6 * rotGMST6;

        } else if (dim == 3) {

            tf = MatrixXd::Zero(3,3);

            MatrixXd rotGd3 = rotY( - M_PI / 2.0 + latRad) * rotZ(-lonRad);
            MatrixXd rotPM3 = polarMotion(xp, yp).inverse();
            MatrixXd rotGMST3 = rotZ(-gmst);

            tf = rotGd3 *rotPM3 * rotGMST3;

        } else {

            tf = MatrixXd::Zero(6,6);
        }

        return tf;
    }

    /**
     * <summary> Converts SEZ to TEME coordinates.</summary>
     * <param name = sez> A VectorXd containing SEZ coordinates (km) </param>
     * <param name = geo> VectorXd containing geodetic coordinates (latitude, longitude (degrees)
     * and altitude (km)) </param>
     * <param name = date> A date_ms structure contatining the date of the transofrmation. </param>
     * <param name = xp> Polar motion coefficient for x axis (radians) </param>
     * <param name = yp> Polar motion coefficient for x axis (radians) </param>
     * <param name = type> Type of the geodetic transformation (not used) </param>
     * @return SEZ to TEME transformation matrix
     */
    MatrixXd SEZToTEME(VectorXd sez, VectorXd geo, date_ms date, const double& xp, const double& yp, const int& type) {
        double jd = getJulianDay(date);
        VectorXd ecef = SEZToECEF(sez, geo, type);
        VectorXd teme = ECEFToTEME(ecef, jd, 0, xp, yp);
        return teme;
    }

    /**
     * <summary> Converts TEME to SEZ coordinates.</summary>
     * <param name = teme> A VectorXd containing TEME coordinates (km) </param>
     * <param name = geo> VectorXd containing geodetic coordinates (latitude, longitude (degrees)
     * and altitude (km)) </param>
     * <param name = date> A date_ms structure contatining the date of the transofrmation. </param>
     * <param name = xp> Polar motion coefficient for x axis (radians) </param>
     * <param name = yp> Polar motion coefficient for x axis (radians) </param>
     * <param name = type> Type of the geodetic transformation (not used) </param>
     * @return SEZ to TEME transformation matrix
     */
    MatrixXd TEMEToSEZ(VectorXd teme, VectorXd geo, date_ms date, const double& xp, const double& yp, const int& type) {
        double jd = getJulianDay(date);
        VectorXd ecef = TEMEToECEF(teme, jd, 0, xp, yp);
        VectorXd sez = ECEFToSEZ(ecef, geo, type);
        return sez;
    }

    /**
     * <summary> Converts SEZ to ECEF coordinates. </summary>
     * <param name = sez> A VectorXd containing SEZ coordinates (km) </param>
     * <param name = geo> VectorXd containing geodetic coordinates (latitude, longitude (degrees)
     * and altitude (km)) </param>
     * <param name = type> Type of the geodetic transformation (not used) </param>
     * @return SEZ to ECEF transformation matrix
     */
    MatrixXd SEZToECEF(VectorXd sez, VectorXd geo, const int& type) {
        VectorXd ecef(6), tempR(3), tempV(3), tempS(6);
        double latRad = geo(0) / 180.0 * M_PI;
        double lonRad = geo(1) / 180.0 * M_PI;
        MatrixXd rot(3,3);
        rot = rotZ(lonRad) * rotY(M_PI / 2.0 - latRad);

        tempR << sez(0), sez(1), sez(2);
        tempV << sez(3), sez(4), sez(5);
        ecef << rot * tempR, rot * tempV;
        ecef += GdToECEF(geo, type);

        return ecef;
    }

    /**
     * <summary> Converts ECEF to SEZ coordinates. </summary>
     * <param name = ecef> A VectorXd containing ECEF coordinates (km) </param>
     * <param name = geo> VectorXd containing geodetic coordinates (latitude, longitude (degrees)
     * and altitude (km)) </param>
     * <param name = type> Type of the geodetic transformation (not used) </param>
     * @return ECEF to SEZ transformation matrix
     */
    MatrixXd ECEFToSEZ(VectorXd ecef, VectorXd geo, const int& type) {
        VectorXd tempR(3), tempV(3), temp(6);
        MatrixXd rot(3, 3);

        double latRad = geo(0) / 180.0 * M_PI;
        double lonRad = geo(1) / 180.0 * M_PI;

        temp = ecef - GdToECEF(geo, type);
        rot = rotY(-(M_PI / 2.0 - latRad)) * rotZ(-lonRad);

        tempR << temp(0), temp(1), temp(2);
        tempV << temp(3), temp(4), temp(5);
        temp << rot * tempR, rot * tempV;

        return temp;
    }

    /**
     * Temporary function for calcluating difference between timestamps (hours, minutes, seconds)
     * @brief diffTimeHMS
     * @param start
     * @param end
     */
    double diffTimeHMS(date_ms start, date_ms end) {
        return (end.hour - start.hour) * 3600 + (end.min - start.min) * 60 + end.sec - start.sec;
    }
}

double getHellingerDistance(VectorXd v1, MatrixXd p1, VectorXd v2, MatrixXd p2) {

	// TODO: LLT, inverse
    VectorXd vDiff = v1 - v2;
    MatrixXd pSum = p1 + p2;
    MatrixXd pProd = p1 * p2;
    double epsilon = (-0.25 * vDiff.transpose() * pSum.inverse() * vDiff)(0,0);
    return 1 - sqrt(sqrt(pProd.determinant())/ (0.5 * pSum).determinant()) * exp(epsilon);


}
