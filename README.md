# Readme #
Implementation of several Random Finite Sets (RFS) filters for processing radar (range, azimuth, elevation) measurements. Data is provided in tracking data message (.tdm) or Network Common Data Form (NetCDF, .nc, .cdf) format.

	[x] Gaussian Mixture Joint Target Detection and Track (Bernoulli) Filter (GM-JoTT)
		[x] beta GM-JoTT Filter (+ probability of detection)
		[ ] robust GM-JoTT Filter (+ clutter rate)
	
	[ ] Gaussian Mixture Probability Hypothesis Density (GM-PHD) Filter - not tested
		[ ] Cardinalized GM-PHD (GM-CPHD) - not fully implemented

# Dependencies #

* Eigen 
* Boost 
* FADBAD 
* dirent.h 