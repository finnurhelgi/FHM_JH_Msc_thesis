This code was created in connection with our master thesis; System Identification of Active Devices for Frequency Domain Studies.
Written with support from Siemens Gamesa Renewable Energy (SGRE).
As there are some aspects of the thesis, namely the contents of chapter 4 of the thesis, which use data from SGRE,
this code is not meant to replicate every step of the thesis, but rather supplement the explaination of the system identification
process outlined in the report.

Note: 	To run the following scripts, Matlab (v.2022B or later) is required, along with the Matrix Fitting Toolbox,
	which is available from the SINTEF website: https://www.sintef.no/en/software/vector-fitting/

The following is a discription of the included files.

	example_system.m: 		Runs the simple example system discribed in Section 3.4.

	black_box_modeling.m: 		Contains the system identification process of the black-box STATCOM and grid 
					described in chapter 5. Note that user input is required to run the script, including 
					selection of vector fitting method, model order, etc.

	accuracy_evaluation.m:		Benchmarks the accuracy of vector fit models. Used in Section 4.2.4.

	cluster_poles.m:		Used in the PC-CF multiport vector fitting method, to identify pole clusters.
	
	Vf_driver_driver.m:		Function used to ease the use of different multiport vector fitting methods.

	VFdriver_non_sym.m:		Modification of VFdriver.m function found in Matrix Fitting Toolbox, to accept
					non-symmetrical matrices inputs.

	VFdriver_multi_SIMO.m: 		Modification of VFdriver.m function found in Matrix Fitting Toolbox, to use the
					multi-SIMO method, as is described in Section 3.5.2

	VFdriver_PCCF.m:		Modification of VFdriver.m function found in Matrix Fitting Toolbox, to use the PC-CF
					method, as is described in Section 3.5.2.

	Participation_factor_analysis_STATCOM.m:	Calculates participation matrix using system eigenvalues

	Participation_factor_calculation.xlxs:		Calculates results of p-factor analysis for the STATCOM model


Finnur Helgi Malmquist, Jeppe Haugaard

Date: 01/08/2024 (d/m/y)

