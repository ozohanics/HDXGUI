//$T indentinput.h GC 1.140 01/23/11 19:16:42

//$6
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//

#ifndef savgolH
#define savgolH
class	SavGol
{	// Translation from VB6
	///;
	///'Savitzky_Golay Smoothing 'If SmoothCurrent is set to true, then the last smoothed data set will be smoothed,
	///'(i.e. The new smoothing will be cumulative over the last smoothing operation. If 'SmoothCurrent is false, then
	///the original Y-data will be smoothed, and the smoothedY array will 'be overwritten 'The Aavitzky-Golay smoothing
	///algorithm essentialy fits the data to a second order polynomial 'within a moving data window. It assumes that
	///the data has a fixed spacing in the x direction, 'but does work even if this is not the case. 'For more info
	///see: '"Smoothing and Differentiation of Data by Simplified Least Squares Procedure", 'Abraham Savitzky and
	///Marcel J. E. Golay, Analytical Chemistry, Vol. 36, No. 8, Page 1627 (1964) 'Degree 2 = 5 point 'Degree 3 = 7
	///point ...etc
	///;
	///needed input data degree of smoothing output data
	///;
	///'Dynamic data arrays

	//
	// -----------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------
	//
public:
	double	*DataYPtr;
	int		peaknr;

	//
	// * The matrix for the Savitzky-Golay Coefficents These are filled in the form load
	// * event
	//
	int		SGCoef[13][14];

	//
	// ===================================================================================================================
	//    Dim SGCoef(1 To 11, 0 To 13) As Integer
	// ===================================================================================================================
	//
	SavGol()
	{
		peaknr = 0;
		DataYPtr = 0;

		//
		// * 'Set the Smoothing Coefficients for Savitzky-Golay 'The zeroth value is the
		// * normalization factor
		//

		int i = 2;

		SGCoef[i][1] = 17;
		SGCoef[i][2] = 12;
		SGCoef[i][3] = -3;
		SGCoef[i][0] = 35;
		i++;

		SGCoef[i][1] = 7;
		SGCoef[i][2] = 6;
		SGCoef[i][3] = 3;
		SGCoef[i][4] = -2;
		SGCoef[i][0] = 21;
		i++;

		SGCoef[i][1] = 59;
		SGCoef[i][2] = 54;
		SGCoef[i][3] = 39;
		SGCoef[i][4] = 14;
		SGCoef[i][5] = -21;
		SGCoef[i][0] = 231;
		i++;

		SGCoef[i][1] = 89;
		SGCoef[i][2] = 84;
		SGCoef[i][3] = 69;
		SGCoef[i][4] = 44;
		SGCoef[i][5] = 9;
		SGCoef[i][6] = -36;
		SGCoef[i][0] = 429;
		i++;

		SGCoef[i][1] = 25;
		SGCoef[i][2] = 24;
		SGCoef[i][3] = 21;
		SGCoef[i][4] = 16;
		SGCoef[i][5] = 9;
		SGCoef[i][6] = 0;
		SGCoef[i][7] = -11;
		SGCoef[i][0] = 143;
		i++;

		SGCoef[i][1] = 167;
		SGCoef[i][2] = 162;
		SGCoef[i][3] = 147;
		SGCoef[i][4] = 122;
		SGCoef[i][5] = 87;
		SGCoef[i][6] = 42;
		SGCoef[i][7] = -13;
		SGCoef[i][8] = -78;
		SGCoef[i][0] = 1105;
		i++;

		SGCoef[i][1] = 43;
		SGCoef[i][2] = 42;
		SGCoef[i][3] = 39;
		SGCoef[i][4] = 34;
		SGCoef[i][5] = 27;
		SGCoef[i][6] = 18;
		SGCoef[i][7] = 7;
		SGCoef[i][8] = -6;
		SGCoef[i][9] = -21;
		SGCoef[i][0] = 323;
		i++;

		SGCoef[i][1] = 269;
		SGCoef[i][2] = 264;
		SGCoef[i][3] = 249;
		SGCoef[i][4] = 224;
		SGCoef[i][5] = 189;
		SGCoef[i][6] = 144;
		SGCoef[i][7] = 89;
		SGCoef[i][8] = 24;
		SGCoef[i][9] = -51;
		SGCoef[i][10] = -136;
		SGCoef[i][0] = 2261;
		i++;

		SGCoef[i][1] = 329;
		SGCoef[i][2] = 324;
		SGCoef[i][3] = 309;
		SGCoef[i][4] = 284;
		SGCoef[i][5] = 249;
		SGCoef[i][6] = 204;
		SGCoef[i][7] = 149;
		SGCoef[i][8] = 84;
		SGCoef[i][9] = 9;
		SGCoef[i][10] = -76;
		SGCoef[i][11] = -171;
		SGCoef[i][0] = 3059;
		i++;

		SGCoef[i][1] = 79;
		SGCoef[i][2] = 78;
		SGCoef[i][3] = 75;
		SGCoef[i][4] = 70;
		SGCoef[i][5] = 63;
		SGCoef[i][6] = 54;
		SGCoef[i][7] = 43;
		SGCoef[i][8] = 30;
		SGCoef[i][9] = 15;
		SGCoef[i][10] = -2;
		SGCoef[i][11] = -21;
		SGCoef[i][12] = -42;
		SGCoef[i][0] = 806;
		i++;

		SGCoef[i][1] = 467;
		SGCoef[i][2] = 462;
		SGCoef[i][3] = 447;
		SGCoef[i][4] = 422;
		SGCoef[i][5] = 387;
		SGCoef[i][6] = 322;
		SGCoef[i][7] = 287;
		SGCoef[i][8] = 222;
		SGCoef[i][9] = 147;
		SGCoef[i][10] = 62;
		SGCoef[i][11] = -33;
		SGCoef[i][12] = -138;
		SGCoef[i][13] = -253;
		SGCoef[i][0] = 5135;
		i++;
	}

	void	DoSmooth(double *y, int pointnr, int Degree, int smoothnr);
	void 	DoSmoothWithOutput(const double *y, double *SmoothedY, int pointnr, int Degree, int smoothnr);
};
#endif
