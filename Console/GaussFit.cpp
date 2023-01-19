//$T indentinput.cpp GC 1.140 01/23/11 19:13:56

//$6
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//

#include <limits>
#include <stdlib.h>
#include <cmath>

#include "lmcurve.h"
#include "GaussFit.h"


double CalcGaussExtern(double t, const double *c)
{
	//if (pGauss != nullptr) { return pGauss->CalcGauss(t, c); }
	double tmp = (-1) * (t - c[GAUSS_PARAM_LOC]) * (t - c[GAUSS_PARAM_LOC]) / (2 * c[GAUSS_PARAM_SIGMA] * c[GAUSS_PARAM_SIGMA]);
	return c[GAUSS_PARAM_MAX] * exp(tmp);
}

double CalcDoubleGaussExtern(double t, const double* c)
{
	//if (pGauss != nullptr) { return pGauss->CalcGauss(t, c); }
	double tmp = (-1) * (t - c[GAUSS_PARAM_LOC]) * (t - c[GAUSS_PARAM_LOC]) / (2 * c[GAUSS_PARAM_SIGMA] * c[GAUSS_PARAM_SIGMA]);
	return c[GAUSS_PARAM_MAX] * exp(tmp);
}


//skewed gauss function
double CalcEMGaussExtern(double t, const double *p)
{//4 param Exponentially modified gaussian
 // H = max
 // L = loc
 // W = sigma
 // S = skewness
 // enum GAUSS_PARAM { GAUSS_PARAM_MAX = 0, GAUSS_PARAM_LOC = 1, GAUSS_PARAM_SIGMA = 2, GAUSS_PARAM_EXTRA = 3 };
	if (p[GAUSS_PARAM_SIGMA] < 0.01)
	{
		return ZERO_CONST.dt;
	}

	if (p[GAUSS_PARAM_EXTRA] < 0.01)
	{
		return ZERO_CONST.dt;
	}

	double S = p[GAUSS_PARAM_EXTRA];
	double W = fabs(p[GAUSS_PARAM_SIGMA])*(4 - S*S); //modositott szelesseg
	if (t < p[GAUSS_PARAM_LOC] - W / 2 / S)
	{
		return 0;
	}

	double value = p[GAUSS_PARAM_MAX];

	double logparam = (2 * S*(t - p[GAUSS_PARAM_LOC]) / W);

	double exponent = (4 / S / S - 1) * (log(1 + logparam) - logparam);


	value = value*exp(exponent);

	return value;
}



// =======================================================================================================================
//	Wrapper for the lm-fit optimization to be able to select which function to use for fitting
// =======================================================================================================================
double TGauss::CalcGauss(double t, const double *p) const
{
	switch (fitting_function_idx)
	{
	case GAUSSFIT:
		return CalcGaussExtern(t, p);
	case SKEWEDGAUSS:
		return CalcEMGaussExtern(t, p);
	case DGAUSS:
		return CalcDoubleGaussExtern(t, p);
	default:
		SOFTWARE_ERROR;
	}
}
//Calculate the area under the curve
double	TGauss::CalcArea(int fitidx)
{
	if (fitidx < 0 || fitidx > fitnr - 1)
	{
		SOFTWARE_ERROR;
	}
	return fitData[fitidx].CalcArea(fitting_function_idx);
}

// =======================================================================================================================
//  Calculate the R2 fit qualifier for data in x_in, y_in using a gaussian with parameters contained in param
//  start and end are the indexes defining a domain for which to calc fit
// =======================================================================================================================
double TGauss::CalcFitError(double *x_in, double *y_in, int start, int end)
{	// calculate R2 fit parameter
	if (start < 0 || end < start || end == start)
	{
		SOFTWARE_ERROR;
	}
	//V688 The 'param' function argument possesses the same name as one of the class members, which can result in a confusion. gaussfit.cpp 215
		//~~~~~~~~~~~~~~
		// calculate average Y
	double	y_avg = 0;
	//~~~~~~~~~~~~~~

	for (int x = start; x < end; x++)
	{
		y_avg = y_avg + y_in[x];
	}

	y_avg = y_avg / (end - start);

	double	SSerr = 0;
	//~~~~~~~~~~~~~~

	for (int i = start; i < end; i++)
	{
		//170309 - V403 - Calc error for the whole EIC and store results
		//if ((x_in[i] > paramVal[GAUSS_PARAM_LOC] - paramVal[GAUSS_PARAM_SIGMA] * GAUSS_ERROR_SIGMA) && (x_in[i] < paramVal[GAUSS_PARAM_LOC] + paramVal[GAUSS_PARAM_SIGMA] * GAUSS_ERROR_SIGMA))
			//ne szamolja a hibat ahol az illesztett gauss 0
		double	value = 0;
		for (int f=0; f<fitnr; f++) //calc all fits
		{
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			value += CalcGauss(x_in[i], &fitData[f].maxi);
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		}

		SSerr = SSerr + (y_in[i] - value) * (y_in[i] - value);
	}

	//~~~~~~~~~~~~~~
	// calculate SStot
	double	SStot = 0;
	//~~~~~~~~~~~~~~

	for (int i = start; i < end; i++)
	{
		SStot = SStot + (y_in[i] - y_avg) * (y_in[i] - y_avg);
	}

	if (SStot < 1e-15)
	{
		return 0.0;
	}
	if (SSerr < 1e-15)
	{
		return 0.0;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double	error = 1 - SSerr / SStot;	// formula for R2
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	return error;
}

// =======================================================================================================================
// find maximum in y_in, substracting a prevoiusly fitted gaussian
// return value and maximum index in max_out and max_pos_out
// =======================================================================================================================
void TGauss::FindCurrMax(double *x_in, double *y_in, int start, int end, double *max_out, double *max_pos_out)
{
	if (start < 0 || (end - start < 1))
	{
		SOFTWARE_ERROR;
	}

	//set starting values
	*max_pos_out = start;
	*max_out = y_in[start];

	double maxI = y_in[start];
	for (int i = start; i < end; i++)
	{
		if (maxI < y_in[i])
		{
			maxI = y_in[i];
			*max_pos_out = x_in[i];
			*max_out = maxI;
			curr_fit_idx = i;
		}
	}
}

// =======================================================================================================================
// find median in y_in, substracting a prevoiusly fitted gaussian
// return value and maximum index in max_out and max_pos_out
// =======================================================================================================================
void TGauss::FindCurrMedian(double *x_in, double *y_in, int start, int end, double *max_out, double *max_pos_out)
{
	if (start < 0 || (end - start < 1))
	{
		SOFTWARE_ERROR;
	}
	//~~~~~~~~~~
	double	currmax = y_in[start];
	//set starting values
	*max_pos_out = start;
	*max_out = currmax;

	//calc median for all points
	double sumIX = 0;
	double sumI = 0;
	double maxI = y_in[start];
	for (int i = start; i < end; i++)
	{
		sumIX += y_in[i] * x_in[i];
		sumI += y_in[i];
		if (maxI < y_in[i])
		{
			maxI = y_in[i];
			curr_fit_idx = i;
		}
	}

	//position of geometric center
	*max_pos_out = sumIX / sumI;
	// really 1/3 from bottom and 2/3 from top
	*max_out = maxI;
}

//write parameters to a text file
void TGaussCalcParam::WriteParamToList(TStringList* List)
{
	if (List == nullptr)
	{
		SOFTWARE_ERROR;
	}
	List->Add("TGaussCalcParam Data");
	List->Add("maxFitPeaks\tfitSigma\tminHeight\tmaxHeight\tfitting_function_idx");
	List->Add(FormatDouble(maxFitPeaks, 0) + std::string("\t") +
		FormatDouble(fitSigma, 2) + std::string("\t") +
		FormatDouble(minHeight, 2) + std::string("\t") +
		FormatDouble(maxHeight, 2) + std::string("\t") +
		FitFunc[fitting_function_idx]);
	List->Add("");
}

void TGaussCalcParam::Assign(const TGaussCalcParam *p)
{
	maxFitPeaks = p->maxFitPeaks;
	fitSigma = p->fitSigma;
	minHeight = p->minHeight;
	maxHeight = p->maxHeight;
	fitting_function_idx = p->fitting_function_idx;
}

void TGaussCalcParam::Get(TGaussCalcParam *p)
{
	p->maxFitPeaks = maxFitPeaks;
	p->fitSigma = fitSigma;
	p->minHeight = minHeight;
	p->maxHeight = maxHeight;
	p->fitting_function_idx = fitting_function_idx;
}

//use binary serialization for param
void TGaussCalcParam::WriteToFile(TBinaryFileWriter *write)
{
	if (write != nullptr)
	{
		//write->Write(mass);
		write->WriteValue(maxFitPeaks);	// number of gausses to fit

		// chromatographic peak parameters, in minutes
		write->WriteValue(fitSigma);

		write->WriteValue(minHeight);		// min XIC peak
		write->WriteValue(maxHeight);

		write->WriteValue(fitting_function_idx); //which fitting function to use

		char spacebytes[PARAM_FILE_PADDING_LENGTH];

		write->WriteValue(spacebytes, PARAM_FILE_PADDING_LENGTH);
	}
}

//fill from binary datafile
void TGaussCalcParam::ReadFromFile(TBinaryFileWriter *read)
{
	if (read != nullptr)
	{
		//max = read->ReadDouble();;

		//number of gausses to fit
		read->ReadValue(&maxFitPeaks);

		// chromatographic peak parameters, in minutes
		read->ReadValue(&fitSigma);
		read->ReadValue(&minHeight);		// min XIC peak
		read->ReadValue(&maxHeight);

		int tmp;
		read->ReadValue(&tmp); //which fitting function to use
		fitting_function_idx = (FITFUNCTYPE)tmp;

		char spacebytes[PARAM_FILE_PADDING_LENGTH + 1];
		read->ReadValue(spacebytes, PARAM_FILE_PADDING_LENGTH);
	}
}

//calc area of several funtions
double TGParam::CalcArea(FITFUNCTYPE fitting_function_idx) //calculate area
{
	// ReSharper disable once CppInitializedValueIsAlwaysRewritten
	double area = 0;
	if (fitting_function_idx == SKEWEDGAUSS)
	{
		double S = extra;
	}

	double gp;
	switch (fitting_function_idx)
	{
	case GAUSSFIT:
		//return CalcOrigGauss(t, p);
		area = (maxi)* sigma * INTEGRATION_PARAM;
		break;
	case SKEWEDGAUSS:
		//area = boost::math::tgamma(gp)*max*extra*sigma / 2 / exp((gp - 1)*(log(gp - 1) - 1));
		gp = 4 / extra / extra;  //gamma func param
		area = tgamma(gp)*maxi*extra*sigma / 2 / exp((gp - 1)*(log(gp - 1) - 1));
		break;
	case DGAUSS:
		area = maxi* (sigma+extra) * INTEGRATION_PARAM/2;
		break;
	default:
		SOFTWARE_ERROR;
	};
	return fabs(area);
}

void TGParam::operator=(const TGParam &input)
{
	maxi = input.maxi;
	loc = input.loc;
	sigma = fabs(input.sigma);
	extra = input.extra;
	fit = input.fit;
	qual = input.qual;
}

//read from binary files
void TGParam::ReadFromFile(TBinaryFileWriter *read)
{
	if (read != nullptr)
	{
		read->ReadValue(&maxi);
		read->ReadValue(&loc);
		read->ReadValue(&sigma);
		read->ReadValue(&extra);
		read->ReadValue(&fit);
		read->ReadValue(&qual);
	}
}

void TGParam::WriteToFile(TBinaryFileWriter *write)
{
	if (write != nullptr)
	{
		write->WriteValue(maxi);
		write->WriteValue(loc);
		write->WriteValue(sigma);
		write->WriteValue(extra);
		write->WriteValue(fit);
		write->WriteValue(qual);
	}
}




// fitting functions using LMFit or ALGLib
int	TGauss::FitLM(double* x_in, double* y_in, double* param_in, int pointnr, int paramnr)
{

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	lm_status_struct   status;
	lm_control_struct   control = lm_control_double;
	//control.stepbound = 1000;
	control.maxcall = 10000;
	control.printflags = 0;         // monitor status (+1) and parameters (+2)
									//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
									// * perform the fit ;

//set global parameters for fit constraint
	minFitLoc = x_in[0];
	maxFitLoc = x_in[pointnr - 1];
//start fitting
	switch (fitting_function_idx)
	{
	case GAUSSFIT:
		lmcurve_fit(paramnr, param_in, pointnr, x_in, y_in, &CalcGaussExtern, &control, &status);
		break;
	case SKEWEDGAUSS:
		lmcurve_fit(paramnr, param_in, pointnr, x_in, y_in, &CalcEMGaussExtern, &control, &status);
		break;
	case DGAUSS:
		lmcurve_fit(paramnr, param_in, pointnr, x_in, y_in, &CalcDoubleGaussExtern, &control, &status);
		break;
	default:
		//SOFTWARE_ERROR;
		WRITE_DEBUG("Unknown fittig function called: ");
		WRITE_DEBUG(to_string(fitting_function_idx));
		break;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	WRITE_DEBUG(lm_infmsg[status.info]);
	return status.info; //information about how the fit ended
}

bool TGauss::DoFit(double* x_in, double* y_in, int pointnr, int fixed_nr)
{

	if (x_in == nullptr || y_in == nullptr)
	{
		return false;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~
	int padding = 7;
	int n_par = GAUSS_PARAM_NR; // number of parameters in model function f
	int m_dat = pointnr + padding;	// number of data pairs

	double* tempX, * tempY;

	

	tempX = new double[m_dat];
	tempY = new double[m_dat];

	for (int a = 0; a <= pointnr; a++)
	{
		if (a==0) 
		{
			tempY[a] = 0.0;
			tempX[a] = x_in[0]-0.1;
		}
		else
		{
			tempY[a] = y_in[a - 1];
			tempX[a] = x_in[a - 1];
		}
		
	}

	for (int a = pointnr+1; a < m_dat; a++)
	{
		
		tempY[a] = 0.0;
		tempX[a] = tempX[a - 1] + 0.1;
	}
	
	
	for (int a = 0; a < m_dat; a++)
	{
		WRITE_DEBUG(FormatDouble(tempX[a], 4) + "\t" + FormatDouble(tempY[a], 4));
	}
	

	param[GAUSS_PARAM_MAX] = y_in[0];
	param[GAUSS_PARAM_LOC] = x_in[0];

	// initialize sigma with default value (1)
	param[GAUSS_PARAM_SIGMA] = GAUSS_SIGMA_INIT;	
	if (fitnr > 0) param[GAUSS_PARAM_SIGMA] = fitData[fitnr - 1].sigma;
	param[GAUSS_PARAM_EXTRA] = GAUSS_ALFA_INIT;
	if (fitnr > 0) param[GAUSS_PARAM_EXTRA] = fitData[fitnr - 1].extra;
	//Try fitting to the top of the peak


	FitLM(tempX, tempY, param, m_dat, n_par - fixed_nr);

	// eredmenyek tarolasa
	fitData[fitnr].maxi = param[GAUSS_PARAM_MAX];
	fitData[fitnr].sigma = param[GAUSS_PARAM_SIGMA];
	fitData[fitnr].loc = param[GAUSS_PARAM_LOC];
	fitData[fitnr].extra = param[GAUSS_PARAM_EXTRA];
	fitData[fitnr].fit = CalcFitError(x_in, y_in, 0, pointnr - 1);	// status.fnorm;

	//std::cout << fitnr << std::endl;
	WRITE_DEBUG(FormatDouble(fitData[fitnr].loc, 2) + "\t" + FormatDouble(fitData[fitnr].maxi, 2) + "\t" + FormatDouble(fitData[fitnr].sigma, 2) + "\t" + FormatDouble(fitData[fitnr].extra, 2) + "\t" + FormatDouble(fitData[fitnr].fit, 2) + "\t");
	fitnr++;
	
	
	delete [] tempX;
	delete [] tempY;
	return true;
}

