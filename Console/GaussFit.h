//$T indentinput.h GC 1.140 01/23/11 19:14:03

//$6
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//

#ifndef GaussFitH
#define GaussFitH

#include "common.h"
#include "TPeakDetect.h"
#include "sortedClasses.h"
#include "TBinaryFileWriter.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// VT
const int		MAXGAUSSIANS = 10;  //maximum number of gaussian fittable
// VT
const double	GAUSS_SIGMA_INIT = 0.5; //atlagos csucsszelesseg
const double	GAUSS_ALFA_INIT = 0.2; //atlagos csucstorzulas


// maximum diff between fit on all points and on top fitSigma
const double	GAUSS_ERROR_SIGMA = 5; //how to calc fit error

const double	MIN_REL_PEAK_HIGHT = 0.025; //don't fit if max int lower than this percentage

enum			GAUSS_PARAM { GAUSS_PARAM_MAX = 0, GAUSS_PARAM_LOC = 1, GAUSS_PARAM_SIGMA = 2, GAUSS_PARAM_EXTRA = 3 };
enum			FITFUNCTYPE { GAUSSFIT = 0, SKEWEDGAUSS = 1, DGAUSS = 2};
enum			FITTYPE {MAXFIT = 0, AREAFIT = 1};
enum			FITALGO { LMMIN_FIT = 0, ALGFIT=1};
const string	FitFunc[4] = {"Gauss","Skewed Gauss","DGauss"}; // DGauss - 2 gauss functions

const double INTEGRATION_PARAM = sqrt(2*M_PI) ; //multiply the area with this to get masslynx like areas SQRT(2*PI)

//const double SG_PARAM_EXTRA_MIN = 0.2; //skewness for Gauss, if too small numerical errors
//const double SG_PARAM_EXTRA_MAX = 0.7; //skewness for Gauss max

const double MAX_BASE_LINE = 0.50; //baseline max height

const int	GAUSS_PARAM_NR = 4; //number of parameters the biggest function has

static double minFitLoc = 0; //min time for gaussian fit, used as constraint by fitting func
static double maxFitLoc = 1000; //max time for gaussian fit, used as constraint by fitting func

enum FIT_RESULTS_MEANING //meaning of fitted func parameters extra
{
	frTail,
	frFront
};




double CalcGaussExtern(double t, const double *p);
double CalcEMGaussExtern(double t, const double *p);

struct TGParam
{
	TGParam()
		:
		maxi(0),
		loc(0),
		sigma(0),
		extra(0),
		fit(0),
		qual(0)
	{ }
	double	maxi;
	double	loc;
	double	sigma;
	double	extra; //[GAUSS_PARAM_NR-3]
	double	fit;
	double	qual;

	void WriteToFile(TBinaryFileWriter *write);

	void ReadFromFile(TBinaryFileWriter *read);

	double	CalcArea(FITFUNCTYPE fitting_function_idx); //calculate area;

	void operator = (const TGParam &input);

};

class	TGaussCalcParam
{
	//
	// -----------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------
	//

public:
	atomic_bool stopped;

	int		maxFitPeaks;	// number of gausses to fit

// chromatographic peak parameters, in minutes
	double	fitSigma;

	double	minHeight;		// min XIC peak
	double	maxHeight;

	FITFUNCTYPE		fitting_function_idx; //which fitting function to use


	// ===================================================================================================================
	// ===================================================================================================================
	//

	void Init()
	{
		maxFitPeaks = 1; // default 1
		fitSigma = 1;			// fit to peak top<>
		minHeight = 100.00;		// min XIC peak
		maxHeight = 1.0e8;
		fitting_function_idx = FITFUNCTYPE::GAUSSFIT;
		stopped = false;
	}

	void WriteParamToList(TStringList* List);
	// ===================================================================================================================
	// ===================================================================================================================
	//
	TGaussCalcParam()	{ Init(); }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	void Assign(const TGaussCalcParam *p);

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	void Get(TGaussCalcParam *p);


	void WriteToFile(TBinaryFileWriter *write);

	void ReadFromFile(TBinaryFileWriter *read);

};

class TGauss :
	public TGaussCalcParam
{
	//
	// -----------------------------------------------------------------------------------------------------------------------
	//    class to fit a simple gaussian function to an intensity array important a*excp(-x2/2c2)
	// -----------------------------------------------------------------------------------------------------------------------
	//
private:

	int 	curr_fit_idx;  //index of the max of the last gaussian in the X array
	int		pointnr;

	double	param[MAXGAUSSIANS * GAUSS_PARAM_NR];
	double	goodparam[MAXGAUSSIANS * GAUSS_PARAM_NR];

	void	FindCurrMax(double *x_in, double *y_in, int start, int end, double *max_out, double *max_pos_out);
	void 	FindCurrMedian(double *x_in, double *y_in, int start, int end, double *max_out, double *max_pos_out);

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//

	void Init()
	{
		for(int i = 0; i < MAXGAUSSIANS * GAUSS_PARAM_NR; i++)
		{
			param[i] = 0;
			//numeric_limits<double>::quiet_NaN();
			goodparam[i] = 0;
			//numeric_limits<double>::quiet_NaN();
		}
		fitnr = 0;
		curr_fit_idx = -1;

		TGaussCalcParam::Init();

		pointnr = ERROR_VALUE;
	}


	// -----------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------
	//
public:
	TGParam	fitData[MAXGAUSSIANS];
	int		fitnr; //number of fitted gaussians

	double	CalcGauss(double t, const double *p) const;
	double	CalcArea(int fitidx); //calculate area for a fit, using the correct function

	double	CalcFitError(double *x_in, double *y_in, int start, int end);

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//

	TGauss()			{ Init(); }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	~		TGauss()	{ }

	//	fitting functions - if no change in algo, they could be called instead of lmfit
		int    FitLM(double* x_in, double* y_in, double* param_in, int pointnr, int paramnr);
		bool   DoFit(double* x_in, double* y_in, int pointnr, int fixed_nr = 0);
};


#endif
