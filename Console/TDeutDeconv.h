#pragma once

#include "solvers.h"
///////////////////////////////////////////////////////////////////////
// OUTPUT PARAMETERS
// Info - return code:
// *-3    matrix is very badly conditioned or exactly singular.
// * -1    N <= 0 was passed
// * 1    task is solved(but matrix A may be ill - conditioned,
// 	check R1 / RInf parameters for condition numbers).
// 	Rep - additional report, following fields are set :
// *rep.r1    condition number in 1 - norm
// * rep.rinf  condition number in inf - norm
// X - array[N], it contains :
// *info > 0 = > solution
// * info = -3 = > filled by zeros
// void alglib::rmatrixsolve(
//		real_2d_array a,
//		ae_int_t n,
//		real_1d_array b,
//		ae_int_t& info,
//		densesolverreport& rep,
//		real_1d_array& x,
//		const xparams _params = alglib::xdefault);
///////////////////////////////////////////////////////////////////////
class TDeutDeconv
{
	public:
		alglib::real_2d_array		aMatrix; //0..N-1, 0..N-1 and contains the parameters
		alglib::real_1d_array		bMatrix; //0..N-1, will contain the deuterated peaks
		alglib::ae_int_t			size; //N
		alglib::ae_int_t			info; //return code for alglib
		alglib::densesolverreport   report;
		alglib::real_1d_array		solution_X;
		double* parameters;

		TDeutDeconv(int sizeN);
		~TDeutDeconv()
		{
			//delete arrays

		}

		void FillArrays(double* isotopes, double* deuteration, int isoNr, int maxD);
		void CalcCurrIsotopeSum();
		double CalcFitDelta() 
		{
			double delta = 0.0;
			for (size_t i = 0; i < size; i++)
			{
				delta += pow(solution_X[i] - bMatrix[i], 2);
			}
			return delta;
		}
		void CalculateDeutIncorp();
		void WriteArrayToStdOut();

};

