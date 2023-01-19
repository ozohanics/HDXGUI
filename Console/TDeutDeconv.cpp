#include "TDeutDeconv.h"
//#include "minmax.h"
#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

TDeutDeconv::TDeutDeconv(int sizeN)
{
	size = sizeN;

	aMatrix.setlength(size, size);
	bMatrix.setlength(size);
	solution_X.setlength(size);
	parameters = new double[size];

	for (size_t i = 0; i < size; i++)
	{
		for (size_t j = 0; j < size; j++)
		{
			aMatrix[i][j] = 0.0;
		}
		bMatrix[i] = 0.0;
		solution_X[i] = 0.0;
		parameters[i] = 0.0;
	}
	info = -1; //error
	parameters[1] = 0.5;
	parameters[2] = 0.5;
}

void TDeutDeconv::FillArrays(double* isotopeComposition, double* deuteratedPeaks, int isoNr, int maxD)
{// use outside info for matrices
	//xy.setlength(iArrayWidth, iArrayHeight);
	for (size_t h = 0; h < size; h++)
	{//cycle for height
		bMatrix[h] = deuteratedPeaks[h];
		std::cout << "B matrix " << h << "\t" << bMatrix[h] << std::endl;
	}

	//fill up the space for normal isotopes, no deut needed
	for (size_t h = 0; h < isoNr; h++)
	{//cycle for height
		//std::cout << "A matrix " << h << "\t";
		for (size_t w = 0; w <=h; w++)
		{//cycle for width
			if (w > maxD) continue;
			if (h>=w)
			{				
				aMatrix[w][h] = isotopeComposition[h-w];
			//	std::cout << aMatrix[w][h] << "\t" ;
			}					
		}
		//std::cout << std::endl;
	}
	
	int isoCount = 0;
	for (size_t h = isoNr; h <isoNr+maxD; h++)
	{//cycle for height
		isoCount++;
		for (size_t w = 1; w <= maxD; w++)
		{//cycle for width		
			if (w+isoCount > maxD ) continue;
			aMatrix[w+isoCount][h] = isotopeComposition[isoNr-w];
		}
	}
	/*
	for (size_t h = maxD; h < isoNr + maxD; h++)
	{//cycle for height
		size_t w = maxD;
		{//cycle for width
			aMatrix[w][h] = isotopeComposition[h-maxD];
		}
	}*/

}

void TDeutDeconv::CalcCurrIsotopeSum()
{
	for (size_t h = 0; h < size; h++)
	{//cycle for height
		for (size_t w = 0; w < size; w++)
		{//cycle for width
			solution_X[h] += parameters[w] * aMatrix[w][h];			
		}
		std::cout << solution_X[h] << "\t";
		std::cout << std::endl;
	}
}

void TDeutDeconv::CalculateDeutIncorp()
{//use ALGLIB to calculate 

	
	std::cout << info << std::endl;
	for (size_t i = 0; i < size; i++)
	{
		std::cout << solution_X[i] << std::endl;
	}

}

void TDeutDeconv::WriteArrayToStdOut()
{
	for (size_t h = 0; h < size; h++)
	{//cycle for height
		for (size_t w = 0; w < size; w++)
		{//cycle for width
			std::cout << aMatrix[w][h] << "\t";
		}
		std::cout << std::endl;
	}
}

