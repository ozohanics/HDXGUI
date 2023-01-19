//$T indentinput.cpp GC 1.140 01/23/11 19:16:37

//$6
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
#include <limits>
#include "math.h"

#include "savgol.h"
#include "common.h"

//
// =======================================================================================================================
// =======================================================================================================================
//
//  Regi simitas, output nelkul
void SavGol::DoSmooth(double *y, int pointnr, int Degree, int smoothnr)
{
	if (pointnr < 5) return;
	if (Degree < 2)
	{
		return;
	}
	if (Degree > 12)
	{
		//Dialogs::ShowMessage(L"Max Smooth Window is 12!");
		WRITE_DEBUG("Max Smooth Window is 12!");
		Degree = 12;
	}
	//~~~~~~~~~~~~
	// Dim i As Integer, J As Integer Dim TempSum As Double On Error Resume Next
	double	TempSum;
	//~~~~~~~~~~~~

	DataYPtr = y;

	// ReDim Temp(1 To NP) As Double
	peaknr = pointnr;

	//~~~~~~~~~~~~~~~
	double	*SmoothedY;
	//~~~~~~~~~~~~~~~

	SmoothedY = new double[peaknr];

	for (int i = 0; i < peaknr; i++)
	{
		SmoothedY[i] = 0;
		//numeric_limits<double>::quiet_NaN();
	}

	// 'we cannot smooth too close to the data bounds
	for (int i = 0; i < Degree; i++)
	{
		SmoothedY[i] = DataYPtr[i];
	}
	// a veget sem simitod
	for (int i = peaknr - Degree; i < peaknr; i++)
	{
		SmoothedY[i] = DataYPtr[i];
	}

	for (int i = Degree; i < (peaknr - Degree); i++)
	{
		TempSum = DataYPtr[i] * SGCoef[Degree][1];

		for (int j = 1; j <= Degree; j++)
		{
			TempSum = TempSum + DataYPtr[i - j] * (SGCoef[Degree][j + 1]);
			TempSum = TempSum + DataYPtr[i + j] * (SGCoef[Degree][j + 1]);
		}

		SmoothedY[i] = TempSum / SGCoef[Degree][0];
	}

	// 'The last smoothed data will be used to create a new smoothed data set,
	// 'therefore the smoothing operations will be additive
	for (int s = 1; s < smoothnr; s++)
	{
		for (int i = Degree; i < (peaknr - Degree); i++)
		{
			TempSum = SmoothedY[i] * SGCoef[Degree][1];

			for (int j = 1; j <= Degree; j++)
			{
				TempSum = TempSum + SmoothedY[i - j] * (SGCoef[Degree][j + 1]);
				TempSum = TempSum + SmoothedY[i + j] * (SGCoef[Degree][j + 1]);
			}

			SmoothedY[i] = TempSum / SGCoef[Degree][0];
		}
	}

	for (int i = Degree; i < (peaknr - Degree); i++)
	{
		DataYPtr[i] = SmoothedY[i];
	}

	delete[] SmoothedY;
}

//
// =======================================================================================================================
// =======================================================================================================================
//
void SavGol::DoSmoothWithOutput(const double *aDataY, double *SmoothedY, int pointnr, int Degree, int smoothnr)
{
	if (Degree < 2)
	{
		return;
	}
	if (Degree > 12)
	{
		//Dialogs::ShowMessage(L"Max Smooth Window is 12!");
		WRITE_DEBUG("Max Smooth Window is 12!");
		Degree = 12;
	}
	//~~~~~~~~~~~~
	// Dim i As Integer, J As Integer Dim TempSum As Double On Error Resume Next
	double	TempSum;
	//~~~~~~~~~~~~

	// ReDim Temp(1 To NP) As Double
	peaknr = pointnr;

	if (SmoothedY == nullptr)
	{
		SOFTWARE_ERROR;
	}

	for (int i = 0; i < peaknr; i++)
	{
		SmoothedY[i] = 0;
		//numeric_limits<double>::quiet_NaN();
	}

	// 'we cannot smooth too close to the data bounds
	for (int i = 0; i < Degree; i++)
	{
		SmoothedY[i] = aDataY[i];
	}
	// a veget sem simitod
	for (int i = peaknr - Degree; i < peaknr; i++)
	{
		SmoothedY[i] = aDataY[i];
	}

	for (int i = Degree; i < (peaknr - Degree); i++)
	{
		TempSum = aDataY[i] * SGCoef[Degree][1];

		for (int j = 1; j <= Degree; j++)
		{
			TempSum = TempSum + aDataY[i - j] * (SGCoef[Degree][j + 1]);
			TempSum = TempSum + aDataY[i + j] * (SGCoef[Degree][j + 1]);
		}

		SmoothedY[i] = TempSum / SGCoef[Degree][0];
	}

	// 'The last smoothed data will be used to create a new smoothed data set,
	// 'therefore the smoothing operations will be additive
	for (int s = 1; s < smoothnr; s++)
	{
		for (int i = Degree; i < (peaknr - Degree); i++)
		{
			TempSum = SmoothedY[i] * SGCoef[Degree][1];

			for (int j = 1; j <= Degree; j++)
			{
				TempSum = TempSum + SmoothedY[i - j] * (SGCoef[Degree][j + 1]);
				TempSum = TempSum + SmoothedY[i + j] * (SGCoef[Degree][j + 1]);
			}

			SmoothedY[i] = TempSum / SGCoef[Degree][0];
		}
	}
}