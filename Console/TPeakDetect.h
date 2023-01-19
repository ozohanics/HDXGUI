#pragma once

#include "common.h"
#include "savgol.h"
//Class for chromatographic peak detection
const int MAX_IDX_PEAKS = 100;

class TPeakDetect : public SavGol
{
private:
	//storage objects
	double* mz; //x scale
	double* DataY; //y scale

	double* fderiv; //first derivative
	double* f2deriv; //first derivative
	int* fpeak; //peak found
	int		peakIdx[MAX_IDX_PEAKS]; //location of peaks, max 100


	int     peaknr;

	int pointNr;
	TStringList fileIO;

public:
	TPeakDetect();
	~TPeakDetect();

	double amplitudeThr; //amplitude threshold
	double slopeThr; //slope threshold

	void ClearData();
	void ReadInput(string fileName);
	void SetData(double* x, double* y, int nr);
	void WriteData(string fileName);
	void CalculateDeriv();
	void DetectPeaks();

	int GetPeakNr() { return peaknr; };
	int GetPeakIDX(int nr); //get the location of peak nr
	int		peakStart[MAX_IDX_PEAKS]; //start of peaks, max 100
	int		peakEnd[MAX_IDX_PEAKS]; //end of peaks, max 100
	double  peakMaxInt[MAX_IDX_PEAKS];
	double  GetPeakRT(int peakID) { return mz[peakIdx[peakID]]; };

	void SumPeaks(TPeakDetect& tpIn); //create a list of peaks that have smallest common area/length
};
