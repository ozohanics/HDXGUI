#pragma once

#include <limits>
#include <time.h>
#include <algorithm>
#include <map>

#include "common.h"
#include "sortedClasses.h"
#include "TBinaryFileWriter.h"



const double		THEO_MASS_ERROR = 0.0005; //mikor egyenlo ket tomega beirtak kozul
const int			MAX_ISOTOPE_SHIFT = 100; //maximum nr of deuterium exchanged
const double		BASELINE = 150.0; //baseline for mass spectra
const int			MAX_COMBINED_SPECTRA_NR = 11;

struct TBinaryFileVersion
{
	TBinaryFileVersion()
		:
		bvMajor(1),
		bvMinor(0)
	{}
	int bvMajor;
	int bvMinor;
};
const TBinaryFileVersion CURRFILEVER;

struct Point2D
{
	double x;
	double y[MAX_COMBINED_SPECTRA_NR + 1]; //last is used as summed spec
	Point2D() { x = 0; for (int i = 0; i <= MAX_COMBINED_SPECTRA_NR; i++) y[i] = 0; }
};

struct TDeutData
{
	double massCenter = 0.0; //mass center
	double deutDelta = 0.0; //mass difference
	double deutTime = 0.0;  //time of labeling
	TDeutData operator = (TDeutData& in)
	{
		TDeutData result;
		result.deutDelta = in.deutDelta;
		result.deutTime = in.deutTime;
		result.massCenter = in.massCenter;
		return result;
	}
};

class TMRC_Iso
{//class to sum up subsequent scans with Maximum Ratio Combining / Goldsmith A: Wireless communications. Cambridge Univ Pr; 2005.
public:
	TMRC_Iso() { lastChargeIdx = ERROR_VALUE; }
	~TMRC_Iso() { Reset(); }

	void AddPoint(double xVal, double yVal, int charge); //add a new point, with conversion to 1 charge

	void CalcMRC(); //calc an averaged spectrum, where all points need to be in at least max-1 places

	void Reset() {
		spec.clear(); lastChargeIdx = ERROR_VALUE; for (int i = 0; i <= MAX_COMBINED_SPECTRA_NR; i++) aPoint.y[i] = 0;
	}

	void Write_Log();

	void SetSpace(int nr) { spec.reserve(nr); }
private:
	Point2D aPoint;
	int lastChargeIdx;
	vector<Point2D> spec; //holds the data of the spectra
};

// -----------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------
class	TExIC	// extracted IC
{
public:
	//order of properties is the same as for the binary file
	double		mass;  //mass of the compound, include in file
	std::string	name;  //name of the compound, include in file
	int			charge;
	int			isotopeNr; //number of added isotopes
	int			maxIsotopes;

	double* intensity;  //intensity for all charge states, include in file
	double* MZ;        //mass for all points
	double* trueRT;	   //peak max RT for each isotope

	vector<double> maxScanInt;
	vector<int> maxScanLoc;

	// write the results and the main data to binary file
	void WriteToFile(TBinaryFileWriter* write);

	//read results from a binary file
	void ReadFromFile(TBinaryFileWriter* read);		
	
	void Copy(TExIC* inPut);

	// ===================================================================================================================
	// ===================================================================================================================
	void Init();

	// ===================================================================================================================
	// ===================================================================================================================
	TExIC() { Init(); }

	// ===================================================================================================================
	// ===================================================================================================================
	~TExIC();

	void WriteToStdOut();
};

// =======================================================================================================================
//    Class to store XIC chromatograms
// =======================================================================================================================
class	TOutput
{
private:
	//variables for iterator
public:
	TExIC*		ICs;
	int			fileXICNr;		// max nr of ICs, number of files

	int			maxIsotopes; // max nr of points per ICs
	TDeutData*  deutResults;
	TDeutData   refMassData; //reference mass center
	TDeutData   isoRefData; //reference mass center, based on isotopes

	void    ZeroIC(std::string name, int icsIndex);

	void WriteToFile(TBinaryFileWriter* write);

	//read output from file
	void ReadFromFile(TBinaryFileWriter* read);

	void Copy(TOutput* inPut);

	TOutput(TOutput* in) //copy constructor
	{
		ICs = nullptr; Reset();
		Copy(in);
		deutResults = in->deutResults;  //toDo: correct it with a copy operator
		refMassData = in->refMassData;
		isoRefData = in->isoRefData;
	}
	TOutput() {
		ICs = nullptr; Reset(); deutResults = nullptr; 
	}
	~TOutput() { Reset(); }

	void			AddICs(int nr, int scannr);		
	void			AddPointToExIC(double mass, double y, double ppm, int icIdx, int isotopeIdx);
	void			SetExICName(std::string name, int idx) { ICs[idx].name = name; }

	void			Reset(); // Delete all data;
	//remove values;					
	void			ReInit();

	void			WriteToStdOut();

	void			CalculateRefMassData();

};

class TIsotopeData
{
public:
	//it is just a storage for isotopic patterns

};

