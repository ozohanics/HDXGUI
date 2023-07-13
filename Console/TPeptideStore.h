#pragma once
#include "rawclass.h"
#include "TPeptideData.h"
#include "TIsotopeData.h"
#include "GaussFit.h"

const int maxIsotopeNr = 7; //number of isotopes to use for calc of mass center 2*maxIsotopeNr+1
const double MIN_ISOTOPE_ABUNDANCE = 5.0e-2; //smallest isotopic peak to calc with
const double MAX_DYNAMIC_RANGE = 500.0; //max difference between highest and smallest peak in a spectrum
const double ISOTOPE_CENTROID_TOP = 0.5; //use 50% of the xEIC to calculate average isotopic distribution
const double DEUTMASSERRORPPM = 50.0;
const int AVGSPECISOTOPENR = 7;

//structure of the peptide list file
enum e_PEPTIDE_INPUT_LIST
{
	inpProteinCOL, //contains the protein ID
	inpPeptideCOL, //contains the peptide sequence
	inpRTsCOL, //contains the start of the retention time window
	inpRTeCOL, //contains the end of the retention time window
	inpZCOL, //contains the charge
	inpTotalCOL //max number of columns
};
//structure of raw data list file 
enum RAWINFO
{
	riFilePath,
	riDeuTime,
	riState,
	riReplicate,
	riSize //number of elements
};

struct TRawDataInfo
{//information about measured data
	std::string rawPath = "";
	double	labelTime = 0.0;
	std::string state = "";
	int replicate = 0;
	bool isRef = false;
};

struct EIC
{
	vector<double> rt;
	vector<double> y;
	double baseMass = 0.0;
};

// properties of an isotopic envelope for a standard peptide
// it assumes a max isotope number of 100
class TEnvelopeDescriptors
{
private:
	double massExp[100]; //experimental mass
//double massTheo[100]; //theoretical mass
	double ppm[100]; //measured mass error
	double intensity[100]; //measured intensity
	int contiguosRegion[100]; //where the ppm is OK
	double minPPM; 
	int lastID; //last filled descriptor
	int start;
	int end;

public:
	TEnvelopeDescriptors() {
		lastID = 0; minPPM = 0;
		start = 0; end = 0;
		for (int i=0; i<100; i++)
		{
			massExp[i] = ppm[i] = intensity[i] = contiguosRegion[i] = -1;
		}
	};
	~TEnvelopeDescriptors() {};

	void addDataPoint(double mass_inn, double ppm_in, double intensity_in)
	{
		if (lastID < 100)
		{
			massExp[lastID] = mass_inn;
			ppm[lastID] = ppm_in;
			intensity[lastID] = intensity_in;
			lastID++;
		}
		else
		{
			throw "too many isotopes error";
		}
	}


	double getMinPPMAvg() {

		if (lastID == 0)
		{
			//throw "Initialize the envelope first error";
			return 0.0;
		}
		
		int minIDX[] = { 0, 0, 0 };

		for (int i=1; i<lastID; i++)
		{
			for (int mn = 2; mn >= 0; mn--) {

				if (ppm[i] < ppm[minIDX[mn]])
				{
					minIDX[mn] = i;
				}
			}
		}

		minPPM = (ppm[minIDX[0]] + ppm[minIDX[1]] + ppm[minIDX[2]]) / 3;
		return minPPM; 

	}; //sort ppm and get best 3 average

	int CalcContiguosRegion(double ppm_error)
	{
		for (int i = 0; i < lastID; i++)
		{
			if (ppm[i] < ppm_error)
			{//if good accuracy, say that it is good?
				contiguosRegion[i] = 1;
			}
			else
			{
				contiguosRegion[i] = -1;
			}
		}
		int regLen = 0;
		int newstart = 0; start = 0;
		int newend = 0; end = 0;
		for (int i = 0; i < lastID; i++)
		{
			if (contiguosRegion[i] == 1)
			{
				if (regLen == 0)
				{
					newstart = i;
				}
				regLen++;
				newend = i;
			}
			else
			{
				if (regLen > 0)
				{
					newend = i-1;
				}
				regLen = 0;
			}
			if ((end - start) < (newend - newstart))
			{
				end = newend;
				start = newstart;
			}
			
		}
		return end - start;
	}

	int getPureEnvelopeStart() {
		return start;
		
	
	}; //find the continuous region of low mass error
	int getPureEnvelopeEnd() {
		return end;
	};

};

////////////////////////////////////////////////////
//	structure is as follows:
//peptides [one for each peptide] . composition
//results [one for each peptide]  . ICs [one for each file] . intensity [ one for each isotope or mass]
//each file is a different labeling experiment
////////////////////////////////////////////////////
class TPeptideStore
{	
public:
	TPeptideData*	peptides; //an array of peptides, only the original data, not calculated
	int				pepNr; //number of peptides

	TRawDataInfo*	dataInfo;
	int				fileNr; //number of raw files

	TRawDataFile*	rawData; //pointer to current raw data
	
	TOutput*		results; //storage for results, one for each peptide, contains all calc values

	//TGauss			gFit; //fit isotope data

	EIC*			xEIC; //extracted ion chrom for one mass only

	bool			useXICOverlap; //overlap calc for isotopes

	double			exp_ppm_error;

	int refFileIdx; //which raw file is the reference

	TPeptideStore() {
		results = nullptr; peptides = nullptr; fileNr = 0; rawData = nullptr; pepNr = 0; dataInfo = nullptr;
		refFileIdx = 0; xEIC = nullptr; useXICOverlap = false; exp_ppm_error = 5.0;
	}
	~TPeptideStore() {
		SAFE_DELETE_ARRAY(results); SAFE_DELETE_ARRAY(peptides);  SAFE_DELETE_ARRAY(dataInfo); SAFE_DELETE_ARRAY(xEIC);
	}

	bool DeletePeptide(int idx);

	int			ReadPeptideFile(std::string fileName); //read file containing peptide input data
	bool		CreateResultStorage(TRawDataFile* rawData, int numFiles); //create storage for data, from all the files to be searched
	bool		SetRawInfo(std::string path, double time, int idx, std::string st="", int repl=0);
	//search the isotopes belonging to a peptide
	bool		SearchAPeptide(TRawDataFile* rawData, int idx, int fileIdx, bool isDeuterated); //search a file for peptide at index idx
	//log results
	void		WriteResultsToStdOut();

	//helper functions to check if the right peaks are used
	double		GetPeptideRT(int pepIdx, int fileIdx, int isoIdx) { return results[pepIdx].ICs[fileIdx].trueRT[isoIdx]; };
	
	//returns mass of the highest peak
	double		CreateAverageIsoPattern(int pepIdx, int fileIdx, int sScan, int eScan, int maxPos, int maxD, double* X, double** Y);
	//create and extracted chromatogram and store results in an internal container xEIC
	void        CreateXIC(int pepIdx, int fileIdx, double mass, double sDelta, TRawDataFile* currRawData);
	
	double		CalcPeptideMassCenter(int pepIdx, int fileIdx); //calculate mass center for an isotopic distribution
	double		CalcReferenceMassCenter(int pepIdx); //calculate mass center for 0 labeling time 
	double		CalcDeuteration(int pepIdx, int fileIdx, int refIdx=ERROR_VALUE); //calculate delta D compared to refIdx
	void		WriteDeuterationDataToStdOut();
	void		WriteDeuterationPercentToStdOut();
	void		WriteDeuterationCentersToStdOut();
	void		WriteDeuterationInMEMHDX();

	int			CalcMaxPepDeuteration(int pepIdx, int minIsoToConsider, bool isDeuterated);
	int 		CheckIsotopeOverlap(int pepIdx, int start, int end, int maxPepD, double** tempY, double* RTIDX, int refIdx);
};

