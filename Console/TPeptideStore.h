#pragma once
#include "rawclass.h"
#include "TPeptideData.h"
#include "TIsotopeData.h"
#include "GaussFit.h"

const int maxIsotopeNr = 7; //number of isotopes to use for calc of mass center 2*maxIsotopeNr+1
const double MIN_ISOTOPE_ABUNDANCE = 5.0e-2; //smallest isotopic peak to calc with
const double MAX_DYNAMIC_RANGE = 100.0; //max difference between highest and smallest peak in a spectrum
const double ISOTOPE_CENTROID_TOP = 0.5; //use 50% of the xEIC to calculate average isotopic distribution
const double DEUTMASSERRORPPM = 50.0;
const int AVGSPECISOTOPENR = 13;

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

	int refFileIdx; //which raw file is the reference

	TPeptideStore() {
		results = nullptr; peptides = nullptr; fileNr = 0; rawData = nullptr; pepNr = 0; dataInfo = nullptr;
		refFileIdx = 0; xEIC = nullptr; useXICOverlap = false;
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
	bool		CheckCoelution(int pepIdx, int fileIdx);
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

