#pragma once

#include <string>
#include "TPeptideStore.h"
#include "rawclass.h"



using namespace std;
class HDXMain
{

	

//private:
public:
	int				peptideNr;
	int				fileNr;
	std::string		peptideListPath;
	std::string		rawFileListPath;
	TStringList		fileList;
	TPeptideStore*	peptideStore;
	
	TRawDataFile*	rawData;
	int				last_file_idx;

	HDXMain() {
		peptideNr = -1; fileNr = -1; peptideListPath = ""; rawFileListPath = ""; last_file_idx = -1;
		fileList.clear();
		peptideStore = new TPeptideStore;
		rawData = new TRawDataFile;
	}
	~HDXMain() {
		delete peptideStore;
		delete rawData;
	}
	//delete old state
	void ReInit()
	{
		peptideNr = -1; fileNr = -1; peptideListPath = ""; rawFileListPath = "";
		fileList.clear();
		delete peptideStore;
		delete rawData;
		peptideStore = new TPeptideStore;
		rawData = new TRawDataFile;
	}

	void SetPPMError(double errorVal) { peptideStore->exp_ppm_error = errorVal; }

	//read a tab delimited text file containing peptide input data, with the format below and no header
	// 	   Start-End Seqeuence	RTstart	RTend	Charge
	//35-70	IEKNETLGGTXLNVGXIPSKALLNNSHYYHMAHGKD	9.3622	9.9622	3
	bool ReadPeptideList(std::string fileName) { peptideNr = peptideStore->ReadPeptideFile(fileName); peptideListPath = fileName; return false; }
	//read a file with the list of datafiles and info about state and no header
	//		Path	DeuterationTime	State	Replicate
	//20210707_E3WT_H2O_MS_02AFAMM-C.mzML	0.0	7.3FAD	1
	bool ReadRawFileList(std::string fileName) { fileList.LoadFromFile(fileName); fileNr = fileList.ItemCount(); rawFileListPath = fileName; return false; }
	//Search for peptide data in all files
	bool SearchDataFiles();

	void SaveResults();

	void SaveModifiedInput();

	//removes a peptide from list
	bool DeletePeptide(int idx); 
	//modify peptide results, by checking a different retention time region
	bool ReEvaluatePeptide(int peptideIdx, int fileIdx, double RTs, double RTe, bool useXICOverlap = true, int charge=-1); //check the peptide with modified parameters, with or without 
																														  //checking for isotopes chromatographic peaks
	//modify isotope inclusion in results, by setting the intensity to 0
	bool ResetIsotope(int peptideIdx, int fileIdx, double isoMin, double isoMax, bool useMax = true);
};

