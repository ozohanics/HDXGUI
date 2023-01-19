//$T indentinput.h GC 1.140 01/23/11 19:16:14

//$6
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//

#ifndef RawClassH
#define RawClassH
#include <limits>
#include <stdio.h>
#include "common.h"

#define _NOSQLITE

#include <ntverp.h>
#if VER_PRODUCTBUILD > 9600
// Windows 10+ SDK code goes here
	#pragma comment(lib, "MSToolkit.lib")
#else
// Windows 8.1- SDK code goes here
	#pragma comment(lib, "MSToolkit-81.lib")
#endif


#include <iostream>
#include <iomanip>
#include <cmath>
#include "MSToolkitTypes.h"
#include "MSReader.h"
#include "MSObject.h"
#include "Spectrum.h"

using namespace std;
using namespace MSToolkit;

const int	MAXSPECNR = 3;
const int	MASS_INT_BYTE_NR = 8;
const int	MASS_INT_BYTE_NR_DRE = 12;
const int	IDX_SPEC_BYTE_NR = 22;
const int	IDX_TIME_LOC = 6;
const int	IDX_START_LOC = 0;
const int	IDX_NR_LOC = 2;
const int	IDX_READ_NR = IDX_SPEC_BYTE_NR / sizeof(unsigned short);
const int	FLOATSIZE = 4;
enum ftType { ftCENTROID = 0, ftDRE = 1, ftMZ = 2, ftEnd = 3 }; //open Waters Centroid/DRE or mzML/XML files
const static string ftExtension[] = { "DAT", "DAT", "MZXML, MZML" };
//redundant
const static int ftTypeNr = 3;

const long	MAX_SCAN_NR = 200000;//number of scans
const long  MAX_RAW_SIZE = 2000000000;
const long	MAX_CHR_TIME = MAX_SCAN_NR;
const double SCANRESOLUTION = 0.5 / 60.0; //kb 1s

const long MAX_SUM_POINT_NR = 1000; //maximum nr of points in a summed spectrum

const int MINSPECDATASIZE = 5; //min 5 point in an MS spectrum
//
// =======================================================================================================================
//    storage struct for IDX file data
// =======================================================================================================================
//
class	IDXData
{
	//
	// -----------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------
	//
public:
	double	time;
	long	specStartPos;	// This means that the spectrum file can be MAX 4GB
	long	specEndPos;
	long	peaknr;

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//

	void	Init() { time = 0; specStartPos = 0; specEndPos = 0; peaknr = 0; }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	IDXData() { Init(); }
};

//
// =======================================================================================================================
//    class for IDX file decoding
// =======================================================================================================================
//
class	TIdxDataFile
{
	//
	// -----------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------
	//
private:
	bool	isIDXRead;

	//
	// -----------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------
	//
public:
	long	specnr;		// number of spectra
	IDXData *mzData;	// storage of data
	double	minRT;
	double	maxRT;

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//

	void	Init() { specnr = 0; isIDXRead = false; mzData = nullptr; minRT = 0; maxRT = 0; }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	TIdxDataFile() { Init(); };

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//

	void	Clear() { if (mzData != nullptr) { SAFE_DELETE_ARRAY(mzData); }Init(); }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	~TIdxDataFile() { Clear(); };

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//

	void	AddSpec(long nr)
	{ // init storage array
		SAFE_DELETE(mzData); mzData = new IDXData[nr];
	}

	void	ReadIDXFile(std::string path, ftType filetype);
	double	DecodeSpecTime(unsigned short *data, long nr);
	long	DecodeStart(unsigned short *data, long nr);
	long	DecodeEnd(unsigned short *data, long nr);
	long	DecodePeakNr(unsigned short *data, long nr);
	long	GetScanAtTime(double time); // only usable after IDX file decode finished
};

//
// =======================================================================================================================
//    class for decoding raw DAT data file, 1 spectrum at a time
// =======================================================================================================================
//
class	RawData
{
	//
	// -----------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------
	//
private:
	double	*Fmz;	// mass array
	double	*Fy;	// intensity array

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//

	void	Init()
	{
		nr = 0; scan = 0; RT = 0; Fmz = nullptr; Fy = nullptr;
	}

	//
	// -----------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------
	//
public:
	long	nr;		// number of peaks in spectrum
	double	RT;		// retention time of current spectrum
	long	scan;	// scan nr of current spectrum

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//

	RawData() { Init(); }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	~RawData() { DeletePeaks(); }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	void DeletePeaks()
	{	
		if (Fmz != nullptr)
		{
			//SAFE_DELETE_ARRAY(Fmz);
			delete[] Fmz;
			Fmz = nullptr;
		}

		if (Fy != nullptr)
		{
			//SAFE_DELETE_ARRAY(Fy);
			delete[] Fy;
			Fy = nullptr;
		}

		nr = 0;
	}

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	void AddPeaks(long currnr)
	{
		DeletePeaks();

		if (currnr > 0)
		{
			Fmz = new double[currnr];
			Fy = new double[currnr];
			nr = currnr;

			for (long i = 0; i < currnr; i++)
			{
				Fmz[i] = numeric_limits<double>::quiet_NaN();
				Fy[i] = numeric_limits<double>::quiet_NaN();
				//numeric_limits<double>::quiet_NaN();
			}
		}
	}

	void CopySpectrumData(MSToolkit::Spectrum& sp)
	{
		AddPeaks(sp.size());
		for (long s = 0; s < sp.size(); s++)
		{
			//read all scans
			Fmz[s] = sp.at((int) s).mz;
			Fy[s] = sp.at((int) s).intensity;
		}
	}

	void				DecodeSpectrum(unsigned short *data, long size, ftType filetype);
	double				GetIntensity(unsigned short *data, long offset);
	double				GetMass(unsigned short *data, long offset);

	const double*	mz() { return Fmz; };
	const double*	y() { return Fy; };
};

//
// =======================================================================================================================
//    combined class for reading IDX and DAT file data
// =======================================================================================================================
//
class	TRawDataFile
{
	//
	// -----------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------
	//
private:

	// variable to store whole raw file data !!! This function does not check if the
	// computer has enough memmory !!
	RawData			*filedata;		// used to store everything in memmory, if enough memmory is available
	TIdxDataFile	idx;			// IDX file data
	FILE			*fp;			// file pointer, automatically closed
	unsigned short	*specbuffer;	// pointer to the raw data read from the file
	ftType			filetype;

	string			fpath;
	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	void			Init();

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	
	long			GetNearestGoodScan(long scan); 
	long			SearchRT(double	RT_in, long startScan, long endScan, double limit); //Get scan for RT_in using binary search
	//
	// -----------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	void OpenRawFile(std::string path);

	long scanNr;

public:

	bool isRef; //reference imdicator for HDX

	long GetSpecNr();
	double			GetRTAtScan(long scan);

	MSToolkit::MSReader *msr; //reader for ms mzML/mzXML files
	MSToolkit::Spectrum spec; //will contain one spectrum

	// the next 3 members are not used anymore, as data will only hold 1 spectrum
	// kept for compatibility reasons
	RawData data;	// spectrum storage, limited nr, used initially, when working file-based
	

	atomic_bool stopped;

	long GetFileScanNr(string fileName, MSToolkit::MSSpectrumType format = MSToolkit::MS2);

	bool GetMSMSScan(long idxIn);

	void ReadMGFFile(string fileName, long scan);

	void Clear();
	//
	// ===================================================================================================================
	// ===================================================================================================================
	//

	TRawDataFile() { Init(); }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	~TRawDataFile() { Clear(); }

	
	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	void				CloseFile() { if (fp != 0) { fclose(fp); fp = 0; }if (specbuffer != nullptr) { SAFE_DELETE_ARRAY(specbuffer); } SAFE_DELETE(msr);  specbuffer = nullptr; msr = nullptr; }
	void				RawToMemory(std::string path); //put all contents of data file in memory struct
	void				XMLToMemory(std::string filePath, double startTime = ERROR_VALUE, double endTime = ERROR_VALUE); //put all contents of data file in memory struct, start time and end time are optional
	void				GetSpectrumString(long start, long peaknr, ftType filetype);	// get raw data from
	void				GetSpectrumByScan(long nr); // get one scan
	void				GetSpectrumByTime(double time);  //fill internal data structure with the spectrum

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//

	long				GetScanAtTime(double time); //get scan number for time

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//V688 The 'filetype' function argument possesses the same name as one of the class members, which can result in a confusion. rawclass.h 423
	void				ReadIDXFile(std::string filename, ftType filetypeIn) { idx.ReadIDXFile(filename, filetypeIn); }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//

	double			GetMaxRT() { return idx.minRT; }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	double			GetMinRT() { return idx.maxRT; }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//

	double GetFileDataMZ(long scannr, long pointnr);

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	double GetFileDataY(long scannr, long pointnr);

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	double	GetFileDataRT(long scannr) { if (scannr < scanNr) { return filedata[scannr].RT; }SOFTWARE_ERROR; return ERROR_VALUE; }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	const double	*GetFileDataMZArray(long scannr) { if (scannr < scanNr) { return filedata[scannr].mz(); }SOFTWARE_ERROR; return nullptr; }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	const double	*GetFileDataYArray(long scannr) { if (scannr < scanNr) { return filedata[scannr].y(); }SOFTWARE_ERROR; return nullptr; }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	long		GetFileDataPeakNr(long scannr) { if (scannr < scanNr) { return filedata[scannr].nr; }SOFTWARE_ERROR; return ERROR_VALUE; }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	void	SetFileType(ftType ftype) { filetype = ftype; }

	// returns number of peaks in generated spectrum, summing up sumNr spectra, and deletening the peaks 
	// that are not present in all spectra
	int GenerateSumSpectra(double *mzOut, double *yOut, long scanIdx, double mass, double massWnd, int sumNr = 3);

	long	BinarySearchClosestSpectrumIndex(long scan_in, double	mass_in, double limit);//get the array idx of mass_in in current spectrum.
};
//needed for Scan RT
long			BinarySearchClosestRT(RawData	*array_in, double	RT_in, unsigned long lo, unsigned long up, double limit);



	

#endif
