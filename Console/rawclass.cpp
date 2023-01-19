//$T indentinput.cpp GC 1.140 01/23/11 19:16:07

//$6
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
#include "rawclass.h"
//#include "pch.h"
//#include "time.h"
//#include "mem.h"
#define _FILE_OFFSET_BITS 64 

//
// =======================================================================================================================
//    func for binary search, to find the closest value to a given one, return its position
// =======================================================================================================================
//
long BinarySearchClosestRT(RawData	*array_in, double	RT_in, unsigned long lo, unsigned long up, double limit)
{
	if (array_in == nullptr)
	{
		return -1;
	}

	if (fabs(RT_in) < 1e-10)
	{
		return -1;
	}

	//~~~~~~~~~~~~~~~~~
	int upLimit = up - 1;
	int loLimit = lo;
	int middle = (up + lo) / 2;
	//~~~~~~~~~~~~~~~~~

	while (lo <= up)
	{
		middle = (up + lo) / 2;

		if (middle < loLimit)
		{
			break;
		}

		if (middle > upLimit)
		{
			break;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		double	delta = fabs(array_in[middle].RT - RT_in);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		if (delta <= limit)
		{
			break; //return middle
		}

		if (RT_in > array_in[middle].RT)
		{
			lo = middle + 1;
		}
		else
		{
			up = middle - 1;
		}
	}

	double delta1 = fabs(RT_in - array_in[min(loLimit, middle - 1)].RT); //check the distance to the one before it
	double delta2 = fabs(RT_in - array_in[max(upLimit, middle + 1)].RT); //check the distance to the one after it
	double delta3 = fabs(RT_in - array_in[middle].RT);

	if (delta1 <= delta2)
	{
		if (delta1 <= delta3)
		{
			return min(loLimit, middle - 1);
		}
		else
		{
			return middle;
		}
	}
	else
	{
		if (delta2 <= delta3)
		{
			return max(upLimit, middle + 1);
		}
		else
		{
			return middle;
		}
	}

	SOFTWARE_ERROR;
	return ERROR_VALUE;
}



//
// =======================================================================================================================
// =======================================================================================================================
//
long TIdxDataFile::GetScanAtTime(double time)
{
	//~~~~~~~~~~~~~~
	double	limit = 0;
	//~~~~~~~~~~~~~~

	if (!isIDXRead)
	{	// file nor read
		return -1;
	}

	if (specnr < 2)
	{
		return -1;
	}

	if ((time < mzData[0].time) || (time > mzData[specnr - 1].time))
	{	// outside chromatogram
		return -1;
	}

	if ((time > mzData[0].time) && (time < mzData[specnr - 1].time))
	{	// inside chromatogram, specify limit
		limit = (mzData[specnr - 1].time - mzData[0].time) / (specnr + 1);
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	long	lo = 0;
	long	up = specnr - 1;
	long	middle = (lo + up + 1) / 2;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	while (lo <= up)
	{
		//~~~~~~~~~~~~~~~~
		double	delta;
		double	elem = time;
		//~~~~~~~~~~~~~~~~

		try
		{
			delta = fabs(mzData[middle].time - elem);
		}
		catch (...)
		{
			middle = middle - 1;
			delta = fabs(mzData[middle].time - elem);
		}

		if (delta < limit)
		{
			return middle;
		}
		else if (elem > mzData[middle].time)
		{
			lo = middle + 1;
		}
		else
		{
			up = middle - 1;
		}

		middle = (up + lo) / 2;
	}

	return -1;
}

//
// =======================================================================================================================
// =======================================================================================================================
//
double TIdxDataFile::DecodeSpecTime(unsigned short *data, long nr)
{
	// get RT for current spec
	if (data == nullptr)
	{
		SOFTWARE_ERROR;
	}

	//~~~~~~~~~~~~~~
	float	time = -1;
	//~~~~~~~~~~~~~~

	memcpy(&time, data + nr, sizeof(float)); //float needs to be 32 bit

	return time;
}

//
// =======================================================================================================================
// =======================================================================================================================
//
long TIdxDataFile::DecodeStart(unsigned short *data, long nr)
{	// find the starting timepoint of chromatogram
	return(data[nr + 1] * 65536 + data[nr]);
}

//
// =======================================================================================================================
// =======================================================================================================================
//
long TIdxDataFile::DecodePeakNr(unsigned short *data, long nr)
{
	// get the number of peaks or scans
	return(data[nr + 1] & 255) * 65536 + data[nr];
}

//
// =======================================================================================================================
// =======================================================================================================================
//
void TIdxDataFile::ReadIDXFile(std::string path, ftType filetype)
{
	if (filetype > ftDRE || !FileExists(path))
	{
		return; //no IDX file to read
	}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Open file and put contents in mzData;
	FILE			*fp;
	long			nr = IDX_READ_NR;
	unsigned short	line[IDX_SPEC_BYTE_NR];
	size_t			filesize;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	isIDXRead = false;

	fopen_s(&fp, path.c_str(), "rb");

	if (fp != 0)
	{
		// read file get filesize
		fseek(fp, 0L, SEEK_END);
		filesize = ftell(fp);

		// go back to the begining
		fseek(fp, 0L, SEEK_SET);

		if (filesize % IDX_SPEC_BYTE_NR != 0)
		{
			// check if the file size is correct
			return;
		}

		specnr = filesize / IDX_SPEC_BYTE_NR;

		AddSpec(specnr + 1);

		//~~~~~~~~~~~~~~~~~~~~~~~
		long maxtimenr = specnr - 1;
		//~~~~~~~~~~~~~~~~~~~~~~~

		for (long i = 0; i < specnr; i++)
		{
			// fseek ( fp , IDX_SPEC_BYTE_NR , SEEK_CUR);
			fread(line, sizeof(unsigned short), nr, fp);	// read end-start elements of size 2
			mzData[i].time = DecodeSpecTime(line, IDX_TIME_LOC);
			mzData[i].specStartPos = DecodeStart(line, IDX_START_LOC);
			mzData[i].peaknr = DecodePeakNr(line, IDX_NR_LOC);

			switch (filetype)
			{
			case ftDRE:
				mzData[i].specEndPos = mzData[i].specStartPos + mzData[i].peaknr * MASS_INT_BYTE_NR_DRE;
				break;

			case ftCENTROID:
				mzData[i].specEndPos = mzData[i].specStartPos + mzData[i].peaknr * MASS_INT_BYTE_NR;
				break;

			default:
				SOFTWARE_ERROR;
			}

			if (i == 1)
			{
				if (mzData[1].specStartPos != mzData[0].specEndPos)
				{
					WRITE_DEBUG("Incorrect file format specified!");
					SOFTWARE_ERROR;
				}
			}

			if (i == 0)
			{
				minRT = mzData[i].time;
			}

			if (i == maxtimenr)
			{
				maxRT = mzData[i].time;
			}
		}

		isIDXRead = true;
	}

	if (fp) fclose(fp);
}

//
// =======================================================================================================================
// =======================================================================================================================
//
double RawData::GetIntensity(unsigned short *data, long offset)
{
	//~~~~~~~~~~~~~
	// transform bytes to double
	double	y;
	long	exponent;
	//~~~~~~~~~~~~~

	// Most az intenzitasokat dekodoljuk shr 6 = /64 and 63 = modulo 64
	exponent = ((data[offset + 1] >> 6) & 63) - 21;
	y = (((data[offset + 1] & 63) << 16) + data[offset + 0]) * pow(2, (double)exponent);

	return y;
}

//
// =======================================================================================================================
// =======================================================================================================================
//
double RawData::GetMass(unsigned short *data, long offset)
{
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// get mass encoded in data
	double	mz;
	long	exponent = (data[offset + 3] >> 11);
	long	upperByte = data[offset + 3] & 2047;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	exponent = exponent - 27;
	mz = (double(upperByte) * 65536 + double(data[offset + 2])) * pow(double(2), double(exponent));

	return mz;
}

//
// =======================================================================================================================
// =======================================================================================================================
//
void RawData::DecodeSpectrum(unsigned short *data, long size, ftType filetype)
{
	// for the current scan, store mass-int pairs
	AddPeaks(size);

	//~~~~~~~~~~~~~~~~~~~~
	long	fileType_offset;
	//~~~~~~~~~~~~~~~~~~~~

	switch (filetype)
	{
	case ftDRE:			fileType_offset = 6;
		break;
	case ftCENTROID:	fileType_offset = FLOATSIZE;	// 4
		break;
	default:			SOFTWARE_ERROR;
		return;
	}
	try {
		for (long i = 0; i < size; i++)
		{
			// read one mass and one intensity first 4 bytes are intensity
			Fmz[i] = GetMass(data, i * fileType_offset);
			Fy[i] = GetIntensity(data, i * fileType_offset);
			if (Fy[i] < 0)
			{
				SOFTWARE_ERROR;
			}
		}
	}
	catch (...)
	{
		//this will occur if thread terminates in the middle of processing
	}
}

//
// =======================================================================================================================
// =======================================================================================================================
//
void TRawDataFile::GetSpectrumString(long start, long peaknr, ftType filetype)
{
	//~~~~~~~~~~~~
	// read a spectrum from file, between start and end file locations
	long	nr = -1;
	//~~~~~~~~~~~~

	if (fp != nullptr)
	{
		// get data
		fseek(fp, 0L, SEEK_SET);  //seek to the end???
		switch (filetype)
		{
		case ftDRE:
			nr = peaknr * MASS_INT_BYTE_NR_DRE / sizeof(unsigned short);	// (end-start)/MASS_INT_BYTE_NR;
			///
			break;

		case ftCENTROID:
			nr = peaknr * MASS_INT_BYTE_NR / sizeof(unsigned short);		// (end-start)/MASS_INT_BYTE_NR;
			///
			break;

		default:
			SOFTWARE_ERROR;
		}

		if (specbuffer != nullptr)
		{
			SAFE_DELETE_ARRAY(specbuffer);
		}

		specbuffer = new unsigned short[nr + 1];			// store results for fread
		fseek(fp, start, SEEK_SET);
		fread(specbuffer, sizeof(unsigned short), nr, fp);	// read end-start elements of size 1
	}
	else
	{
		WRITE_DEBUG("File is not open! Please open file before reading!");
	}
}

//
// =======================================================================================================================
// =======================================================================================================================
//
void TRawDataFile::GetSpectrumByScan(long nr)
{
	// get spectrum of scan nr;
	if (idx.specnr <= 0)
	{
		// NO IDX FILE ERROR
		return;
	}

	if (nr > idx.specnr - 1)
	{
		return;
	}

	//GetSpectrumString(idx.mzData[nr].specStartPos, idx.mzData[nr].specEndPos, idx.mzData[nr].peaknr, filetype);
	GetSpectrumString(idx.mzData[nr].specStartPos, idx.mzData[nr].peaknr, filetype);
	data.DecodeSpectrum(specbuffer, idx.mzData[nr].peaknr, filetype);
	data.scan = nr;
	data.RT = idx.mzData[nr].time;
}

//
// =======================================================================================================================
// =======================================================================================================================
//
void TRawDataFile::GetSpectrumByTime(double time)
{
	//~~~~~~~~~~~~~~~~~
	// get the closest spectrum to time
	double	t1, t2;
	double	diff1, diff2;
	long	index = -1;
	//~~~~~~~~~~~~~~~~~

	if (idx.specnr <= 0)
	{
		// NO IDX FILE ERROR
		return;
	}

	for (long i = 0; i < idx.specnr - 1; i++)
	{
		t1 = idx.mzData[i].time;
		t2 = idx.mzData[i + 1].time;
		if (t2 > time)
		{
			diff2 = t2 - time;
			diff1 = time - t1;
			if (diff2 <= diff1)
			{
				// closer to t2
				index = i + 1;
				break;
			}
			else
			{
				index = i;
				break;
			}
		}
	}

	if (index == -1)
	{
		return; // time not found error
	}

	//GetSpectrumString(idx.mzData[index].specStartPos, idx.mzData[index].specEndPos, idx.mzData[index].peaknr, filetype);
	GetSpectrumString(idx.mzData[index].specStartPos, idx.mzData[index].peaknr, filetype);
	data.DecodeSpectrum(specbuffer, index, filetype);
	data.scan = index;
	data.RT = idx.mzData[index].time;
}

//
// =======================================================================================================================
// =======================================================================================================================
//
long TRawDataFile::GetScanAtTime(double time_in)
{
	if (filetype != ftMZ)
	{
		return idx.GetScanAtTime(time_in);
	}


	long scanAtRT = ERROR_VALUE;
	if (scanNr > 0 && scanNr < MAX_CHR_TIME)
	{
		scanAtRT = BinarySearchClosestRT(filedata, time_in, 0, scanNr, 0.01);
	}

	return scanAtRT;
}

//
// =======================================================================================================================
// =======================================================================================================================
//
void TRawDataFile::XMLToMemory(std::string filePath, double startTime, double endTime)
{
	//filePath = ReplaceSubStr(filePath, "/", "\\");
	
	if (!FileExists(filePath)) return;
	if (filetype <= ftDRE) return;
	try
	{
		SAFE_DELETE_ARRAY(filedata);
		filedata = nullptr;

		//use MSReader
		OpenRawFile(filePath);
		//Clear();
		

		if (scanNr < 1)
		{
			return;
		}

		long startScan = 1;
		long endScan = scanNr;

		if (startTime > 0) startScan = SearchRT(startTime, startScan, endScan, SCANRESOLUTION);
		if (endTime > 0) endScan = SearchRT(endTime, startScan, endScan, SCANRESOLUTION);

		//std::cout << endScan << "__" << startScan << std::endl;
		//return;

		if (endScan > startScan)
		{
			
			filedata = new RawData[endScan - startScan]; //3 scan / second?
		}
		else
		{
			WRITE_DEBUG("WRONG RT Window or bad MZXML file!!!");
			SOFTWARE_ERROR;
		}

		size_t currScan = 0;

		//spec.clear();
		//msr->readFile(nullptr, spec, 0);//try to read first spectrum.

		for (long i = startScan; i < endScan; i++)
		{
			if (stopped) { return; }
			spec.clear();

			msr->readFile(nullptr, spec, i);//try to read first spectrum.
			if (spec.size() > MINSPECDATASIZE) //needs to have some actual data, not just one noise
			{
				filedata[currScan].CopySpectrumData(spec);
				filedata[currScan].scan = i;
				filedata[currScan].RT = spec.getRTime();
#ifdef LOG_DEBUG
				//LOG read data
				//WRITE_INFO(FormatDouble(filedata[currScan].RT, 4));
				/*for (int n=0; n < filedata[currScan].nr; n++ )
				{
					if (filedata[currScan].mz()[n] > 1200 && filedata[currScan].mz()[n] < 1250)
					{
						WRITE_INFO(FormatDouble(filedata[currScan].mz()[n], 4) + "\t" + FormatDouble(filedata[currScan].y()[n], 4));
					}

				}*/
				//WRITE_INFO("\n");
#endif // LOG_DEBUG
				currScan++;
			}
		}

		scanNr = currScan;
	}
	catch (...)
	{
		// say that you couldnt open the file
		//Dialogs::ShowMessage(u"ERROR!!! Could not open the RAW file!");
		WRITE_DEBUG("ERROR!!! Could not read MZXML into memory!");
	}
}

long TRawDataFile::GetFileScanNr(string fileName, MSToolkit::MSSpectrumType format /*= MSToolkit::MS2*/)
{
	SOFTWARE_ERROR;
	scanNr = ERROR_VALUE;
	fpath = fileName;
	msr->setFilter(format);
	msr->addFilter(MSToolkit::MS1);
	msr->addFilter(MSToolkit::MS2);
	msr->checkFileFormat(fpath.c_str());
	//msr->createIndex();
	msr->readFile(fpath.c_str(), spec);

	return msr->getLastScan();

	//return scanNr;
}

bool TRawDataFile::GetMSMSScan(long idxIn)
{
	//V688 The 'idx' function argument possesses the same name as one of the class members, which can result in a confusion. rawclass.h 328
	spec.clear();
	if (idxIn > ERROR_VALUE && idxIn < scanNr)
	{
		msr->setFilter(MSToolkit::MS2);
		msr->readFile(nullptr, spec, idxIn);
		if (spec.size() > 0) return  true;
	}
	return false;
}

void TRawDataFile::ReadMGFFile(string fileName, long scan)
{
	//ToDo: implement function
	//WRITE_DEBUG(fileName); //to prevent an unused variable warning
	//WRITE_DEBUG(to_string(scan)); //to prevent an unused variable warning
	//WRITE_DEBUG("Unimplemented ReadMGFFile");
	//start implementing 171006
	try
	{
		long scanIDX = 0;
		bool mgfOK = false;
		if (scan == -1)
		{//load a new file
			if (!FileExists(fileName)) return;
			CloseFile();

			msr = new MSToolkit::MSReader;
			msr->addFilter(MSToolkit::MS2);
			//use MSReader
			spec.clear();
			msr->checkFileFormat(fileName.c_str());
			//msr->createIndex();
			msr->setHighResMGF(true);
			msr->setOnePlusMGF(true);


			spec.clear();
			mgfOK = msr->readMGFFile(fileName.c_str(), spec);//try to read first spectrum
			if (mgfOK)	scanIDX++;
		}
		else
		{
			//set spec to the desired spectrum, file already open
			//while (mgfOK)
			{
				spec.clear();
				mgfOK = msr->readMGFFile("", spec, true);//try to read first spectrum
				if (mgfOK)	scanIDX++;
				if (scanIDX == scan+1 && mgfOK)
				{
					return;
				}
			}
		}
		
	}
	catch (...)
	{
		WRITE_DEBUG("MGF file read error");
		return;
	}

}

void TRawDataFile::Clear()
{
	CloseFile();
	SAFE_DELETE_ARRAY(filedata);
	SAFE_DELETE_ARRAY(specbuffer);
	spec.clear();
	scanNr = ERROR_VALUE;
	SAFE_DELETE(msr);
	filedata = nullptr;
	specbuffer = nullptr;
	msr = nullptr;
}

void TRawDataFile::OpenRawFile(std::string path)
{
	if (!FileExists(path)) return;
	try
	{
		CloseFile();

		msr = new MSToolkit::MSReader;
		msr->addFilter(MSToolkit::MS1);

		fpath = path;
#if defined (_MSC_VER) //for microsoft compiler
		fopen_s(&fp, path.c_str(), "rb"); // open file for read in binary mode
#else
		fp = fopen(path.c_str(), "rb"); // open file for read in binary mode
#endif

		if (filetype > ftDRE)
		{
			//use MSReader
			spec.clear();
			msr->checkFileFormat(path.c_str());
			//msr->createIndex();
			msr->readFile(fpath.c_str(), spec);//try to read first spectrum
			scanNr = msr->getLastScan();

			spec.clear();
			if (scanNr < 1)
			{
				WRITE_DEBUG("Wrong file format! " + path);
			}
		}
	}
	catch (...)
	{
		// say that you couldnt open the file
		//Dialogs::ShowMessage(u"ERROR!!! Could not open the RAW file!");
		fclose(fp);
		WRITE_DEBUG("ERROR!!! Could not open the RAW file!");
	}
}

//
// =======================================================================================================================
// =======================================================================================================================
//
void TRawDataFile::RawToMemory(std::string path)
{
	// read all file to filedata
	//Clear();
	SAFE_DELETE_ARRAY(filedata);
	filedata = nullptr;

	OpenRawFile(path);

	//in measurement file to big, one should not put it in RAM, limit 2GB on 32bits
	if (idx.specnr > MAX_SCAN_NR)
	{
		SOFTWARE_ERROR;
	}

	filedata = new RawData[idx.specnr];

	for (long i = 0; i < idx.specnr - 1; i++)
	{
		if (stopped) { return; }
		GetSpectrumString(idx.mzData[i].specStartPos, idx.mzData[i].peaknr, filetype);
		filedata[i].DecodeSpectrum(specbuffer, idx.mzData[i].peaknr, filetype);
		filedata[i].scan = i;
		filedata[i].RT = idx.mzData[i].time;
	}

	scanNr = GetSpecNr();

	if (fp != 0)
	{
		fclose(fp);
		fp = 0;
	}
}

double TRawDataFile::GetFileDataMZ(long scannr, long pointnr)
{
	if (scannr < scanNr)
	{
		if (filedata[scannr].nr < 1)
		{
			return -1;
		}

		if (pointnr < filedata[scannr].nr)
		{
			return filedata[scannr].mz()[pointnr];
		}
	}
	else
	{
		SOFTWARE_ERROR;
	}

	return ERROR_VALUE;
}

double TRawDataFile::GetFileDataY(long scannr, long pointnr)
{
	if (scannr < scanNr)
	{
		if (filedata[scannr].nr < 1)
		{
			return -1;
		}

		if (pointnr < filedata[scannr].nr)
		{
			return filedata[scannr].y()[pointnr];
		}
	}
	else
	{
		SOFTWARE_ERROR;
	}
	return ERROR_VALUE;
}

long TRawDataFile::SearchRT(double RT_in, long startScan, long endScan, double limit)
{
	if (!FileExists(fpath)) return ERROR_VALUE;
	if (filetype != ftMZ) return ERROR_VALUE;

	if (fabs(RT_in) <= SCANRESOLUTION) //nothing in the first second
	{
		return ERROR_VALUE;
	}

	//~~~~~~~~~~~~~~~~~
	long up = endScan;
	long lo = startScan;
	long middle = ERROR_VALUE;
	//~~~~~~~~~~~~~~~~~

	while (lo <= up)
	{
		middle = (up + lo) / 2;

		if (middle < startScan)
		{
			break;
		}

		if (middle > endScan)
		{
			break;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		double	delta = fabs(GetRTAtScan(middle) - RT_in);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		if (delta <= limit)
		{
			break; //return middle
		}

		if (RT_in > GetRTAtScan(middle))
		{
			lo = middle + 1;
		}
		else
		{
			up = middle - 1;
		}
	}

	return GetNearestGoodScan(middle);
}

//checks if the datapoints from one scan are present in the neighbouring ones. If not, sets their intensity to 0;
int TRawDataFile::GenerateSumSpectra(double *mzOut, double *yOut, long scanIdx, double mass, double massWnd, int sumNr)
{
	long windowStart = ERROR_VALUE, windowEnd = ERROR_VALUE;
	long massIdx = SearchPPM(filedata[scanIdx].mz(), mass, 0, filedata[scanIdx].nr, 15); //5ppm ket scan kozott
	if (massIdx < 0) return ERROR_VALUE;

	for (int i = massIdx - 1; i > 0; i--)
	{
		if (filedata[scanIdx].mz()[i] < mass - massWnd / 2)
		{
			windowStart = i;
			break;
		}
	}
	for (int i = massIdx + 1; i < filedata[scanIdx].nr; i++)
	{
		if (filedata[scanIdx].mz()[i] > mass + massWnd / 2)
		{
			windowEnd = i;
			break;
		}
	}

	if (windowEnd < 0  || windowStart < 0) return ERROR_VALUE;

	if (scanIdx < sumNr / 2 + 1 || scanIdx > scanNr - sumNr / 2 - 1) return ERROR_VALUE;
	try
	{
		int count = 0;
		for (int p = windowStart; p < windowEnd; p++)
		{
			if (p-windowStart >= MAX_SUM_POINT_NR) break;
			mzOut[count] = filedata[scanIdx].mz()[p];
			yOut[count] = filedata[scanIdx].y()[p];
			count++;
		}

		return count;
	}
	catch (...)
	{
		return ERROR_VALUE;
	}
	return ERROR_VALUE;
}

long TRawDataFile::BinarySearchClosestSpectrumIndex(long scan_in, double mass_in, double limit)
{//find mz in scan
	auto array_in = GetFileDataMZArray(scan_in);
	unsigned long lo = 1;
	unsigned long up = this->GetFileDataPeakNr(scan_in) - 1;

	if (array_in == nullptr)
	{
		return ERROR_VALUE;
	}

	if (fabs(mass_in) < 100)
	{
		return ERROR_VALUE;
	}

	//~~~~~~~~~~~~~~~~~
	int upLimit = up - 1;
	int loLimit = lo;
	int middle = (up + lo) / 2;
	//~~~~~~~~~~~~~~~~~

	while (lo <= up)
	{
		middle = (up + lo) / 2;

		if (middle < loLimit)
		{
			break;
		}

		if (middle > upLimit)
		{
			break;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		double	delta = fabs(array_in[middle] - mass_in);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		if (delta <= limit)
		{
			break; //return middle
		}

		if (mass_in > array_in[middle])
		{
			lo = middle + 1;
		}
		else
		{
			up = middle - 1;
		}
	}

	double delta1 = fabs(mass_in - array_in[min(loLimit, middle - 1)]); //check the distance to the one before it
	double delta2 = fabs(mass_in - array_in[max(upLimit, middle + 1)]); //check the distance to the one after it
	double delta3 = fabs(mass_in - array_in[middle]);

	if (delta1 <= delta2)
	{
		if (delta1 <= delta3)
		{
			return min(loLimit, middle - 1);
		}
		else
		{
			return middle;
		}
	}
	else
	{
		if (delta2 <= delta3)
		{
			return max(upLimit, middle + 1);
		}
		else
		{
			return middle;
		}
	}

	SOFTWARE_ERROR;
	return ERROR_VALUE;
}

long TRawDataFile::GetSpecNr()
{
	if (filetype <= ftDRE) return idx.specnr; //if Waters, use IDX
	return scanNr;
}

void TRawDataFile::Init()
{
	specbuffer = nullptr;
	fp = nullptr;
	filetype = ftMZ;
	filedata = nullptr;

	spec.clear();
	fpath.clear();
	scanNr = 0;
	stopped = false;

	msr = nullptr;
	isRef = false;
}

double TRawDataFile::GetRTAtScan(long scan)
{
	if (!FileExists(fpath)) return ERROR_VALUE;
	if (filetype != ftMZ) return ERROR_VALUE;

	try
	{
		if (scan < 1 || scan > scanNr) return ERROR_VALUE;

		for (int sc = scan; sc > 0; sc--)
		{
			spec.clear();
			msr->readFile(nullptr, spec, sc);//try to read last spectrum.
			if (spec.size() > 0)
			{
				return spec.getRTime();
			}
		}
	}
	catch (...)
	{
		return ERROR_VALUE;
	}
	return ERROR_VALUE;
}

long TRawDataFile::GetNearestGoodScan(long scan)
{
	if (!FileExists(fpath)) return ERROR_VALUE;
	if (filetype != ftMZ) return ERROR_VALUE;

	try
	{
		//use MSReader

		if (scan < 1 || scan > scanNr) return ERROR_VALUE;

		for (int sc = scan; sc > 0; sc--)
		{
			spec.clear();
			msr->readFile(nullptr, spec, sc);//try to read last spectrum.
			if (spec.size() > 0)
			{
				return sc;
			}
		}
	}
	catch (...)
	{
		return ERROR_VALUE;
	}
	return ERROR_VALUE;
}
