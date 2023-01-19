#include "TIsotopeData.h"



void TMRC_Iso::AddPoint(double xVal, double yVal, int charge)
{//create a new point
	double new_x = (xVal - masses::Hp) * charge; //0 charge
	double delta = new_x * PPM_DIFF / 1.0e6; //use 15 ppm for binning
	int idx = ERROR_VALUE;
	for (int i = 0; i < spec.size(); ++i)
	{//search if x exists
		if (fabs(spec[i].x - new_x < delta))
		{
			idx = i;
			break;
		}
	}
	if (idx == ERROR_VALUE)
	{//set a new point
		aPoint.x = new_x;
		aPoint.y[0] = yVal;
		spec.push_back(aPoint);
		lastChargeIdx = 0;
	}
	else
	{//add to existing
		spec[idx].x = (new_x + spec[idx].x) / 2; //average
		if (lastChargeIdx < MAX_COMBINED_SPECTRA_NR)
		{
			for (int c = 0; c < MAX_COMBINED_SPECTRA_NR; c++)
			{
				if (spec[idx].y[c] < 1)
				{
					lastChargeIdx = c;
					break;
				}
			}
			spec[idx].y[lastChargeIdx] += yVal;

		}
		else
		{
			WRITE_DEBUG("To many charge states added")
				SOFTWARE_ERROR;
		}
	}
}

void TMRC_Iso::CalcMRC()
{//calculate a sum for all charge states
	for (int i = 0; i < spec.size(); ++i)
	{
		double ysum = 0; //sum intensity
		int count = 0; //number of peaks
		for (int c = 0; c <= lastChargeIdx; c++)
		{
			if (spec[i].y[c] > 0)
			{
				count++;
				ysum += spec[i].y[c];
			}
		}
		if (count < lastChargeIdx - 1) ysum = 0;
		spec[i].y[MAX_COMBINED_SPECTRA_NR] = ysum; //use last item as summed intensity
	}
}

void TMRC_Iso::Write_Log()
{
	for (int i = 0; i < spec.size(); ++i)
	{
		string tmp = FormatDouble(spec[i].x, 4) + "\t";
		for (int c = 0; c <= MAX_COMBINED_SPECTRA_NR; c++)
		{
			tmp += FormatDouble(spec[i].y[c], 2) + "\t";
		}
		WRITE_INFO(tmp);
	}
}


#pragma region TEXIC

void TExIC::WriteToFile(TBinaryFileWriter* write)
{
	if (write != nullptr)
	{
		write->WriteValue(mass);
		write->WriteValue((int)name.length());
		write->WriteValue(name);
		write->WriteValue(charge);

		write->WriteValue(isotopeNr);

		for (int c = 0; c < isotopeNr; c++)
		{
			write->WriteValue(intensity[c]);
			write->WriteValue(MZ[c]);
		}
	}
}

void TExIC::ReadFromFile(TBinaryFileWriter* read)
{
	//SOFTWARE_ERROR
	if (read != nullptr)
	{
		read->ReadValue(&mass);
		int name_Length;
		read->ReadValue(&name_Length);
		read->ReadValue(&name);
		read->ReadValue(&charge);
		read->ReadValue(&isotopeNr);
		int array_size = isotopeNr;
		SAFE_DELETE_ARRAY(intensity);
		SAFE_DELETE_ARRAY(MZ);
		intensity = new double[array_size];
		MZ = new double[array_size];

		for (int c = 0; c < isotopeNr; c++)
		{	
			read->ReadValue(&intensity[c]);
			read->ReadValue(&MZ[c]);		
		}
	}
}

void TExIC::Copy(TExIC* inPut)
{
	//copy data from another ExIC
	if (inPut != nullptr)
	{
		mass = inPut->mass;
		name = inPut->name;
		charge = inPut->charge;
		isotopeNr = inPut->isotopeNr;
		maxIsotopes = inPut->maxIsotopes;
		delete[] intensity;
		intensity = new double[maxIsotopes];
		delete[] MZ;
		MZ = new double[maxIsotopes];
		delete[] trueRT;
		trueRT = new double[maxIsotopes];
		for (int i = 0; i<maxIsotopes; i++)
		{
			MZ[i] = inPut->MZ[i];
			intensity[i] = inPut->intensity[i];
			trueRT[i] = inPut->trueRT[i];
		}
		
		
		maxScanInt.assign(inPut->maxScanInt.begin(), inPut->maxScanInt.end());
		maxScanLoc.assign(inPut->maxScanLoc.begin(), inPut->maxScanLoc.end());
	}
}

void TExIC::Init()
{
	mass = numeric_limits<double>::quiet_NaN();
	name = "";
	charge = 0;
	//ToDo : make possible setting this value
	isotopeNr = 0;
	maxIsotopes = 0;

	intensity = nullptr;
	MZ = nullptr;
	trueRT = nullptr;

	isotopeNr = ZERO_CONST.it;
}

TExIC::~TExIC()
{
	SAFE_DELETE_ARRAY(intensity);
	SAFE_DELETE_ARRAY(MZ);		
	SAFE_DELETE_ARRAY(trueRT);
}

void TExIC::WriteToStdOut()
{//std::cout << peptides[pepIdx].composition.isoMass[idx] << "\t" << foundmz << "\t" << RT[s] << "\t" << intensity << std::endl;
	for (size_t i = 0; i < isotopeNr; i++)
	{
		std::cout << "Results for scan index: " << i << std::endl;
		for (size_t x = 0; x < isotopeNr; x++)
		{
			if (intensity[x] > BASELINE)
			{
				std::cout << name << "\t" << x << "\t" << MZ[x] << "\t" << intensity[x] << std::endl;
			}
			
		}
	}
}


#pragma endregion TEXIC
#pragma region TOutput
// =======================================================================================================================
//    Init new array for XIC-s input: nr-maxim number of XICs, scans-nr of spectra in raw file
// =======================================================================================================================
void TOutput::AddICs(int fileNr, int maxD) // files x isotopeNr
{
	//create an Extracted Ion Chromatogram ExIC for each file. 
	//will contain all the isotopic data
	SAFE_DELETE_ARRAY(ICs);
	SAFE_DELETE_ARRAY(deutResults);
	ICs = new TExIC[fileNr];
	deutResults = new TDeutData[fileNr];
	
	fileXICNr = fileNr;

	maxIsotopes = maxD;
	
	for (size_t i = 0; i < fileNr; i++)
	{
		ICs[i].intensity = new double[maxD];
		ICs[i].MZ = new double[maxD];
		ICs[i].trueRT = new double[maxD];
		ICs[i].maxIsotopes = maxD;
		for (int n = 0; n < maxIsotopes; n++)
		{
			ICs[i].intensity[n] = 0.0;
			ICs[i].MZ[n] = 0.0;
			ICs[i].trueRT[n] = 0.0;
		}
	}
	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Function:  Initialize vales for an extracted ion chromatogram
// Result:    None
// Parameters: double mass - mass of the ICs that we want to initialize
//		       std::string name - name of the ICs
//			   double *RT - pointer to an array containing the retention time values
//			   int start - starting retention time
//			   int end - ending retention time
//			   int charge - charge state of the compound the ICs will contain
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void TOutput::ZeroIC(string name, int icsIndex)
{
	//~~~~~~~~~~~~~~~~~
	// add a new point to the correct ICs;
	ICs[icsIndex].name = name;

	for (int IsoIdx = 0; IsoIdx < maxIsotopes; IsoIdx++ )
	{	
		ICs[icsIndex].intensity[IsoIdx] = 0.0000;
		ICs[icsIndex].MZ[IsoIdx] = 0.0000;
	}
}

// =======================================================================================================================
// Store an intensity "y" The found intensity has "ppm" mass uncertainty.
// The new point is stored in the extracted ion chromatogram specified by index "index" and has charge state - charge
// RT is the retention time value where the new point should be stored, mass holds the mass to be set
// =======================================================================================================================
void TOutput::AddPointToExIC(double mass, double y, double mz, int icsIndex, int IsoIdx)
{
	if (icsIndex >= fileXICNr)
	{
		SOFTWARE_ERROR;
	}

	if (IsoIdx == -1)
	{
		WRITE_DEBUG("Isotope number error");
		SOFTWARE_ERROR;
	}
	if (ICs[icsIndex].isotopeNr >= maxIsotopes )
	{
		WRITE_DEBUG("Isotope too many; error");
		SOFTWARE_ERROR;
	}

	ICs[icsIndex].mass = mass;
	ICs[icsIndex].intensity[IsoIdx] = y;
	ICs[icsIndex].MZ[IsoIdx] = mz;
	if (ICs[icsIndex].isotopeNr < IsoIdx + 1)
	{
		ICs[icsIndex].isotopeNr = IsoIdx + 1;
	}
	
}


void TOutput::WriteToFile(TBinaryFileWriter* write)
{
	if (write != nullptr)
	{
		write->WriteValue(CURRFILEVER.bvMajor); // file type 1.0 contains the nr of ICs and then the ICS structure
		write->WriteValue(CURRFILEVER.bvMinor); // file type 1.0 contains the nr of ICs and then the ICS structure
		write->WriteValue(fileXICNr);
		if (ICs != nullptr)
		{
			for (int i = 0; i < fileXICNr; i++)
			{
				ICs[i].WriteToFile(write);
			}
		}
	}
}

void TOutput::ReadFromFile(TBinaryFileWriter* read)
{
	//SOFTWARE_ERROR
	if (read != nullptr)
	{
		//
		int bvMajor = -1; read->ReadValue(&bvMajor);
		int bvMinor = -1; read->ReadValue(&bvMinor);
		if (CURRFILEVER.bvMajor != bvMajor || CURRFILEVER.bvMinor != bvMinor)
		{//check if binary file version is OK
			SOFTWARE_ERROR;
		}

		read->ReadValue(&fileXICNr);
		if (fileXICNr > 0)
		{
			SAFE_DELETE_ARRAY(ICs);
			ICs = new TExIC[fileXICNr];
			for (int i = 0; i < fileXICNr; i++)
			{
				ICs[i].ReadFromFile(read);
			}
		}
	}
}

void TOutput::Copy(TOutput* inPut)
{
	//copy one data to another
	fileXICNr = inPut->fileXICNr;

	maxIsotopes = inPut->maxIsotopes;

	SAFE_DELETE_ARRAY(ICs);
	ICs = new TExIC[fileXICNr];
	for (int i = 0; i < fileXICNr; i++)
	{
		ICs[i].Copy(inPut->ICs);
	}
}

void TOutput::Reset() // Delete all data
{
	SAFE_DELETE_ARRAY(ICs);
	fileXICNr = 0; maxIsotopes = 0; 
}

void TOutput::ReInit()
{
	

	refMassData.deutDelta = ERROR_VALUE;
	refMassData.massCenter = ERROR_VALUE;
	for (int f = 0; f < fileXICNr; f++)
	{
		deutResults[f].deutDelta = ERROR_VALUE;
		deutResults[f].massCenter = ERROR_VALUE;
		for (int n = 0; n < maxIsotopes; n++)
		{
			ICs[f].intensity[n] = 0.0;
			ICs[f].MZ[n] = 0.0;
			ICs[f].trueRT[n] = 0.0;
		}
	}
	
}

void TOutput::WriteToStdOut()
{
	//std::cout << peptides[pepIdx].composition.isoMass[idx] << "\t" << foundmz << "\t" << RT[s] << "\t" << intensity << std::endl;
	for (size_t i = 0; i < fileXICNr; i++)
	{
		std::cout << "Results for file index: " << i << std::endl;
		ICs[i].WriteToStdOut();
	}
}

void TOutput::CalculateRefMassData()
{//check for all data with 0 sec deut time

	int count = 0;
	refMassData.massCenter = 0.0;
	for (size_t i = 0; i < fileXICNr; i++)
	{
		if (deutResults[i].deutTime <= NUMERICAL_ZERO_DOUBLE)
		{
			refMassData.massCenter += deutResults[i].massCenter;
			//std::cout << i << "\t" << deutResults[i].deutTime << "\t" << refMassData.massCenter << std::endl;
			count++;
		}
	}
	if (count>0)
	{
		refMassData.massCenter = refMassData.massCenter / count;
	}
	

	//std::cout << refMassData.massCenter  <<std::endl;
}

#pragma endregion TOutput
