#include "TPeptideStore.h"
#include "peptide.h"
#include <math.h>       /* isnan, sqrt */
#include <vector>
#include <numeric>
#include <iostream>
#include "TPeakDetect.h"
#include <algorithm>

float average(std::vector<double> const& v) {
	if (v.empty()) {
		return 0;
	}

	auto const count = static_cast<float>(v.size());
	return std::reduce(v.begin(), v.end()) / count;
}

bool TPeptideStore::DeletePeptide(int idx)
{
	if (idx < pepNr && idx >= 0)
	{
		//peptides[p].ClearData();
		results[idx].ReInit();

		return true;
	}

	return false;
}

int TPeptideStore::ReadPeptideFile(std::string fileName)
{
	TStringList dataList;
	TStringList splitList;
	TIsoCalc	calcIsotopes;
	splitList.separator = '\t';
	if (FileExists(fileName))
	{
		dataList.LoadFromFile(fileName);
	}
	else
	{
		pepNr = ERROR_VALUE;
		return pepNr;
	}

	//ClearData();
	SAFE_DELETE_ARRAY(peptides);
	SAFE_DELETE_ARRAY(results);
	peptides = new TPeptideData[dataList.size()];
	results = new TOutput[dataList.size()];

	pepNr = 0;

	double		mass = 0, sRT, eRT;
	int			charge;
	std::string tmp;
	TPeptide aPeptide;

	for (unsigned int i = 0; i < dataList.size(); ++i)
	{
		splitList.SetText(dataList[i]);
		peptides[i].proteinID = trim(splitList[inpProteinCOL]);
		peptides[i].sequence = trim(splitList[inpPeptideCOL]);
		peptides[i].mass = aPeptide.GetMass(peptides[i].sequence);
		peptides[i].composition.SetComposition(aPeptide.composition, ATOMTYPES);
		
		calcIsotopes.CopyCompositionData(&peptides[i].composition);		
		calcIsotopes.CalcIsoCluster();
		peptides[i].composition.CopyComposition(&calcIsotopes);
		peptides[i].CalcMaxDeuteration();
		
		try //try to convert proteinID to seq position
		{
			auto sepPos = peptides[i].proteinID.find("-");
			peptides[i].seqStart = stoi(peptides[i].proteinID.substr(0, sepPos));
			peptides[i].seqEnd = stoi(peptides[i].proteinID.substr(sepPos + 1));
		}
		catch (...)
		{
			//do nothing
		}

		try
		{
			SafeConvert(splitList[inpRTsCOL], &sRT);
			SafeConvert(splitList[inpRTeCOL], &eRT);
			charge = stoi(splitList[inpZCOL]);
			peptides[i].endRT = eRT;
			peptides[i].startRT = sRT;
			peptides[i].z = charge;
		}
		catch (...)
		{
			//Dialogs::ShowMessage(L"Error in RT input list!");
			WRITE_DEBUG("Error in input list!");
			SOFTWARE_ERROR;
		}

		pepNr++;
	}
	return pepNr;
}

bool TPeptideStore::CreateResultStorage(TRawDataFile* rawData, int numFiles)
{
	if (results != nullptr && rawData != nullptr && numFiles > 0 && pepNr > 0)
	{
		fileNr = numFiles;
		for (size_t i = 0; i < pepNr; i++)
		{
			results[i].AddICs(numFiles, MAX_ISOTOPE_SHIFT);
		}		
		SAFE_DELETE_ARRAY(dataInfo);
		dataInfo = new TRawDataInfo[numFiles];
		SAFE_DELETE_ARRAY(xEIC);
		xEIC = new EIC[numFiles * pepNr + 1 ];  //avoid init with 0
		for (size_t i = 0; i < numFiles * pepNr + 1; i++)
		{
			xEIC[i].rt.clear();
			xEIC[i].y.clear();
		}
		return true;	
	}
	
	return false;
}

bool TPeptideStore::SetRawInfo(std::string path, double time, int idx, std::string st/*=""*/, int repl/*=0*/)
{
	if (dataInfo != nullptr && idx < fileNr)
	{
		dataInfo[idx].labelTime = time;
		dataInfo[idx].rawPath = path;
		dataInfo[idx].isRef = time <= 0.01;
		dataInfo[idx].replicate = repl;
		dataInfo[idx].state = st;
		if (dataInfo[idx].isRef) refFileIdx = idx;
		return true;
	}

	return false;
}

bool TPeptideStore::SearchAPeptide(TRawDataFile* rawData, int pepIdx, int fileIdx, bool isDeuterated)
{
	//read the file data and check the corresponding isotopic peaks;
	//do search
	results[pepIdx].deutResults[fileIdx].deutTime = dataInfo[fileIdx].labelTime;
	results[pepIdx].ZeroIC(peptides[pepIdx].sequence, fileIdx);
	this->rawData = rawData;

	double	mz, ppm = DEUTMASSERRORPPM;
	long	start = ERROR_VALUE, end = ERROR_VALUE;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Check if there is anything to search for
	long		rawSpecNr = rawData->GetSpecNr();
	double*		RT = new double[rawSpecNr];
	double*		RTIDX = new double[rawSpecNr];
	double		maxInt = -1.0;  //maximum intensity in all xEIC
	size_t		maxIntPos = 0;

	for (long i = 0; i < rawSpecNr; i++)
	{
		RT[i] = rawData->GetFileDataRT(i);
	}
	
	//check if valid retention times
	try
	{
		//set start and end scans
		start = rawData->GetScanAtTime(peptides[pepIdx].startRT);
		end   = rawData->GetScanAtTime(peptides[pepIdx].endRT);
	}
	catch (exception * e)
	{
		WRITE_DEBUG(e->what());
		return false;
	}

	if (start == -1) start = 1;
	if (end < start) end = rawData->GetSpecNr();
	
	int maxPepD = peptides[pepIdx].composition.isoDataNr;
	int minIsoToConsider = 0; //calculate the number of isotopes to use minimum

	for (size_t idx = 0; idx < peptides[pepIdx].composition.isoDataNr; idx++)
	{
		if (peptides[pepIdx].composition.isoData[idx]  > MIN_ISOTOPE_ABUNDANCE) minIsoToConsider++;
		maxPepD = minIsoToConsider;
	}

	if (isDeuterated) maxPepD = peptides[pepIdx].CalcMaxDeuteration() + minIsoToConsider + 1;
	
	//create temp storage for isotopes to be searched maxPepD in number
	double** tempY = new double* [maxPepD];
	double* tempMass = new double[maxPepD];
	
	//fill up with data
	for (size_t i = 0; i < maxPepD; i++)
	{
		//isDeuterated = true;
		tempY[i] = new double[end - start + 1];
		for (size_t x = 0; x < end-start + 1; x++)
		{
			tempY[i][x] = 0.0;
		}
		if (isDeuterated && i > peptides[pepIdx].composition.maxIdx)
		{
			tempMass[i] = peptides[pepIdx].composition.isoMass[peptides[pepIdx].composition.maxIdx] + (masses::D - masses::H) * (i- peptides[pepIdx].composition.maxIdx);
		}
		else
		{
			tempMass[i] = peptides[pepIdx].composition.isoMass[i];
		}
	}

	//vectors to hold data about the found isotopic peaks
	vector<int> maxLocVec(maxPepD,0); //location of maximum
	vector<double> maxIntVec(maxPepD,0.0); //intensity of maximum 
	//auto p = make_unique<int*[]>(5);
	//p[0] = new int[5];
	double sDelta = 0.0; //mass error
	//search all scans for the calculated isotope masses in tempMass

	for (size_t s = start; s < end; s++)
	{//check all scans	
		RTIDX[s] = s;
		double sumInt = 0.0; //sum of all peaks found in a scan
		int isoCount = 0;
		int foundCases = 0; //how many isotopes found
		//isoInScan[s] = make_unique<double *>(new double[maxPepD]);
		
		//use max intensity of a spectrum to cut the small peaks
		double scanMaxInt = 0;
		
		for (size_t idx = 0; idx < maxPepD; idx++)
		{// cycle for all isotopes	

			tempY[idx][s - start] = 0.0;
			mz = tempMass[idx] / peptides[pepIdx].z + masses::Hp;
			sDelta = ppm * mz / 1.0e6;
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			int found = BinarySearchClosest((double*)rawData->GetFileDataMZArray(s), (double*)rawData->GetFileDataYArray(s), mz, 0, rawData->GetFileDataPeakNr(s), BASELINE, sDelta, false);
			double intensity = 0.0;
			double foundmz = 0.0;
			
			if (found > 0)
			{	//found the peaks
				intensity = rawData->GetFileDataY(s, found);
				foundmz = rawData->GetFileDataMZ(s, found);
				isoCount++;
				sumInt += intensity;
				foundCases++;
				if (intensity > scanMaxInt)
				{
					scanMaxInt = intensity;
				}
				if (intensity > maxIntVec[idx])
				{
					maxIntVec[idx] = intensity;
					maxLocVec[idx] = s;
				}
			}
			//isoInScan[s][idx] = intensity;
			
			tempY[idx][s - start] = intensity;
		}//end isotopes



		for (size_t idx = 0; idx < maxPepD; idx++)
		{
			if (tempY[idx][s - start] < scanMaxInt / MAX_DYNAMIC_RANGE) tempY[idx][s - start] = 0.0;
		}
		
		//check if this is the highest spectrum
		if (sumInt > maxInt && (isoCount >= minIsoToConsider)) //find at least some % of isotopes
		{
			maxInt = sumInt;
			maxIntPos = s-start; //collation of the best spectrum
		}
		results[pepIdx].ICs[fileIdx].isotopeNr = max(isoCount, results[pepIdx].ICs[fileIdx].isotopeNr);
	}//end scans

	
	//set retention data of the isotopic peaks; what is the location of their maximum and what is the intensity
	//use this to filter out overlaping peaks
	results[pepIdx].ICs[fileIdx].maxScanInt.assign(maxIntVec.begin(), maxIntVec.end());
	results[pepIdx].ICs[fileIdx].maxScanLoc.assign(maxLocVec.begin(), maxLocVec.end());

#if 1
	if (useXICOverlap) {
		// if desired, we can detect chroma peaks for all isotopes	
		//create storage for results
#include "savgol.h"
		SavGol* sg = new SavGol;

		struct peakData
		{
			double loc[3]; //peak location
			double y[3]; //peak intensity
			int size = -1; //number of peaks
		};
		vector<peakData> peaks(maxPepD);
		//create storage for peak param, common
		double commonRT = 0.0, commonStart = 0.0, commonEnd = 0.0, prevInt = 0.0;
		//ToDo: current algorithm assumes that the isotope 0 is found and is the target isotope
		//		and also that the most intense chromatographic peak will be the target
		double* smoothY = new double[end - start];
		int pwidth = 0;
		for (int iso = 0; iso < maxPepD; iso++)
		{//for all isotopes
			peaks[iso].size = 0;
			TPeakDetect tp;
			sg->DoSmoothWithOutput(tempY[iso], smoothY, end - start, 7, 1);
			tp.SetData(&RTIDX[start], smoothY, end - start);
			tp.DetectPeaks();
			//ToDo: here do something with the detected peaks
			//tp.WriteData(to_string(pepIdx) +"_" +to_string(fileIdx) + "_" + to_string(iso));
			int pn = tp.GetPeakNr();
			peaks[iso].size = min(pn, 3);
			for (int n = 0; n < peaks[iso].size; ++n)
			{	//check peak RT						
				(peaks[iso].loc[n] = tp.GetPeakRT(n));

				(peaks[iso].y[n] = tp.peakMaxInt[n]);
				pwidth = tp.peakEnd[n] - tp.peakStart[n];
			}
		}
		SAFE_DELETE_ARRAY(smoothY);
		SAFE_DELETE(sg);

		double rtDelta = pwidth*0.5; //difference that we count as 0
		// make it dependent on the width of the current chromatographic peak
		//find the first isotope that is not 0 intensity and find its highest peak???
		int maxIntLoc = 0; //max intensity location
		bool foundMax = false;
		double maxLocalInt = BASELINE;
		for (int iso = peptides[pepIdx].composition.maxIdx; iso < maxPepD/2.0; iso++)
		{
			//int iso = peptides[pepIdx].composition.maxIdx;
			for (int s = 0; s < peaks[iso].size; s++)
			{
				if (peaks[iso].y[s] > maxLocalInt)
				{
					maxLocalInt = peaks[iso].y[s];
					maxIntLoc = (long)peaks[iso].loc[s];
					foundMax = true;
				}
			}
			if (foundMax) break; //found a good isotope and its location. Use that as a reference
		}

		//std::cout << start <<"__" << maxIntPos << "__" << maxIntLoc << std::endl;

		for (int iso = maxPepD/2.0-1; iso < maxPepD; iso++)
		{//for all isotopes
			//std::cout << iso << "__";
			bool isoFound = false;
			for (int s = 0; s < peaks[iso].size; s++)
			{
				//std::cout << s << "__";
				//if such a peak is found for current isotope, then OK, else erase
				//std::cout << peaks[iso].loc[s] << "__" << peaks[iso].y[s] << std::endl;
				if (fabs(peaks[iso].loc[s] - maxIntLoc) < rtDelta)
				{
					isoFound = true;
					//std::cout << iso << "___";
				}
			}
			//std::cout << std::endl;
			if (!isoFound)
			{
				for (int p = start; p < end; p++)
				{
					tempY[iso][p - start] = 0.0;
				}
			}
			//std::cout << std::endl;
		}
		maxIntPos = max(maxIntLoc - (int)start, 0);
		//this->peptides[pepIdx].			maxIntLoc = maxIntPos + start;

		//WRITE_INFO(maxIntPos+start);
		//WRITE_INFO(RT[maxIntPos + start]);
		//std::cout << maxIntPos << std::endl;
	}
#endif	


	//store the position of peptide maximum
	for (size_t i = 0; i < maxPepD; i++)
	{
		double isoMaxInt = 0.0;
		for (size_t s = start; s < end; s++)
		{
			//always get the max int;
			if (isoMaxInt < tempY[i][s - start])
			{
				isoMaxInt = tempY[i][s - start];
				results[pepIdx].ICs[fileIdx].trueRT[i] = RT[s]; //store as seconds and not scans;
			}
		}
	}

	
	//create average spectra for isotopic pattern comparison
	double scanMaxMass = CreateAverageIsoPattern(pepIdx, fileIdx, start, end, maxIntPos, maxPepD, tempMass, tempY);

	//CheckCoelution(pepIdx, fileIdx);
	//if not found, set data to theoretical value
	if (scanMaxMass < 100) scanMaxMass = peptides[pepIdx].composition.isoMass[0] / peptides[pepIdx].z + masses::Hp;
	//set data for extracted ion chrom
	CreateXIC(pepIdx, fileIdx, scanMaxMass, sDelta, rawData);

	//return true;
	for (size_t i = 0; i < maxPepD; i++)
	{
		SAFE_DELETE_ARRAY(tempY[i]);
	}
	SAFE_DELETE_ARRAY(tempMass);
	SAFE_DELETE_ARRAY(tempY);
	SAFE_DELETE_ARRAY(RT);
	SAFE_DELETE_ARRAY(RTIDX);

	//WriteResultsToStdOut();
	//return true;

	return false;
}


void TPeptideStore::WriteResultsToStdOut()
{
	for (size_t i = 0; i < pepNr; i++)
	{
		std::cout << "Writing results for peptide: " << peptides[i].proteinID << "\t" << peptides[i].sequence << std::endl;
		results[i].WriteToStdOut();

	}
}

bool TPeptideStore::CheckCoelution(int pepIdx, int fileIdx)
{//needs to return bool
	//calc the max loc and set bad ones to 0
	auto ret = freopen("dbg.txt", "w", stdout);
	std::vector<int> v;
	/*for (int x = 0; x < results[pepIdx].ICs[fileIdx].maxScanLoc.size(); x++)
	{
		v.push_back(results[pepIdx].ICs[fileIdx].maxScanLoc[x]);
	}

	int mElem = v.size() / 2;
	std::cout << v.size() <<"  " << mElem << "  " << results[pepIdx].ICs[fileIdx].isotopeNr << std::endl;
	if (v.size() < 2) return false;

	auto m = v.begin() + v.size() / 2;
	std::nth_element(v.begin(), m, v.end());
	
	
	auto median = v[mElem];
	//std::cout << "\nThe median is " << v[mElem] << '\n';
	int prevZero = 0; //number of prev bad ones
	for (int i = 0; i < results[pepIdx].ICs[fileIdx].isotopeNr; i++)
	{
		//set to 0 the bad ones
		std::cout << results[pepIdx].ICs[fileIdx].intensity[i] << "\t" << results[pepIdx].ICs[fileIdx].maxScanLoc[i] << std::endl;
		
		if (results[pepIdx].ICs[fileIdx].maxScanLoc[i] > median + 4 || results[pepIdx].ICs[fileIdx].maxScanLoc[i] < median - 4)
		{
			prevZero++;
			results[pepIdx].ICs[fileIdx].intensity[i] = 0.0;
		}
		if (prevZero > 2)
		{
			results[pepIdx].ICs[fileIdx].intensity[i] = 0.0;
		}
	}*/
	SearchAPeptide(rawData, pepIdx, fileIdx, true);
	CalcPeptideMassCenter(pepIdx, fileIdx);
	CalcDeuteration(pepIdx, fileIdx);

	fclose(stdout);
	return false;
}

double TPeptideStore::CreateAverageIsoPattern(int pepIdx, int fileIdx, int sScan, int eScan, int maxPos, int maxPepD, double* Mass, double** Y)
{
	//create average spectra for isotopic pattern comparison
	double resultIsotopeDistrib = 0.0;
	double maxIsoInt = 0.0;
	double scanMaxMass = 0.0;
	
	//std::cout << sScan << std::endl;
	//std::cout << eScan << std::endl;
	//auto ret = freopen("dbg.txt", "w", stdout);
	for (size_t i = 0; i < maxPepD; i++)
	{
		//toDo: boundary checks
		//std::cout << Mass[i] / peptides[pepIdx].z + masses::Hp << "\t";
		resultIsotopeDistrib = 0.0;
		int isoCounter = 0;
		for (size_t r = max(0,maxPos - AVGSPECISOTOPENR); r < min(maxPos+ AVGSPECISOTOPENR,eScan) ; r++)
			//ToDo: create average of spectra
		{
			if (Y[0][r] > ISOTOPE_CENTROID_TOP * Y[0][maxPos])
			{
				//std::cout << Y[i][r] << "\t";
				resultIsotopeDistrib += Y[i][r];
				isoCounter++;
			}

		}
		//std::cout << std::endl;
		//std::cout << isoCounter << std::endl;
		if (isoCounter < 1) isoCounter = 1;
		if (maxIsoInt < resultIsotopeDistrib / isoCounter)
		{
			maxIsoInt = resultIsotopeDistrib / isoCounter;
			scanMaxMass = Mass[i] / peptides[pepIdx].z + masses::Hp;
		}
		results[pepIdx].AddPointToExIC
		(
			peptides[pepIdx].mass,
			resultIsotopeDistrib / isoCounter,
			Mass[i] / peptides[pepIdx].z + masses::Hp,
			fileIdx,
			i
		);
	}
	//fclose(stdout);
	return scanMaxMass;
}

void TPeptideStore::CreateXIC(int pepIdx, int fileIdx, double mass, double sDelta, TRawDataFile* currRawData)
{
	
	if (!currRawData) return;
	
	xEIC[pepIdx + fileIdx * pepNr].baseMass = mass;
	xEIC[pepIdx + fileIdx * pepNr].rt.clear();
	xEIC[pepIdx + fileIdx * pepNr].y.clear();
	//extract ion chrom
	//sDelta = 0.15; //get XIC with high mass window
	for (long s = 1; s < currRawData->GetSpecNr() - 1; s++)
	{
		double fY = 0.0; //result
		int found = BinarySearchClosest((double*)currRawData->GetFileDataMZArray(s), (double*)currRawData->GetFileDataYArray(s), mass, 0, currRawData->GetFileDataPeakNr(s), 0, sDelta, false);
		if (found > 0)
		{
			fY = rawData->GetFileDataY(s, found);
		}
		//set RT
		xEIC[pepIdx + fileIdx * pepNr].rt.push_back(currRawData->GetFileDataRT(s));
		//set intensity
		xEIC[pepIdx + fileIdx * pepNr].y.push_back(fY);
	}
}

double TPeptideStore::CalcPeptideMassCenter(int pepIdx, int fileIdx)
{//use the following formula
	double WeightedSumInt = 0.0;
	double SumInt = 0.0;
	double maxInt = 0.0; //base peak intensity
	int maxLoc = 0; //location of base peak

	int zeroCount = 0; //number of consecutive zeros
	int endPos = results[pepIdx].ICs[fileIdx].isotopeNr; //end of isotopes

	for (int i = 0; i < results[pepIdx].ICs[fileIdx].isotopeNr; i++)
	{//sum intensities
		if (maxInt > BASELINE && results[pepIdx].ICs[fileIdx].intensity[i] < maxInt/ MAX_DYNAMIC_RANGE)
		{
			zeroCount++;
		}
		if (zeroCount > 2) {
			endPos = i;
			break;
		}
		if (results[pepIdx].ICs[fileIdx].intensity[i] > maxInt)
		{
			maxInt = results[pepIdx].ICs[fileIdx].intensity[i];
			maxLoc = i;
		}
	}

	int start = max(0, maxLoc - maxIsotopeNr);
	int end = min(endPos, maxLoc + maxIsotopeNr);

	for (size_t i = start; i < end; i++)
	{//sum intensities

			SumInt += results[pepIdx].ICs[fileIdx].intensity[i] ;
			WeightedSumInt += results[pepIdx].ICs[fileIdx].intensity[i] * results[pepIdx].ICs[fileIdx].MZ[i];
	}

	double massCenter = peptides[pepIdx].z * (WeightedSumInt / SumInt);
	results[pepIdx].deutResults[fileIdx].massCenter = massCenter - (peptides[pepIdx].z * masses::Hp);

	return massCenter;
}

double TPeptideStore::CalcReferenceMassCenter(int pepIdx)
{
	//create an average for reference data and calculate its center
	results[pepIdx].CalculateRefMassData();
	//calculate theoretical mass center
	results[pepIdx].isoRefData.deutTime = 0.0;
	results[pepIdx].isoRefData.deutDelta = 0.0;
	results[pepIdx].isoRefData.massCenter = 0;
	int pointNr = peptides[pepIdx].composition.isoDataNr;
	//for each file check if it is reference
	double aCenter = 0.0, intensityWeight = 0.0;
	for (int idx = 0; idx < pointNr; idx++)
	{
		aCenter += peptides[pepIdx].composition.isoMass[idx] * peptides[pepIdx].composition.isoData[idx];
		intensityWeight += peptides[pepIdx].composition.isoData[idx];
	}
	return results[pepIdx].isoRefData.massCenter = aCenter / intensityWeight;
}

double TPeptideStore::CalcDeuteration(int pepIdx, int fileIdx, int refIdx)
{
	double deut = 0.0;
	if (refIdx == ERROR_VALUE)
	{
		deut = (results[pepIdx].deutResults[fileIdx].massCenter - results[pepIdx].isoRefData.massCenter);
	}
	else
	{
		deut = (results[pepIdx].deutResults[fileIdx].massCenter - results[pepIdx].deutResults[refIdx].massCenter);
	}
	deut = (results[pepIdx].deutResults[fileIdx].massCenter - results[pepIdx].isoRefData.massCenter);
	// (result 100% - results[pepIdx].deutResults[refIdx].massCenter)) * peptides[pepIdx].maxD;
	//toDo: incorporate fully deuterated data 
	results[pepIdx].deutResults[fileIdx].deutDelta = deut; // / (peptides[pepIdx].maxD * (masses::D - masses::AH);
	return deut;
}

void TPeptideStore::WriteDeuterationDataToStdOut()
{
	std::cout << "File\tTime\t";
	for (size_t p = 0; p < pepNr; p++)
	{
		bool isBad = false;
		for (size_t i = 0; i < fileNr; i++)
		{
			if (isnan(results[p].deutResults[i].massCenter)) isBad = true;
		}
		if (isBad) continue;
		std::cout << peptides[p].sequence << "\t";
	}
	std::cout << std::endl;

	for (size_t f = 0; f < fileNr; f++)
	{
		std::cout << dataInfo[f].rawPath << "\t" << dataInfo[f].labelTime;
		for (size_t p = 0; p < pepNr; p++)
		{
			bool isBad = false;
			for (size_t i = 0; i < fileNr; i++)
			{
				if (isnan(results[p].deutResults[i].massCenter)) isBad = true;
			}
			if (isBad) continue;

			 std::cout << "\t" << results[p].deutResults[f].deutDelta;
		}
		std::cout << std::endl;
	}
	//<< "\t" << peptides[p].maxD << "\t" << results[p].deutResults[f].deutDelta/(peptides[p].maxD * (masses::D-masses::H))
}

void TPeptideStore::WriteDeuterationPercentToStdOut()
{
	std::cout << "File\tTime\t";
	for (size_t p = 0; p < pepNr; p++)
	{
		bool isBad = false;
		for (size_t i = 0; i < fileNr; i++)
		{
			if (isnan(results[p].deutResults[i].massCenter)) isBad = true;
		}
		if (isBad) continue;
		std::cout << " " << peptides[p].proteinID << "\t"; //put space before, because Excel needs it to display correctly
	}
	std::cout << std::endl;

	for (size_t f = 0; f < fileNr; f++)
	{
		std::cout << dataInfo[f].rawPath << "\t" << dataInfo[f].labelTime;
		for (size_t p = 0; p < pepNr; p++)
		{
			bool isBad = false;
			for (size_t i = 0; i < fileNr; i++)
			{
				if (isnan(results[p].deutResults[i].massCenter)) isBad = true;
			}
			if (isBad) continue;

			std::cout << "\t" << results[p].deutResults[f].deutDelta / (peptides[p].maxD * (masses::D - masses::H));
		}
		std::cout << std::endl;
	}
	//<< "\t" << peptides[p].maxD << "\t" << results[p].deutResults[f].deutDelta/(peptides[p].maxD * (masses::D-masses::H))
}

void TPeptideStore::WriteDeuterationCentersToStdOut()
{//write only centers
	std::cout << "File\tTime\t";
	for (size_t p = 0; p < pepNr; p++)
	{
		bool isBad = false;
		for (size_t i = 0; i < fileNr; i++)
		{
			if (isnan(results[p].deutResults[i].massCenter)) isBad = true;
		}
		if (isBad) continue;
		std::cout << (peptides[p].endRT + peptides[p].startRT)/2 << "\t";
	}
	std::cout << std::endl;

	for (size_t f = 0; f < fileNr; f++)
	{
		std::cout << dataInfo[f].rawPath << "\t" << dataInfo[f].labelTime;
		for (size_t p = 0; p < pepNr; p++)
		{
			bool isBad = false;
			for (size_t i = 0; i < fileNr; i++)
			{
				if (isnan(results[p].deutResults[i].massCenter)) isBad = true;
			}
			if (isBad) continue;

			std::cout << "\t" << results[p].deutResults[f].massCenter;
		}
		std::cout << std::endl;
	}
}

void TPeptideStore::WriteDeuterationInMEMHDX()
{
	//store specific format to MEMHDX
	//Start	End	Sequence	MaxUptake	State	Exposure	z	Replicate	Center
	std::cout << "Start,End\tSequence\tMaxUptake\tState\tExposure\tz\tReplicate\tCenter" << std::endl;
	for (size_t p = 0; p < pepNr; p++)
	{
		bool isBad = false;
		for (size_t i = 0; i < fileNr; i++)
		{
			if (isnan(results[p].deutResults[i].massCenter)) isBad = true;
			if ((results[p].deutResults[i].massCenter) < 0) isBad = true;
		}
		if (isBad) continue;
		
		for (size_t f = 0; f < fileNr; f++)
		{


		//write peptide data
		std::cout << peptides[p].proteinID << "\t" << peptides[p].sequence << "\t" << peptides[p].CalcMaxDeuteration() << "\t";

		std::string path = dataInfo[f].rawPath;
		std::size_t botDirPos = path.find_last_of("\\");
		string state = dataInfo[f].state;
		//path.substr(botDirPos, path.length());
		std::cout << state << "\t" << dataInfo[f].labelTime << "\t";
		string z = "1";
		string replicate = to_string(dataInfo[f].replicate);

		std::cout << z << "\t" << replicate << "\t";

		std::cout << results[p].deutResults[f].massCenter;
		std::cout << std::endl;
		}
	}
}


