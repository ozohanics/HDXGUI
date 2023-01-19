#include "HDXMain.h"

bool HDXMain::SearchDataFiles()
{
	TStringList splitList;
	splitList.separator = '\t';
	TStringList stateList;
	TStringList replicateList;

	for (size_t f = 0; f < fileNr; f++)
	{
		splitList.SetText(fileList.Item(f));
		std::string filePath = splitList.Item(RAWINFO::riFilePath);
		std::string tmp = splitList.Item(RAWINFO::riDeuTime);
		if (tmp.length() < 2) { break; }//error
		double deuTime = stod(tmp);
		std::string State = splitList.Item(RAWINFO::riState); 
		stateList.push_back(State);
		std::string Replicate = splitList.Item(RAWINFO::riReplicate);	
		replicateList.push_back(Replicate);

		if (!FileExists(filePath))
		{
			std::cout << "File does not exist: " << filePath << std::endl;
			continue;
		}
		//read first file
		rawData->XMLToMemory(filePath);
		//create storage for data, assuming all files share the same length
		if (f==0) peptideStore->CreateResultStorage(rawData, fileNr);
		//search one file for all peptides
		peptideStore->SetRawInfo(filePath, deuTime, f, State, stoi(Replicate));

		
		for (size_t p = 0; p < peptideNr; p++)
		{			
			//std::cout << fileList.Item(f) << "\t" << "Searching peptide: " << p << std::endl;
			if (deuTime < NUMERICAL_ZERO_DOUBLE)
			{
				peptideStore->SearchAPeptide(rawData, p, f, false);
				rawData->isRef = true;
			}
			else
			{
				peptideStore->SearchAPeptide(rawData, p, f, true);
			}
			peptideStore->CalcPeptideMassCenter(p, f);
		}
	}

	for (size_t p = 0; p < peptideNr; p++)
	{
		peptideStore->CalcReferenceMassCenter(p);
		for (size_t f = 0; f < fileNr; f++)
		{
			peptideStore->CalcDeuteration(p, f);
		}		
	}

	SaveResults();
	


	return false;
}

void HDXMain::SaveResults()
{
	auto ret = freopen("stdoutput.txt", "w", stdout);

	//peptideStore.WriteResultsToStdOut();
	peptideStore->WriteDeuterationDataToStdOut();
	peptideStore->WriteDeuterationPercentToStdOut();
	peptideStore->WriteDeuterationCentersToStdOut();
	auto ret2 = fclose(stdout);

	ret = freopen("MEMHDX_out.txt", "w", stdout);
	peptideStore->WriteDeuterationInMEMHDX();
	ret2 = fclose(stdout);
	
}

void HDXMain::SaveModifiedInput()
{
	auto ret = freopen("modified_pepfile.txt", "w", stdout);
	for (int p = 0; p < peptideNr; p++)
	{
		//eliminate deleted peptides
		if (peptideStore->results[p].deutResults[0].massCenter < 0) continue;

		std::cout << peptideStore->peptides[p].proteinID << "\t" << peptideStore->peptides[p].sequence << "\t" <<
			peptideStore->peptides[p].startRT << "\t" << peptideStore->peptides[p].endRT << "\t" <<
			peptideStore->peptides[p].z << std::endl;
	}
	auto ret2 = fclose(stdout);
}

bool HDXMain::DeletePeptide(int idx)
{
	//peptideNr--;
	return peptideStore->DeletePeptide(idx);
	
}

bool HDXMain::ReEvaluatePeptide(int p, int f, double RTs, double RTe, bool useXICOverlap, int charge)
{//analog with search, but only for one peptide
		//search one file for all peptides
		peptideStore->useXICOverlap = useXICOverlap;
		peptideStore->peptides[p].endRT = RTe;
		peptideStore->peptides[p].startRT = RTs;
		if (charge > 0)
		{
			peptideStore->peptides[p].z = charge;
		}
		if (last_file_idx != f)
		{
			rawData->XMLToMemory(peptideStore->dataInfo[f].rawPath);
			last_file_idx = f;
		}

		if (peptideStore->dataInfo[f].labelTime <= NUMERICAL_ZERO_DOUBLE)
		{
			//not deuterated
			peptideStore->SearchAPeptide(rawData, p, f, false);
			rawData->isRef = true;
		}
		else
		{
			//deuterated
			peptideStore->SearchAPeptide(rawData, p, f, true);
		}
		peptideStore->CalcPeptideMassCenter(p, f);
		
		for (size_t i = 0; i < fileNr; i++)
		{//recalc for all files
			peptideStore->CalcDeuteration(p, i);
		}
		//peptideStore->CalcDeuteration(p, f);
	
		return true;
}

bool HDXMain::ResetIsotope(int peptideIdx, int fileIdx, double isoMin, double isoMax, bool useMax)
{
	bool isoFound = false;
	for (int i = 0; i < peptideStore->results[peptideIdx].ICs[fileIdx].isotopeNr; i++)
	{
		double aMass = peptideStore->results[peptideIdx].ICs[fileIdx].MZ[i];
		if (aMass < isoMax && aMass > isoMin)
		{
			//mass found, set to zero 
			peptideStore->results[peptideIdx].ICs[fileIdx].intensity[i] = 0;
			isoFound = true;
			break;
		}	
	}
	peptideStore->CalcPeptideMassCenter(peptideIdx, fileIdx);
	peptideStore->CalcDeuteration(peptideIdx, fileIdx);
	double scanMaxMass = 0.0;
	double scanMaxInt = 0.0;
	for (int m = 0; m < peptideStore->results[peptideIdx].ICs[fileIdx].isotopeNr; m++)
	{
		if (scanMaxInt < peptideStore->results[peptideIdx].ICs[fileIdx].intensity[m])
		{
			scanMaxInt = peptideStore->results[peptideIdx].ICs[fileIdx].intensity[m];
			scanMaxMass = peptideStore->results[peptideIdx].ICs[fileIdx].MZ[m];
		}
	}
	if (last_file_idx != fileIdx)
	{
		last_file_idx = fileIdx;
		rawData->XMLToMemory(peptideStore->dataInfo[fileIdx].rawPath);
	}
	if (!useMax)
	{//use isotope maximum
		int maxidx = peptideStore->peptides[peptideIdx].composition.maxIdx;
		scanMaxMass = peptideStore->peptides[peptideIdx].composition.isoMass[maxidx] / peptideStore->peptides[peptideIdx].z + masses::Hp;
	}
	
	double sDelta = scanMaxMass * DEUTMASSERRORPPM / 1e6; //scan with DEUTMASSERRORPPM ppm error
	peptideStore->CreateXIC(peptideIdx, fileIdx, scanMaxMass, sDelta, rawData);

	return isoFound;
}

