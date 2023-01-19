//$T indentinput.cpp GC 1.140 01/23/11 19:14:44

//$6
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//

#include "isotopeCalc.h"
#include "globals.h"
#include "constants.h"
// =======================================================================================================================
//	Fill IsoCalc with results from IPC calc method
// =======================================================================================================================
void TIsoCalc::StoreIPCResults()
{

	string tmpstr = "";

	for (unsigned int i = 0; i < IsoType::atNr; i++)
	{
		tmpstr += isoNames[i] + to_string(composition[i]);
	}

	def.sp = ms::computeIsotopePatternThr(tmpstr, 0.00001);
	def.sp.sortByMass();

	isoDataNr = 0;
	double mz_resolution = def.sp.masses[0] / MAXMASSRESOLUTION; //use 50000 resolution

	double yVal = -1;
	maxIdx = -1;

	for (int i = 0; i < def.sp.size(); ++i)
	{
		i = CalcIsoPoint(i, mz_resolution);
		isoMass[isoDataNr] = iPoint.iMass;
		isoData[isoDataNr] = iPoint.iAbun;

		isoDataNr++;
		//if ( isoDataNr > 5) break;
	}

	CalcMaxIntIsoIdx();
}

int TIsoCalc::CalcIsoPoint(int idx, double resolution)
{
	int lastIdx = idx;

	for (int i = idx + 1; i < def.sp.size(); ++i)
	{
		if (def.sp.masses[i] - def.sp.masses[i - 1] < resolution)
		{
			//set last idx;
			lastIdx = i;
		}
		else
		{
			break;
		}
	}

	iPoint.iAbun = 0;
	iPoint.iMass = 0;
	for (int i = idx; i <= lastIdx; ++i)
	{
		iPoint.iAbun += def.sp.intensities[i];
		iPoint.iMass += def.sp.masses[i] * def.sp.intensities[i];
	}

	iPoint.iMass = iPoint.iMass / iPoint.iAbun;
	return lastIdx;
}

// =======================================================================================================================
//	Calculate isotopic cluster for input Mass
// =======================================================================================================================
int TIsoCalc::CalcIsoCluster(double Mass)
{
	ClearCluster();
	clusterStart = 0;
	clusterEnd = 1;

	StoreIPCResults();
	/*string st = "dbg" + to_string(Mass) + to_string(0) + ".txt";
	auto ret = freopen(st.c_str(), "w", stdout);
	std::cout << "Data follows" << std::endl;
	PrintIsoCluster();
	fclose(stdout);
	*/
	return isoDataNr;
}

void TIsoCalc::PrintIsoCluster()
{
	//just a debug function, to see if results are correct
	TStringList *clusterList = new TStringList;
	PrintCompositionData();

	string tmp = "";
	for (int i = 0; i < IsoType::atNr; i++)
	{
		tmp += isoNames[i] + to_string(composition[i]);
	}
	clusterList->Add(tmp);

	for (int i = 0; i < isoDataNr; i++)
	{
		cout << FormatDouble(isoData[i], 4) + "\t" + FormatDouble(isoMass[i], 4) << endl;
		clusterList->Add(FormatDouble(isoData[i], 4) + "\t" + FormatDouble(isoMass[i], 4));
	}
	clusterList->SaveToFile(FormatDouble(isoMass[0], 4) + ".txt");//"tmpIso.txt");
}

void TIsoCalc::ClearCluster()
{
	clusterStart = 0;
	clusterEnd = 0;
	for (int i = 0; i < MAXCALCISO; i++)
	{
		IsoCluster[i] = 0.0;
	}
	isoDataNr = 0;
}

void TIsoCalc::SetFormula(int* atoms, int atomNr)
{
	//set formula to calc with
	if (atoms == nullptr)
	{
		SOFTWARE_ERROR;
		return;
	}
	SetComposition(atoms, atomNr);

	nAtom = IsoType::atNr;
	//valid = true;
}

void TIsoCalc::SetFormula(int atC, int atH, int atN, int atO, int atS, int atP)
{
	//set formula to calc with
	composition[IsoType::atC] = atC;
	composition[IsoType::atH] = atH;
	composition[IsoType::atN] = atN;
	composition[IsoType::atO] = atO;
	composition[IsoType::atS] = atS;
	composition[IsoType::atP] = atP;
	nAtom = IsoType::atNr;
	valid = true;
}


void TIsoCalc::SetDefaults()
{
	//set default abundances using only
	nAtom = -1;
	// set abundances and atom numbers
	composition[0] = 1;
	for (int i = 1; i < ATOMTYPES; i++)
	{
		composition[i] = 0;
	}
	// set isotopic abundance for ATOM 1 (glcnac)
	Abund[0][0] = 100;
	Abund[0][1] = 9.6509;
	Abund[0][2] = 1.4176;

	// set the number of isotopes
	nPeak[0] = 3;
	//valid = false;
}

void TIsoCalc::SetFormula(TCompositionData* inputComp) //specify composition
{
	InputComposition(inputComp);
	nAtom = IsoType::atNr;
}

TIsoCalc::TIsoCalc()
{
	EraseArray(IsoCluster, MAXCALCISO);
	SetDefaults();
	ClearCluster();
	//UseAtoms();
}

void TIsoCalc::EraseArray(double *arr, int nr)
{
	for (int i = 0; i < nr; i++){ arr[i] = 0; }
}

void TComposition::InputComposition(const TCompositionData* ind)
{
	if (ind == nullptr) return;
	SetComposition(ind->composition, ATOMTYPES);

	//valid = ind->isValid();
}

void TComposition::InputData(const TComposition* ind)
{
	if (ind == nullptr) return;
	for (int i = 0; i < ATOMTYPES; i++)
	{
		composition[i] = ind->composition[i];
	}
	for (int i = 0; i < MAXMOLISONR; i++)
	{
		isoData[i] = ind->isoData[i];
		isoMass[i] = ind->isoMass[i];
	}

	isoDataNr = ind->isoDataNr;
	CalcMaxIntIsoIdx();
	valid = ind->valid;
}

void TComposition::OutputData(TComposition* out)
{
	if (out == nullptr) return;
	for (int i = 0; i < ATOMTYPES; i++)
	{
		out->composition[i] = composition[i];
	}
	for (int i = 0; i < MAXMOLISONR; i++)
	{
		out->isoData[i] = isoData[i];
		out->isoMass[i] = isoMass[i];
	}
	out->isoDataNr = isoDataNr;
	CalcMaxIntIsoIdx();
	out->valid = valid;
}

void TComposition::CalcMaxIntIsoIdx()
{
	double yVal = -1;
	maxIdx = -1;
	for (int i = 0; i < isoDataNr; i++)
	{
		if (yVal < isoData[i])
		{
			maxIdx = i;
			yVal = isoData[i];
		}
	}
}

TComposition::TComposition()
{
	for (int i = 0; i < MAXMOLISONR; i++)
	{
		isoData[i] = 0.0;
		isoMass[i] = 0.0;
	}
	isoDataNr = -1;
	maxIdx = -1;
}

double TComposition::CalcCompDelta(TComposition *incomp)
{ //scale from 0 to 1?
	if (incomp == nullptr) return 1e3;
	double sum = 0;
	for (int c = 0; c < incomp->isoDataNr; c++)
	{
		sum += (isoData[c] - incomp->isoData[c])*(isoData[c] - incomp->isoData[c]);
	}
	return sqrt(sum);
}

void TComposition::NormalizeData()
{
	void CalcMaxIntIsoIdx();
	for (int i = 0; i < isoDataNr; i++)
	{
		isoData[i] = isoData[i] / isoData[maxIdx];
	}
}

void TComposition::CopyComposition(TComposition *incomp)
{
	for (int i = 0; i < MAXMOLISONR; i++)
	{
		isoData[i] = incomp->isoData[i];//isotopic distribution for current composition
		isoMass[i] = incomp->isoMass[i];//mass for current composition
	}

	isoDataNr = incomp->isoDataNr; //number of isotopes bigger than ISOPREC
	maxIdx = incomp->maxIdx;
	CopyCompositionData(incomp);
}

void TCompositionData::SetComposition(const int* atoms, int atomNr)
{
	if (atoms != nullptr)
	{
		valid = false;
		for (int i = 0; i < min(atomNr, ATOMTYPES); i++)
		{
			composition[i] = atoms[i];
			if (atoms[i] > 0 && atoms[i] < MAXATOMPERMOL)
			{
				valid = true;
			}
		}
	}
}

void TCompositionData::ValidateComposition()
{
	valid = false;
	int sum = 0;
	for (int t = 0; t < IsoType::atNr; t++)
	{
		sum += composition[t];
	}
	if (sum > 0 && sum < MAXATOMPERMOL) valid = true;
}

void TCompositionData::Clear()
{
	valid = false;
	for (int i = 0; i < ATOMTYPES; i++)
	{
		composition[i] = 0;
	}
}

TCompositionData::TCompositionData()
{
	Clear();
}

void TCompositionData::PrintCompositionData()
{
	string tmp = "";
	for (int i = 0; i < IsoType::atNr; i++)
	{
		tmp += isoNames[i] + to_string(composition[i]);
	}
	cout << tmp << endl;
}

std::string TCompositionData::CompositionToString()
{
	string tmp = "";
	for (int i = 0; i < IsoType::atNr; i++)
	{
		tmp += isoNames[i] + to_string(composition[i]);
	}
	return tmp;
}

void TCompositionData::TextToComposition(string compText) //use data stored in text format CxHyNzOw ...
{
	Clear();
	for (int i = 0; i < IsoType::atNr; i++)
	{
		composition[i] = GetNumberAfterLetter(compText, isoNames[i]);
	}
	valid = true;
}

void TCompositionData::CopyCompositionData(TCompositionData *acomp)
{
	Clear();
	//WRITE_INFO("Copy composition");
	for (int i = 0; i < IsoType::atNr; i++)
	{
		composition[i] = acomp->composition[i];
		//WRITE_INFO(isoNames[i] + " " + to_string(composition[i]));
	}
	valid = acomp->valid;
}

void TCompositionData::AddComposition(TCompositionData *acomp)
{
	//WRITE_INFO("Add composition");
	for (int i = 0; i < IsoType::atNr; i++)
	{
		composition[i] += acomp->composition[i];
		//WRITE_INFO(isoNames[i] + " " + to_string(composition[i])+" " + to_string(acomp->composition[i]));
	}
	valid = true;
}
