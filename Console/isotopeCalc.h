//$T indentinput.h GC 1.140 01/23/11 19:14:48

//$6
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	This unit contains the TIsoCalc class, that can be used to calculate the mass and isotopic pattern for molecules up to
//	1000 C, H, N, O, S atoms.
//	2014-07-10: v2.0
//
//	Usage example for TIsoCalc
//		TIsoCalc isos;
//			//CHNOS - YGGFL, leucine enkephaline
//		isos.SetFormula(28, 38, 5, 7, 0);
//		isos.CalcIsoCluster(); //calculates the data, but stores result in isoData/isoMass arrays, isoDataNr - number of isotopes
//		isos.PrintIsoCluster();
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//

#ifndef isotopeCalcH
#define isotopeCalcH

#include "common.h"

enum IsoType
{//the different atoms considered CHNOS, up to 5 isotopes
	atC=0,
	atH=1,
	atN=2,
	atO=3,
	atS=4,
	atP=5,
	atNr=6
};

const string isoNames[IsoType::atNr] = { "C", "H", "N", "O", "S", "P" };

enum IsoNr
{
	baseIso,
	firstIso,
	secondIso,
	thirdIso,
	fourthIso,
	fifthIso,
	endIso
};

const int		ATOMTYPES = IsoType::atNr; //max 8 atom types can be used
const int		MAXISOTOPESFORATOM = IsoNr::endIso;
const int		MAXCALCISO = 250; 
const double	ISOPREC = 1; // min isotope abundance to calc, only OK for deuteration
const int		MAXATOMPERMOL = 1000; //maximum number of atoms in a molecule, we still want to use
//ne foglalkozzunk a 0.1-2%nal kisebb csucsokkal

const int MAXISOTYPES = 15; //max number of isotopes an atom type has
const int MAXMOLISONR = MAXISOTYPES * 3; //max number of isotopes a molecule has

const int IsoAbundanceNr[IsoType::atNr] = { 2, 2, 2, 3, 5, 1}; //nr of isotopes per atom type

const double MAXMASSRESOLUTION = 50000.0; //MS resolution, for isotope calc

class TCompositionData
{
protected:
	bool valid;
public:
	int composition[ATOMTYPES];//number of each atom type
	void SetComposition(const int* atoms, int atomNr);
	//check if there is a valid composition stored
	void ValidateComposition();

	void Clear();
	bool isValid(){ return valid; }
	TCompositionData();
	~TCompositionData(){};

	void PrintCompositionData();

	std::string CompositionToString();

	void TextToComposition(string compText);//use data stored in text format CxHyNzOw ...;

	void CopyCompositionData(TCompositionData *acomp);
	void AddComposition(TCompositionData *acomp);

	void DebugDataOut() {
		//write debug information to file

		WRITE_INFO(CompositionToString());
		WRITE_INFO("\n");
	}
};

class TComposition : public TCompositionData
{
public:	//TComposition

	double isoData[MAXMOLISONR];//isotopic distribution for current composition
	double isoMass[MAXMOLISONR];//mass for current composition
	int isoDataNr; //number of isotopes bigger than ISOPREC
	int maxIdx;

	TComposition();

	void CopyComposition(TComposition *incomp);

	void InputComposition(const TCompositionData* ind);

	void InputData(const TComposition* ind);

	void OutputData(TComposition* out);

	void CalcMaxIntIsoIdx(); //get the index of the biggest isotopic peak

	double CalcCompDelta(TComposition *incomp); //compare current composition to another

	void NormalizeData();
};

class			TIsoCalc :public TComposition
{
	//
	// -----------------------------------------------------------------------------------------------------------------------
	//    specific isotope cluster calculator for glycopeptide isotope abundance based on Kubinyi 1991
	//		Uses GlcNAc as base atom, and generates all isotopic patters like they are build up of this atom.
	// -----------------------------------------------------------------------------------------------------------------------
	//
private:
	double	Abund[ATOMTYPES][MAXISOTOPESFORATOM];
	int		nPeak[ATOMTYPES];
	int		clusterStart; //from which isotope do we begin filling the cluster
	int		clusterEnd; //how many isotopes we use

	double	IsoCluster[MAXCALCISO];
	int		nAtom;			// number of atoms max to consider < ATOMTYPES
	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	struct IsoPoint
	{
		double iMass = 0.0;
		double iAbun = 0.0;
	} iPoint;

	void	EraseArray(double *arr, int nr);

	//
	// -----------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------
	//
	

	//sum up isotopic intensities and return last index summed up
	int CalcIsoPoint(int idx, double resolution);


public:

	void SetDefaults(); //set default values

	void SetFormula(int atC=0, int atH=0, int atN=0, int atO=0, int atS=0, int atP=0); //specify composition
	void SetFormula(int* atoms, int atomNr); //specify composition
	void SetFormula(TCompositionData* inputComp); //specify composition

	void StoreIPCResults();

	void ClearCluster(); //clear previous data

	void PrintIsoCluster();

	TIsoCalc();

	int CalcIsoCluster(double Mass = -1);
};
#endif
