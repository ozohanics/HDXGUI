//$T indentinput.h GC 1.140 01/23/11 19:15:42

//$6
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//

#ifndef peptideH
#define peptideH

//#include "common.h"
//#include "TProject.h"
#include "isotopecalc.h"
const int LOSSLENGTH = 10;


class	MyList
{
	//
	// -----------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------
	//
private:
	std::string	Text;

	//
	// -----------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------
	//
public:
	std::string											List[100];
	int													nr;

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//

	MyList() : Text("")
	{
		nr = 0;  // List=new std::string[100];
	}

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	~MyList() {
		try {
			nr = 0;
		}
		catch (...) {} // nothing to delete
	}

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	void												Clear() {
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		int	  i; for (i = 0; i < nr; i++) { List[i].clear(); }Text.clear(); nr = 0;
	}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	void												Add(const std::string &str) { if (nr < 100) { List[nr] = str; nr++; } }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	void												AddtoIndex(const std::string &str, int index) { if (index < 100) { List[index] = str; nr = index + 1; } }

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	int GetIndex(const std::string &Name, int ini)
	{
		std::string compStr = Name + "=";
		for (int i = 0; i < nr; i++)
		{			
			if (List[i].compare(0, compStr.length(), compStr) == 0 && ini == 1)
			{
				return i;
			}
			else if ((List[i] == Name) && (ini == 0))
			{
				return i;
			}
		}

		return -1;
	}

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	std::string GetValue(const std::string &Name)
	{
		std::string compStr = Name + "=";
		for (int i = 0; i < nr; i++)
		{
			if (List[i].compare(0, compStr.length(), compStr) == 0)
			{
				return(List[i].substr(Name.size() + 2, List[i].size() - Name.size()));
			}
		}

		return "";
	}

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	void SetValue(const std::string &Name, const std::string &Value)
	{
		int i = GetIndex(Name, 1);
		if ((i > -1) && (i < 100))
		{
			List[i] = Name + "=" + Value;
		}
		else
		{
			Add(Name + "=" + Value);
		}
	}

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	std::string						GetText() {
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		int	  i; Text.clear(); for (i = 0; i < nr; i++) { Text = Text + List[i] + "\r\n"; }return Text;
	}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	void							Delete(int index) { if (index < 100) { List[index].clear(); } }
};

// Peptide Class aminosavak sz?ma
static const int		MAX_AMINO_ACID = 24;
static const int		MAX_PEPTIDE_LENGTH = 1500;

// aminosavak indexel?se
static const std::string AMINO_ACID = "GAVLIMFWPSTCNQYHDEKRcX_";   //include modified aminoacids

// atomt?megek
static const double		MT_H = 1.00782, MT_O = 15.9949;

static const double		MT_ACID[MAX_AMINO_ACID] =
{
	57.02146,
	71.03711,
	99.06841,
	113.08406,
	113.08406,
	131.04049,
	147.06841,
	186.07931,
	97.05276,
	87.03203,
	101.04768,
	103.00919,
	114.04293,
	128.05858,
	163.06333,
	137.05891,
	115.02694,
	129.04259,
	128.09496,
	156.10111,
	103.00919 + masses::CARBAMIDOMETH,
	103.00919 - 1.00785 //X, Cys-Cys bonded peptide
};

static const int COMP_ACID_C[MAX_AMINO_ACID] =
{
	2,
	3,
	5,
	6,
	6,
	5,
	9,
	11,
	5,
	3,
	4,
	3,
	4,
	5,
	9,
	6,
	4,
	5,
	6,
	6,
	5,
	3
};
static const int COMP_ACID_H[MAX_AMINO_ACID] =
{
	3,
	5,
	9,
	11,
	11,
	9,
	9,
	10,
	7,
	5,
	7,
	5,
	6,
	8,
	9,
	7,
	5,
	7,
	12,
	12,
	8,
	4
};
static const int COMP_ACID_N[MAX_AMINO_ACID] =
{
	1,
	1,
	1,
	1,
	1,
	1,
	1,
	2,
	1,
	1,
	1,
	1,
	2,
	2,
	1,
	3,
	1,
	1,
	2,
	4,
	2,
	1
};
static const int COMP_ACID_O[MAX_AMINO_ACID] =
{
	1,
	1,
	1,
	1,
	1,
	1,
	1,
	1,
	1,
	2,
	2,
	1,
	2,
	2,
	2,
	1,
	3,
	3,
	1,
	1,
	2,
	1
};
static const int COMP_ACID_S[MAX_AMINO_ACID] =
{
	0,
	0,
	0,
	0,
	0,
	1,
	0,
	0,
	0,
	0,
	0,
	1,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	1,
	1
};

class					TPeptide
{	// calculate mass of a peptide, will include modifications too
	int			sequence[MAX_PEPTIDE_LENGTH];	// aa. sequence
	size_t		length; // length of peptide

	double		mass;	// mass of peptide
	int			charge; // charge of peptide
	std::string	MySequence;

	
	//
	// ===================================================================================================================
	// ===================================================================================================================
	//

	void calculate_sequence(const std::string &asequence)
	{
		for (int x = 0; x < MAX_PEPTIDE_LENGTH; x++)
		{
			sequence[x] = -1;
		}
		for (int s = 0; s < IsoType::atNr; s++)
		{
			composition[s] = 0;
		}

		length = asequence.size();
		for (int i = 0; i < length; i++)
		{
			int found = -1;
			char aaVal = asequence[i];
			found = (int) AMINO_ACID.find(aaVal);
			/*
			// 			if (aaVal == 'c')
			// 			{
			// 				int xx = 0;
			// 			}
			*/
			if (found > -1)
			{
				sequence[i] = found;
				composition[IsoType::atC] += COMP_ACID_C[found];
				composition[IsoType::atH] += COMP_ACID_H[found];
				composition[IsoType::atN] += COMP_ACID_N[found];
				composition[IsoType::atO] += COMP_ACID_O[found];
				composition[IsoType::atS] += COMP_ACID_S[found];
			}
			else
			{
				if (asequence[i] == 'J')
				{//Leu or Ile in sequence
					sequence[i] = 4;
				}
				else
				{
					WRITE_DEBUG("No valid peptide sequence found! Please make sure that no B, X, Z is used in fasta file!");
					SOFTWARE_ERROR;
				}
			}
		}
	}

	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	void calculate_mass()
	{
		mass = MT_H;
		for (int i = 0; i < length; i++)
		{
			try
			{
				mass += MT_ACID[sequence[i]];
			}
			catch (...)
			{
				WRITE_DEBUG("Calculating peptide mass failed!");
				SOFTWARE_ERROR;
			}
		}

		mass += (MT_O + MT_H);
		composition[IsoType::atH] += 2;
		composition[IsoType::atO] += 1;
	}

	//
	// -----------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------
	//
public:
	//
	// ===================================================================================================================
	// ===================================================================================================================
	//
	int 	composition[IsoType::atNr];

	TPeptide() : MySequence("")
	{//constructor
		
		length = 0; // length of peptide
		mass = 0;	// mass of peptide
		charge = 0; // charge of peptide


		//V730 Not all members of a class are initialized inside the constructor. Consider inspecting: sequence, composition. peptide.h 498
		for (int i = 0; i < MAX_PEPTIDE_LENGTH; i++)
		{
			sequence[i] = ERROR_VALUE;
		}

		for (int i = 0; i < IsoType::atNr; i++)
		{
			composition[i] = ERROR_VALUE;
		}

	}

	
	double	GetMass(const std::string &seq) { MySequence = seq; calculate_sequence(seq); calculate_mass(); return mass; }
};

// end of Peptide class


#endif
