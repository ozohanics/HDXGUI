#pragma once
#include "isotopeCalc.h"
#include "common.h"

/*
	This class holds the data of one peptide
	it contains isotopic profiles as well as search results?
*/
class TPeptideData
{
public:
	//isotopic pattern with base isotopes
	TComposition	composition;
		//isoBase; 
	//isotopic pattern with max number of deuterium 
	//TIsotopeData	isoDeut;
	//peptide sequence
	std::string		sequence;
	//protein from which the peptide comes from
	std::string		proteinID;
	//peptide mass 
	double			mass;
	//peptide expected retention time window
	double			startRT;
	double			endRT;
	//info about the place of the sequence position in the protein
	int				seqStart;
	int				seqEnd;

	//peptide charge
	int				z;
	//max deuterium uptake possible
	int				maxD;
	//location of maximum int 
	long			maxIntLoc;

	TPeptideData& operator=(const TPeptideData& other)
	{
		if (this == &other)
			return *this;
		this->composition = other.composition;
		this->sequence = other.sequence;
		this->proteinID = other.proteinID;
		this->mass = other.mass;
		this->startRT = other.startRT;
		this->endRT = other.endRT;
		this->seqStart = other.seqStart;
		this->seqEnd = other.seqEnd;
		this->z = other.z;
		this->maxD = other.maxD;
		this->maxIntLoc = other.maxIntLoc;
		return *this;
	}

	TPeptideData() {
		sequence = ""; proteinID = ""; mass = ERROR_VALUE; startRT = ERROR_VALUE; endRT = ERROR_VALUE; z = ERROR_VALUE; maxD = ERROR_VALUE; seqStart
			= ERROR_VALUE; seqEnd = seqStart; maxIntLoc = ERROR_VALUE;
	}
	~TPeptideData() {}
	//debug
	void PrintComposition() { composition.PrintCompositionData(); }
	void PrintData() { std::cout << proteinID << "\t" << sequence << "\t" << mass << "\t" << startRT << "\t" << endRT << "\t" << z << std::endl; }
	int  CalcMaxDeuteration();
};

