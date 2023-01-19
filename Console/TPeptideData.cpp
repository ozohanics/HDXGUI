#include "TPeptideData.h"

int TPeptideData::CalcMaxDeuteration()
{
	//count number of prolines
	int proNr = std::count(sequence.begin(), sequence.end(), 'P');
	maxD = sequence.length() - proNr - 1;
	return maxD;
}
