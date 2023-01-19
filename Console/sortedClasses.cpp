//---------------------------------------------------------------------------

#include "sortedClasses.h"

//---------------------------------------------------------------------------

void FinalResults::SortResultsByI(FinalResult* arr, int left, int right)
{
	int i = left, j = right;
	double tmp;
	double pivot = arr[(left + right) / 2].finalIntensity;
	
	/* partition */
	while (i <= j) {
		while (arr[i].finalIntensity < pivot)
			i++;
		while (arr[j].finalIntensity > pivot)
			j--;
		if (i <= j) {
			tmp = arr[i].finalIntensity;
			arr[i].finalIntensity = arr[j].finalIntensity;
			arr[j].finalIntensity = tmp;
			i++;
			j--;
		}
	};

	/* recursion */
	if (left < j)
		SortResultsByI(arr, left, j);
	if (i < right)
		SortResultsByI(arr, i, right);
}

//constructor
FinalResults::FinalResults()
{
	nr = 0;
}

//copy constructor
void FinalResults::Copy(FinalResults *inPut)
{
	if (inPut == nullptr)
	{
		return;
	}
	nr = inPut->nr;
	for (int i = 0; i < inPut->nr; i++)
	{
		results[i].Copy(&inPut->results[i]);
	}
}

void FinalResult::Copy(FinalResult* inPut)
{
	finalIntensity = inPut->finalIntensity;
	finalRT = inPut->finalRT;
	finalSigma = inPut->finalSigma;

	finalChargeIndex = inPut->finalChargeIndex;
	finalFitIndex = inPut->finalFitIndex;

	peakRTScore = inPut->peakRTScore;
}

// =========================================================================
// Clear values in the storage
// =========================================================================
void FinalResults::Clear()
{
	if (nr >= MAX_THEO_CHARGE*MAXFITNR)
	{
		SOFTWARE_ERROR;
	}
	for (int i = 0; i < nr; i++)
	{
		results[i].Init();
	}
	nr = 0;
}

void FinalResults::ClearScores()
{
	for (int i = 0; i < nr; i++)
	{
		results[i].ScoreInit();
	}
}

// =========================================================================
// Add new values to the storage
// =========================================================================
int FinalResults::AddElem(double I, double RT, double sigma, int charge, int fit)
{
	if (nr >= MAX_THEO_CHARGE*MAXFITNR - 1)
	{
		SOFTWARE_ERROR;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	FinalResult  *temp = new FinalResult[nr + 1];
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	for (int i = 0; i < nr; i++)
	{
		temp[i] = results[i];
	}

	//~~~~~~~~~~~
	int found = -1;
	//~~~~~~~~~~~

	for (int i = nr - 1; i > -1; i--)
	{
		if (I > temp[i].finalIntensity)
		{
			found = i + 1;
			results[found].Set(I, RT, sigma, charge, fit);

			for (int x = found + 1; x < nr + 1; x++)
			{
				results[x] = temp[x - 1];
			}
			break;
		}
	}

	if (found == -1)
	{
		results[0].Set(I, RT, sigma, charge, fit);
		found = 0;
		for (int x = found + 1; x < nr + 1; x++)
		{
			results[x] = temp[x - 1];
		}
	}

	nr++;

	SAFE_DELETE_ARRAY(temp);

	return found;
};

// ===========================================================================
// Delete a "row" from storage space, it will move the rest in the empty space
// ===========================================================================
void FinalResults::DeleteIndex(int idx)
{
	if (idx > nr - 1 || idx < 0)
	{
		SOFTWARE_ERROR;
		return;
	}

	for (int i = idx; i < nr - 1; i++)
	{
		results[i] = results[i + 1];
	}

	nr--;
}

void FinalResults::DeleteZeroRT()
{
	for (int i = 0; i < nr; i++)
	{
		if (results[i].peakRTScore < FLT_EPSILON)
		{
			DeleteIndex(i);
			i--;
		}
	}
}

void FinalResults::Sort()
{//sort by intensity, by default
	if (nr <= 0) return;
	SortResultsByI(results, 0, nr-1);
}

void FinalResults::Swap(int id1, int id2)
{//swap 2 elements
	if (id1 < 0 || id1 > nr - 1 || id2 < 0 || id2 > nr - 1) SOFTWARE_ERROR;

	FinalResult tmp;
	tmp.Assign(results[id2]);
	results[id2].Assign(results[id1]);
	results[id1].Assign(tmp);
}

void FinalResult::Init()
{
	finalIntensity = 0;
	finalRT = 0;
	finalSigma = 0;
	finalChargeIndex = -1;
	finalFitIndex = -1;
	peakRTScore = 0;
}

void FinalResult::ScoreInit()
{
	peakRTScore = 0;
}

void FinalResult::Add(const FinalResult& f2)
{
	finalIntensity += f2.finalIntensity;
	finalRT += f2.finalRT;
	finalSigma += f2.finalSigma;
	peakRTScore = 0;
}

void FinalResult::Assign(const FinalResult& f2)
{
	finalIntensity = f2.finalIntensity;
	finalRT = f2.finalRT;
	finalSigma = f2.finalSigma;
	finalChargeIndex = f2.finalChargeIndex;
	finalFitIndex = f2.finalFitIndex;
	peakRTScore = f2.peakRTScore;
}

void FinalResult::Set(double I, double RT, double sigma, int charge, int fit)
{
	finalIntensity = I;
	finalRT = RT;
	finalSigma = sigma;
	finalChargeIndex = charge;
	finalFitIndex = fit;
}

void FinalResult::Score(double rts)
{
	peakRTScore = rts;
}

FinalResult& FinalResult::operator= (const FinalResult& f2)
{
	if (this == &f2) return *this;
	Assign(f2);
	return *this;
};

FinalResult::FinalResult()
{
	Init();
}

FinalResult::FinalResult(FinalResult &inF)
{
	Assign(inF);
}
