//---------------------------------------------------------------------------

#ifndef sortedClassesH
#define sortedClassesH
//---------------------------------------------------------------------------
#include "common.h"
#include <limits>

const int			MAX_THEO_CHARGE = 6;
const int			MAXFITNR = 5;

template<class T>
class SortedClass
{
private:
	std::string *name;
	T	*data;
	int *pCharge;
	int *pFit;
	size_t maxnr;
	

public:
	size_t	nr;
	// ===================================================================================================================
	// Constructor, set default values
	// ===================================================================================================================
	void Clear()
	{
		for (size_t i = 0; i < maxnr; i++)
		{
			data[i] = NULL;//numeric_limits<T>::quiet_NaN();
			pCharge[i] = numeric_limits<int>::quiet_NaN();
			pFit[i] = numeric_limits<int>::quiet_NaN();
			name[i].clear();
		}
		nr = 0;
	}

	SortedClass<T>(size_t n)
	{
		if (n <= 0) n = 1;
		data = new T[n];
		if (data == nullptr) SOFTWARE_ERROR;
		pCharge = new int[n];
		pFit = new int[n];
		name = new std::string[n];
		maxnr = n;
		Clear();
	}

	// ===================================================================================================================
	// ===================================================================================================================
	~SortedClass()	{ delete[] data; delete[] (pCharge); SAFE_DELETE_ARRAY(pFit); delete[] name; }

	// ===================================================================================================================
	// ===================================================================================================================
	T &operator [](size_t i)	{ return data[i]; }

	// ===================================================================================================================
	// ===================================================================================================================
	T	GetData(size_t i)	{ if (i < maxnr){ return data[i]; } else{ SOFTWARE_ERROR; } return NULL; };

	// ===================================================================================================================
	// ===================================================================================================================
	std::string GetName(size_t i){ if (i < maxnr){ return name[i]; }return ""; };

	// ===================================================================================================================
	// ===================================================================================================================
	int GetCharge(size_t i){ if (i < maxnr){ return pCharge[i]; }return ERROR_VALUE; };

	// ===================================================================================================================
	// ===================================================================================================================
	int GetFit(size_t i)	{ if (i < maxnr){ return pFit[i]; }return ERROR_VALUE; };

	void SetElem(size_t idx, T elem, int charge = ERROR_VALUE, int fit = ERROR_VALUE, std::string inName = "") //always store the biggest MAXNR elements
	{ //don't use this, because it acn break the sorting, only usefull for copying from another sorted source
		if (idx < maxnr)
		{
			data[idx] = elem;
			pCharge[idx] = charge;
			pFit[idx] = fit;
			name[idx] = inName;
			//_nr = min(idx,maxnr);
		}
	}

	void SetFit(size_t idx, int fit = ERROR_VALUE) //always store the biggest MAXNR elements
	{ //don't use this, because it acn break the sorting, only usefull for copying from another sorted source
		if (idx < maxnr)
		{
			pFit[idx] = fit;
		}
	}

	// ===================================================================================================================
	//	add a new element to the storage. order it ascending based on elem. If no more space, replace the smallest
	// ===================================================================================================================
	void AddElem(T elem, int charge = ERROR_VALUE, int fit = ERROR_VALUE, std::string inName = "")
	{
		if (nr == maxnr)
		{
			//SOFTWARE_ERROR;
			PushElem(elem, charge, fit, inName);
			return;
		}

		//~~~~~~~~~~~
		int found = nr;
		//~~~~~~~~~~~

		for (int i = 0; i < nr; i++)
		{
			if (elem < data[i])
			{
				found = i;
				for (int x = nr - 1; x >= found; x--)
				{
					data[x + 1] = data[x];
					pCharge[x + 1] = pCharge[x];
					pFit[x + 1] = pFit[x];
					name[x + 1] = name[x];
				}
				break;
			}
		}

		data[found] = elem;
		pCharge[found] = charge;
		pFit[found] = fit;
		name[found] = inName;

		nr++;
	};
private:
	// ===================================================================================================================
	//	Add a new element by replacing the smallest if no more storage available
	// ===================================================================================================================
	void PushElem(T elem, int charge = ERROR_VALUE, int fit = ERROR_VALUE, std::string inName = "") //always store the biggest MAXNR elements
	{
		//SOFTWARE_ERROR;
		//hibas, nem jol rakja be, az ADD sem
		// maradnak ures helyek
		//~~~~~~~~~~~
		size_t found = nr;
		//~~~~~~~~~~~
		for (size_t i = 0; i < nr; i++)
		{
			if (elem < data[i])   //should always give false at first
			{
				found = i;
				break;
			}
		}

		if (found == 0) return;

		//toljuk le a regieket
		for (size_t i = 1; i < found; i++)
		{
			data[i - 1] = data[i];
			pCharge[i - 1] = pCharge[i];
			pFit[i - 1] = pFit[i];
			name[i - 1] = name[i];
		}
		//berakni az uj elemet a regi ala
		data[found - 1] = elem;
		pCharge[found - 1] = charge;
		pFit[found - 1] = fit;
		name[found - 1] = inName;
	}
};

// Class to hold results from groupping
class FinalResult
{
public:
	double		finalIntensity;
	double		finalRT;
	double		finalSigma;
	int 		finalChargeIndex;
	int 		finalFitIndex;
	//scoring after chem filters
	double  	peakRTScore;

	void Copy(FinalResult* inPut);
	void Init();
	void ScoreInit();
	void Add(const FinalResult& f2);
	void Assign(const FinalResult& f2);
	void Set(double I, double RT, double sigma, int charge, int fit);
	void Score(double rts);

	FinalResult& operator= (const FinalResult& f2);
	FinalResult(FinalResult &inF);
	FinalResult();
};

class FinalResults
{
private:
	void SortResultsByI(FinalResult* arr, int left, int right);
public:
	FinalResult results[MAX_THEO_CHARGE*MAXFITNR];
	size_t nr;

	FinalResults();
	void Copy(FinalResults *inPut);

	// ===================================================================================================================
	// Clear values in the storage
	// ===================================================================================================================
	void Clear();

	void ClearScores();

	// ===================================================================================================================
	// Add new values to the storage
	// ===================================================================================================================
	int AddElem(double I, double RT, double sigma, int charge, int fit);

	// ===================================================================================================================
	// Delete a "row" from storage space, it will move the rest in the empty space
	// ===================================================================================================================
	void DeleteIndex(int idx);
	void DeleteZeroRT();

	void Sort(); //sort by Intensity

	void Swap(int id1, int id2); //swap 2 elements

};

#endif
