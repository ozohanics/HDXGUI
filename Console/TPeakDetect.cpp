#include "TPeakDetect.h"


TPeakDetect::TPeakDetect()
{
	mz = nullptr;
	DataY = nullptr;
	fderiv = nullptr;
	f2deriv = nullptr;
	fpeak = nullptr;
	pointNr = 0;
	amplitudeThr = 10;
	slopeThr = 0.5;
	peaknr = 0;
	//V730 Not all members of a class are initialized inside the constructor. Consider inspecting: peakIdx. tpeakdetect.cpp 4
	for (int i = 0; i < MAX_IDX_PEAKS; i++)
	{
		peakIdx[i] = ERROR_VALUE;
		peakStart[i] = ERROR_VALUE;
		peakEnd[i] = ERROR_VALUE;
		peakMaxInt[i] = ERROR_VALUE;
	}
	
}


TPeakDetect::~TPeakDetect()
{
	ClearData();
}

void TPeakDetect::ReadInput(string fileName)
{
	if (FileExists(fileName))
	{
		fileIO.LoadFromFile(fileName);
		ClearData();
		pointNr = fileIO.ItemCount();
		mz = new double[pointNr];
		DataY = new double[pointNr];

		for (int i = 0; i < pointNr; i++)
		{
			mz[i] = 0;
			DataY[i] = 0;
			string tmp = GetPart(fileIO[i],1);
			SafeConvert(tmp, &mz[i]);
			tmp = GetPart(fileIO[i], 2);
			SafeConvert(tmp, &DataY[i]);				
		}
	}
}

void TPeakDetect::SetData(double* x, double* yin, int nr)
{
	ClearData();
	if (x != nullptr && yin != nullptr && nr > 0)
	{
		pointNr = nr;
		mz = new double[pointNr];
		DataY = new double[pointNr];
		for (int i = 0; i < nr; i++)
		{
			mz[i] = x[i];
			//DataY[i] = yin[i];
		}

		//do a very basic smoothing
		DoSmoothWithOutput(yin, DataY, nr, 11, 1);
	}
}

void TPeakDetect::CalculateDeriv()
{
	if (pointNr <= 0) return;
	int nr = 9;
	int dist[9] = { -4, -3, -2, -1, 0, 1, 2, 3, 4 };
	int factor[9] = { -3, -3, -2, -1, 0, 1, 2, 3, 3 };

	SAFE_DELETE_ARRAY(fderiv);
	fderiv = new double[pointNr];

	for (int i = 0; i < pointNr; i++)
	{
		if ((i < (nr - 1) / 2 + 1) || (i>pointNr-(nr-1)/2))
		{
			fderiv[i] = 0;
			continue;
		}

		double sum = 0;
		for (int n = 0; n < nr; n++)
		{
			sum += DataY[i + dist[n]] * factor[n];
		}
		fderiv[i] = sum;
	}

}

void TPeakDetect::ClearData()
{
	SAFE_DELETE_ARRAY(mz);
	SAFE_DELETE_ARRAY(DataY);
	SAFE_DELETE_ARRAY(fderiv);
	SAFE_DELETE_ARRAY(f2deriv);
	SAFE_DELETE_ARRAY(fpeak);
	pointNr = 0;
	peaknr = 0;
}

//TPeakDetect tp;
//tp.SetData(output.ICs[i].RT[n], smoothed_intensity, array_size);
//tp.CalculateDeriv();
//tp.DetectPeaks();

void TPeakDetect::DetectPeaks()
{
	/*Signal greater than amplitude threshold?	Down zero crossing in first derivative?	Crossing slope greater than slope threshold?*/
	if (pointNr <= 0) return;

	CalculateDeriv();

	SAFE_DELETE_ARRAY(fpeak);
	fpeak = new int [pointNr];

	fpeak[0] = 0;
	peaknr = 0;
	peakStart[0] = 0;

	for (int i = 1; i < pointNr; i++)
	{
		if (peaknr == 99) break; //too many peaks
		bool signalGood = (DataY[i] > amplitudeThr); //magassag, ne nezzuk a <0 csucokat
		bool fderivGood = (Sign(fderiv[i - 1]) > Sign(fderiv[i])); //irany valtas
		bool slopeGood = ((fderiv[i - 1] - fderiv[i]) > slopeThr); //ne osszon nullaval
		
		if (signalGood && fderivGood && slopeGood)
		{
			fpeak[i] = fpeak[i - 1]++;
			peakIdx[peaknr] = i;
			
			for (int s = i-1; s >= 1; s--)
			{	//ha a derivalt eleg kicsi, akkor azt kinevezzuk kezdetnek
				if ((Sign(fderiv[s -1]) != Sign(fderiv[s])))
				{
					peakStart[peaknr] = s;
					break;
				}
			}
			for (int e = i + 1; e < pointNr-1; e++)
			{	//ha a derivalt eleg kicsi, akkor azt kinevezzuk vegnek
				if ((Sign(fderiv[e + 1]) != Sign(fderiv[e])))
				{
					peakEnd[peaknr] = e;
					break;
				}
			}

			peaknr++;
		}
		else
		{
			fpeak[i] = fpeak[i - 1];
		}
	}
	
	//*max_element(arr, arr + n);
	for (int n = 0; n < peaknr; n++)
	{
		if (peakStart[n] < 0) continue;
		if (peakEnd[n] < 0) continue;
		if (peakStart[n] > peakEnd[n]) continue;
		peakMaxInt[n] = *max_element(&DataY[peakStart[n]], &DataY[peakEnd[n]]);
	}
		
}

void TPeakDetect::WriteData(string fileName)
{
	if (peaknr < 1 || pointNr < 1) return;
	
	TStringList save;
	save.Add("Number of peaks: " + to_string(peaknr) + "\n");
	

	for (int n = 0; n < peaknr; n++)
	{
		save.Add(FormatDouble(peakIdx[n], 2) + "\t" + FormatDouble(mz[peakIdx[n]], 2) + "\t" + FormatDouble(peakStart[n], 2) + "\t" + FormatDouble(peakEnd[n], 2) + "\t" + FormatDouble(peakMaxInt[n],2));
		save.Add("\n");
	}
	for (int i = 0; i < pointNr; i++)
	{
		save.Add(FormatDouble(mz[i], 4) + "\t" + FormatDouble(DataY[i], 4) + "\t" + FormatDouble(fderiv[i], 4) );
	}
	save.SaveToFile(fileName);

}

int TPeakDetect::GetPeakIDX(int nr)
{
	if (nr >= 0 && nr < peaknr)
	{
		return peakIdx[nr];
	}
	SOFTWARE_ERROR;
	return -1;
}

void TPeakDetect::SumPeaks(TPeakDetect& tpIn)
{
	//first copy peaks
	//std::cout << tpIn.peaknr <<"_|_"<< this->peaknr << std::endl;
	if (this->peaknr < 1)
	{//copy data
		for (int i=0; i<tpIn.peaknr; i++)
		{
			this->peakEnd[i] = tpIn.peakEnd[i];
			this->peakStart[i] = tpIn.peakStart[i];
			this->peakIdx[i] = tpIn.GetPeakIDX(i);
			this->peakMaxInt[i] = tpIn.peakMaxInt[i];
		}
		this->peaknr = tpIn.peaknr;
		return;
	}
	
	for (int t = 0; t < this->peaknr; t++)
	{//for all present peaks
		bool found = false;
		for (int i = 0; i < tpIn.peaknr; i++)
		{//check if overlaps are present
			bool overlap1 = tpIn.peakStart[i] <= this->peakStart[t] && tpIn.peakEnd[i] <= this->peakEnd[t] && tpIn.peakEnd[i] >= this->peakStart[t];
			bool overlap2 = tpIn.peakStart[i] >= this->peakStart[t] && tpIn.peakStart[i] <= this->peakEnd[t] && tpIn.peakEnd[i] >= this->peakEnd[t];
			bool overlap3 = tpIn.peakStart[i] <= this->peakStart[t] && tpIn.peakEnd[i] >= this->peakEnd[t];
			bool overlap4 = tpIn.peakStart[i] >= this->peakStart[t] && tpIn.peakEnd[i] <= this->peakEnd[t];
			if (overlap3 || overlap2 || overlap1 || overlap4)
			{//peaks overlap, do something
				this->peakStart[t] = max(this->peakStart[t], tpIn.peakStart[i]);
				this->peakEnd[t] = min(tpIn.peakEnd[i], this->peakEnd[t]);
				found = true;
				break;
			}
		}
		if (!found)
		{

			this->peakStart[t] = 0;
			this->peakEnd[t] = 0;
		}
	}

	int skipCount = 0;
	int finalPeakNr = peaknr;
	for (int i = 0; i < this->peaknr; i++)
	{
		skipCount = i;
		if (this->peakEnd[i] == 0)
		{
			finalPeakNr--;
			for	(int s=i+1; s<peaknr; s++)
			{
				if (this->peakEnd[s] == 0)
				{ 
					finalPeakNr--;
					skipCount = s;
					break;
				}
			}
			this->peakEnd[i] = peakEnd[skipCount];
			this->peakStart[i] = peakStart[skipCount];
			this->peakIdx[i] = GetPeakIDX(skipCount);
			this->peakMaxInt[i] = peakMaxInt[skipCount];
			i = skipCount;			
		}
	}
	peaknr = finalPeakNr;

}


