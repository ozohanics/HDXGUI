//$T indentinput.cpp GC 1.140 01/23/11 19:13:21

//$6
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//	Ebben a fileban levo fuggvenyek nagy resze a GlycoMiner eseti feladatainak a megoldasara keszult.
//	Nem emlekszem a legtobb esetben miert csinaltam oket, de tobbseguk nagyon egyszeru
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "common.h"


using namespace std;
using namespace std::chrono;

// =======================================================================================================================
//    search for value in input array  that falls close to elem. Return index
//	  lo is start index, up is end index in array
// =======================================================================================================================
int BinSearch(const double *arrayin, double elem, int lo, int up)
{
	int		last = -1;
	double	limit = 0.25;
	int		middle = (lo + up) / 2;
	int		upLimit = up - 1;
	int		loLimit = lo;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~

	while (lo <= up)
	{
		if (middle < loLimit)
		{
			break;
		}

		if (middle > upLimit)
		{
			break;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		double	delta = fabs(arrayin[middle] - elem);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		last = middle;

		if (delta < limit)
		{
			return middle;
		}
		else if (elem > arrayin[middle])
		{
			lo = middle + 1;
		}
		else if (elem < arrayin[middle])
		{
			up = middle - 1;
		}

		middle = (up + lo) / 2;
	}

	return last;
}

// =======================================================================================================================
//  Use binary search for int array
// =======================================================================================================================
int Search(const int *array, int elem, int lo, int up)
{
	//~~~~~~~~~~~~~~~~~~~~~~~
	int middle = (lo + up) / 2;
	int upLimit = up - 1;
	int loLimit = lo;
	//~~~~~~~~~~~~~~~~~~~~~~~

	while (lo <= up)
	{
		if (middle < loLimit)
		{
			break;
		}

		if (middle > upLimit)
		{
			break;
		}

		if (array[middle] == elem)
		{
			return middle;
		}
		else if (elem > array[middle])
		{
			lo = middle + 1;
		}
		else if (elem < array[middle])
		{
			up = middle - 1;
		}

		middle = (up + lo) / 2;
	}

	return -1;
}

// =======================================================================================================================
// =======================================================================================================================
int Search(const double *array, const double *intensity, double elem, int lo, int up, double limit, double precision)
{
	return BinarySearchBase(array, intensity, elem, lo, up, limit, precision, false);
}

// =======================================================================================================================
// =======================================================================================================================
int Search(const double *array, double elem, int lo, int up, double limit)
{
	return BinarySearchBase(array, NULL, elem, lo, up, NULL, limit, false);
}

// =======================================================================================================================
// =======================================================================================================================
int SearchPPM(const double *array, double elem, int lo, int up, double MAX_PPM)
{
	return BinarySearchBase(array, NULL, elem, lo, up, NULL, MAX_PPM, true);
}

// =======================================================================================================================
// =======================================================================================================================
int SearchPPM(const double *array, const double *intensity, double elem, int lo, int up, double limit, double MAX_PPM)
{
	return BinarySearchBase(array, intensity, elem, lo, up, limit, MAX_PPM, true);
}

// =======================================================================================================================
//    template func for binary search
// =======================================================================================================================
int BinarySearchBase
(
	const double	*array,
	const double	*intensity,
	double	elem,
	int		lo,
	int		up,
	double	limit,
	double	MAX_PPM,
	bool	isppm
)
{
	if (array == NULL)
	{
		return -1;
	}

	if (fabs(elem) < 1e-10)
	{
		return -1;
	}

	//~~~~~~~~~~~~~~~~~
	int upLimit = up - 1;
	int loLimit = lo;
	//~~~~~~~~~~~~~~~~~

	if (!isppm)
	{	// if it is not PPM based, calculate the correct PPM
		MAX_PPM = MAX_PPM / elem * 1.0e6;
	}

	while (lo <= up)
	{
		int middle = (up + lo) / 2;
		if (middle < loLimit)
		{
			break;
		}

		if (middle > upLimit)
		{
			break;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		double	delta = fabs(array[middle] - elem) / elem * 1.0e6;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		if (delta <= MAX_PPM)
		{
			if (intensity == NULL)
			{
				return middle;
			}
			else if (intensity[middle] > limit)
			{
				return middle;
			}
		}

		if (elem > array[middle])
		{
			lo = middle + 1;
		}
		else
		{
			up = middle - 1;
		}
	}

	return ERROR_VALUE;
}

//
// =======================================================================================================================
//    func for binary search, to find the closest value to a given one, return its position
// =======================================================================================================================
//
int BinarySearchClosest
(
	const double	*array,
	const double	*intensity,
	double	elem,
	int		lo,
	int		up,
	double	limit,
	double	MAX_PPM,
	bool	isppm
)
{
	if (array == NULL)
	{
		SOFTWARE_ERROR;
		return -1;
	}

	if (intensity == NULL)
	{
		SOFTWARE_ERROR;
		return -1;
	}


	if (fabs(elem) < 1e-10)
	{
		return -1;
	}

	//~~~~~~~~~~~~~~~~~
	int upLimit = up - 1;
	int loLimit = lo;
	int middle = (up + lo) / 2;
	//~~~~~~~~~~~~~~~~~

	if (!isppm)
	{	// if it is not PPM based, calculate the correct PPM
		MAX_PPM = MAX_PPM / elem * 1.0e6;
	}

	while (lo <= up)
	{
		middle = (up + lo) / 2;

		if (middle < loLimit)
		{
			break;
		}

		if (middle > upLimit)
		{
			break;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		double	delta = fabs(array[middle] - elem) / elem * 1.0e6;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		if (delta <= MAX_PPM)
		{
			if (intensity == NULL)
			{
				break; //return middle
			}
			else if (intensity[middle] > limit)
			{
				break; //return middle
			}
		}

		if (elem > array[middle])
		{
			lo = middle + 1;
		}
		else
		{
			up = middle - 1;
		}
	}

	//this part is not working well
	double delta1 = fabs(elem - array[max(loLimit, middle - 1)]) / elem * 1.0e6;; //check the distance to the one before it
	double delta2 = fabs(elem - array[min(upLimit, middle + 1)]) / elem * 1.0e6;; //check the distance to the one after it
	double delta3 = fabs(elem - array[middle]) / elem * 1.0e6;;

	if (delta3 > MAX_PPM) return ERROR_VALUE;

	int found_point = middle;

	if (delta1 < delta3)
	{
		if (delta1 < delta2)
		{		
			found_point =  max(loLimit, middle - 1);
		}
		else if (delta2 < delta3)
		{
			found_point = min(upLimit, middle + 1);
		}
	}
	else if (delta2 < delta3)
	{
		found_point = min(upLimit, middle + 1);
	}
	//else {	return middle;	}

	if (intensity == nullptr) return found_point;
	else if (intensity[found_point] > limit)
	{
		return found_point;
	}
	else return ERROR_VALUE;


	SOFTWARE_ERROR;
	return ERROR_VALUE;
}

//
//
// =======================================================================================================================
//    func for binary search, to find the closest value to a given one, return the difference
// =======================================================================================================================
//
double BinarySearchClosestDiff
(
	const double	*array,
	const double	*intensity,
	double	elem,
	int		lo,
	int		up,
	double	limit,
	double	MAX_PPM,
	bool	isppm
)
{
	if (array == NULL)
	{
		return -1;
	}

	if (fabs(elem) < 1e-10)
	{
		return -1;
	}

	//~~~~~~~~~~~~~~~~~
	int upLimit = up - 1;
	int loLimit = lo;
	int middle = (up + lo) / 2;
	//~~~~~~~~~~~~~~~~~

	if (!isppm)
	{	// if it is not PPM based, calculate the correct PPM
		MAX_PPM = MAX_PPM / elem * 1.0e6;
	}

	while (lo <= up)
	{
		middle = (up + lo) / 2;

		if (middle < loLimit)
		{
			break;
		}

		if (middle > upLimit)
		{
			break;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		double	delta = fabs(array[middle] - elem) / elem * 1.0e6;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		if (delta <= MAX_PPM)
		{
			if (intensity == NULL)
			{
				break; //return middle
			}
			else if (intensity[middle] > limit)
			{
				break; //return middle
			}
		}

		if (elem > array[middle])
		{
			lo = middle + 1;
		}
		else
		{
			up = middle - 1;
		}
	}

	double delta1 = fabs(elem - array[min(loLimit, middle - 1)]); //check the distance to the one before it
	double delta2 = fabs(elem - array[max(upLimit, middle + 1)]); //check the distance to the one after it
	double delta3 = fabs(elem - array[middle]);

	if (delta1 <= delta2)
	{
		if (delta1 <= delta3)
		{
			return delta1;
		}
		else
		{
			return delta3;
		}
	}
	else
	{
		if (delta2 <= delta3)
		{
			return delta2;
		}
		else
		{
			return delta3;
		}
	}

	//soha nem erhetne ide
	SOFTWARE_ERROR;
	return ERROR_VALUE;
}

// =======================================================================================================================
// =======================================================================================================================
int SearchLin(const double *array, double elem, int lo, int up)
{
	for (int i = lo; i < up; i++)
	{
		double delta = fabs(array[i] - elem);
		if (delta < MAX_DIFF)
		{
			return i;
		}
	}

	return -1;
}

// =======================================================================================================================
// =======================================================================================================================
int SearchLin(const int *array, int elem, int lo, int up)
{
	//~~
	int i;
	//~~

	for (i = lo; i < up; i++)
	{
		if (array[i] == elem)
		{
			return i;
		}
	}

	return -1;
}

// =======================================================================================================================
// =======================================================================================================================
int SearchLin(const double *array, double elem, int lo, int up, double limit)
{
	for (int i = lo; i < up; i++)
	{
		double delta = fabs(array[i] - elem);
		if (delta < limit)
		{
			return i;
		}
	}

	return -1;
}

int SearchLin(const int *arrayin, int elem, int up)
{
	for (int i = 0; i < up; i++)
	{
		if (arrayin[i] == elem)
		{
			return i;
		}
	}

	return ERROR_VALUE;
}

// =======================================================================================================================
// =======================================================================================================================
int SearchPPMLin(const double *array, double elem, int lo, int up, double limit)
{
	//V550 An odd precise comparison: elem == 0.0. It's probably better to use a comparison with defined precision: fabs(A - B) < Epsilon. common.cpp 512

	if (elem <= FLT_EPSILON) return -1;

	for (int i = lo; i < up; i++)
	{
		double delta = fabs(array[i] - elem) * 1.0e6 / elem;
		if (delta < limit)
		{
			return i;
		}
	}

	return -1;
}

// =======================================================================================================================
// =======================================================================================================================
std::string FormatDouble(double in, int nr)
{
	//~~~~~~~~~~~~~~~~~~~~~
	std::string	temp;
	char buffer[64];
	temp = "%." + std::to_string(nr) + "f";
#ifdef __GNUC__
	sprintf(buffer, temp.c_str(), in);
#else
	sprintf_s(buffer, sizeof(buffer), temp.c_str(), in);
#endif // __GNUC__
	return buffer;
}

// =======================================================================================================================
// =======================================================================================================================
std::string FormatDoubleSci(double in, int nr)
{
	//~~~~~~~~~~~~~~~~~~~~~
	std::string	temp;
	char buffer[64];
	//~~~~~~~~~~~~~~~~~~~~~
	temp = "%." + std::to_string(nr) + "E";
#ifdef __GNUC__
	sprintf(buffer, temp.c_str(), in);
#else
	sprintf_s(buffer, sizeof(buffer), temp.c_str(), in);
#endif // __GNUC__
	return buffer;
}

std::string FormatBool(bool in)
{
	if (in) return "true";
	return "false";
}

// =======================================================================================================================
// =======================================================================================================================
long GetFileSize(std::string filename)
{
	long fsize = ERROR_VALUE;
	try {
		std::ifstream in(filename, std::ifstream::in | std::ifstream::binary);
		in.seekg(0, std::ifstream::end);
		fsize = (int)in.tellg();
		in.close();
	}
	catch (...)
	{
		WRITE_DEBUG("Unable to open file: " + filename);
	}
	return fsize;
}

void OutputSpectrumArray(const double* x, const double* y, long nr)
{
	std::ofstream outfile;
	outfile.open(DebugPath, std::ios_base::app);
	
	for (int idx = 0; idx < nr; idx++)
	{
		outfile << x[idx] << "\t" << y[idx] << "\n";
	}

	outfile.close();
}

// =======================================================================================================================
//	read any file and return its text
// =======================================================================================================================
std::string ReadFiles(std::string fileName)
{
	if (!FileExists(fileName))
	{
		throw invalid_argument("File does not exist: " + fileName);
	}
	//use streams to read file
	std::ifstream inData(fileName.c_str(), std::ios::in | std::ios::binary);
	if (inData)
	{
		std::string contents;
		inData.seekg(0, std::ios::end);
		contents.resize((size_t)inData.tellg());
		inData.seekg(0, std::ios::beg);
		inData.read(&contents[0], contents.size());
		inData.close();
		return(contents);
	}
	WRITE_DEBUG("File not open " + fileName);
	throw(errno);

	return "";
}