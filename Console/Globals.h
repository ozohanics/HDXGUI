/*$T indentinput.h GC 1.140 05/09/14 21:56:58 */

/*$6
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/
#pragma once

#ifndef GLOBALSH
#define GLOBALSH

#include <string.h>
#include <ctype.h>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <functional>
#include <cctype>
#include <atomic>
#include <fstream>
#include <memory>

#include "TStringList.h"

#include "constants.h"
#include "ms_isocalc.hpp"



const std::string DebugPath = "debugLog.txt";
const std::string DebugInfo = "debugInfo.txt";

#define WRITE_INFO(x)			{ std::ofstream outfile; outfile.open(DebugInfo, std::ios_base::app); \
	outfile << x << "\n"; outfile.close(); }


#define WRITE_DEBUG(x)			{ std::ofstream outfile; outfile.open(DebugPath, std::ios_base::app); \
	outfile << x << " in File " << __FILE__ << " at Line: " << std::to_string(__LINE__) << ", In Function: " << __FUNCTION__ << "\n"; outfile.close(); }

#ifndef SOFTWARE_ERROR

#define SOFTWARE_ERROR			{std::ofstream outfile; outfile.open(DebugPath, std::ios_base::app); \
	outfile << "Error in File " << __FILE__ << " at Line: "	<< std::to_string(__LINE__) << ", In Function: " << __FUNCTION__ ; outfile.close(); \
	throw std::runtime_error(std::string("\n\tin File: ") + __FILE__ + ",\n\tat Line: " + std::to_string(__LINE__) + ",\n\tIn Function: " + std::string(__FUNCTION__)); }

#endif // !SOFTWARE_ERROR

#define SAFE_DELETE(x)			{ if(x != nullptr) { delete x; x = nullptr; } }
#define SAFE_DELETE_ARRAY(x)	{ if(x != nullptr) { delete [] x; x = nullptr; } }
#define SAFE_DELETE_ARRAYS(x)	{ if(x != nullptr) { delete [] x; x = nullptr; } }

//void					SAFE_DELETE_ARRAYS(std::string *x);

template <typename T>  T tmax(T a, T b) { return (((a) > (b)) ? (a) : (b)); }
template <typename T>  T tmin(T a, T b) { return (((a) < (b)) ? (a) : (b)); }


#include <stdio.h>  /* defines FILENAME_MAX */
#include <direct.h>
#define GetCurrentDir _getcwd

const double DIV_MIN = sqrt(FLT_MIN); //minimum number for which you might want a division operation

/*
#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif
*/
#define USE_DIGITS 5


enum ftMSMSResultsTable
{
	ftmCounter,
	ftmCharge,
	ftmMass,
	ftmIntensity,
	ftmRetentionTime,
	ftmQualityScore,
	ftmMarkerScore,
	ftmGlycoScore,
	ftmGlycan,
	ftmPeptide,
	ftmPPM,
	ftmEnd
};


typedef struct isotope{
	struct isotope *next, *previous;
	double p;
	double mass;
	isotope(): p(0)
	{
		next = nullptr; previous = nullptr; p = 0.0; mass = 0.0;
	}
	//~isotope() { SAFE_DELETE(next); SAFE_DELETE(previous); }
} isotope;

typedef struct compound {
	struct compound *next, *previous;
	short amount;
	struct isotope *isotopes;
	compound(): amount(0)
	{ next = nullptr; previous = nullptr; isotopes = nullptr; }
	//~compound() { SAFE_DELETE(next); SAFE_DELETE(previous); SAFE_DELETE(isotopes); }
} compound;

typedef struct element {
	struct element *next, *previous;
	char *symbol;
	struct isotope *isotopes;
	element(): symbol(nullptr)
	{ next = nullptr; previous = nullptr; isotopes = nullptr; }
	//	~element() { SAFE_DELETE(next); SAFE_DELETE(previous); SAFE_DELETE(isotopes); }
} element;

/*Extern parameters can move to own class*/
class TCompounds
{
public:
	compound *verbindung;
	isotope *peaks;
	element *elements;
	//isotope calc param
	int fast_calc;
	double FWHM_limit;
	char atC, atH, atO, atN, atS, atP;

	TCompounds()
	{
		verbindung = nullptr;
		peaks = nullptr;
		elements = nullptr;
		//isotope calc param
		fast_calc = 80;
		FWHM_limit = 0.02; //ezeket olvassza ossze
		atC = 'C';
		atH = 'H';
		atO = 'O';
		atN = 'N';
		atS = 'S';
		atP = 'P';
	}
	~TCompounds()
	{
		SAFE_DELETE(peaks);
		SAFE_DELETE(elements);
		SAFE_DELETE(verbindung);
	}

	void ReinitCalc()
	{
		SAFE_DELETE(peaks);
		SAFE_DELETE(verbindung);
	};

	int pars_peptid(char *formel);
	int pars_chem_form(char *formel, int len);

	int is_symbol(char *probe);
	int add_component(char *symbol, int number);
	int pars_amino_acid(char *formel);
	int add_amino_acid(char acid);
	int print_sum(void);

	int init_elements(void);
	FILE* open_file(char *filename);
	void add_element(element* ce);
	void add_isotope(isotope* ci, element* ce);
	element* parse_element(FILE* data, char* linebuffer);
	isotope* parse_isotope(FILE* data, char* linebuffer);

	int calculate_peaks(void);
	void print_result();
	void free_list(isotope *target);
	void summarize_peaks(void);
	isotope *add_peak(isotope *base, isotope *peak);
	void cut_peaks(isotope *spectrum);

	int StorePeaks(double* massOut, double* intOut);
};

//#define MIN_DIF 0.0009
#define MIN_INT 0.009

#define MAX_LINE (80)
#define NUM_ELEMENTS (86)
#define ELEMENTFILE  "elemente"

#define MAX_DIGITS 10
#define MAX_PEP_LINE 81

std::string GetExePath();

// trim from start
std::string &ltrim(std::string &s);
// trim from end
std::string &rtrim(std::string &s);
// trim from both ends
std::string &trim(std::string &s);

void		SafeConvert(std::string aText, double *out);
void		SafeConvert(std::string aText, double *out, double min, double max);

class		Globals
{
	/*
	-----------------------------------------------------------------------------------------------------------------------
	-----------------------------------------------------------------------------------------------------------------------
	*/
public:
	std::string GetPath()
	{
		char cCurrentPath[FILENAME_MAX];

		if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
		{
			return "";
		}
		return cCurrentPath;
	}

	void Init()
	{
		passcount = 0;

		//V730 It is possible that not all members of a class are initialized inside the constructor. Consider inspecting:
		//current_glycan_masses, MAXCHRPEAKS, MAX_SIGMA_DIFF, MAX_NEIGHBOUR_PEAKS. globals.h 164

		MIN_PEAK_SIGMA = 1e-2;
		MAX_PEAK_SIGMA = 3.5;
		MINREGIONHEIGHT = 5;
		MINGOODPEAKS = 4;
		SN_MASS_WINDOW = 2.5;
		ISO_MASS_WINDOW = 5.5;
		BASELINE_PERC = 0.5; //median
		BASELINE_MIN_POINTS = 15; //don't calc baseline if less than nr of point
		MIN_POINTS_IN_WINDOW = 5;
		SN_POINT_NR = 5;
		DELTARRT = 6;
		DELTART = 0.5;
		SEARCHWINDOWSIZE = 2.5;
		RTSCOREWINDOW = 0.15;
		RTSCOREVALUE = 0.3333333;
		GOODSCOREVAL = 10.0;
		GOODSCOREINCREMENT = 1.0;
		PENALTY_WEIGHT_NR = 0;
		PENALTY_WEIGHT_VAL = nullptr;


		GAUSS_FIT_SIGMA = 3.5;

		progressCount = -1;

		OVERLAP_FACTOR = 0.05;
		PREBASEISO_PERC = 0.5;

		tc.init_elements();

	
	}
	Globals(void): AppPath("c:\\"), DebugPath("debugLog.txt"), FUNCIDX(".idx"),	FUNCFILE(".DAT")
	{
		//WRITE_DEBUG(GetPath());
		AppPath = GetPath();
		DebugPath = AppPath + "\\debugLog.txt";
		IniPath = AppPath + "\\patternIni.ini";
		Init();
	}

	~Globals(void) { SAFE_DELETE_ARRAY(PENALTY_WEIGHT_VAL); }

	int		progressCount; //value to check for operations in a thread

	std::string AppPath;
	std::string DebugPath;
	std::string IniPath;

	double		*current_glycan_masses;

	int			passcount;

	/* constants */

	std::string FUNCIDX;
	std::string FUNCFILE;

	std::string DEF_RAW_FILE;
	std::string DEF_GLC_FILE;
	std::string DEF_PEP_FILE;
	std::string DEF_LIST_FILE;

	double		MIN_PEAK_SIGMA;
	double		MAX_PEAK_SIGMA;
	double		SN_MASS_WINDOW;			/* mass window for baseline calc */
	double		ISO_MASS_WINDOW;		/* mass window for an isotopic distribution */
	double		BASELINE_PERC;			/*median limit for baseline, percentage of points giving the baseline*/
	int			BASELINE_MIN_POINTS; //don't calc baseline if less than nr of point
	int			MIN_POINTS_IN_WINDOW;
	int			SN_POINT_NR;

	//int			MAX_NEIGHBOUR_PEAKS;	/* the peaks to check if there is anything similar to searched mass */

	double		RTSCOREWINDOW;			/* minutes */
	double		RTSCOREVALUE;			/* 3 windows to 0, number of windows to consider */

	double		GOODSCOREVAL;			/* what is the score recieved manually */
	double		GOODSCOREINCREMENT;		/* what is the score recieved manually */

	/* VT */
	double		MINREGIONHEIGHT;		/* minimum intenzitas egy regioban */
	int			MINGOODPEAKS;			/* how many points to request in a chrom peak minimum */

	/* ToDo: make this setting available on the UI */
	int			DELTARRT;			/* how close can two peaks be in the RRT table, for their places to be uncertain, 1
									* Fuc unit delta is 4 */
	double		DELTART;

	double		SEARCHWINDOWSIZE;	/* minutes */

	int			PENALTY_WEIGHT_NR;
	double		*PENALTY_WEIGHT_VAL;

	double		GAUSS_FIT_SIGMA;

	double		OVERLAP_FACTOR; //factor to multiply the calculated overlap value, for comparison with input
	double		PREBASEISO_PERC; //how big can the peak of same charge state be, before the base isotope

	TCompounds  tc; //compund types for isotope calc
	ms::Spectrum sp;

};

std::string GetPart(std::string text, int first);
bool		IsDelimiter(std::string in, std::string test, int pos);
int			GetNumberAfterLetter(std::string inputText, std::string Letter);
std::string ReplaceSubStr(std::string source, const std::string &substr, const std::string &replstr);
std::string GetDateTimeStr(); //create a string of the current date and time for individual file name
int Sign(double in); //get sign

extern Globals def;

#endif
