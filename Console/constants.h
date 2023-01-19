#pragma once
#ifndef CONSTANTSH
#define CONSTANTSH

//DONT_TOUCH_NEXT_LINE
const std::string Version = "0.1_b00";

//# define LOG_DEBUG

union ZERO
{
	double dt;
	int it;
	void* pt;
	char ct;
	ZERO() { it = 0; dt = it; pt = nullptr; ct = 0; }
};

const ZERO ZERO_CONST;

const std::string NumberDelimiter = "0123456789";
const std::string ExtendedNumberDelimiter = "0123456789,.-";

const int		PrgBar_UPDATE_INTERVAL = 200;
const double	PrgBar_UPDATE_STEP = 0.05;

const int				MAX_MATCH = 100;
const double			MAX_DIFF = 0.15;		// in Da OR PPM
const double			MAX_MARKER_PREC = 150;	// only in PPM
const double			PPM_DIFF = 15;
const double			MAX_DELTA = 0.25;		// maximum diff in Da between the match and search

const  double		MIN_GLC = 120.0f;		// min glc mass
//static double			MIN_MASS_LIMIT = 300;
const int				MAX_SERIES = 100;

const double		INSTR_ACCURACY = 0.0001;
const double		PPM_ACC_GAIN = 10.0;

const int		ION_SERES_WEIGHT = 10;
const double		ION_SERIES[] = { 0.0, 0.1, 0.5, 0.7, 3.0, 5.0, 11.0, 13.0, 17.0, 19.0 };
const int		ION_SERIES_NR = 9;

const int 		NO_DATA_PASSED = -1234; //for functions that have predefined value parameters

const int		MAX_RESULT_DATA = 10;
const int		MAX_CANDIDATE_PEAK_NR = 10;

const double	NOISE_SERIES[] = { 0.0, 1.0, 5.0, 15.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0, 200.0, 250.0, 300.0 };
const int	NOISE_SERIES_NR = 14;

const double	MARKER_WEIGHT = 1;
const int	MARKER_SIGNIF = 2;
const int	MARKER_NORM = 9;

const int	MAXDIGPEPTIDENR = 10000;

const int	MAXPEAKSPERSPECTRUM = 100000;

const int	 numel = 8;	// GLC_SEQ.Length();
const int    ANTENNANR = 6;

const int	HEADER_LINES = 4;
const int	PATH_LINE = 0;
const int	NR_LINE = 2;

const double MINRT = 1;// minimum retention time considered valid
const double MINLEN = 2; //minimum length for glycan struct
const int MAXRTDATA = 30;

const int MAXISOSTORED = 10; //number of isotopes stored for each compound we look for

const int HONR = 6; //number of O in Hex
const int NONR = 6; //number of O in HexNAc
const int FONR = 4; //number of O in Fuc-1

const std::string NAMESTR = "NAME";
const std::string RRTSTR = "RRT";
const std::string MASSSTR = "MASS";

const int MINFILENAMELEN = 4;

const double NUMERICAL_ZERO_DOUBLE = 0.000000000001; //anything smaller we consider 0

#ifndef MAX_PATH
const int	MAX_PATH = FILENAME_MAX;
#endif

//defaults neeeded everywhere + INI parameters

const int PARAM_FILE_PADDING_LENGTH = sizeof(double) * 24 - sizeof(bool);;

const char dec_point = std::use_facet< std::numpunct<char> >(std::cout.getloc()).decimal_point();
#endif
