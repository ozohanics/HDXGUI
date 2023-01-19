//$T indentinput.h GC 1.140 01/23/11 19:13:27

//$6
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//

#pragma once

#ifndef commonH
#define commonH

#include <math.h>
#include <float.h>
#include <iostream>
#include <stdio.h>
#include <string>
#include "Globals.h"
#include "masses.h"
#include "constants.h"

#include <algorithm>
#include <functional>
#include <cctype>
#include "time.h"
//#include <locale>

#ifdef __GNUC__ //gives error if func not in header
// =======================================================================================================================
//    replaces all substr in source with replstr
// =======================================================================================================================
#endif
#include <chrono>

//#define DEBUGISOTOPES 1

long GetTimeMs64();


int			BinSearch(const double *array, double elem, int lo, int up);
int			Search(const double array[], const  double *intesity, double elem, int lo, int up, double limit, double precision);
int			Search(const double *array, double elem, int lo, int up, double limit);
// =======================================================================================================================
// =======================================================================================================================
int			SearchLin(const double *array, double elem, int lo, int up);
int			SearchLin(const double *array, double elem, int lo, int up, double limit);
int			SearchLin(const int *arrayin, int elem, int up);
// =======================================================================================================================
// =======================================================================================================================

int			SearchPPMLin(const double *array, double elem, int lo, int up, double MAX_PPM);
int			SearchPPM(const double *array, double elem, int lo, int up, double MAX_PPM);
int			SearchPPM(const double *array, const double *intensity, double elem, int lo, int up, double limit, double MAX_PPM);

int			BinarySearchBase(const double	*array, const double	*intensity, double	elem, int		lo, int		up, double	limit, double	MAX_PPM, bool	isppm = false);
double		BinarySearchClosestDiff(const double	*array, const double	*intensity, double	elem, int		lo, int		up, double	limit, double	MAX_PPM, bool	isppm);
int			BinarySearchClosest(const double	*array, const double	*intensity, double	elem, int		lo, int		up, double	limit, double	MAX_PPM, bool isppm = false);

std::string	FormatDouble(double in, int nr);
std::string FormatDoubleSci(double in, int nr);
std::string FormatBool(bool in);

std::string	ReadFiles(std::string filename);
long		GetFileSize(std::string filename);

void		OutputSpectrumArray(const double* x, const double* y, long nr);
#endif
