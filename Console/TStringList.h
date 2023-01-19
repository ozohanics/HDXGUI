/*$T indentinput.h GC 1.140 04/28/16 09:27:04 */


/*$6
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */


#pragma once
#include <vector>
#include <list>
#include <string>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/stat.h>

const int							ERROR_VALUE = -1;
using namespace std;


#ifndef SOFTWARE_ERROR

#define SOFTWARE_ERROR			{std::ofstream outfile; outfile.open("c:\\", std::ios_base::app); \
	outfile << "Error in File " << __FILE__ << " at Line: "	<< std::to_string(__LINE__) << ", In Function: " << __FUNCTION__ ; outfile.close(); \
	throw std::runtime_error(std::string("\n\tin File: ") + __FILE__ + ",\n\tat Line: " + std::to_string(__LINE__) + ",\n\tIn Function: " + std::string(__FUNCTION__)); }
#endif // !SOFTWARE_ERROR

typedef std::vector<std::string>			StringContainer;
typedef StringContainer::iterator	StringIterator;

/*
 =======================================================================================================================
 =======================================================================================================================
 */

inline bool FileExists(const std::string &name)
{
	/*~~~~~~~~~~~~~~~*/
	struct _stat64 buffer;
	/*~~~~~~~~~~~~~~~*/

	return( _stat64(name.c_str(), &buffer) == 0);
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
inline int CompareStrings(const std::string &a, const std::string &b)
{
	return strcmp(a.c_str(), b.c_str());
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
inline int CompareStringsReverse(const std::string &a, const std::string &b)
{
	return -1 * strcmp(a.c_str(), b.c_str());
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
inline bool CompareStringsBool(const std::string &a, const std::string &b)
{
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	int ret = strcmp(a.c_str(), b.c_str());
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	return(ret >= 0) ? false : true;
}	/* compare two strings, return false if equal or greater */

/*
 =======================================================================================================================
 =======================================================================================================================
 */
inline bool CompareStringsReverseBool(const std::string &a, const std::string &b)
{
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	int ret = strcmp(a.c_str(), b.c_str());
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	return(ret >= 0) ? true : false;
}	/* compare two strings, return true if equal or greater */

/*
 * split string by delim, and return idx part;
 * 1-based
 */
std::string SplitString(std::string in, char delim, int idx);

/*
 =======================================================================================================================
    Error message definitions ;
    ! Class Name: TStringList Description: This class will emulate the TStringList class of VCL, to be used when
    translating projects from borland c++ to portable c++. Only STL types should be used
 =======================================================================================================================
 */
class TStringList :
	public std::vector <std::string>
{
/*
 -----------------------------------------------------------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------------------
 */
public:
	/*
	 ===================================================================================================================
	 ===================================================================================================================
	 */
	TStringList(void)	{ separator = '\n'; /* read new content line by line */ Sorted = false; }

	/* TStringList(void){} */
	char		separator;

	/* add a new string to list end, return index */
	int			Add(const std::string &stringIn);

	/* add a new string to the specified place, return index */
	int			Insert(const std::string &stringIn, int index);

	/* load data into the list, return size */
	int			LoadFromFile(std::string fileName);

	bool		SaveToFile(const std::string &fileName);

	/* set vector content to text. split in the process. return new size */
	int			SetText(const std::string &text);

	/*
	 * items are considered as Name=Value pairs ;
	 * Get Name part
	 */
	std::string Item(unsigned int idx);
	std::string NameFromIndex(unsigned int idx);

	/* Get Value part */
	std::string ValueOfName(const std::string &name);	/* return value for first occurence of name */
	std::string ValueFromIndex(unsigned int idx);

	/* search for a stringh in vector */
	int			IndexOf(const std::string &text);

	/* sort vector */
	void		Sort(bool ascending = true);
	bool		Sorted;

	std::string Text(); /* get stringlist content */

	/*
	 ===================================================================================================================
	 ===================================================================================================================
	 */
	size_t		ItemCount()		{ return size(); }

/*
 -----------------------------------------------------------------------------------------------------------------------
    return the number of items
 -----------------------------------------------------------------------------------------------------------------------
 */
private:

	/* vector based functions to split a string. Highly effective */
	std::vector<std::string>	&split(const std::string & s, char delim, std::vector < std::string > &elems);
	std::vector<std::string>	split(const std::string &s, char delim);
	std::vector<std::string>	split(const std::string &s, const std::string &delim);

	bool						IsName(const std::string &name);
};
