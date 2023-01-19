/*$T indentinput.cpp GC 1.140 04/28/16 09:27:59 */


/*$6
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */


#include "pch.h"
#include "TStringList.h"

/*
 =======================================================================================================================
 =======================================================================================================================
 */

std::vector<std::string> &TStringList::split(const std::string &s, char delim, std::vector<std::string> &elems)
{
	/*~~~~~~~~~~~~~~~~*/
	using namespace std;
	/*~~~~~~~~~~~~~~~~*/

	std::stringstream ss(s);

	/*~~~~~~~~~~~~~*/
	std::string item;
	/*~~~~~~~~~~~~~*/

	while(std::getline(ss, item, delim))
	{
		elems.push_back(item);
	}

	elems.pop_back();
	elems.push_back(item);

	return elems;
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
std::vector<std::string> TStringList::split(const std::string &s, char delim)
{
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	std::vector<std::string>	elems;
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	return split(s, delim, elems);
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
std::vector<std::string> TStringList::split(const std::string &s, const std::string &delim)
{
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	using namespace std;
	std::vector<std::string> elems;
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/*
	 * std::vector<std::string> copys;
	 * *set vector to first string
	 */
	elems.push_back(s);

	for(unsigned int i = 1; i < delim.length(); i++)
	{
		/*~~~~~~~~~~~~~~~~~~~~~~~~*/
		/* for all delimiters split string */
		int cur_size = elems.size();
		/*~~~~~~~~~~~~~~~~~~~~~~~~*/

		for(int ji = 0; ji < cur_size; ji++)
		{
			split(elems[ji], delim[i], elems);
		}

		/* erase strings that were not splitted */
		elems.erase(elems.begin(), elems.begin() + cur_size);
	}

	return elems;
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
int TStringList::LoadFromFile(std::string fileName)
{

	/* Error check */
	this->clear();
	if(!FileExists(fileName))
	{
		throw std::invalid_argument("File does not exist: " + fileName);
		return ERROR_VALUE;
	}

	/*~~~~~~~~~~~~*/
	/* use streams to read file */
	std::ifstream inData;
	/*~~~~~~~~~~~~*/

	inData.open(fileName.c_str());

	/*~~~~~~*/
	std::string ss;
	/*~~~~~~*/

	if((!inData.bad()) && (inData.is_open()))
	{
		inData.seekg(0, ifstream::end);

		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		size_t fsize = (size_t) inData.tellg();
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

		inData.seekg(0, ifstream::beg);
		clear();	/* clear previous contents */
		this->reserve(fsize + 1);

		/* read lines */
		while(getline(inData, ss, separator))
		{
			if(!inData.fail()) push_back(ss);
		}
	}
	else
	{
		throw invalid_argument("File could not be opened!" + fileName);
	}

	if(inData.is_open())
	{
		inData.close();
	}

	/* return size of list */
	return (int) size();
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
bool TStringList::SaveToFile(const std::string &fileName)
{
	/*~~~~~~~~~~~~~~~~~~~*/
	bool ret = false;
	std::ofstream	output;
	/*~~~~~~~~~~~~~~~~~~~*/

	output.open(fileName.c_str());
	if(output.good())
	{
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		std::vector<std::string>::iterator	it;
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

		for(it = begin(); it < end(); ++it)
		{
			output << it[0] << this->separator;
		}

		output.close();
		ret = true;
	}

	return ret;
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
int TStringList::SetText(const std::string &text)
{
	clear();
	if(text.size() <= 0) return 0;

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	std::vector<std::string>	elems;
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	elems = split(text, separator);

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	std::vector<std::string>::iterator	it;
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	it = elems.begin();
	assign(it, elems.end());
	return size();
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
void TStringList::Sort(bool ascending)
{
	if(ascending)
	{
		std::sort(begin(), end(), CompareStringsBool);
	}
	else
	{
		std::sort(begin(), end(), CompareStringsReverseBool);
	}
}

/*
 =======================================================================================================================
    add a new string to list end, return index
 =======================================================================================================================
 */
int TStringList::Add(const std::string &stringIn)
{
	push_back(stringIn);
	if(Sorted)
	{

		/* generate a sorted list */
		std::sort(this->begin(), this->end());
	}

	return size() - 1;
}

/*
 =======================================================================================================================
    add a new string to the specified place, return index
 =======================================================================================================================
 */
int TStringList::Insert(const std::string &stringIn, int index)
{
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	std::vector<std::string>::iterator	it;
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	it = begin() + index;
	insert(it, stringIn);
	return index;
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
int TStringList::IndexOf(const std::string &text)
{
	/*~~~~~~~~~~~~~~*/
#pragma warning(disable : 244)
	StringIterator it;
	/*~~~~~~~~~~~~~~*/

	it = std::find(begin(), end(), text);
	if(it != end())
	{
		return it - begin();
	}

	return ERROR_VALUE;
#pragma warning(default : 244)
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
std::string TStringList::Text()
{

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	/* return vector contents */
	std::string							tmp;
	std::vector<std::string>::iterator	it;
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	for(it = begin(); it < end(); ++it)
	{
		tmp += it[0];
	}

	return tmp;
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
std::string TStringList::Item(unsigned int idx) /* return string at position idx */
{
	if(idx < 0 || idx > this->size()) SOFTWARE_ERROR;

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	std::vector<std::string>::iterator	it = begin() + idx;
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	//return it[0];
	return this->data()[idx];
}

/*
 =======================================================================================================================
    split string by delim, and return idx part;
    1-based
 =======================================================================================================================
 */
std::string SplitString(std::string in, char delim, int idx)
{
	std::stringstream ss(in);

	/*~~~~~~~~~~~~~~~~~~*/
	std::string item = "";
	int count = 0;
	/*~~~~~~~~~~~~~~~~~~*/

	if(in.find(delim) == string::npos) return item;

	while(std::getline(ss, item, delim))
	{
		count++;
		if(count == idx)	/* found the desired part */
		{
			return item;
		}
	}

	/* if (count < idx) - if count too much just return last piece */
	return item;
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
std::string TStringList::NameFromIndex(unsigned int idx)
{	/* get data at index and split it according to "="
	 * ;
	 * return first part */
	if(idx < 0 || idx > this->size()) SOFTWARE_ERROR;

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	std::string in = *this[idx].data();
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	return SplitString(in, '=', 1);
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
std::string TStringList::ValueFromIndex(unsigned int idx)
{
	if(idx < 0 || idx > this->size()) SOFTWARE_ERROR;

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	std::string in = *this[idx].data();
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	return SplitString(in, '=', 2);
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
bool TStringList::IsName(const std::string &name)
{
	if(name.length() > 5) return true;
	return false;
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
std::string TStringList::ValueOfName(const std::string &name)	/* return value for first occurence of name */
{
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	std::vector<std::string>::iterator	it;
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	for(it = begin(); it < end(); ++it)
	{
		if(SplitString(it[0], '=', 1) == name)
		{
			return SplitString(it[0], '=', 2);
		}
	}

	SOFTWARE_ERROR return "";	/* nothing found */
}
