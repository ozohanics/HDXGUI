#include "pch.h"
#include "Globals.h"
#include <locale>
#include <iostream>
#include <string>

static inline std::string& ltrim(std::string& s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int c) {return !std::isspace(c); }));
	return s;
}


std::string& rtrim(std::string& s)
{
	s.erase(std::find_if(s.rbegin(), s.rend(), [](int c) {return !std::isspace(c); }).base(), s.end());
	//s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(pointer_to_unary_functionn<int, int>(std::isspace))).base(), s.end());
	return s;
}

std::string & trim(std::string &s)
{
	return ltrim(rtrim(s));
}

//
// =======================================================================================================================
// =======================================================================================================================
//
void SafeConvert(std::string aText, double *out)
{
	if (out == nullptr)
	{
		SOFTWARE_ERROR;
	}

	//numpunct
	aText = trim(aText);
	if (aText.length() < 1) return;
	char point = std::use_facet< std::numpunct<char> >(std::cout.getloc()).decimal_point();
	if ('.' != point) aText = ReplaceSubStr(aText, ".", to_string(point));

	try
	{
		*out = stod(aText);
	}
	catch (...)
	{
		SOFTWARE_ERROR;
	}
}

//
// =======================================================================================================================
// =======================================================================================================================
//
void SafeConvert(std::string aText, double *out, double min, double max)
{
	if (out == nullptr)
	{
		SOFTWARE_ERROR;
	}

	//numpunct
	char point = std::use_facet< std::numpunct<char> >(std::cout.getloc()).decimal_point();
	aText = trim(aText);
	if ('.' != point) aText = ReplaceSubStr(aText, ".", to_string(point));
	if (aText.length() == 0) return;

	try
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~
		double	val = stod(aText);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~

		if (min <= val && val <= max)
		{
			*out = val;
		}
	}
	catch (...)
	{
		SOFTWARE_ERROR;
	}
}



#ifdef __GNUC__
//do nothing
#else
// =======================================================================================================================
//    replaces all substr in source with replstr
// =======================================================================================================================
std::string ReplaceSubStr(std::string source, const std::string &substr, const std::string &replstr)
{
	int i = source.find(substr);
	if ((i != 0) && (substr != replstr))
	{
		do
		{
			source.replace(i, substr.size(), "");
			source.insert(i, replstr);
			i = source.find(substr);
		} while (i > 0);
	}

	return source;
}
#endif
// =======================================================================================================================
// =======================================================================================================================
std::string GetPart(std::string text, int first)
{  //return parts of delimited text by \t
	TStringList split;
	split.separator = '\t';
	int nr = split.SetText(text);

	if (nr < first - 1 || first < 1)
	{
		WRITE_DEBUG("GetPart function called with invalid parameter!")
			return "";
	}

	return split[first - 1];
}

//check if there is a character from test string in the input
bool IsDelimiter(std::string in, std::string test, int pos)
{
	for (unsigned int i = 0; i < test.length(); i++)
	{
		int currpos = in.find_first_of(test.substr(i, 1), pos);
		if (currpos == pos) return true;
	}

	return false;
}


int GetNumberAfterLetter(std::string inputText, std::string Letter) //return 0 if not found
{
	int LetterPos = inputText.find(Letter);
	int NumberAfter = 0;

	if (LetterPos < 0)
	{
		//cout << "Letter not found" << endl;
		WRITE_DEBUG("Letter not found");
		return NumberAfter;
	}
	std::string tmp = "";

	const std::string NumberDelimiter = "-1234567890";

	for (unsigned int i = LetterPos; i < inputText.length(); i++)
	{
		if (IsDelimiter(inputText, NumberDelimiter, i) == true)
		{
			cout << "delimiter found" << endl;
			int len = 1;
			for (unsigned int x = i + 1; x < inputText.length(); x++)
			{
				if (IsDelimiter(inputText, NumberDelimiter, x) == false || x == inputText.length())
				{
					len = x - i;
					break;
				}
			}
			tmp = inputText.substr(LetterPos + 1, len);
			break;
		}
	}
	try
	{
		NumberAfter = stoi(tmp);
	}
	catch (...)
	{
		//cout << "No number after letter" << endl;
		WRITE_DEBUG("No number after letter! ");
		SOFTWARE_ERROR;
	}

	return  NumberAfter;
}

//create a string of the current date and time for individual file name
std::string GetDateTimeStr()
{
	std::string out = "";
	struct tm newtime;
	char am_pm[] = "AM";
	__time64_t long_time;
	char timebuf[26];
	errno_t err;

	// Get time as 64-bit integer.
	_time64(&long_time);
	// Convert to local time.
	err = _localtime64_s(&newtime, &long_time);
	if (err)
	{
		printf("Invalid argument to _localtime64_s.");
		exit(1);
	}
	if (newtime.tm_hour > 12)        // Set up extension.
		strcpy_s(am_pm, sizeof(am_pm), "PM");
	if (newtime.tm_hour > 12)        // Convert from 24-hour
		newtime.tm_hour -= 12;        // to 12-hour clock.
	if (newtime.tm_hour == 0)        // Set hour to 12 if midnight.
		newtime.tm_hour = 12;

	// Convert to an ASCII representation.
	err = asctime_s(timebuf, 26, &newtime);
	if (err)
	{
		printf("Invalid argument to asctime_s.");
		exit(1);
	}
	//printf("%.19s %s\n", timebuf, am_pm);
	out += timebuf;
	return out;
}

int Sign(double in)
{
	return in >= 0 ? 1 : -1;
}

Globals def;
