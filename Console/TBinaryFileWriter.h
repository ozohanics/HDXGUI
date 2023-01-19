#pragma once

#include <iostream>
#include "common.h"

static const int MAX_BYTES_TO_READ = 1024;
#pragma warning (disable : 4250) //disable message related to inheritance
class TBinaryFileWriter :
	public fstream
{
public:
	TBinaryFileWriter(void);

	~TBinaryFileWriter(void);

	fstream* OpenFile(std::string fileName);
	fstream* OpenTextFile(std::string fileName);
	void	  CloseFile();

	fstream* WriteValue(const double &value);
	fstream* WriteValue(const float &value);
	fstream* WriteValue(const int &value);
	fstream* WriteValue(const bool &value);
	fstream* WriteValue(const long &value);
	fstream* WriteValue(const string &value);
	fstream* WriteValue(const char* value, int size);

	fstream* ReadValue(double* out);
	fstream* ReadValue(float* out);
	fstream* ReadValue(int* out);
	fstream* ReadValue(bool* out);
	fstream* ReadValue(long* out);
	fstream* ReadValue(string *out);
	fstream* ReadValue(char* out, int size);
};
