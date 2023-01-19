#include "TBinaryFileWriter.h"

TBinaryFileWriter::TBinaryFileWriter(void)
{
	//default constructor
}

TBinaryFileWriter::~TBinaryFileWriter(void)
{
	CloseFile();
}

fstream* TBinaryFileWriter::OpenFile(std::string fileName)
{
	if (!FileExists(fileName))
	{//create file
		FILE *fp;
		fopen_s(&fp, fileName.c_str(), "wb");
		fclose(fp);
	}
	open(fileName, ios::binary | ios::out | ios::in);

	if (bad()) WRITE_DEBUG("Could not open/create: " + fileName);

	return this;
}

fstream* TBinaryFileWriter::OpenTextFile(std::string fileName)
{
	if (!FileExists(fileName))
	{//create file
		FILE *fp;
		fopen_s(&fp, fileName.c_str(), "wb");
		fclose(fp);
	}
	open(fileName, ios::out | ios::in);

	if (bad()) WRITE_DEBUG("Could not open/create: " + fileName);

	return this;
}

void	  TBinaryFileWriter::CloseFile()
{
	if (good()) {close();}
}

fstream* TBinaryFileWriter::WriteValue(const double &value) //write T to a file in binary mode
{
	int v_size = sizeof(value);
	WriteValue(reinterpret_cast<const char *>(&value), v_size);
	return this;
}

fstream* TBinaryFileWriter::WriteValue(const float &value)
{
	int v_size = sizeof(value);
	WriteValue(reinterpret_cast<const char *>(&value), v_size);
	return this;
}

fstream* TBinaryFileWriter::WriteValue(const int &value)
{
	int v_size = sizeof(value);
	WriteValue(reinterpret_cast<const char *>(&value), v_size);
	return this;
}

fstream* TBinaryFileWriter::WriteValue(const bool &value)
{
	int v_size = sizeof(value);
	char c[2] = { '\1', '\0' };
	WriteValue(reinterpret_cast<const char *>(&value), v_size);
	return this;
}

fstream* TBinaryFileWriter::WriteValue(const long &value)
{
	int v_size = sizeof(value);
	WriteValue(reinterpret_cast<const char *>(&value), v_size);
	return this;
}

fstream* TBinaryFileWriter::WriteValue(const string &value)
{
	int v_size = value.length()*sizeof(value[0]);
	WriteValue(value.data(), v_size);
	WriteValue(0); //use 0 char to end string data
	return this;
}

fstream* TBinaryFileWriter::WriteValue(const char* value, int size)
{
	write(value,size);
	return this;
}

fstream* TBinaryFileWriter::ReadValue(double* out)
{
	return ReadValue(reinterpret_cast<char *>(out), sizeof(*out));
}

fstream* TBinaryFileWriter::ReadValue(float* out)
{
	return ReadValue(reinterpret_cast<char *>(out), sizeof(*out));
}

fstream* TBinaryFileWriter::ReadValue(int* out)
{
	return ReadValue(reinterpret_cast<char *>(out), sizeof(*out));
}

fstream* TBinaryFileWriter::ReadValue(bool* out)
{
	return ReadValue(reinterpret_cast<char *>(out), sizeof(*out));
}

fstream* TBinaryFileWriter::ReadValue(long* out)
{
	return ReadValue(reinterpret_cast<char *>(out), sizeof(*out));
}

fstream* TBinaryFileWriter::ReadValue(string *out)
{
	char* buffer = new char[MAX_BYTES_TO_READ];
	getline(buffer,MAX_BYTES_TO_READ,'\0');
	seekg(3,ios_base::cur); //move one char sizeof(int)-1
	*out = buffer;
	return this;
}

fstream* TBinaryFileWriter::ReadValue(char* out, int size)
{
	read(out, size);
	return this;
}