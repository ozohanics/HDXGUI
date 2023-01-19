// HDXConsole.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include "HDXMain.h"
#include "TDeutDeconv.h"

int main()
{
    //std::cout << "Hello World!\n"; 
    HDXMain HDXApp;

    HDXApp.ReadPeptideList("pepfile.txt");
    HDXApp.ReadRawFileList("rawfile.txt");
    HDXApp.SearchDataFiles();
    //TDeutDeconv tdd(8);
	//double iso[] = { 1, 0.5, 0.25, 0.125, 0.0625};
	//double peaks[] = { 0, 0.5, 0.75, 0.375, 0.1875, 0.09375, 0.03125, 0	};
    //tdd.FillArrays(iso, peaks, 5, 3);
    //tdd.CalcCurrIsotopeSum();
    //tdd.WriteArrayToStdOut();
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
