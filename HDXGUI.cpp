// Dear ImGui: standalone example application for SDL2 + OpenGL
// (SDL is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan/Metal graphics context creation, etc.)
// If you are new to Dear ImGui, read documentation from the docs/ folder + read the top of imgui.cpp.
// Read online: https://github.com/ocornut/imgui/tree/master/docs

#include "imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"
#include <stdio.h>
#include <string>
#include <SDL.h>
#include <SDL_image.h>
#include "implot.h" 

#include "Console/TPeptideStore.h"
#include "Console/rawclass.h"
#include "Console/HDXMain.h"
#include <iostream>
//#include "imgui_demo.cpp"
#include "ImGuiFileDialog/ImGuiFileDialog.h"

#define __STDC_WANT_LIB_EXT1__ 1

#if defined(IMGUI_IMPL_OPENGL_ES2)
#include <SDL_opengles2.h>
#else
#include <SDL_opengl.h>
#endif
#include <stdlib.h>

#ifdef _MSC_VER
#define sprintf sprintf_s
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

template<typename T>
T FindMax(T* arr, size_t n)
{
	return *std::max_element(arr, arr + n);
}


// Helper to display a little (?) mark which shows a tooltip when hovered.
// In your own code you may want to display an actual icon if you are using a merged icon fonts (see docs/FONTS.md)
static void HelpMarker(const char* desc)
{
	ImGui::TextDisabled("(?)");
	if (ImGui::IsItemHovered())
	{
		ImGui::BeginTooltip();
		ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
		ImGui::TextUnformatted(desc);
		ImGui::PopTextWrapPos();
		ImGui::EndTooltip();
	}
}

enum class MouseButons
{
	R, M, L, RR, LL, MM, None
	//right, middle, left, double right, double left
};


static struct UIMessageData
{
	float xCoord = 0.0;
	float yCoord = 0.0;
	std::string info = "";
	MouseButons btn;
} currentUIMsg; 


void ButtonSelector(const char* label, ImGuiMouseButton* b) {
	ImGui::PushID(label);
	if (ImGui::RadioButton("LMB", *b == ImGuiMouseButton_Left))
		*b = ImGuiMouseButton_Left;
	ImGui::SameLine();
	if (ImGui::RadioButton("RMB", *b == ImGuiMouseButton_Right))
		*b = ImGuiMouseButton_Right;
	ImGui::SameLine();
	if (ImGui::RadioButton("MMB", *b == ImGuiMouseButton_Middle))
		*b = ImGuiMouseButton_Middle;
	ImGui::PopID();
}

void ModSelector(const char* label, ImGuiKeyModFlags* k) {
	ImGui::PushID(label);
	ImGui::CheckboxFlags("Ctrl", (unsigned int*)k, ImGuiKeyModFlags_Ctrl); ImGui::SameLine();
	ImGui::CheckboxFlags("Shift", (unsigned int*)k, ImGuiKeyModFlags_Shift); ImGui::SameLine();
	ImGui::CheckboxFlags("Alt", (unsigned int*)k, ImGuiKeyModFlags_Alt); ImGui::SameLine();
	ImGui::CheckboxFlags("Super", (unsigned int*)k, ImGuiKeyModFlags_Super);
	ImGui::PopID();
}

void InputMapping(const char* label, ImGuiMouseButton* b, ImGuiKeyModFlags* k) {
	ImGui::LabelText("##", "%s", label);
	if (b != NULL) {
		ImGui::SameLine(100);
		ButtonSelector(label, b);
	}
	if (k != NULL) {
		ImGui::SameLine(300);
		ModSelector(label, k);
	}
}


template <typename T>
inline T RandomRange(T min, T max) {
	T scale = rand() / (T)RAND_MAX;
	return min + scale * (max - min);
}


bool ShowSpectrumPlot(std::string title, float* x1, float*y1, int size1, float center1)
{
	bool isMouseInside = false;
	if (ImPlot::BeginPlot(title.c_str())) {
		ImPlot::SetupAxes("mz", "Intensity", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
		ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 3);
		ImPlot::PlotLine("isotopes", x1, y1, size1);
		//plot centers
		float centX[4], centY[4];
		float vertPanningVal = FindMax(y1, size1) * 1.1;
		centX[0] = center1 - 0.025; //possible to set error bars for center accuracy
		centX[1] = center1 - 0.025;
		centX[2] = center1 + 0.025;
		centX[3] = center1 + 0.025;
		centY[0] = 0.0;
		centY[1] = vertPanningVal; //vertical panning
		centY[2] = vertPanningVal;
		centY[3] = 0.0;

		ImPlot::PlotShaded("center", centX, centY, 4);

		auto mousePos = ImPlot::GetPlotMousePos(IMPLOT_AUTO, IMPLOT_AUTO);		
		currentUIMsg.xCoord = mousePos.x;
		currentUIMsg.yCoord = mousePos.y;

		float sxMax = ImPlot::GetPlotLimits().X.Max;
		float sxMin = ImPlot::GetPlotLimits().X.Min;
		float syMax = ImPlot::GetPlotLimits().Y.Max;
		float syMin = ImPlot::GetPlotLimits().Y.Min;
		//check if mouse is inside
		if (currentUIMsg.xCoord < sxMax && currentUIMsg.xCoord > sxMin && currentUIMsg.yCoord< syMax && currentUIMsg.yCoord > syMin) isMouseInside = true;

		ImPlot::EndPlot();
	}
	return isMouseInside;
}

//x1,y1 - isotopes, x2, y2 - reference, center is loaded separately
bool ShowSpectrumPlotWithRef(std::string title, float* x1, float* y1, float* x2, float* y2, int size1, int size2, float center1)
{
	bool isMouseInside = false;
	if (ImPlot::BeginPlot(title.c_str(), ImVec2(-1, -1))) {
		ImPlot::SetNextLineStyle(ImVec4(0, 0, 0, -1), 2);
		ImPlot::SetupAxes("mz", "Intensity", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
		ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_Outside);
		ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 3);
		
		//ImPlot::SetupFinish();
		//plot current isotopes
		if (size1 > 0) ImPlot::PlotLine("isotopes", x1, y1, size1);
		//plot reference
		if (size2 > 0) ImPlot::PlotLine("reference", x2, y2, size2);		
		//plot centers
		float centX[4], centY[4];
		float vertPanningVal = FindMax(y1, size1) * 1.1;
		centX[0] = center1 - 0.025; //possible to set error bars for center accuracy
		centX[1] = center1 - 0.025;
		centX[2] = center1 + 0.025;
		centX[3] = center1 + 0.025;
		centY[0] = 0.0;
		centY[1] = vertPanningVal; //vertical panning
		centY[2] = vertPanningVal;
		centY[3] = 0.0;

		ImPlot::PlotShaded("center", centX, centY, 4);

		auto mousePos = ImPlot::GetPlotMousePos(IMPLOT_AUTO, IMPLOT_AUTO);
		currentUIMsg.xCoord = mousePos.x;
		currentUIMsg.yCoord = mousePos.y;

		float sxMax = ImPlot::GetPlotLimits().X.Max;
		float sxMin = ImPlot::GetPlotLimits().X.Min;
		float syMax = ImPlot::GetPlotLimits().Y.Max;
		float syMin = ImPlot::GetPlotLimits().Y.Min;
		//check if mouse is inside
		if (currentUIMsg.xCoord < sxMax && currentUIMsg.xCoord > sxMin && currentUIMsg.yCoord< syMax && currentUIMsg.yCoord > syMin) isMouseInside = true;
		
		ImPlot::EndPlot();
	}
	return isMouseInside;
}


//x1,y1 - isotopes, x2, y2 - reference, center is loaded separately
bool ShowChromatogram(HDXMain &HDXApp, int item_current_idx, int file_list_idx)
{
	bool myMouse = false;
	int pNr = HDXApp.peptideStore->pepNr;
	float icX[4], icY[4];
	icX[0] = HDXApp.peptideStore->peptides[item_current_idx].startRT;
	icX[1] = HDXApp.peptideStore->peptides[item_current_idx].startRT;
	icX[2] = HDXApp.peptideStore->peptides[item_current_idx].endRT;
	icX[3] = HDXApp.peptideStore->peptides[item_current_idx].endRT;
	icY[0] = 0;
	icY[3] = 0;
	icY[1] = FindMax(HDXApp.peptideStore->xEIC[item_current_idx + pNr * file_list_idx].y.data(), HDXApp.peptideStore->xEIC[item_current_idx + pNr * file_list_idx].y.size());
	icY[2] = icY[1];

	//ImGui::Begin("Chromatogram");
	if (ImPlot::BeginPlot((HDXApp.peptideStore->peptides[item_current_idx].sequence + "_" + FormatDouble(HDXApp.peptideStore->xEIC[item_current_idx + pNr * file_list_idx].baseMass, 2)).c_str(), ImVec2(-1, -1))) {
		ImPlot::SetupAxes("RT", "Intensity", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
		//ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 3);
		ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_Horizontal);
		ImPlot::PlotShaded("SeachWindow", icX, icY, 4);
		ImPlot::PlotLine("EIC", HDXApp.peptideStore->xEIC[item_current_idx + pNr * file_list_idx].rt.data(), HDXApp.peptideStore->xEIC[item_current_idx + pNr * file_list_idx].y.data(), HDXApp.peptideStore->xEIC[item_current_idx + pNr * file_list_idx].rt.size());
		auto mousePos = ImPlot::GetPlotMousePos(IMPLOT_AUTO, IMPLOT_AUTO);
		currentUIMsg.xCoord = mousePos.x;
		currentUIMsg.yCoord = mousePos.y;
		float sxMax = ImPlot::GetPlotLimits().X.Max;
		float sxMin = ImPlot::GetPlotLimits().X.Min;
		float syMax = ImPlot::GetPlotLimits().Y.Max;
		float syMin = ImPlot::GetPlotLimits().Y.Min;
		//check if mouse is inside
		if (currentUIMsg.xCoord < sxMax && currentUIMsg.xCoord > sxMin && currentUIMsg.yCoord< syMax && currentUIMsg.yCoord > syMin) myMouse = true;

		ImPlot::EndPlot();
	}
	return myMouse;
}


//x1,y1 - isotopes, x2, y2 - reference, center is loaded separately
bool ShowRAWSpectrum(std::string title, float* x1, float* y1, int size1)
{
	bool isMouseInside = false;
	if (ImPlot::BeginPlot(title.c_str(), ImVec2(-1, -1))) {
		ImPlot::SetNextLineStyle(ImVec4(0, 0, 0, -1), 2);
		ImPlot::SetupAxes("mz", "Intensity", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
		ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);
		ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 3);

		//ImPlot::SetupFinish();
		//plot current isotopes
		if (size1 > 0) ImPlot::PlotLine(title.c_str(), x1, y1, size1);
		
		ImPlot::EndPlot();
	}
	return isMouseInside;
}

void ShowMouseSetupPanel(bool* showMe)
{
	ImGui::Begin("Mouse Panel", showMe);
	ImPlotInputMap& map = ImPlot::GetInputMap();
	InputMapping("Pan", &map.Pan, &map.PanMod);
	InputMapping("Fit", &map.Fit, NULL);
	InputMapping("Select", &map.Select, &map.SelectMod);
	InputMapping("SelectHorzMod", NULL, &map.SelectHorzMod);
	InputMapping("SelectVertMod", NULL, &map.SelectVertMod);
	InputMapping("SelectCancel", &map.SelectCancel, NULL);
	InputMapping("Menu", &map.Menu, NULL);
	InputMapping("OverrideMod", NULL, &map.OverrideMod);
	InputMapping("ZoomMod", NULL, &map.ZoomMod);
	ImGui::SliderFloat("ZoomRate", &map.ZoomRate, -1, 1);
	ImGui::End();
}

bool DisplayFileDialog(std::string& filePathName, std::string dialogKey)
{
	// display and action if ok
	if (ImGuiFileDialog::Instance()->Display(dialogKey))
	{
		if (ImGuiFileDialog::Instance()->IsOk())
		{
			filePathName = ImGuiFileDialog::Instance()->GetFilePathName();
			//filePath = ImGuiFileDialog::Instance()->GetCurrentPath();
			//std::string filter = ImGuiFileDialog::Instance()->GetCurrentFilter();
			// here convert from string because a string was passed as a userDatas, but it can be what you want
			std::string userDatas;
			if (ImGuiFileDialog::Instance()->GetUserDatas())
				userDatas = std::string((const char*)ImGuiFileDialog::Instance()->GetUserDatas());
			auto selection = ImGuiFileDialog::Instance()->GetSelection(); // multiselection
			// action
		}
		// close
		ImGuiFileDialog::Instance()->Close();
		return true;
	}

	return false;
}


static HDXMain HDXApp;
static bool is_HDX_calc_done = false;

bool StartHDXCalc(std::string pepPath = "pepfile.txt", std::string rawPath = "rawfile.txt", bool useOverlap = true)
{
	HDXApp.ReInit();
	HDXApp.peptideStore->useXICOverlap = useOverlap;

	HDXApp.rawFileListPath = rawPath;
	HDXApp.peptideListPath = pepPath;
	
	HDXApp.ReadPeptideList(HDXApp.peptideListPath);
	HDXApp.ReadRawFileList(HDXApp.rawFileListPath);
	HDXApp.SearchDataFiles();
	is_HDX_calc_done = true;
	return true;
}

bool ModellAA()
{
	int pepmin = 0;
	int pepmax = 0;
	for (int i = 0; i < HDXApp.peptideNr; i++)
	{
		if (HDXApp.peptideStore->peptides[i].seqEnd > pepmax) pepmax = HDXApp.peptideStore->peptides[i].seqEnd;
	}

	struct AAData
	{
		double D; //deut incorporation
		double Dmax; //total D for peptide
		int aaStart; //smallest peptide data
		int aaEnd; //smallest peptide data
		int aaLoc; //location in protein
	} tmpAA;

	vector<double> proteinD(pepmax + 1, -1.0);
	vector<string> protein(pepmax + 1, "X");
	vector<AAData> proteinAA(pepmax + 1, tmpAA);

#define MIN(a,b) (((a)<=(b)) ? (a) : (b))
#define MAX(a,b) (((a)>=(b)) ? (a) : (b))
	for (int f = 0; f < HDXApp.fileNr; f++)
	{
		//if (HDXApp.rawData->isRef) continue;
		if (f != 1) continue;

		for (int i = 0; i < HDXApp.peptideNr; i++)
		{
			double aaDeutPercent = HDXApp.peptideStore->results[i].deutResults[f].deutDelta / HDXApp.peptideStore->peptides[i].maxD;
			for (int d = HDXApp.peptideStore->peptides[i].seqStart; d < HDXApp.peptideStore->peptides[i].seqEnd; d++)
			{
				if (proteinD[d - 1] < 0) //not set
				{
					protein[d - 1] = HDXApp.peptideStore->peptides[i].sequence.substr(d - HDXApp.peptideStore->peptides[i].seqStart, 1);
					if (string(protein[d - 1]) != "P") proteinD[d - 1] = aaDeutPercent;
					proteinAA[d - 1].aaStart = HDXApp.peptideStore->peptides[i].seqStart;
					proteinAA[d - 1].aaEnd = HDXApp.peptideStore->peptides[i].seqEnd;
					proteinAA[d - 1].aaLoc = d - 1;
					proteinAA[d - 1].D = proteinD[d - 1];
					proteinAA[d - 1].Dmax = HDXApp.peptideStore->results[i].deutResults[f].deutDelta;
				}
				else
				{
					//if set 
					if (string(protein[d - 1]) != "P")
					{
						double maxD = HDXApp.peptideStore->results[i].deutResults[f].deutDelta;
						//maximum D for peptide
						if (proteinAA[d - 1].aaStart <= HDXApp.peptideStore->peptides[i].seqStart && proteinAA[d - 1].aaEnd >= HDXApp.peptideStore->peptides[i].seqEnd)
							proteinAA[d - 1].Dmax = maxD;
						proteinAA[d - 1].aaStart = MAX(HDXApp.peptideStore->peptides[i].seqStart, proteinAA[d - 1].aaStart);
						proteinAA[d - 1].aaEnd = MIN(HDXApp.peptideStore->peptides[i].seqEnd, proteinAA[d - 1].aaEnd);
						proteinAA[d - 1].D = MIN(aaDeutPercent, proteinAA[d - 1].D);

					}
				}
			}
		}

		for (int p = 0; p < pepmax; p++)
		{
			int s = proteinAA[p].aaStart;
			int e = proteinAA[p].aaEnd;
			int seqStart = p;
			int seqEnd = p;
			bool seqFound = false;
			for (int x = p + 1; x < pepmax; x++)
			{
				if (seqFound) break;
				if (s < proteinAA[x].aaStart)
				{
					seqStart = x;
					for (size_t i = 0; i < pepmax; i++)
					{
						if (e == proteinAA[i].aaEnd)
						{
							seqEnd = i;
						}
						else
						{
							seqFound = true;
							break;
						}
					}

				}
			}
			if (seqFound)
			{
				double newD_before = MAX(0.0, proteinAA[seqStart].D - proteinAA[p].D);
				double newD_after = proteinAA[seqStart].D - newD_before;
				for (size_t b = p; b < seqStart; b++)
				{
					proteinAA[b].D = newD_before / (seqStart - p);
				}
				for (size_t b = p; b < seqStart; b++)
				{

				}
			}
		}
	}
	return false;
}

static float xs1[1001], ys1[1001]; //spectra
static float xs2[1001], ys2[1001]; //reference
static float xs3[1001], ys3[1001]; //raw data

int FillSpectrumArray(int item_current_idx, int file_list_idx)
{
	int arrCounter = 0;
	for (int c = 0; c < HDXApp.peptideStore->results[item_current_idx].ICs[file_list_idx].isotopeNr; c++)
	{
		if (c == 0)
		{
			xs1[arrCounter] = HDXApp.peptideStore->results[item_current_idx].ICs[file_list_idx].MZ[c] - 0.5;
			ys1[arrCounter] = 0;
			arrCounter++;
		}
		xs1[arrCounter] = HDXApp.peptideStore->results[item_current_idx].ICs[file_list_idx].MZ[c];
		ys1[arrCounter] = 0;
		arrCounter++;
		xs1[arrCounter] = HDXApp.peptideStore->results[item_current_idx].ICs[file_list_idx].MZ[c];
		ys1[arrCounter] = HDXApp.peptideStore->results[item_current_idx].ICs[file_list_idx].intensity[c];
		arrCounter++;
		xs1[arrCounter] = HDXApp.peptideStore->results[item_current_idx].ICs[file_list_idx].MZ[c];
		ys1[arrCounter] = 0;
		arrCounter++;
		if (c == HDXApp.peptideStore->results[item_current_idx].ICs[file_list_idx].isotopeNr - 1)
		{
			xs1[arrCounter] = HDXApp.peptideStore->results[item_current_idx].ICs[file_list_idx].MZ[c] + 0.5;
			ys1[arrCounter] = 0;
			arrCounter++;
		}
	}
	return arrCounter;
}

int FillRefArray(int item_current_idx, int refIdx)
{
	int refCounter = 0;
	for (int c = 0; c < HDXApp.peptideStore->results[item_current_idx].ICs[refIdx].isotopeNr; c++)
	{
		if (c == 0)
		{
			xs2[refCounter] = HDXApp.peptideStore->results[item_current_idx].ICs[refIdx].MZ[c] - 0.5;
			ys2[refCounter] = 0;
			refCounter++;
		}
		xs2[refCounter] = HDXApp.peptideStore->results[item_current_idx].ICs[refIdx].MZ[c];
		ys2[refCounter] = 0;
		refCounter++;
		xs2[refCounter] = HDXApp.peptideStore->results[item_current_idx].ICs[refIdx].MZ[c];
		ys2[refCounter] = HDXApp.peptideStore->results[item_current_idx].ICs[refIdx].intensity[c];
		refCounter++;
		xs2[refCounter] = HDXApp.peptideStore->results[item_current_idx].ICs[refIdx].MZ[c];
		ys2[refCounter] = 0;
		refCounter++;
		if (c == HDXApp.peptideStore->results[item_current_idx].ICs[refIdx].isotopeNr - 1)
		{
			xs2[refCounter] = HDXApp.peptideStore->results[item_current_idx].ICs[refIdx].MZ[c] + 0.5;
			ys2[refCounter] = 0;
			refCounter++;
		}
	}
	return refCounter;
}

int FillIsoRefArray(int item_current_idx, int refIdx , int ys1Count)
{

	int pointNr = HDXApp.peptideStore->peptides[item_current_idx].composition.isoDataNr; //see all isotopes
	//pointNr = HDXApp.peptideStore->results[item_current_idx].ICs[refIdx].isotopeNr; //see only points in the data

	int refCounter = 0;
	double maxI = FindMax(ys1, ys1Count);

	for (int c = 0; c < pointNr; c++)
	{
		if (c == 0)
		{
			xs2[refCounter] = HDXApp.peptideStore->peptides[item_current_idx].composition.isoMass[c] / HDXApp.peptideStore->peptides[item_current_idx].z + masses::Hp - 0.5;
			ys2[refCounter] = 0;
			refCounter++;
		}
		xs2[refCounter] = HDXApp.peptideStore->peptides[item_current_idx].composition.isoMass[c] / HDXApp.peptideStore->peptides[item_current_idx].z + masses::Hp;
		ys2[refCounter] = 0;
		refCounter++;
		xs2[refCounter] = HDXApp.peptideStore->peptides[item_current_idx].composition.isoMass[c] / HDXApp.peptideStore->peptides[item_current_idx].z + masses::Hp;
		ys2[refCounter] = HDXApp.peptideStore->peptides[item_current_idx].composition.isoData[c] * maxI / HDXApp.peptideStore->peptides[item_current_idx].composition.isoData[HDXApp.peptideStore->peptides[item_current_idx].composition.maxIdx];
		refCounter++;
		xs2[refCounter] = HDXApp.peptideStore->peptides[item_current_idx].composition.isoMass[c] / HDXApp.peptideStore->peptides[item_current_idx].z + masses::Hp;
		ys2[refCounter] = 0;
		refCounter++;
		if (c == pointNr - 1)
		{
			xs2[refCounter] = HDXApp.peptideStore->peptides[item_current_idx].composition.isoMass[c] / HDXApp.peptideStore->peptides[item_current_idx].z + masses::Hp + 0.5;
			ys2[refCounter] = 0;
			refCounter++;
		}
	}
	return refCounter;
}


int FillRAWSpectrumArray(int item_current_idx, int file_list_idx, double window)
{
	int arrCounter = 0;
	
	if (HDXApp.last_file_idx != file_list_idx)
	{
		HDXApp.last_file_idx = file_list_idx;
		HDXApp.rawData->XMLToMemory(HDXApp.peptideStore->dataInfo[file_list_idx].rawPath);
	}
	double curr_mass = HDXApp.peptideStore->peptides[item_current_idx].mass / HDXApp.peptideStore->peptides[item_current_idx].z + masses::Hp;
	long   curr_scan = window; //HDXApp.peptideStore->peptides[item_current_idx].maxIntLoc + 
	//find index of mass in current scan, used for extracting a part of the spectrum 
	long current_scan__low_idx = HDXApp.rawData->BinarySearchClosestSpectrumIndex(curr_scan, curr_mass, 0.25) - 5; 
	long current_scan__high_idx = current_scan__low_idx + 17;
		//HDXApp.rawData->BinarySearchClosestSpectrumIndex(curr_scan, curr_mass + 30, 0.25);

	for (int c = current_scan__low_idx; c < current_scan__high_idx; c++)
	{
		if (c == current_scan__low_idx)
		{
			xs3[arrCounter] = HDXApp.rawData->GetFileDataMZArray(curr_scan)[c] - 0.5;
			ys3[arrCounter] = 0;
			arrCounter++;
		}
		xs3[arrCounter] = HDXApp.rawData->GetFileDataMZArray(curr_scan)[c];
		ys3[arrCounter] = 0;
		arrCounter++;
		xs3[arrCounter] = HDXApp.rawData->GetFileDataMZArray(curr_scan)[c];
		ys3[arrCounter] = HDXApp.rawData->GetFileDataYArray(curr_scan)[c];
		arrCounter++;
		xs3[arrCounter] = HDXApp.rawData->GetFileDataMZArray(curr_scan)[c];
		ys3[arrCounter] = 0;
		arrCounter++;
		if (c == current_scan__high_idx - 1)
		{
			xs3[arrCounter] = HDXApp.rawData->GetFileDataMZArray(curr_scan)[c] + 0.5;
			ys3[arrCounter] = 0;
			arrCounter++;
		}
	}
	return arrCounter;
}


void SaveSpectralDataToFile()
{
	TStringList fileOut;
	for (size_t i = 0; i < 100; i++)
	{
		if (ys1[i] > 0.5)//
		{
			fileOut.Add(FormatDouble(xs1[i], 4) + "\t" + FormatDouble(ys1[i], 4));
		}
		
	}
	
	fileOut.SaveToFile("SpectrumX1.txt");
	fileOut.clear();
	//fileOut.Add();
	for (size_t i = 0; i < 100; i++)
	{
		if (ys2[i] > 0.5)
		{
			fileOut.Add(FormatDouble(xs2[i], 4) + "\t" + FormatDouble(ys2[i], 4));
		}
	}
	fileOut.SaveToFile("SpectrumX2.txt");
}


bool saveScreenshotBMP(std::string filepath, SDL_Window* SDLWindow, SDL_Renderer* SDLRenderer) {
	SDL_Surface* saveSurface = NULL;
	SDL_Surface* infoSurface = NULL;
	infoSurface = SDL_GetWindowSurface(SDLWindow);
	if (infoSurface == NULL) {
		std::cerr << "Failed to create info surface from window in saveScreenshotBMP(string), SDL_GetError() - " << SDL_GetError() << "\n";
	}
	else {
		unsigned char* pixels = new (std::nothrow) unsigned char[infoSurface->w * infoSurface->h * infoSurface->format->BytesPerPixel];
		if (pixels == 0) {
			std::cerr << "Unable to allocate memory for screenshot pixel data buffer!\n";
			return false;
		}
		else {
			if (SDL_RenderReadPixels(SDLRenderer, &infoSurface->clip_rect, infoSurface->format->format, pixels, infoSurface->w * infoSurface->format->BytesPerPixel) != 0) {
				std::cerr << "Failed to read pixel data from SDL_Renderer object. SDL_GetError() - " << SDL_GetError() << "\n";
				delete[] pixels;
				return false;
			}
			else {
				saveSurface = SDL_CreateRGBSurfaceFrom(pixels, infoSurface->w, infoSurface->h, infoSurface->format->BitsPerPixel, infoSurface->w * infoSurface->format->BytesPerPixel, infoSurface->format->Rmask, infoSurface->format->Gmask, infoSurface->format->Bmask, infoSurface->format->Amask);
				if (saveSurface == NULL) {
					std::cerr << "Couldn't create SDL_Surface from renderer pixel data. SDL_GetError() - " << SDL_GetError() << "\n";
					delete[] pixels;
					return false;
				}
				SDL_SaveBMP(saveSurface, filepath.c_str());
				SDL_FreeSurface(saveSurface);
				saveSurface = NULL;
			}
			delete[] pixels;
		}
		SDL_FreeSurface(infoSurface);
		infoSurface = NULL;
	}
	return true;
}

void Screenshot(int x, int y, int w, int h, const char* filename)
{
	unsigned char* pixels = new unsigned char[w * h * 4]; // 4 bytes for RGBA
	glReadPixels(x, y, w, h, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

	SDL_Surface* surf = SDL_CreateRGBSurfaceFrom(pixels, w, h, 8 * 4, w * 4, 0, 0, 0, 0);
	SDL_SaveBMP(surf, filename);

	SDL_FreeSurface(surf);
	delete[] pixels;
}

// Main code
int main(int, char**)
{
	
	// Setup SDL
	// (Some versions of SDL before <2.0.10 appears to have performance/stalling issues on a minority of Windows systems,
	// depending on whether SDL_INIT_GAMECONTROLLER is enabled or disabled.. updating to latest version of SDL is recommended!)
	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) != 0)
	{
		printf("Error: %s\n", SDL_GetError());
		return -1;
	}

	// GL 3.0 + GLSL 130
	const char* glsl_version = "#version 130";
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, 0);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);

	// Create window with graphics context
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
	SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
	SDL_WindowFlags window_flags = (SDL_WindowFlags)(SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI );
	SDL_Window* window = SDL_CreateWindow("HDXGUI SDL2+OpenGL3", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1400, 1000, window_flags);
	SDL_GLContext gl_context = SDL_GL_CreateContext(window);
	SDL_GL_MakeCurrent(window, gl_context);
	SDL_GL_SetSwapInterval(1); // Enable vsync
	//SDL_Surface
	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImPlot::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();

	// Setup Platform/Renderer backends
	ImGui_ImplSDL2_InitForOpenGL(window, gl_context);
	ImGui_ImplOpenGL3_Init(glsl_version);
	SDL_Surface* sshot = SDL_GetWindowSurface(window);

	// Our state
	bool show_mouse_window = false;
	bool show_another_window = false;
	bool useOverlap = true; //used for isotope extracted chromatogram overlap check

	ImVec4 clear_color = ImVec4(0.145f, 0.55f, 0.60f, 0.70f);
	
	//datafiles for searching, needed in two forms, because of Input
	char rawPath[255] = "rawfile.txt";
	char pepPath[255] = "pepfile.txt";
	std::string rawStr = "rawfile.txt";
	std::string pepStr = "pepfile.txt";
	
	// Main loop
	static bool done = false;
	std::string resTitle = "Results ";
	static int item_current_idx = 0; // Here we store our selection data as an index.
	static int file_list_idx = 0; // Here we store our selection data as an index.
	static double refError = 0.25;
	static double deutProtLimit = -0.05;
	static double retTimeWindow = 0.25;
	static double rawScanWindow = 1;
	static bool useIsoMax = true; //use isotope data instead of scan basepeak
	int rawCounter = ERROR_VALUE;

	while (!done)
	{
		static float f = 0.0f;
		static int counter = 0;
		
		/*for (int i = 1; i < 100; i++)
		{//check which key is which
		//delete is 76
		//enter is 40
			auto keys = io.KeysDown[i];
			if (keys) resTitle = "Results " + to_string(i);
		}*/
		if (is_HDX_calc_done)
		{
			auto keys = io.KeysDown[76];
			if (keys) // delete the selected peptide
			{
				//delete not working correctly yet
				HDXApp.DeletePeptide(item_current_idx);
				//item_current_idx = 0;
			}
		}
		Sleep(25);
		// Poll and handle events (inputs, window resize, etc.)
		// You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
		// - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application.
		// - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application.
		// Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
		SDL_Event event;
		while (SDL_PollEvent(&event))
		{
			ImGui_ImplSDL2_ProcessEvent(&event);
			if (event.type == SDL_QUIT)
				done = true;
			if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_CLOSE && event.window.windowID == SDL_GetWindowID(window))
				done = true;
		}

		// Start the Dear ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplSDL2_NewFrame();
		ImGui::NewFrame();
		{//frame begin

			int nr = HDXApp.peptideStore->pepNr;
			std::string defVal = "";
			vector<std::string> peptideNames(2,defVal);

			if (nr > 0) peptideNames.clear();
			for (int i = 0; i < nr; i++)
			{
				peptideNames.push_back(to_string(i)+"_" + HDXApp.peptideStore->peptides[i].sequence + "_" + to_string(HDXApp.peptideStore->peptides[i].z) + "+");
			}


			//Control Panel
			ImGui::Begin("Control Panel");
			ImGui::Checkbox("Show Parameter Window", &show_another_window); ImGui::SameLine();
			ImGui::Checkbox("Show Mouse Settings", &show_mouse_window); ImGui::SameLine();
			//ImGui::Checkbox("Use isotopes for XIC", &useIsoMax); ImGui::SameLine();
			ImGui::Checkbox("Use isotopes overlap", &useOverlap); ImGui::SameLine();
			HDXApp.peptideStore->useXICOverlap = useOverlap;

			ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(7 / 7.0f, 0.6f, 0.6f));
			ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(1 / 7.0f, 0.7f, 0.7f));
			ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(5 / 7.0f, 0.8f, 0.8f));
			
			ImGui::NewLine();
			if (ImGui::Button("Start Calculation")) {
				is_HDX_calc_done = StartHDXCalc(pepStr, rawStr, useOverlap);
			}
			ImGui::SameLine(); 
			if (ImGui::Button("Save Spectrum")) { SaveSpectralDataToFile(); }
			ImGui::SameLine();
			if (ImGui::Button("Save Results")) 
			{
				if (is_HDX_calc_done)
				{
					HDXApp.SaveResults();
					HDXApp.SaveModifiedInput();
				}
			}

			if (is_HDX_calc_done)
			{//show setings for peptide charge
				double eRT = HDXApp.peptideStore->peptides[item_current_idx].endRT;
				double sRT = HDXApp.peptideStore->peptides[item_current_idx].startRT;
				int myCharge = HDXApp.peptideStore->peptides[item_current_idx].z;
				if (ImGui::Button("Increase charge"))
				{
					HDXApp.ReEvaluatePeptide(item_current_idx, file_list_idx, sRT, eRT, useOverlap, myCharge + 1);
				}
				ImGui::SameLine();
				if (ImGui::Button("Decrease charge"))
				{
					HDXApp.ReEvaluatePeptide(item_current_idx, file_list_idx, sRT, eRT, useOverlap, myCharge - 1);
				}
				ImGui::SameLine();
				if (ImGui::Button("ReCalc peptide"))
				{
					HDXApp.ReEvaluatePeptide(item_current_idx, file_list_idx, sRT, eRT, useOverlap, myCharge);
				}


			}
			ImGui::NewLine();
			ImGui::PushItemWidth(450);

			const char* combo_preview_value = peptideNames[item_current_idx].c_str();  // Pass in the preview value visible before opening the combo (it could be anything)
			if (ImGui::BeginCombo("Peptide List", combo_preview_value))//,ImVec2(0, 15 * ImGui::GetTextLineHeightWithSpacing())
			{
				for (int n = 0; n < nr; n++)
				{
					const bool is_selected = (item_current_idx == n);
					if (ImGui::Selectable(peptideNames[n].c_str(), is_selected))
						item_current_idx = n;
					// Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
					if (is_selected)
					{
						ImGui::SetItemDefaultFocus();
					}
				}
				ImGui::EndCombo();
			}
			/*ImGui::SameLine();
			ImGui::BeginChildFrame(1, ImVec2(150, 150));
			ImGui::Button("TestBtn");
			ImGui::EndChildFrame();*/
			ImGui::NewLine();

			if (ImGui::BeginListBox("Raw File List",ImVec2(0,75)))
			{
				for (int n = 0; n < HDXApp.fileNr; n++)
				{
					const bool is_selected = (file_list_idx == n);
					if (ImGui::Selectable(HDXApp.fileList[n].c_str(), is_selected))
						file_list_idx = n;
					// Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
					if (is_selected)
					{
						ImGui::SetItemDefaultFocus();
					}
				}
				ImGui::EndListBox();
			}
			ImGui::PopItemWidth();
			ImGui::NewLine();
			ImGui::End(); //End Control panel

			//Plot data
			if (is_HDX_calc_done)
			{
				
				int refIdx = HDXApp.peptideStore->refFileIdx;
				//current file data
				int arrCounter = FillSpectrumArray(item_current_idx, file_list_idx);
				//reference file data
				//int refCounter = FillRefArray(item_current_idx, refIdx);
				int refCounter = FillIsoRefArray(item_current_idx, refIdx, arrCounter);

				double center = ((HDXApp.peptideStore->results[item_current_idx].deutResults[file_list_idx].massCenter) / HDXApp.peptideStore->peptides[item_current_idx].z + masses::Hp);
				double refcenter = ((HDXApp.peptideStore->results[item_current_idx].isoRefData.massCenter) / HDXApp.peptideStore->peptides[item_current_idx].z + masses::Hp);
				//store a raw spectrum
				rawCounter = FillRAWSpectrumArray(item_current_idx, file_list_idx, rawScanWindow);

				//Plot
				ImGui::Begin("Spectra Panel");
				ImGui::BulletText("See the plot's context menu for setting. \n Spectrum: Double left click for XIC, 	double middle to remove isotope. \n Double middle click inside chromatogram to set window");

				bool isMouseIn = ShowSpectrumPlotWithRef(("Spectrum: " + HDXApp.peptideStore->peptides[item_current_idx].sequence + "_" + to_string(HDXApp.peptideStore->peptides[item_current_idx].z) + "_deltaD_" + FormatDouble(HDXApp.peptideStore->results[item_current_idx].deutResults[file_list_idx].deutDelta, 2)).c_str(), xs1, ys1, xs2, ys2, arrCounter, refCounter, center);
				//spectrum mouse events
				if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Middle) && isMouseIn)
				{
					//currentUIMsg.btn = MouseButons::MM;
					HDXApp.ResetIsotope(item_current_idx, file_list_idx, currentUIMsg.xCoord - 0.07, currentUIMsg.xCoord + 0.07, useIsoMax);
				}
				if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left) && isMouseIn)
				{

					if (HDXApp.last_file_idx != file_list_idx)
					{
						HDXApp.last_file_idx = file_list_idx;
						HDXApp.rawData->XMLToMemory(HDXApp.peptideStore->dataInfo[file_list_idx].rawPath);
					}
					double scanMaxMass = currentUIMsg.xCoord;
					double sDelta = scanMaxMass* DEUTMASSERRORPPM /1e6; //scan with high inacuracy
					HDXApp.peptideStore->CreateXIC(item_current_idx, file_list_idx, scanMaxMass, sDelta, HDXApp.rawData);
				}

				ImGui::End();
				
				ImGui::Begin("Chroma Panel");
				bool chromMouse = false;
				chromMouse = ShowChromatogram(HDXApp,item_current_idx, file_list_idx);
				//chromatogram mouse events
				if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Middle) && chromMouse)
				{
					//currentUIMsg.btn = MouseButons::MM;
					HDXApp.ReEvaluatePeptide(item_current_idx, file_list_idx, currentUIMsg.xCoord - retTimeWindow, currentUIMsg.xCoord + retTimeWindow, useOverlap);
				}
				
				ImGui::End(); //Plot end


				ImGui::Begin(resTitle.c_str());
				static ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable;
				//const float TEXT_BASE_WIDTH = ImGui::CalcTextSize("A").x;
				

				if (ImGui::BeginTable("ResultTable", 5, flags))
				{
					
					// Display headers so we can inspect their interaction with borders.
					// (Headers are not the main purpose of this section of the demo, so we are not elaborating on them too much. See other sections for details)
					if (true)
					{
						ImGui::TableSetupColumn("Nr.");
						ImGui::TableSetupColumn("Start-End");
						ImGui::TableSetupColumn("Sequence");
						ImGui::TableSetupColumn("deltaD");
						ImGui::TableSetupColumn("percentD");
						ImGui::TableHeadersRow();
					}

					bool isColored = false;
					for (int row = 0; row < HDXApp.peptideStore->pepNr; row++)
					{//peptides in a row
																	
						ImGui::TableNextRow();
						if (isColored)
						{
							ImGui::PopStyleColor(2);
						}
						
						const bool is_selected = (item_current_idx == row);
						for (int column = 0; column < 5; column++)
						{
							ImGui::TableSetColumnIndex(column);
							switch (column)
							{
							case 0:
								ImGui::TextUnformatted(to_string(row).c_str());
								break;
							case 1:
								ImGui::TextUnformatted(HDXApp.peptideStore->peptides[row].proteinID.c_str());
								break;
							case 2:							
								if (ImGui::Selectable(HDXApp.peptideStore->peptides[row].sequence.c_str(), is_selected)) {
									item_current_idx = row;
									rawScanWindow = HDXApp.peptideStore->results[row].ICs[file_list_idx].maxScanLoc[0];
								}
								break;
							case 3:
									ImGui::TextUnformatted(FormatDouble(HDXApp.peptideStore->results[row].deutResults[file_list_idx].deutDelta, 2).c_str());								
								break;
							case 4:
								ImGui::TextUnformatted(FormatDouble(HDXApp.peptideStore->results[row].deutResults[file_list_idx].deutDelta/ HDXApp.peptideStore->peptides[row].maxD, 2).c_str());
							default:
								break;
							}
							

						}
						//ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(1.0f, 0.0f, 0.0f));
						//ImGui::PushStyleColor(ImGuiCol_TextSelectedBg, (ImVec4)ImColor::HSV(1.0f, 0.0f, 0.0f));
						//ImGui::PushStyleColor(ImGuiCol_NavHighlight, (ImVec4)ImColor::HSV(1.0f, 0.0f, 0.0f));
						//ImGui::PushStyleColor(ImGuiCol_FrameBgHovered, (ImVec4)ImColor::HSV(1.0f, 0.0f, 0.0f));
						//check peptides where the measured data does not fit with the theory
						if (fabsf(HDXApp.peptideStore->results[row].deutResults[file_list_idx].deutDelta) > refError && HDXApp.peptideStore->dataInfo[file_list_idx].labelTime < 0.05) //check for ref file
						{
							ImGui::PushStyleColor(ImGuiCol_TableRowBg, (ImVec4)ImColor::HSV(1.0/7.0f, 0.1f, 0.7f));
							ImGui::PushStyleColor(ImGuiCol_TableRowBgAlt, (ImVec4)ImColor::HSV(2.0 / 7.0f, 0.2f, 0.8f));
							//ImGui::Text(FormatDouble(HDXApp.peptideStore->results[row].deutResults[file_list_idx].deutDelta, 2).c_str());
							isColored = true;
						}
						//check peptides that are protected
						else if ((HDXApp.peptideStore->results[row].deutResults[file_list_idx].deutDelta) < deutProtLimit && HDXApp.peptideStore->dataInfo[file_list_idx].labelTime > 0.05)//check for ref file
						{
							ImGui::PushStyleColor(ImGuiCol_TableRowBg, (ImVec4)ImColor::HSV(6.0 / 7.0f, 0.5f, 0.5f));
							ImGui::PushStyleColor(ImGuiCol_TableRowBgAlt, (ImVec4)ImColor::HSV(7.0 / 7.0f, 0.2f, 0.4f));
							//ImGui::Text(FormatDouble(HDXApp.peptideStore->results[row].deutResults[file_list_idx].deutDelta, 2).c_str());
							isColored = true;
						}
						else
						{
							isColored = false;
						}
						
						
					}
					if (isColored)
					{
						ImGui::PopStyleColor(2);
					}
					//ImGui::Selectable(HDXApp.peptideStore->peptides[item_current_idx].sequence.c_str(), &item_current_idx);
					//ImGui::SetKeyboardFocusHere(item_current_idx);
					//ImGui::SetItemDefaultFocus();
					
					
					ImGui::EndTable();
				}
				
				ImGui::End();


			}

		}//frame ending

		// 3. Show another simple window.
		if (show_another_window)
		{
			ImGui::Begin("Parameter Window", &show_another_window);   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
			ImGui::Text("Please set correct parameters");
			ImGui::NewLine();
			if (ImGui::Button("Open Raw File List")) // open Dialog Simple
				ImGuiFileDialog::Instance()->OpenDialog("RawKey", "Choose Raw File", ".txt", ".", "");
			strcpy(rawPath, rawStr.c_str());
			strcpy(pepPath, pepStr.c_str());
			ImGui::PushItemWidth(250);
			ImGui::InputText("Raw file list path", rawPath, IM_ARRAYSIZE(rawPath));
			ImGui::NewLine();
			if (ImGui::Button("Open Peptide List File")) // open Dialog Simple
				ImGuiFileDialog::Instance()->OpenDialog("PepKey", "Choose Peptide File", ".txt", ".", "");
			DisplayFileDialog(rawStr, "RawKey");
			DisplayFileDialog(pepStr, "PepKey");
			ImGui::InputText("Peptide file list path", pepPath, IM_ARRAYSIZE(pepPath));
			ImGui::PopItemWidth();
			ImGui::NewLine();
			ImGui::PushItemWidth(50);
			
			ImGui::InputDouble("Set Reference Max Error", &refError, 0.0,0.0,"%.2f",ImGuiInputTextFlags_EnterReturnsTrue);
			
			ImGui::InputDouble("Set Deuteration Min Limit", &deutProtLimit, 0.0, 0.0, "%.2f", ImGuiInputTextFlags_EnterReturnsTrue);
				//Text("", pepPath, IM_ARRAYSIZE(pepPath));
			ImGui::NewLine();
			
			ImGui::InputDouble("Set XIC Window", &retTimeWindow, 0.0, 0.0, "%.2f", ImGuiInputTextFlags_EnterReturnsTrue);
			//Text("", pepPath, IM_ARRAYSIZE(pepPath));
			ImGui::NewLine();
			ImGui::InputDouble("Set RAW Window", &rawScanWindow, 0.0, 0.0, "%.2f", ImGuiInputTextFlags_EnterReturnsTrue);
			int rawScanWindowInt = rawScanWindow; ImGui::SameLine();
			ImGui::PushItemWidth(250);
			ImGui::SliderInt("RAW spectra", &rawScanWindowInt, 1, HDXApp.rawData->GetSpecNr()-1);
			rawScanWindow = rawScanWindowInt;

			ImGui::PopItemWidth();

			//ImGui::Begin("Parameter Window", &show_another_window);
			auto title = HDXApp.rawData->GetFileDataRT(rawScanWindow);
			ShowRAWSpectrum(to_string(title), xs3, ys3, rawCounter);


			ImGui::End();
			

		}

		if (show_mouse_window)
		{
			ShowMouseSetupPanel(&show_mouse_window);
		}
		//debuging
	
	
		// Rendering
		ImGui::Render();
		glViewport(0, 0, (int)io.DisplaySize.x, (int)io.DisplaySize.y);
		glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
		glClear(GL_COLOR_BUFFER_BIT);
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		SDL_GL_SwapWindow(window);



		//if (currentUIMsg.btn != MouseButons::None) Sleep(500);
		
	}

	// Cleanup
	try
	{
		SDL_FreeSurface(sshot);

		ImGui_ImplOpenGL3_Shutdown();
		ImGui_ImplSDL2_Shutdown();
		ImPlot::DestroyContext();
		ImGui::DestroyContext();

		SDL_GL_DeleteContext(gl_context);
		SDL_DestroyWindow(window);
		SDL_Quit();
	}
	catch (const std::exception&)
	{

	}
	

	return 0;
}