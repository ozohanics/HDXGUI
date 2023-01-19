/***************************************************************************
 *  file                 :  ipc.c                                          *
 *  copyright            : (C) 2001-2005 by Dirk Nolting                   *
 *  email                : nolting@uni-duesseldorf.de                      *
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "globals.h"
#include <signal.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <ppl.h>

#include "constants.h"

int TCompounds::is_symbol(char *probe)
{
	element *cur;

	cur = elements;
	while (cur)
	{
		if (!strcmp(probe, cur->symbol))
			return 1;
		cur = cur->next;
	}
	return 0;
}

int TCompounds::add_component(char *symbol, int number)
{
	element *el;
	compound *newco = NULL, *co;

	el = elements;
	while (el)
	{
		if (!strcmp(symbol, el->symbol))
		{
			newco = (compound *)malloc(sizeof(compound));
			newco->isotopes = el->isotopes;
			newco->amount = number;
			newco->next = NULL;
		}
		el = el->next;
	}

	if (!newco)
	{
		printf("Unknown element: %s. Check input or file elemente\n", symbol);
		return 0;
	}

	co = verbindung;
	if (!verbindung)
	{
		newco->previous = NULL;
		verbindung = newco;
		return 1;
	}

	while (co)
	{
		if (co->isotopes == newco->isotopes)
		{
			co->amount += newco->amount;
			free(newco);
			newco = nullptr;
			return 1;
		}
		co = co->next;
	}

	co = verbindung;
	while (co->next)
		co = co->next;
	co->next = newco;
	newco->previous = co;

	return 1;
}

int TCompounds::pars_chem_form(char *formel, int len)
{
	char par[MAX_DIGITS], par1[MAX_DIGITS];
	int m = 0, number = 0;

	if (!formel) return 0;

	int t = 0;
	while (t < len)		
	{
		int a = 0;
		while (isalpha(formel[t]))
		{			
			par[a] = formel[t];
			par[a + 1] = '\0';
			a++;
			t++;
		}

		m = 0;
		while (isdigit(formel[t]))
		{
			par1[m] = formel[t];
			m++;
			t++;		
		}
		
		//store data
		par1[m] = '\0';
		number = atoi(par1);

		//WRITE_INFO(par);
		//WRITE_INFO(to_string(number));
		if (!is_symbol(par)) return 0;

		/* Adding omitted last 1, e.g. CH3Cl -> CH3Cl1 */
		if (!number) number = 1;
		if (!add_component(par, number)) return 0;
	};
	return 1;
}

int TCompounds::add_amino_acid(char acid)
{
	switch (acid)
	{
	case 'A':
		add_component(&atC, 3);
		add_component(&atH, 5);
		add_component(&atN, 1);
		add_component(&atO, 1);
		break;
	case 'R':
		add_component(&atC, 6);
		add_component(&atH, 12);
		add_component(&atN, 4);
		add_component(&atO, 1);
		break;
	case 'N':
		add_component(&atC, 4);
		add_component(&atH, 6);
		add_component(&atN, 2);
		add_component(&atO, 2);
		break;
	case 'B':
	case 'D':
		add_component(&atC, 4);
		add_component(&atH, 5);
		add_component(&atN, 1);
		add_component(&atO, 3);
		break;
	case 'C':
		add_component(&atC, 3);
		add_component(&atH, 5);
		add_component(&atN, 1);
		add_component(&atO, 1);
		add_component(&atS, 1);
		break;
	case 'E':
		add_component(&atC, 5);
		add_component(&atH, 7);
		add_component(&atN, 1);
		add_component(&atO, 3);
		break;
	case 'Z':
	case 'Q':
		add_component(&atC, 5);
		add_component(&atH, 8);
		add_component(&atN, 2);
		add_component(&atO, 2);
		break;
	case 'X':
	case 'G':
		add_component(&atC, 2);
		add_component(&atH, 3);
		add_component(&atN, 1);
		add_component(&atO, 1);
		break;
	case 'H':
		add_component(&atC, 6);
		add_component(&atH, 7);
		add_component(&atN, 3);
		add_component(&atO, 1);
		break;
	case 'I':
	case 'L':
		add_component(&atC, 6);
		add_component(&atH, 11);
		add_component(&atN, 1);
		add_component(&atO, 1);
		break;
	case 'K':
		add_component(&atC, 6);
		add_component(&atH, 12);
		add_component(&atN, 2);
		add_component(&atO, 1);
		break;
	case 'M':
		add_component(&atC, 5);
		add_component(&atH, 9);
		add_component(&atN, 1);
		add_component(&atO, 1);
		add_component(&atS, 1);
		break;
	case 'F':
		add_component(&atC, 9);
		add_component(&atH, 9);
		add_component(&atN, 1);
		add_component(&atO, 1);
		break;
	case 'P':
		add_component(&atC, 5);
		add_component(&atH, 7);
		add_component(&atN, 1);
		add_component(&atO, 1);
		break;
	case 'S':
		add_component(&atC, 3);
		add_component(&atH, 5);
		add_component(&atN, 1);
		add_component(&atO, 2);
		break;
	case 'T':
		add_component(&atC, 4);
		add_component(&atH, 7);
		add_component(&atN, 1);
		add_component(&atO, 2);
		break;
	case 'W':
		add_component(&atC, 11);
		add_component(&atH, 10);
		add_component(&atN, 2);
		add_component(&atO, 1);
		break;
	case 'Y':
		add_component(&atC, 9);
		add_component(&atH, 9);
		add_component(&atN, 1);
		add_component(&atO, 2);
		break;
	case 'V':
		add_component(&atC, 5);
		add_component(&atH, 9);
		add_component(&atN, 1);
		add_component(&atO, 1);
		break;
	default:
		printf("Unknown symbol: %c\n", acid);
		return 0;
	}
	return 1;
}

int TCompounds::pars_amino_acid(char *formel)
{
	while (*formel)
	{
		if (!add_amino_acid(*formel)) return 0;
		//V532 160614 OO added ()
		(*formel)++;
	}
	add_component(&atO, 1);
	add_component(&atH, 2);

	return 1;
}

int TCompounds::pars_peptid(char *formel)
{
	FILE *peptid_file;
	char linebuffer[MAX_PEP_LINE];
	int index = 0;

	if (!(peptid_file = fopen(formel, "r")))
	{
		printf("Can´t open file: %s\n", formel);
		return 0;
	}

	while (fgets(linebuffer, MAX_PEP_LINE, peptid_file))
	{
		while (linebuffer[index] != '\n')
		{
			add_amino_acid(linebuffer[index]);
			index++;
		}
		index = 0;
	}

	add_component(&atO, 1);
	add_component(&atH, 2);

	fclose(peptid_file);

	return 1;
}

int TCompounds::print_sum()
{
	element *el;
	compound *co;

	printf("\nChemical formula: ");
	co = verbindung;

	while (co)
	{
		el = elements;
		while (el && (co->isotopes != el->isotopes))
			el = el->next;

		if (!el) return 0;
		printf("%s%i ", el->symbol, co->amount);
		co = co->next;
	}

	printf("\n");
	return 1;
}

FILE* TCompounds::open_file(char *filename)
{
	FILE *datenfile;
	if (!(datenfile = fopen(filename, "r")))
	{
		printf("Ca&atNt open file: %s.\n", filename);
		return 0;
	}

	return datenfile;
}

void TCompounds::add_element(element *ce)
{
	element *cur;
	ce->next = NULL;
	if (!(elements))
	{
		elements = ce;
		elements->previous = NULL;
		return;
	}

	cur = elements;
	while (cur->next)
	{
		cur = cur->next;
		/*printf("Element %d\n", nr);
		if (nr > NUM_ELEMENTS) break;*/
	}

	/*	if (cur->previous) printf("Element %s after %s\n", cur->symbol, cur->previous->symbol);*/
	ce->previous = cur;
	cur->next = ce;
}

void TCompounds::add_isotope(isotope *ci, element *ce)
{
	isotope *cur;
	/*if (ci) printf("Element %s with Isotope: %4f\n", ce->symbol,ci->mass);*/
	ci->next = NULL;
	if (!(ce->isotopes))
	{
		ce->isotopes = ci;
		ci->previous = NULL;
		return;
	}

	cur = ce->isotopes;

	while (cur->next)
	{
		cur = cur->next;
	}

	/*if (cur->previous) printf("Isotope %4f after %4f\n", cur->mass, cur->previous->mass);*/

	cur->next = ci;
	ci->previous = cur;

	return;
}

int TCompounds::init_elements()
{
	FILE* data;
	char linebuffer[MAX_LINE];
	int count;
	element *ce;
	isotope *ci;
	data = fopen(ELEMENTFILE, "r");
	if (!data)
	{
		printf("Unable to open elements file! \n");
		fclose(data);
		return 0;
	}

	fgets(linebuffer, MAX_LINE, data);  /* read table headings*/
	fgets(linebuffer, MAX_LINE, data);
	fgets(linebuffer, MAX_LINE, data);  /* get first barred line*/

	for (count = 0; count < NUM_ELEMENTS; count++) {
		ce = parse_element(data, linebuffer);
		add_element(ce);
		ci = parse_isotope(data, linebuffer);
		/* do&atNt add isotopes with 0 probability since they do&atNt affect results
		 but increase computing time.*/
		if (ci->p > 0.0) {
			add_isotope(ci, ce);
		}

		while (fgetc(data) != '_') {  /*if the first char is&atNt a bar then this is an isotope line*/
			ci = parse_isotope(data, linebuffer);
			/*printf("isotope %4f ... %4f\n", ci->p, ci->mass );*/
			add_isotope(ci, ce);
		}
		fgets(linebuffer, MAX_LINE, data);  /* otherwise read the rest of the bar line*/
	}

	fclose(data);
	return 1;
}

element* TCompounds::parse_element(FILE* data, char* linebuffer) {
	char* token;
	int number;
	element* ce = (element*)malloc(sizeof(element));
	ce->isotopes = NULL;

	/* read enough characters to get the element number and symbol*/
	fgets(linebuffer, 9, data);

	token = (char*)strtok(linebuffer, " ");
	number = atoi(token);
	token = (char*)strtok(NULL, " ");
	ce->symbol = _strdup(token);

	/*printf("Element symbol: %s\tElement number: %d\n",ce->symbol,number);*/
	return ce;
}

isotope* TCompounds::parse_isotope(FILE* data, char* linebuffer) {
	char* token;
	int number;
	isotope* ci = (isotope *)malloc(sizeof(isotope));

	fgets(linebuffer, MAX_LINE, data);

	token = (char*)strtok(linebuffer, " ");   /*get the number in 3rd column*/
	number = atoi(token);

	token = (char*)strtok(NULL, "( ");  /* get the mass*/
	ci->mass = atof(token);

	token = (char*)strtok(NULL, ")# "); /* get the number in paretheses and remove*/
	/* the ')', a '#',and any trailing spaces*/

	token = (char*)strtok(NULL, "(");   /*get the isotopic composition*/
	ci->p = (atof(token)) / 100;

	/*printf("\tisotope mass: %f\tisotope percent: %f\n",ci->mass,ci->p);*/
	return ci;
}

void TCompounds::free_list(isotope *target)
{
	while (target->next)
	{
		target = target->next;
		free(target->previous);
	}
	free(target);
}

void TCompounds::cut_peaks(isotope *spectrum)
{
	int dummy = 1;

	while ((spectrum->next) && (dummy < fast_calc))
	{
		dummy++;
		spectrum = spectrum->next;
	}

	if (spectrum->next)
	{
		free_list(spectrum->next);
		spectrum->next = NULL;
	}
}

void TCompounds::summarize_peaks()
{
	isotope *dummy, *d2;

	for (dummy = peaks; dummy; dummy = dummy->next)
		/* Differenz wegen Rundungsfehlern */
		while (dummy->next && (dummy->next->mass - dummy->mass < FWHM_limit))
			/*Felbontas allitast megirni*/
		{
			d2 = dummy->next;
			dummy->next = d2->next;
			if (dummy->next)
				dummy->next->previous = dummy;
			dummy->mass = (dummy->mass*dummy->p + d2->mass*d2->p) / (dummy->p + d2->p);
			dummy->p += d2->p;

			free(d2);
		}
}

isotope* TCompounds::add_peak(isotope *base, isotope *peak)
{
	static isotope *reiter;

	if (!(base->mass))
	{
		peak->next = NULL;
		peak->previous = NULL;
		reiter = peak;
		return peak;
	}

	if (peak->mass >= reiter->mass)
		while ((reiter->next) && (reiter->mass < peak->mass))
			reiter = reiter->next;
	else
	{
		while ((reiter->previous) && (reiter->mass > peak->mass))
			reiter = reiter->previous;
		reiter = reiter->next;
	}

	if ((reiter->mass) >= (peak->mass))
	{
		peak->next = reiter;
		peak->previous = reiter->previous;
		peak->previous->next = peak;
		reiter->previous = peak;
		return base;
	}
	else
	{
		reiter->next = peak;
		peak->next = NULL;
		peak->previous = reiter;
		return base;
	}
	return 0;
}

int TCompounds::calculate_peaks(){
	compound *c;
	isotope *npeaks, *p, *i, *np1;
	size_t anzahl;

	if (!(peaks = (isotope*)malloc(sizeof(isotope))))
		return 0;
	peaks->mass = 0;
	peaks->p = 1;
	peaks->previous = NULL;
	peaks->next = NULL;


	for (c = verbindung; c; c = c->next){
		//parallel_for (size_t(0), c->amount, [&](size_t anzahl)
		for (anzahl = 0; anzahl < c->amount; anzahl++)
		{
			if (!(npeaks = (isotope*)malloc(sizeof(isotope))))
				return 0;
			npeaks->mass = 0;
			for (p = peaks; p; p = p->next){
				for (i = c->isotopes; i; i = i->next){
					if (!(np1 = (isotope*)malloc(sizeof(isotope))))
						return 0;

					np1->mass = p->mass + i->mass;
					np1->p = p->p * i->p;
					if (!(npeaks = add_peak(npeaks, np1)))
						return 0;
				}
			}
			free_list(peaks);
			peaks = npeaks;
			summarize_peaks();
			if (fast_calc)
				cut_peaks(peaks);
		}
		//);
	}
	return 1;
}

void TCompounds::print_result(){
	isotope *d;
	double maxp = 0, relint = 0, sump = 0;
	int permutationen = 0;

	printf("\n");

	for (d = peaks; d; d = d->next)
	{
		permutationen++;
		sump += d->p;
	}

	summarize_peaks();
	for (d = peaks; d; d = d->next)
		if (d->p > maxp)
			maxp = d->p;

	for (d = peaks; d; d = d->next)
	{
		if ((relint = (d->p / maxp) * 100) > MIN_INT)
			printf("M= %f, p= %e, rel. Int.= %f%%\n",
			d->mass, d->p, relint);
	}
	if (!(fast_calc))
		printf("\nNumber of  permutations: %i\n", permutationen);
	else
	{
		sump = (rint(sump * 10000) / 100);
		printf("\nCovered Intensity: %2.2f%% \n", sump);
	}
}

int TCompounds::StorePeaks(double* massOut, double* intOut)
{
	isotope *d;
	double maxp = 0, relint = 0, sump = 0;
	int permutationen = 0;

	for (d = peaks; d; d = d->next)
	{
		permutationen++;
		sump += d->p;
	}

	summarize_peaks();
	for (d = peaks; d; d = d->next)
		if (d->p > maxp)
			maxp = d->p;

	int peaknr = 0;
	for (d = peaks; d; d = d->next)
	{
		if ((relint = (d->p / maxp)) > MIN_INT / 100)
		{
			massOut[peaknr] = d->mass;
			intOut[peaknr] = relint;

			peaknr++;
		}
		//modification to store just 7 isotopes
		if (peaknr > 7) break;
		//			d->mass, d->p, relint;
	}
	if (fast_calc)
	{
		sump = (rint(sump * 10000) / 100);
	}
	return peaknr;
}
