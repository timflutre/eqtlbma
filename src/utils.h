/** \file utils.h
 *
 *  `utils.cpp' gathers functions useful for any program.
 *  Copyright (C) 2011-2012  T. Flutre
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef UTILS_H
#define UTILS_H

#include <ctime>
#include <sys/stat.h>

#include <vector>
#include <map>
#include <string>
#include <fstream>
using namespace std;

#include <gsl/gsl_sort.h>
#include <gsl/gsl_cdf.h>

#include "zlib.h"

vector<string> & split (const string & s, char delim, vector<string> & tokens);

vector<string> split (const string & s, char delim);

vector<string> & split (const string & s, const char * delim,
			vector<string> & tokens);

vector<string> split (const string & s, const char * delim);

string split (const string & s, const char * delim, const size_t & idx);

string elapsedTime (const time_t & startRawTime, const time_t & endRawTime);

string time2string (const time_t & inTime);

void openFile (const string & pathToFile, ifstream & fileStream);

void openFile (const string & pathToFile, ofstream & fileStream);

void openFile (const string & pathToFile, gzFile & fileStream,
	       const char * mode);

void closeFile (const string & pathToFile, ifstream & fileStream);

void closeFile (const string & pathToFile, ofstream & fileStream);

void closeFile (const string & pathToFile, gzFile & fileStream);

int getline (gzFile & fileStream, string & line);

void gzwriteLine (gzFile & fileStream, const string & line,
		  const string & pathToFile, const size_t & lineId);

vector<string> loadOneColumnFile (const string & inFile,
				  const int & verbose);

map<string, string> loadTwoColumnFile (const string & inFile,
				       const int & verbose);

void loadTwoColumnFile (const string & inFile, map<string, string> & mItems,
			vector<string> & vKeys, const int & verbose);

vector<size_t> loadOneColumnFileAsNumbers (const string & inFile,
					   const int & verbose);

vector<size_t> getCounters (const size_t & nbIterations,
			    const size_t & nbSteps);

void printCounter (const size_t & currentIter,
		   const vector<size_t> & vCounters);

void progressBar (string msg, double currentIter, double nbIterations);

template <class T> inline string toString (const T & t);

string copyString (const string input);

void replaceAll (string & str, const string & from, const string & to);

double round (double x);

bool doesFileExist (const string & filename);

vector<string> scanInputDirectory (const string & inDir, const int & verbose);

bool isDirectory (const char path[]);

void createDirectory (const string & dirName);

void changeDirectory (const string & dirName);

string getCurrentDirectory (void);

int removeDir (string path);

void removeFiles (const vector<string> & vFileNames);

size_t getSeed (void);

void qqnorm (double * ptData, const size_t & n);

double log10_weighted_sum (const double * vec, const size_t size);

double log10_weighted_sum (const double * vec, const double * weights,
			   const size_t size);

bool isNonZero (size_t i);

bool isNonNpos (size_t i);

bool isNan (double i);

double getMaxMemUsedByProcess (void);

string getMaxMemUsedByProcess2Str (void);

void printCmdLine (ostream & os, int argc, char ** argv);

#endif
