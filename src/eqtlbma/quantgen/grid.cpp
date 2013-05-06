/** \file grid.cpp
 *
 *  `Grid' is a class 
 *  Copyright (C) 2013 Timothee Flutre
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

#include "quantgen/grid.hpp"

using namespace std;

using namespace utils;

namespace quantgen {

  Grid::Grid(const string & gridFile, const bool & makeFixMaxh,
	     const int & verbose)
  {
    if (! gridFile.empty()) {
      if (verbose > 0)
	cout << "load grid in " << gridFile << " ..." << endl << flush;
    
      gzFile gridStream;
      vector<string> tokens;
      string line;
      openFile(gridFile, gridStream, "rb");
      while (getline(gridStream, line)) {
	split(line, " \t", tokens);
	if (tokens.size() != 2) {
	  cerr << "ERROR: format of file " << gridFile
	       << " should be phi2<space/tab>oma2" << endl;
	  exit (1);
	}
	phi2s.push_back(atof(tokens[0].c_str()));
	oma2s.push_back(atof(tokens[1].c_str()));
	if (makeFixMaxh) {
	  phi2s_fix.push_back(0.0);
	  oma2s_fix.push_back(atof(tokens[0].c_str())
			      + atof(tokens[1].c_str()));
	  phi2s_maxh.push_back(atof(tokens[0].c_str())
			       + atof(tokens[1].c_str()));
	  oma2s_maxh.push_back(0.0);
	}
      }
      if (! gzeof(gridStream)) {
	cerr << "ERROR: can't read successfully file " << gridFile
	     << " up to the end" << endl;
	exit(1);
      }
      closeFile(gridFile, gridStream);
    
      if (verbose > 0)
	cout << "grid size: " << phi2s.size() << endl;
    }
  }

} //namespace quantgen
