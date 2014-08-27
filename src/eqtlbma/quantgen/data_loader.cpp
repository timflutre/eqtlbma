/** \file data_loader.cpp
 *
 *  `data_loader' gathers functions useful to load gen-etics/-omics data
 *  Copyright (C) 2011-2014 Timothée Flutre
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

#include "quantgen/data_loader.hpp"

using namespace std;

using namespace utils;

namespace quantgen {
  
  map<string,string> loadTwoColumnFile(const string & file,
				       const int & verbose)
  {
    map<string,string> mItems;
    if(file.empty())
      return mItems;
    
    string line;
    gzFile stream;
    vector<string> tokens;
    size_t nb_lines = 0;
    
    openFile(file, stream, "rb");
    if(verbose > 0)
      cout <<"load file " << file << " ..." << endl;
    
    while(getline(stream, line)){
      nb_lines++;
      split(line, " \t,", tokens);
      if(tokens.size() != 2){
	cerr << "ERROR: file " << file << " should have only two columns"
	     << " at line " << nb_lines << endl;
	exit(EXIT_FAILURE);
      }
      if(tokens[0][0] == '#')
	continue;
      if(mItems.find(tokens[0]) == mItems.end())
	mItems.insert(make_pair(tokens[0], tokens[1]));
    }
    
    if(! gzeof(stream)){
      cerr << "ERROR: can't read successfully file "
	   << file << " up to the end" << endl;
      exit (1);
    }
    closeFile(file, stream);
    
    if(verbose > 0)
      cout << "items loaded: " << mItems.size() << endl;
    
    return mItems;
  }
  
} // namespace quantgen
