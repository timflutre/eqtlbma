/** \file utils_io.hpp
 *
 *  `utils_io' gathers functions useful for any program.
 *  Copyright (C) 2011-2013 Timothee Flutre
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

#ifndef UTILS_UTILS_IO_HPP
#define UTILS_UTILS_IO_HPP

#include <cstdlib>
#include <ctime>

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>

#include "zlib.h"

namespace utils {

  std::vector<std::string> & split (const std::string & s, char delim, std::vector<std::string> & tokens);

  std::vector<std::string> split (const std::string & s, char delim);

  std::vector<std::string> & split (const std::string & s, const char * delim,
			  std::vector<std::string> & tokens);

  std::vector<std::string> split (const std::string & s, const char * delim);

  std::string split (const std::string & s, const char * delim, const size_t & idx);

  double getElapsedTime (const clock_t & startTime);

  std::string getElapsedTime (const time_t & startRawTime, const time_t & endRawTime);

  std::string getDateTime (const time_t & inTime);

  void openFile (const std::string & pathToFile, std::ifstream & fileStream);

  void openFile (const std::string & pathToFile, std::ofstream & fileStream);

  void openFile (const std::string & pathToFile, gzFile & fileStream,
		 const char * mode);

  void closeFile (const std::string & pathToFile, std::ifstream & fileStream);

  void closeFile (const std::string & pathToFile, std::ofstream & fileStream);

  void closeFile (const std::string & pathToFile, gzFile & fileStream);

  int getline (gzFile & fileStream, std::string & line);

  void gzwriteLine (gzFile & fileStream, const std::string & line,
		    const std::string & pathToFile, const size_t & lineId);

  std::vector<size_t> getCounters (const size_t & nbIterations,
			      const size_t & nbSteps);

  void printCounter (const size_t & currentIter,
		     const std::vector<size_t> & vCounters);

  void progressBar (std::string msg, double currentIter, double nbIterations);

  template <class T> inline std::string toString (const T & t);

  std::string copyString (const std::string input);

  void replaceAll (std::string & str, const std::string & from, const std::string & to);

  bool doesFileExist (const std::string & filename);

  std::vector<std::string> scanInputDirectory (const std::string & inDir, const int & verbose);

  bool isDirectory (const char path[]);

  void createDirectory (const std::string & dirName);

  void changeDirectory (const std::string & dirName);

  std::string getCurrentDirectory (void);

  int removeDir (std::string path);

  void removeFiles (const std::vector<std::string> & vFileNames);

  double getMaxMemUsedByProcess (void);

  std::string getMaxMemUsedByProcess2Str (void);

  std::string getCmdLine (int argc, char ** argv);

  /** \brief Fill a vector with the keys of a map
   *  \note http://stackoverflow.com/a/771463/597069
   *  \note http://stackoverflow.com/a/10632266/597069
   */
  template <typename M, typename V>
  void keys2vec(const M & m, V & v)
  {
    for(typename M::const_iterator it = m.begin(); it != m.end(); ++it)
      v.push_back(it->first);
  }

} // namespace utils

#endif // UTILS_UTILS_IO_HPP
