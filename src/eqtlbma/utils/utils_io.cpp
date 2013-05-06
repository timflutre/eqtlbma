/** \file utils_io.cpp
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

#include <cstdlib>
#include <cmath>
#include <cerrno>
#include <sys/stat.h>
#include <sys/time.h>
#include <dirent.h>

#include <iomanip>
#include <algorithm>
#include <sstream>

#include "utils/utils_io.hpp"

using namespace std;

namespace utils {

// http://stackoverflow.com/questions/1644868/c-define-macro-for-debug-printing/1644898#1644898
#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif
#define debug_print(fmt, ...)						\
  do { if (DEBUG_TEST) fprintf(stderr, fmt, __VA_ARGS__); } while (0)

/** \brief Split a string with one delimiter.
 *  \note http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c/236803#236803
 */
  vector<string> &
  split (
    const string & s,
    char delim,
    vector<string> & tokens)
  {
    tokens.clear();
    stringstream ss(s);
    string item;
    while(getline(ss, item, delim)) {
      tokens.push_back(item);
    }
    return tokens;
  }

/** \brief Split a string with one delimiter.
 */
  vector<string>
  split (
    const string & s,
    char delim)
  {
    vector<string> tokens;
    return split (s, delim, tokens);
  }

/** \brief Split a string with several delimiters.
 */
  vector<string> &
  split (
    const string & s,
    const char * delim,
    vector<string> & tokens)
  {
    tokens.clear();
    char * pch;
    pch = strtok ((char *) s.c_str(), delim);
    while (pch != NULL)
    {
      tokens.push_back (string(pch));
      pch = strtok (NULL, delim);
    }
    return tokens;
  }

/** \brief Split a string with several delimiters.
 */
  vector<string>
  split (
    const string & s,
    const char * delim)
  {
    vector<string> tokens;
    char * pch;
    pch = strtok ((char *) s.c_str(), delim);
    while (pch != NULL)
    {
      tokens.push_back (string(pch));
      pch = strtok (NULL, delim);
    }
    return tokens;
  }

/** \brief Split a string with several delimiters and return only the content
 *  of one token.
 */
  string
  split (
    const string & s,
    const char * delim,
    const size_t & idx)
  {
    vector<string> tokens = split (s, delim);
    if (tokens.size() < idx)
    {
      cerr << "ERROR: not enough tokens after splitting string" << endl;
      exit (1);
    }
    return (tokens[idx]);
  }

/** \brief Return the processor time in seconds consumed by the program
 *  since 'startTime'.
 */
  double
  getElapsedTime (
    const clock_t & startTime)
  {
    return ((clock () - startTime) / double(CLOCKS_PER_SEC));
  }

/** \brief Return a string with the elapsed time in d, h, m and s.
 *  \note http://stackoverflow.com/a/2419597/597069
 */
  string
  getElapsedTime (
    const time_t & startRawTime,
    const time_t & endRawTime)
  {
    double elapsed = difftime (endRawTime, startRawTime); // in sec
    float nbDays = floor (elapsed / 86400),
      nbHours = floor (fmod ((elapsed / 3600), 24)),
      nbMins = floor (fmod ((elapsed / 60), 60)),
      nbSecs = floor (fmod (elapsed, 60));
    char buffer[100];
    snprintf (buffer, 100, "%01.0fd %01.0fh %01.0fm %01.0fs",
	      nbDays, nbHours, nbMins, nbSecs);
    return string(buffer);
  }

/** \brief Return a string with Y-M-D H:M:S without end-of-line.
 */
  string
  getDateTime (
    const time_t & inTime)
  {
    struct tm * ptTm = localtime (&inTime);
    char buffer[100];
    snprintf (buffer, 126, "%i-%02i-%02i %02i:%02i:%02i",
	      1900 + ptTm->tm_year, ptTm->tm_mon + 1, ptTm->tm_mday,
	      ptTm->tm_hour, ptTm->tm_min, ptTm->tm_sec);
    return string(buffer);
  }

  void
  openFile (
    const string & pathToFile,
    ifstream & fileStream)
  {
    fileStream.open (pathToFile.c_str());
    if (! fileStream.is_open())
    {
      cerr << "ERROR: can't open file " << pathToFile << " to read ("
	   << boolalpha
	   << "fail=" << fileStream.fail()
	   << ", bad=" << fileStream.bad()
	   << noboolalpha
	   << ")" << endl;
      exit (1);
    }
  }

  void
  openFile (
    const string & pathToFile,
    ofstream & fileStream)
  {
    fileStream.open (pathToFile.c_str());
    if (! fileStream.is_open())
    {
      cerr << "ERROR: can't open file " << pathToFile << " to write ("
	   << boolalpha
	   << "fail=" << fileStream.fail()
	   << ", bad=" << fileStream.bad()
	   << noboolalpha
	   << ")" << endl;
      exit (1);
    }
  }

  void
  openFile (
    const string & pathToFile,
    gzFile & fileStream,
    const char * mode)
  {
    fileStream = gzopen (pathToFile.c_str(), mode);
    if (fileStream == NULL)
    {
      cerr << "ERROR: can't open file " << pathToFile
	   << " with mode " << *mode
	   << " (errno=" << errno << ")" << endl;
      exit (1);
    }
  }

  void
  closeFile (
    const string & pathToFile,
    ifstream & fileStream)
  {
    // http://gehrcke.de/2011/06/reading-files-in-c-using-ifstream-dealing-correctly-with-badbit-failbit-eofbit-and-perror/comment-page-1/#comment-6060
    if (fileStream.bad())
    {
      cerr << "ERROR: stream of file " << pathToFile
	   << " has badbit=true before closing" << endl;
      exit (1);
    }
    fileStream.close ();
  }

  void
  closeFile (
    const string & pathToFile,
    ofstream & fileStream)
  {
    if (! fileStream.good())
    {
      cerr << "ERROR: stream of file " << pathToFile
	   << " returns good()=false before closing" << endl
	   << boolalpha
	   << "fail()=" << fileStream.fail()
	   << " bad()=" << fileStream.bad()
	   << " eof()=" << fileStream.eof()
	   << noboolalpha << endl;
      exit (1);
    }
    fileStream.close ();
  }

  void
  closeFile (
    const string & pathToFile,
    gzFile & fileStream)
  {
    int ret = gzclose (fileStream);
    if (ret != Z_OK)
    {
      cerr << "ERROR: can't close the file " << pathToFile
	   << ", gzclose() returned " << ret << endl;
      exit (1);
    }
  }

  int
  getline (
    gzFile & fileStream,
    string & line)
  {
    int res = 1, c;
    line.clear ();
    while (true)
    {
      c = gzgetc (fileStream);
      if (c == -1) // eof or error
      {
	res = 0;
	break;
      }
      else if (c == 10) // 10 is ASCII code for '\n'
	break;
      else
	line.push_back (c);
    }
    return res;
  }

  void
  gzwriteLine (
    gzFile & fileStream,
    const string & line,
    const string & pathToFile,
    const size_t & lineId)
  {
    // if (gzprintf (fileStream, "%s", line.c_str()) <= 0)
    if (gzputs (fileStream, line.c_str()) < 0)
    {
      cerr << "ERROR: can't write line " << lineId
	   << " in file " << pathToFile << endl;
      exit (1);
    }
  }

/** \brief Used by scandir.
 *  \note unused parameter, see http://stackoverflow.com/q/1486904/597069
 */
  static int dummy_selector (const struct dirent * /*dir_entry*/)
  {
    return 1;
  }

/** \brief Return a vector with the iterations corresponding to nbSteps.
 *  \note Useful with verbose to print at which iteration a loop is.
 */
  vector<size_t>
  getCounters (
    const size_t & nbIterations,
    const size_t & nbSteps = 5)
  {
    vector<size_t> vCounters;
    size_t step = (size_t) floor (nbIterations / nbSteps);
    for (size_t i = 1; i < nbSteps; ++i)
      vCounters.push_back (i * step);
    vCounters.push_back (nbIterations);
    return vCounters;
  }

/** \brief Print the nb of iterations already complete in percentage of
 *  the total loop size.
 */
  void
  printCounter (
    const size_t & currentIter,
    const vector<size_t> & vCounters)
  {
    size_t i = 0;
    while (i < vCounters.size())
    {
      if (currentIter == vCounters[i])
      {
	printf ("%.0f%%\n", (float) 100 * currentIter / vCounters[vCounters.size()-1]);
	fflush (stdout);
	break;
      }
      ++i;
    } 
  }

/** \brief Display a progress bar on stdout.
 *  \note adapted from the GEMMA package by Xiang Zhou
 */
  void
  progressBar (
    string msg,
    double currentIter,
    double nbIterations)
  {
    double progress = (100.0 * currentIter / nbIterations);
    int barsize = (int) (progress / 2.0);
    char bar[51];
  
    cout << "\r" << msg;
    for (int i = 0; i < 50; i++)
    {
      if (i < barsize)
	bar[i] = '=';
      else
	bar[i] = ' ';
      cout << bar[i];
    }
    cout << setprecision(2) << fixed << progress << "%" << flush;
  }

/** \brief Convert int, float, etc into a string.
 *  \note http://notfaq.wordpress.com/2006/08/30/c-convert-int-to-string/
 */
  template <class T>
  inline string toString (const T & t)
  {
    stringstream ss;
    ss << t;
    return ss.str();
  }

/** \brief Copy a string into another.
 */
  string
  copyString (
    const string & input)
  {
    string output;
    for (string::const_iterator it = input.begin();
	 it != input.end();
	 ++it)
    {
      output += *it;
    }
    return output;
  }

/** \brief Replace part of a string with another string.
 *  \note http://stackoverflow.com/a/3418285/597069
 */
  void
  replaceAll (
    string & str,
    const string & from,
    const string & to)
  {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != string::npos)
    {
      str.replace(start_pos, from.length(), to);
      start_pos += to.length();  // in case 'to' contains 'from', eg. replacing 'x' with 'yx'
    }
  }

/** Return true if file exists.
 */
  bool
  doesFileExist (
    const string & filename)
  {
    bool fexists = false;
    struct stat buffer;
    fexists = ( stat(filename.c_str(), &buffer) == 0);
    return fexists;
  }

/** \brief List a given directory.
 */
  vector<string>
  scanInputDirectory (
    const string & inDir,
    const int & verbose)
  {
    vector<string> vInFiles;
    struct dirent ** inFiles = NULL;
    int nbInFiles;
    if (verbose > 0)
    {
      cout << "scan directory " << inDir << " ..." << endl;
    }
    nbInFiles = scandir(inDir.c_str(), &inFiles, dummy_selector, alphasort);
    if (nbInFiles == -1)
    {
      cerr << "ERROR: can't scan " << inDir << endl;
      exit (1);
    }
    else if (nbInFiles == 0)
    {
      cerr << "ERROR: " << inDir << " contains no file" << endl;
      exit (1);
    }
    else
    {
      for (int s = 0; s < nbInFiles; ++s)
      {
	if (string(inFiles[s]->d_name) == "." ||
	    string(inFiles[s]->d_name) == "..")
	{
	  free (inFiles[s]);
	  continue;
	}
	char path[1024];
	int nbChar;
	if (inDir[inDir.size()-1] != '/')
	  nbChar = sprintf (path, "%s/%s", inDir.c_str(), inFiles[s]->d_name);
	else
	  nbChar = sprintf (path, "%s%s", inDir.c_str(), inFiles[s]->d_name);
	if (nbChar < 0)
	{
	  cerr << "ERROR: variable 'path' is not big enough" << endl;
	}
	vInFiles.push_back (string(path));
	free (inFiles[s]);
      }
      if (verbose > 0)
	cout << "nb of files: " << vInFiles.size() << endl;
    }
    free (inFiles);
    return vInFiles;
  }

/** \brief Return true if the given path is a directory.
 */
  bool
  isDirectory (
    const char path[])
  {
    bool res = false;
    if (strlen (path) > 0)
    {
      struct stat st;
      if (stat(path, &st) != 0)
      {
	fprintf (stderr, "ERROR: stat failed for path %s\n", path);
	fprintf (stderr, "errno=%i %s\n", errno, strerror(errno));
	exit (1);
      }
      if (S_ISDIR(st.st_mode))
	res = true;
    }
    return res;
  }

  void
  createDirectory (
    const string & dirName)
  {
    if (mkdir (dirName.c_str(), 0774) != 0) // u=rwx g=rwx o=r--
    {
      cerr << "ERROR: can't create directory " << dirName
	   << " (errno=" << errno << ")" << endl;
      exit (1);
    }
  }

  void
  changeDirectory (
    const string & dirName)
  {
    if (chdir (dirName.c_str()) != 0) // u=rwx g=rwx o=r--
    {
      cerr << "ERROR: can't change directory to " << dirName
	   << " (errno=" << errno << ")" << endl;
      exit (1);
    }
  }

  string
  getCurrentDirectory (
    void)
  {
    char buf[FILENAME_MAX];
    if (getcwd (buf, sizeof(buf)) == NULL)
    {
      cerr << "ERROR: can't get current working directory (errno="
	   << errno << ")" << endl;
      exit (1);
    }
    string cwd (buf);
    return cwd;
  }

/** \brief Remove a directory even if it is not empty.
 *  \note http://stackoverflow.com/a/1149769/597069
 *  \note Don't do anything if the supplied path is empty
 *  or if the directory doesn't exist.
 */
  int
  removeDir(
    string path)
  {
    if (path.empty())
      return 0;
  
    if (path[path.size()] == '.')
      return 0;
  
    if (path[path.length()-1] != '/')
      path += "/";
  
    // create a pointer to a directory
    DIR *pdir = NULL;
    pdir = opendir (path.c_str());
    if (pdir == NULL)
    {
      if (errno == 2) // No such file or directory
	return 0;
      else
      {
	cerr << "ERROR: opendir returned NULL for path " << path << endl;
	fprintf (stderr, "errno=%i %s\n", errno, strerror(errno));
	return errno;
      }
    }
  
    struct dirent *pent = NULL;
    char file[1024];
    int counter = 1; // use this to skip the first TWO which cause an infinite loop (and eventually, stack overflow)
    while (true)
    {
      pent = readdir (pdir); // while there is still something in the directory
      if (pent == NULL)
      {
	if (errno != 0) // if pent has not been initialised correctly
	{
	  cerr << "ERROR: readdir returned NULL for path " << path << endl;
	  fprintf (stderr, "errno=%i %s\n", errno, strerror(errno));
	  return errno; // we couldn't do it
	}
	else // if the directory is empty
	  break;
      }
      if (counter > 2)
      {
	for (int i = 0; i < 256; i++)
	  file[i] = '\0';
	strcat(file, path.c_str());
	// otherwise, it was initialised correctly, so let's delete the file~
	strcat(file, pent->d_name); // concatenate the strings to get the complete path
	if (isDirectory(file) == true)
	  removeDir(file);
	else // it's a file, we can use remove
	  remove(file);
      }
      counter++;
    }
  
    // finally, let's clean up
    closedir (pdir); // close the directory
    if (rmdir(path.c_str()) != 0)
    {
      if (errno != 0)
      {
	cerr << "ERROR: rmdir returned an error" << endl;
	fprintf (stderr, "errno=%i %s\n", errno, strerror(errno));
	return errno;
      }
    }
  
    return 0;
  }

  void
  removeFiles (
    const vector<string> & vFileNames)
  {
    for (size_t i = 0; i < vFileNames.size(); ++i)
    {
      if (remove (vFileNames[i].c_str()) != 0)
      {
	cerr << "ERROR: can't remove file" << vFileNames[i] << endl;
	exit (1);
      }
    }
  }

  double
  getMaxMemUsedByProcess (void)
  {
    double vmHWM = 0.0;
    string pathToFile = "/proc/self/status";
  
    if (! doesFileExist (pathToFile))
    {
      cerr << "WARNING: " << pathToFile << " doesn't exist,"
	   << " can't track memory usage" << endl << flush;
    }
  
    string line;
    ifstream stream;
    vector<string> tokens;
    openFile (pathToFile, stream);
    while (getline (stream, line))
    {
      if (line.find("VmHWM") != string::npos)
      {
	split (line, ":", tokens);
	if (tokens.size() != 2)
	{
	  cerr << "ERROR: file " << pathToFile
	       << " has a different format" << endl;
	  exit (1);
	}
	replaceAll (tokens[1], " ", "");
	replaceAll (tokens[1], "kB", "");
	vmHWM = atof (tokens[1].c_str());
	break;
      }
    }
    closeFile (pathToFile, stream);
  
    return vmHWM;
  }

  string getMaxMemUsedByProcess2Str (void)
  {
    char str[128];
    double maxMem = getMaxMemUsedByProcess (); // in kB
    snprintf (str, 127, "%.0f kB", maxMem);
    return string(str);
  }

  string
  getCmdLine (
    int argc,
    char ** argv)
  {
    ostringstream oss;
    oss << argv[0];
    for(int i = 1; i < argc; ++i)
      oss << " " << argv[i];
    return (oss.str());
  }

} // namespace utils
