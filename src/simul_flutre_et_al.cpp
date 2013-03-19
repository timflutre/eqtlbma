/** \file simul.cpp
 *
 *  `simul' generates data to test eQtlBma.
 *  Copyright (C) 2012-2013 Timothee Flutre
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
 *
 *  g++ -Wall -Wextra -g utils.cpp simul_flutre_et_al.cpp -lgsl -lgslcblas -lz -o simul_flutre_et_al
 */

#include <cmath>
#include <ctime>
#include <getopt.h>
#include <sys/stat.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <limits>
using namespace std;

#include "utils.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifndef VERSION
#define VERSION "0.0"
#endif

/** \brief Display the help on stdout.
 */
void help (char ** argv)
{
  cout << "`" << argv[0] << "'"
       << " generates data to test eQtlBma." << endl
       << endl
       << "Usage: " << argv[0] << " [OPTIONS] ..." << endl
       << endl
       << "Options:" << endl
       << "  -h, --help\tdisplay the help and exit" << endl
       << "  -V, --version\toutput version information and exit" << endl
       << "  -v, --verbose\tverbosity level (0/default=1/2/3)" << endl
       << "      --name\tname of the scenario (optional)" << endl
       << "      --tissues\tnumber of tissues (default=1)" << endl
       << "      --scenario\twhich set of configurations the data are simulated" << endl
       << "\t\t1 (default): from global null to fully consistent (eg. if S=3: 000,100,110,111)" << endl
       << "\t\t2: null, consistent, each tissue-specific, 30% activity (eg. if S=5: 11000)" << endl
       << "\t\t and complementary 70% activity (00111), but requires --tissues >= 5" << endl
       << "      --repl\tnumber of replicates (default=1)" << endl
       << "      --seed\tmain seed (default=1859)" << endl
       << "      --err-var\tchoice of the standard deviations of the errors" << endl
       << "\t\t1: default, they are all equal to 1" << endl
       << "\t\t2: they are slightly different (uniformly 1, 1.5 or 2)" << endl
       << "\t\t3: they are 'very' different (uniformly 1, 2 or 4)" << endl
       << "      --het\tchoice of the heterogeneity phi^2 / (phi^2 + omega^2)" << endl
       << "\t\t1: default, no heterogeneity (0)" << endl
       << "\t\t2: low heterogeneity (0.2)" << endl
       << "\t\t3: medium heterogeneity (0.5)" << endl
       << "\t\t4: high heterogeneity (0.8)" << endl
       << "      --pve\tproportion of variance explained (default=0.2)" << endl
       << "      --maf\tminor allele frequency (default=0.3)" << endl
       << "      --inds\tnumber of individuals (default=100)" << endl
       << "      --genes\tnumber of genes per configuration (default=1000)" << endl
       << endl;
}

/** \brief Display version and license information on stdout.
 */
void version (char ** argv)
{
  cout << argv[0] << " " << __DATE__ << " " << __TIME__ << endl
       << endl
       << "Copyright (C) 2012-2013 Timothee Flutre." << endl
       << "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>" << endl
       << "This is free software; see the source for copying conditions.  There is NO" << endl
       << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl
       << endl
       << "Written by Timothee Flutre." << endl;
}

/** \brief Parse the command-line arguments and check the values of the 
 *  compulsory ones.
 */
void
parseArgs (
  int argc,
  char ** argv,
  string & simName,
  size_t & nbTissues,
  size_t & scenario,
  size_t & nbReplicates,
  size_t & mainSeed,
  size_t & choiceErrVar,
  size_t & choiceHet,
  double & pve,
  double & maf, 
  size_t & nbInds,
  size_t & nbGenesPerConfig,
  int & verbose)
{
  int c = 0;
  while (true)
  {
    static struct option long_options[] =
    {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'V'},
      {"verbose", required_argument, 0, 'v'},
      {"name", required_argument, 0, 0},
      {"scenario", required_argument, 0, 0},
      {"tissues", required_argument, 0, 0},
      {"repl", required_argument, 0, 0},
      {"seed", required_argument, 0, 0},
      {"err-var", required_argument, 0, 0},
      {"het", required_argument, 0, 0},
      {"pve", required_argument, 0, 0},
      {"maf", required_argument, 0, 0},
      {"inds", required_argument, 0, 0},
      {"genes", required_argument, 0, 0},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    c = getopt_long (argc, argv, "hVv:",
                     long_options, &option_index);
    if (c == -1)
      break;
    switch (c)
    {
    case 0:
      if (long_options[option_index].flag != 0)
        break;
      if (strcmp(long_options[option_index].name, "name") == 0)
      {
	simName = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "scenario") == 0)
      {
	scenario = atol (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "tissues") == 0)
      {
	nbTissues = atol (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "repl") == 0)
      {
	nbReplicates = atol (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "seed") == 0)
      {
	mainSeed = atol (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "err-var") == 0)
      {
	choiceErrVar = atol (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "het") == 0)
      {
	choiceHet = atol (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "pve") == 0)
      {
	pve = atof (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "maf") == 0)
      {
	maf = atof (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "inds") == 0)
      {
	nbInds = atol (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "genes") == 0)
      {
	nbGenesPerConfig = atol (optarg);
	break;
      }
    case 'h':
      help (argv);
      exit (0);
    case 'V':
      version (argv);
      exit (0);
    case 'v':
      verbose = atoi(optarg);
      break;
    case '?':
      printf ("\n"); help (argv);
      abort ();
    default:
      printf ("\n"); help (argv);
      abort ();
    }
  }
  if (scenario == 2 & nbTissues < 5)
  {
    cerr << "ERROR: --configs 2 requires --tissues >= 5" << endl;
    exit (1);
  }
}

class Params
{
  string _nameOfProjectDir; // will contain each replicate directory
  string _simName; // optional
  size_t _nbTissues;
  size_t _scenario;
  vector<vector<size_t> > _vvConfigs;
  size_t _nbReplicates;
  size_t _mainSeed;
  vector<size_t> _vRepSeeds;
  size_t _choiceErrVar;
  vector<double> _vSigmas;
  size_t _choiceHet;
  double _het;
  double _pve;
  double _maf;
  size_t _nbInds;
  size_t _nbGenesPerConfig;
  size_t _nbGenes;
  size_t _nbSnps;
  vector<vector<size_t> > _vvGenotypes;
public:
  Params (const string & simName,
	  const size_t & nbTissues,
	  const size_t & scenario,
	  const size_t & nbReplicates,
	  const size_t & mainSeed,
	  const size_t & choiceErrVar, const size_t & choiceHet,
	  const double & pve, const double & maf, const size_t & nbInds,
	  const size_t & nbGenesPerConfig);
  void drawSigma (const size_t & s, const gsl_rng * rng);
  void write (const size_t & r);
  void simulAndWriteGeneCoords (const int & verbose);
  void simulAndWriteSnpCoords (const int & verbose);
  void simulAndWriteOrLoadGenotypes (const gsl_rng * rng, const int & verbose);
  void simulAndWritePhenotypes (const gsl_rng * rng, const int & verbose);
  void simulAndWritePhenoForGivenRep (const string & truthFileName,
				      gzFile & truthFileStream,
				      stringstream & ssTruth,
				      const vector<string> & vPhenoFileNames,
				      const gsl_rng * rng);
};

Params::Params (
  const string & simName,
  const size_t & nbTissues,
  const size_t & scenario,
  const size_t & nbReplicates,
  const size_t & mainSeed,
  const size_t & choiceErrVar,
  const size_t & choiceHet,
  const double & pve,
  const double & maf, 
  const size_t & nbInds,
  const size_t & nbGenesPerConfig)
{
  _simName = simName;
  _nbTissues = nbTissues;
  
  _scenario = scenario;
  if (_scenario == 1)
  {
    _vvConfigs.resize (_nbTissues, vector<size_t> (_nbTissues, 0));
    for (size_t i = 0; i < _vvConfigs.size(); ++i)
      for (size_t j = i; j < _vvConfigs[i].size(); ++j)
	_vvConfigs[i][j] = 1;
  }
  else if (_scenario == 2)
  {
    _vvConfigs.resize (1 + _nbTissues + 2,
		       vector<size_t> (_nbTissues, 0));
    // _vvConfigs[0] is the fully consistent config
    // then each tissue-specific,
    // and finally the complementary 30% and 70%
    size_t boundary = (size_t) ceil ((double) _nbTissues / 3.0);
    for (size_t s = 0; s < _nbTissues; ++s)
    {
      _vvConfigs[0][s] = 1;
      _vvConfigs[1+s][s] = 1;
      if (s < boundary)
	_vvConfigs[1+_nbTissues][s] = 1;
      else
	_vvConfigs[1+_nbTissues+1][s] = 1;
    }
  }
  
  _nbReplicates = nbReplicates;
  _mainSeed = mainSeed;
  _vRepSeeds = vector<size_t> (_nbReplicates, 0);
  
  _choiceErrVar = choiceErrVar;
  _vSigmas = vector<double> (_nbTissues, 1);
  
  _choiceHet = choiceHet;
  _het = 0;
  if (_choiceHet == 2)
    _het = 0.2;
  else if (_choiceHet == 3)
    _het = 0.5;
  else if (_choiceHet == 4)
    _het = 0.8;
  
  _pve = pve;
  _maf = maf;
  _nbInds = nbInds;
  
  _nbGenesPerConfig = nbGenesPerConfig;
  _nbGenes = _nbGenesPerConfig * (_vvConfigs.size() + 1);
  _nbSnps = _nbGenes;
  
  stringstream ssNameOfProjectDir;
  ssNameOfProjectDir << "simul";
  if (! _simName.empty())
    ssNameOfProjectDir << "_" << _simName;
  ssNameOfProjectDir << "_tissues-" << _nbTissues
		     << "_scenario-" << _scenario;
  if (_choiceErrVar == 1)
    ssNameOfProjectDir << "_err-var-same";
  else if (_choiceErrVar == 2)
    ssNameOfProjectDir << "_err-var-diff-small";
  else if (_choiceErrVar == 3)
    ssNameOfProjectDir << "_err-var-diff-large";
  if (_choiceHet == 1)
    ssNameOfProjectDir << "_het-no";
  if (_choiceHet == 2)
    ssNameOfProjectDir << "_het-low";
  if (_choiceHet == 3)
    ssNameOfProjectDir << "_het-med";
  if (_choiceHet == 4)
    ssNameOfProjectDir << "_het-high";
  ssNameOfProjectDir << "_pve-" << (size_t) floor (_pve * 10)
		     << "_maf-" << (size_t) floor (_maf * 10);
  _nameOfProjectDir = ssNameOfProjectDir.str();
  if (doesFileExist (_nameOfProjectDir))
  {
    cerr << "ERROR: directory " << _nameOfProjectDir << " already exists"
	 << endl;
    exit (1);
  }
}

void
Params::drawSigma (
  const size_t & s,
  const gsl_rng * rng)
{
  if (_choiceErrVar == 2)
  {
    // _vSigmas[s] = 1 / gsl_ran_gamma (rng, 3, 0.5);
    double x = gsl_rng_uniform (rng);
    if (x < 0.3)
      _vSigmas[s] = 1;
    else if (x >= 0.3 & x < 0.6)
      _vSigmas[s] = 1.5;
    else
      _vSigmas[s] = 2;
  }
  else if (_choiceErrVar == 3)
  {
    // _vSigmas[s] = 1 / gsl_ran_gamma (rng, 5, 0.1);
    double x = gsl_rng_uniform (rng);
    if (x < 0.3)
      _vSigmas[s] = 1;
    else if (x >= 0.3 & x < 0.6)
      _vSigmas[s] = 2;
    else
      _vSigmas[s] = 4;
  }
}

void
Params::write (
  const size_t & r)
{
  string fileName = "params.txt";
  ofstream fileStream;
  openFile (fileName, fileStream);
  
  if (! _simName.empty())
    fileStream << "sim.name " << _simName << endl;
  fileStream << "nb.tissues " << _nbTissues
	     << endl
	     << "set.configs " << _scenario
	     << endl
	     << "choice.errvar " << _choiceErrVar
	     << endl
	     << "sigmas";
  for (size_t s = 0; s < _vSigmas.size(); ++s)
    fileStream << " " << _vSigmas[s];
  fileStream << endl;
  fileStream << "choice.het " << _choiceHet
	     << endl
	     << "het " << _het
	     << endl
	     << "pve " << _pve
	     << endl
	     << "maf " << _maf
	     << endl
	     << "nb.inds " << _nbInds
	     << endl
	     << "nb.genes.per.config " << _nbGenesPerConfig
	     << endl
	     << "nb.genes " << _nbGenes
	     << endl
	     << "nb.snps " << _nbSnps
	     << endl
	     << "rep.id " << (r+1)
	     << endl
	     << "rep.seed " << _vRepSeeds[r]
	     << endl;
  
  closeFile (fileName, fileStream);
}

void
Params::simulAndWriteGeneCoords (const int & verbose)
{
  stringstream ssFileName;
  ssFileName << "gene_coords";
  if (! _simName.empty())
    ssFileName << "_" << _simName;
  ssFileName << "_tissues-" << _nbTissues
	     << "_scenario-" << _scenario
	     << ".bed.gz";
  if (! doesFileExist (ssFileName.str()))
  {
    if (verbose > 0)
      cout << "write file " << ssFileName.str() << endl;
    gzFile fileStream;
    openFile (ssFileName.str(), fileStream, "wb");
    
    stringstream ssTxt;
    for (size_t g = 0; g < _nbGenes; ++g)
    {
      ssTxt.str("");
      ssTxt << "chr1"
	    << "\t" << g * 50
	    << "\t" << g * 50 + 11
	    << "\tgene" << (g+1)
	    << endl;
      gzwriteLine (fileStream, ssTxt.str(), ssFileName.str(), (g+1));
    }
    
    closeFile (ssFileName.str(), fileStream);
  }
  else
  {
    if (verbose > 0)
      cout << "read file " << ssFileName.str() << endl;
    gzFile fileStream;
    openFile (ssFileName.str(), fileStream, "rb");
    
    string line;
    size_t nbLines = 0;
    while (getline (fileStream, line))
      ++nbLines;
    if (nbLines != _nbGenes)
    {
      cerr << "ERROR: different nb of genes in file " << ssFileName.str() << endl;
      exit (1);
    }
    
    closeFile (ssFileName.str(), fileStream);
  }
}

void
Params::simulAndWriteSnpCoords (const int & verbose)
{
  stringstream ssFileName;
  ssFileName << "snp_coords";
  if (! _simName.empty())
    ssFileName << "_" << _simName;
  ssFileName << "_tissues-" << _nbTissues
	     << "_scenario-" << _scenario
	     << ".bed.gz";
  if (! doesFileExist (ssFileName.str()))
  {
    if (verbose > 0)
      cout << "write file " << ssFileName.str() << endl;
    gzFile fileStream;
    openFile (ssFileName.str(), fileStream, "wb");
    
    stringstream ssTxt;
    for (size_t v = 0; v < _nbSnps; ++v)
    {
      ssTxt.str("");
      ssTxt << "chr1"
	    << "\t" << v * 50
	    << "\t" << v * 50 + 1
	    << "\tsnp" << (v+1)
	    << endl;
      gzwriteLine (fileStream, ssTxt.str(), ssFileName.str(), (v+1));
    }
    
    closeFile (ssFileName.str(), fileStream);
  }
  else
  {
    if (verbose > 0)
      cout << "read file " << ssFileName.str() << endl;
    gzFile fileStream;
    openFile (ssFileName.str(), fileStream, "rb");
    
    string line;
    size_t nbLines = 0;
    while (getline (fileStream, line))
      ++nbLines;
    if (nbLines != _nbGenes)
    {
      cerr << "ERROR: different nb of genes in file " << ssFileName.str()
	   << endl;
      exit (1);
    }
    
    closeFile (ssFileName.str(), fileStream);
  }
}

void
Params::simulAndWriteOrLoadGenotypes (
  const gsl_rng * rng,
  const int & verbose)
{
  stringstream ssFileName;
  ssFileName << "genotypes";
  if (! _simName.empty())
    ssFileName << "_" << _simName;
  ssFileName << "_tissues-" << _nbTissues
	     << "_scenario-" << _scenario
	     << ".txt.gz";
  gzFile fileStream;
  
  if (! doesFileExist (ssFileName.str()))
  {
    if (verbose > 0)
      cout << "simulate genotypes and write in file " << ssFileName.str()
	   << endl;
    vector<vector<size_t> > vvGenotypes (_nbSnps,
					 vector<size_t> (_nbInds, 0));
    openFile (ssFileName.str(), fileStream, "wb");
    
    stringstream ssTxt;
    ssTxt << "ind1";
    for (size_t n = 1; n < _nbInds; ++n)
      ssTxt << " ind" << (n+1);
    ssTxt << endl;
    gzwriteLine (fileStream, ssTxt.str(), ssFileName.str(), 1);
    
    gsl_rng_set (rng, _mainSeed);
    for (size_t v = 0; v < _nbSnps; ++v)
    {
      ssTxt.str("");
      ssTxt << "snp" << (v+1);
      for (size_t n = 0; n < _nbInds; ++n)
      {
	vvGenotypes[v][n] = gsl_ran_binomial (rng, _maf, 2); // assume HWE
	ssTxt << " " << vvGenotypes[v][n];
      }
      ssTxt << endl;
      gzwriteLine (fileStream, ssTxt.str(), ssFileName.str(), (v+1));
    }
    _vvGenotypes = vvGenotypes;
    
    stringstream ssGenoListName;
    ssGenoListName << "list_genotypes";
    if (! _simName.empty())
      ssGenoListName << "_" << _simName;
    ssGenoListName << "_tissues-" << _nbTissues
		   << "_scenario-" << _scenario
		   << ".txt";
    ofstream genoListStream;
    openFile (ssGenoListName.str(), genoListStream);
    vector<string> vGenoFileNames;
    for (size_t s = 0; s < _nbTissues; ++s)
      genoListStream << "tissue" << (s+1) << " " << getCurrentDirectory ()
		     << "/" << ssFileName.str() << endl;
    closeFile (ssGenoListName.str(), genoListStream);
  }
  
  else
  {
    if (verbose > 0)
      cout << "load genotypes from file " << ssFileName.str() << endl;
    vector<size_t> vGenotypes;
    openFile (ssFileName.str(), fileStream, "rb");
    
    string line;
    if (! getline (fileStream, line))
    {
      cerr << "ERROR: problem with the header of file " << ssFileName.str()
	   << endl;
      exit (1);
    }
    if (line.empty())
    {
      cerr << "ERROR: file " << ssFileName.str() << " is empty" << endl;
      exit (1);
    }
    vector<string> tokens;
    split (line, " \t", tokens);
    if (tokens.size() != _nbInds)
    {
      cerr << "ERROR: different nb of individuals in file " << ssFileName.str()
	   << endl;
      exit (1);
    }
    vGenotypes.resize (tokens.size(), 0);
    while (getline (fileStream, line))
    {
      split (line, " \t", tokens);
      for (size_t n = 1; n < tokens.size(); ++n)
	vGenotypes[(n-1)] = atol (tokens[n].c_str());
      _vvGenotypes.push_back (vGenotypes);
    }
    if (! gzeof (fileStream))
    {
      cerr << "ERROR: can't read file " << ssFileName.str()
	   << " successfully up to the end" << endl;
      exit (1);
    }
    if (_vvGenotypes.size() != _nbGenes)
    {
      cerr << "ERROR: different nb of genes in file " << ssFileName.str()
	   << endl;
      exit (1);
    }
  }
  
  closeFile (ssFileName.str(), fileStream);
}

void
Params::simulAndWritePhenoForGivenRep (
  const string & truthFileName,
  gzFile & truthFileStream,
  stringstream & ssTruth,
  const vector<string> & vPhenoFileNames,
  const gsl_rng * rng)
{
  vector<gzFile> vPhenoFileStreams;
  for (size_t s = 0; s < _nbTissues; ++s)
  {
    gzFile phenoFileStream;
    vPhenoFileStreams.push_back (phenoFileStream);
    openFile (vPhenoFileNames[s], vPhenoFileStreams[s], "wb");
  }
  
  stringstream ssPheno;
  ssPheno << "ind1";
  for (size_t n = 1; n < _nbInds; ++n)
    ssPheno << " ind" << (n+1);
  ssPheno << endl;
  for (size_t s = 0; s < _nbTissues; ++s)
    gzwriteLine (vPhenoFileStreams[s], ssPheno.str(),
		 vPhenoFileNames[s], 1);
  
  double oma2_plus_phi2 = _pve / ((1 - _pve) * 2 * _maf * (1 - _maf));
  double phi2 = _het * oma2_plus_phi2; // variance of the b's
  double oma2 = oma2_plus_phi2 - phi2; // variance of \bar{b}
  double bbar; // average effect size across tissues
  vector<double> vBs (_nbTissues, 0); // vector of standardized effect sizes
  
  for (size_t g = 0; g < _nbGenes; ++g)
  {
    if (g < _nbGenesPerConfig) // null config
    {
      ssTruth.str("");
      ssTruth << "gene" << (g+1) << " snp" << (g+1)
	      << " 0 NA NA "; // pve phi omega
      for (size_t i = 0; i < _nbTissues; ++i)
	ssTruth << "0"; // config
      ssTruth << " 0"; // bbar
      for (size_t i = 0; i < vBs.size(); ++i)
	ssTruth << " " << vBs[i];
      ssTruth << endl;
      gzwriteLine (truthFileStream, ssTruth.str(), truthFileName, (g+1));
      
      for (size_t s = 0; s < _nbTissues; ++s)
      {
	ssPheno.str("");
	ssPheno << "gene" << (g+1);
	for (size_t n = 0; n < _nbInds; ++n)
	  ssPheno << " " << gsl_ran_gaussian_ziggurat (rng, _vSigmas[s]);
	ssPheno << endl;
	gzwriteLine (vPhenoFileStreams[s], ssPheno.str(),
		     vPhenoFileNames[s], (g+1));
      }
    }
    else // simulate all other configs
    {
      ssTruth.str("");
      ssTruth << "gene" << (g+1) << " snp" << (g+1) << " " << _pve
	      << " " << sqrt(phi2) << " " << sqrt(oma2) << " ";
      size_t idxConfig = 0;
      while (g >= _nbGenesPerConfig + (idxConfig + 1) * _nbGenesPerConfig)
	++idxConfig;
      for (size_t s = 0; s < _nbTissues; ++s)
	ssTruth << _vvConfigs[idxConfig][s];
      bbar = gsl_ran_gaussian_ziggurat (rng, sqrt(oma2));
      ssTruth << " " << bbar;
      for (size_t s = 0; s < _nbTissues; ++s)
      {
	if (_vvConfigs[idxConfig][s] == 0)
	  vBs[s] = 0;
	else
	  vBs[s] = bbar + gsl_ran_gaussian_ziggurat (rng, sqrt(phi2));
	ssTruth << " " << vBs[s];
      }
      ssTruth << endl;
      gzwriteLine (truthFileStream, ssTruth.str(), truthFileName, (g+1));
      
      for (size_t s = 0; s < _nbTissues; ++s)
      {
	ssPheno.str("");
	ssPheno << "gene" << (g+1);
	for (size_t n = 0; n < _nbInds; ++n)
	  ssPheno << " " << vBs[s] * _vSigmas[s] * _vvGenotypes[g][n]
	    + gsl_ran_gaussian_ziggurat (rng, _vSigmas[s]);
	ssPheno << endl;
	gzwriteLine (vPhenoFileStreams[s], ssPheno.str(),
		     vPhenoFileNames[s], (g+1));
      }
    }
  }
  
  for (size_t s = 0; s < _nbTissues; ++s)
    closeFile (vPhenoFileNames[s], vPhenoFileStreams[s]);
}

void
Params::simulAndWritePhenotypes (
  const gsl_rng * rng,
  const int & verbose)
{
  if (verbose > 0)
    cout << "simulate and write phenotypes ..." << endl;
  
  createDirectory (_nameOfProjectDir);
  changeDirectory (_nameOfProjectDir);
  
  gsl_rng_set (rng, _mainSeed);
  for(size_t r = 0; r < _nbReplicates; ++r)
    _vRepSeeds[r] = gsl_rng_uniform_int (rng, 1000000);
  
  for (size_t r = 0; r < _nbReplicates; ++r)
  {
    if (verbose > 0)
      cout << "replicate " << (r+1) << "/" << _nbReplicates << endl << flush;
    stringstream ssDirName;
    ssDirName << "rep-" << setfill('0')
	      << setw ((int) ceil (log10 (_nbReplicates+1))) << (r+1);
    if (doesFileExist (ssDirName.str()))
    {
      cerr << "ERROR: directory " << ssDirName.str() << " already exists"
	   << endl;
      exit (1);
    }
    createDirectory (ssDirName.str());
    changeDirectory (ssDirName.str());
    
    string truthFileName = "truth.txt.gz";
    gzFile truthFileStream;
    openFile (truthFileName, truthFileStream, "wb");
    
    stringstream ssTruth;
    ssTruth << "gene snp pve phi oma config bbar";
    for (size_t s = 0; s < _nbTissues; ++s)
      ssTruth << " bg.s" << (s+1);
    ssTruth << endl;
    gzwriteLine (truthFileStream, ssTruth.str(), truthFileName, 1);
    
    gsl_rng_set (rng, _vRepSeeds[r]);
    
    stringstream ssPhenoListName;
    ssPhenoListName << "list_phenotypes.txt";
    ofstream phenoListStream;
    openFile (ssPhenoListName.str(), phenoListStream);
    vector<string> vPhenoFileNames;
    for (size_t s = 0; s < _nbTissues; ++s)
    {
      stringstream ssPhenoFileName;
      ssPhenoFileName << "phenotypes_tissue" << (s+1) << ".txt.gz";
      phenoListStream << "tissue" << (s+1) << " " << ssPhenoFileName.str()
		      << endl;
      vPhenoFileNames.push_back (ssPhenoFileName.str());
      if (_choiceErrVar != 1)
	Params::drawSigma (s, rng);
    }
    closeFile (ssPhenoListName.str(), phenoListStream);
    simulAndWritePhenoForGivenRep (truthFileName, truthFileStream, ssTruth,
				   vPhenoFileNames, rng);
    
    closeFile (truthFileName, truthFileStream);
    
    Params::write (r);
    
    changeDirectory ("..");
  }
  
  changeDirectory ("..");
}

vector<vector<float> >
getGrid (
  const string & gridType,
  const bool & noHet
  )
{
  static const float oma2_plus_phi2[] = {0.1*0.1, 0.2*0.2, 0.4*0.4, 0.8*0.8, 1.6*1.6}; // average effect size
  vector<float> oma2_over_oma2_plus_phi2; // homogeneity
  if (gridType.compare("general") == 0)
  {
    oma2_over_oma2_plus_phi2.push_back (0);
    oma2_over_oma2_plus_phi2.push_back (0.25);
    oma2_over_oma2_plus_phi2.push_back (0.5);
    oma2_over_oma2_plus_phi2.push_back (0.75);
    oma2_over_oma2_plus_phi2.push_back (1);
  }
  else
  {
    if (noHet)
      oma2_over_oma2_plus_phi2.push_back (1);
    else
    {
      oma2_over_oma2_plus_phi2.push_back (0.75);
      oma2_over_oma2_plus_phi2.push_back (1);
    }
  }
  
  vector<vector<float> > grid (5 * oma2_over_oma2_plus_phi2.size(),
			       vector<float> (2, 0));
  size_t k = 0;
  for (size_t i = 0; i < 5; ++i)
  {
    for (size_t j = 0; j < oma2_over_oma2_plus_phi2.size(); ++j)
    {
      grid[k][0] = oma2_plus_phi2[i] * (1 - oma2_over_oma2_plus_phi2[j]);
      grid[k][1] = oma2_plus_phi2[i] * oma2_over_oma2_plus_phi2[j];
      ++k;
    }
  }
  return (grid);
}

void
writeGrids (const int & verbose)
{
  gzFile fileStream;
  stringstream ssTxt;
  string fileName;
  
  fileName = "grid_phi2_oma2_general.txt.gz";
  if (! doesFileExist (fileName))
  {
    if (verbose > 0)
      cout << "write large grid" << endl;
    vector<vector<float> > grid = getGrid ("general", false);
    openFile (fileName, fileStream, "wb");
    for (size_t i = 0; i < grid.size(); ++i)
    {
      ssTxt.str("");
      // ssTxt << setprecision(7) << tmp;
      ssTxt << grid[i][0] << "\t" << grid[i][1] << endl;
      gzwriteLine (fileStream, ssTxt.str(), fileName, i+1);
    }
    closeFile (fileName, fileStream);
  }
  
  fileName = "grid_phi2_oma2_with-configs.txt.gz";
  if (! doesFileExist (fileName))
  {
    if (verbose > 0)
      cout << "write small grid" << endl;
    vector<vector<float> > grid = getGrid ("with-configs", false);
    openFile (fileName, fileStream, "wb");
    for (size_t i = 0; i < grid.size(); ++i)
    {
      ssTxt.str("");
      ssTxt << grid[i][0] << "\t" << grid[i][1] << endl;
      gzwriteLine (fileStream, ssTxt.str(), fileName, i+1);
    }
    closeFile (fileName, fileStream);
  }
  
  fileName = "grid_phi2_oma2_with-configs_no-het.txt.gz";
  if (! doesFileExist (fileName))
  {
    if (verbose > 0)
      cout << "write very small grid" << endl;
    vector<vector<float> > grid = getGrid ("with-configs", true);
    openFile (fileName, fileStream, "wb");
    for (size_t i = 0; i < grid.size(); ++i)
    {
      ssTxt.str("");
      ssTxt << grid[i][0] << "\t" << grid[i][1] << endl;
      gzwriteLine (fileStream, ssTxt.str(), fileName, i+1);
    }
    closeFile (fileName, fileStream);
  }
}

void
run (
  const string & simName,
  const size_t & nbTissues,
  const size_t & scenario,
  const size_t & nbReplicates,
  const size_t & mainSeed,
  const size_t & choiceErrVar,
  const size_t & choiceHet,
  const double & pve,
  const double & maf, 
  const size_t & nbInds,
  const size_t & nbGenesPerConfig,
  const int & verbose)
{
  gsl_rng * rng = NULL;
  gsl_rng_env_setup();
  rng = gsl_rng_alloc (gsl_rng_default);
  if (rng == NULL)
  {
    cerr << "ERROR: can't allocate memory for the RNG" << endl;
    exit (1);
  }
  
  writeGrids (verbose);
  
  Params iParams (simName, nbTissues, scenario, nbReplicates, mainSeed,
		  choiceErrVar, choiceHet, pve, maf, nbInds,
		  nbGenesPerConfig);
  iParams.simulAndWriteGeneCoords (verbose);
  iParams.simulAndWriteSnpCoords (verbose);
  iParams.simulAndWriteOrLoadGenotypes (rng, verbose);
  iParams.simulAndWritePhenotypes (rng, verbose);
  
  gsl_rng_free (rng);
}

int main (int argc, char ** argv)
{
  int verbose = 1;
  string simName;
  size_t nbTissues = 1, scenario = 1, nbReplicates = 1, mainSeed = 1859,
    choiceErrVar = 1, choiceHet = 1, nbInds = 100, nbGenesPerConfig = 1000;
  double pve = 0.2, maf = 0.3;
  
  parseArgs (argc, argv, simName, nbTissues, scenario, nbReplicates,
	     mainSeed, choiceErrVar, choiceHet, pve, maf, nbInds,
	     nbGenesPerConfig, verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << basename(argv[0])
	 << " " << getDateTime (startRawTime) << endl
	 << "version " << VERSION << " compiled " << __DATE__
	 << " " << __TIME__ << endl
	 << "cmd-line: " << getCmdLine (argc, argv) << endl
	 << "cwd: " << getCurrentDirectory() << endl;
    cout << flush;
  }
  
  run (simName, nbTissues, scenario, nbReplicates, mainSeed, choiceErrVar,
       choiceHet, pve, maf, nbInds, nbGenesPerConfig, verbose);
  
  if (verbose > 0)
  {
    time (&endRawTime);
    cout << "END " << basename(argv[0])
	 << " " << getDateTime (endRawTime) << endl
	 << "elapsed -> " << getElapsedTime(startRawTime, endRawTime) << endl
	 << "max.mem -> " << getMaxMemUsedByProcess2Str () << endl;
  }
  
  return EXIT_SUCCESS;
}
