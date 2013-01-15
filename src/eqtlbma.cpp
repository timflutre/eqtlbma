/** \file eqtlbma.cpp
 *
 *  `eqtlbma' performs eQTL mapping in multiple subgroups via a Bayesian model.
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
 *  g++ -Wall -Wextra -g -DEQTLBMA_MAIN -DLIB_MVLR MVLR.cc utils.cpp eqtlbma.cpp -lgsl -lgslcblas -lz -o eqtlbma
 */

#include <cmath>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <getopt.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>
#include <sstream>
#include <numeric>
using namespace std;

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_gamma.h>

#include "utils.h"

#ifdef LIB_MVLR
#include "MVLR.h"
#endif

/** \brief Display the help on stdout.
 */
void help (char ** argv)
{
  cout << "`" << argv[0] << "'"
       << " performs eQTL mapping in multiple subgroups via a Bayesian model." << endl
       << endl
       << "Usage: " << argv[0] << " [OPTIONS] ..." << endl
       << endl
       << "Options:" << endl
       << "  -h, --help\tdisplay the help and exit" << endl
       << "  -V, --version\toutput version information and exit" << endl
       << "  -v, --verbose\tverbosity level (0/default=1/2/3)" << endl
       << "  -g, --geno\tfile with absolute paths to genotype files" << endl
       << "\t\ttwo columns: subgroup identifier<space/tab>path to file" << endl
       << "\t\tadd '#' at the beginning of a line to comment it" << endl
       << "\t\tsubgroup file: can be in three formats (VCF/IMPUTE/custom)" << endl
       << "\t\tVCF: see specifications on 1kG website" << endl
       << "\t\tIMPUTE: row 1 is a header chr<del>name<del>coord<del>a1<del>a2" << endl
       << "\t\t followed by >sample1_a1a1<del>sample1_a1a2<del>sample1_a2a2<del>..." << endl
       << "\t\tcustom: row 1 for sample names, column 1 for SNP names, genotypes as allele dose" << endl
       << "      --scoord\tfile with the SNP coordinates (compulsory if custom genotype format)" << endl
       << "\t\tshould be in the BED format (delimiter: tab)" << endl
       << "\t\tSNPs in the genotype files without coordinate are skipped" << endl
       << "  -p, --pheno\tfile with absolute paths to phenotype files" << endl
       << "\t\ttwo columns: subgroup identifier<space/tab>path to file" << endl
       << "\t\tcan be a single line (single subgroup)" << endl
       << "\t\tadd '#' at the beginning of a line to comment it" << endl
       << "\t\tsubgroup file: row 1 for sample names, column 1 for feature names" << endl
       << "\t\tsubgroups can have different features" << endl
       << "\t\tall features should be in the --fcoord file" << endl
       << "      --fcoord\tfile with the feature coordinates" << endl
       << "\t\tshould be in the BED format (delimiter: tab)" << endl
       << "\t\tfeatures in the phenotype files without coordinate are skipped" << endl
       << "      --anchor\tfeature boundary(ies) for the cis region" << endl
       << "\t\tdefault=FSS, can also be FSS+FES" << endl
       << "      --cis\tlength of half of the cis region (in bp)" << endl
       << "\t\tapart from the anchor(s), default=100000" << endl
       << "  -o, --out\tprefix for the output files" << endl
       << "\t\tall output files are gzipped and have a header line" << endl
       << "      --step\tstep of the analysis to perform" << endl
       << "\t\t1: only separate analysis of each subgroup, without permutation" << endl
       << "\t\t2: only separate analysis of each subgroup, with permutation" << endl
       << "\t\t3: both separate and joint analysis, without permutation" << endl
       << "\t\t4: both separate and joint analysis, with permutation for joint only" << endl
       << "\t\t5: both separate and joint analysis, with permutation for both" << endl
       << "      --outraw\twrite the output file with all raw ABFs" << endl
       << "\t\tby default, only the ABFs averaged over the grids are written" << endl
       << "\t\twriting all raw ABFs can be too much if the number of subgroups is large" << endl
       << "      --qnorm\tquantile-normalize the phenotypes" << endl
       << "      --maf\tminimum minor allele frequency (default=0.0)" << endl
       << "      --covar\tfile with absolute paths to covariate files" << endl
       << "\t\ttwo columns: subgroup identifier<space/tab>path to file" << endl
       << "\t\tcan be a single line (single subgroup)" << endl
       << "\t\tadd '#' at the beginning of a line to comment it" << endl
       << "\t\tsubgroup file: row 1 is a header sample<space/tab>covariate1 ..." << endl
       << "\t\tall sample names should be in the respective genotype and phenotype files" << endl
       << "\t\tthe covariates should be numbers, no missing value is allowed" << endl
       << "\t\tsubgroups can have different covariates" << endl
       << "\t\tthe order of rows is not important" << endl
       << "      --gridL\tfile with a 'large' grid for prior variances in standardized effect sizes" << endl
       << "\t\tfirst column is phi^2 and second column is omega^2, no header" << endl
       << "\t\tthis grid is used with model 1 ('general alternative') trying to capture" << endl
       << "\t\t all sorts of types of heterogeneity" << endl
       << "\t\trequired with --step 3, 4 or 5, whatever --bfs is" << endl
       << "      --gridS\tfile with a 'small' grid of values for phi^2 and omega^2" << endl
       << "\t\tsame format as --gridL, but usually a single value for phi^2 (eg. 0.01)" << endl
       << "\t\trequired with --step 3, 4 or 5 if --bfs is subset or all" << endl
       << "      --bfs\twhich Bayes Factors to compute for the joint analysis" << endl
       << "\t\teach BF for a given configuration is the average of the BFs over one of the grids, with equal weights" << endl
       << "\t\tonly the Laplace-approximated BF from Wen and Stephens is implemented" << endl
       << "\t\t'gen' (default): general way to capture any level of heterogeneity" << endl
       << "\t\t correspond to the consistent configuration with the large grid" << endl
       << "\t\t fixed-effect and maximum-heterogeneity BFs are also calculated" << endl
       << "\t\t'sin': compute also the BF for each singleton (subgroup-specific configuration)" << endl
       << "\t\t they use the small grid (BFgen-sin is also reported)" << endl
       << "\t\t'all': compute also the BFs for all configurations (costly if many subgroups)" << endl
       << "\t\t all BFs use the small grid" << endl
#ifdef LIB_MVLR
       << "      --mvlr\tuse the multivariate version of the ABF" << endl
       << "\t\tallows for correlation between samples in the errors" << endl
       << "\t\tespecially useful when subgroups come from same individuals" << endl
       << "\t\tcaution, no separate analysis can be performed with this option" << endl
       << "      --fitsig\tparam used when estimating the variance of the errors with --mvlr" << endl
       << "\t\tdefault=0 (null model), but can be between 0 and 1" << endl
#endif
       << "      --nperm\tnumber of permutations" << endl
       << "\t\tdefault=0, recommended=10000" << endl
       << "      --seed\tseed for the two random number generators" << endl
       << "\t\tone for the permutations, another for the trick" << endl
       << "\t\tby default, both are initialized via microseconds from epoch" << endl
       << "\t\tthe RNGs are re-seeded before each subgroup and before the joint analysis" << endl
       << "\t\tthis, along with --trick 2, allows for proper comparison of separate and joint analyzes" << endl
       << "      --trick\tapply trick to speed-up permutations" << endl
       << "\t\tstop after the tenth permutation for which the test statistic" << endl
       << "\t\t is better than or equal to the true value, and sample from" << endl
       << "\t\t a uniform between 11/(nbPermsSoFar+2) and 11/(nbPermsSoFar+1)" << endl
       << "\t\tif '1', the permutations really stops" << endl
       << "\t\tif '2', all permutations are done but the test statistics are not computed" << endl
       << "\t\tallows to compare different test statistics on the same permutations" << endl
       << "      --permsep\twhich permutation procedure for the separate analysis" << endl
       << "\t\t1 (default): use the minimum P-value over SNPs and subgroups as a test statistic" << endl
       << "\t\t2: use the minimum P-value over SNPs but in each subgroup separately (breaks correlations)" << endl
       << "      --pbf\twhich BF to use as the test statistic for the joint-analysis permutations" << endl
       << "\t\t'gen' (default): general BF (see --bfs above)" << endl
       << "\t\t'sin': average only over the singletons" << endl
       << "\t\t'gen-sin': 0.5 BFgen + 0.5 BFsin" << endl
       << "\t\t'all': average over all configurations" << endl
       << "      --maxbf\tuse the maximum ABF over SNPs as test statistic for permutations" << endl
       << "\t\totherwise the average ABF over SNPs is used (more Bayesian)" << endl
       << "  -f, --ftr\tfile with a list of features to analyze" << endl
       << "\t\tone feature name per line" << endl
       << "\t\tallows to easily parallelize a whole analyzis" << endl
       << "  -s, --snp\tfile with a list of SNPs to analyze" << endl
       << "\t\tone SNP name per line" << endl
       << "      --sbgrp\tidentifier of the subgroup to analyze" << endl
       << "\t\tuseful for quick analysis and debugging" << endl
       << endl;
}

/** \brief Display version and license information on stdout.
 */
void version (char ** argv)
{
  cout << argv[0] << " " << __DATE__ << " " << __TIME__ << endl
       << endl
       << "Copyright (C) 2012-2013 Timothee Flutre and Xiaoquan Wen." << endl
       << "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>" << endl
       << "This is free software; see the source for copying conditions.  There is NO" << endl
       << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl
       << endl
       << "Written by Timothee Flutre and Xiaoquan Wen." << endl;
}

/** \brief Parse the command-line arguments and check the values of the 
 *  compulsory ones.
 */
void
parseArgs (
  int argc,
  char ** argv,
  string & genoPathsFile,
  string & snpCoordsFile,
  string & phenoPathsFile,
  string & ftrCoordsFile,
  string & anchor,
  size_t & lenCis,
  string & outPrefix,
  bool & outRaw,
  int & whichStep,
  bool & needQnorm,
  float & minMaf,
  string & covarFile,
  string & largeGridFile,
  string & smallGridFile,
  string & whichBfs,
  bool & mvlr,
  float & propFitSigma,
  size_t & nbPerms,
  size_t & seed,
  int & trick,
  int & whichPermSep,
  string & whichPermBf,
  bool & useMaxBfOverSnps,
  string & ftrsToKeepFile,
  string & snpsToKeepFile,
  string & sbgrpToKeep,
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
	{"geno", required_argument, 0, 'g'},
	{"scoord", required_argument, 0, 0},
	{"pheno", required_argument, 0, 'p'},
	{"fcoord", required_argument, 0, 0},
	{"anchor", required_argument, 0, 0},
	{"cis", required_argument, 0, 0},
	{"out", required_argument, 0, 'o'},
	{"outraw", no_argument, 0, 0},
	{"step", required_argument, 0, 0},
	{"qnorm", no_argument, 0, 0},
	{"maf", required_argument, 0, 0},
	{"covar", required_argument, 0, 0},
	{"gridL", required_argument, 0, 0},
	{"gridS", required_argument, 0, 0},
	{"bfs", required_argument, 0, 0},
	{"mvlr", no_argument, 0, 0},
	{"fitsig", required_argument, 0, 0},
	{"nperm", required_argument, 0, 0},
	{"seed", required_argument, 0, 0},
	{"trick", required_argument, 0, 0},
	{"permsep", required_argument, 0, 0},
	{"pbf", required_argument, 0, 0},
	{"maxbf", no_argument, 0, 0},
	{"ftr", required_argument, 0, 'f'},
	{"snp", required_argument, 0, 's'},
	{"sbgrp", required_argument, 0, 0},
	{0, 0, 0, 0}
      };
    int option_index = 0;
    c = getopt_long (argc, argv, "hVv:g:p:o:f:s:",
		     long_options, &option_index);
    if (c == -1)
      break;
    switch (c)
    {
    case 0:
      if (long_options[option_index].flag != 0)
	break;
      if (strcmp(long_options[option_index].name, "scoord") == 0)
      {
	snpCoordsFile = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "fcoord") == 0)
      {
	ftrCoordsFile = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "anchor") == 0)
      {
	anchor = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "cis") == 0)
      {
	lenCis = atol (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "outraw") == 0)
      {
	outRaw = true;
	break;
      }
      if (strcmp(long_options[option_index].name, "step") == 0)
      {
	whichStep = atoi (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "qnorm") == 0)
      {
	needQnorm = true;
	break;
      }
      if (strcmp(long_options[option_index].name, "maf") == 0)
      {
	minMaf = atof (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "covar") == 0)
      {
	covarFile = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "gridL") == 0)
      {
	largeGridFile = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "gridS") == 0)
      {
	smallGridFile = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "bfs") == 0)
      {
	whichBfs = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "mvlr") == 0)
      {
	mvlr = true;
	break;
      }
      if (strcmp(long_options[option_index].name, "fitsig") == 0)
      {
	propFitSigma = atof (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "nperm") == 0)
      {
	nbPerms = atol (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "seed") == 0)
      {
	seed = atol (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "trick") == 0)
      {
	trick = atoi (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "permsep") == 0)
      {
	whichPermSep = atoi (optarg);
	break;
      }
      if (strcmp(long_options[option_index].name, "pbf") == 0)
      {
	whichPermBf = optarg;
	break;
      }
      if (strcmp(long_options[option_index].name, "maxbf") == 0)
      {
	useMaxBfOverSnps = true;
	break;
      }
      if (strcmp(long_options[option_index].name, "sbgrp") == 0)
      {
	sbgrpToKeep = optarg;
	break;
      }
    case 'h':
      help (argv);
      exit (0);
    case 'V':
      version (argv);
      exit (0);
    case 'v':
      verbose = atoi (optarg);
      break;
    case 'g':
      genoPathsFile = optarg;
      break;
    case 'p':
      phenoPathsFile = optarg;
      break;
    case 'o':
      outPrefix = optarg;
      break;
    case 'f':
      ftrsToKeepFile = optarg;
      break;
    case 's':
      snpsToKeepFile = optarg;
      break;
    case '?':
      printf ("\n"); help (argv);
      abort ();
    default:
      printf ("\n"); help (argv);
      abort ();
    }
  }
  if (genoPathsFile.empty())
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: missing compulsory option -g\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (genoPathsFile))
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: can't file '%s'\n\n", genoPathsFile.c_str());
    help (argv);
    exit (1);
  }
  if (! snpCoordsFile.empty() && ! doesFileExist (snpCoordsFile))
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: can't file '%s'\n\n", snpCoordsFile.c_str());
    help (argv);
    exit (1);
  }
  if (phenoPathsFile.empty())
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: missing compulsory option -p\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (phenoPathsFile))
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: can't find '%s'\n\n", phenoPathsFile.c_str());
    help (argv);
    exit (1);
  }
  if (ftrCoordsFile.empty())
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: missing compulsory option --fcoord\n\n");
    help (argv);
    exit (1);
  }
  if (! doesFileExist (ftrCoordsFile))
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: can't find '%s'\n\n", ftrCoordsFile.c_str());
    help (argv);
    exit (1);
  }
  if (anchor.empty())
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: SNPs in trans not yet implemented, see --anchor and --cis\n\n");
    help (argv);
    exit (1);
  }
  if (outPrefix.empty())
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: missing compulsory option -o\n\n");
    help (argv);
    exit (1);
  }
  // check writing files is possible
  stringstream ssOutFile;
  ssOutFile << outPrefix << "_test.txt.gz";
  gzFile outStream;
  openFile (ssOutFile.str(), outStream, "wb");
  closeFile (ssOutFile.str(), outStream);
  remove (ssOutFile.str().c_str());
  if (whichStep != 1 && whichStep != 2 && whichStep != 3 && whichStep != 4
      && whichStep != 5)
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: --step should be 1, 2, 3, 4 or 5\n\n");
    help (argv);
    exit (1);
  }
  if ((whichStep == 3 || whichStep == 4 || whichStep == 5)
      && largeGridFile.empty())
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: missing compulsory option --gridL with --step 3, 4 or 5\n\n");
    help (argv);
    exit (1);
  }
  if (! largeGridFile.empty() && ! doesFileExist (largeGridFile))
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: can't find '%s'\n\n", largeGridFile.c_str());
    help (argv);
    exit (1);
  }
  if ((whichStep == 3 || whichStep == 4 || whichStep == 5)
      && (whichBfs.compare("sin") == 0 || whichBfs.compare("all") == 0)
      && smallGridFile.empty())
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: missing compulsory option --gridS with --step 3, 4 or 5 and --bfs sin or all\n\n");
    help (argv);
    exit (1);
  }
  if (! smallGridFile.empty() && ! doesFileExist (smallGridFile))
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: can't find '%s'\n\n", smallGridFile.c_str());
    help (argv);
    exit (1);
  }
  if (whichBfs.compare("gen") != 0 && whichBfs.compare("sin") != 0
      && whichBfs.compare("all") != 0)
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: --bf should be 'gen', 'sin' or 'all'\n\n");
    help (argv);
    exit (1);
  }
  if (mvlr && whichStep == 5)
    fprintf (stderr, "WARNING: separate analysis won't be performed with --mvlr\n\n");
  if ((whichStep == 2 || whichStep == 4 || whichStep == 5) && nbPerms == 0)
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: --step %i but nbPerms = 0, see --nperm\n\n", whichStep);
    help (argv);
    exit (1);
  }
  if (trick != 0 && trick != 1 && trick != 2)
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: --trick should be 0, 1 or 2\n\n");
    help (argv);
    exit (1);
  }
  if ((whichStep == 2 || whichStep == 5) &&
      whichPermSep != 1 && whichPermSep != 2)
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: --permsep should be 1 or 2\n\n");
    help (argv);
    exit (1);
  }
  if ((whichStep == 4 || whichStep == 5) && whichBfs.compare("gen") == 0
      && whichPermBf.compare("gen") != 0)
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: if --bfs gen, then --pbf should be 'gen' too\n\n");
    help (argv);
    exit (1);
  }
  if ((whichStep == 4 || whichStep == 5) && whichBfs.compare("sin") == 0
      && whichPermBf.compare("all") == 0)
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: if --bfs sin, then --pbf should be 'gen', 'sin' or 'gen-sin'\n\n");
    help (argv);
    exit (1);
  }
  if (! ftrsToKeepFile.empty() && ! doesFileExist (ftrsToKeepFile))
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: can't find '%s'\n\n", ftrsToKeepFile.c_str());
    help (argv);
    exit (1);
  }
  if (! snpsToKeepFile.empty() && ! doesFileExist (snpsToKeepFile))
  {
    printCmdLine (cerr, argc, argv);
    fprintf (stderr, "ERROR: can't find '%s'\n\n", snpsToKeepFile.c_str());
    help (argv);
    exit (1);
  }
  if (seed == string::npos)
    seed = getSeed ();
}

class Grid
{
public:
  vector<double> phi2s;
  vector<double> oma2s;
  vector<double> phi2s_fix;
  vector<double> oma2s_fix;
  vector<double> phi2s_maxh;
  vector<double> oma2s_maxh;
  Grid (const string & gridFile,
	const bool & makeFixMaxh,
	const int & verbose);
  size_t size (void) const { return phi2s.size(); }
};

Grid::Grid (
  const string & gridFile,
  const bool & makeFixMaxh,
  const int & verbose)
{
  if (! gridFile.empty())
  {
    if (verbose > 0)
      cout << "load grid in " << gridFile << " ..." << endl << flush;
    
    gzFile gridStream;
    vector<string> tokens;
    string line;
    openFile (gridFile, gridStream, "rb");
    while (getline (gridStream, line))
    {
      split (line, " \t", tokens);
      if (tokens.size() != 2)
      {
	cerr << "ERROR: format of file " << gridFile
	     << " should be phi2<space/tab>oma2" << endl;
	exit (1);
      }
      phi2s.push_back (atof (tokens[0].c_str()));
      oma2s.push_back (atof (tokens[1].c_str()));
      if (makeFixMaxh)
      {
	phi2s_fix.push_back (0.0);
	oma2s_fix.push_back (atof (tokens[0].c_str())
			     + atof (tokens[1].c_str()));
	phi2s_maxh.push_back (atof (tokens[0].c_str())
			      + atof (tokens[1].c_str()));
	oma2s_maxh.push_back (0.0);
      }
    }
    if (! gzeof (gridStream))
    {
      cerr << "ERROR: can't read successfully file "
	   << gridFile
	   << " up to the end" << endl;
      exit (1);
    }
    closeFile (gridFile, gridStream);
    
    if (verbose > 0)
      cout << "grid size: " << phi2s.size() << endl;
  }
}

/** \brief Fit a simple linear regression model.
 *  \note y_i = mu + g_i * beta + e_i with e_i ~ N(0,sigma^2)
 *  \note missing values should have been already filtered out
 *  \note g and y should have the same size, strictly larger than 2
 */
void
fitSimpleLinearRegression (
  const vector<double> & g,
  const vector<double> & y,
  double & pve,
  double & sigmahat,
  double & betahat,
  double & sebetahat,
  double & pval)
{
  if (g.size() == y.size() && g.size() > 2)
  {
    size_t i = 0, n = g.size();
    double ym = 0, gm = 0, yty = 0, gtg = 0, gty = 0;
    for(i=0; i<n; ++i)
    {
      ym += y[i];
      gm += g[i];
      yty += y[i] * y[i];
      gtg += g[i] * g[i];
      gty += g[i] * y[i];
    }
    ym /= n;
    gm /= n;
    double vg = gtg - n * gm * gm;  // variance of g
#ifdef DEBUG
    printf ("n=%zu ym=%f gm=%f yty=%f gtg=%f gty=%f vg=%f\n",
	    n, ym, gm, yty, gtg, gty, vg);
#endif
    if(vg > 1e-8)
    {
      betahat = (gty - n * gm * ym) / vg;
      double rss1 = yty - 1/vg * (n*ym*(gtg*ym - gm*gty) - gty*(n*gm*ym - gty));
      if (fabs(betahat) > 1e-8)
	sigmahat = sqrt(rss1 / (n-2));
      else  // case where y is not variable enough among samples
	sigmahat = sqrt((yty - n * ym * ym) / (n-2));  // sqrt(rss0/(n-2))
      sebetahat = sigmahat / sqrt(gtg - n*gm*gm);
      double muhat = (ym*gtg - gm*gty) / (gtg - n*gm*gm);
      double mss = 0;
      for(i=0; i<n; ++i)
	mss += pow(muhat + betahat * g[i] - ym, 2);
      pval = gsl_cdf_fdist_Q (mss/pow(sigmahat,2), 1, n-2);
      pve = mss / (mss + rss1);
    }
    else
    {
      betahat = 0;
      sebetahat = numeric_limits<double>::infinity();
      sigmahat = sqrt((yty - n * ym * ym) / (n-2));  // sqrt(rss0/(n-2))
      pval = 1; // or should be drawn from U([0,1])?
      pve = 0;
    }
  }
}

/** \brief Fit a multiple linear regression model.
 *  \note y_i = beta_0 + beta_1 * g_i + beta_2 * c_i + ... + e_i with e_i ~ N(0,sigma^2)
 *  \note missing values should have been already filtered out
 */
void
fitMultipleLinearRegression (
  const vector<double> & vG,
  const vector<double> & vY,
  const vector<vector<double> > & vvCovars,
  double & pve,
  double & sigmahat,
  vector<vector<double> > & vvResPredictors)
{
  size_t N = vY.size(), P = 2 + vvCovars.size();
  gsl_vector * Y = gsl_vector_alloc (N), * Bhat = gsl_vector_alloc (P);
  gsl_matrix * X = gsl_matrix_alloc (N, P),
    * covBhat = gsl_matrix_alloc (P, P);
  gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (N, P);
  
  for (size_t n = 0; n < N; ++n)
  {
    gsl_vector_set (Y, n, vY[n]);
    gsl_matrix_set (X, n, 0, 1.0);
    gsl_matrix_set (X, n, 1, vG[n]);
    for (size_t p = 2; p < P; ++p)
      gsl_matrix_set (X, n, p, vvCovars[(p-2)][n]);
  }
  
  double rss;
  size_t rank;
  gsl_multifit_linear_svd (X, Y, GSL_DBL_EPSILON, &rank, Bhat, covBhat,
			   &rss, work);
  pve = 1 - rss / gsl_stats_tss (Y->data, Y->stride,
				 Y->size);
  sigmahat = sqrt (rss / (N-rank));
  
  vvResPredictors[0][0] = gsl_vector_get (Bhat, 1);
  vvResPredictors[0][1] = sqrt (gsl_matrix_get (covBhat, 1, 1));
  vvResPredictors[0][2] = 2 * gsl_cdf_tdist_Q (
    fabs(vvResPredictors[0][0] / vvResPredictors[0][1]), N - rank);
  
  for (size_t p = 2; p < P; ++p)
  {
    vvResPredictors[(p-1)][0] = gsl_vector_get (Bhat, p);
    vvResPredictors[(p-1)][1] = sqrt (gsl_matrix_get (covBhat, p, p));
    vvResPredictors[(p-1)][2] = 2 * gsl_cdf_tdist_Q (
      fabs(vvResPredictors[(p-1)][0] / vvResPredictors[(p-1)][1]), N - rank);
  }
  
  gsl_vector_free (Y);
  gsl_vector_free (Bhat);
  gsl_matrix_free (X);
  gsl_matrix_free (covBhat);
  gsl_multifit_linear_free (work);
}

/** \brief Return the approximate Bayes Factor from Wen and Stephens (2011)
 *  \note univariate version of the ABF
 *  \param vvStdSstats vector of vectors, one per subgroup, with betahat as
 *  first entry, sebetahat as second and the frequentist test stat t as third
 */
double
getAbfFromStdSumStats (
  const vector<size_t> & vNs,
  const vector<vector<double> > & vvStdSstats,
  const double & phi2,
  const double & oma2)
{
  double l10AbfAll = 0, bhat = 0, varbhat = 0, t = 0, bbarhat_num = 0,
    bbarhat_denom = 0, varbbarhat = 0;
  vector<double> l10AbfsSingleSbgrp;
  
  for (size_t s = 0; s < vNs.size(); ++s)
  {
    if (vNs[s] > 2)
    {
      bhat = vvStdSstats[s][0];
      varbhat = pow (vvStdSstats[s][1], 2);
      t = vvStdSstats[s][2];
      double lABF_single;
      if (fabs(t) < 1e-8)
      {
	lABF_single = 0;
      }
      else
      {
	bbarhat_num += bhat / (varbhat + phi2);
	bbarhat_denom += 1 / (varbhat + phi2);
	varbbarhat += 1 / (varbhat + phi2);
	lABF_single = 0.5 * log10(varbhat)
	  - 0.5 * log10(varbhat + phi2)
	  + (0.5 * pow(t,2) * phi2 / (varbhat + phi2)) / log(10);
      }
      l10AbfsSingleSbgrp.push_back (lABF_single);
#ifdef DEBUG
      printf ("l10AbfsSingleSbgrp[%ld]=%e\n", s+1, l10AbfsSingleSbgrp[s]);
#endif
    }
  }
  
  double bbarhat = (bbarhat_denom != 0) ?
    bbarhat_num / bbarhat_denom
    : 0;
  varbbarhat = (varbbarhat != 0) ?
    1 / varbbarhat
    : numeric_limits<double>::infinity();
  if (bbarhat != 0 && varbbarhat < numeric_limits<double>::infinity())
  {
    double T2 = pow(bbarhat, 2.0) / varbbarhat;
    double lABF_bbar = (T2 != 0) ?
      0.5 * log10(varbbarhat) - 0.5 * log10(varbbarhat + oma2)
      + (0.5 * T2 * oma2 / (varbbarhat + oma2)) / log(10)
      : 0;
#ifdef DEBUG
    printf ("bbarhat=%e varbbarhat=%e T2=%e lABF_bbar=%e\n",
	    bbarhat, varbbarhat, T2, lABF_bbar);
#endif
    l10AbfAll = lABF_bbar;
    for (size_t i = 0; i < l10AbfsSingleSbgrp.size(); ++i)
      l10AbfAll += l10AbfsSingleSbgrp[i];
  }
  else // MAF very close, or equal, to 0
    l10AbfAll = 0;
  
  return l10AbfAll;
}

struct Snp
{
  string name; // eg. rs7263524
  string chr; // eg. chr21
  size_t coord; // 1-based coordinate
  vector<vector<double> > vvGenos; // genotypes of samples per subgroup
  vector<vector<bool> > vvIsNa; // missing values per subgroup
  vector<double> vMafs; // minor allele frequencies per subgroup
};

struct ResFtrSnp
{
  string snp; // name of the SNP
  vector<size_t> vNs; // sample sizes per subgroup
  vector<double> vPves; // proportions of variance explained (R2)
  vector<double> vSigmahats; // MLEs of the error variance per subgroup
  
  // summary stats for predictors (genotype + covariates) per subgroup
  // key: name of predictors ; values: summary stats
  // 0:betahat ; 1:sebetahat ; 2:betaPval
  vector<map<string, vector<double> > > vMapPredictors;
  
  // raw ABFs
  // keys: 'gen', 'gen-fix', 'gen-maxh', '1-2-3', '1', '2', etc
  map<string, vector<double> > mUnweightedAbfs;
  
  // averaged ABFs (over a grid for all, over configurations for some)
  // keys: same as above
  map<string, double> mWeightedAbfs;
};

struct Ftr
{
  string name; // eg. ENSG00000182816
  string chr; // eg. chr21
  size_t start; // 1-based coordinate
  size_t end; // idem
  vector<vector<double> > vvPhenos; // phenotypes of samples per subgroup
  vector<vector<bool> > vvIsNa; // missing values per subgroup
  
  vector<Snp*> vPtCisSnps; // pointers to SNP(s) in cis
  vector<ResFtrSnp> vResFtrSnps; // for each SNP in vPtCisSnps
  
  // test statistic: min P-value over SNPs in each subgroup
  vector<double> vPermPvalsSep; // permutation P-values in each subgroup
  vector<size_t> vNbPermsSoFar; // nb of permutations in each subgroup
  vector<double> vMinTruePvals; // min true P-value over SNPs in each subgroup
  
  // test statistic: min P-value over SNPs across subgroups
  double sepPermPval;
  size_t nbPermsSoFarSep;
  double minTruePval;
  
  // test statistic: ABF (max or avg over SNPs)
  double jointPermPval;
  size_t nbPermsSoFarJoint;
  double maxL10TrueAbf;
  double avgL10TrueAbf;
};

void
Snp_init (
  Snp & iSnp,
  const string & name,
  const size_t & nbSubgroups)
{
  iSnp.name = name;
  iSnp.chr.clear();
  iSnp.coord = string::npos;
  iSnp.vvGenos.resize (nbSubgroups);
  iSnp.vvIsNa.resize (nbSubgroups);
  iSnp.vMafs = (vector<double> (nbSubgroups,
				numeric_limits<double>::quiet_NaN()));
  for (size_t s = 0; s < nbSubgroups; ++s)
  {
    iSnp.vvGenos[s] = (vector<double> ());
    iSnp.vvIsNa[s] = (vector<bool> ());
  }
}

// assume both features are on the same chromosome
bool Snp_compByCoord (
  const Snp* pt_iSnp1,
  const Snp* pt_iSnp2)
{
  bool res = false;
  if (pt_iSnp1->coord < pt_iSnp2->coord)
    res = true;
  return res;
}

int Snp_isInCis (
  const Snp & iSnp,
  const size_t & ftrStart,
  const size_t & ftrEnd,
  const string & anchor,
  const size_t & lenCis)
{
  int res = -1;
  if (anchor.compare("FSS+FES") == 0)
  {
    if (((ftrStart >= lenCis &&
	 iSnp.coord >= ftrStart - lenCis) ||
	(ftrStart < lenCis)) &&
	iSnp.coord <= ftrEnd + lenCis)
      res = 0;
    else if (iSnp.coord > ftrEnd + lenCis)
      res = 1;
  }
  else if (anchor.compare("FSS") == 0)
  {
    if (((ftrStart >= lenCis &&
	 iSnp.coord >= ftrStart - lenCis) ||
	(ftrStart < lenCis)) &&
	iSnp.coord <= ftrStart + lenCis)
      res = 0;
    else if (iSnp.coord > ftrStart + lenCis)
      res = 1;
  }
  return res;
}

void
ResFtrSnp_init (
  ResFtrSnp & iResFtrSnp,
  const string & snpName,
  const size_t & nbSubgroups)
{
  iResFtrSnp.snp = snpName;
  iResFtrSnp.vNs.assign (nbSubgroups, 0);
  iResFtrSnp.vMapPredictors.assign (nbSubgroups,
  				    map<string, vector<double> > ());
  for (size_t s = 0; s < nbSubgroups; ++s)
    iResFtrSnp.vMapPredictors[s]["genotype"] =
      vector<double>  (3, numeric_limits<double>::quiet_NaN());
  iResFtrSnp.vSigmahats.assign (nbSubgroups,
  				numeric_limits<double>::quiet_NaN());
  iResFtrSnp.vPves.assign (nbSubgroups,
			   numeric_limits<double>::quiet_NaN());
}

/** \brief Compute the summary statistics for a given feature-snp pair
 *  in a given subgroup by simple linear regression
 */
void
ResFtrSnp_getSstatsOneSbgrp (
  ResFtrSnp & iResFtrSnp,
  const Ftr & iFtr,
  const Snp & iSnp,
  const size_t & s,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const bool & needQnorm)
{
  vector<double> y, g;
  size_t idxPheno, idxGeno;
  for (size_t i = 0; i < vvSampleIdxPhenos[s].size(); ++i)
  {
    idxPheno = vvSampleIdxPhenos[s][i];
    idxGeno = vvSampleIdxGenos[s][i];
    if (idxPheno != string::npos
	&& idxGeno != string::npos
	&& ! iFtr.vvIsNa[s][idxPheno]
	&& ! iSnp.vvIsNa[s][idxGeno])
    {
      y.push_back (iFtr.vvPhenos[s][idxPheno]);
      g.push_back (iSnp.vvGenos[s][idxGeno]);
    }
  }
  
  iResFtrSnp.vNs[s] = y.size();
  
  if (iResFtrSnp.vNs[s] > 2)
  {
    if (needQnorm)
      qqnorm (&y[0], y.size());
    
    fitSimpleLinearRegression (g, y, iResFtrSnp.vPves[s],
			       iResFtrSnp.vSigmahats[s],
			       iResFtrSnp.vMapPredictors[s]["genotype"][0],
			       iResFtrSnp.vMapPredictors[s]["genotype"][1],
			       iResFtrSnp.vMapPredictors[s]["genotype"][2]);
  }
}

/** \brief Compute the summary statistics for a given feature-snp pair
 *  in a given subgroup by multiple linear regression
 */
void
ResFtrSnp_getSstatsOneSbgrp (
  ResFtrSnp & iResFtrSnp,
  const Ftr & iFtr,
  const Snp & iSnp,
  const size_t & s,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const bool & needQnorm,
  const vector<map<string, vector<double> > > & vSbgrp2Covars)
{
  vector<double> vY, vG;
  size_t idxPheno, idxGeno, j;
  vector<vector<double> > vvCovars (vSbgrp2Covars[s].size(),
				    vector<double> ());
  for (size_t i = 0; i < vvSampleIdxPhenos[s].size(); ++i)
  {
    idxPheno = vvSampleIdxPhenos[s][i];
    idxGeno = vvSampleIdxGenos[s][i];
    if (idxPheno != string::npos
	&& idxGeno != string::npos
	&& ! iFtr.vvIsNa[s][idxPheno]
	&& ! iSnp.vvIsNa[s][idxGeno])
    {
      vY.push_back (iFtr.vvPhenos[s][idxPheno]);
      vG.push_back (iSnp.vvGenos[s][idxGeno]);
      j = 0;
      for (map<string, vector<double> >::const_iterator it =
	     vSbgrp2Covars[s].begin(); it != vSbgrp2Covars[s].end(); ++it)
      {
	vvCovars[j].push_back (it->second[i]);
	++j;
      }
    }
  }
  
  iResFtrSnp.vNs[s] = vY.size();
  
  if (iResFtrSnp.vNs[s] > 2)
  {
    if (needQnorm)
      qqnorm (&vY[0], vY.size());
    
    vector<vector<double> > vvResPredictors (
      1 + vvCovars.size(),
      vector<double> (3, numeric_limits<double>::quiet_NaN()));
    fitMultipleLinearRegression (vG, vY, vvCovars,
				 iResFtrSnp.vPves[s],
				 iResFtrSnp.vSigmahats[s],
				 vvResPredictors);
    iResFtrSnp.vMapPredictors[s]["genotype"] = vvResPredictors[0];
    j = 1;
    for (map<string, vector<double> >::const_iterator it =
	   vSbgrp2Covars[s].begin(); it != vSbgrp2Covars[s].end(); ++it)
    {
      iResFtrSnp.vMapPredictors[s][it->first] = vvResPredictors[j];
	++j;
    }
  }
}

/** \brief Compute the summary statistics for a given feature-snp pair
 *  in a given subgroup by simple linear regression after permuting
 *  the phenotypes according to the given permutation
 */
void
ResFtrSnp_getSstatsPermOneSbgrp (
  ResFtrSnp & iResFtrSnp,
  const Ftr & iFtr,
  const Snp & iSnp,
  const size_t & s,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const bool & needQnorm,
  const gsl_permutation * perm)
{
  vector<double> y, g;
  size_t idxPheno, idxGeno, p;
  for (size_t i = 0; i < vvSampleIdxPhenos[s].size(); ++i)
  {
    p = gsl_permutation_get (perm, i);
    idxPheno = vvSampleIdxPhenos[s][p];
    idxGeno = vvSampleIdxGenos[s][i];
    if (idxPheno != string::npos
	&& idxGeno != string::npos
	&& ! iFtr.vvIsNa[s][idxPheno]
	&& ! iSnp.vvIsNa[s][idxGeno])
    {
      y.push_back (iFtr.vvPhenos[s][idxPheno]);
      g.push_back (iSnp.vvGenos[s][idxGeno]);
    }
  }
  
  iResFtrSnp.vNs[s] = y.size();
  
  if (iResFtrSnp.vNs[s] > 2)
  {
    if (needQnorm)
      qqnorm (&y[0], y.size());
    
    fitSimpleLinearRegression (g, y, iResFtrSnp.vPves[s],
			       iResFtrSnp.vSigmahats[s],
			       iResFtrSnp.vMapPredictors[s]["genotype"][0],
			       iResFtrSnp.vMapPredictors[s]["genotype"][1],
			       iResFtrSnp.vMapPredictors[s]["genotype"][2]);
  }
}

/** \brief Compute the summary statistics for a given feature-snp pair
 *  in a given subgroup by multiple linear regression after permuting
 *  the phenotypes according to the given permutation
 */
void
ResFtrSnp_getSstatsPermOneSbgrp (
  ResFtrSnp & iResFtrSnp,
  const Ftr & iFtr,
  const Snp & iSnp,
  const size_t & s,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const bool & needQnorm,
  const vector<map<string, vector<double> > > & vSbgrp2Covars,
  const gsl_permutation * perm)
{
  vector<double> vY, vG;
  size_t idxPheno, idxGeno, j, p;
  vector<vector<double> > vvCovars (vSbgrp2Covars[s].size(),
				    vector<double> ());
  for (size_t i = 0; i < vvSampleIdxPhenos[s].size(); ++i)
  {
    p = gsl_permutation_get (perm, i);
    idxPheno = vvSampleIdxPhenos[s][p];
    idxGeno = vvSampleIdxGenos[s][i];
    if (idxPheno != string::npos
	&& idxGeno != string::npos
	&& ! iFtr.vvIsNa[s][idxPheno]
	&& ! iSnp.vvIsNa[s][idxGeno])
    {
      vY.push_back (iFtr.vvPhenos[s][idxPheno]);
      vG.push_back (iSnp.vvGenos[s][idxGeno]);
      j = 0;
      for (map<string, vector<double> >::const_iterator it =
	     vSbgrp2Covars[s].begin(); it != vSbgrp2Covars[s].end(); ++it)
      {
	vvCovars[j].push_back (it->second[i]);
	++j;
      }
    }
  }
  
  iResFtrSnp.vNs[s] = vY.size();
  
  if (iResFtrSnp.vNs[s] > 2)
  {
    if (needQnorm)
      qqnorm (&vY[0], vY.size());
    
    vector<vector<double> > vvResPredictors (
      1 + vvCovars.size(),
      vector<double> (3, numeric_limits<double>::quiet_NaN()));
    fitMultipleLinearRegression (vG, vY, vvCovars,
				 iResFtrSnp.vPves[s],
				 iResFtrSnp.vSigmahats[s],
				 vvResPredictors);
    iResFtrSnp.vMapPredictors[s]["genotype"] = vvResPredictors[0];
    j = 1;
    for (map<string, vector<double> >::const_iterator it =
	   vSbgrp2Covars[s].begin(); it != vSbgrp2Covars[s].end(); ++it)
    {
      iResFtrSnp.vMapPredictors[s][it->first] = vvResPredictors[j];
	++j;
    }
  }
}

/** \brief Standardize the summary statistics and correct for small 
 *  sample size if requested
 */
void
ResFtrSnp_getStdStatsAndCorrSmallSampleSize (
  const ResFtrSnp & iResFtrSnp,
  vector<vector<double> > & vvStdSstats)
{
  for (size_t s = 0; s < iResFtrSnp.vNs.size(); ++s)
  {
    if (iResFtrSnp.vNs[s] > 2)
    {
      double n = iResFtrSnp.vNs[s],
	bhat = iResFtrSnp.vMapPredictors[s].find("genotype")->second[0] /
	iResFtrSnp.vSigmahats[s],
	sebhat = iResFtrSnp.vMapPredictors[s].find("genotype")->second[1] /
	iResFtrSnp.vSigmahats[s],
	t = bhat / sebhat;
      // apply quantile-quantile transformation
      t = gsl_cdf_gaussian_Pinv (gsl_cdf_tdist_P (-fabs(bhat/sebhat),
						  n-2), 1.0);
      vector<double> vStdSstats;
      if (fabs(t) > 1e-8)
      {
	double sigmahat = fabs (iResFtrSnp.vMapPredictors[s].find("genotype")
				->second[0]) /
	  (fabs (t) * sebhat);
	bhat = iResFtrSnp.vMapPredictors[s].find("genotype")->second[0] /
	  sigmahat;
	sebhat = bhat / t;
      }
      else
      {
	bhat = 0;
	sebhat = numeric_limits<double>::infinity();
      }
      vStdSstats.push_back (bhat);
      vStdSstats.push_back (sebhat);
      vStdSstats.push_back (t);
      vvStdSstats.push_back (vStdSstats);
    }
    else // iResFtrSnp.vNs[s] <= 2
      vvStdSstats.push_back (vector<double> (
			       3, numeric_limits<double>::quiet_NaN()));
  }
}

/** \brief Calculate the ABF for the consistent configuration,
 *  averaged over the large grid for the general model, the fixed-effect
 *  model and the maximum-heterogeneity model
 *  \note the configuration is represented by 'gen', 'gen-fix'
 *  and 'gen-maxh'
 */
void
ResFtrSnp_calcAbfsConstLargeGrid (
  ResFtrSnp & iResFtrSnp,
  const Grid & iGridL,
  const vector<vector<double> > & vvStdSstats)
{
  vector<double> vL10AbfsConst (iGridL.size(), 0.0),
    vL10AbfsConstFix (iGridL.size(), 0.0),
    vL10AbfsConstMaxh (iGridL.size(), 0.0);
  
  for (size_t gridIdx = 0; gridIdx < iGridL.size(); ++gridIdx)
  {
    vL10AbfsConst[gridIdx] = (getAbfFromStdSumStats (
			       iResFtrSnp.vNs,
			       vvStdSstats,
			       iGridL.phi2s[gridIdx],
			       iGridL.oma2s[gridIdx]));
    vL10AbfsConstFix[gridIdx] = (getAbfFromStdSumStats (
				   iResFtrSnp.vNs,
				   vvStdSstats,
				   0,
				   iGridL.phi2s[gridIdx]
				   + iGridL.oma2s[gridIdx]));
    vL10AbfsConstMaxh[gridIdx] = (getAbfFromStdSumStats (
				    iResFtrSnp.vNs,
				    vvStdSstats,
				    iGridL.phi2s[gridIdx]
				    + iGridL.oma2s[gridIdx],
				    0));
  }
  iResFtrSnp.mUnweightedAbfs.insert (make_pair ("gen",
						vL10AbfsConst));
  iResFtrSnp.mUnweightedAbfs.insert (make_pair ("gen-fix",
						vL10AbfsConstFix));
  iResFtrSnp.mUnweightedAbfs.insert (make_pair ("gen-maxh",
						vL10AbfsConstMaxh));
  
  iResFtrSnp.mWeightedAbfs.insert (make_pair ("gen",
					      log10_weighted_sum (
						&(vL10AbfsConst[0]),
						vL10AbfsConst.size())));
  iResFtrSnp.mWeightedAbfs.insert (make_pair ("gen-fix",
					      log10_weighted_sum (
						&(vL10AbfsConstFix[0]),
						vL10AbfsConstFix.size())));
  iResFtrSnp.mWeightedAbfs.insert (make_pair ("gen-maxh",
					      log10_weighted_sum (
						&(vL10AbfsConstMaxh[0]),
						vL10AbfsConstMaxh.size())));
}

/** \brief Calculate the ABF for each singleton, averaged over the small grid
 *  \note the configuration is represented by the index of the given subgroup
 *  (eg. '1' if in the first subgroup)
 */
void
ResFtrSnp_calcAbfsSingletons (
  ResFtrSnp & iResFtrSnp,
  const Grid & iGridS,
  const vector<vector<double> > & vvStdSstatsAll)
{
  stringstream ssConfig;
  vector<size_t> vNs (iResFtrSnp.vNs.size(), 0);
  vector< vector<double> > vvStdSstatsSingletons;
  vector<double> vL10Abfs;
  
  for (size_t s = 0; s < iResFtrSnp.vNs.size(); ++s)
  {
    ssConfig.str("");
    ssConfig << (s+1);
    
    if (iResFtrSnp.vNs[s] > 2)
    {
      vvStdSstatsSingletons.clear ();
      for (size_t i = 0; i < iResFtrSnp.vNs.size(); ++i)
      {
	if (s == i)
	{
	  vNs[i] = iResFtrSnp.vNs[i];
	  vvStdSstatsSingletons.push_back (vvStdSstatsAll[i]);
	}
	else
	{
	  vNs[i] = 0;
	  vvStdSstatsSingletons.push_back (vector<double> (3, 0));
	}
      }
      
      vL10Abfs.assign (iGridS.size(), 0);
      for (size_t gridIdx = 0; gridIdx < iGridS.size(); ++gridIdx)
	vL10Abfs[gridIdx] = getAbfFromStdSumStats (vNs,
						   vvStdSstatsSingletons,
						   iGridS.phi2s[gridIdx],
						   iGridS.oma2s[gridIdx]);
      iResFtrSnp.mUnweightedAbfs.insert (make_pair (ssConfig.str(),
						    vL10Abfs));
      iResFtrSnp.mWeightedAbfs.insert (make_pair(ssConfig.str(),
						 log10_weighted_sum (
						   &(vL10Abfs[0]),
						   vL10Abfs.size())));
    }
    else // iResFtrSnp.vNs[s] <= 2
    {
      iResFtrSnp.mUnweightedAbfs.insert (
	make_pair (ssConfig.str(),
		   vector<double> (
		     iGridS.size(),
		     numeric_limits<double>::quiet_NaN())));
      iResFtrSnp.mWeightedAbfs.insert (
	make_pair (ssConfig.str(),
		   numeric_limits<double>::quiet_NaN()));
    }
  }
}

/** \brief Calculate the two averaged ABFs 'sin' and 'gen-sin'
 */
void
ResFtrSnp_calcAbfsAvgSinAndGenSin (
  ResFtrSnp & iResFtrSnp)
{
  stringstream ssConfig;
  vector<double> vL10Abfs, vWeights;
  for (size_t s = 0; s < iResFtrSnp.vNs.size(); ++s)
  {
    ssConfig.str("");
    ssConfig << (s+1);
    if (! isNan (iResFtrSnp.mWeightedAbfs[ssConfig.str()]))
    {
      vL10Abfs.push_back (iResFtrSnp.mWeightedAbfs[ssConfig.str()]);
      vWeights.push_back ((1.0 / 2.0) *
			  (1.0 / gsl_sf_choose (iResFtrSnp.vNs.size(), 1)));
    }
  }
  
  if (vL10Abfs.size() > 0)
  {
    iResFtrSnp.mWeightedAbfs.insert (
      make_pair ("sin",
		 log10_weighted_sum (&(vL10Abfs[0]), vL10Abfs.size())));
    
    vL10Abfs.push_back (iResFtrSnp.mWeightedAbfs["gen"]);
    vWeights.push_back (1.0 / 2.0);
    iResFtrSnp.mWeightedAbfs.insert (
      make_pair ("gen-sin",
		 log10_weighted_sum (&(vL10Abfs[0]), &(vWeights[0]),
				     vL10Abfs.size())));
  }
  else
  {
    iResFtrSnp.mWeightedAbfs.insert (
      make_pair ("sin", numeric_limits<double>::quiet_NaN()));
    iResFtrSnp.mWeightedAbfs.insert (
      make_pair ("gen-sin", numeric_limits<double>::quiet_NaN()));
  }
}

/** \brief The configuration is represented only by the index/indices
 *  of the subgroup(s) in which the eQTL is active (eg. '2' or '6-7-11')
 */
static void
prepareConfig (
  stringstream & ssConfig,
  vector<bool> & vIsEqtlInConfig,
  const gsl_combination * comb)
{
  ssConfig.str("");
  vIsEqtlInConfig.clear();
  
  vIsEqtlInConfig.assign (comb->n, false);
  
  ssConfig << gsl_combination_get (comb, 0) + 1;
  vIsEqtlInConfig[gsl_combination_get (comb, 0)] = true;
  if (comb->k > 1)
  {
    for (size_t i = 1; i < comb->k; ++i)
    {
      ssConfig << "-" << gsl_combination_get (comb, i) + 1;
      vIsEqtlInConfig[gsl_combination_get (comb, i)] = true;
    }
  }
}

/** \brief Calculate the ABF for all configurations, averaged over
 *  the small grid
 */
void
ResFtrSnp_calcAbfsAllConfigs (
  ResFtrSnp & iResFtrSnp,
  const Grid & iGridS,
  const vector<vector<double> > & vvStdSstatsAll)
{
  gsl_combination * comb;
  stringstream ssConfig;
  vector<bool> vIsEqtlInConfig; // T,T,F if S=3 and config="1-2"
  vector<size_t> vNs (iResFtrSnp.vNs.size(), 0);
  vector<vector<double> > vvStdSstatsSubset;
  vector<double> vL10Abfs;
  
  for (size_t k = 1; k <= iResFtrSnp.vNs.size(); ++k)
  {
    comb = gsl_combination_calloc (iResFtrSnp.vNs.size(), k);
    if (comb == NULL)
    {
      cerr << "ERROR: can't allocate memory for the combination" << endl;
      exit (1);
    }
    while (true)
    {
      prepareConfig (ssConfig, vIsEqtlInConfig, comb);
      vvStdSstatsSubset.clear();
      for (size_t s = 0; s < iResFtrSnp.vNs.size(); ++s)
      {
	if (iResFtrSnp.vNs[s] > 2 && vIsEqtlInConfig[s])
	{
	  vNs[s] = iResFtrSnp.vNs[s];
	  vvStdSstatsSubset.push_back (vvStdSstatsAll[s]);
	}
	else
	{
	  vNs[s] = 0;
	  vvStdSstatsSubset.push_back (vector<double> (3, 0));
	}
      }
      if (accumulate (vNs.begin(), vNs.end(), 0) > 2)
      {
	vL10Abfs.assign (iGridS.size(), 0);
	for (size_t gridIdx = 0; gridIdx < iGridS.size(); ++gridIdx)
	  vL10Abfs[gridIdx] = getAbfFromStdSumStats (vNs,
						     vvStdSstatsSubset,
						     iGridS.phi2s[gridIdx],
						     iGridS.oma2s[gridIdx]);
	iResFtrSnp.mUnweightedAbfs.insert (make_pair (ssConfig.str(),
						      vL10Abfs));
	iResFtrSnp.mWeightedAbfs.insert (make_pair(ssConfig.str(),
						   log10_weighted_sum (
						     &(vL10Abfs[0]),
						     vL10Abfs.size())));
      }
      else
      {
	iResFtrSnp.mUnweightedAbfs.insert (
	  make_pair(ssConfig.str(), vector<double> (
		      iGridS.size(), numeric_limits<double>::quiet_NaN())));
	iResFtrSnp.mWeightedAbfs.insert (
	  make_pair(ssConfig.str(), numeric_limits<double>::quiet_NaN()));
      }
      if (gsl_combination_next (comb) != GSL_SUCCESS)
	break;
    }
    gsl_combination_free (comb);
  }
}

/** \brief Calculate the ABF averaged over all configurations,
 *  each being itself averaged over the small grid
 */
void
ResFtrSnp_calcAbfAvgAll (
  ResFtrSnp & iResFtrSnp)
{
  gsl_combination * comb;
  stringstream ssConfig;
  vector<bool> vIsEqtlInConfig; // T,T,F if S=3 and config="1-2"
  vector<double> vL10Abfs, vWeights;
  
  for (size_t k = 1; k <= iResFtrSnp.vNs.size(); ++k)
  {
    comb = gsl_combination_calloc (iResFtrSnp.vNs.size(), k);
    if (comb == NULL)
    {
      cerr << "ERROR: can't allocate memory for the combination" << endl;
      exit (1);
    }
    while (true)
    {
      prepareConfig (ssConfig, vIsEqtlInConfig, comb);
      if (! isNan (iResFtrSnp.mWeightedAbfs[ssConfig.str()]))
      {
	vL10Abfs.push_back (iResFtrSnp.mWeightedAbfs[ssConfig.str()]);
	vWeights.push_back ((1.0 / (double) iResFtrSnp.vNs.size()) *
			    (1.0 / gsl_sf_choose (iResFtrSnp.vNs.size(), k)));
      }
      if (gsl_combination_next (comb) != GSL_SUCCESS)
	break;
    }
    gsl_combination_free (comb);
  }
  
  if (vL10Abfs.size() > 0)
    iResFtrSnp.mWeightedAbfs.insert (
      make_pair ("all",
		 log10_weighted_sum (&(vL10Abfs[0]), &(vWeights[0]),
				     vL10Abfs.size())));
  else
    iResFtrSnp.mWeightedAbfs.insert (
      make_pair ("all", numeric_limits<double>::quiet_NaN()));
}

/** \brief Calculate the ABFs for a given feature-SNP pair
 *  \note depending on whichBfs, this function calculates only the ABF
 *  for the consistent configuration, or that one as well as the ABFs
 *  for each subgroup-specific configuration, or the ABFs for all
 *  \note whatever the configuration, the ABFs for each grid value are 
 *  averaged with equal weights
 *  \note the ABF for the consistent configuration is always calculated
 *  on the large grid, and sometimes also on the small one
 */
void
ResFtrSnp_calcAbfs (
  ResFtrSnp & iResFtrSnp,
  const string & whichBfs,
  const Grid & iGridL,
  const Grid & iGridS)
{
  vector<vector<double> > vvStdSstats;
  ResFtrSnp_getStdStatsAndCorrSmallSampleSize (iResFtrSnp,
					       vvStdSstats);
  ResFtrSnp_calcAbfsConstLargeGrid (iResFtrSnp, iGridL, vvStdSstats);
  if (whichBfs.compare("sin") == 0)
  {
    ResFtrSnp_calcAbfsSingletons (iResFtrSnp, iGridS, vvStdSstats);
    ResFtrSnp_calcAbfsAvgSinAndGenSin (iResFtrSnp);
  }
  else if (whichBfs.compare("all") == 0)
  {
    ResFtrSnp_calcAbfsAllConfigs (iResFtrSnp, iGridS, vvStdSstats);
    ResFtrSnp_calcAbfsAvgSinAndGenSin (iResFtrSnp);
    ResFtrSnp_calcAbfAvgAll (iResFtrSnp);
  }
}

void
ResFtrSnp_prepareDataForMvlr (
  ResFtrSnp & iResFtrSnp,
  const Ftr & iFtr,
  const Snp & iSnp,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const vector<map<string, vector<double> > > & vSbgrp2Covars,
  const bool & needQnorm,
  vector<vector<double> > & Y,
  vector<vector<double> > & Xg,
  vector<vector<double> > & Xc)
{
  // don't initialize Y yet as the feature may not be present in all subgroups
  Xg.assign (1, vector<double> ());
  Xc.assign (vSbgrp2Covars[0].size(), vector<double> ());
  for (size_t s = 0; s < iFtr.vvPhenos.size(); ++s)
  {
    vector<double> tmpY;
    size_t idxPheno, idxGeno;
    for (size_t i = 0; i < vvSampleIdxPhenos[s].size(); ++i)
    {
      idxPheno = vvSampleIdxPhenos[s][i];
      idxGeno = vvSampleIdxGenos[s][i];
      if (idxPheno != string::npos
	  && idxGeno != string::npos
	  && ! iFtr.vvIsNa[s][idxPheno]
	  && ! iSnp.vvIsNa[s][idxGeno])
      {
	tmpY.push_back (iFtr.vvPhenos[s][idxPheno]);
	Xg[0].push_back (iSnp.vvGenos[s][idxGeno]);
	for (map<string, vector<double> >::const_iterator it =
	       vSbgrp2Covars[s].begin(); it != vSbgrp2Covars[s].end(); ++it)
	  Xc[s].push_back (it->second[i]);
      }
    }
    if (tmpY.size() != 0)
    {
      iResFtrSnp.vNs[s] = tmpY.size();
      if (needQnorm)
	qqnorm (&tmpY[0], tmpY.size());
      Y.push_back (tmpY);
    }
  }
}

void
ResFtrSnp_prepareDataForMvlrPerm (
  ResFtrSnp & iResFtrSnp,
  const Ftr & iFtr,
  const Snp & iSnp,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const vector<map<string, vector<double> > > & vSbgrp2Covars,
  const bool & needQnorm,
  const gsl_permutation * perm,
  vector<vector<double> > & Y,
  vector<vector<double> > & Xg,
  vector<vector<double> > & Xc)
{
  // don't initialize Y yet as the feature may not be present in all subgroups
  Xg.assign (1, vector<double> ());
  Xc.assign (vSbgrp2Covars[0].size(), vector<double> ());
  for (size_t s = 0; s < iFtr.vvPhenos.size(); ++s)
  {
    vector<double> tmpY;
    size_t idxPheno, idxGeno, p;
    for (size_t i = 0; i < vvSampleIdxPhenos[s].size(); ++i)
    {
      p = gsl_permutation_get (perm, i);
      idxPheno = vvSampleIdxPhenos[s][p];
      idxGeno = vvSampleIdxGenos[s][i];
      if (idxPheno != string::npos
	  && idxGeno != string::npos
	  && ! iFtr.vvIsNa[s][idxPheno]
	  && ! iSnp.vvIsNa[s][idxGeno])
      {
	tmpY.push_back (iFtr.vvPhenos[s][idxPheno]);
	Xg[0].push_back (iSnp.vvGenos[s][idxGeno]);
	for (map<string, vector<double> >::const_iterator it =
	       vSbgrp2Covars[s].begin(); it != vSbgrp2Covars[s].end(); ++it)
	  Xc[s].push_back (it->second[i]);
      }
    }
    if (tmpY.size() != 0)
    {
      iResFtrSnp.vNs[s] = tmpY.size();
      if (needQnorm)
	qqnorm (&tmpY[0], tmpY.size());
      Y.push_back (tmpY);
    }
  }
}

void
ResFtrSnp_calcAbfsConstLargeGridMvlr (
  ResFtrSnp & iResFtrSnp,
  const Grid & iGridL,
  const float & propFitSigma,
  vector<vector<double> > & Y,
  vector<vector<double> > & Xg,
  vector<vector<double> > & Xc)
{
#ifdef LIB_MVLR
  MVLR iMvlr;
  iMvlr.set_sigma_option (propFitSigma);
  iMvlr.init (Y, Xg, Xc);
  vector<vector<int> > vvGamma (1, vector<int> (iResFtrSnp.vNs.size(), 1));
  iMvlr.set_effect_vec (iGridL.phi2s, iGridL.oma2s);
  vector<double> vL10AbfsConst = iMvlr.compute_log10_ABF_vec (vvGamma);
  iMvlr.set_effect_vec (iGridL.phi2s_fix, iGridL.oma2s_fix);
  vector<double> vL10AbfsConstFix = iMvlr.compute_log10_ABF_vec (vvGamma);
  iMvlr.set_effect_vec (iGridL.phi2s_maxh, iGridL.oma2s_maxh);
  vector<double> vL10AbfsConstMaxh = iMvlr.compute_log10_ABF_vec (vvGamma);
#endif
  iResFtrSnp.mUnweightedAbfs.insert (make_pair ("gen",
						vL10AbfsConst));
  iResFtrSnp.mUnweightedAbfs.insert (make_pair ("gen-fix",
						vL10AbfsConstFix));
  iResFtrSnp.mUnweightedAbfs.insert (make_pair ("gen-maxh",
						vL10AbfsConstMaxh));
  
  iResFtrSnp.mWeightedAbfs.insert (make_pair ("gen",
					      log10_weighted_sum (
						&(vL10AbfsConst[0]),
						vL10AbfsConst.size())));
  iResFtrSnp.mWeightedAbfs.insert (make_pair ("gen-fix",
					      log10_weighted_sum (
						&(vL10AbfsConstFix[0]),
						vL10AbfsConstFix.size())));
  iResFtrSnp.mWeightedAbfs.insert (make_pair ("gen-maxh",
					      log10_weighted_sum (
						&(vL10AbfsConstMaxh[0]),
						vL10AbfsConstMaxh.size())));
}

void
ResFtrSnp_calcAbfsSingletonsMvlr (
  ResFtrSnp & iResFtrSnp,
  const Grid & iGridS,
  const float & propFitSigma,
  vector<vector<double> > & Y,
  vector<vector<double> > & Xg,
  vector<vector<double> > & Xc)
{
  stringstream ssConfig;
  vector<vector<int> > vvGamma (1, vector<int> ());
  
  size_t i = 0; // for cases where ftr is absent in some subgroups
  for (size_t s = 0; s < iResFtrSnp.vNs.size(); ++s)
  {
    if (iResFtrSnp.vNs[s] < 2)
      continue;
    
    ssConfig.str("");
    ssConfig << (s+1);
    
    vvGamma[0].assign (Y.size(), 0);
    vvGamma[0][i] = 1;
    ++i;
    
#ifdef LIB_MVLR
    MVLR iMvlr;
    iMvlr.set_sigma_option (propFitSigma);
    iMvlr.init (Y, Xg, Xc);
    iMvlr.set_effect_vec (iGridS.phi2s, iGridS.oma2s);
    vector<double> vL10Abfs = iMvlr.compute_log10_ABF_vec (vvGamma);
    iResFtrSnp.mUnweightedAbfs.insert (make_pair (ssConfig.str(),
						  vL10Abfs));
    iResFtrSnp.mWeightedAbfs.insert (make_pair(ssConfig.str(),
					       log10_weighted_sum (
						 &(vL10Abfs[0]),
						 vL10Abfs.size())));
#endif
  }
}

void
ResFtrSnp_calcAbfsAllConfigsMvlr (
  ResFtrSnp & iResFtrSnp,
  const Grid & iGridS,
  const float & propFitSigma,
  vector<vector<double> > & Y,
  vector<vector<double> > & Xg,
  vector<vector<double> > & Xc)
{
  gsl_combination * comb;
  stringstream ssConfig;
  vector<bool> vIsEqtlInConfig; // T,T,F if S=3 and config="1-2"
  vector<vector<int> > vvGamma (1, vector<int> ());
  
  for (size_t k = 1; k <= iResFtrSnp.vNs.size(); ++k)
  {
    comb = gsl_combination_calloc (iResFtrSnp.vNs.size(), k);
    if (comb == NULL)
    {
      cerr << "ERROR: can't allocate memory for the combination" << endl;
      exit (1);
    }
    
    while (true)
    {
      bool isFtrAbsentInOneSubgroup = false; // skip config if true
      for (size_t i = 0; i < comb->k; ++i)
	if (iResFtrSnp.vNs[gsl_combination_get (comb, i)] < 2)
	{
	  isFtrAbsentInOneSubgroup = true;
	  break;
	}
      if (isFtrAbsentInOneSubgroup)
      {
	if (gsl_combination_next (comb) != GSL_SUCCESS)
	  break;
	continue;
      }
      
      prepareConfig (ssConfig, vIsEqtlInConfig, comb);
      vvGamma[0].assign (Y.size(), 0);
      size_t i = 0;
      for (size_t s = 0; s < vIsEqtlInConfig.size(); ++s)
	if (vIsEqtlInConfig[s])
	{
	  vvGamma[0][i] = 1;
	  ++i;
	}
#ifdef LIB_MVLR
      MVLR iMvlr;
      iMvlr.set_sigma_option (propFitSigma);
      iMvlr.init (Y, Xg, Xc);
      iMvlr.set_effect_vec (iGridS.phi2s, iGridS.oma2s);
      vector<double> vL10Abfs = iMvlr.compute_log10_ABF_vec (vvGamma);
      iResFtrSnp.mUnweightedAbfs.insert (make_pair (ssConfig.str(),
						    vL10Abfs));
      iResFtrSnp.mWeightedAbfs.insert (make_pair(ssConfig.str(),
						 log10_weighted_sum (
						   &(vL10Abfs[0]),
						   vL10Abfs.size())));
#endif
      if (gsl_combination_next (comb) != GSL_SUCCESS)
	break;
    }
    gsl_combination_free (comb);
  }
}

void
ResFtrSnp_calcAbfsMvlr (
  ResFtrSnp & iResFtrSnp,
  vector<vector<double> > & Y,
  vector<vector<double> > & Xg,
  vector<vector<double> > & Xc,
  const string & whichBfs,
  const Grid & iGridL,
  const Grid & iGridS,
  const float & propFitSigma)
{
  ResFtrSnp_calcAbfsConstLargeGridMvlr (iResFtrSnp, iGridL, propFitSigma,
					Y, Xg, Xc);
  if (whichBfs.compare("sin") == 0)
  {
    ResFtrSnp_calcAbfsSingletonsMvlr (iResFtrSnp, iGridS, propFitSigma,
				      Y, Xg, Xc);
    ResFtrSnp_calcAbfsAvgSinAndGenSin (iResFtrSnp);
  }
  else if (whichBfs.compare("all") == 0)
  {
    ResFtrSnp_calcAbfsAllConfigsMvlr (iResFtrSnp, iGridS, propFitSigma,
				      Y, Xg, Xc);
    ResFtrSnp_calcAbfsAvgSinAndGenSin (iResFtrSnp);
    ResFtrSnp_calcAbfAvgAll (iResFtrSnp);
  }
}

void
Ftr_init (
  Ftr & iFtr,
  const string & name,
  const size_t & nbSubgroups)
{
  iFtr.name = name;
  iFtr.chr.clear();
  iFtr.start = string::npos;
  iFtr.end = string::npos;
  iFtr.vvPhenos.resize (nbSubgroups);
  iFtr.vvIsNa.resize (nbSubgroups);
  iFtr.vPermPvalsSep.assign (nbSubgroups, numeric_limits<double>::quiet_NaN());
  iFtr.vNbPermsSoFar.assign (nbSubgroups, 0);
  iFtr.vMinTruePvals.assign (nbSubgroups, numeric_limits<double>::quiet_NaN());
  iFtr.sepPermPval = numeric_limits<double>::quiet_NaN();
  iFtr.nbPermsSoFarSep = 0;
  iFtr.minTruePval = numeric_limits<double>::quiet_NaN();
  iFtr.jointPermPval = numeric_limits<double>::quiet_NaN();
  iFtr.nbPermsSoFarJoint = 0;
  iFtr.maxL10TrueAbf = numeric_limits<double>::quiet_NaN();
  iFtr.avgL10TrueAbf = numeric_limits<double>::quiet_NaN();
}

// assume both features are on the same chromosome
bool Ftr_compByCoord (
  const Ftr* pt_iFtr1,
  const Ftr* pt_iFtr2)
{
  bool res = false;
  if ((pt_iFtr1->start < pt_iFtr2->start) ||
      (pt_iFtr1->start == pt_iFtr2->start && 
       pt_iFtr1->end < pt_iFtr2->end ))
      res = true;
  return res;
}

void
Ftr_getCisSnps (
  Ftr & iFtr,
  const map<string, vector<Snp*> > & mChr2VecPtSnps,
  const string & anchor,
  const size_t & lenCis)
{
  map<string, vector<Snp*> >::const_iterator itVecPtSnps =
    mChr2VecPtSnps.find(iFtr.chr);
  if (itVecPtSnps == mChr2VecPtSnps.end())
    cerr << "WARNING: feature " << iFtr.name << " is on chr '" << iFtr.chr
	 << "', is it encoded identically in the SNP coords file?" << endl;
  else
  {
    for (size_t snpId = 0; snpId < itVecPtSnps->second.size(); ++snpId)
    {
      Snp * ptSnp = (itVecPtSnps->second)[snpId];
      int inCis = Snp_isInCis (*ptSnp, iFtr.start, iFtr.end,
			       anchor, lenCis);
      if (inCis == 1)
	break;
      else if (inCis == -1)
	continue;
      iFtr.vPtCisSnps.push_back ((itVecPtSnps->second)[snpId]);
    }
  }
}

size_t
Ftr_getNbCisSnpsInGivenSubgroup (
  const Ftr & iFtr,
  const size_t & s)
{
  size_t nbCisSnps = 0;
  for (vector<ResFtrSnp>::const_iterator itP = iFtr.vResFtrSnps.begin();
       itP != iFtr.vResFtrSnps.end(); ++itP)
    if (itP->vNs[s] > 2)
      ++nbCisSnps;
  return nbCisSnps;
}

/** \brief Infer associations between a given feature and its SNPs in cis
 */
void
Ftr_inferAssos (
  Ftr & iFtr,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const int & whichStep,
  const bool & needQnorm,
  const vector<map<string, vector<double> > > & vSbgrp2Covars,
  const Grid & iGridL,
  const Grid & iGridS,
  const string & whichBfs,
  const bool & mvlr,
  const float & propFitSigma,
  const int & verbose)
{
  size_t nbSubgroups = iFtr.vvPhenos.size();
  for (size_t snpId = 0; snpId < iFtr.vPtCisSnps.size(); ++snpId)
  {
    if (verbose > 0)
      cout << iFtr.name << " " << iFtr.vPtCisSnps[snpId]->name << endl;
    ResFtrSnp iResFtrSnp;
    ResFtrSnp_init (iResFtrSnp, iFtr.vPtCisSnps[snpId]->name, nbSubgroups);
    if (whichStep == 1 || whichStep == 2 ||
	(! mvlr && (whichStep == 3 || whichStep == 4 || whichStep == 5)))
    {
      for (size_t s = 0; s < nbSubgroups; ++s)
      {
	if (iFtr.vvPhenos[s].size() > 0 &&
	    iFtr.vPtCisSnps[snpId]->vvGenos[s].size() > 0)
	{
	  ResFtrSnp_getSstatsOneSbgrp (iResFtrSnp, iFtr,
				       *(iFtr.vPtCisSnps[snpId]), s,
				       vvSampleIdxPhenos, vvSampleIdxGenos,
				       needQnorm, vSbgrp2Covars);
	}
      }
      if (whichStep == 3 || whichStep == 4 || whichStep == 5)
	ResFtrSnp_calcAbfs (iResFtrSnp, whichBfs, iGridL, iGridS);
    }
    else // multivariate model with BFs only (no summary stats)
    {
      vector<vector<double> > Y, Xg, Xc;
      ResFtrSnp_prepareDataForMvlr (iResFtrSnp, iFtr, *(iFtr.vPtCisSnps[snpId]),
				    vvSampleIdxPhenos, vvSampleIdxGenos,
				    vSbgrp2Covars, needQnorm, Y, Xg, Xc);
      ResFtrSnp_calcAbfsMvlr (iResFtrSnp, Y, Xg, Xc, whichBfs, iGridL, iGridS,
			      propFitSigma);
    }
    iFtr.vResFtrSnps.push_back (iResFtrSnp);
  }
}

/** \brief Retrieve the lowest genotype P-value over SNPs of the given feature
 *  for the given subgroup
 */
void
Ftr_findMinTrueGenoPvals (
  Ftr & iFtr,
  const size_t & s)
{
  iFtr.vMinTruePvals[s] = 2;
  for (size_t snpId = 0; snpId < iFtr.vResFtrSnps.size(); ++snpId)
  {
    if (iFtr.vResFtrSnps[snpId].vNs[s] > 2 &&
	iFtr.vResFtrSnps[snpId].vMapPredictors[s]["genotype"][2]
	< iFtr.vMinTruePvals[s])
      iFtr.vMinTruePvals[s] =
	iFtr.vResFtrSnps[snpId].vMapPredictors[s]["genotype"][2];
  }
  if (iFtr.vMinTruePvals[s] == 2)
    iFtr.vMinTruePvals[s] = numeric_limits<double>::quiet_NaN();
}

/** \brief Make permutations for the separate analysis for a given feature
 *  for a given subgroup (index s)
 */
void
Ftr_makePermsSepOneSubgrp (
  Ftr & iFtr,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const bool & needQnorm,
  const vector<map<string, vector<double> > > & vSbgrp2Covars,
  const size_t & nbPerms,
  const int & trick,
  const size_t & s,
  const gsl_rng * & rngPerm,
  const gsl_rng * & rngTrick)
{
  gsl_permutation * perm = NULL;
  double minPermBetaPval;
  bool shuffleOnly = false;
  Snp iSnp;
  
  Ftr_findMinTrueGenoPvals (iFtr, s);
  if (! isNan (iFtr.vMinTruePvals[s]))
  {
    iFtr.vPermPvalsSep[s] = 1;
    
    perm = gsl_permutation_calloc (vvSampleIdxPhenos[s].size());
    if (perm == NULL)
    {
      cerr << "ERROR: can't allocate memory for the permutation" << endl;
      exit (1);
    }
    
    for(size_t permId = 0; permId < nbPerms; ++permId)
    {
      gsl_ran_shuffle (rngPerm, perm->data, perm->size, sizeof(size_t));
      if (shuffleOnly)
	continue;
      
      ++(iFtr.vNbPermsSoFar[s]);
      minPermBetaPval = 1;
      
      for (size_t snpId = 0; snpId < iFtr.vPtCisSnps.size(); ++snpId)
      {
	iSnp = *(iFtr.vPtCisSnps[snpId]);
	ResFtrSnp iResFtrSnp;
	ResFtrSnp_init (iResFtrSnp, iSnp.name, 1);
	if (iFtr.vvPhenos[s].size() > 0 &&
	    iFtr.vPtCisSnps[snpId]->vvGenos[s].size() > 0)
	  ResFtrSnp_getSstatsPermOneSbgrp (iResFtrSnp, iFtr, iSnp, 0,
					   vvSampleIdxPhenos,
					   vvSampleIdxGenos,
					   needQnorm, vSbgrp2Covars, perm);
	if (iResFtrSnp.vNs[0] > 2 &&
	    iResFtrSnp.vMapPredictors[0]["genotype"][2] < minPermBetaPval)
	  minPermBetaPval = iResFtrSnp.vMapPredictors[0]["genotype"][2];
      }
      
      if (minPermBetaPval <= iFtr.vMinTruePvals[s])
	++(iFtr.vPermPvalsSep[s]);
      if (trick != 0 && iFtr.vPermPvalsSep[s] == 11)
      {
	if (trick == 1)
	  break;
	else if (trick == 2)
	  shuffleOnly = true;
      }
    }
    
    if (iFtr.vNbPermsSoFar[s] == nbPerms)
      iFtr.vPermPvalsSep[s] /= (nbPerms + 1);
    else
      iFtr.vPermPvalsSep[s] = gsl_ran_flat (rngTrick,
					    (11 / ((double) (iFtr.vNbPermsSoFar[s] + 2))),
					    (11 / ((double) (iFtr.vNbPermsSoFar[s] + 1))));
    
    gsl_permutation_free (perm);
  }
}

/** \brief Retrieve the lowest P-value over SNPs and subgroups
 *  of the given feature
 */
void
Ftr_findMinPvalOverSnpsAndSubgroups (
  Ftr & iFtr)
{
  size_t nbSubgroups = iFtr.vvPhenos.size();
  iFtr.minTruePval = 1;
  for (size_t snpId = 0; snpId < iFtr.vResFtrSnps.size(); ++snpId)
    for (size_t s = 0; s < nbSubgroups; ++s)
      if (iFtr.vResFtrSnps[snpId].vMapPredictors[s]["genotype"][2] < iFtr.minTruePval)
	iFtr.minTruePval = iFtr.vResFtrSnps[snpId].vMapPredictors[s]["genotype"][2];
}

/** \brief Make permutations for the separate analysis for a given feature
 *  using the minimum P-value over SNPs and subgroups as a test statistic
 */
void
Ftr_makePermsSepAllSubgrps (
  Ftr & iFtr,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const bool & needQnorm,
  const vector<map<string, vector<double> > > & vSbgrp2Covars,
  const size_t & nbPerms,
  const int & trick,
  const gsl_rng * & rngPerm,
  const gsl_rng * & rngTrick)
{
  gsl_permutation * perm = NULL;
  
  perm = gsl_permutation_calloc (vvSampleIdxPhenos[0].size());
  if (perm == NULL)
  {
    cerr << "ERROR: can't allocate memory for the permutation" << endl;
    exit (1);
  }
  
  iFtr.sepPermPval = 1;
  iFtr.nbPermsSoFarSep = 0;
  
  size_t nbSubgroups = iFtr.vvPhenos.size();
  double minPermPval, permPval;
  bool shuffleOnly = false;
  Snp iSnp;
  
  Ftr_findMinPvalOverSnpsAndSubgroups (iFtr);
  
  for(size_t permId = 0; permId < nbPerms; ++permId)
  {
    gsl_ran_shuffle (rngPerm, perm->data, perm->size, sizeof(size_t));
    if (shuffleOnly)
      continue;
    
    ++iFtr.nbPermsSoFarSep;
    minPermPval = 1;
    
    for (size_t snpId = 0; snpId < iFtr.vPtCisSnps.size(); ++snpId)
    {
      iSnp = *(iFtr.vPtCisSnps[snpId]);
      ResFtrSnp iResFtrSnp;
      ResFtrSnp_init (iResFtrSnp, iSnp.name, nbSubgroups);
      for (size_t s = 0; s < nbSubgroups; ++s)
      {
	if (iFtr.vvPhenos[s].size() > 0 &&
	    iFtr.vPtCisSnps[snpId]->vvGenos[s].size() > 0)
	  ResFtrSnp_getSstatsPermOneSbgrp (iResFtrSnp, iFtr, iSnp, s,
					   vvSampleIdxPhenos,
					   vvSampleIdxGenos,
					   needQnorm, vSbgrp2Covars, perm);
	permPval = iResFtrSnp.vMapPredictors[s]["genotype"][2];
	if (permPval < minPermPval)
	  minPermPval = permPval;
      }
    }
    
    if (minPermPval <= iFtr.minTruePval)
      ++iFtr.sepPermPval;
    if (trick != 0 && iFtr.sepPermPval == 11)
    {
      if (trick == 1)
	break;
      else if (trick == 2)
	shuffleOnly = true;
    }
  }
  
  if (iFtr.nbPermsSoFarSep == nbPerms)
    iFtr.sepPermPval /= (nbPerms + 1);
  else
    iFtr.sepPermPval = gsl_ran_flat (rngTrick,
				     (11 / ((double) (iFtr.nbPermsSoFarSep + 2))),
				     (11 / ((double) (iFtr.nbPermsSoFarSep + 1))));
  
  gsl_permutation_free (perm);
}

/** \brief Retrieve the highest log10(ABF) over SNPs of the given feature
 *  \note whichPermBf is 'gen', 'sin', 'gen-sin' or 'all'
 */
void
Ftr_findMaxL10TrueAbf (
 Ftr & iFtr,
  const string & whichPermBf)
{
  iFtr.maxL10TrueAbf = - numeric_limits<double>::infinity();
  for (size_t snpId = 0; snpId < iFtr.vResFtrSnps.size(); ++snpId)
    if (iFtr.vResFtrSnps[snpId].mWeightedAbfs[whichPermBf] > iFtr.maxL10TrueAbf)
      iFtr.maxL10TrueAbf = iFtr.vResFtrSnps[snpId].mWeightedAbfs[whichPermBf];
}

/** \brief Average the log10(ABF) over SNPs of the given feature
 *  \note whichPermBf is 'gen', 'sin', 'gen-sin' or 'all'
 */
void
Ftr_avgL10TrueAbfs (
  Ftr & iFtr,
  const string & whichPermBf)
{
  vector<double> vL10Abfs;
  for (size_t snpId = 0; snpId < iFtr.vResFtrSnps.size(); ++snpId)
    if (! isNan (iFtr.vResFtrSnps[snpId].mWeightedAbfs[whichPermBf]))
      vL10Abfs.push_back (iFtr.vResFtrSnps[snpId].mWeightedAbfs[whichPermBf]);
  iFtr.avgL10TrueAbf = log10_weighted_sum (&(vL10Abfs[0]), vL10Abfs.size());
}

/** \brief Make permutations for the joint analysis for a given feature
 *  using as a test statistic the BF precised by whichPermBf
 */
void
Ftr_makePermsJoint (
  Ftr & iFtr,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const bool & needQnorm,
  const vector<map<string, vector<double> > > & vSbgrp2Covars,
  const Grid & iGridL,
  const Grid & iGridS,
  const bool & mvlr,
  const float & propFitSigma,
  const size_t & nbPerms,
  const int & trick,
  const string & whichPermBf,
  const bool & useMaxBfOverSnps,
  const gsl_rng * & rngPerm,
  const gsl_rng * & rngTrick)
{
  gsl_permutation * perm = NULL;
  
  perm = gsl_permutation_calloc (vvSampleIdxPhenos[0].size());
  if (perm == NULL)
  {
    cerr << "ERROR: can't allocate memory for the permutation" << endl;
    exit (1);
  }
  
  iFtr.jointPermPval = 1;
  iFtr.nbPermsSoFarJoint = 0;
  
  size_t nbSubgroups = iFtr.vvPhenos.size();
  double maxL10PermAbf, l10PermAbf, avgL10PermAbf;
  vector<double> vL10PermAbfs;
  bool shuffleOnly = false;
  Snp iSnp;
  
  if (useMaxBfOverSnps)
    Ftr_findMaxL10TrueAbf (iFtr, whichPermBf);
  else
    Ftr_avgL10TrueAbfs (iFtr, whichPermBf);
  
  for(size_t permId = 0; permId < nbPerms; ++permId)
  {
    gsl_ran_shuffle (rngPerm, perm->data, perm->size, sizeof(size_t));
    if (shuffleOnly)
      continue;
    
    ++iFtr.nbPermsSoFarJoint;
    if (useMaxBfOverSnps)
      maxL10PermAbf = - numeric_limits<double>::infinity();
    else // if avg over SNPs
    {
      vL10PermAbfs.clear ();
      avgL10PermAbf = numeric_limits<double>::quiet_NaN();
    }
    
    for (size_t snpId = 0; snpId < iFtr.vPtCisSnps.size(); ++snpId)
    {
      iSnp = *(iFtr.vPtCisSnps[snpId]);
      ResFtrSnp iResFtrSnp;
      ResFtrSnp_init (iResFtrSnp, iSnp.name, nbSubgroups);
      if (! mvlr) // univariate model
      {
	for (size_t s = 0; s < nbSubgroups; ++s)
	  if (iFtr.vvPhenos[s].size() > 0 &&
	      iFtr.vPtCisSnps[snpId]->vvGenos[s].size() > 0)
	    ResFtrSnp_getSstatsPermOneSbgrp (iResFtrSnp, iFtr, iSnp, s,
					     vvSampleIdxPhenos,
					     vvSampleIdxGenos,
					     needQnorm, vSbgrp2Covars, perm);
	vector<vector<double> > vvStdSstats;
	ResFtrSnp_getStdStatsAndCorrSmallSampleSize (iResFtrSnp,
						     vvStdSstats);
	if (whichPermBf.compare("gen") == 0)
	  ResFtrSnp_calcAbfsConstLargeGrid (iResFtrSnp, iGridL, vvStdSstats);
	else if (whichPermBf.compare("sin") == 0)
	{
	  ResFtrSnp_calcAbfsSingletons (iResFtrSnp, iGridS, vvStdSstats);
	  ResFtrSnp_calcAbfsAvgSinAndGenSin (iResFtrSnp);
	}
	else if (whichPermBf.compare("gen-sin") == 0)
	{
	  ResFtrSnp_calcAbfsConstLargeGrid (iResFtrSnp, iGridL, vvStdSstats);
	  ResFtrSnp_calcAbfsSingletons (iResFtrSnp, iGridS, vvStdSstats);
	  ResFtrSnp_calcAbfsAvgSinAndGenSin (iResFtrSnp);
	}
	else if (whichPermBf.compare("all") == 0)
	{
	  ResFtrSnp_calcAbfsAllConfigs (iResFtrSnp, iGridS, vvStdSstats);
	  ResFtrSnp_calcAbfAvgAll (iResFtrSnp);
	}
      }
      else // multivariate model
      {
	vector<vector<double> > Y, Xg, Xc;
	ResFtrSnp_prepareDataForMvlrPerm (iResFtrSnp, iFtr, iSnp,
					  vvSampleIdxPhenos, vvSampleIdxGenos,
					  vSbgrp2Covars, needQnorm, perm,
					  Y, Xg, Xc);
	if (whichPermBf.compare("gen") == 0)
	  ResFtrSnp_calcAbfsConstLargeGridMvlr (iResFtrSnp, iGridL, propFitSigma,
						Y, Xg, Xc);
	else if (whichPermBf.compare("sin") == 0)
	{
	  ResFtrSnp_calcAbfsSingletonsMvlr (iResFtrSnp, iGridS, propFitSigma,
					    Y, Xg, Xc);
	  ResFtrSnp_calcAbfsAvgSinAndGenSin (iResFtrSnp);
	}
	else if (whichPermBf.compare("gen-sin") == 0)
	{
	  ResFtrSnp_calcAbfsConstLargeGridMvlr (iResFtrSnp, iGridL, propFitSigma,
						Y, Xg, Xc);
	  ResFtrSnp_calcAbfsSingletonsMvlr (iResFtrSnp, iGridS, propFitSigma,
					    Y, Xg, Xc);
	  ResFtrSnp_calcAbfsAvgSinAndGenSin (iResFtrSnp);
	}
	else if (whichPermBf.compare("all") == 0)
	{
	  ResFtrSnp_calcAbfsAllConfigsMvlr (iResFtrSnp, iGridS, propFitSigma,
					    Y, Xg, Xc);
	  ResFtrSnp_calcAbfAvgAll (iResFtrSnp);
	}
      }
      if (useMaxBfOverSnps)
      {
	l10PermAbf = iResFtrSnp.mWeightedAbfs[whichPermBf];
	if (l10PermAbf > maxL10PermAbf)
	  maxL10PermAbf = l10PermAbf;
      }
      else // if avg over SNPs
	if (! isNan (iResFtrSnp.mWeightedAbfs[whichPermBf]))
	  vL10PermAbfs.push_back (iResFtrSnp.mWeightedAbfs[whichPermBf]);
    }
    
    if (useMaxBfOverSnps)
    {
      if (maxL10PermAbf >= iFtr.maxL10TrueAbf)
	++iFtr.jointPermPval;
    }
    else // if avg over SNPs
    {
      avgL10PermAbf = log10_weighted_sum (&(vL10PermAbfs[0]), vL10PermAbfs.size());
      if (avgL10PermAbf >= iFtr.avgL10TrueAbf)
	++iFtr.jointPermPval;
    }
    if (trick != 0 && iFtr.jointPermPval == 11)
    {
      if (trick == 1)
	break;
      else if (trick == 2)
	shuffleOnly = true;
    }
  }
  
  if (iFtr.nbPermsSoFarJoint == nbPerms)
    iFtr.jointPermPval /= (nbPerms + 1);
  else
    iFtr.jointPermPval = gsl_ran_flat (rngTrick,
				       (11 / ((double) (iFtr.nbPermsSoFarJoint + 2))),
				       (11 / ((double) (iFtr.nbPermsSoFarJoint + 1))));
  
  gsl_permutation_free (perm);
}

void
loadListsGenoAndPhenoFiles (
  const string & genoPathsFile,
  const string & phenoPathsFile,
  const string & sbgrpToKeep,
  map<string, string> & mGenoPaths,
  map<string, string> & mPhenoPaths,
  vector<string> & vSubgroups,
  const bool & mvlr,
  const int & verbose)
{
  mPhenoPaths = loadTwoColumnFile (phenoPathsFile, verbose);
  for (map<string, string>::const_iterator it = mPhenoPaths.begin();
       it != mPhenoPaths.end(); ++it)
    if (sbgrpToKeep.empty() || sbgrpToKeep.compare(it->first) == 0)
      vSubgroups.push_back (it->first);
  
  mGenoPaths = loadTwoColumnFile (genoPathsFile, verbose);
  if (mGenoPaths.size() > 1 && mGenoPaths.size() != mPhenoPaths.size())
  {
    cerr << "ERROR: there should be only one genotype file"
	 << " or as many as phenotype files" << endl;
    exit (1);
  }
  if (mvlr)
    for (map<string, string>::const_iterator it = mGenoPaths.begin();
	 it != mGenoPaths.end(); ++it)
      if (it->second.compare(mGenoPaths.begin()->second) != 0)
      {
	cerr << "ERROR: --mvlr requires the same genotypes in all subgroups"
	     << endl;
	exit (1);
      }
  
  if (! sbgrpToKeep.empty())
    if (verbose > 0)
      cout << "analyze a single subgroup: " << sbgrpToKeep << endl;
}

void
loadSamplesAllPhenos (
  const map<string, string> & mPhenoPaths,
  const vector<string> & vSubgroups,
  vector<string> & vSamples,
  vector<vector<string> > & vvSamples,
  const int & verbose)
{
  gzFile fileStream;
  string line;
  for (size_t s = 0; s < vSubgroups.size(); ++s)
  {
    openFile (mPhenoPaths.find(vSubgroups[s])->second, fileStream, "rb");
    if (! getline (fileStream, line))
    {
      cerr << "ERROR: problem with the header of file "
	   << mPhenoPaths.find(vSubgroups[s])->second << endl;
      exit (1);
    }
    if (line.empty())
    {
      cerr << "ERROR: file " << mPhenoPaths.find(vSubgroups[s])->second
	   << " is empty" << endl;
      exit (1);
    }
    closeFile (mPhenoPaths.find(vSubgroups[s])->second, fileStream);
    if (s == 0)
    {
      split (line, " \t", vSamples);
      if (vSamples[0].compare("Id") == 0)
	vSamples.erase (vSamples.begin());
      vvSamples.push_back (vSamples);
    }
    else
    {
      vector<string> tokens;
      split (line, " \t", tokens);
      if (tokens[0].compare("Id") == 0)
	tokens.erase (tokens.begin());
      vvSamples.push_back (tokens);
      for (size_t i = 0; i < tokens.size(); ++i)
	if (find (vSamples.begin(), vSamples.end(), tokens[i])
	    == vSamples.end())
	  vSamples.push_back (tokens[i]);
    }
  }
  if (verbose > 0)
  {
    cout << "nb of samples (phenotypes): "
	 << vSamples.size() << endl << flush;
    for (size_t s = 0; s < vSubgroups.size(); ++s)
    {
      cout << "s" << (s+1) << " (" << vSubgroups[s] << "): "
	   << vvSamples[s].size() << " samples" << endl << flush;
      if (verbose > 1)
      {
	for (size_t i = 0; i < vvSamples[s].size(); ++i)
	  cout << vvSamples[s][i] << endl;
      }
    }
  }
}

void
loadSamplesAllGenos (
  const map<string, string> & mGenoPaths,
  const vector<string> & vSubgroups,
  vector<string> & vSamples,
  vector<vector<string> > & vvSamples,
  const int & verbose)
{
  gzFile fileStream;
  string line, sample;
  vector<string> tokens, tokens2;
  size_t i;
  for (size_t s = 0; s < vSubgroups.size(); ++s)
  {
    if (mGenoPaths.find(vSubgroups[s]) == mGenoPaths.end())
    {
      cerr << "ERROR: can't find genotypes for subgroup " << vSubgroups[s]
	   << endl;
      exit (1);
    }
    openFile (mGenoPaths.find(vSubgroups[s])->second, fileStream, "rb");
    if (! getline (fileStream, line))
    {
      cerr << "ERROR: problem with the header of file "
	   << mGenoPaths.find(vSubgroups[s])->second << endl;
      exit (1);
    }
    if (line.empty())
    {
      cerr << "ERROR: file " << mGenoPaths.find(vSubgroups[s])->second
	   << " is empty" << endl;
      exit (1);
    }
    
    if (line.find("##fileformat=VCF") != string::npos) // VCF format
    {
      while (getline (fileStream, line))
      {
	if (line.find("#CHROM") == string::npos)
	  continue;
	split (line, " \t", tokens);
	closeFile (mGenoPaths.find(vSubgroups[s])->second, fileStream);
	vvSamples.push_back (vector<string> (tokens.size()-9));
	for (size_t i = 9; i < tokens.size(); ++i)
	{
	  vvSamples[s][i-9] = tokens[i];
	  if (find (vSamples.begin(), vSamples.end(), tokens[i])
	      == vSamples.end())
	    vSamples.push_back (tokens[i]);
	}
	break;
      }
    }
    else // not VCF
    {
      split (line, " \t", tokens);
      closeFile (mGenoPaths.find(vSubgroups[s])->second, fileStream);
      
      if (tokens[0].compare("chr") == 0) // IMPUTE format
      {
	if ((tokens.size() - 5) % 3 != 0)
	{
	  cerr << "ERROR: the header of IMPUTE file "
	       << mGenoPaths.find(vSubgroups[s])->second
	       << " is badly formatted" << endl;
	  exit (1);
	}
	tokens2.clear();
	i = 5;
	while (i < tokens.size())
	{
	  sample = split (tokens[i], "_a", 0); // sampleX_a1a1, sampleX_a1a2 or sampleX_a2a2
	  tokens2.push_back (sample);
	  i = i + 3;
	}
	vvSamples.push_back (tokens2);
	for (i = 0; i < tokens2.size(); ++i)
	  if (find (vSamples.begin(), vSamples.end(), tokens2[i])
	      == vSamples.end())
	    vSamples.push_back (tokens2[i]);
      }
      else // allele dosage
      {
	if (tokens[0].compare("Id") == 0)
	  tokens.erase (tokens.begin());
	vvSamples.push_back (tokens);
	for (i = 0; i < tokens.size(); ++i)
	  if (find (vSamples.begin(), vSamples.end(), tokens[i])
	      == vSamples.end())
	    vSamples.push_back (tokens[i]);
      }
    }
  }
  
  if (verbose > 0)
  {
    cout << "nb of samples (genotypes): "
	 << vSamples.size() << endl << flush;
    for (size_t s = 0; s < vSubgroups.size(); ++s)
    {
      cout << "s" << (s+1) << " (" << vSubgroups[s] << "): "
	   << vvSamples[s].size() << " samples" << endl << flush;
      if (verbose > 1)
      {
	for (size_t i = 0; i < vvSamples[s].size(); ++i)
	  cout << vvSamples[s][i] << endl;
      }
    }
  }
}

void
loadSamples (
  const map<string, string> & mGenoPaths,
  const map<string, string> & mPhenoPaths,
  const vector<string> & vSubgroups,
  vector<string> & vSamples,
  vector<vector<size_t> > & vvSampleIdxGenos,
  vector<vector<size_t> > & vvSampleIdxPhenos,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load samples ..." << endl << flush;
  
  vector<string> vAllSamplesPhenos;
  vector<vector<string> > vvSamplesPhenos;
  loadSamplesAllPhenos (mPhenoPaths, vSubgroups, vAllSamplesPhenos,
			vvSamplesPhenos, verbose);
  
  vector<string> vAllSamplesGenos;
  vector<vector<string> > vvSamplesGenos;
  loadSamplesAllGenos (mGenoPaths, vSubgroups, vAllSamplesGenos,
		       vvSamplesGenos, verbose);
  
  // fill vSamples by merging vAllSamplesPhenos and vAllSamplesGenos
  vector<string>::iterator it;
  for (it = vAllSamplesPhenos.begin();
       it != vAllSamplesPhenos.end(); ++it)
    if (find(vSamples.begin(), vSamples.end(), *it) == vSamples.end())
      vSamples.push_back (*it);
  for (it = vAllSamplesGenos.begin();
       it != vAllSamplesGenos.end(); ++it)
    if (find(vSamples.begin(), vSamples.end(), *it) == vSamples.end())
      vSamples.push_back (*it);
  if (verbose > 0)
    cout << "total nb of samples: "
	 << vSamples.size() << endl << flush;
  
  // fill vvSampleIdxPhenos
  // vvSampleIdxPhenos[3][0] = 5 means that the 1st sample in vSamples
  // corresponds to the 6th sample in the 4th subgroup
  // it is npos if the sample is absent from this subgroup
  for (size_t s = 0; s < vvSamplesPhenos.size(); ++s)
  {
    vector<size_t> vSampleIdxPhenos (vSamples.size(), string::npos);
    for (size_t i = 0; i < vSamples.size(); ++i)
    {
      vector<string>::iterator it = find(vvSamplesPhenos[s].begin(),
					 vvSamplesPhenos[s].end(),
					 vSamples[i]);
      if (it != vvSamplesPhenos[s].end())
	vSampleIdxPhenos[i] = it - vvSamplesPhenos[s].begin();
    }
    vvSampleIdxPhenos.push_back (vSampleIdxPhenos);
  }
  
  // fill vvSampleIdxGenos
  for (size_t s = 0; s < vvSamplesGenos.size(); ++s)
  {
    vector<size_t> vSampleIdxGenos (vSamples.size(), string::npos);
    for (size_t i = 0; i < vSamples.size(); ++i)
    {
      vector<string>::iterator it = find(vvSamplesGenos[s].begin(),
					 vvSamplesGenos[s].end(),
					 vSamples[i]);
      if (it != vvSamplesGenos[s].end())
	vSampleIdxGenos[i] = it - vvSamplesGenos[s].begin();
    }
    vvSampleIdxGenos.push_back (vSampleIdxGenos);
  }
}

void
loadPhenos (
  const map<string, string> & mPhenoPaths,
  const vector<string> & vSubgroups,
  const vector<string> & vFtrsToKeep,
  map<string, Ftr> & mFtrs,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load phenotypes ..." << endl << flush;
  
  gzFile phenoStream;
  string line;
  vector<string> tokens;
  size_t nbSamples, nbLines, nbFtrsToKeepPerSubgroup;
  
  for (size_t s = 0; s < vSubgroups.size(); ++s)
  {
    nbFtrsToKeepPerSubgroup = 0;
    nbLines = 0;
    openFile (mPhenoPaths.find(vSubgroups[s])->second, phenoStream, "rb");
    if (! getline (phenoStream, line))
    {
      cerr << "ERROR: problem with the header of file "
	   << mPhenoPaths.find(vSubgroups[s])->second << endl;
      exit (1);
    }
    ++nbLines;
    split (line, " \t", tokens);
    if (tokens[0].compare("Id") == 0)
      nbSamples = tokens.size() - 1;
    else
      nbSamples = tokens.size();
    
    while (getline (phenoStream, line))
    {
      ++nbLines;
      split (line, " \t", tokens);
      if (tokens.size() != nbSamples + 1)
      {
	cerr << "ERROR: not enough columns on line " << nbLines
	     << " of file " << mPhenoPaths.find(vSubgroups[s])->second
	     << " (" << tokens.size() << " != " << nbSamples + 1 << ")"
	     << endl;
	exit (1);
      }
      if (! vFtrsToKeep.empty()
	  && find (vFtrsToKeep.begin(), vFtrsToKeep.end(), tokens[0])
	  == vFtrsToKeep.end())
	continue;
      ++nbFtrsToKeepPerSubgroup;
      
      if (mFtrs.find(tokens[0]) == mFtrs.end())
      {
	Ftr iFtr;
	Ftr_init (iFtr, tokens[0], mPhenoPaths.size());
	iFtr.vvPhenos[s].assign (nbSamples, numeric_limits<double>::quiet_NaN());
	iFtr.vvIsNa[s].assign (nbSamples, true);
	for (size_t i = 1; i < tokens.size(); ++i)
	  if (tokens[i].compare("NA") != 0)
	  {
	    iFtr.vvPhenos[s][i-1] = atof (tokens[i].c_str());
	    iFtr.vvIsNa[s][i-1] = false;
	  }
	mFtrs.insert (make_pair (tokens[0], iFtr));
      }
      else
      {
	mFtrs[tokens[0]].vvPhenos[s].assign (nbSamples, numeric_limits<double>::quiet_NaN());
	mFtrs[tokens[0]].vvIsNa[s].assign (nbSamples, true);
	for (size_t i = 1; i < tokens.size() ; ++i)
	  if (tokens[i].compare("NA") != 0)
	  {
	    mFtrs[tokens[0]].vvPhenos[s][i-1] = atof (tokens[i].c_str());
	    mFtrs[tokens[0]].vvIsNa[s][i-1] = false;
	  }
      }
    }
    if (! gzeof (phenoStream))
    {
      cerr << "ERROR: can't read successfully file "
	   << mPhenoPaths.find(vSubgroups[s])->second
	   << " up to the end" << endl;
      exit (1);
    }
    
    closeFile (mPhenoPaths.find(vSubgroups[s])->second, phenoStream);
    if (verbose > 0)
      cout << "s" << (s+1) << " (" << vSubgroups[s] << "): " << (nbLines-1)
	   << " features (to keep: " << nbFtrsToKeepPerSubgroup << ")"
	   << endl << flush;
  }
  
  if (mFtrs.size() == 0)
  {
    cerr << "ERROR: no feature to analyze" << endl;
    exit (1);
  }
  if (verbose > 0)
    cout << "total nb of features with phenotypes: " << mFtrs.size() << endl;
}

void
loadFtrInfo (
  const string & ftrCoordsFile,
  map<string, Ftr> & mFtrs,
  map<string, vector<Ftr*> > & mChr2VecPtFtrs,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load feature coordinates ..." << endl << flush;
  
  // parse the BED file
  gzFile ftrCoordsStream;
  openFile (ftrCoordsFile, ftrCoordsStream, "rb");
  string line;
  vector<string> tokens;
  while (getline (ftrCoordsStream, line))
  {
    split (line, " \t", tokens);
    if (mFtrs.find(tokens[3]) == mFtrs.end())
      continue;
    mFtrs[tokens[3]].chr = tokens[0];
    mFtrs[tokens[3]].start = atol (tokens[1].c_str()) + 1;
    mFtrs[tokens[3]].end = atol (tokens[2].c_str());
    
    if (mChr2VecPtFtrs.find(tokens[0]) == mChr2VecPtFtrs.end())
      mChr2VecPtFtrs.insert (make_pair (tokens[0],
					vector<Ftr*> ()));
    mChr2VecPtFtrs[tokens[0]].push_back (&(mFtrs[tokens[3]]));
  }
  if (! gzeof (ftrCoordsStream))
  {
    cerr << "ERROR: can't read successfully file "
	 << ftrCoordsFile
	 << " up to the end" << endl;
    exit (1);
  }
  closeFile (ftrCoordsFile, ftrCoordsStream);
  
  // check that all features have coordinates
  map<string, Ftr>::iterator it = mFtrs.begin();
  while (it != mFtrs.end())
  {
    if (it->second.chr.empty())
    {
      cerr << "WARNING: skip feature " << it->second.name
	   << " because it has no coordinate" << endl;
      mFtrs.erase (it++);
    }
    else
      ++it;
  }
  
  // sort the features per chr
  for (map<string, vector<Ftr*> >::iterator it = mChr2VecPtFtrs.begin();
       it != mChr2VecPtFtrs.end(); ++it)
    sort (it->second.begin(), it->second.end(), Ftr_compByCoord);
  
  if (mFtrs.size() == 0)
  {
    cerr << "ERROR: no feature to analyze" << endl;
    exit (1);
  }
  if (verbose > 0)
    cout << "total nb of features to analyze: " << mFtrs.size() << endl;
}

void
loadGenosAndSnpInfoFromImpute (
  gzFile & genoStream,
  string & line,
  size_t & nbLines,
  const map<string, string> & mGenoPaths,
  const float & minMaf,
  const vector<string> & vSubgroups,
  const size_t & s,
  const vector<string> & vSnpsToKeep,
  const map<string, vector<Ftr*> > & mChr2VecPtFtrs,
  size_t & nbSnpsToKeepPerSubgroup,
  map<string, Snp> & mSnps,
  map<string, vector<Snp*> > & mChr2VecPtSnps)
{
  vector<string> tokens;
  size_t nbSamples;
  double maf, AA, AB, BB;
  
  split (line, " \t", tokens);
  nbSamples = (size_t) (tokens.size() - 5) / 3;
  
  while (getline (genoStream, line))
  {
    ++nbLines;
    split (line, " \t", tokens);
    if (tokens.size() != (size_t) (3 * nbSamples + 5))
    {
      cerr << "ERROR: not enough columns on line " << nbLines
	   << " of file " << mGenoPaths.find(vSubgroups[s])->second
	   << " (" << tokens.size() << " != " << (3 * nbSamples + 5 ) << ")"
	   << endl;
      exit (1);
    }
    if (mChr2VecPtFtrs.find (tokens[0]) == mChr2VecPtFtrs.end())
      continue;
    if (! vSnpsToKeep.empty()
	&& find (vSnpsToKeep.begin(), vSnpsToKeep.end(), tokens[1])
	== vSnpsToKeep.end())
      continue;
    
    maf = 0;
    if (mSnps.find(tokens[1]) == mSnps.end())
    {
      Snp iSnp;
      Snp_init (iSnp, tokens[1], mGenoPaths.size());
      iSnp.vvIsNa[s].resize (nbSamples, false);
      iSnp.vvGenos[s].resize (nbSamples,
			      numeric_limits<double>::quiet_NaN());
      for (size_t i = 0; i < nbSamples; ++i)
      {
	AA = atof(tokens[5+3*i].c_str());
	AB = atof(tokens[5+3*i+1].c_str());
	BB = atof(tokens[5+3*i+2].c_str());
	if (AA == 0 && AB == 0 && BB == 0)
	  iSnp.vvIsNa[s][i] = true;
	else
	{
	  iSnp.vvGenos[s][i] = 0 * AA + 1 * AB + 2 * BB;
	  maf += iSnp.vvGenos[s][i];
	}
      }
      maf /= 2 * (nbSamples
		  - count (iSnp.vvIsNa[s].begin(),
			   iSnp.vvIsNa[s].end(),
			   true));
      if ((maf <= 0.5 && maf < minMaf) || (maf > 0.5 && 1-maf < minMaf))
	continue;
      ++nbSnpsToKeepPerSubgroup;
      iSnp.vMafs[s] = maf <= 0.5 ? maf : (1 - maf);
      iSnp.chr = tokens[0];
      iSnp.coord = atol (tokens[2].c_str());
      mSnps.insert (make_pair (tokens[1], iSnp));
      
      if (mChr2VecPtSnps.find(tokens[0]) == mChr2VecPtSnps.end())
	mChr2VecPtSnps.insert (make_pair (tokens[0],
					  vector<Snp*> ()));
      mChr2VecPtSnps[tokens[0]].push_back (&(mSnps[tokens[1]]));
    }
    else
    {
      mSnps[tokens[1]].vvIsNa[s].resize (nbSamples, false);
      mSnps[tokens[1]].vvGenos[s].resize (nbSamples,
					  numeric_limits<double>::quiet_NaN());
      for (size_t i = 0; i < nbSamples; ++i)
      {
	AA = atof(tokens[5+3*i].c_str());
	AB = atof(tokens[5+3*i+1].c_str());
	BB = atof(tokens[5+3*i+2].c_str());
	if (AA == 0 && AB == 0 && BB == 0)
	  mSnps[tokens[1]].vvIsNa[s][i] = true;
	else
	{
	  mSnps[tokens[1]].vvGenos[s][i] = 0 * AA + 1 * AB + 2 * BB;
	  maf += mSnps[tokens[1]].vvGenos[s][i];
	}
      }
      maf /= 2 * (nbSamples
		  - count (mSnps[tokens[1]].vvIsNa[s].begin(),
			   mSnps[tokens[1]].vvIsNa[s].end(),
			   true));
      if ((maf <= 0.5 && maf < minMaf) || (maf > 0.5 && 1-maf < minMaf))
	continue;
      ++nbSnpsToKeepPerSubgroup;
      mSnps[tokens[1]].vMafs[s] = maf <= 0.5 ? maf : (1 - maf);
    }
  }
}

void
loadGenosAndSnpInfoFromVcf (
  gzFile & genoStream,
  string & line,
  size_t & nbLines,
  const map<string, string> & mGenoPaths,
  const float & minMaf,
  const vector<string> & vSubgroups,
  const size_t & s,
  const vector<string> & vSnpsToKeep,
  const map<string, vector<Ftr*> > & mChr2VecPtFtrs,
  size_t & nbSnpsToKeepPerSubgroup,
  map<string, Snp> & mSnps,
  map<string, vector<Snp*> > & mChr2VecPtSnps)
{
  vector<string> tokens, tokens2, tokens3;
  size_t nbSamples;
  double maf;
  
  while (getline (genoStream, line))
  {
    ++nbLines;
    if (line.find("#CHROM") == string::npos)
      continue;
    split (line, " \t", tokens);
    nbSamples = tokens.size() - 9;
    break;
  }
  
  while (getline (genoStream, line))
  {
    ++nbLines;
    split (line, " \t", tokens);
    if (tokens.size() != nbSamples + 9)
    {
      cerr << "ERROR: not enough columns on line " << nbLines
	   << " of file " << mGenoPaths.find(vSubgroups[s])->second
	   << " (" << tokens.size() << " != " << nbSamples + 9 << ")"
	   << endl;
      exit (1);
    }
    if (tokens[8].find("GT") == string::npos)
    {
      cerr << "ERROR: missing GT in 9-th field on line " << nbLines
	   << " of file " << mGenoPaths.find(vSubgroups[s])->second
	   << endl;
      exit (1);
    }
    if (mChr2VecPtFtrs.find (tokens[0]) == mChr2VecPtFtrs.end())
      continue;
    if (! vSnpsToKeep.empty()
	&& find (vSnpsToKeep.begin(), vSnpsToKeep.end(), tokens[2])
	== vSnpsToKeep.end())
      continue;
    
    maf = 0;
    if (mSnps.find(tokens[2]) == mSnps.end())
    {
      Snp iSnp;
      Snp_init (iSnp, tokens[2], mGenoPaths.size());
      iSnp.vvIsNa[s].resize (nbSamples, false);
      iSnp.vvGenos[s].resize (nbSamples,
			      numeric_limits<double>::quiet_NaN());
      for (size_t i = 0; i < nbSamples; ++i)
      {
	split (tokens[9+i], ":", tokens2); // if several subfields
	if (tokens2[0].compare(".") == 0)
	  iSnp.vvIsNa[s][i] = true;
	else
	{
	  split (tokens2[0], "|/", tokens3);
	  iSnp.vvGenos[s][i] = 0;
	  if (tokens3[0].compare("1") == 0)
	    iSnp.vvGenos[s][i] += 1;
	  if (tokens3[1].compare("1") == 0)
	    iSnp.vvGenos[s][i] += 1;
	  maf += iSnp.vvGenos[s][i];
	}
      }
      maf /= 2 * (nbSamples
		  - count (iSnp.vvIsNa[s].begin(),
			   iSnp.vvIsNa[s].end(),
			   true));
      if ((maf <= 0.5 && maf < minMaf) || (maf > 0.5 && 1-maf < minMaf))
	continue;
      ++nbSnpsToKeepPerSubgroup;
      iSnp.vMafs[s] = maf <= 0.5 ? maf : (1 - maf);
      iSnp.chr = tokens[0];
      iSnp.coord = atol (tokens[1].c_str());
      mSnps.insert (make_pair (tokens[2], iSnp));
      
      if (mChr2VecPtSnps.find(tokens[0]) == mChr2VecPtSnps.end())
	mChr2VecPtSnps.insert (make_pair (tokens[0],
					  vector<Snp*> ()));
      mChr2VecPtSnps[tokens[0]].push_back (&(mSnps[tokens[2]]));
    }
    else
    {
      mSnps[tokens[2]].vvIsNa[s].resize (nbSamples, false);
      mSnps[tokens[2]].vvGenos[s].resize (nbSamples,
					  numeric_limits<double>::quiet_NaN());
      for (size_t i = 0; i < nbSamples; ++i)
      {
	split (tokens[9+i], ":", tokens2); // if several subfields
	if (tokens2[0].compare(".") == 0)
	  mSnps[tokens[2]].vvIsNa[s][i] = true;
	else
	{
	  split (tokens2[0], "|/", tokens3);
	  mSnps[tokens[2]].vvGenos[s][i] = 0;
	  if (tokens3[0].compare("1") == 0)
	    mSnps[tokens[2]].vvGenos[s][i] += 1;
	  if (tokens3[1].compare("1") == 0)
	    mSnps[tokens[2]].vvGenos[s][i] += 1;
	  maf += mSnps[tokens[2]].vvGenos[s][i];
	}
      }
      maf /= 2 * (nbSamples
		  - count (mSnps[tokens[2]].vvIsNa[s].begin(),
			   mSnps[tokens[2]].vvIsNa[s].end(),
			   true));
      if ((maf <= 0.5 && maf < minMaf) || (maf > 0.5 && 1-maf < minMaf))
	continue;
      ++nbSnpsToKeepPerSubgroup;
      mSnps[tokens[2]].vMafs[s] = maf <= 0.5 ? maf : (1 - maf);
    }
  }
}

void
loadGenosAndSnpInfo (
  const map<string, string> & mGenoPaths,
  const float & minMaf,
  const vector<string> & vSubgroups,
  const vector<string> & vSnpsToKeep,
  const map<string, vector<Ftr*> > & mChr2VecPtFtrs,
  map<string, Snp> & mSnps,
  map<string, vector<Snp*> > & mChr2VecPtSnps,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load genotypes and SNP coordinates ..." << endl << flush;
  
  gzFile genoStream;
  string line;
  size_t nbLines, nbSnpsToKeepPerSubgroup;
  
  // load each file
  for (size_t s = 0; s < vSubgroups.size(); ++s)
  {
    clock_t timeBegin = clock();
    nbSnpsToKeepPerSubgroup = 0;
    nbLines = 0;
    openFile (mGenoPaths.find(vSubgroups[s])->second, genoStream, "rb");
    if (! getline (genoStream, line))
    {
      cerr << "ERROR: problem with the header of file "
	   << mGenoPaths.find(vSubgroups[s])->second << endl;
      exit (1);
    }
    ++nbLines;
    
    if (line.find("##fileformat=VCF") != string::npos) // VCF format
      loadGenosAndSnpInfoFromVcf (genoStream, line, nbLines, mGenoPaths,
				  minMaf, vSubgroups, s, vSnpsToKeep,
				  mChr2VecPtFtrs, nbSnpsToKeepPerSubgroup,
				  mSnps, mChr2VecPtSnps);
    else // IMPUTE format
      loadGenosAndSnpInfoFromImpute (genoStream, line, nbLines, mGenoPaths,
				     minMaf, vSubgroups, s, vSnpsToKeep,
				     mChr2VecPtFtrs, nbSnpsToKeepPerSubgroup,
				     mSnps, mChr2VecPtSnps);
    
    if (! gzeof (genoStream))
    {
      cerr << "ERROR: can't read successfully file "
	   << mGenoPaths.find(vSubgroups[s])->second
	   << " up to the end" << endl;
      exit (1);
    }
    
    closeFile (mGenoPaths.find(vSubgroups[s])->second , genoStream);
    if (verbose > 0)
      cout << "s" << (s+1) << " (" << vSubgroups[s] << "): " << (nbLines-1)
	   << " SNPs (to keep: " << nbSnpsToKeepPerSubgroup
	   << ", loaded in " << fixed << setprecision(2)
	   << (clock() - timeBegin) / (double(CLOCKS_PER_SEC)*60.0)
	   << " min)" << endl << flush;
  }
  
  // sort the SNPs per chr
  for (map<string, vector<Snp*> >::iterator it = mChr2VecPtSnps.begin();
       it != mChr2VecPtSnps.end(); ++it)
    sort (it->second.begin(), it->second.end(), Snp_compByCoord);
  
  if (verbose > 0)
    cout << "nb of SNPs: " << mSnps.size() << endl;
}

/** \brief Load genotypes from each subgroup file.
 *  \note format: row 1 for sample names, column 1 for SNP names,
 *  genotypes as allele dose
 */
void
loadGenos (
  const map<string, string> & mGenoPaths,
  const float & minMaf,
  const vector<string> & vSubgroups,
  const vector<string> & vSnpsToKeep,
  map<string, Snp> & mSnps,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load genotypes ..." << endl << flush;
  
  gzFile genoStream;
  string line;
  vector<string> tokens;
  size_t nbSamples, nbLines;
  double maf;
  
  for (size_t s = 0; s < vSubgroups.size(); ++s)
  {
    clock_t timeBegin = clock();
    nbLines = 0;
    openFile (mGenoPaths.find(vSubgroups[s])->second, genoStream, "rb");
    if (! getline (genoStream, line))
    {
      cerr << "ERROR: problem with the header of file "
	   << mGenoPaths.find(vSubgroups[s])->second << endl;
      exit (1);
    }
    ++nbLines;
    split (line, " \t", tokens);
    if (tokens[0].compare("Id") == 0)
      nbSamples = tokens.size() - 1;
    else
      nbSamples = tokens.size();
    
    while (getline (genoStream, line))
    {
      ++nbLines;
      split (line, " \t", tokens);
      if (tokens.size() != nbSamples + 1)
      {
	cerr << "ERROR: not enough columns on line " << nbLines
	     << " of file " << mGenoPaths.find(vSubgroups[s])->second
	     << " (" << tokens.size() << " != " << nbSamples + 1 << ")"
	     << endl;
	exit (1);
      }
      if (! vSnpsToKeep.empty()
	  && find (vSnpsToKeep.begin(), vSnpsToKeep.end(), tokens[0])
	  == vSnpsToKeep.end())
	continue;
      
      maf = 0;
      if (mSnps.find(tokens[0]) == mSnps.end())
      {
	Snp iSnp;
	Snp_init (iSnp, tokens[0], mGenoPaths.size());
	iSnp.vvIsNa[s].resize (nbSamples, false);
	iSnp.vvGenos[s].resize (nbSamples,
				numeric_limits<double>::quiet_NaN());
	for (size_t i = 1; i < tokens.size(); ++i)
	{
	  if (tokens[i].compare("NA") == 0)
	    iSnp.vvIsNa[s][i-1] = true;
	  else
	  {
	    iSnp.vvGenos[s][i-1] = atof (tokens[i].c_str());
	    maf += iSnp.vvGenos[s][i-1];
	  }
	}
	maf /= 2 * (nbSamples
		    - count (iSnp.vvIsNa[s].begin(),
			     iSnp.vvIsNa[s].end(),
			     true));
	if ((maf <= 0.5 && maf < minMaf) || (maf > 0.5 && 1-maf < minMaf))
	  continue;
	iSnp.vMafs[s] = maf <= 0.5 ? maf : (1 - maf);
	mSnps.insert (make_pair (tokens[0], iSnp));
      }
      else
      {
	mSnps[tokens[0]].vvIsNa[s].resize (nbSamples, false);
	mSnps[tokens[0]].vvGenos[s].resize (nbSamples,
					    numeric_limits<double>::quiet_NaN());
	for (size_t i = 1; i < tokens.size() ; ++i)
	{
	  if (tokens[i].compare("NA") == 0)
	    mSnps[tokens[0]].vvIsNa[s][i-1] = true;
	  else
	  {
	    mSnps[tokens[0]].vvGenos[s][i-1] = atof (tokens[i].c_str());
	    maf += mSnps[tokens[0]].vvGenos[s][i-1];
	  }
	}
	maf /= 2 * (nbSamples
		    - count (mSnps[tokens[0]].vvIsNa[s].begin(),
			     mSnps[tokens[0]].vvIsNa[s].end(),
			     true));
	if ((maf <= 0.5 && maf < minMaf) || (maf > 0.5 && 1-maf < minMaf))
	  continue;
	mSnps[tokens[0]].vMafs[s] = maf <= 0.5 ? maf : (1 - maf);
      }
    }
    if (! gzeof (genoStream))
    {
      cerr << "ERROR: can't read successfully file "
	   << mGenoPaths.find(vSubgroups[s])->second
	   << " up to the end" << endl;
      exit (1);
    }
    
    closeFile (mGenoPaths.find(vSubgroups[s])->second , genoStream);
    if (verbose > 0)
      cout << "s" << (s+1) << " (" << vSubgroups[s] << "): " << (nbLines-1)
	   << " SNPs (loaded in " << fixed << setprecision(2)
	   << (clock() - timeBegin) / (double(CLOCKS_PER_SEC)*60.0)
	   << " min)" << endl << flush;
  }
  
  if (mSnps.size() == 0)
  {
    cerr << "ERROR: no SNP to analyze" << endl;
    exit (1);
  }
  if (verbose > 0)
    cout << "total nb of SNPs with genotypes: " << mSnps.size() << endl;
}

void
loadSnpInfo (
  const string & snpCoordsFile,
  map<string, Snp> & mSnps,
  map<string, vector<Snp*> > & mChr2VecPtSnps,
  const int & verbose)
{
  if (verbose > 0)
    cout << "load SNP coordinates ..." << endl << flush;
  
  // parse the BED file
  gzFile snpCoordsStream;
  openFile (snpCoordsFile, snpCoordsStream, "rb");
  string line;
  vector<string> tokens;
  while (getline (snpCoordsStream, line))
  {
    split (line, " \t", tokens);
    if (mSnps.find(tokens[3]) == mSnps.end())
      continue;
    mSnps[tokens[3]].chr = tokens[0];
    mSnps[tokens[3]].coord = atol (tokens[1].c_str()) + 1;
    
    if (mChr2VecPtSnps.find(tokens[0]) == mChr2VecPtSnps.end())
      mChr2VecPtSnps.insert (make_pair (tokens[0],
					vector<Snp*> ()));
    mChr2VecPtSnps[tokens[0]].push_back (&(mSnps[tokens[3]]));
  }
  if (! gzeof (snpCoordsStream))
  {
    cerr << "ERROR: can't read successfully file "
	 << snpCoordsFile
	 << " up to the end" << endl;
    exit (1);
  }
  closeFile (snpCoordsFile, snpCoordsStream);
  
  // check that all SNPs have coordinates
  map<string, Snp>::iterator it = mSnps.begin();
  while (it != mSnps.end())
  {
    if (it->second.chr.empty())
    {
      cerr << "WARNING: skip SNP " << it->second.name
	   << " because it has no coordinate" << endl;
      mSnps.erase (it++);
    }
    else
      ++it;
  }
  
  // sort the SNPs per chr
  for (map<string, vector<Snp*> >::iterator it = mChr2VecPtSnps.begin();
       it != mChr2VecPtSnps.end(); ++it)
    sort (it->second.begin(), it->second.end(), Snp_compByCoord);
  
  if (mSnps.size() == 0)
  {
    cerr << "ERROR: no SNP to analyze" << endl;
    exit (1);
  }
  if (verbose > 0)
    cout << "total nb of SNPs to analyze: " << mSnps.size() << endl;
}

void
loadListCovarFiles (
  const string & covarPathsFile,
  const string & sbgrpToKeep,
  const vector<string> & vSubgroups,
  map<string, string> & mCovarPaths,
  const int & verbose)
{
  mCovarPaths = loadTwoColumnFile (covarPathsFile, verbose);
  
  map<string, string>::iterator it = mCovarPaths.begin();
  while (it != mCovarPaths.end())
  {
    if (! sbgrpToKeep.empty() && sbgrpToKeep.compare(it->first) != 0)
    {
      mCovarPaths.erase (it++);
      continue;
    }
    if (find (vSubgroups.begin(), vSubgroups.end(), it->first)
	== vSubgroups.end())
    {
      cerr << "WARNING: skip covariates of subgroup " << it->first
	   << " as there is no corresponding genotype nor phenotype files"
	   << endl;
      mCovarPaths.erase (it++);
    }
    else
      ++it;
  }
}

void
loadCovarsFromFiles (
  const vector<string> & vSubgroups,
  const vector<string> & vSamples,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const map<string, string> & mCovarPaths,
  vector<map<string, vector<double> > > & vSbgrp2Covars)
{
  gzFile covarStream;
  vector<string> tokens;
  string line;
  size_t nbLines = 0;

  for (size_t s = 0; s < vSubgroups.size(); ++s)
  {
    if (mCovarPaths.find(vSubgroups[s]) == mCovarPaths.end())
      continue;
    
    nbLines = 0;
    openFile (mCovarPaths.find(vSubgroups[s])->second, covarStream, "rb");
    
    // parse the header line to get covar names and order
    if (! getline (covarStream, line))
    {
      cerr << "ERROR: problem with the header of file "
	   << mCovarPaths.find(vSubgroups[s])->second << endl;
      exit (1);
    }
    if (line.empty())
    {
      cerr << "ERROR: file " << mCovarPaths.find(vSubgroups[s])->second
	   << " is empty" << endl;
      exit (1);
    }
    ++nbLines;
    vector<string> vCovars;
    split (line, " \t", vCovars);
    if (vCovars[0].compare ("sample") == 0)
      vCovars.erase (vCovars.begin());
    else
    {
      cerr << "ERROR: file " << mCovarPaths.find(vSubgroups[s])->second
	   << " should have a header line starting with 'sample'" << endl;
      exit (1);
    }
    
    // parse the rest of the file into a temporary container
    map<string, vector<double> > mSample2Covars;
    while (getline (covarStream, line))
    {
      ++nbLines;
      split (line, " \t", tokens);
      if (tokens.size() != vCovars.size() + 1)
      {
	cerr << "ERROR: not enough columns on line " << nbLines
	     << " of file " << mCovarPaths.find(vSubgroups[s])->second
	     << " (" << tokens.size() << " != " << vCovars.size() + 1 << ")"
	     << endl;
	exit (1);
      }
      if (find (vSamples.begin(), vSamples.end(), tokens[0])
	  == vSamples.end())
      {
	cerr << "WARNING: skip sample " << tokens[0]
	     << " in line " << nbLines
	     << " of file " << mCovarPaths.find(vSubgroups[s])->second
	     << " because it is absent from all genotype and phenotype files"
	     << endl;
	continue;
      }
      mSample2Covars.insert (make_pair (tokens[0], vector<double> (
					  vCovars.size(),
					  numeric_limits<double>::quiet_NaN())));
      for (size_t c = 0; c < vCovars.size(); ++c)
	mSample2Covars[tokens[0]][c] = atof (tokens[(c+1)].c_str());
    }
    if (! gzeof (covarStream))
    {
      cerr << "ERROR: can't read successfully file "
	   << mCovarPaths.find(vSubgroups[s])->second
	   << " up to the end" << endl;
      exit (1);
    }
    closeFile (mCovarPaths.find(vSubgroups[s])->second, covarStream);
    
    // check that all samples with covar from the given subgroup also have
    // a genotype and a phenotype
    if (vvSampleIdxGenos.size() == 1)
    {
      for (size_t i = 0; i < vSamples.size(); ++i)
	if (vvSampleIdxGenos[s][i] != string::npos
	    && vvSampleIdxPhenos[s][i] != string::npos
	    && mSample2Covars.find(vSamples[i]) == mSample2Covars.end())
	{
	  cerr << "ERROR: sample " << vSamples[i] << " has genotype"
	       << " and phenotype in subgroup " << (s+1)
	       << " but no covariate" << endl;
	  exit (1);
	}
    }
    else
    {
      for (size_t i = 0; i < vSamples.size(); ++i)
	if (vvSampleIdxGenos[s][i] != string::npos
	    && vvSampleIdxPhenos[s][i] != string::npos
	    && mSample2Covars.find(vSamples[i]) == mSample2Covars.end())
	{
	  cerr << "ERROR: sample " << vSamples[i] << " has genotype"
	       << " and phenotype in subgroup " << (s+1)
	       << " but no covariate" << endl;
	  exit (1);
	}
    }
    
    // fill the final containers with samples in the right order
    map<string, vector<double> > mCovars;
    for (size_t c = 0; c < vCovars.size(); ++c)
    {
      mCovars.insert (make_pair (vCovars[c], vector<double> (
				   vSamples.size(),
				   numeric_limits<double>::quiet_NaN())));
      for (size_t i = 0; i < vSamples.size(); ++i)
	if (mSample2Covars.find(vSamples[i]) != mSample2Covars.end())
	  mCovars[vCovars[c]][i] = mSample2Covars[vSamples[i]][c];
    }
    vSbgrp2Covars.push_back (mCovars);
  }
}

void
loadCovariates (
  const string & covarPathsFile,
  const string & sbgrpToKeep,
  const vector<string> & vSubgroups,
  const vector<string> & vSamples,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  vector<map<string, vector<double> > > & vSbgrp2Covars,
  const int & verbose)
{
  if (covarPathsFile.empty())
    vSbgrp2Covars = vector<map<string, vector<double> > >
      (vvSampleIdxPhenos.size(), map<string, vector<double> > ());
  else
  {
    if (verbose > 0)
      cout << "load covariates ..." << endl << flush;
    map<string, string> mCovarPaths;
    loadListCovarFiles (covarPathsFile, sbgrpToKeep, vSubgroups, mCovarPaths,
			verbose);
    loadCovarsFromFiles (vSubgroups, vSamples, vvSampleIdxGenos,
			 vvSampleIdxPhenos, mCovarPaths, vSbgrp2Covars);
    if (verbose > 0)
      for (size_t s = 0; s < vSbgrp2Covars.size(); ++s)
	cout << "s" << (s+1) << " (" << vSubgroups[s] << "): "
	     << vSbgrp2Covars[s].size() << " covariates" << endl;
  }
}

/** \brief Infer associations between each feature and its SNPs in cis
 *  via a loop over all features
 */
void
inferAssos (
  map<string, Ftr> & mFtrs,
  const map<string, vector<Snp*> > & mChr2VecPtSnps,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const string & anchor,
  const size_t & lenCis,
  const int & whichStep,
  const bool & needQnorm,
  const vector<map<string, vector<double> > > & vSbgrp2Covars,
  const Grid & iGridL,
  const Grid & iGridS,
  const string & whichBfs,
  const bool & mvlr,
  const float & propFitSigma,
  const int & verbose)
{
  if (verbose > 0)
  {
    cout << "look for association between each pair feature-SNP ..." << endl
	 << "anchor=" << anchor << " lenCis=" << lenCis;
    if (whichStep == 3 || whichStep == 4 || whichStep == 5)
      cout << " multivariate=" << boolalpha << mvlr
	   << " propFitSigma=" << propFitSigma;
    cout << endl << flush;
  }
  
  clock_t timeBegin = clock();
  size_t nbAnalyzedPairs = 0;
  size_t countFtrs = 0;
  
  for (map<string, Ftr>::iterator itF = mFtrs.begin();
       itF != mFtrs.end(); ++itF)
  {
    ++countFtrs;
    if (verbose == 1)
      progressBar ("", countFtrs, mFtrs.size());
    Ftr_getCisSnps (itF->second, mChr2VecPtSnps, anchor, lenCis);
    if (verbose > 1)
      cout << itF->second.name << ": " << itF->second.vPtCisSnps.size()
	   << " SNPs in cis" << endl << flush;
    if (itF->second.vPtCisSnps.size() > 0)
    {
      Ftr_inferAssos (itF->second, vvSampleIdxPhenos, vvSampleIdxGenos,
		      whichStep, needQnorm, vSbgrp2Covars, iGridL, iGridS,
		      whichBfs, mvlr, propFitSigma, verbose-1);
      nbAnalyzedPairs += itF->second.vResFtrSnps.size();
    }
  }
  
  if (verbose > 0)
  {
    if(verbose == 1)
      cout << " (" << fixed << setprecision(2) << (clock() - timeBegin) /
	(double(CLOCKS_PER_SEC)*60.0) << " min)" << endl << flush;
    cout << "nb of analyzed feature-SNP pairs: " << nbAnalyzedPairs << endl;
  }
}

/** \brief Make permutations for the separate analysis
 *  via a loop over all features
 */
void
makePermsSep (
  map<string, Ftr> & mFtrs,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const bool & needQnorm,
  const vector<map<string, vector<double> > > & vSbgrp2Covars,
  const size_t & nbPerms,
  const size_t & seed,
  const int & trick,
  const int & whichPermSep,
  const gsl_rng * rngPerm,
  const gsl_rng * rngTrick,
  const int & verbose)
{
  if (whichPermSep == 1)
  {
    clock_t timeBegin = clock();
    gsl_rng_set (rngPerm, seed);
    if (trick != 0)
      gsl_rng_set (rngTrick, seed);
    size_t countFtrs = 0;
    for (map<string, Ftr>::iterator itF = mFtrs.begin();
	 itF != mFtrs.end(); ++itF)
    {
      ++countFtrs;
      if (verbose == 1)
	progressBar ("sep", countFtrs, mFtrs.size());
      if (itF->second.vResFtrSnps.size() > 0)
	Ftr_makePermsSepAllSubgrps (itF->second, vvSampleIdxPhenos,
				    vvSampleIdxGenos, needQnorm, vSbgrp2Covars,
				    nbPerms, trick, rngPerm, rngTrick);
    }
    if (verbose == 1)
      cout << " (" << fixed << setprecision(2)<< (clock() - timeBegin) /
	(double(CLOCKS_PER_SEC)*60.0) << " min)" << endl << flush;
  }
  else
  {
    size_t nbSubgroups = vvSampleIdxPhenos.size();
    stringstream ss;
    for (size_t s = 0; s < nbSubgroups; ++s)
    {
      clock_t timeBegin = clock();
      gsl_rng_set (rngPerm, seed);
      if (trick != 0)
	gsl_rng_set (rngTrick, seed);
      ss.str("");
      ss << "s" << (s+1);
      size_t countFtrs = 0;
      for (map<string, Ftr>::iterator itF = mFtrs.begin();
	   itF != mFtrs.end(); ++itF)
      {
	++countFtrs;
	if (verbose == 1)
	  progressBar (ss.str(), countFtrs, mFtrs.size());
	if (itF->second.vPtCisSnps.size() > 0)
	  Ftr_makePermsSepOneSubgrp (itF->second, vvSampleIdxPhenos,
				     vvSampleIdxGenos, needQnorm, vSbgrp2Covars,
				     nbPerms, trick, s, rngPerm, rngTrick);
      }
      if (verbose == 1)
	cout << " (" << fixed << setprecision(2) << (clock() - timeBegin) /
	  (double(CLOCKS_PER_SEC)*60.0) << " min)" << endl << flush;
    }
  }
}

/** \brief Make permutations for the joint analysis
 *  via a loop over all features
 */
void
makePermsJoint (
  map<string, Ftr> & mFtrs,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const vector<vector<size_t> > & vvSampleIdxGenos,
  const bool & needQnorm,
  const vector<map<string, vector<double> > > & vSbgrp2Covars,
  const Grid & iGridL,
  const Grid & iGridS,
  const bool & mvlr,
  const float & propFitSigma,
  const size_t & nbPerms,
  const size_t & seed,
  const int & trick,
  const string & whichPermBf,
  const bool & useMaxBfOverSnps,
  const gsl_rng * rngPerm,
  const gsl_rng * rngTrick,
  const int & verbose)
{
  clock_t timeBegin = clock();
  gsl_rng_set (rngPerm, seed);
  if (trick != 0)
    gsl_rng_set (rngTrick, seed);
  size_t countFtrs = 0;
  for (map<string, Ftr>::iterator itF = mFtrs.begin();
       itF != mFtrs.end(); ++itF)
  {
    ++countFtrs;
    if (verbose == 1)
      progressBar ("joint", countFtrs, mFtrs.size());
    if (itF->second.vResFtrSnps.size() > 0)
      Ftr_makePermsJoint (itF->second, vvSampleIdxPhenos,
			  vvSampleIdxGenos, needQnorm, vSbgrp2Covars, iGridL,
			  iGridS, mvlr, propFitSigma, nbPerms, trick, whichPermBf,
			  useMaxBfOverSnps, rngPerm, rngTrick);
  }
  if (verbose == 1)
    cout << " (" << fixed << setprecision(2) << (clock() - timeBegin) /
      (double(CLOCKS_PER_SEC)*60.0) << " min)" << endl << flush;
}

/** \brief Make permutations for the separate and/or the joint analysis
 *  depending on the --step option
 */
void
makePerms (
  map<string, Ftr> & mFtrs,
  const vector<vector<size_t> > & vvSampleIdxPhenos,
  const int & whichStep,
  const bool & needQnorm,
  const vector<map<string, vector<double> > > & vSbgrp2Covars,
  const Grid & iGridL,
  const Grid & iGridS,
  const bool & mvlr,
  const float & propFitSigma,
  const size_t & nbPerms,
  const size_t & seed,
  const int & trick,
  const int & whichPermSep,
  const string & whichPermBf,
  const bool & useMaxBfOverSnps,
  const int & verbose)
{
  if (verbose > 0)
  {
    cout << "get feature-level P-values by permuting phenotypes ..." << endl
	 << "permutation"<< (nbPerms > 1 ? "s=" : "=") << nbPerms
	 << ", seed=" << seed
	 << ", trick=" << trick;
    if (whichStep == 2 || whichStep == 5)
      cout << ", permSep=" << whichPermSep;
    if (whichStep == 4 || whichStep == 5)
      cout << ", permBf=" << whichPermBf;
    cout << endl << flush;
  }
  
  gsl_rng * rngPerm = NULL, * rngTrick = NULL;
  gsl_rng_env_setup();
  rngPerm = gsl_rng_alloc (gsl_rng_default);
  if (rngPerm == NULL)
  {
    cerr << "ERROR: can't allocate memory for the RNG" << endl;
    exit (1);
  }
  if (trick != 0)
  {
    rngTrick = gsl_rng_alloc (gsl_rng_default);
    if (rngTrick == NULL)
    {
      cerr << "ERROR: can't allocate memory for the RNG" << endl;
      exit (1);
    }
  }
  
  if (whichStep == 2 || (whichStep == 5 && ! mvlr))
    makePermsSep (mFtrs, vvSampleIdxPhenos, vvSampleIdxPhenos, needQnorm,
		  vSbgrp2Covars, nbPerms, seed, trick, whichPermSep,rngPerm,
		  rngTrick, verbose);
  
  if (whichStep == 4 || whichStep == 5)
    makePermsJoint (mFtrs, vvSampleIdxPhenos, vvSampleIdxPhenos,
		    needQnorm, vSbgrp2Covars, iGridL, iGridS, mvlr, propFitSigma,
		    nbPerms, seed, trick, whichPermBf, useMaxBfOverSnps,
		    rngPerm, rngTrick, verbose);
  
  gsl_rng_free (rngPerm);
  if (trick != 0)
    gsl_rng_free (rngTrick);
}

void
writeResSstats (
  const string & outPrefix,
  const map<string, Ftr> & mFtrs,
  const map<string, Snp> & mSnps,
  const vector<string> & vSubgroups,
  const vector<map<string, vector<double> > > & vSbgrp2Covars,
  const int & verbose)
{
  if (verbose > 0)
    cout << "write results of summary statistics in each subgroup ..."
	 << endl << flush;
  
  for (size_t s = 0; s < vSubgroups.size(); ++s)
  {
    stringstream ssOutFile, ssTxt;
    ssOutFile << outPrefix << "_sumstats_" << vSubgroups[s] << ".txt.gz";
    if (verbose > 0)
      cout << "file " << ssOutFile.str() << endl << flush;
    gzFile outStream;
    openFile (ssOutFile.str(), outStream, "wb");
    
    ssTxt << "ftr snp maf n pve sigmahat"
	  << " betahat.geno sebetahat.geno betapval.geno";
    for (map<string, vector<double> >::const_iterator it =
	   vSbgrp2Covars[s].begin(); it != vSbgrp2Covars[s].end(); ++it)
      ssTxt << " betahat." << it->first
	    << " sebetahat." << it->first
	    << " betapval." << it->first;
    ssTxt << endl;
    size_t lineId = 1;
    gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
    
    ssTxt.precision (7);
    ssTxt.setf (ios::scientific);
    for (map<string, Ftr>::const_iterator itF = mFtrs.begin();
	 itF != mFtrs.end(); ++itF)
    {
      const Ftr * ptF = &(itF->second);
      for (vector<ResFtrSnp>::const_iterator itP = ptF->vResFtrSnps.begin();
	   itP != ptF->vResFtrSnps.end(); ++itP)
	if (itP->vNs[s] > 2)
	{
	  ssTxt.str("");
	  ++lineId;
	  ssTxt << ptF->name
		<< " " << itP->snp
		<< " " << mSnps.find(itP->snp)->second.vMafs[s]
		<< " " << itP->vNs[s]
		<< " " << itP->vPves[s]
		<< " " << itP->vSigmahats[s]
		<< " " << (itP->vMapPredictors[s].find("genotype")->second)[0]
		<< " " << (itP->vMapPredictors[s].find("genotype")->second)[1]
		<< " " << (itP->vMapPredictors[s].find("genotype")->second)[2];
	  for (map<string, vector<double> >::const_iterator itC =
		 vSbgrp2Covars[s].begin(); itC != vSbgrp2Covars[s].end();
	       ++itC)
	    ssTxt << " " << (itP->vMapPredictors[s].find(itC->first)->
			     second)[0]
		  << " " << (itP->vMapPredictors[s].find(itC->first)->
			     second)[1]
		  << " " << (itP->vMapPredictors[s].find(itC->first)->
			     second)[2];
	  ssTxt << endl;
	  gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
	}
    }
    
    closeFile (ssOutFile.str(), outStream);
  }
}

void
writeResSepPermPval (
  const string & outPrefix,
  const map<string, Ftr> & mFtrs,
  const vector<string> & vSubgroups,
  const size_t & seed,
  const int & verbose)
{
  if (verbose > 0)
    cout << "write results of feature-level P-values in each subgroup separately ..."
	 << endl << flush;
  
  for (size_t s = 0; s < vSubgroups.size(); ++s)
  {
    stringstream ssOutFile, ssTxt;
    ssOutFile << outPrefix << "_permPval_" << vSubgroups[s] << ".txt.gz";
    if (verbose > 0)
      cout << "file " << ssOutFile.str() << endl << flush;
    gzFile outStream;
    openFile (ssOutFile.str(), outStream, "wb");
    size_t lineId = 0;
    
    ssTxt << "# seed=" << seed << endl;
    ++lineId;
    gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
    
    ssTxt.str("");
    ssTxt << "ftr nbSnps permPval nbPerms minTruePval" << endl;
    ++lineId;
    gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
    
    for (map<string, Ftr>::const_iterator itF = mFtrs.begin();
	 itF != mFtrs.end(); ++itF)
    {
      const Ftr * ptF = &(itF->second);
      size_t nbCisSnps = Ftr_getNbCisSnpsInGivenSubgroup (*ptF, s);
      if (nbCisSnps > 0)
      {
	ssTxt.str("");
	++lineId;
	ssTxt << ptF->name
	      << " " << nbCisSnps
	      << " " << ptF->vPermPvalsSep[s]
	      << " " << ptF->vNbPermsSoFar[s]
	      << " " << ptF->vMinTruePvals[s]
	      << endl;
	gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
      }
    }
    
    closeFile (ssOutFile.str(), outStream);
  }
}

void
writeResSepPermPval (
  const string & outPrefix,
  const map<string, Ftr> & mFtrs,
  const size_t & seed,
  const int & verbose)
{
  if (verbose > 0)
    cout << "write results of feature-level P-values, each subgroup ..."
	 << endl << flush;
  
  stringstream ssOutFile, ssTxt;
  ssOutFile << outPrefix << "_sepPermPvals.txt.gz";
  if (verbose > 0)
    cout << "file " << ssOutFile.str() << endl << flush;
  gzFile outStream;
  openFile (ssOutFile.str(), outStream, "wb");
  size_t lineId = 0;
  
  ssTxt << "# seed=" << seed << endl;
  ++lineId;
  gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
  
  ssTxt.str("");
  ssTxt << "ftr nbSnps sepPermPval nbPerms minTruePval" << endl;
  ++lineId;
  gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
  
  for (map<string, Ftr>::const_iterator itF = mFtrs.begin();
       itF != mFtrs.end(); ++itF)
  {
    if (itF->second.vPtCisSnps.size() > 0)
    {
      ssTxt.str("");
      ++lineId;
      ssTxt << itF->second.name
	    << " " << itF->second.vPtCisSnps.size()
	    << " " << itF->second.sepPermPval
	    << " " << itF->second.nbPermsSoFarSep
	    << " " << itF->second.minTruePval
	    << endl;
      gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
    }
  }
  
  closeFile (ssOutFile.str(), outStream);
}

/** \brief 
 *  \note the 'const' BFs use the large grid while the other use the small one,
 *  but the header line lists the large grid only
 */
void
writeResAbfsRaw (
  const string & outPrefix,
  const map<string, Ftr> & mFtrs,
  const size_t & nbSubgroups,
  const Grid & iGridL,
  const Grid & iGridS,
  const string & whichBfs,
  const int & verbose)
{
  if (verbose > 0)
    cout << "write results of Bayes Factors (one per grid value) ..."
	 << endl << flush;
  
  gsl_combination * comb;
  stringstream ssOutFile, ssConfig, ssTxt;
  ssOutFile << outPrefix << "_l10abfs_raw.txt.gz";
  if (verbose > 0)
    cout << "file " << ssOutFile.str() << endl << flush;
  gzFile outStream;
  openFile (ssOutFile.str(), outStream, "wb");
  
  // write header line
  ssTxt << "ftr snp config";
  for (size_t i = 0; i < iGridL.size(); ++i)
    ssTxt << " l10abf.grid" << (i+1);
  ssTxt << endl;
  size_t lineId = 1;
  gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
  
  // write results
  ssTxt.precision (7);
  ssTxt.setf (ios::scientific);
  for (map<string, Ftr>::const_iterator itF = mFtrs.begin();
       itF != mFtrs.end(); ++itF)
  {
    const Ftr * ptF = &(itF->second);
    
    for (vector<ResFtrSnp>::const_iterator itP = ptF->vResFtrSnps.begin();
	 itP != ptF->vResFtrSnps.end(); ++itP)
    {
      // write gen BFs (large grid)
      ssTxt.str("");
      ssTxt << ptF->name
	    << " " << itP->snp
	    << " gen";
      for (size_t i = 0; i < iGridL.size(); ++i)
	ssTxt << " " << itP->mUnweightedAbfs.find("gen")->second[i];
      ssTxt << endl;
      ++lineId;
      gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
      
      // write gen-fix BFs (large grid)
      ssTxt.str("");
      ssTxt << ptF->name
	    << " " << itP->snp
	    << " gen-fix";
      for (size_t i = 0; i < iGridL.size(); ++i)
	ssTxt << " " << itP->mUnweightedAbfs.find("gen-fix")->second[i];
      ssTxt << endl;
      ++lineId;
      gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
      
      // write gen-maxh BFs (large grid)
      ssTxt.str("");
      ssTxt << ptF->name
	    << " " << itP->snp
	    << " gen-maxh";
      for (size_t i = 0; i < iGridL.size(); ++i)
	ssTxt << " " << itP->mUnweightedAbfs.find("gen-maxh")->second[i];
      ssTxt << endl;
      ++lineId;
      gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
      
      if (whichBfs.compare("gen") != 0)
      {
	// write the BFs for each config (small grid)
	for (size_t k = 1; k <= nbSubgroups; ++k)
	{
	  comb = gsl_combination_calloc (nbSubgroups, k);
	  if (comb == NULL)
	  {
	    cerr << "ERROR: can't allocate memory for the combination"
		 << endl;
	    exit (1);
	  }
	  while (true)
	  {
	    ssTxt.str("");
	    ssTxt << ptF->name
		  << " " << itP->snp;
	    ssConfig.str("");
	    ssConfig << gsl_combination_get (comb, 0) + 1;
	    if (comb->k > 1)
	      for (size_t i = 1; i < k; ++i)
		ssConfig << "-" << gsl_combination_get (comb, i) + 1;
	    ssTxt << " " << ssConfig.str();
	    for (size_t i = 0; i < iGridL.size(); ++i)
	    {
	      if (i < iGridS.size()
		  && itP->mUnweightedAbfs.find(ssConfig.str()) !=
		  itP->mUnweightedAbfs.end()) // for ftr not present in all subgroups
		ssTxt << " " << itP->mUnweightedAbfs.find(ssConfig.str())->
		  second[i];
	      else
		ssTxt << " " << numeric_limits<double>::quiet_NaN();
	    }
	    ssTxt << endl;
	    ++lineId;
	    gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
	    if (gsl_combination_next (comb) != GSL_SUCCESS)
	      break;
	  }
	  if (whichBfs.compare("sin") == 0)
	    break;
	  gsl_combination_free (comb);
	}
      }
    }
  }
  
  closeFile (ssOutFile.str(), outStream);
}

void
writeResAbfsAvgGrids (
  const string & outPrefix,
  const map<string, Ftr> & mFtrs,
  const size_t & nbSubgroups,
  const string & whichBfs,
  const int & verbose)
{
  if (verbose > 0)
    cout << "write results of Bayes Factors (one per configuration) ..."
	 << endl << flush;
  
  gsl_combination * comb;
  stringstream ssOutFile, ssConfig, ssTxt;
  ssOutFile << outPrefix << "_l10abfs_avg-grids.txt.gz";
  if (verbose > 0)
    cout << "file " << ssOutFile.str() << endl << flush;
  gzFile outStream;
  openFile (ssOutFile.str(), outStream, "wb");
  
  // write header line
  ssTxt << "ftr snp nb.subgroups nb.samples"
	<< " l10abf.gen l10abf.gen.fix l10abf.gen.maxh";
  if (whichBfs.compare("sin") == 0 || whichBfs.compare("all") == 0)
    ssTxt << " l10abf.sin l10abf.gen.sin";
  if (whichBfs.compare("all") == 0)
    ssTxt << " l10abf.all";
  if (whichBfs.compare("gen") != 0)
  {
    for (size_t k = 1; k <= nbSubgroups; ++k)
    {
      comb = gsl_combination_calloc (nbSubgroups, k);
      if (comb == NULL)
      {
	cerr << "ERROR: can't allocate memory for the combination" << endl;
	exit (1);
      }
      while (true)
      {
	ssTxt << " l10abf." << gsl_combination_get (comb, 0) + 1;
	if (comb->k > 1)
	  for (size_t i = 1; i < k; ++i)
	    ssTxt << "-" << gsl_combination_get (comb, i) + 1;
	if (gsl_combination_next (comb) != GSL_SUCCESS)
	  break;
      }
      gsl_combination_free (comb);
      if (whichBfs.compare("sin") == 0)
	break;
    }
  }
  ssTxt << endl;
  size_t lineId = 1;
  gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
  
  // write results
  size_t n;
  ssTxt.precision (7);
  ssTxt.setf (ios::scientific);
  for (map<string, Ftr>::const_iterator itF = mFtrs.begin();
       itF != mFtrs.end(); ++itF)
  {
    const Ftr * ptF = &(itF->second);
    for (vector<ResFtrSnp>::const_iterator itP = ptF->vResFtrSnps.begin();
	 itP != ptF->vResFtrSnps.end(); ++itP)
    {
      n = accumulate (itP->vNs.begin(), itP->vNs.end(), 0);
      if (n > 2)
      {
	ssTxt.str("");
	ssTxt << ptF->name
	      << " " << itP->snp
	      << " " << count_if (itP->vNs.begin(), itP->vNs.end(), isNonZero)
	      << " " << n
	      << " " << itP->mWeightedAbfs.find("gen")->second
	      << " " << itP->mWeightedAbfs.find("gen-fix")->second
	      << " " << itP->mWeightedAbfs.find("gen-maxh")->second;
	if (whichBfs.compare("sin") == 0 || whichBfs.compare("all") == 0)
	  ssTxt << " " << itP->mWeightedAbfs.find("sin")->second
		<< " " << itP->mWeightedAbfs.find("gen-sin")->second;
	if (whichBfs.compare("all") == 0)
	  ssTxt << " " << itP->mWeightedAbfs.find("all")->second;
	if (whichBfs.compare("gen") != 0)
	{
	  for (size_t k = 1; k <= nbSubgroups; ++k)
	  {
	    comb = gsl_combination_calloc (nbSubgroups, k);
	    if (comb == NULL)
	    {
	      cerr << "ERROR: can't allocate memory for the combination"
		   << endl;
	      exit (1);
	    }
	    while (true)
	    {
	      ssConfig.str("");
	      ssConfig << gsl_combination_get (comb, 0) + 1;
	      if (comb->k > 1)
		for (size_t i = 1; i < k; ++i)
		  ssConfig << "-" << gsl_combination_get (comb, i) + 1;
	      ssTxt << " " << itP->mWeightedAbfs.find(ssConfig.str())->
		second;
	      if (gsl_combination_next (comb) != GSL_SUCCESS)
		break;
	    }
	    gsl_combination_free (comb);
	    if (whichBfs.compare("sin") == 0)
	      break;
	  }
	}
	ssTxt << endl;
	++lineId;
	gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
      }
    }
  }
  
  closeFile (ssOutFile.str(), outStream);
}

void
writeResJointPermPval (
  const string & outPrefix,
  const map<string, Ftr> & mFtrs,
  const size_t & seed,
  const string & whichPermBf,
  const int & verbose)
{
  if (verbose > 0)
    cout << "write results of feature-level P-values, all subgroups jointly ..."
	 << endl << flush;
  
  stringstream ssOutFile, ssTxt;
  ssOutFile << outPrefix << "_jointPermPvals.txt.gz";
  if (verbose > 0)
    cout << "file " << ssOutFile.str() << endl << flush;
  gzFile outStream;
  openFile (ssOutFile.str(), outStream, "wb");
  size_t lineId = 0;
  
  ssTxt << "# permBf=" << whichPermBf << " seed=" << seed << endl;
  ++lineId;
  gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
  
  ssTxt.str("");
  ssTxt << "ftr nbSnps jointPermPval nbPerms maxL10TrueAbf avgL10TrueAbf" << endl;
  ++lineId;
  gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
  
  for (map<string, Ftr>::const_iterator itF = mFtrs.begin();
       itF != mFtrs.end(); ++itF)
  {
    if (itF->second.vPtCisSnps.size() > 0)
    {
      ssTxt.str("");
      ++lineId;
      ssTxt << itF->second.name
	    << " " << itF->second.vPtCisSnps.size()
	    << " " << itF->second.jointPermPval
	    << " " << itF->second.nbPermsSoFarJoint
	    << " " << itF->second.maxL10TrueAbf
	    << " " << itF->second.avgL10TrueAbf
	    << endl;
      gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), lineId);
    }
  }
  
  closeFile (ssOutFile.str(), outStream);
}

void
writeRes (
  const string & outPrefix,
  const bool & outRaw,
  const map<string, Ftr> & mFtrs,
  const map<string, Snp> & mSnps,
  const vector<string> & vSubgroups,
  const vector<map<string, vector<double> > > & vSbgrp2Covars,
  const int & whichStep,
  const Grid & iGridL,
  const Grid & iGridS,
  const string & whichBfs,
  const int & whichPermSep,
  const bool & mvlr,
  const size_t & seed,
  const string & whichPermBf,
  const int & verbose)
{
  if (whichStep == 1 || whichStep == 2 ||
      (! mvlr && (whichStep == 3 || whichStep == 4 || whichStep == 5)))
    writeResSstats (outPrefix, mFtrs, mSnps, vSubgroups, vSbgrp2Covars,
		    verbose);
  
  if (whichStep == 2 || (whichStep == 5 && ! mvlr))
  {
    if (whichPermSep == 1)
      writeResSepPermPval (outPrefix, mFtrs, seed, verbose);
    else
      writeResSepPermPval (outPrefix, mFtrs, vSubgroups, seed, verbose);;
  }
  
  if (whichStep == 3 || whichStep == 4 || whichStep == 5)
  {
    if (outRaw)
      writeResAbfsRaw (outPrefix, mFtrs, vSubgroups.size(), iGridL, iGridS,
		       whichBfs, verbose);
    writeResAbfsAvgGrids (outPrefix, mFtrs, vSubgroups.size(), whichBfs,
			  verbose);
  }
  
  if (whichStep == 4 || whichStep == 5)
    writeResJointPermPval (outPrefix, mFtrs, seed, whichPermBf, verbose);
}

void
run (
  const string & genoPathsFile,
  const string & snpCoordFile,
  const string & phenoPathsFile,
  const string & ftrCoordsFile,
  const string & anchor,
  const size_t & lenCis,
  const string & outPrefix,
  const bool & outRaw,
  const int & whichStep,
  const bool & needQnorm,
  const float & minMaf,
  const string & covarPathsFile,
  const string & largeGridFile,
  const string & smallGridFile,
  const string & whichBfs,
  const bool & mvlr,
  const float & propFitSigma,
  const size_t & nbPerms,
  const size_t & seed,
  const int & trick,
  const int & whichPermSep,
  const string & whichPermBf,
  const bool & useMaxBfOverSnps,
  const string & ftrsToKeepFile,
  const string & snpsToKeepFile,
  const string & sbgrpToKeep,
  const int & verbose)
{
  vector<string> vFtrsToKeep = loadOneColumnFile (ftrsToKeepFile, verbose);
  vector<string> vSnpsToKeep = loadOneColumnFile (snpsToKeepFile, verbose);
  Grid iGridL (largeGridFile, true, verbose);
  Grid iGridS (smallGridFile, false, verbose);
  
  map<string, string> mGenoPaths, mPhenoPaths;
  vector<string> vSubgroups;
  loadListsGenoAndPhenoFiles (genoPathsFile, phenoPathsFile, sbgrpToKeep,
			      mGenoPaths, mPhenoPaths, vSubgroups, mvlr,
			      verbose);
  
  vector<string> vSamples;
  vector<vector<size_t> > vvSampleIdxGenos, vvSampleIdxPhenos;
  loadSamples (mGenoPaths, mPhenoPaths, vSubgroups, vSamples,
	        vvSampleIdxGenos, vvSampleIdxPhenos, verbose);
  
  vector<map<string, vector<double> > > vSbgrp2Covars;
  loadCovariates (covarPathsFile, sbgrpToKeep, vSubgroups, vSamples,
		  vvSampleIdxGenos, vvSampleIdxPhenos, vSbgrp2Covars,
		  verbose);
  
  map<string, Ftr> mFtrs;
  map<string, vector<Ftr*> > mChr2VecPtFtrs;
  loadPhenos (mPhenoPaths, vSubgroups, vFtrsToKeep, mFtrs, verbose);
  loadFtrInfo (ftrCoordsFile, mFtrs, mChr2VecPtFtrs, verbose);
  
  map<string, Snp> mSnps;
  map<string, vector<Snp*> > mChr2VecPtSnps;
  if (snpCoordFile.empty())
    loadGenosAndSnpInfo (mGenoPaths, minMaf, vSubgroups, vSnpsToKeep,
			 mChr2VecPtFtrs, mSnps, mChr2VecPtSnps, verbose);
  else
  {
    loadGenos (mGenoPaths, minMaf, vSubgroups, vSnpsToKeep, mSnps, verbose);
    loadSnpInfo (snpCoordFile, mSnps, mChr2VecPtSnps, verbose);
  }
  
  inferAssos (mFtrs, mChr2VecPtSnps, vvSampleIdxPhenos, vvSampleIdxGenos,
	      anchor, lenCis, whichStep, needQnorm, vSbgrp2Covars, iGridL,
	      iGridS, whichBfs, mvlr, propFitSigma, verbose);
  if (whichStep == 2 || whichStep == 4 || whichStep == 5)
    makePerms (mFtrs, vvSampleIdxPhenos, whichStep, needQnorm, vSbgrp2Covars,
	       iGridL, iGridS, mvlr, propFitSigma, nbPerms, seed, trick,
	       whichPermSep, whichPermBf, useMaxBfOverSnps, verbose);
  
  writeRes (outPrefix, outRaw, mFtrs, mSnps, vSubgroups, vSbgrp2Covars,
	    whichStep, iGridL, iGridS, whichBfs, whichPermSep, mvlr,
	    seed, whichPermBf, verbose);
}

#ifdef EQTLBMA_MAIN

int main (int argc, char ** argv)
{
  int verbose = 1, whichStep = 0, trick = 0, whichPermSep = 1;
  size_t lenCis = 100000, nbPerms = 0, seed = string::npos;
  float minMaf = 0.0, propFitSigma = 0.0;
  bool outRaw = false, needQnorm = false, mvlr = false,
    useMaxBfOverSnps = false;
  string genoPathsFile, snpCoordFile, phenoPathsFile, ftrCoordsFile,
    anchor = "FSS", outPrefix, covarPathsFile, largeGridFile, smallGridFile,
    whichBfs = "gen", whichPermBf = "gen", ftrsToKeepFile, snpsToKeepFile,
    sbgrpToKeep;
  
  parseArgs (argc, argv, genoPathsFile, snpCoordFile, phenoPathsFile,
	     ftrCoordsFile, anchor, lenCis, outPrefix, outRaw, whichStep,
	     needQnorm, minMaf, covarPathsFile, largeGridFile, smallGridFile,
	     whichBfs, mvlr, propFitSigma, nbPerms, seed, trick, whichPermSep,
	     whichPermBf, useMaxBfOverSnps, ftrsToKeepFile, snpsToKeepFile, 
	     sbgrpToKeep, verbose);
  
  time_t startRawTime, endRawTime;
  if (verbose > 0)
  {
    time (&startRawTime);
    cout << "START " << argv[0] << " (" << time2string (startRawTime) << ")"
	 << endl
	 << "compiled -> " << __DATE__ << " " << __TIME__
	 << endl << flush;
    printCmdLine (cout, argc, argv);
  }
  
  run (genoPathsFile, snpCoordFile, phenoPathsFile, ftrCoordsFile, anchor,
       lenCis, outPrefix, outRaw, whichStep, needQnorm, minMaf, covarPathsFile,
       largeGridFile, smallGridFile, whichBfs, mvlr, propFitSigma, nbPerms, seed,
       trick, whichPermSep, whichPermBf, useMaxBfOverSnps, ftrsToKeepFile,
       snpsToKeepFile, sbgrpToKeep, verbose);
  
  if (verbose > 0)
  {
    time (&endRawTime);
    cout << "END " << argv[0] << " (" << time2string (endRawTime) << ")"
	 << endl
	 << "elapsed -> " << elapsedTime(startRawTime, endRawTime)
	 << endl
	 << "max.mem -> " << getMaxMemUsedByProcess2Str ()
	 << endl;
  }
  
  return EXIT_SUCCESS;
}

#endif
