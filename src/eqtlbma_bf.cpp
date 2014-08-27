/** \file eqtlbma_bf.cpp
 *
 *  `eqtlbma_bf' performs eQTL mapping in multiple subgroups via a Bayesian model.
 *  Copyright (C) 2012-2013 Timothée Flutre
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

#include <cmath>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <getopt.h>
#include <libgen.h>
#include <sys/stat.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <set>
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

#include <omp.h>

#include "utils/utils_io.hpp"
#include "utils/utils_math.hpp"
using namespace utils;

#include "quantgen/data_loader.hpp"
#include "quantgen/gene.hpp"
#include "quantgen/snp.hpp"
#include "quantgen/gene_snp_pair.hpp"
using namespace quantgen;

#include "tabix/bgzf.h"
#include "tabix/tabix.h"

/** \brief Display the help on stdout.
 */
void help(char ** argv)
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
       << "      --geno\tfile with absolute paths to genotype files" << endl
       << "\t\ttwo columns: subgroup identifier<space/tab>path to file" << endl
       << "\t\tadd '#' at the beginning of a line to comment it" << endl
       << "\t\tsubgroup file: can be in three formats (VCF/IMPUTE/custom)" << endl
       << "\t\tVCF: see specifications on 1kG website" << endl
       << "\t\tIMPUTE: row 1 is a header chr<del>name<del>coord<del>a1<del>a2" << endl
       << "\t\t followed by >sample1_a1a1<del>sample1_a1a2<del>sample1_a2a2<del>..." << endl
       << "\t\tcustom: genotypes as allele dose, same as for MatrixEQTL" << endl
       << "\t\t and missing data can be NA or -1 (as used by vcftools --012)" << endl
       << "      --scoord\tfile with the SNP coordinates" << endl
       << "\t\tcompulsory if custom genotype format; forbidden otherwise" << endl
       << "\t\tshould be in the BED format (delimiter: tab)" << endl
       << "\t\tSNPs in the genotype files without coordinate are skipped (see also --snp)" << endl
       << "\t\tif a tabix-indexed file is also present, it will be used" << endl
       << "      --exp\tfile with absolute paths to expression level files" << endl
       << "\t\ttwo columns: subgroup identifier<space/tab>path to file" << endl
       << "\t\tadd '#' at the beginning of a line to comment it" << endl
       << "\t\tsubgroup file: custom format, same as for MatrixEQTL" << endl
       << "\t\t row 1 for sample names, column 1 for gene names" << endl
       << "\t\tsubgroups can have different genes" << endl
       << "\t\tall genes should be in the --gcoord file" << endl
       << "      --gcoord\tfile with the gene coordinates" << endl
       << "\t\tshould be in the BED format (delimiter: tab)" << endl
       << "\t\tgenes in the exp level files without coordinates are skipped" << endl
       << "      --anchor\tgene boundary(ies) for the cis region" << endl
       << "\t\tdefault=TSS (assumed to be start in BED file)" << endl
       // << "\t\tdefault=TSS, can also be TSS+TES" << endl
       << "      --cis\tlength of half of the cis region (radius, in bp)" << endl
       << "\t\tapart from the anchor(s), default=100000" << endl
       << "      --inss\tfile with absolute paths to files with summary statistics" << endl
       << "\t\ttwo columns: subgroup identifier<space/tab>path to file" << endl
       << "\t\tadd '#' at the beginning of a line to comment it" << endl
       << "\t\tsstats file: custom format, similar to the one from --outss (see below)" << endl
       << "\t\t header should have gene, snp, n, sigmahat, betahat.geno and sebetahat.geno" << endl
       << "\t\t order doesn't matter" << endl
       << "      --out\tprefix for the output files" << endl
       << "\t\tall output files are gzipped and have a header line" << endl
       << "      --lik\tlikelihood to use" << endl
       << "\t\t'normal' (default)" << endl
       << "\t\t'poisson' or 'quasipoisson'" << endl
       << "      --analys\tanalysis to perform" << endl
       << "\t\t'sep': separate analysis of each subgroup" << endl
       << "\t\t'join': joint analysis of all subgroups" << endl
       << "      --outss\twrite the output file with all summary statistics" << endl
       << "      --outw\twrite the output file with the ABFs averaged over the grid" << endl
       << "\t\tgrid weights are uniformly equal" << endl
       << "      --qnorm\tquantile-normalize the exp levels to a N(0,1)" << endl
       << "      --maf\tminimum minor allele frequency (default=0.0)" << endl
       << "      --covar\tfile with absolute paths to covariate files" << endl
       << "\t\ttwo columns: subgroup identifier<space/tab>path to file" << endl
       << "\t\tcan be a single line (single subgroup)" << endl
       << "\t\tadd '#' at the beginning of a line to comment it" << endl
       << "\t\tsubgroup file: row 1 is a header sample<space/tab>covariate1 ..." << endl
       << "\t\tall sample names should be in the respective genotype and exp level files" << endl
       << "\t\tthe covariates should be numbers, no missing value is allowed" << endl
       << "\t\tsubgroups can have different covariates" << endl
       << "\t\tthe order of rows is not important" << endl
       << "      --gridL\tfile with a 'large' grid for prior variances in standardized effect sizes" << endl
       << "\t\tfirst column is phi^2 and second column is omega^2, no header" << endl
       << "\t\tthis grid is used with model 1 ('general alternative') trying to capture" << endl
       << "\t\t all sorts of heterogeneity" << endl
       << "\t\trequired with --analys join" << endl
       << "      --gridS\tfile with a 'small' grid of values for phi^2 and omega^2" << endl
       << "\t\tsame format as --gridL" << endl
       << "\t\trequired with --analyis join if --bfs is 'sin' or 'all'" << endl
       << "      --bfs\twhich Bayes Factors to compute for the joint analysis" << endl
       << "\t\tonly the Laplace-approximated BF from Wen and Stephens (AoAS 2013) is implemented" << endl
       << "\t\tif --outw, each BF for a given configuration is the average of the BFs over one of the grids, with equal weights" << endl
       << "\t\t'gen' (default): general way to capture any level of heterogeneity" << endl
       << "\t\t correspond to the consistent configuration with the large grid" << endl
       << "\t\t fixed-effect and maximum-heterogeneity BFs are also calculated" << endl
       << "\t\t'sin': compute also the BF for each singleton (subgroup-specific configuration)" << endl
       << "\t\t they use the small grid (BF_BMAlite is also reported)" << endl
       << "\t\t'all': compute also the BFs for all configurations (costly if many subgroups)" << endl
       << "\t\t all BFs use the small grid (BF_BMA is also reported)" << endl
       << "      --error\tmodel for the errors (if --analys join)" << endl
       << "\t\t'uvlr': default, errors are not correlated between subgroups (different individuals)" << endl
       << "\t\t'mvlr': errors can be correlated between subgroups (same individuals)" << endl
       << "\t\t'hybrid': errors can be correlated between pairs of subgroups (common individuals)" << endl
       << "      --fiterr\tparam used when estimating the variance of the errors (if --analys join, only with 'mvlr' or 'hybrid')" << endl
       << "\t\tdefault=0.5 but can be between 0 (null model) and 1 (full model)" << endl
       << "      --nperm\tnumber of permutations" << endl
       << "\t\tdefault=0, otherwise 10000 is recommended" << endl
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
       << "      --tricut\tcutoff for the trick (default=10)" << endl
       << "\t\tstop permutations once the nb of permutations for which permTestStat is more extreme" << endl
       << "\t\t than trueTestStat equals this cutoff" << endl
       << "      --permsep\twhich permutation procedure for the separate analysis" << endl
       << "\t\t0 (default): no permutations are done for the separate analysis" << endl
       << "\t\t1: use the minimum P-value over SNPs and subgroups as a test statistic (keeps correlations)" << endl
       << "\t\t2: use the minimum P-value over SNPs but in each subgroup separately (breaks correlations)" << endl
       << "      --pbf\twhich BF to use as the test statistic for the joint-analysis permutations" << endl
       << "\t\t'none' (default): no permutations are done for the joint analysis" << endl
       << "\t\t'gen': general BF (see --bfs above)" << endl
       << "\t\t'gen-sin': 0.5 BFgen + 0.5 BFsin (also called BF_BMAlite)" << endl
       << "\t\t'all': average over all configurations (also called BF_BMA)" << endl
       << "      --maxbf\tuse the maximum ABF over SNPs as test statistic for permutations" << endl
       << "\t\totherwise the average ABF over SNPs is used (more Bayesian)" << endl
       << "      --thread\tnumber of threads (default=1, parallelize over SNPs)" << endl
       << "      --snp\tfile with a list of SNPs to analyze" << endl
       << "\t\tone SNP name per line, useful when launched in parallel" << endl
       << "\t\tprogram exits if an empty file is given" << endl
       << "      --sbgrp\tidentifier of the subgroup to analyze" << endl
       << "\t\tuseful for quick analysis and debugging" << endl
       << "\t\tcan be 'sbgrp1+sbgrp3' for instance" << endl
       << endl;
}

/** \brief Display version and license information on stdout.
 */
void version(char ** argv)
{
  cout << argv[0] << " " << VERSION << endl
       << endl
       << "Copyright (C) 2012-2014 Timothée Flutre and Xiaoquan Wen." << endl
       << "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>" << endl
       << "This is free software; see the source for copying conditions.  There is NO" << endl
       << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl
       << endl
       << "Written by Timothée Flutre and Xiaoquan Wen." << endl;
}

/** \brief Parse the command-line arguments and check the values of the 
 *  compulsory ones.
 */
void
parseCmdLine(
  int argc,
  char ** argv,
  string & file_genopaths,
  string & snpCoordsFile,
  string & file_exppaths,
  string & file_genecoords,
  string & anchor,
  size_t & radius,
  string & file_sstats,
  string & out_prefix,
  bool & save_sstats,
  bool & save_weighted_abfs,
  string & likelihood,
  string & analysis,
  bool & need_qnorm,
  float & min_maf,
  string & covarFile,
  string & file_largegrid,
  string & file_smallgrid,
  string & bfs,
  string & error_model,
  float & prop_cov_errors,
  size_t & nb_types,
  size_t & nb_permutations,
  size_t & seed,
  int & trick,
  size_t & trick_cutoff,
  int & perm_sep,
  string & permbf,
  bool & use_max_bf,
  int & nb_threads,
  string & file_snpstokeep,
  vector<string> & subgroups_tokeep,
  int & verbose)
{
  int c = 0;
  while(true){
    static struct option long_options[] = {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'V'},
      {"verbose", required_argument, 0, 'v'},
      {"geno", required_argument, 0, 0},
      {"scoord", required_argument, 0, 0},
      {"exp", required_argument, 0, 0},
      {"gcoord", required_argument, 0, 0},
      {"anchor", required_argument, 0, 0},
      {"cis", required_argument, 0, 0},
      {"inss", required_argument, 0, 0},
      {"out", required_argument, 0, 0},
      {"outss", no_argument, 0, 0},
      {"outm", no_argument, 0, 0},
      {"outw", no_argument, 0, 0},
      {"lik", required_argument, 0, 0},
      {"analys", required_argument, 0, 0},
      {"qnorm", no_argument, 0, 0},
      {"maf", required_argument, 0, 0},
      {"covar", required_argument, 0, 0},
      {"gridL", required_argument, 0, 0},
      {"gridS", required_argument, 0, 0},
      {"bfs", required_argument, 0, 0},
      {"error", required_argument, 0, 0},
      {"fiterr", required_argument, 0, 0},
      {"nperm", required_argument, 0, 0},
      {"seed", required_argument, 0, 0},
      {"trick", required_argument, 0, 0},
      {"tricut", required_argument, 0, 0},
      {"permsep", required_argument, 0, 0},
      {"pbf", required_argument, 0, 0},
      {"maxbf", no_argument, 0, 0},
      {"thread", required_argument, 0, 0},
      {"snp", required_argument, 0, 0},
      {"sbgrp", required_argument, 0, 0},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    c = getopt_long(argc, argv, "hVv:",
		    long_options, &option_index);
    if(c == -1)
      break;
    switch (c){
    case 0:
      if(long_options[option_index].flag != 0)
	break;
      if(strcmp(long_options[option_index].name, "help") == 0){
	help(argv);
	exit(0);
      }
      if(strcmp(long_options[option_index].name, "version") == 0){
	version(argv);
	exit(0);
      }
      if(strcmp(long_options[option_index].name, "verbose") == 0){
	verbose = atoi(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "geno") == 0){
	file_genopaths = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "scoord") == 0){
	snpCoordsFile = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "exp") == 0){
	file_exppaths = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "gcoord") == 0){
	file_genecoords = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "anchor") == 0){
	anchor = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "cis") == 0){
	radius = atol(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "inss") == 0){
	file_sstats = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "out") == 0){
	out_prefix = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "outss") == 0){
	save_sstats = true;
	break;
      }
      if(strcmp(long_options[option_index].name, "outw") == 0){
	save_weighted_abfs = true;
	break;
      }
      if(strcmp(long_options[option_index].name, "lik") == 0){
	likelihood = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "analys") == 0){
	analysis = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "qnorm") == 0){
	need_qnorm = true;
	break;
      }
      if(strcmp(long_options[option_index].name, "maf") == 0){
	min_maf = atof(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "covar") == 0){
	covarFile = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "gridL") == 0){
	file_largegrid = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "gridS") == 0){
	file_smallgrid = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "bfs") == 0){
	bfs = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "error") == 0){
	error_model = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "fiterr") == 0){
	prop_cov_errors = atof(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "nperm") == 0){
	nb_permutations = atol(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "seed") == 0){
	seed = atol(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "trick") == 0){
	trick = atoi(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "tricut") == 0){
	trick_cutoff = atol(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "permsep") == 0){
	perm_sep = atoi(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "pbf") == 0){
	permbf = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "maxbf") == 0){
	use_max_bf = true;
	break;
      }
      if(strcmp(long_options[option_index].name, "thread") == 0){
	nb_threads = atoi(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "snp") == 0){
	file_snpstokeep = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "sbgrp") == 0){
	split(optarg, "+", subgroups_tokeep);
	break;
      }
    case 'h':
      help(argv);
      exit(0);
    case 'V':
      version(argv);
      exit(0);
    case 'v':
      verbose = atoi(optarg);
      break;
    case '?':
      printf("\n"); help(argv);
      abort();
    default:
      printf("\n"); help(argv);
      abort();
    }
  }
  
  if(file_sstats.empty()){
    if(file_genopaths.empty()){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: missing compulsory option --geno" << endl << endl;
      help(argv);
      exit(EXIT_FAILURE);
    }
    if(! doesFileExist(file_genopaths)){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: can't find " << file_genopaths << endl << endl;
      help(argv);
      exit(EXIT_FAILURE);
    }
    if(! snpCoordsFile.empty() && ! doesFileExist(snpCoordsFile)){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: can't find " << snpCoordsFile << endl << endl;
      help(argv);
      exit(EXIT_FAILURE);
    }
    if(file_exppaths.empty()){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: missing compulsory option --exp" << endl << endl;
      help(argv);
      exit(EXIT_FAILURE);
    }
    if(! doesFileExist(file_exppaths)){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: can't find " << file_exppaths << endl << endl;
      help(argv);
      exit(EXIT_FAILURE);
    }
    if(file_genecoords.empty()){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: missing compulsory option --gcoord" << endl << endl;
      help(argv);
      exit(EXIT_FAILURE);
    }
    if(! doesFileExist(file_genecoords)){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: can't find " << file_genecoords << endl << endl;
      help(argv);
      exit(EXIT_FAILURE);
    }
    if(anchor.empty()){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: SNPs in trans not yet implemented, see --anchor and --cis" << endl << endl;
      help(argv);
      exit(EXIT_FAILURE);
    }
    if(anchor != "TSS" && anchor != "TSS+TES"){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: --anchor should be TSS or TSS+TES" << endl << endl;
      help(argv);
      exit(EXIT_FAILURE);
    }
  } else{ // ! file_stats.empty()
    if(! doesFileExist(file_sstats)){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: can't find " << file_sstats << endl << endl;
      help(argv);
      exit(EXIT_FAILURE);
    }
    if(analysis != "join"){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: --inss requires --analys join" << endl << endl;
      help(argv);
      exit(EXIT_FAILURE);
    }
    if(error_model != "uvlr"){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: --inss requires --error uvlr" << endl << endl;
      help(argv);
      exit(EXIT_FAILURE);
    }
  }
  if(out_prefix.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --out" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  // check writing files is possible
  stringstream ssOutFile;
  ssOutFile << out_prefix << "_test.txt.gz";
  gzFile outStream;
  openFile(ssOutFile.str(), outStream, "wb");
  closeFile(ssOutFile.str(), outStream);
  remove(ssOutFile.str().c_str());
  if(likelihood.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --lik" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(likelihood != "normal" && likelihood != "poisson"
     && likelihood != "quasipoisson"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --lik " << likelihood << " is not valid" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(likelihood != "normal" && analysis != "join"
     && error_model != "uvlr"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --lik " << likelihood << " is not valid with --analys join"
	 << " and --error " << error_model << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --analys" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis != "sep" && analysis != "join"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --analys " << analysis << " is not valid" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis != "join" && file_largegrid.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --gridL with --analys join" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(! file_largegrid.empty() && ! doesFileExist(file_largegrid)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find " << file_largegrid << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis == "join" &&
     (bfs == "sin" || bfs == "all") &&
     file_smallgrid.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --gridS is required with --analys join and --bfs " << bfs << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(bfs!= "gen" && bfs != "sin" && bfs != "all"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --bfs " << bfs << " is not valid" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(error_model != "uvlr" && error_model != "mvlr" &&
     error_model != "hybrid"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --error " << error_model << " is not valid" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis == "join" && error_model == "mvlr")
    cerr << "WARNING: summary statistics per subgroup won't be saved with --error mvlr" << endl;
  if(trick != 0 && trick != 1 && trick != 2){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --trick " << trick << " is not valid" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(trick != 0 && trick_cutoff > nb_permutations){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --tricut " << trick_cutoff << " is larger than --nperm "
	 << nb_permutations << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(perm_sep != 0 && perm_sep != 1 && perm_sep != 2){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --permsep " << perm_sep << " is not valid" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis == "sep" && nb_permutations > 0 &&
     perm_sep != 1 && perm_sep != 2){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: if --type sep --nperm > 0, --permsep should be '1' or '2'" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(permbf != "none" && permbf != "gen" &&
     permbf != "gen-sin" && permbf != "all"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --pbf " << permbf << " is unvalid" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis == "join" && nb_permutations > 0 && permbf == "none"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: if --analys join --nperm > 0, --pbf should be different than 'none'" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis == "join" && nb_permutations > 0 && permbf != "none" &&
     bfs == "gen" && permbf != "gen"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: if --analys join --bfs gen --nperm > 0, --pbf should be 'gen'" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis == "join" && nb_permutations > 0 && permbf != "none" &&
     bfs == "sin" && permbf == "all"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: if --analys join --bfs sin --nperm > 0, --pbf should be 'gen' or 'gen-sin'" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(nb_threads <= 0){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --thread " << nb_threads << " is invalid" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(! file_snpstokeep.empty() && ! doesFileExist(file_snpstokeep)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find " << file_snpstokeep << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(seed == string::npos)
    seed = getSeed();
}

void loadListsGenoExplevelAndCovarFiles(
  const string & file_genopaths,
  const string & file_exppaths,
  const string & file_covarpaths,
  const vector<string> & subgroups_tokeep,
  const string & error_model,
  const int & verbose,
  map<string,string> & subgroup2genofile,
  map<string,string> & subgroup2explevelfile,
  map<string,string> & subgroup2covarfile,
  vector<string> & subgroups)
{
  map<string,string>::iterator it;
  
  // load paths to exp levels files
  subgroup2explevelfile = loadTwoColumnFile(file_exppaths, verbose);
  
  // erase subgroups with exp levels but not to keep
  it = subgroup2explevelfile.begin();
  while(it != subgroup2explevelfile.end()){
    if(! subgroups_tokeep.empty()
	&& find(subgroups_tokeep.begin(), subgroups_tokeep.end(), it->first)
       == subgroups_tokeep.end())
      subgroup2explevelfile.erase(it++);
    else
      ++it;
  }
  
  // load paths to geno files
  subgroup2genofile = loadTwoColumnFile(file_genopaths, verbose);
  
  // erase subgroups with geno but no exp levels
  it = subgroup2genofile.begin();
  while(it != subgroup2genofile.end()){
    if(subgroup2explevelfile.find(it->first) == subgroup2explevelfile.end())
      subgroup2genofile.erase(it++);
    else
      ++it;
  }
  
  // erase subgroups with exp levels but no geno
  it = subgroup2explevelfile.begin();
  while(it != subgroup2explevelfile.end()){
    if(subgroup2genofile.find(it->first) == subgroup2genofile.end())
      subgroup2explevelfile.erase(it++);
    else
      ++it;
  }
  
  // for correlated errors, check that, at least, all subgroups 
  // have genotypes in the same file
  if(error_model != "uvlr")
    for(it = subgroup2genofile.begin(); it != subgroup2genofile.end(); ++it)
      if(it->second.compare(subgroup2genofile.begin()->second) != 0){
	cerr << "ERROR: --error mvlr/hybrid requires the same genotypes in a single file for all subgroups"
	     << endl;
	exit(EXIT_FAILURE);
      }
  
  // identify the names of all remaining subgroups to be considered
  for(it = subgroup2explevelfile.begin(); it != subgroup2explevelfile.end();
       ++it)
    subgroups.push_back(it->first);
  
  // load paths to covariate files
  subgroup2covarfile = loadTwoColumnFile(file_covarpaths, verbose);
  
  // erase subgroups with covariates but neither geno nor exp levels
  it = subgroup2covarfile.begin();
  while(it != subgroup2covarfile.end()){
    if(find(subgroups.begin(), subgroups.end(), it->first) == subgroups.end())
      subgroup2covarfile.erase(it++);
    else++it;
  }
  
  if(verbose > 0)
    cout << "analyze " << subgroups.size() << " subgroup"
	 << (subgroups.size() > 1 ? "s" : "")
	 <<" (identifier):" << endl;
  for(vector<string>::const_iterator it = subgroups.begin();
       it != subgroups.end(); ++it)
    cout << *it << " (" << it - subgroups.begin() + 1 << ")" << endl;
}

void loadSamplesFromMatrixEqtl(
  const map<string, string> & subgroup2file,
  const string & data_type,
  const int & verbose,
  vector<string> & samples,
  map<string,vector<string> > & subgroup2samples)
{
  gzFile fileStream;
  string line;
  for(map<string,string>::const_iterator it = subgroup2file.begin();
       it != subgroup2file.end(); ++it){
    openFile(it->second, fileStream, "rb");
    if(! getline(fileStream, line)){
      cerr << "ERROR: problem with the header of file " << it->second << endl;
      exit(EXIT_FAILURE);
    }
    if(line.empty()){
      cerr << "ERROR: file " << it->second << " is empty" << endl;
      exit(EXIT_FAILURE);
    }
    closeFile(it->second, fileStream);
    if(it == subgroup2file.begin()){
      split(line, " \t", samples);
      if(samples[0] == "Id" || samples[0] == "id" || samples[0] == "ID")
	samples.erase(samples.begin());
      subgroup2samples.insert(make_pair(it->first, samples));
    }
    else{
      vector<string> tokens;
      split(line, " \t", tokens);
      if(tokens[0] == "Id" || tokens[0] == "id" || tokens[0] == "ID")
	tokens.erase(tokens.begin());
      subgroup2samples.insert(make_pair(it->first, tokens));
      for(vector<string>::const_iterator itT = tokens.begin();
	   itT != tokens.end(); ++itT)
	if(find(samples.begin(), samples.end(), *itT) == samples.end())
	  samples.push_back(*itT);
    }
  }
  if(verbose > 0){
    cout << "nb of samples (" << data_type << "): " << samples.size() << endl
	 << flush;
    for(map<string,vector<string> >::const_iterator it =
	   subgroup2samples.begin(); it != subgroup2samples.end(); ++it){
      cout << it->first << ": " << it->second.size() << " samples" << endl;
      // if(verbose > 1){
      // 	for(vector<string>::const_iterator itS = it->second.begin();
      // 	     itS != it->second.end(); ++itS)
      // 	  cout << *itS << endl;
      // }
    }
  }
}

void loadSamplesFromGenotypes(
  const map<string, string> & subgroup2genofile,
  const int & verbose,
  vector<string> & samples,
  map<string,vector<string> > & subgroup2samples_genotypes)
{
  gzFile fileStream;
  string line, sample;
  vector<string> tokens, tokens2;
  size_t i;
  for(map<string,string>::const_iterator it = subgroup2genofile.begin();
       it != subgroup2genofile.end(); ++it){
    openFile(it->second, fileStream, "rb");
    if(! getline(fileStream, line)){
      cerr << "ERROR: problem with the header of file " << it->second << endl;
      exit(EXIT_FAILURE);
    }
    if(line.empty()){
      cerr << "ERROR: file " << it->second << " is empty" << endl;
      exit(EXIT_FAILURE);
    }
    
    // if file in VCF format
    if(line.find("##fileformat=VCF") != string::npos){
      while(getline(fileStream, line)){
	if(line.find("#CHROM") == string::npos)
	  continue;
	split(line, " \t", tokens);
	closeFile(it->second, fileStream);
	subgroup2samples_genotypes.insert(
	  make_pair(it->first, vector<string>(tokens.begin()+9, tokens.end())));
	break;
      }
    }
    else{ // not VCF
      split(line, " \t", tokens);
      closeFile(it->second, fileStream);
      
      // if file in IMPUTE format
      if(tokens[0] == "chr"
	 && (tokens[1] == "name" || tokens[1] == "id")
	 && tokens[2] == "coord"
	 && tokens[3] == "a1"
	 && tokens[4] == "a2"){
	if((tokens.size() - 5) % 3 != 0){
	  cerr << "ERROR: the header of IMPUTE file " << it->second
	       << " is badly formatted" << endl;
	  exit(EXIT_FAILURE);
	}
	tokens2.clear();
	i = 5;
	while(i < tokens.size()){
	  sample = split(tokens[i], "_a", 0); // sampleX_a1a1, sampleX_a1a2 or sampleX_a2a2
	  tokens2.push_back(sample);
	  i = i + 3;
	}
	subgroup2samples_genotypes.insert(make_pair(it->first, tokens2));
	for(vector<string>::const_iterator itT = tokens2.begin();
	     itT != tokens2.end(); ++itT)
	  if(find(samples.begin(), samples.end(), *itT) == samples.end())
	    samples.push_back(*itT);
      }
      else{ // if file in MatrixEQTL format (allele dosage)
	if(tokens[0] == "Id" || tokens[0] == "id" || tokens[0] == "ID")
	  tokens.erase(tokens.begin());
	subgroup2samples_genotypes.insert(make_pair(it->first, tokens));
      }
    } // if file is not VCF
    
    for(vector<string>::const_iterator itS =
	   subgroup2samples_genotypes[it->first].begin();
	 itS != subgroup2samples_genotypes[it->first].end(); ++itS){
      if(find(samples.begin(), samples.end(), *itS) == samples.end())
	samples.push_back(*itS);
    }
  } // end for loop over subgroups
  
  if(verbose > 0){
    cout << "nb of samples (genotypes): " << samples.size() << endl << flush;
    for(map<string,vector<string> >::const_iterator it =
	   subgroup2samples_genotypes.begin();
	 it != subgroup2samples_genotypes.end(); ++it){
      cout << it->first << ": " << it->second.size() << " samples" << endl;
      // if(verbose > 1){
      // 	for(vector<string>::const_iterator itS = it->second.begin();
      // 	     itS != it->second.end(); ++itS)
      // 	  cout << *itS << endl;
      // }
    }
  }
}

void loadSamples(const map<string,string> & subgroup2genofile,
		 const map<string,string> & subgroup2explevelfile,
		 const map<string,string> & subgroup2covarfile,
		 const string & error_model,
		 const int & verbose,
		 Samples & samples)
{
  if(verbose > 0)
    cout << "load samples ..." << endl << flush;
  
  vector<string> samples_explevels;
  map<string,vector<string> > subgroup2samples_explevels;
  loadSamplesFromMatrixEqtl(subgroup2explevelfile, "explevels", verbose,
			    samples_explevels, subgroup2samples_explevels);
  
  vector<string> samples_genotypes;
  map<string,vector<string> > subgroup2samples_genotypes;
  loadSamplesFromGenotypes(subgroup2genofile, verbose, samples_genotypes,
			   subgroup2samples_genotypes);
  
  // fill samples by merging samples_explevels and samples_genotypes
  samples.AddSamplesIfNew(samples_explevels);
  samples.AddSamplesIfNew(samples_genotypes);
  if(verbose > 0)
    cout << "total nb of samples: " << samples.GetTotalNbSamples()
	 << endl << flush;
  
  vector<string> samples_covariates;
  map<string,vector<string> > subgroup2samples_covariates;
  loadSamplesFromMatrixEqtl(subgroup2covarfile, "covariates", verbose,
			    samples_covariates, subgroup2samples_covariates);
  
  // check that each covariate sample has both geno and pheno
  for(vector<string>::const_iterator it = samples_covariates.begin();
       it != samples_covariates.end(); ++it)
    if(samples.IsAbsent(*it)){
      cerr << "ERROR: sample " << *it << " has covariates but neither expression levels nor genotypes" << endl;
      exit(EXIT_FAILURE);
    }
  
  // do the mapping between the vector of all samples and each subgroup
  samples.AddSamplesFromData(subgroup2samples_genotypes, "genotype");
  samples.AddSamplesFromData(subgroup2samples_explevels, "explevel");
  samples.AddSamplesFromData(subgroup2samples_covariates, "covariate");
  
  if(error_model == "hybrid"){
    if(verbose > 0){
      cout << "nb of samples with genotype and exp level (pairs of subgroups):" << endl;
      samples.ShowPairs(cout);
    }
    if(verbose > 1)
      samples.ShowAllMappings(cerr);
  }
}

/** \brief Parse the BED file
 */
void loadGeneInfo(const string & file_genecoords, const int & verbose,
		  map<string,Gene> & gene2object,
		  map<string,vector<Gene*> > & mChr2VecPtGenes)
{
  if(verbose > 0)
    cout << "load gene coordinates ..." << endl << flush;
  
  vector<string> lines;
  readFile(file_genecoords, lines);
  
  vector<string> tokens;
  for(size_t i = 0; i < lines.size(); ++i){
    split(lines[i], " \t", tokens);
    if(gene2object.find(tokens[3]) != gene2object.end())
      continue; // in case of redundancy
    if(tokens[1] == tokens[2]){
      cerr << "ERROR: start and end coordinates of " << tokens[3]
	   << " should be different (at least 1 bp)" << endl;
      exit(1);
    }
    Gene gene(tokens[3], tokens[0], tokens[1], tokens[2]);
    if(tokens.size() >= 6)
      gene.SetStrand(tokens[5]);
    gene2object.insert(make_pair(gene.GetName(), gene));
    if(mChr2VecPtGenes.find(gene.GetChromosome()) == mChr2VecPtGenes.end())
      mChr2VecPtGenes.insert(make_pair(gene.GetChromosome(),
				       vector<Gene*>()));
    mChr2VecPtGenes[gene.GetChromosome()].push_back(
      &(gene2object[gene.GetName()]));
  }
  
  // sort the genes per chr
  for(map<string,vector<Gene*> >::iterator it = mChr2VecPtGenes.begin();
       it != mChr2VecPtGenes.end(); ++it)
    sort(it->second.begin(), it->second.end(), pt_gene_lt_pt_gene);
  
  if(verbose > 0)
    cout << "total nb of genes with coordinates: " << gene2object.size()
	 << endl;
}

void loadExplevels(const map<string,string> & subgroup2explevelfile,
		   const int & verbose,
		   map<string,Gene> & gene2object)
{
  if(gene2object.empty())
    return;
  
  if(verbose > 0)
    cout << "load gene expression levels ..." << endl << flush;
  
  gzFile explevelstream;
  string subgroup, explevelfile, line;
  vector<string> tokens;
  size_t nb_samples, nb_lines, nb_genes_tokeep_per_subgroup;
  
  for(map<string,string>::const_iterator it = subgroup2explevelfile.begin();
       it != subgroup2explevelfile.end(); ++it){
    nb_genes_tokeep_per_subgroup = 0;
    nb_lines = 0;
    subgroup = it->first;
    explevelfile = it->second;
    openFile(explevelfile, explevelstream, "rb");
    if(! getline(explevelstream, line)){
      cerr << "ERROR: problem with the header of file " << explevelfile
	   << endl;
      exit(EXIT_FAILURE);
    }
    ++nb_lines;
    split(line, " \t", tokens);
    if(tokens[0] == "Id" || tokens[0] == "id" || tokens[0] == "ID")
      nb_samples = tokens.size() - 1;
    else
      nb_samples = tokens.size();
    
    while(getline(explevelstream, line)){
      ++nb_lines;
      split(line, " \t", tokens);
      if(tokens.size() != nb_samples + 1){
	cerr << "ERROR: not enough columns on line " << nb_lines << " of file "
	     << explevelfile << " (" << tokens.size() << " != "
	     << nb_samples + 1 << ")" << endl;
	exit(EXIT_FAILURE);
      }
      if(gene2object.find(tokens[0]) == gene2object.end())
	continue; // skip if no coordinate
      gene2object[tokens[0]].AddSubgroup(subgroup, tokens.begin()+1,
					 tokens.end());
      ++nb_genes_tokeep_per_subgroup;
    }
    if(! gzeof(explevelstream)){
      cerr << "ERROR: can't read successfully file " << explevelfile
	   << " up to the end" << endl;
      exit(EXIT_FAILURE);
    }
    
    closeFile(explevelfile, explevelstream);
    if(verbose > 0)
      cout << subgroup << " (" << explevelfile << "): "
	   << (nb_lines-1) << " genes (to keep: "
	   << nb_genes_tokeep_per_subgroup << ")"
	   << endl << flush;
  }
  
  map<string,Gene>::iterator it = gene2object.begin();
  while(it != gene2object.end()){
    if(! it->second.HasExplevelsInAtLeastOneSubgroup()){
      cerr << "WARNING: skip gene " << it->second.GetName()
      	   << " because it has no expression level in any subgroup" << endl;
      gene2object.erase(it++);
    }
    else
      ++it;
  }
  
  if(verbose > 0)
    cout << "total nb of genes to analyze: " << gene2object.size() << endl;
}

void loadSnpsToKeep(const string & file_snpstokeep, const int & verbose,
		    set<string> & sSnpsToKeep)
{
  if(! file_snpstokeep.empty())
 {  
    string line;
    gzFile stream;
    vector<string> tokens;
    size_t nb_lines = 0;
    
    openFile(file_snpstokeep, stream, "rb");
    if(verbose > 0)
      cout <<"load file " << file_snpstokeep << " ..." << endl;
    
    while(getline(stream, line)){
      nb_lines++;
      split(line, " \t,", tokens);
      if(tokens.size() != 1){
	cerr << "ERROR: file " << file_snpstokeep
	     << " should have only one column"
	     << " at line " << nb_lines << endl;
	exit(EXIT_FAILURE);
      }
      if(tokens[0][0] == '#')
	continue;
      if(sSnpsToKeep.find(tokens[0]) == sSnpsToKeep.end())
	sSnpsToKeep.insert(tokens[0]);
    }
    
    if(! gzeof(stream)){
      cerr << "ERROR: can't read successfully file "
	   << file_snpstokeep << " up to the end" << endl;
      exit(EXIT_FAILURE);
    }
    closeFile(file_snpstokeep, stream);
    
    if(verbose > 0)
      cout << "nb of SNPs to keep: " << sSnpsToKeep.size() << endl;
  }
}

void loadGenosAndSnpInfoFromImpute(
  const map<string,string>::const_iterator & it_subgroup2genofile,
  const set<string> & sSnpsToKeep,
  const map<string, vector<Gene*> > & mChr2VecPtGenes,
  gzFile & genoStream,
  string & line,
  size_t & nb_lines,
  size_t & nb_snps_tokeep_per_subgroup,
  map<string, Snp> & snp2object,
  map<string, vector<Snp*> > & mChr2VecPtSnps)
{
  string subgroup = it_subgroup2genofile->first,
    genofile = it_subgroup2genofile->second;
  vector<string> tokens;
  split(line, " \t", tokens); // header line
  if((tokens.size() - 5) % 3 != 0){
    cerr << "ERROR: wrong number of columns on line " << nb_lines
	 << " of file " << genofile << endl;
    exit(EXIT_FAILURE);
  }
  size_t nb_samples = static_cast<size_t>((tokens.size() - 5) / 3);
  
  while(getline(genoStream, line)){
    ++nb_lines;
    split(line, " \t", tokens);
    if(tokens.size() != 3 * nb_samples + 5){
      cerr << "ERROR: not enough columns on line " << nb_lines << " of file "
	   << genofile << " (" << tokens.size() << " != "
	   << (3 * nb_samples + 5) << ")" << endl;
      exit(EXIT_FAILURE);
    }
    if(mChr2VecPtGenes.find(tokens[0]) == mChr2VecPtGenes.end())
      continue; // no gene on this chromosome
    if(! sSnpsToKeep.empty() && sSnpsToKeep.find(tokens[1])
	== sSnpsToKeep.end())
      continue;
    if(snp2object.find(tokens[1]) == snp2object.end()){
      Snp snp(tokens[1], tokens[0], tokens[2]);
      snp2object.insert(make_pair(snp.GetName(), snp));
    }
    snp2object[tokens[1]].AddSubgroup(subgroup, tokens.begin()+5,
				      tokens.end(), "impute");
    ++nb_snps_tokeep_per_subgroup;
  }
}

void
loadGenosAndSnpInfoFromVcf (
  const map<string,string>::const_iterator & it_subgroup2genofile,
  const set<string> & sSnpsToKeep,
  const map<string, vector<Gene*> > & mChr2VecPtGenes,
  gzFile & genoStream,
  string & line,
  size_t & nb_lines,
  size_t & nb_snps_tokeep_per_subgroup,
  map<string, Snp> & snp2object,
  map<string, vector<Snp*> > & mChr2VecPtSnps)
{
  string subgroup = it_subgroup2genofile->first,
    genofile = it_subgroup2genofile->second;
  vector<string> tokens, tokens2, tokens3;
  
  // skip the first lines of meta-data
  while(getline(genoStream, line)){
    ++nb_lines;
    if(line.find("#CHROM") != string::npos)
      break;
  }
  split(line, " \t", tokens);
  size_t nb_samples = tokens.size() - 9;
  
  while(getline(genoStream, line)){
    ++nb_lines;
    split(line, " \t", tokens);
    if(tokens.size() != nb_samples + 9){
      cerr << "ERROR: not enough columns on line " << nb_lines
	   << " of file " << genofile << " (" << tokens.size() << " != "
	   << nb_samples + 9 << ")" << endl;
      exit(EXIT_FAILURE);
    }
    if(tokens[8].find("GT") == string::npos){
      cerr << "ERROR: missing GT in 9-th field on line " << nb_lines
	   << " of file " << genofile << endl;
      exit(EXIT_FAILURE);
    }
    if(mChr2VecPtGenes.find(tokens[0]) == mChr2VecPtGenes.end())
      continue; // no gene on this chromosome
    if(! sSnpsToKeep.empty() && sSnpsToKeep.find(tokens[2])
	== sSnpsToKeep.end())
      continue;
    if(snp2object.find(tokens[2]) == snp2object.end()){
      Snp snp(tokens[2], tokens[0], tokens[1]);
      snp2object.insert(make_pair(snp.GetName(), snp));
    }
    snp2object[tokens[2]].AddSubgroup(subgroup, tokens.begin()+9,
				      tokens.end(), "vcf");
    ++nb_snps_tokeep_per_subgroup;
  }
}

void duplicateGenosPerSnpInAllSubgroups(
  const map<string,string> & subgroup2genofile,
  const int & verbose,
  map<string, Snp> & snp2object)
{
  if(verbose > 0)
    cout << "duplicate genotypes for other subgroups (same file) ..."
	 << flush << endl;
  map<string,string>::const_iterator it = subgroup2genofile.begin();
  string first_subgroup = it->first;
  ++it;
  while(it != subgroup2genofile.end()){
    clock_t startTime = clock();
    for(map<string,Snp>::iterator it_snp = snp2object.begin();
	 it_snp != snp2object.end(); ++it_snp)
      it_snp->second.DuplicateGenotypesFromFirstSubgroup(first_subgroup,
							 it->first);
    if(verbose > 0)
      cout << it->first << " (" << it->second << "): "
	   << snp2object.size() << " SNPs duplicated in "
	   << fixed << setprecision(2) << getElapsedTime(startTime) << " sec"
	   << endl << flush;
    ++it;
  }
}

void loadGenosAndSnpInfo(
  const map<string, string> & subgroup2genofile,
  const float & min_maf,
  const set<string> & sSnpsToKeep,
  const map<string, vector<Gene*> > & mChr2VecPtGenes,
  const int & verbose,
  map<string, Snp> & snp2object,
  map<string, vector<Snp*> > & mChr2VecPtSnps)
{
  if(verbose > 0)
    cout << "load genotypes and SNP coordinates ..." << endl << flush;
  
  gzFile genoStream;
  string line;
  size_t nb_lines, nb_snps_tokeep_per_subgroup;
  
  bool same_files = false;
  for(map<string,string>::const_iterator it = subgroup2genofile.begin();
       it != subgroup2genofile.end(); ++it){
    if(it != subgroup2genofile.begin() && it->second.compare(
	  subgroup2genofile.begin()->second) == 0){
      same_files = true;
      break; // avoid loading same file several times
    }
    
    clock_t startTime = clock();
    nb_snps_tokeep_per_subgroup = 0;
    nb_lines = 0;
    openFile(it->second, genoStream, "rb");
    if(! getline(genoStream, line)){
      cerr << "ERROR: problem with the header of file " << it->second << endl;
      exit(EXIT_FAILURE);
    }
    ++nb_lines;
    
    if(line.find("##fileformat=VCF") != string::npos) // VCF format
      loadGenosAndSnpInfoFromVcf(it, sSnpsToKeep, mChr2VecPtGenes,
				 genoStream, line, nb_lines,
				 nb_snps_tokeep_per_subgroup, snp2object,
				 mChr2VecPtSnps);
    else if(line.find("chr") != string::npos
	    && (line.find("name") != string::npos
		|| line.find("id") != string::npos)
	    && line.find("coord") != string::npos
	    && line.find("a1") != string::npos
	    && line.find("a2") != string::npos) // IMPUTE format
      loadGenosAndSnpInfoFromImpute(it, sSnpsToKeep, mChr2VecPtGenes,
				    genoStream, line, nb_lines,
				    nb_snps_tokeep_per_subgroup,
				    snp2object, mChr2VecPtSnps);
    else if(line.compare(0, 2, "id") == 0 || line.compare(0, 2, "Id") == 0
	     || line.compare(0, 2, "ID") == 0){
      cerr << "ERROR: file " << it->second
	   << " seems to be in the custom format but --scoord is missing"
	   << endl;
      exit(EXIT_FAILURE);
    }
    else{
      cerr << "ERROR: can't recognize the format of file " << it->second
	   << endl;
      exit(EXIT_FAILURE);
    }
    
    if(! gzeof(genoStream)){
      cerr << "ERROR: can't read successfully file " << it->second
	   << " up to the end" << endl;
      exit(EXIT_FAILURE);
    }
    
    closeFile(it->second , genoStream);
    if(verbose > 0)
      cout << it->first << " (" << it->second << "): " << (nb_lines-1)
	   << " SNPs (to keep: " << nb_snps_tokeep_per_subgroup << ", loaded in "
	   << fixed << setprecision(2) << getElapsedTime(startTime) << " sec)"
	   << endl << flush;
  }
  
  if(verbose > 0)
    cout << "discard SNPs with missing values ..." << endl << flush;
  map<string,Snp>::iterator it = snp2object.begin();
  while(it != snp2object.end()){
    it->second.EraseIfMissingValuesPerSubgroup();
    if(! it->second.HasGenotypesInAtLeastOneSubgroup()){
      if(verbose > 0)
	cerr << "WARNING: skip SNP " << it->second.GetName()
	     << " because it has missing values in each subgroup" << endl;
      snp2object.erase(it++);
    } else
      ++it;
  }
  
  if(min_maf > 0){
    if(verbose > 0)
      cout << "filter SNPs with MAF < " << min_maf << " ..." << endl << flush;
    map<string,Snp>::iterator it = snp2object.begin();
    while(it != snp2object.end()){
      it->second.EraseIfLowMafPerSubgroup(min_maf);
      if(! it->second.HasGenotypesInAtLeastOneSubgroup()){
	if(verbose > 0)
	  cerr << "WARNING: skip SNP " << it->second.GetName()
	       << " because it has a low MAF in each subgroup" << endl;
	snp2object.erase(it++);
      } else
	++it;
    }
  }
  
  if(same_files)
    duplicateGenosPerSnpInAllSubgroups(subgroup2genofile, verbose,
				       snp2object);
  
  for(map<string, Snp>::const_iterator it = snp2object.begin();
      it != snp2object.end(); ++it){
    if(mChr2VecPtSnps.find(it->second.GetChromosome()) == mChr2VecPtSnps.end())
      mChr2VecPtSnps.insert(make_pair(it->second.GetChromosome(), vector<Snp*>()));
    mChr2VecPtSnps[it->second.GetChromosome()].push_back(&(snp2object[it->second.GetName()]));
  }
  
  // sort the SNPs per chr
  for(map<string, vector<Snp*> >::iterator it = mChr2VecPtSnps.begin();
       it != mChr2VecPtSnps.end(); ++it)
    sort(it->second.begin(), it->second.end(), pt_snp_lt_pt_snp);
  
  if(verbose > 0)
    cout << "nb of SNPs: " << snp2object.size() << endl;
}

/** \brief Parse the tabix-indexed BED file
 */
void loadSnpInfo(const string & file_snpcoords,
		 const string & file_snpcoords_idx,
		 const set<string> & sSnpsToKeep,
		 const map<string,Gene> gene2object,
		 const string & anchor,
		 const size_t & radius,
		 const int & verbose,
		 map<string, Snp> & snp2object,
		 map<string, vector<Snp*> > & mChr2VecPtSnps)
{
  struct stat stat_bed, stat_tbi;
  stat(file_snpcoords.c_str(), &stat_bed);
  stat(file_snpcoords_idx.c_str(), &stat_tbi);
  if(stat_bed.st_mtime > stat_tbi.st_mtime){
    cerr << "ERROR: index file (" << file_snpcoords_idx
	 << ") is older than data file (" << file_snpcoords << ")" << endl;
    exit(EXIT_FAILURE);
  }
  
  tabix_t * t;
  if((t = ti_open(file_snpcoords.c_str(), 0)) == 0){
    cerr << "ERROR: fail to open the data file (tabix)" << endl;
    exit(EXIT_FAILURE);
  }
  if(ti_lazy_index_load(t) < 0){
    cerr << "ERROR: failed to load the index file (tabix)" << endl;
    exit(EXIT_FAILURE);
  }
  
  const char *s;
  ti_iter_t iter;
  vector<string> tokens;
  for(map<string,Gene>::const_iterator it = gene2object.begin();
      it != gene2object.end(); ++it){
    int t_id, t_beg, t_end, len;
    if(ti_parse_region(t->idx, it->second.GetRegionInTabixFormat(anchor, radius).c_str(),
		       &t_id, &t_beg, &t_end) == 0){
      iter = ti_queryi(t, t_id, t_beg, t_end);
      while((s = ti_read(t, iter, &len)) != 0){
	
	split(string(s), "\t", tokens);
	if(! sSnpsToKeep.empty() && sSnpsToKeep.find(tokens[3])
	   == sSnpsToKeep.end())
	  continue;
	if(snp2object.find(tokens[3]) != snp2object.end())
	  continue; // in case of redundancy
	Snp snp(tokens[3], tokens[0], tokens[2]);
	snp2object.insert(make_pair(snp.GetName(), snp));
	
      }
      ti_iter_destroy(iter);
    }
  }
  ti_close(t);
}

/** \brief Parse the unindexed BED file
 */
void loadSnpInfo(const string & snpCoordsFile,
		 const set<string> & sSnpsToKeep,
		 const int & verbose,
		 map<string, Snp> & snp2object,
		 map<string, vector<Snp*> > & mChr2VecPtSnps)
{
  gzFile snpCoordsStream;
  openFile(snpCoordsFile, snpCoordsStream, "rb");
  string line;
  vector<string> tokens;
  while(getline(snpCoordsStream, line)){
    split (line, " \t", tokens);
    if(! sSnpsToKeep.empty() && sSnpsToKeep.find(tokens[3])
	== sSnpsToKeep.end())
      continue;
    if(snp2object.find(tokens[3]) != snp2object.end())
      continue; // in case of redundancy
    if(tokens[1] == tokens[2]){
      cerr << "ERROR: start and end coordinates of " << tokens[3]
	   << " should be different (at least 1 bp)" << endl;
      exit(1);
    }
    Snp snp(tokens[3], tokens[0], tokens[2]);
    snp2object.insert(make_pair(snp.GetName(), snp));
  }
  if(! gzeof(snpCoordsStream))
 {
    cerr << "ERROR: can't read successfully file " << snpCoordsFile
	 << " up to the end" << endl;
    exit(EXIT_FAILURE);
  }
  closeFile(snpCoordsFile, snpCoordsStream);
}

/** \brief Parse the BED file (indexed or not)
 */
void loadSnpInfo(const string & file_snpcoords,
		 const set<string> & sSnpsToKeep,
		 const map<string,Gene> gene2object,
		 const string & anchor,
		 const size_t & radius,
		 const int & verbose,
		 map<string, Snp> & snp2object,
		 map<string, vector<Snp*> > & mChr2VecPtSnps)
{
  if(verbose > 0)
    cout << "load SNP coordinates";
  clock_t startTime = clock();
  
  stringstream file_snpcoords_idx;
  file_snpcoords_idx << file_snpcoords << ".tbi";
  if(doesFileExist(file_snpcoords_idx.str())){
    cout << " (tabix-indexed BED file) ..." << endl << flush;
    loadSnpInfo(file_snpcoords, file_snpcoords_idx.str(), sSnpsToKeep,
		gene2object, anchor, radius, verbose, snp2object,
		mChr2VecPtSnps);
  }
  else{
    cout << " (unindexed BED file) ..." << endl << flush;
    loadSnpInfo(file_snpcoords, sSnpsToKeep, verbose,
		snp2object, mChr2VecPtSnps);
  }
  
  for(map<string, Snp>::const_iterator it = snp2object.begin();
      it != snp2object.end(); ++it){
    if(mChr2VecPtSnps.find(it->second.GetChromosome()) == mChr2VecPtSnps.end())
      mChr2VecPtSnps.insert(make_pair(it->second.GetChromosome(), vector<Snp*>()));
    mChr2VecPtSnps[it->second.GetChromosome()].push_back(&(snp2object[it->second.GetName()]));
  }
  
  // sort the SNPs per chr
  for(map<string,vector<Snp*> >::iterator it = mChr2VecPtSnps.begin();
       it != mChr2VecPtSnps.end(); ++it)
    sort(it->second.begin(), it->second.end(), pt_snp_lt_pt_snp);
  
  if(verbose > 0)
    cout << "total nb of SNPs with coordinates: " << snp2object.size()
	 << " (loaded in " << fixed << setprecision(2)
	 << getElapsedTime(startTime) << " sec)" << endl;
}

void loadGenos(const map<string, string> & subgroup2genofile,
	       const float & min_maf, const int & verbose,
	       map<string, Snp> & snp2object)
{
  if(snp2object.empty())
    return;
  
  if(verbose > 0)
    cout << "load genotypes (custom format) ..." << endl << flush;
  
  gzFile genoStream;
  string line;
  vector<string> tokens;
  size_t nb_samples, nb_lines, nb_snps_tokeep_per_subgroup;
  
  bool same_files = false;
  for(map<string,string>::const_iterator it = subgroup2genofile.begin();
       it != subgroup2genofile.end(); ++it){
    if(it != subgroup2genofile.begin() && it->second.compare(
	  subgroup2genofile.begin()->second) == 0){
      same_files = true;
      break; // avoid loading same file several times
    }
    
    clock_t startTime = clock();
    nb_snps_tokeep_per_subgroup = 0;
    nb_lines = 0;
    openFile(it->second, genoStream, "rb");
    if(! getline(genoStream, line)){
      cerr << "ERROR: problem with the header of file " << it->second << endl;
      exit(EXIT_FAILURE);
    }
    ++nb_lines;
    split(line, " \t", tokens);
    if(tokens[0] == "chr"){
      cerr << "ERROR: don't use --scoord if genotypes in IMPUTE format" << endl;
      exit(1);
    }
    if(tokens[0].find("##fileformat=VCF") != string::npos){
      cerr << "ERROR: don't use --scoord if genotypes in VCF format" << endl;
      exit(1);
    }
    if(tokens[0] == "Id" || tokens[0] == "id" || tokens[0] == "ID")
      nb_samples = tokens.size() - 1;
    else
      nb_samples = tokens.size();
    
    while(getline(genoStream, line)){
      ++nb_lines;
      split(line, " \t", tokens);
      if(tokens.size() != nb_samples + 1){
	cerr << "ERROR: not enough columns on line " << nb_lines << " of file "
	     << it->second << " (" << tokens.size() << " != "
	     << nb_samples + 1 << ")" << endl;
	exit(EXIT_FAILURE);
      }
      if(snp2object.find(tokens[0]) == snp2object.end())
	continue; // skip if no coordinate
      snp2object[tokens[0]].AddSubgroup(it->first, tokens.begin()+1,
					tokens.end(), "dose");
      ++nb_snps_tokeep_per_subgroup;
    }
    if(! gzeof(genoStream)){
      cerr << "ERROR: can't read successfully file " << it->second
	   << " up to the end" << endl;
      exit(EXIT_FAILURE);
    }
    
    closeFile(it->second, genoStream);
    if(verbose > 0)
      cout << it->first << " (" << it->second << "): " << (nb_lines-1)
	   << " SNPs (to keep: " << nb_snps_tokeep_per_subgroup << ", loaded in "
	   << fixed << setprecision(2) << getElapsedTime(startTime) << " sec)"
	   << endl << flush;
  }
  
  if(min_maf > 0){
    if(verbose > 0)
      cout << "filter SNPs with MAF < " << min_maf << " ..." << endl << flush;
    map<string,Snp>::iterator it = snp2object.begin();
    while(it != snp2object.end()){
      it->second.EraseIfLowMafPerSubgroup(min_maf);
      if(verbose > 1 &&
	 ! it->second.HasGenotypesInAtLeastOneSubgroup()){
	cerr << "WARNING: skip SNP " << it->second.GetName()
	     << " because it has a low MAF in each subgroup" << endl;
	snp2object.erase(it++);
      } else
	++it;
    }
  }
  
  if(same_files)
    duplicateGenosPerSnpInAllSubgroups(subgroup2genofile, verbose,
				       snp2object);
  
  if(verbose > 0)
    cout << "total nb of SNPs to analyze: " << snp2object.size() << endl;
}

void
loadListCovarFiles (
  const string & file_covarpaths,
  const string & sbgrpToKeep,
  const vector<string> & subgroups,
  map<string, string> & mCovarPaths,
  const int & verbose)
{
  mCovarPaths = loadTwoColumnFile (file_covarpaths, verbose);
  
  map<string, string>::iterator it = mCovarPaths.begin();
  while(it != mCovarPaths.end())
 {
    if(! sbgrpToKeep.empty() && sbgrpToKeep.compare(it->first) != 0)
   {
      mCovarPaths.erase (it++);
      continue;
    }
    if(find (subgroups.begin(), subgroups.end(), it->first)
	== subgroups.end())
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

/** \brief Allow to load several covariate files per subgroup.
 */
void
loadListCovarFiles (
  const string & file_covarpaths,
  const string & sbgrpToKeep,
  const vector<string> & subgroups,
  map<string, vector<string> > & mCovarPaths,
  const int & verbose)
{
  string line;
  gzFile stream;
  vector<string> tokens;
  size_t nb_lines = 0, nbLoadedCovarFiles = 0;
  
  openFile(file_covarpaths, stream, "rb");
  if(verbose > 0)
    cout <<"load file " << file_covarpaths << " ..." << endl;
  
  while(getline(stream, line)){
    nb_lines;
    split(line, " \t,", tokens);
    if(tokens.size() != 2){
      cerr << "ERROR: file " << file_covarpaths
	   << " should have only two columns at line " << nb_lines << endl;
      exit(EXIT_FAILURE);
    }
    if(tokens[0][0] == '#')
      continue;
    if(! sbgrpToKeep.empty() && sbgrpToKeep.compare(tokens[0]) != 0)
      continue;
    if(find(subgroups.begin(), subgroups.end(), tokens[0])
       == subgroups.end()){
      cerr << "WARNING: skip covariates of subgroup " << tokens[0]
	   << " as there is no corresponding genotype / phenotype files"
	   << endl;
      continue;
    }
    if(mCovarPaths.find(tokens[0]) == mCovarPaths.end())
      mCovarPaths.insert(make_pair(tokens[0], vector<string> ()));
    mCovarPaths[tokens[0]].push_back(tokens[1]);
    ++nbLoadedCovarFiles;
  }
  
  if(! gzeof(stream)){
    cerr << "ERROR: can't read successfully file "
	 << file_covarpaths << " up to the end" << endl;
    exit(EXIT_FAILURE);
  }
  closeFile(file_covarpaths, stream);
  
  if(verbose > 0)
    cout << "items loaded: " << nbLoadedCovarFiles << " files for "
	 << mCovarPaths.size() << " subgroups" << endl;
}

void loadCovariates(const map<string,string> subgroup2covarfile,
		    const int & verbose, Covariates & covariates)
{
  if(subgroup2covarfile.empty())
    return;
  
  if(verbose > 0)
    cout << "load covariates ..." << endl << flush;
  
  gzFile covarstream;
  vector<string> tokens;
  string subgroup, covarfile, line;
  size_t nb_samples, nb_lines;
  
  for(map<string,string>::const_iterator it = subgroup2covarfile.begin();
       it != subgroup2covarfile.end(); ++it){
    map<string,vector<string> > covariate2values;
    nb_lines = 0;
    subgroup = it->first;
    covarfile = it->second;
    openFile(covarfile, covarstream, "rb");
    if(! getline(covarstream, line)){
      cerr << "ERROR: problem with the header of file " << covarfile << endl;
      exit(EXIT_FAILURE);
    }
    ++nb_lines;
    split(line, " \t", tokens);
    if(tokens[0] == "Id" || tokens[0] == "id" || tokens[0] == "ID")
      nb_samples = tokens.size() - 1;
    else
      nb_samples = tokens.size();
    
    while(getline(covarstream, line)){
      ++nb_lines;
      split(line, " \t", tokens);
      if(tokens.size() != nb_samples + 1){
	cerr << "ERROR: not enough columns on line " << nb_lines << " of file "
	     << covarfile << " (" << tokens.size() << " != "
	     << nb_samples + 1 << ")" << endl;
	exit(EXIT_FAILURE);
      }
      covariate2values.insert(
	make_pair(tokens[0], vector<string>(tokens.begin()+1, tokens.end())));
    }
    if(! gzeof(covarstream)){
      cerr << "ERROR: can't read successfully file " << covarfile
	   << " up to the end" << endl;
      exit(EXIT_FAILURE);
    }
    
    closeFile(covarfile, covarstream);
    covariates.AddSubgroup(subgroup, covariate2values);
  }
  
  if(verbose > 0)
    for(map<string,string>::const_iterator it = subgroup2covarfile.begin();
	it != subgroup2covarfile.end(); ++it)
      cout << it->first << " (" << it->second << "): "
	   << covariates.GetNbCovariates(it->first)
	   << " covariates" << endl;
}

void loadRawInputData(
  const string & file_genopaths,
  const string & file_snpcoords,
  const string & file_exppaths,
  const string & file_genecoords,
  const string & anchor,
  const size_t & radius,
  const float & min_maf,
  const string & file_covarpaths,
  const string & error_model,
  const vector<string> & subgroups_tokeep,
  const set<string> & sSnpsToKeep,
  const int & verbose,
  vector<string> & subgroups,
  Samples & samples,
  map<string,Snp> & snp2object,
  map<string,vector<Snp*> > & mChr2VecPtSnps,
  Covariates & covariates,
  map<string,Gene> & gene2object)
{
  map<string,string> subgroup2genofile, subgroup2explevelfile,
    subgroup2covarfile;
  loadListsGenoExplevelAndCovarFiles(file_genopaths, file_exppaths,
				     file_covarpaths, subgroups_tokeep,
				     error_model, verbose, subgroup2genofile,
				     subgroup2explevelfile,
				     subgroup2covarfile, subgroups);
  
  loadSamples(subgroup2genofile, subgroup2explevelfile, subgroup2covarfile,
	      error_model, verbose, samples);
  
  loadCovariates(subgroup2covarfile, verbose, covariates);
  
  map<string,vector<Gene*> > mChr2VecPtGenes;
  loadGeneInfo(file_genecoords, verbose, gene2object, mChr2VecPtGenes);
  loadExplevels(subgroup2explevelfile, verbose, gene2object);
  if(gene2object.empty())
    return;
  
  if(file_snpcoords.empty())
    loadGenosAndSnpInfo(subgroup2genofile, min_maf, sSnpsToKeep,
  			mChr2VecPtGenes, verbose, snp2object, mChr2VecPtSnps);
  else{
    loadSnpInfo(file_snpcoords, sSnpsToKeep, gene2object, anchor, radius,
  		verbose, snp2object, mChr2VecPtSnps);
    loadGenos(subgroup2genofile, min_maf, verbose, snp2object);
  }
  if(snp2object.empty())
    return;
}

void loadListSstatsFile(
  const string & file_sstats,
  const int & verbose,
  map<string,string> & subgroup2sstatsfile)
{
  string line;
  gzFile stream;
  vector<string> tokens;
  size_t nb_lines = 0;
  
  openFile(file_sstats, stream, "rb");
  if(verbose > 0)
    cout <<"load file " << file_sstats << " ..." << endl;
  
  while(getline(stream, line)){
    ++nb_lines;
    split(line, " \t,", tokens);
    if(tokens.size() != 2){
      cerr << "ERROR: file " << file_sstats 
	   << " should have only two columns at line " << nb_lines << endl;
      exit(EXIT_FAILURE);
    }
    if(tokens[0][0] == '#')
      continue;
    if(subgroup2sstatsfile.find(tokens[0]) == subgroup2sstatsfile.end())
      subgroup2sstatsfile.insert(make_pair(tokens[0], tokens[1]));
  }
  
  if(! gzeof(stream)){
    cerr << "ERROR: can't read successfully file "
	 << file_sstats << " up to the end" << endl;
    exit(EXIT_FAILURE);
  }
  closeFile(file_sstats, stream);
  
  if(verbose > 0)
    cout << "items loaded: " << subgroup2sstatsfile.size() << " files" << endl;
}

void fillGeneSnpPairsWithSstats(
  const map<string,string> & subgroup2sstatsfile,
  const int & verbose,
  map<string,Gene> & gene2object,
  map<string,Snp> & snp2object)
{
  vector<string> lines, tokens;
  map<string,size_t> col2idx;
  col2idx["gene"] = string::npos;
  col2idx["snp"] = string::npos;
  col2idx["n"] = string::npos;
  col2idx["sigmahat"] = string::npos;
  col2idx["betahat.geno"] = string::npos;
  col2idx["sebetahat.geno"] = string::npos;
  Gene * pt_gene = NULL;
  Snp * pt_snp = NULL;
  vector<GeneSnpPair>::iterator it_gsp;
  size_t idx_snp = string::npos;
  
  // loop over subgroups
  for(map<string,string>::const_iterator it_sf = subgroup2sstatsfile.begin();
      it_sf != subgroup2sstatsfile.end(); ++it_sf){
    if(verbose > 0)
      cout << "load summary statistics for subgroup "
	   << it_sf->first << " ..." << endl;
    readFile(it_sf->second, lines);
    
    // parse file header
    for(map<string,size_t>::iterator it_c = col2idx.begin();
	it_c != col2idx.end(); ++it_c)
      it_c->second = string::npos;
    split(lines[0], "\t", tokens);
    for(size_t col_id = 0; col_id < tokens.size(); ++col_id)
      if(col2idx.find(tokens[col_id]) != col2idx.end())
	col2idx[tokens[col_id]] = col_id;
    for(map<string,size_t>::const_iterator it_c = col2idx.begin();
	it_c != col2idx.end(); ++it_c)
      if(it_c->second == string::npos){
	cerr << "ERROR: missing " << it_c->second << " in header of "
	     << it_sf->second << endl;
	exit(EXIT_FAILURE);
      }
    
    // parse file content
    for(size_t line_id = 1; line_id < lines.size(); ++line_id){
      split(lines[line_id], "\t", tokens);
      
      // get gene (create it if necessary)
      if(gene2object.find(tokens[col2idx["gene"]]) == gene2object.end())
	gene2object.insert(make_pair(tokens[col2idx["gene"]],
				     Gene(tokens[col2idx["gene"]])));
      pt_gene = &(gene2object[tokens[col2idx["gene"]]]);
      
      // get snp (create it if necessary)
      if(snp2object.find(tokens[col2idx["snp"]]) == snp2object.end())
	snp2object.insert(make_pair(tokens[col2idx["snp"]],
				    Snp(tokens[col2idx["snp"]])));
      pt_snp = &(snp2object[tokens[col2idx["snp"]]]);
      
      // get gene-snp pair (create it if necessary)
      idx_snp = pt_gene->FindIdxSnp(pt_snp);
      if(idx_snp == string::npos){
	pt_gene->AddCisSnp(pt_snp);
	it_gsp = pt_gene->AddGeneSnpPair(pt_snp->GetName(), "uvlr");
      } else
	it_gsp = pt_gene->FindGeneSnpPair(idx_snp);
      
      it_gsp->SetSstats(it_sf->first,
			atol(tokens[col2idx["n"]].c_str()),
			atof(tokens[col2idx["sigmahat"]].c_str()),
			atof(tokens[col2idx["betahat.geno"]].c_str()),
			atof(tokens[col2idx["sebetahat.geno"]].c_str()));
      
    } // end of loop over lines
    
  } // end of loop over subgroups
}

void loadSummaryStats(
  const string & file_sstats,
  const int & verbose,
  vector<string> & subgroups,
  map<string,Gene> & gene2object,
  map<string,Snp> & snp2object)
{
  map<string,string> subgroup2sstatsfile;
  loadListSstatsFile(file_sstats, verbose, subgroup2sstatsfile);
  
  keys2vec(subgroup2sstatsfile, subgroups);
  
  fillGeneSnpPairsWithSstats(subgroup2sstatsfile, verbose, gene2object,
  			     snp2object);
}

void testForAssociations(
  const bool & hasDataNotSstats,
  const map<string,vector<Snp*> > & mChr2VecPtSnps,
  const string & anchor,
  const size_t & radius,
  const vector<string> & subgroups,
  const Samples & samples,
  const string & likelihood,
  const string & analysis,
  const bool & need_qnorm,
  const Covariates & covariates,
  const Grid & iGridL,
  const Grid & iGridS,
  const string & bfs,
  const string & error_model,
  const float & prop_cov_errors,
  const int & verbose,
  map<string,Gene> & gene2object)
{
  if(verbose > 0){
    cout << "test for association between each pair gene-SNP ..." << endl
	 << "analysis=" << analysis
	 << " likelihood=" << likelihood
	 << " error_model=" << error_model;
    if(error_model != "uvlr") // i.e. if 'mvlr' or 'hybrid'
      cout << " prop_cov_errors=" << prop_cov_errors;
    if(hasDataNotSstats)
      cout << " anchor=" << anchor << " radius=" << radius;
    // cout << " threads=1";
    cout << endl << flush;
  }
  
  clock_t startTime = clock();
  size_t nbAnalyzedGenes = 0, nbAnalyzedPairs = 0;
  size_t countGenes = 0;
  
  for(map<string,Gene>::iterator itG = gene2object.begin();
       itG != gene2object.end(); ++itG){
    ++countGenes;
    if(verbose == 1)
      progressBar("", countGenes, gene2object.size());
    if(verbose > 1)
      cerr << "gene " << itG->first << endl;
    
    if(hasDataNotSstats){
      itG->second.SetCisSnps(mChr2VecPtSnps, anchor, radius);
      if(! itG->second.HasAtLeastOneCisSnpInAtLeastOneSubgroup()){
	if(verbose > 1)
	  cerr << "WARNING: skip gene " << itG->second.GetName()
	       << " because it has no SNP in cis" << endl;
	continue;
      }
    }
    
    if(analysis == "join" && error_model != "uvlr"
       && ! itG->second.HasExplevelsInAllSubgroups(subgroups)){
      if(verbose > 1)
	cerr << "WARNING: skip gene " << itG->second.GetName()
	     << " because option --error " << error_model
	     << " requires expression levels in all subgroups" << endl;
      continue;
    }
    itG->second.TestForAssociations(hasDataNotSstats, subgroups, samples,
				    likelihood, analysis, need_qnorm,
				    covariates, iGridL, iGridS, bfs,
				    error_model, prop_cov_errors, verbose-1);
    ++nbAnalyzedGenes;
    nbAnalyzedPairs += itG->second.GetNbGeneSnpPairs();
  }
  
  if(verbose > 0){
    if(verbose == 1)
      cout << " (" << fixed << setprecision(2) << getElapsedTime(startTime)
	   << " sec)" << endl << flush;
    cout << "nb of analyzed gene-SNP pairs: " << nbAnalyzedPairs
	 << " (" << nbAnalyzedGenes << " genes)" << endl;
  }
}

void makePermutationsSep(
  const vector<string> & subgroups,
  const Samples & samples,
  const string & likelihood,
  const bool & need_qnorm,
  const Covariates & covariates,
  const size_t & nb_permutations,
  const size_t & seed,
  const int & trick,
  const size_t & trick_cutoff,
  const int & perm_sep,
  const gsl_rng * rngPerm,
  const gsl_rng * rngTrick,
  const int & verbose,
  map<string, Gene> & gene2object)
{
  if(perm_sep == 1){
    clock_t startTime = clock();
    gsl_rng_set(rngPerm, seed);
    if(trick != 0)
      gsl_rng_set(rngTrick, seed);
    size_t countGenes = 0;
    for(map<string,Gene>::iterator itG = gene2object.begin();
	itG != gene2object.end(); ++itG){
      ++countGenes;
      if(verbose == 1)
	progressBar("sep", countGenes, gene2object.size());
      if(verbose > 1)
	cerr << "gene " << itG->first << endl;
      if(! itG->second.HasAtLeastOneCisSnpInAtLeastOneSubgroup())
	continue;
      itG->second.MakePermutationsSepAllSubgroups(subgroups, samples,
						  likelihood,
						  need_qnorm, covariates,
						  nb_permutations, trick,
						  trick_cutoff,
						  rngPerm, rngTrick);
    }
    if(verbose == 1)
      cout << " (" << fixed << setprecision(2) << getElapsedTime(startTime)
	   << " sec)" << endl << flush;
  }
  else if(perm_sep == 2){
    stringstream ss;
    for(vector<string>::const_iterator it_sbgrp = subgroups.begin();
	it_sbgrp != subgroups.end(); ++it_sbgrp){
      clock_t startTime = clock();
      gsl_rng_set(rngPerm, seed);
      if(trick != 0)
	gsl_rng_set(rngTrick, seed);
      size_t countGenes = 0;
      ss.str("");
      ss << "sep" << (it_sbgrp - subgroups.begin()) + 1;
      for(map<string,Gene>::iterator itG = gene2object.begin();
	  itG != gene2object.end(); ++itG){
	++countGenes;
	if(verbose == 1)
	  progressBar(ss.str(), countGenes, gene2object.size());
	if(verbose > 1)
	  cerr << "gene " << itG->first << endl;
	if(! itG->second.HasAtLeastOneCisSnpInAtLeastOneSubgroup())
	  continue;
	itG->second.MakePermutationsSepPerSubgroup(*it_sbgrp, samples,
						   likelihood,
						   need_qnorm, covariates,
						   nb_permutations, trick,
						   trick_cutoff,
						   rngPerm, rngTrick);
      }
      if(verbose == 1)
	cout << " (" << fixed << setprecision(2) << getElapsedTime(startTime)
	     << " sec)" << endl << flush;
    }
  }
}

void makePermutationsJoin(
  const vector<string> & subgroups,
  const Samples & samples,
  const string & likelihood,
  const bool & need_qnorm,
  const Covariates & covariates,
  const Grid & iGridL,
  const Grid & iGridS,
  const string & error_model,
  const float & prop_cov_errors,
  const size_t & nb_permutations,
  const size_t & seed,
  const int & trick,
  const size_t & trick_cutoff,
  const string & permbf,
  const bool & use_max_bf,
  const gsl_rng * rngPerm,
  const gsl_rng * rngTrick,
  const int & verbose,
  map<string, Gene> & gene2object)
{
  clock_t startTime = clock();
  gsl_rng_set(rngPerm, seed);
  if(trick != 0)
    gsl_rng_set(rngTrick, seed);
  size_t countGenes = 0;
  
  for(map<string,Gene>::iterator itG = gene2object.begin();
       itG != gene2object.end(); ++itG){
    ++countGenes;
    if(verbose == 1)
      progressBar("join", countGenes, gene2object.size());
    if(verbose > 1)
      cerr << "gene " << itG->first << endl;
    if(! itG->second.HasAtLeastOneCisSnpInAtLeastOneSubgroup())
      continue;
    if(error_model != "uvlr" &&
       ! itG->second.HasExplevelsInAllSubgroups(subgroups))
      continue;
    itG->second.MakePermutationsJoin(subgroups, samples, likelihood, need_qnorm,
				     covariates, iGridL, iGridS, error_model,
				     prop_cov_errors, nb_permutations, trick,
				     trick_cutoff, permbf, use_max_bf,
				     rngPerm, rngTrick);
  }
  
  if(verbose == 1)
    cout << " (" << fixed << setprecision(2) << getElapsedTime(startTime)
	 << " sec)" << endl << flush;
}

void
makePermutations(
  const vector<string> & subgroups,
  const Samples & samples,
  const string & likelihood,
  const string & analysis,
  const bool & need_qnorm,
  const Covariates & covariates,
  const Grid & iGridL,
  const Grid & iGridS,
  const string & error_model,
  const float & prop_cov_errors,
  const size_t & nb_permutations,
  const size_t & seed,
  const int & trick,
  const size_t & trick_cutoff,
  const int & perm_sep,
  const string & permbf,
  const bool & use_max_bf,
  const int & verbose,
  map<string, Gene> & gene2object)
{
  if(verbose > 0){
    cout << "get gene-level p-values by permuting expression levels ..." << endl
	 << "permutation"<< (nb_permutations > 1 ? "s=" : "=") << nb_permutations
	 << " seed=" << seed;
    if(trick != 0){
      cout << " trick=" << trick
	   << " trick_cutoff=" << trick_cutoff;
    }
    if(analysis == "sep")
      cout << " perm_sep=" << perm_sep;
    else if(analysis == "join")
      cout << " perm_bf=" << permbf;
    cout << " threads=" << omp_get_max_threads();
    cout << endl << flush;
  }
  
  gsl_rng * rngPerm = NULL, * rngTrick = NULL;
  gsl_rng_env_setup();
  rngPerm = gsl_rng_alloc (gsl_rng_default);
  if(rngPerm == NULL){
    cerr << "ERROR: can't allocate memory for the RNG" << endl;
    exit(EXIT_FAILURE);
  }
  if(trick != 0){
    rngTrick = gsl_rng_alloc (gsl_rng_default);
    if(rngTrick == NULL){
      cerr << "ERROR: can't allocate memory for the RNG" << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if(analysis == "sep" && perm_sep != 0)
    makePermutationsSep(subgroups, samples, likelihood, need_qnorm, covariates,
			nb_permutations, seed, trick, trick_cutoff,
			perm_sep, rngPerm, rngTrick,
			verbose, gene2object);
  if(analysis == "join" && permbf != "none")
    makePermutationsJoin(subgroups, samples, likelihood, need_qnorm, covariates, iGridL,
			 iGridS, error_model, prop_cov_errors, nb_permutations,
			 seed, trick, trick_cutoff, permbf, use_max_bf,
			 rngPerm, rngTrick, verbose, gene2object);
  
  gsl_rng_free(rngPerm);
  if(trick != 0)
    gsl_rng_free(rngTrick);
}

void writeResSstats(
  const string & out_prefix,
  const vector<string> & subgroups,
  const map<string,Gene> & gene2object,
  const map<string,Snp> & snp2object,
  const int & verbose)
{
  if(verbose > 0)
    cout << "write results of summary statistics in each subgroup ..."
	 << endl << flush;
  
  string subgroup, sep = "\t";
  
  for(size_t s = 0; s < subgroups.size(); ++s){
    subgroup = subgroups[s];
    stringstream ssOutFile, ssTxt;
    ssOutFile << out_prefix << "_sumstats_" << subgroup << ".txt.gz";
    if(verbose > 0)
      cout << "file " << ssOutFile.str() << endl << flush;
    gzFile outStream;
    openFile(ssOutFile.str(), outStream, "wb");
    
    ssTxt << "gene" << sep << "snp" << sep << "maf" << sep << "n" << sep
	  << "pve" << sep << "sigmahat" << sep << "betahat.geno" << sep
	  << "sebetahat.geno" << sep << "betapval.geno";
    ssTxt << endl;
    size_t nb_lines = 1;
    gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
    
    ssTxt.precision(6);
    ssTxt.setf(ios::scientific);
    for(map<string,Gene>::const_iterator it_gene = gene2object.begin();
	 it_gene != gene2object.end(); ++it_gene){
      if(! it_gene->second.HasAtLeastOneCisSnp(subgroup))
	continue;
      for(vector<GeneSnpPair>::const_iterator it_pair
	     = it_gene->second.BeginPair();
	   it_pair != it_gene->second.EndPair(); ++it_pair){
	if(! it_pair->HasResults(subgroup))
	  continue;
	ssTxt.str("");
	++nb_lines;
	ssTxt << it_gene->second.GetName()
	      << sep << it_pair->GetSnpName()
	      << sep << snp2object.find(it_pair->GetSnpName())->second.GetMinorAlleleFreq(subgroup)
	      << sep << it_pair->GetSampleSize(subgroup)
	      << sep << it_pair->GetPve(subgroup)
	      << sep << it_pair->GetSigmahat(subgroup)
	      << sep << it_pair->GetBetahatGeno(subgroup)
	      << sep << it_pair->GetSebetahatGeno(subgroup)
	      << sep << it_pair->GetBetapvalGeno(subgroup);
	ssTxt << "\n";
	gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
      } // end of loop over cis snps
      
    } // end of loop over genes
    
    closeFile(ssOutFile.str(), outStream);
    
  } // end of loop over subgroups
}

void writeResSepPermPval(
  const string & out_prefix,
  const map<string, Gene> & gene2object,
  const vector<string> & subgroups,
  const size_t & seed,
  const int & verbose)
{
  if(verbose > 0)
    cout << "write results of gene-level P-values"
	 << " in each subgroup separately"
	 << " (perm_sep=2) ..."
	 << endl << flush;
  
  string sep = "\t";
  
  for(vector<string>::const_iterator it_sbgrp  = subgroups.begin();
      it_sbgrp != subgroups.end(); ++it_sbgrp){
    stringstream ssOutFile, ssTxt;
    ssOutFile << out_prefix << "_sepPermPvals_" << *it_sbgrp << ".txt.gz";
    if(verbose > 0)
      cout << "file " << ssOutFile.str() << endl << flush;
    gzFile outStream;
    openFile(ssOutFile.str(), outStream, "wb");
    size_t nb_lines = 0;
    
    ssTxt << "# seed=" << seed << endl;
    ++nb_lines;
    gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
    
    ssTxt.str("");
    ssTxt << "gene" << sep << "nb.snps" << sep << "sep.perm.pval"
	  << sep << "nb.permutations" << sep << "true.min.pval" << endl;
    ++nb_lines;
    gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
    
    for(map<string,Gene>::const_iterator it_gene = gene2object.begin();
	 it_gene != gene2object.end(); ++it_gene){
      if(it_gene->second.GetNbGeneSnpPairs() > 0){
	ssTxt.str("");
	++nb_lines;
	ssTxt << it_gene->second.GetName()
	      << sep << it_gene->second.GetNbGeneSnpPairs(*it_sbgrp)
	      << sep << it_gene->second.GetPermutationPvalueSep(*it_sbgrp)
	      << sep << it_gene->second.GetNbPermutationsSep(*it_sbgrp);
	ssTxt.setf(ios::scientific);
	ssTxt << sep << setprecision(6) << it_gene->second.GetTrueMinPval(*it_sbgrp);
	ssTxt.unsetf(ios::scientific);
	ssTxt << "\n";
	gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
      }
    }
    
    closeFile(ssOutFile.str(), outStream);
  }
}

void writeResSepPermPval(
  const string & out_prefix,
  const map<string, Gene> & gene2object,
  const size_t & seed,
  const int & verbose)
{
  if(verbose > 0)
    cout << "write results of gene-level P-values"
	 << " in each subgroup separately"
	 << " (perm_sep=1) ..."
	 << endl << flush;
  
  string sep = "\t";
  
  stringstream ssOutFile, ssTxt;
  ssOutFile << out_prefix << "_sepPermPvals.txt.gz";
  if(verbose > 0)
    cout << "file " << ssOutFile.str() << endl << flush;
  gzFile outStream;
  openFile(ssOutFile.str(), outStream, "wb");
  size_t nb_lines = 0;
  
  ssTxt << "# seed=" << seed << endl;
  ++nb_lines;
  gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
  
  ssTxt.str("");
  ssTxt << "gene" << sep << "nb.snps" << sep << "sep.perm.pval"
	<< sep << "nb.permutations" << sep << "true.min.pval" << endl;
  ++nb_lines;
  gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
  
  for(map<string,Gene>::const_iterator it_gene = gene2object.begin();
      it_gene != gene2object.end(); ++it_gene){
    if(it_gene->second.GetNbGeneSnpPairs() > 0){
      ssTxt.str("");
      ++nb_lines;
      ssTxt << it_gene->second.GetName()
	    << sep << it_gene->second.GetNbGeneSnpPairs()
	    << sep << it_gene->second.GetPermutationPvalueSep()
	    << sep << it_gene->second.GetNbPermutationsSep();
      ssTxt.setf(ios::scientific);
      ssTxt << sep << setprecision(6) << it_gene->second.GetTrueMinPval();
      ssTxt.unsetf(ios::scientific);
      ssTxt << "\n";
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
    }
  }
  
  closeFile(ssOutFile.str(), outStream);
}

/** \brief 
 *  \note the 'const' BFs use the large grid while the other use the small one,
 *  but the header line lists the large grid only
 */
void writeResAbfsRaw(
  const string & out_prefix,
  const map<string,Gene> & gene2object,
  const size_t & nb_subgroups,
  const Grid & iGridL,
  const Grid & iGridS,
  const string & bfs,
  const int & verbose)
{
  if(verbose > 0)
    cout << "write results of Bayes Factors (one per grid value) ..."
	 << endl << flush;
  
  string sep = "\t";
  
  gsl_combination * comb;
  stringstream ssOutFile, ssConfig, ssTxt;
  ssOutFile << out_prefix << "_l10abfs_raw.txt.gz";
  if(verbose > 0)
    cout << "file " << ssOutFile.str() << endl << flush;
  gzFile outStream;
  openFile(ssOutFile.str(), outStream, "wb");
  
  // write header line
  ssTxt << "gene" << sep << "snp" << sep << "config";
  for(size_t i = 0; i < iGridL.size(); ++i)
    ssTxt << sep << "l10abf.grid" << (i+1);
  ssTxt << endl;
  size_t nb_lines = 1;
  gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
  
  // write results
  ssTxt.precision(6);
  ssTxt.setf(ios::scientific);
  for(map<string,Gene>::const_iterator it_gene = gene2object.begin();
       it_gene != gene2object.end(); ++it_gene){
    for(vector<GeneSnpPair>::const_iterator it_pair
	   = it_gene->second.BeginPair(); it_pair != it_gene->second.EndPair();
	 ++it_pair){
      
      // write gen BFs (large grid)
      ssTxt.str("");
      ssTxt << it_gene->first
	    << sep << it_pair->GetSnpName()
	    << sep << "gen";
      for(vector<double>::const_iterator it
	    = it_pair->BeginUnweightedAbf("gen");
	  it != it_pair->EndUnweightedAbf("gen"); ++it)
	ssTxt << sep << *it;
      ssTxt << "\n";
      ++nb_lines;
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
      
      // write gen-fix BFs (large grid)
      ssTxt.str("");
      ssTxt << it_gene->first
	    << sep << it_pair->GetSnpName()
	    << sep << "gen-fix";
      for(vector<double>::const_iterator it
	    = it_pair->BeginUnweightedAbf("gen-fix");
	  it != it_pair->EndUnweightedAbf("gen-fix"); ++it)
	ssTxt << sep << *it;
      ssTxt << "\n";
      ++nb_lines;
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
      
      // write gen-maxh BFs (large grid)
      ssTxt.str("");
      ssTxt << it_gene->first
	    << sep << it_pair->GetSnpName()
	    << sep << "gen-maxh";
      for(vector<double>::const_iterator it
	    = it_pair->BeginUnweightedAbf("gen-maxh");
	  it != it_pair->EndUnweightedAbf("gen-maxh"); ++it)
	ssTxt << sep << *it;
      ssTxt << "\n";
      ++nb_lines;
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
      
      // write the BFs for each config (small grid)
      if(bfs != "gen"){
      	for(size_t k = 1; k <= nb_subgroups; ++k){
      	  comb = gsl_combination_calloc(nb_subgroups, k);
      	  if(comb == NULL){
      	    cerr << "ERROR: can't allocate memory for the combination"
      		 << endl;
      	    exit(EXIT_FAILURE);
      	  }
      	  while(true){
      	    ssTxt.str("");
      	    ssTxt << it_gene->first
      		  << sep << it_pair->GetSnpName();
      	    ssConfig.str("");
      	    ssConfig << gsl_combination_get(comb, 0) + 1;
      	    if(comb->k > 1)
      	      for(size_t i = 1; i < k; ++i)
      		ssConfig << "-" << gsl_combination_get(comb, i) + 1;
      	    ssTxt << sep << ssConfig.str();
      	    for(size_t j = 0; j < iGridL.size(); ++j){
      	      if(j < iGridS.size())
      		ssTxt << sep << *(it_pair->BeginUnweightedAbf(ssConfig.str()) + j);
      	      else
      		ssTxt << sep << NaN;
      	    }
      	    ssTxt << "\n";
      	    ++nb_lines;
      	    gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
      	    if(gsl_combination_next(comb) != GSL_SUCCESS)
      	      break;
      	  }
      	  if(bfs == "sin")
      	    break;
      	  gsl_combination_free(comb);
      	}
      }
    }
  }
  
  closeFile(ssOutFile.str(), outStream);
}

void writeResAbfsAvgGrids(
  const string & out_prefix,
  const map<string, Gene> & gene2object,
  const size_t & nb_subgroups,
  const string & bfs,
  const string & error_model,
  const int & verbose)
{
  if(verbose > 0)
    cout << "write results of Bayes Factors (one per configuration) ..."
	 << endl << flush;
  
  string sep = "\t";
  
  gsl_combination * comb;
  stringstream ssOutFile, ssConfig, ssTxt;
  ssOutFile << out_prefix << "_l10abfs_avg-grids.txt.gz";
  if(verbose > 0)
    cout << "file " << ssOutFile.str() << endl << flush;
  gzFile outStream;
  openFile(ssOutFile.str(), outStream, "wb");
  
  // write header line
  ssTxt << "gene" << sep << "snp" << sep << "nb.subgroups";
  ssTxt << sep << "l10abf.gen"
	<< sep << "l10abf.gen.fix"
	<< sep << "l10abf.gen.maxh";
  if(bfs == "sin" || bfs == "all")
    ssTxt << sep << "l10abf.gen.sin";
  if(bfs == "all")
    ssTxt << sep << "l10abf.all";
  if(bfs != "gen"){
    for(size_t k = 1; k <= nb_subgroups; ++k){
      comb = gsl_combination_calloc(nb_subgroups, k);
      if(comb == NULL){
	cerr << "ERROR: can't allocate memory for the combination" << endl;
	exit(EXIT_FAILURE);
      }
      while(true){
	ssTxt << sep << "l10abf." << gsl_combination_get(comb, 0) + 1;
	if(comb->k > 1)
	  for(size_t i = 1; i < k; ++i)
	    ssTxt << "-" << gsl_combination_get(comb, i) + 1;
	if(gsl_combination_next(comb) != GSL_SUCCESS)
	  break;
      }
      gsl_combination_free(comb);
      if(bfs == "sin")
	break;
    }
  }
  ssTxt << endl;
  size_t nb_lines = 1;
  gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
  
  // write results
  ssTxt.precision(6);
  ssTxt.setf(ios::scientific);
  for(map<string, Gene>::const_iterator it_gene = gene2object.begin();
       it_gene != gene2object.end(); ++it_gene){
    for(vector<GeneSnpPair>::const_iterator it_pair
	   = it_gene->second.BeginPair(); it_pair != it_gene->second.EndPair();
	 ++it_pair){
      ssTxt.str("");
      ssTxt << it_gene->first
	    << sep << it_pair->GetSnpName()
	    << sep << it_pair->GetNbSubgroups()
	    << sep << it_pair->GetWeightedAbf("gen")
	    << sep << it_pair->GetWeightedAbf("gen-fix")
	    << sep << it_pair->GetWeightedAbf("gen-maxh");
      if(bfs == "sin" || bfs == "all")
	ssTxt << sep << it_pair->GetWeightedAbf("gen-sin");
      if(bfs == "all")
	ssTxt << sep << it_pair->GetWeightedAbf("all");
      if(bfs != "gen"){
	for(size_t k = 1; k <= nb_subgroups; ++k){
	  comb = gsl_combination_calloc(nb_subgroups, k);
	  if(comb == NULL){
	    cerr << "ERROR: can't allocate memory for the combination"
		 << endl;
	    exit(EXIT_FAILURE);
	  }
	  while(true){
	    ssConfig.str("");
	    ssConfig << gsl_combination_get(comb, 0) + 1;
	    if(comb->k > 1)
	      for(size_t i = 1; i < k; ++i)
		ssConfig << "-" << gsl_combination_get(comb, i) + 1;
	    ssTxt << sep << it_pair->GetWeightedAbf(ssConfig.str());
	    if(gsl_combination_next(comb) != GSL_SUCCESS)
	      break;
	  }
	  gsl_combination_free(comb);
	  if(bfs == "sin")
	    break;
	}
      }
      ssTxt << "\n";
      ++nb_lines;
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
    }
  }
  
  closeFile (ssOutFile.str(), outStream);
}

void writeResJoinPermPval(
  const string & out_prefix,
  const map<string,Gene> & gene2object,
  const size_t & seed,
  const string & permbf,
  const bool & use_max_bf,
  const int & verbose)
{
  if(verbose > 0)
    cout << "write results of gene-level P-values, all subgroups jointly ..."
	 << endl << flush;
  
  string sep = "\t";
  
  stringstream ssOutFile, ssTxt;
  ssOutFile << out_prefix << "_joinPermPvals.txt.gz";
  if(verbose > 0)
    cout << "file " << ssOutFile.str() << endl << flush;
  gzFile outStream;
  openFile (ssOutFile.str(), outStream, "wb");
  size_t nb_lines = 0;
  
  ssTxt << "# perm.bf=" << permbf << " seed=" << seed << endl;
  ++nb_lines;
  gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
  
  ssTxt.str("");
  ssTxt << "gene" << sep << "nb.snps" << sep << "join.perm.pval"
	<< sep << "nb.permutations" << sep << "true.l10abf"
	<< sep << "med.perm.l10abf" << endl;
  ++nb_lines;
  gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
  
  for(map<string,Gene>::const_iterator it_gene = gene2object.begin();
      it_gene != gene2object.end(); ++it_gene){
    if(it_gene->second.GetNbGeneSnpPairs() > 0){
      ssTxt.str("");
      ++nb_lines;
      ssTxt << it_gene->second.GetName()
	    << sep << it_gene->second.GetNbGeneSnpPairs()
	    << sep << it_gene->second.GetPermutationPvalueJoin()
	    << sep << it_gene->second.GetNbPermutationsJoin();
      ssTxt.setf(ios::scientific);
      ssTxt << setprecision(6) << sep << it_gene->second.GetTrueL10Abf(use_max_bf)
	    << sep << it_gene->second.GetMedianPermL10Abf();
      ssTxt.unsetf(ios::scientific);
      ssTxt << "\n";
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
    }
  }
  
  closeFile(ssOutFile.str(), outStream);
}

void writeRes(
  const string & out_prefix,
  const bool & save_sstats,
  const bool & save_weighted_abfs,
  const vector<string> & subgroups,
  const map<string, Gene> & gene2object,
  const map<string, Snp> & snp2object,
  const string & analysis,
  const Grid & iGridL,
  const Grid & iGridS,
  const string & bfs,
  const size_t & nb_permutations,
  const int & perm_sep,
  const string & error_model,
  const size_t & seed,
  const string & permbf,
  const bool & use_max_bf,
  const int & verbose)
{
  if(analysis == "sep" ||
     (analysis == "join" && save_sstats &&
      error_model != "mvlr"))
    writeResSstats(out_prefix, subgroups, gene2object, snp2object, verbose);
  
  if(analysis == "sep" && nb_permutations > 0 &&
     perm_sep != 0){
    if(perm_sep == 1)
      writeResSepPermPval(out_prefix, gene2object, seed, verbose);
    else if(perm_sep == 2)
      writeResSepPermPval(out_prefix, gene2object, subgroups, seed, verbose);;
  }
  
  if(analysis == "join"){
    writeResAbfsRaw(out_prefix, gene2object, subgroups.size(), iGridL, iGridS,
		    bfs, verbose);
    if(save_weighted_abfs)
      writeResAbfsAvgGrids(out_prefix, gene2object, subgroups.size(), bfs,
			   error_model, verbose);
  }
  
  if(analysis == "join" && nb_permutations > 0)
    writeResJoinPermPval(out_prefix, gene2object, seed, permbf,
			 use_max_bf, verbose);
}

void run(const string & file_genopaths,
	 const string & file_snpcoords,
	 const string & file_exppaths,
	 const string & file_genecoords,
	 const string & anchor,
	 const size_t & radius,
	 const string & file_sstats,
	 const string & out_prefix,
	 const bool & save_sstats,
	 const bool & save_weighted_abfs,
	 const string & likelihood,
	 const string & analysis,
	 const bool & need_qnorm,
	 const float & min_maf,
	 const string & file_covarpaths,
	 const string & file_largegrid,
	 const string & file_smallgrid,
	 const string & bfs,
	 const string & error_model,
	 const float & prop_cov_errors,
	 const size_t & nb_permutations,
	 const size_t & seed,
	 const int & trick,
	 const size_t & trick_cutoff,
	 const int & perm_sep,
	 const string & permbf,
	 const bool & use_max_bf,
	 const int & nb_threads,
	 const string & file_snpstokeep,
	 const vector<string> & subgroups_tokeep,
	 const int & verbose)
{
  set<string> sSnpsToKeep;
  if(! file_snpstokeep.empty()){
    loadSnpsToKeep(file_snpstokeep, verbose, sSnpsToKeep);
    if(sSnpsToKeep.empty())
      return;
  }
  
  vector<string> subgroups;
  Samples samples;
  map<string,Snp> snp2object;
  map<string,vector<Snp*> > mChr2VecPtSnps;
  Covariates covariates;
  map<string,Gene> gene2object;
  
  if(file_sstats.empty())
    loadRawInputData(file_genopaths, file_snpcoords, file_exppaths,
		     file_genecoords, anchor, radius, min_maf, file_covarpaths,
		     error_model, subgroups_tokeep, sSnpsToKeep, verbose,
		     subgroups, samples, snp2object, mChr2VecPtSnps,
		     covariates, gene2object);
  else
    loadSummaryStats(file_sstats, verbose, subgroups, gene2object, snp2object);
  if(gene2object.empty() || snp2object.empty())
    return;
  
  Grid iGridL(file_largegrid, true, verbose);
  Grid iGridS(file_smallgrid, false, verbose);
  
  testForAssociations(file_sstats.empty(), mChr2VecPtSnps, anchor, radius,
		      subgroups, samples, likelihood, analysis,
		      need_qnorm, covariates, iGridL, iGridS, bfs,
		      error_model, prop_cov_errors, verbose, gene2object);
  if(nb_permutations > 0 &&
     (perm_sep != 0 || permbf != "none")){
    omp_set_num_threads(nb_threads);
    makePermutations(subgroups, samples, likelihood, analysis, need_qnorm,
		     covariates, iGridL, iGridS, error_model, prop_cov_errors,
		     nb_permutations, seed, trick, trick_cutoff, perm_sep,
		     permbf, use_max_bf, verbose, gene2object);
  }
  
  writeRes(out_prefix, save_sstats, save_weighted_abfs, subgroups, gene2object,
	   snp2object, analysis, iGridL, iGridS, bfs, nb_permutations,
	   perm_sep, error_model, seed, permbf, use_max_bf, verbose);
}

int main(int argc, char ** argv)
{
  int verbose = 1, trick = 0, perm_sep = 0, nb_threads = 1;
  size_t radius = 100000, nb_types = string::npos, nb_permutations = 0,
    seed = string::npos, trick_cutoff = 10;
  float min_maf = 0.0, prop_cov_errors = 0.5;
  bool save_sstats = false, save_weighted_abfs = false, need_qnorm = false,
    use_max_bf = false;
  string file_genopaths, file_snpcoords, file_exppaths, file_genecoords,
    anchor = "TSS", file_sstats, out_prefix, likelihood = "normal",
    analysis, file_covarpaths, file_largegrid, file_smallgrid,
    bfs = "gen", error_model = "uvlr", permbf = "none",
    file_snpstokeep;
  vector<string> subgroups_tokeep;
  
  parseCmdLine(argc, argv, file_genopaths, file_snpcoords, file_exppaths,
	       file_genecoords, anchor, radius, file_sstats, out_prefix,
	       save_sstats, save_weighted_abfs, likelihood, analysis,
	       need_qnorm, min_maf, file_covarpaths, file_largegrid,
	       file_smallgrid, bfs, error_model, prop_cov_errors, nb_types,
	       nb_permutations, seed, trick, trick_cutoff, perm_sep,
	       permbf, use_max_bf, nb_threads, file_snpstokeep,
	       subgroups_tokeep, verbose);
  
  time_t time_start, time_end;
  if(verbose > 0){
    time(&time_start);
    cout << "START " << basename(argv[0])
	 << " " << getDateTime(time_start) << endl
	 << "version " << VERSION << " compiled " << __DATE__
	 << " " << __TIME__ << endl
	 << "cmd-line: " << getCmdLine(argc, argv) << endl
	 << "cwd: " << getCurrentDirectory() << endl << flush;
  }
  
  run(file_genopaths, file_snpcoords, file_exppaths, file_genecoords, anchor,
      radius, file_sstats, out_prefix, save_sstats, save_weighted_abfs,
      likelihood, analysis, need_qnorm, min_maf, file_covarpaths,
      file_largegrid, file_smallgrid, bfs, error_model, prop_cov_errors,
      nb_permutations, seed, trick, trick_cutoff, perm_sep, permbf,
      use_max_bf, nb_threads, file_snpstokeep, subgroups_tokeep, verbose);
  
  if(verbose > 0){
    time (&time_end);
    cout << "END " << basename(argv[0])
	 << " " << getDateTime(time_end) << endl
	 << "elapsed -> " << getElapsedTime(time_start, time_end) << endl
	 << "max.mem -> " << getMaxMemUsedByProcess2Str() << endl;
  }
  
  return EXIT_SUCCESS;
}
