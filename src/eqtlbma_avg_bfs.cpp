/** \file eqtlbma_avg_bfs.cpp
 *
 *  `eqtlbma_avg_bfs' averages the raw BFs, over the grid only, or over both
 *  the grid and the configurations, and can also compute posteriors.
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

#include <ctime>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <getopt.h>
#include <cmath>
#include <libgen.h>

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
#include <gsl/gsl_cdf.h>

#include <omp.h>

#include "utils/utils_io.hpp"
#include "utils/utils_math.hpp"
using namespace utils;

#ifndef VERSION
#define VERSION "0.0"
#endif

/** \brief Display the help on stdout.
 */
void help(char ** argv)
{
  cout << "`" << argv[0] << "'"
       << " averages the raw BFs over the grid only," << endl
       << "or over both the grid and the configurations," << endl
       << "and it can also compute posteriors." << endl
       << endl
       << "Usage: " << argv[0] << " [OPTIONS] ..." << endl
       << endl
       << "Options:" << endl
       << "  -h, --help\tdisplay the help and exit" << endl
       << "  -V, --version\toutput version information and exit" << endl
       << "  -v, --verbose\tverbosity level (0/default=1/2/3)" << endl
       << "      --in\tpattern to glob '_l10abfs_raw' files from 'eqtlbma_bf'" << endl
       << "      --gwts\tfile with grid weights (one per line, only the value)" << endl
       << "      --gtk\tind-ex/icies of grid weights to keep (all by default)" << endl
       << "\t\te.g. '1+3+5+7+9' to keep only those with no heterogeneity" << endl
       << "      --model\twhich model (default=configs/types)" << endl
       << "      --nsubgrp\tnumber of subgroups" << endl
       << "      --dim\tdimension of the model (nb of configs or types)" << endl
       << "      --cwts\tfile with configuration weights (one per line, name<sep>value)" << endl
       << "\t\tonly a subset of the configs can be given, in agreement with --nsubgrp and --dim" << endl
       << "      --tswts\tfile with type and subgroup weights (one per line, name<sep>value)" << endl
       << "      --save\tprecise what to save (bf/post/bf+post)" << endl
       << "\t\t'post' requires also options --pi0 and --post" << endl
       << "      --pi0\tproba for a gene to have no eQTL in any subgroup" << endl
       << "\t\tif not provided, BFs will be saved instead of posterior probability" << endl
       << "      --post\tsave various kinds of posterior probabilities (e.g. 'a+b')" << endl
       << "\t\ta: the gene has at least one eQTL" << endl
       << "\t\tb: the SNP is 'the' eQTL for the gene, in at least one subgroup, given that the gene has exactly one eQTL,\n\t\tassuming all cis SNPs are equally likely and a single eQTL per gene" << endl
       << "\t\tc: the SNP is 'an' eQTL for the gene, in at least one subgroup, given that the gene contains at least one eQTL\n\t\tand that SNPs are independent" << endl
       << "\t\td: the SNP is an eQTL in subgroup s, given that it is 'the' eQTL for the gene, the configs/types being marginalized" << endl
       << "      --gene\tfile with subset of gene(s) to keep (one per line)" << endl
       << "      --snp\tfile with subset of snp(s) to keep (one per line)" << endl
       << "\t\tcaution, it can change the gene-level BFs and posteriors" << endl
       << "      --gene-snp\tfile with subset of gene-snp pai(s) to keep (gene<tab>snp, one per line)" << endl
       << "\t\tcaution, it can change the gene-level BFs and posteriors" << endl
       << "      --bestsnp\treport the best SNP(s) per gene" << endl
       << "\t\t0: report all SNPs (default)" << endl
       << "\t\t1: report only the single best SNP (pick one if tie)" << endl
       << "\t\t2: report the best SNP(s) listed in decreasing order of their probability of being the eQTL (conditional on the gene containing an eQTL), such that the sum of these probabilities exceeds 0.95" << endl
       << "      --bestdim\treport the best config/type per SNP (and its posterior)" << endl
       << "      --alldim\treport also BF and/or posterior for all dimensions (configs or types)" << endl
       << "\t\tcaution, the number of configurations can be big" << endl
       << "      --out\tname of the output file (gzipped)" << endl
       << "\t\tif --cwts is not provided, the output file will be used as input for 'eqtlbma_hm'" << endl
       << "      --thread\tnumber of threads (default=1)" << endl
       << endl;
}

/** \brief Display version and license information on stdout.
 */
void version(char ** argv)
{
  cout << argv[0] << " " << VERSION << endl
       << endl
       << "Copyright (C) 2013 Timothee Flutre." << endl
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
parseCmdLine(
  int argc,
  char ** argv,
  int & verbose,
  string & file_pattern,
  string & file_grid_weights,
  vector<size_t> & grid_idx_to_keep,
  string & model,
  size_t & nb_subgroups,
  size_t & dim,
  string & file_config_weights,
  string & file_type_subgroup_weights,
  vector<string> & quantities_to_save,
  double & pi0,
  vector<string> & post_probas_to_save,
  string & file_genes_to_keep,
  string & file_snps_to_keep,
  string & file_gene_snp_pairs_to_keep,
  int & save_best_snps,
  bool & save_best_dim,
  bool & save_all_dims,
  string & file_hm,
  int & nb_threads)
{
  int c = 0;
  while(true){
    static struct option long_options[] = {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'V'},
      {"verbose", required_argument, 0, 'v'},
      {"in", required_argument, 0, 0},
      {"gwts", required_argument, 0, 0},
      {"gtk", required_argument, 0, 0},
      {"model", required_argument, 0, 0},
      {"nsubgrp", required_argument, 0, 0},
      {"dim", required_argument, 0, 0},
      {"cwts", required_argument, 0, 0},
      {"tswts", required_argument, 0, 0},
      {"save", required_argument, 0, 0},
      {"pi0", required_argument, 0, 0},
      {"post", required_argument, 0, 0},
      {"gene", required_argument, 0, 0},
      {"snp", required_argument, 0, 0},
      {"gene-snp", required_argument, 0, 0},
      {"bestsnp", required_argument, 0, 0},
      {"bestdim", no_argument, 0, 0},
      {"alldim", no_argument, 0, 0},
      {"out", required_argument, 0, 0},
      {"thread", required_argument, 0, 0},
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
      if(strcmp(long_options[option_index].name, "in") == 0){
	file_pattern = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "gwts") == 0){
	file_grid_weights = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "gtk") == 0){
	vector<string> tmp;
	split(optarg, "+", tmp);
	for(size_t i = 0; i < tmp.size(); ++i)
	  grid_idx_to_keep.push_back(atoi(tmp[i].c_str()));
	break;
      }
      if(strcmp(long_options[option_index].name, "model") == 0){
	model = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "nsubgrp") == 0){
	nb_subgroups = atol(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "dim") == 0){
	dim = atoi(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "cwts") == 0){
	file_config_weights = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "tswts") == 0){
	file_type_subgroup_weights = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "save") == 0){
	split(optarg, "+", quantities_to_save);
	break;
      }
      if(strcmp(long_options[option_index].name, "pi0") == 0){
	pi0 = atof(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "post") == 0){
	split(optarg, "+", post_probas_to_save);
	break;
      }
      if(strcmp(long_options[option_index].name, "gene") == 0){
	file_genes_to_keep = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "snp") == 0){
	file_snps_to_keep = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "gene-snp") == 0){
	file_gene_snp_pairs_to_keep = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "bestsnp") == 0){
	save_best_snps = atoi(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "bestdim") == 0){
	save_best_dim = true;
	break;
      }
      if(strcmp(long_options[option_index].name, "alldim") == 0){
	save_all_dims = true;
	break;
      }
      if(strcmp(long_options[option_index].name, "out") == 0){
	file_hm = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "thread") == 0){
	nb_threads = atoi(optarg);
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
  if(file_pattern.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --bf-out" << endl << endl;
    help(argv);
    exit(1);
  }
  if(file_grid_weights.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --gwts" << endl << endl;
    help(argv);
    exit(1);
  }
  if(! doesFileExist(file_grid_weights)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find file " << file_grid_weights << endl << endl;
    help(argv);
    exit(1);
  }
  if(! model.empty() && model != "configs" && model != "types"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --model " << model << " is unvalid" << endl << endl;
    help(argv);
    exit(1);
  }
  if(! model.empty() && nb_subgroups == string::npos){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --nsubgrp" << endl << endl;
    help(argv);
    exit(1);
  }
  if(! model.empty() && dim == string::npos){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --dim" << endl << endl;
    help(argv);
    exit(1);
  }
  if(! file_config_weights.empty() && ! doesFileExist(file_config_weights)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find file " << file_config_weights << endl << endl;
    help(argv);
    exit(1);
  }
  if(! file_type_subgroup_weights.empty() && ! doesFileExist(file_type_subgroup_weights)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find file " << file_type_subgroup_weights << endl << endl;
    help(argv);
    exit(1);
  }
  if(file_hm.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --hm-in" << endl << endl;
    help(argv);
    exit(1);
  }
  if(find(quantities_to_save.begin(), quantities_to_save.end(), "post")
     != quantities_to_save.end() && post_probas_to_save.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --post with --save post" << endl << endl;
    help(argv);
    exit(1);
  }
  if(! isNan(pi0) && post_probas_to_save.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --post with --pi0" << endl << endl;
    help(argv);
    exit(1);
  }
  if(isNan(pi0) && ! post_probas_to_save.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --pi0 with --post" << endl << endl;
    help(argv);
    exit(1);
  }
  if(nb_threads <= 0){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --thread " << nb_threads << " is invalid" << endl << endl;
    help(argv);
    exit(1);
  }
  if(save_best_snps != 0 && save_best_snps != 1 && save_best_snps != 2){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --bestsnp " << save_best_snps << " is invalid" << endl << endl;
    help(argv);
    exit(1);
  }
  if(save_best_snps == 2
     && find(quantities_to_save.begin(), quantities_to_save.end(), "post")
     == quantities_to_save.end()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --bestsnp 2 can only be used with --save post" << endl << endl;
    help(argv);
    exit(1);
  }
  if(! file_genes_to_keep.empty() && ! doesFileExist(file_genes_to_keep)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find file " << file_genes_to_keep << endl << endl;
    help(argv);
    exit(1);
  }
  if(! file_snps_to_keep.empty() && ! doesFileExist(file_snps_to_keep)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find file " << file_snps_to_keep << endl << endl;
    help(argv);
    exit(1);
  }
  if(! file_snps_to_keep.empty())
    cout << "WARNING: using --snp can change the gene-level BFs and posteriors" << endl;
  if(! file_gene_snp_pairs_to_keep.empty()
     && ! doesFileExist(file_gene_snp_pairs_to_keep)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find file " << file_gene_snp_pairs_to_keep << endl << endl;
    help(argv);
    exit(1);
  }
  if(! file_gene_snp_pairs_to_keep.empty())
    cout << "WARNING: using --gene-snp can change the gene-level BFs and posteriors" << endl;
}

class BayesFactor
{
public:
  string gene;
  string snp;
  string config;
  double log10_val;
  BayesFactor(void);
  BayesFactor(const string & g, const string & s, const string & c);
  void avg_raw_bfs(const vector<string> & tokens,
		   const vector<double> & grid_weights,
		   const vector<size_t> & grid_idx_to_keep);
};

BayesFactor::BayesFactor(void)
{
}

BayesFactor::BayesFactor(const string & g, const string & s, const string & c)
{
  gene = g;
  snp = s;
  config = c;
}

void BayesFactor::avg_raw_bfs(const vector<string> & tokens,
			      const vector<double> & grid_weights,
			      const vector<size_t> & grid_idx_to_keep)
{
  vector<double> log10_values(grid_idx_to_keep.size(), NaN);
  vector<double> weights(grid_idx_to_keep.size(), NaN);
  for(size_t i = 0; i < grid_idx_to_keep.size(); ++i){
    log10_values[i] = atof(tokens[grid_idx_to_keep[i]].c_str());
    weights[i] = grid_weights[grid_idx_to_keep[i]];
  }
  log10_val = log10_weighted_sum(&(log10_values[0]), &(weights[0]),
				 log10_values.size());
}

class Snp
{
public:
  string name_;
  vector<vector<double> > raw_log10_bfs_; // dim1 is config/subgroup; dim2 is grid
  vector<double> dim_log10_bfs_; // raw_log10_bfs avg over grid (and subgroups if "types" model)
  double log10_bf_;
  double post_the_eQTL_;
  double post_an_eQTL_;
  vector<double> post_dims_;
  vector<double> post_subgroups_;
  size_t best_dim_idx_;
  Snp(void);
  Snp(const string & n, const vector<vector<double> > & vv);
  void avg_raw_bfs(const vector<double> & grid_weights,
		   const vector<size_t> & grid_idx_to_keep,
		   const vector<double> & config_weights);
  void avg_raw_bfs(const vector<double> & grid_weights,
		   const vector<size_t> & grid_idx_to_keep,
		   const vector<double> & type_weights,
		   const vector<vector<double> > & subgroup_weights);
};

Snp::Snp(void)
{
}

Snp::Snp(const string & n, const vector<vector<double> > & vv)
{
  name_ = n;
  raw_log10_bfs_ = vv;
  log10_bf_ = NaN;
  best_dim_idx_ = string::npos;
}

void Snp::avg_raw_bfs(const vector<double> & grid_weights,
		      const vector<size_t> & grid_idx_to_keep,
		      const vector<double> & config_weights)
{
  vector<double> log10_grid_values_tmp, grid_weights_tmp;
  dim_log10_bfs_ = vector<double>(config_weights.size(), NaN);
  
  // for each config, average BFs over grid
  for(size_t j = 0; j < config_weights.size(); ++j){
    log10_grid_values_tmp.assign(grid_idx_to_keep.size(), NaN);
    grid_weights_tmp.assign(grid_idx_to_keep.size(), NaN);
    for(size_t l = 0; l < grid_idx_to_keep.size(); ++l){
      log10_grid_values_tmp[l] = raw_log10_bfs_[j][grid_idx_to_keep[l]];
      grid_weights_tmp[l] = grid_weights[grid_idx_to_keep[l]];
    }
    dim_log10_bfs_[j] = log10_weighted_sum(&(log10_grid_values_tmp[0]),
					   &(grid_weights_tmp[0]),
					   log10_grid_values_tmp.size());
    if(best_dim_idx_ == string::npos
       || dim_log10_bfs_[j] > dim_log10_bfs_[best_dim_idx_])
      best_dim_idx_ = j;
  }
  
  // average BFs over configs
  log10_bf_ = log10_weighted_sum(&(dim_log10_bfs_[0]),
				 &(config_weights[0]),
				 dim_log10_bfs_.size());
}

void Snp::avg_raw_bfs(const vector<double> & grid_weights,
		      const vector<size_t> & grid_idx_to_keep,
		      const vector<double> & type_weights,
		      const vector<vector<double> > & subgroup_weights)
{
  vector<double> log10_grid_values_tmp, grid_weights_tmp;
  dim_log10_bfs_ = vector<double>(type_weights.size(), NaN);
  
  // for each type, average BFs over grid and subgroups
  for(size_t k = 0; k < type_weights.size(); ++k){
    log10_grid_values_tmp.assign(grid_idx_to_keep.size(), NaN);
    grid_weights_tmp.assign(grid_idx_to_keep.size(), NaN);
    for(size_t l = 0; l < grid_idx_to_keep.size(); ++l){
      log10_grid_values_tmp[l] = 0.0;
      for(size_t s = 0; s < subgroup_weights[0].size(); ++s)
	log10_grid_values_tmp[l] += log10(subgroup_weights[k][s]
					  * pow(10, raw_log10_bfs_[s][grid_idx_to_keep[l]])
					  + 1 - subgroup_weights[k][s]);
      grid_weights_tmp[l] = grid_weights[grid_idx_to_keep[l]];
    }
    dim_log10_bfs_[k] = log10_weighted_sum(&(log10_grid_values_tmp[0]),
					   &(grid_weights_tmp[0]),
					   log10_grid_values_tmp.size());
  }
  
  // average BFs over types
  log10_bf_ = log10_weighted_sum(&(dim_log10_bfs_[0]),
				 &(type_weights[0]),
				 dim_log10_bfs_.size());
}

bool operator< (const Snp& lhs, const Snp& rhs){
  return(lhs.log10_bf_ < rhs.log10_bf_);
}

bool operator> (const Snp& lhs, const Snp& rhs){
  return operator<(rhs,lhs);
}

bool operator<=(const Snp& lhs, const Snp& rhs){
  return !operator>(lhs,rhs);
}

bool operator>=(const Snp& lhs, const Snp& rhs){
  return !operator<(lhs,rhs);
}


class Gene
{
public:
  string name_;
  vector<Snp> snps_;
  double log10_bf_; // averaged over SNPs
  double log10_obs_lik_; // requires pi0
  double post_; // requires pi0
  size_t idx_last_best_snp_; // once snps_ is sorted wrt their log10_bf_
  Gene(void);
  Gene(const string & n);
  void avg_raw_bfs(const vector<double> & grid_weights,
		   const vector<size_t> & grid_idx_to_keep,
		   const vector<double> & config_weights);
  void avg_raw_bfs(const vector<double> & grid_weights,
		   const vector<size_t> & grid_idx_to_keep,
		   const vector<double> & type_weights,
		   const vector<vector<double> > & subgroup_weights);
  void calc_log10_obs_lik(const double & pi0);
  void calc_posterior(const double & pi0);
  void calc_cond_snp_posteriors_the_eQTL(void);
  void calc_cond_snp_posteriors_an_eQTL(void);
  void calc_cond_snp_posteriors_config(const vector<double> & config_weights);
  void calc_cond_snp_posteriors_type(const vector<double> & type_weights);
  void calc_cond_snp_posteriors_subgroup(
    const vector<string> & config_names,
    const vector<vector<string> > & config2subgroups);
  void calc_cond_snp_posteriors_subgroup(const vector<double> & type_weights,
					 const vector<vector<double> > & subgroup_weights,
					 const vector<double> & grid_weights,
					 const vector<size_t> & grid_idx_to_keep);
  void identify_best_snps(const int & save_best_snps);
};

Gene::Gene(void)
{
}

Gene::Gene(const string & n)
{
  name_ = n;
  log10_bf_ = NaN;
  post_ = NaN;
}

void Gene::avg_raw_bfs(const vector<double> & grid_weights,
		       const vector<size_t> & grid_idx_to_keep,
		       const vector<double> & config_weights)
{
  vector<double> snp_log10_bfs(snps_.size(), NaN);
  for(size_t p = 0; p < snps_.size(); ++p){
    snps_[p].avg_raw_bfs(grid_weights, grid_idx_to_keep, config_weights);
    snp_log10_bfs[p] = snps_[p].log10_bf_;
  }
  log10_bf_ = log10_weighted_sum(&(snp_log10_bfs[0]), snp_log10_bfs.size());
}

void Gene::avg_raw_bfs(const vector<double> & grid_weights,
		       const vector<size_t> & grid_idx_to_keep,
		       const vector<double> & type_weights,
		       const vector<vector<double> > & subgroup_weights)
{
  vector<double> snp_log10_bfs(snps_.size(), NaN);
  for(size_t p = 0; p < snps_.size(); ++p){
#ifdef DEBUG
    fprintf(stderr, "average BFs of %s\n", snps_[p].name_.c_str());
#endif
    snps_[p].avg_raw_bfs(grid_weights, grid_idx_to_keep, type_weights, subgroup_weights);
    snp_log10_bfs[p] = snps_[p].log10_bf_;
  }
  log10_bf_ = log10_weighted_sum(&(snp_log10_bfs[0]), snp_log10_bfs.size());
}

void Gene::calc_log10_obs_lik(const double & pi0)
{
  double vec[2], wts[2];
  wts[0] = pi0;
  wts[1] = 1 - pi0;
  vec[0] = 0;
  vec[1] = log10_bf_;
  log10_obs_lik_ = log10_weighted_sum(vec, wts, 2);
}

// Pr(z_g = 1 | Y, X, Theta)
void Gene::calc_posterior(const double & pi0)
{
  if(! isNan(pi0)){
    calc_log10_obs_lik(pi0);
    post_ = pow(10, log10(1.0-pi0) + log10_bf_ - log10_obs_lik_);
    if(post_ > 1.0)
      post_ = 1.0;
  }
}

// Pr(v_gp = 1 | Y, X, Theta, z_g = 1)
// assuming same prior for all SNPs (1/m_g)
void Gene::calc_cond_snp_posteriors_the_eQTL(void)
{
  vector<double> snps_log10_bf(snps_.size(), NaN),
    snps_one(snps_.size(), 1.0);
  for(size_t p = 0; p < snps_.size(); ++p)
    snps_log10_bf[p] = snps_[p].log10_bf_;
  double log10_sum_snps_bf = log10_weighted_sum(&(snps_log10_bf[0]),
						&(snps_one[0]),
						snps_.size());
  
  for(size_t p = 0; p < snps_.size(); ++p){
    snps_[p].post_the_eQTL_ = pow(10, snps_[p].log10_bf_ - log10_sum_snps_bf);
    if(snps_[p].post_the_eQTL_ > 1.0)
      snps_[p].post_the_eQTL_ = 1.0;
  }
}

void Gene::calc_cond_snp_posteriors_an_eQTL(void)
{
  double vec[2], wts[2];
  wts[0] = 1.0 - 1/(double)snps_.size();
  wts[1] = 1/(double)snps_.size();
  vec[0] = 0;
  
  for(size_t p = 0; p < snps_.size(); ++p){
    vec[1] = snps_[p].log10_bf_;
    snps_[p].post_an_eQTL_ = pow(10, -log10((double)snps_.size())
				 + snps_[p].log10_bf_
				 - log10_weighted_sum(vec, wts, 2));
    if(snps_[p].post_the_eQTL_ > 1.0)
      snps_[p].post_an_eQTL_ = 1.0;
  }
}

// Pr(c_gpj = 1 | Y, X, Theta, z_g = 1, v_gp = 1)
void Gene::calc_cond_snp_posteriors_config(const vector<double> & config_weights)
{
  for(size_t p = 0; p < snps_.size(); ++p){
    snps_[p].post_dims_ = vector<double>(config_weights.size(), NaN);
    for(size_t j = 0; j < config_weights.size(); ++j){
      snps_[p].post_dims_[j] = pow(10, log10(config_weights[j])
				   + snps_[p].dim_log10_bfs_[j]
				   - snps_[p].log10_bf_);
      if(snps_[p].post_dims_[j] > 1.0)
	snps_[p].post_dims_[j] = 1.0;
      if(j == 0
	 || snps_[p].post_dims_[j] > snps_[p].post_dims_[snps_[p].best_dim_idx_])
	snps_[p].best_dim_idx_ = j;
    }
  }
}

// Pr(t_gpk = 1 | Y, X, Theta, z_g = 1, v_gp = 1)
void Gene::calc_cond_snp_posteriors_type(const vector<double> & type_weights)
{
  for(size_t p = 0; p < snps_.size(); ++p){
    snps_[p].post_dims_ = vector<double>(type_weights.size(), NaN);
    for(size_t k = 0; k < type_weights.size(); ++k){
      snps_[p].post_dims_[k] = pow(10, log10(type_weights[k])
				   + snps_[p].dim_log10_bfs_[k]
				   - snps_[p].log10_bf_);
      if(snps_[p].post_dims_[k] > 1.0)
	snps_[p].post_dims_[k] = 1.0;
      if(k == 0
	 || snps_[p].post_dims_[k] > snps_[p].post_dims_[snps_[p].best_dim_idx_])
	snps_[p].best_dim_idx_ = k;
    }
  }
}

// Pr(gamma_gps = 1 | Y, X, Theta, z_g = 1, v_gp = 1)
// averaging over configs
void Gene::calc_cond_snp_posteriors_subgroup(
  const vector<string> & config_names,
  const vector<vector<string> > & config2subgroups)
{
  size_t nb_subgroups = (size_t) log2(config_names.size() + 1);
  
  stringstream subgroup_id;
  for(size_t p = 0; p < snps_.size(); ++p){
    snps_[p].post_subgroups_ = vector<double>(nb_subgroups, NaN);
    
    for(size_t s = 0; s < nb_subgroups; ++s){
      subgroup_id.str("");
      subgroup_id << s+1;
      snps_[p].post_subgroups_[s] = 0.0;
      for(size_t j = 0; j < config_names.size(); ++j)
      	if(find(config2subgroups[j].begin(), config2subgroups[j].end(), subgroup_id.str())
      	   != config2subgroups[j].end())
      	  snps_[p].post_subgroups_[s] += snps_[p].post_dims_[j];
      if(snps_[p].post_subgroups_[s] > 1.0)
	snps_[p].post_subgroups_[s] = 1.0;
    } // end of "for each subgroup"
  } // end of "for each SNP"
}

// Pr(gamma_gps = 1 | Y, X, Theta, z_g = 1, v_gp = 1)
// averaging over types
void Gene::calc_cond_snp_posteriors_subgroup(const vector<double> & type_weights,
					     const vector<vector<double> > & subgroup_weights,
					     const vector<double> & grid_weights,
					     const vector<size_t> & grid_idx_to_keep)
{
  size_t nb_subgroups = subgroup_weights[0].size(),
    dim = type_weights.size(),
    grid_size = grid_idx_to_keep.size();
  vector<double> log10_grid_values_tmp(grid_size, NaN),
    grid_weights_tmp(grid_size, NaN);
  
  for(size_t p = 0; p < snps_.size(); ++p){
    snps_[p].post_subgroups_.assign(nb_subgroups, NaN);
    
    for(size_t s = 0; s < nb_subgroups; ++s){
      snps_[p].post_subgroups_[s] = 0.0;
      
      for(size_t k = 0; k < dim; ++k){
	double log10_type_k_tmp = log10(type_weights[k])
	  + log10(subgroup_weights[k][s])
	  - snps_[p].log10_bf_;
	
	for(size_t l = 0; l < grid_size; ++l){
	  log10_grid_values_tmp[l] = snps_[p].raw_log10_bfs_[s][grid_idx_to_keep[l]];
	  for(size_t s2 = 0; s2 < nb_subgroups; ++s2)
	    if(s2 != s)
	      log10_grid_values_tmp[l] +=
		log10(subgroup_weights[k][s2]
		      * pow(10, snps_[p].raw_log10_bfs_[s2][grid_idx_to_keep[l]])
		      + 1 - subgroup_weights[k][s2]);
	  grid_weights_tmp[l] = grid_weights[grid_idx_to_keep[l]];
	}
	
	log10_type_k_tmp +=
	  log10_weighted_sum(&(log10_grid_values_tmp[0]),
			     &(grid_weights_tmp[0]),
			     grid_size);
	
	snps_[p].post_subgroups_[s] += pow(10, log10_type_k_tmp);
      }
      
      if(snps_[p].post_subgroups_[s] > 1.0)
	snps_[p].post_subgroups_[s] = 1.0;
      
    } // end of "for each subgroup"
  } // end of "for each SNP"
}

void Gene::identify_best_snps(const int & save_best_snps){
  if(save_best_snps == 0)
    idx_last_best_snp_ = snps_.size() - 1;
  else{
    // sort vector "snps" decreasingly w.r.t. their log10(BF)
    sort(snps_.rbegin(), snps_.rend());
    
    if(save_best_snps == 1)
      idx_last_best_snp_ = 0;
    else{
      // find the min nb of SNP(s) such that sum of their proba 
      // of being the eQTL exceeds 0.95
      double cum_sum = 0.0;
      for(size_t p = 0; p < snps_.size(); ++p){
	cum_sum += snps_[p].post_the_eQTL_;
	if(cum_sum >= 0.95){
	  idx_last_best_snp_ = p;
	  break;
	}
      }
    }
  }
}


void loadFileGridWeights(const string & file_grid_weights,
			 const int & verbose,
			 vector<double> & grid_weights,
			 vector<size_t> & grid_idx_to_keep)
{
  if(verbose > 0)
    cout <<"load grid weights ..." << endl << flush;
  
  string line;
  gzFile stream;
  vector<string> tokens;
  size_t nb_lines = 0;
  
  openFile(file_grid_weights, stream, "rb");
  
  while(getline(stream, line)){
    nb_lines++;
    split(line, " \t", tokens);
    grid_weights.push_back(atof(tokens[0].c_str()));
  }
  if(! gzeof(stream)){
    cerr << "ERROR: can't read successfully file "
	 << file_grid_weights << " up to the end" << endl;
    exit(1);
  }
  
  closeFile(file_grid_weights, stream);
  
  if(verbose > 0)
    for(size_t i = 0; i < grid_weights.size(); ++i)
      cout << "grid weight " << i+1 << ": "
	   << setprecision(4) << scientific
	   << grid_weights[i] << endl;
  
  if(grid_idx_to_keep.empty()){
    if(verbose > 0)
      cout << "use all grid weights" << endl;
    for(size_t i = 0; i < grid_weights.size(); ++i)
      grid_idx_to_keep.push_back(i);
  }
  else{
    if(*(max_element(grid_idx_to_keep.begin(), grid_idx_to_keep.end()))
       >= grid_weights.size()){
      cerr << "ERROR: --gtk doesn't correspond to --gwts" << endl;
      exit(1);
    }
    if(verbose > 0)
      cout << "use only " << grid_idx_to_keep.size() << " grid weights" << endl;
  }
}

void loadFileConfigWeights(const string & file_config_weights,
			   const size_t & nb_subgroups,
			   const size_t & dim,
			   const int & verbose,
			   vector<string> & config_names,
			   vector<double> & config_weights,
			   vector<vector<string> > & config2subgroups)
{
  if(! file_config_weights.empty()){
    if(verbose > 0)
      cout <<"load config weights ..." << endl << flush;
    
    string line;
    gzFile stream;
    vector<string> tokens;
    size_t nb_lines = 0;
    
    openFile(file_config_weights, stream, "rb");
    
    while(getline(stream, line)){
      nb_lines++;
      split(line, " \t,", tokens);
      config_names.push_back(tokens[0]);
      config_weights.push_back(atof(tokens[1].c_str()));
      config2subgroups.push_back(vector<string>());
      split(tokens[0], '-', config2subgroups[config2subgroups.size()-1]);
    }
    if(! gzeof(stream)){
      cerr << "ERROR: can't read successfully file "
	   << file_config_weights << " up to the end" << endl;
      exit(1);
    }
    
    closeFile(file_config_weights, stream);
    
    if(config_weights.size() != dim){
      cerr << "ERROR: file " << file_config_weights << " should contain "
	   << dim << " config weights" << endl;
      exit(1);
    }
    
    if(verbose > 0)
      for(size_t k = 0; k < dim; ++k)
	cout << "config weight " << k+1 << " (" << config_names[k] << "): "
	     << setprecision(4) << scientific
	     << config_weights[k] << endl;
  }
  else{ // if "types" model
    for(size_t s = 0; s < nb_subgroups; ++s)
      config_names.push_back(toString(s+1));
  }
}

void loadFileTypeSubgroupWeights(const string & file_type_subgroup_weights,
				 const size_t & nb_subgroups,
				 const size_t & dim,
				 const int & verbose,
				 vector<double> & type_weights,
				 vector<vector<double> > & subgroup_weights)
{
  if(! file_type_subgroup_weights.empty()){
    if(verbose > 0)
      cout <<"load type and subgroup weights ..." << endl << flush;
    
    string line;
    gzFile stream;
    vector<string> tokens;
    size_t nb_lines = 0;
    
    openFile(file_type_subgroup_weights, stream, "rb");
    
    size_t idx_type = 0, idx_subgrp = 0;
    while(getline(stream, line)){
      nb_lines++;
      split(line, " \t,", tokens);
      if(tokens[0].find("type") != string::npos)
	type_weights.push_back(atof(tokens[1].c_str()));
      else if(tokens[0].find("subgroup") != string::npos){
	if(idx_subgrp == 0)
	  subgroup_weights.push_back(vector<double>());
	subgroup_weights[idx_type].push_back(atof(tokens[1].c_str()));
	++idx_subgrp;
	if(idx_subgrp == nb_subgroups){
	  ++idx_type;
	  idx_subgrp = 0;
	}
      }
    }
    if(! gzeof(stream)){
      cerr << "ERROR: can't read successfully file "
	   << file_type_subgroup_weights << " up to the end" << endl;
      exit(1);
    }
    
    closeFile(file_type_subgroup_weights, stream);
    
    if(type_weights.size() != dim){
      cerr << "ERROR: file " << file_type_subgroup_weights << " should contain "
	   << dim << " type weights" << endl;
      exit(1);
    }
    if(subgroup_weights.size() != dim){
      cerr << "ERROR: file " << file_type_subgroup_weights << " should contain "
	   << dim << " types for subgroups-per-type weights" << endl;
      exit(1);
    }
    for(size_t k = 0; k < dim; ++k)
      if(subgroup_weights[k].size() != nb_subgroups){
	cerr << "ERROR: file " << file_type_subgroup_weights << " should contain "
	     << nb_subgroups << " subgroups for subgroups-per-type "
	     << k+1 << " weights" << endl;
	exit(1);
      }
    
    if(verbose > 0)
      for(size_t k = 0; k < dim; ++k){
	cout << "type weight " << k+1 << ": "
	     << setprecision(4) << scientific
	     << type_weights[k] << endl;
	for(size_t s = 0; s < nb_subgroups; ++s)
	  cout << "subgroup " << s+1 << " per type " << k+1 << ": "
	       << subgroup_weights[k][s] << endl;
      }
  }
}

void loadFileGeneSnpPairsToKeep(
  const string & file_gene_snp_pairs_to_keep,
  const int & verbose,
  map<string,vector<string> > & gene_snp_pairs_to_keep)
{
  if(! file_gene_snp_pairs_to_keep.empty()){
    vector<string> lines;
    readFile(file_gene_snp_pairs_to_keep, lines);
    
    size_t nb_gene_snp_pairs_to_keep = 0;
    vector<string> tokens;
    for(vector<string>::const_iterator it = lines.begin();
	it != lines.end(); ++it){
      split(*it, "\t", tokens);
      if(tokens.size() < 2){
	cerr << "ERROR: filw with gene-snp pairs to keep should have"
	     << " two columns gene<tab>snp" << endl;
	exit(1);
      }
      if(gene_snp_pairs_to_keep.find(tokens[0]) == gene_snp_pairs_to_keep.end())
	gene_snp_pairs_to_keep.insert(make_pair(tokens[0], vector<string>()));
      gene_snp_pairs_to_keep[tokens[0]].push_back(tokens[1]);
      ++nb_gene_snp_pairs_to_keep;
    }
    
    if(verbose > 0 && ! gene_snp_pairs_to_keep.empty())
      cout << "nb of gene-snp pairs to keep: " << nb_gene_snp_pairs_to_keep
	   << " (" << gene_snp_pairs_to_keep.size() << " genes)" << endl;
  }
}

void averageBFs(const string & file_bf_out,
		const vector<string> & lines,
		const vector<double> & grid_weights,
		const vector<size_t> & grid_idx_to_keep,
		vector<BayesFactor> & bfs)
{
  vector<string> tokens;
  
  split(lines[0], " \t,", tokens);
  if(tokens[0].compare("gene") != 0 && tokens[1].compare("snp") != 0
     && tokens[2].compare("config") != 0){
    cerr << "ERROR: wrong header line in file " << file_bf_out << endl;
    exit(1);
  }
  
  for(size_t i = 1; i < lines.size(); ++i){
    split(lines[i], " \t,", tokens);
    BayesFactor bf(tokens[0], tokens[1], tokens[2]);
    bf.avg_raw_bfs(vector<string> (tokens.begin()+3,
				   tokens.begin()+3+grid_weights.size()),
		   grid_weights, grid_idx_to_keep);
    bfs.push_back(bf);
  }
}

void saveAvgBFs(const vector<BayesFactor> & bfs,
		const string & file_hm,
		gzFile & stream_hm,
		stringstream & txt,
		const size_t & nb_lines_hm)
{
  for(size_t i = 0; i < bfs.size(); ++i){
    txt.str("");
    txt << bfs[i].gene << "\t" << bfs[i].snp
	<< "\t" << bfs[i].config << "\t" << bfs[i].log10_val << "\n";
    gzwriteLine(stream_hm, txt.str(), file_hm, nb_lines_hm);
  }
}

void averageRawBFsOverGrid(const vector<string> & files_bf_out,
			   const vector<double> & grid_weights,
			   const vector<size_t> & grid_idx_to_keep,
			   const string & file_hm,
			   const int & verbose)
{
  if(verbose > 0)
    cout << "average raw BFs over grid only and save them ..." << endl;
  
  gzFile stream_hm;
  openFile(file_hm, stream_hm, "wb");
  stringstream txt;
  txt.precision(6);
  txt.setf(ios::scientific);
  size_t nb_lines_hm = 0;
  
  txt << "gene\tsnp\tconfig\tl10abf.grid.avg" << endl;
  nb_lines_hm = 1;
  gzwriteLine(stream_hm, txt.str(), file_hm, nb_lines_hm);
  
  vector<string> lines;
  vector<BayesFactor> bfs;
  for(size_t f = 0; f < files_bf_out.size(); ++f){
    if(verbose == 1)
      progressBar("", f+1, files_bf_out.size());
    
    if(verbose > 1)
      cout << "read file " << files_bf_out[f] << endl << flush;
    clock_t startTime = clock();
    readFile(files_bf_out[f], lines);
    if(verbose > 1)
      cout << "nb of lines: " << lines.size() << " (loaded in "
	   << fixed << setprecision(2) << getElapsedTime(startTime) << " sec)"
	   << endl << flush;
    
    averageBFs(files_bf_out[f], lines, grid_weights, grid_idx_to_keep, bfs);
    
    saveAvgBFs(bfs, file_hm, stream_hm, txt, nb_lines_hm);
    
    lines.clear();
    bfs.clear();
  }
  if(verbose == 1)
    cerr << endl << flush;
  
  closeFile(file_hm, stream_hm);
}


void parseLines(const string & file_bf_out,
		const vector<string> & lines,
		const vector<string> & genes_to_keep,
		const vector<string> & snps_to_keep,
		const map<string,vector<string> > & gene_snp_pairs_to_keep,
		const vector<double> & grid_weights,
		const vector<string> & config_names,
		vector<Gene> & genes)
{
  vector<string> tokens;
  
  split(lines[0], " \t,", tokens);
  if(tokens[0].compare("gene") != 0 && tokens[1].compare("snp") != 0
     && tokens[2].compare("config") != 0){
    cerr << "ERROR: wrong header line in file " << file_bf_out << endl;
    exit(1);
  }
  
  string curr_gene, curr_snp;
  Gene gene;
  vector<vector<double> > log10_bfs; // dim1 is config, dim2 is grid
  map<string,vector<string> >::const_iterator it_g;
  vector<string>::const_iterator it_s;
  
  for(size_t i = 1; i < lines.size(); ++i){
    split(lines[i], " \t,", tokens);
    
    // skip if "meta-analysis BF"
    // and, for "types" model, only load subgroup-specific configs
    if(find(config_names.begin(), config_names.end(), tokens[2])
       == config_names.end())
      continue;
    
    // skip genes, snps or gene-snp pairs if necessary
    if(! genes_to_keep.empty() 
       && find(genes_to_keep.begin(), genes_to_keep.end(), tokens[0])
       == genes_to_keep.end())
      continue;
    if(! snps_to_keep.empty() 
       && find(snps_to_keep.begin(), snps_to_keep.end(), tokens[1])
       == snps_to_keep.end())
      continue;
    if(! gene_snp_pairs_to_keep.empty()){
      it_g = gene_snp_pairs_to_keep.find(tokens[0]);
      if(it_g == gene_snp_pairs_to_keep.end())
	continue;
      else{
	it_s = find(it_g->second.begin(), it_g->second.end(), tokens[1]);
	if(it_s == it_g->second.end())
	  continue;
      }
    }
    
    if(tokens[0].compare(curr_gene) != 0){ // if first or new gene
      if(! gene.name_.empty()){ // if new gene
    	gene.snps_.push_back(Snp(curr_snp, log10_bfs));
    	genes.push_back(gene);
      }
      curr_gene = tokens[0];
      gene = Gene(curr_gene);
      log10_bfs.clear();
    }
    
    if(tokens[1].compare(curr_snp) != 0){ // if first or new snp
      if(! log10_bfs.empty()) // if new snp
    	gene.snps_.push_back(Snp(curr_snp, log10_bfs));
      curr_snp = tokens[1];
      log10_bfs.clear();
    }
    
    size_t config_idx = log10_bfs.size();
    log10_bfs.push_back(vector<double> (grid_weights.size(), NaN));
    for(size_t j = 0; j < grid_weights.size(); ++j)
      log10_bfs[config_idx][j] = atof(tokens[3+j].c_str());
  } // end of "for each line"
  
  if(! gene.name_.empty() && ! curr_snp.empty()){
    gene.snps_.push_back(Snp(curr_snp, log10_bfs));
    genes.push_back(gene);
  }
}

void averageBFs(const vector<double> & grid_weights,
		const vector<size_t> & grid_idx_to_keep,
		const string & model,
		const vector<string> & config_names,
		const vector<double> & config_weights,
		const vector<vector<string> > & config2subgroups,
		const vector<double> & type_weights,
		const vector<vector<double> > & subgroup_weights,
		const vector<string> & quantities_to_save,
		const double & pi0,
		const vector<string> & post_probas_to_save,
		const int & save_best_snps,
		const bool & save_best_dim,
		const bool & save_all_dims,
		const int & nb_threads,
		vector<Gene> & genes)
{
#pragma omp parallel for num_threads(nb_threads)
  for(int g = 0; g < (int) genes.size(); ++g){
#ifdef DEBUG
    fprintf(stderr, "average BFs of %s\n", genes[g].name_.c_str());
#endif
    
    // average BFs
    if(model == "configs")
      genes[g].avg_raw_bfs(grid_weights, grid_idx_to_keep, config_weights);
    else if(model == "types")
      genes[g].avg_raw_bfs(grid_weights, grid_idx_to_keep, type_weights, subgroup_weights);
    
    // calculate posteriors
    if(find(quantities_to_save.begin(), quantities_to_save.end(), "post")
       != quantities_to_save.end()){
      
      if(find(post_probas_to_save.begin(), post_probas_to_save.end(), "a")
	 != post_probas_to_save.end())
	genes[g].calc_posterior(pi0);
      
      if(find(post_probas_to_save.begin(), post_probas_to_save.end(), "b")
	 != post_probas_to_save.end())
	genes[g].calc_cond_snp_posteriors_the_eQTL();
      
      if(find(post_probas_to_save.begin(), post_probas_to_save.end(), "c")
	 != post_probas_to_save.end())
	genes[g].calc_cond_snp_posteriors_an_eQTL();
      
      if(save_best_dim || save_all_dims){
	if(model == "configs")
	  genes[g].calc_cond_snp_posteriors_config(config_weights);
	else if(model == "types")
	  genes[g].calc_cond_snp_posteriors_type(type_weights);
      }
      
      if(find(post_probas_to_save.begin(), post_probas_to_save.end(), "d")
	 != post_probas_to_save.end()){
	if(model == "configs"){
	  if(! (save_best_dim || save_all_dims))
	    genes[g].calc_cond_snp_posteriors_config(config_weights);
	  genes[g].calc_cond_snp_posteriors_subgroup(config_names,
						     config2subgroups);
	}
	else if(model == "types")
	  genes[g].calc_cond_snp_posteriors_subgroup(type_weights,
						     subgroup_weights,
						     grid_weights,
						     grid_idx_to_keep);
      }
    } // end of "if save posteriors"
    
    genes[g].identify_best_snps(save_best_snps);
    
  } // end of "for each gene"
}

void saveAvgBFs(const size_t & nb_subgroups,
		const vector<Gene> & genes,
		const vector<string> & config_names,
		const vector<string> & quantities_to_save,
		const vector<string> & post_probas_to_save,
		const string & model,
		const bool & save_best_dim,
		const bool & save_all_dims,
		const string & file_hm,
		gzFile & stream_hm,
		stringstream & txt,
		size_t & nb_lines_hm)
{
  for(size_t g = 0; g < genes.size(); ++g){
    for(size_t p = 0; p <= genes[g].idx_last_best_snp_; ++p){
      txt.str("");
      txt << setprecision(4) << scientific
	  << genes[g].name_
	  << "\t" << genes[g].snps_[p].name_;
      
      // save averaged BFs
      if(find(quantities_to_save.begin(), quantities_to_save.end(), "bf")
      	 != quantities_to_save.end()){
      	txt << "\t" << genes[g].log10_bf_
      	    << "\t" << genes[g].snps_[p].log10_bf_;
      	if(save_all_dims)
      	  for(size_t k = 0; k < genes[g].snps_[p].dim_log10_bfs_.size(); ++k)
      	    txt << "\t" << genes[g].snps_[p].dim_log10_bfs_[k];
      }
      
      // save posteriors
      if(find(quantities_to_save.begin(), quantities_to_save.end(), "post")
      	 != quantities_to_save.end()){
      	if(find(post_probas_to_save.begin(), post_probas_to_save.end(), "a")
      	   != post_probas_to_save.end())
      	  txt << "\t" << genes[g].post_;
      	if(find(post_probas_to_save.begin(), post_probas_to_save.end(), "b")
      	   != post_probas_to_save.end())
      	  txt << "\t" << genes[g].snps_[p].post_the_eQTL_;
      	if(find(post_probas_to_save.begin(), post_probas_to_save.end(), "c")
      	   != post_probas_to_save.end())
      	  txt << "\t" << genes[g].snps_[p].post_an_eQTL_;
      	if(find(post_probas_to_save.begin(), post_probas_to_save.end(), "d")
      	   != post_probas_to_save.end()){
      	  for(size_t s = 0; s < nb_subgroups; ++s)
      	    txt << "\t" << genes[g].snps_[p].post_subgroups_[s];
      	}
      	if(save_all_dims)
      	  for(size_t k = 0; k < genes[g].snps_[p].post_dims_.size(); ++k)
      	    txt << "\t" << genes[g].snps_[p].post_dims_[k];
      }
      
      if(save_best_dim){
	if(model == "configs")
	  txt << "\t" << config_names[genes[g].snps_[p].best_dim_idx_];
	else if(model == "types")
	  txt << "\t" << genes[g].snps_[p].best_dim_idx_ + 1;
	txt << "\t" << genes[g].snps_[p].post_dims_[genes[g].snps_[p].best_dim_idx_];
      }
      
      txt << "\n";
      ++nb_lines_hm;
      gzwriteLine(stream_hm, txt.str(), file_hm, nb_lines_hm);
      
    } // end of "for each SNP"
  } // end of "for each gene"
}

void averageRawBFsOverGridAndOthers(
  const vector<string> & files_bf_out,
  const vector<double> & grid_weights,
  const vector<size_t> & grid_idx_to_keep,
  const string & model,
  const size_t & nb_subgroups,
  const size_t & dim,
  const vector<string> & config_names,
  const vector<double> & config_weights,
  const vector<vector<string> > & config2subgroups,
  const vector<double> & type_weights,
  const vector<vector<double> > & subgroup_weights,
  const vector<string> & quantities_to_save,
  const double & pi0,
  const vector<string> & post_probas_to_save,
  const vector<string> & genes_to_keep,
  const vector<string> & snps_to_keep,
  const map<string,vector<string> > & gene_snp_pairs_to_keep,
  const int & save_best_snps,
  const bool & save_best_dim,
  const bool & save_all_dims,
  const string & file_hm,
  const int & nb_threads,
  const int & verbose)
{
  if(verbose > 0)
    cout << "average raw BFs over hyperparameters and save them ..." << endl;
  
  gzFile stream_hm;
  openFile(file_hm, stream_hm, "wb");
  stringstream txt;
  size_t nb_lines_hm = 0;
  
  // write pi0 if posteriors are saved
  if(find(quantities_to_save.begin(), quantities_to_save.end(), "post")
     != quantities_to_save.end()){
    txt << "#pi0=" << setprecision(4) << scientific << pi0 << endl;
    ++nb_lines_hm;
    gzwriteLine(stream_hm, txt.str(), file_hm, nb_lines_hm);
  }
  
  // write header line
  txt.str("");
  txt << "gene\tsnp";
  if(find(quantities_to_save.begin(), quantities_to_save.end(), "bf")
     != quantities_to_save.end()){
    txt << "\tgene.log10.bf\tsnp.log10.bf";
    if(save_all_dims){
      if(model == "configs")
	for(size_t j = 0; j < dim; ++j)
	  txt << "\tlog10.bf." << config_names[j];
      else if(model == "types")
	for(size_t k = 0; k < dim; ++k)
	  txt << "\tlog10.bf.type." << k+1;
    }
  } // end of "if save bf"
  if(find(quantities_to_save.begin(), quantities_to_save.end(), "post")
     != quantities_to_save.end()){
    if(find(post_probas_to_save.begin(), post_probas_to_save.end(), "a")
       != post_probas_to_save.end())
      txt << "\tgene.post";
    if(find(post_probas_to_save.begin(), post_probas_to_save.end(), "b")
       != post_probas_to_save.end())
      txt << "\tsnp.post.the";
    if(find(post_probas_to_save.begin(), post_probas_to_save.end(), "c")
       != post_probas_to_save.end())
      txt << "\tsnp.post.an";
    if(find(post_probas_to_save.begin(), post_probas_to_save.end(), "d")
       != post_probas_to_save.end()){
      for(size_t s = 0; s < nb_subgroups; ++s)
	txt << "\tsnp.post.s" << s+1;
    }
    if(save_all_dims){
      if(model == "configs")
	for(size_t j = 0; j < dim; ++j)
	  txt << "\tpost.config." << config_names[j];
      else if(model == "types")
	for(size_t k = 0; k < dim; ++k)
	  txt << "\tpost.type." << k+1;
    }
  } // end of "if save post"
  if(save_best_dim){
    if(model == "configs")
      txt << "\tbest.config"
	  << "\tpost.best.config";
    else if(model == "types")
      txt << "\tbest.type"
	  << "\tpost.best.type";
  }
  txt << endl;
  ++nb_lines_hm;
  gzwriteLine(stream_hm, txt.str(), file_hm, nb_lines_hm);
  
  // read each input file, process lines and save results
  vector<string> lines;
  vector<Gene> genes;
  for(size_t f = 0; f < files_bf_out.size(); ++f){
    if(verbose == 1)
      progressBar("", f+1, files_bf_out.size());
    
    if(verbose > 1)
      cout << "read file " << files_bf_out[f] << endl << flush;
    clock_t startTime = clock();
    readFile(files_bf_out[f], lines); // quicker to read whole file at once
    if(verbose > 1)
      cout << "nb of lines: " << lines.size() << " (loaded in "
	   << fixed << setprecision(2) << getElapsedTime(startTime) << " sec)"
	   << endl << flush;
    
    startTime = clock();
    parseLines(files_bf_out[f], lines, genes_to_keep, snps_to_keep,
	       gene_snp_pairs_to_keep, grid_weights, config_names, genes);
    if(verbose > 1){
      size_t nb_snps = 0;
      for(size_t g = 0; g < genes.size(); ++g)
	nb_snps += genes[g].snps_.size();
      cout << "nb of genes: " << genes.size() << " (" << nb_snps << " gene-SNP pairs"
	   << ", parsed in " << fixed << setprecision(2)
	   << getElapsedTime(startTime) << " sec)" << endl << flush;
    }
    
    averageBFs(grid_weights, grid_idx_to_keep, model,
	       config_names, config_weights,
	       config2subgroups, type_weights, subgroup_weights,
	       quantities_to_save,
	       pi0, post_probas_to_save, save_best_snps,
	       save_best_dim, save_all_dims, nb_threads, genes);
    
    saveAvgBFs(nb_subgroups, genes, config_names, quantities_to_save,
  	       post_probas_to_save, model, save_best_dim, save_all_dims,
	       file_hm, stream_hm, txt, nb_lines_hm);
    
    lines.clear();
    genes.clear();
  } // enf of "for each raw BF file"
  if(verbose == 1)
    cerr << endl << flush;
  
  closeFile(file_hm, stream_hm);
}

void run(const string & file_pattern,
	 const string & file_grid_weights,
	 vector<size_t> & grid_idx_to_keep,
	 const string & model,
	 const size_t & nb_subgroups,
	 const size_t & dim,
	 const string & file_config_weights,
	 const string & file_type_subgroup_weights,
	 const vector<string> & quantities_to_save,
	 const double & pi0,
	 const vector<string> & post_probas_to_save,
	 const string & file_genes_to_keep,
	 const string & file_snps_to_keep,
	 const string & file_gene_snp_pairs_to_keep,
	 const int & save_best_snps,
	 const bool & save_best_dim,
	 const bool & save_all_dims,
	 const string & file_hm,
	 const int & nb_threads,
	 const int & verbose)
{
  if(verbose > 0 && !isNan(pi0))
    cout << "pi0: " << setprecision(4) << scientific << pi0 << endl;
  
  vector<double> grid_weights;
  loadFileGridWeights(file_grid_weights, verbose,
		      grid_weights, grid_idx_to_keep);
  
  vector<string> config_names;
  vector<double> config_weights;
  vector<vector<string> > config2subgroups;
  loadFileConfigWeights(file_config_weights, nb_subgroups, dim, verbose,
			config_names, config_weights, config2subgroups);
  
  vector<double> type_weights;
  vector<vector<double> > subgroup_weights;
  loadFileTypeSubgroupWeights(file_type_subgroup_weights, nb_subgroups, dim,
			      verbose, type_weights, subgroup_weights);
  
  vector<string> genes_to_keep;
  if(! file_genes_to_keep.empty()){
    readFile(file_genes_to_keep, genes_to_keep);
    if(! genes_to_keep.empty())
      cout << "genes to keep: " << genes_to_keep.size() << endl;
  }
  
  vector<string> snps_to_keep;
  if(! file_snps_to_keep.empty()){
    readFile(file_snps_to_keep, snps_to_keep);
    if(! snps_to_keep.empty())
      cout << "snps to keep: " << snps_to_keep.size() << endl;
  }
  
  map<string,vector<string> > gene_snp_pairs_to_keep;
  loadFileGeneSnpPairsToKeep(file_gene_snp_pairs_to_keep, verbose,
			     gene_snp_pairs_to_keep);
  
  vector<string> files_bf_out = glob(file_pattern);
  if(verbose > 0)
    cerr << "nb of files with raw BFs: " << files_bf_out.size()
	 << endl << flush;
  
  // reduce input for 'eqtlbma_hm'
  if(file_config_weights.empty() && file_type_subgroup_weights.empty())
    averageRawBFsOverGrid(files_bf_out, grid_weights, grid_idx_to_keep,
			  file_hm, verbose);
  else
    averageRawBFsOverGridAndOthers(files_bf_out, grid_weights,
				   grid_idx_to_keep, model,
				   nb_subgroups, dim,
				   config_names, config_weights,
				   config2subgroups, type_weights,
				   subgroup_weights,
				   quantities_to_save, pi0,
				   post_probas_to_save, genes_to_keep,
				   snps_to_keep, gene_snp_pairs_to_keep,
				   save_best_snps, save_best_dim,
				   save_all_dims, file_hm, nb_threads,
				   verbose);
}

int main(int argc, char ** argv)
{
#ifdef DEBUG
  fprintf(stderr, "DEBUG\n");
#endif
  int verbose = 1, nb_threads = 1, save_best_snps = 0;
  size_t nb_subgroups = string::npos, dim = string::npos;
  string file_pattern, file_grid_weights, model = "configs",
    file_config_weights, file_type_subgroup_weights, file_hm,
    file_genes_to_keep, file_snps_to_keep, file_gene_snp_pairs_to_keep;
  double pi0 = NaN;
  vector<size_t> grid_idx_to_keep;
  vector<string> quantities_to_save, post_probas_to_save;
  bool save_all_dims = false, save_best_dim = false;
  
  parseCmdLine(argc, argv, verbose, file_pattern, file_grid_weights,
	       grid_idx_to_keep, model, nb_subgroups, dim,
	       file_config_weights, file_type_subgroup_weights,
	       quantities_to_save, pi0, post_probas_to_save,
	       file_genes_to_keep, file_snps_to_keep, file_gene_snp_pairs_to_keep,
	       save_best_snps, save_best_dim, save_all_dims, file_hm, nb_threads);
  
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
  
  run(file_pattern, file_grid_weights, grid_idx_to_keep, model,
      nb_subgroups, dim, file_config_weights, file_type_subgroup_weights,
      quantities_to_save, pi0, post_probas_to_save, file_genes_to_keep,
      file_snps_to_keep, file_gene_snp_pairs_to_keep, save_best_snps,
      save_best_dim, save_all_dims, file_hm, nb_threads, verbose);
  
  if(verbose > 0){
    time(&time_end);
    cout << "END " << basename(argv[0])
	 << " " << getDateTime(time_end) << endl
	 << "elapsed -> " << getElapsedTime(time_start, time_end) << endl
	 << "max.mem -> " << getMaxMemUsedByProcess2Str() << endl;
  }
  
  return EXIT_SUCCESS;
}
