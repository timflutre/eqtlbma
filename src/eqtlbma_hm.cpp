/** \file eqtlbma_hm.cpp
 *
 *  `eqtlbma_hm' implements the EM algorithm to fit the hierarchical model from eQtlBma.
 *  Copyright (C) 2012-2015 Xiaoquan Wen, Timothée Flutre
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
 * g++ -Wall -Wextra -g eqtlbma_hm.cpp hm_methods.cpp ../eqtlbma/utils/utils_io.cpp ../eqtlbma/utils/utils_math.cpp -I../eqtlbma -lgsl -lgslcblas -lz -fopenmp -o eqtlbma_hm
 */

#include <cmath>
#include <cstring>
#include <getopt.h>
#include <libgen.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "utils/utils_io.hpp"
#include "utils/utils_math.hpp"
using namespace utils;

#include <omp.h>

#include "hm_classes.hpp"

#ifndef VERSION
#define VERSION "0.0"
#endif

//-----------------------------------------------------------------------------

class Controller
{
  void load_data_one_file(const string & file,
                          const vector<string> & configs_tokeep,
                          const bool & keep_gen_abfs,
                          vector<string> & lines);
  
public:
  vector<gene_eQTL> genes_;
  size_t nb_subgroups_;
  size_t grid_size_;
  string model_; // configs or types
  size_t dim_; // nb of configs or nb of types
  vector<string> config_names_;
  double thresh_; // on log-lik, to stop EM
  size_t max_nb_iters_; // to stop EM
  bool squarem_; // to speed-up EM
  map<string,bool> param2fixed_;
  int nb_threads_; 
  int verbose_;
  
  // data structures for the likelihood and parameters to estimate
  double log10_obs_lik_;
  double log10_obs_lik_0;
  double log10_obs_lik_1;
  double log10_obs_lik_2;
    double new_log10_obs_lik_;
  double pi0_;
  double pi0_0;
  double pi0_1;
  double pi0_2;
    double new_pi0_;
  vector<double> grid_wts_; // lambda_l
  vector<double> grid_wts_0;
  vector<double> grid_wts_1;
  vector<double> grid_wts_2;
    vector<double> new_grid_wts_;
  vector<double> config_prior_; // eta_j
  vector<double> config_prior_0;
  vector<double> config_prior_1;
  vector<double> config_prior_2;
    vector<double> new_config_prior_;
  vector<double> type_prior_; // pi_k
  vector<double> type_prior_0;
  vector<double> type_prior_1;
  vector<double> type_prior_2;
    vector<double> new_type_prior_;
  vector<vector<double> > subgroup_prior_; // q_ks
  vector<vector<double> > subgroup_prior_0;
  vector<vector<double> > subgroup_prior_1;
  vector<vector<double> > subgroup_prior_2;
    vector<vector<double> > new_subgroup_prior_;
  
  // data structures for the confidence intervals
  double left_pi0_, right_pi0_;
  vector<double> left_grids_, right_grids_;
  vector<double> left_configs_, right_configs_;
  vector<double> left_types_, right_types_;
  vector<vector<double> > left_subgroups_, right_subgroups_;
  
  // utilitary data structures for EM updates
  vector<double> gene_wts_ones_; // vector of 1's per gene
  vector<double> grid_wts_ones_; // idem per grid
  vector<double> config_wts_ones_; // idem per config
  vector<double> type_wts_ones_; // idem per type
  vector<vector<double> > grid_genes_; // log10 values per grid per gene
  vector<vector<double> > config_genes_; // idem per config per gene
  vector<vector<double> > type_genes_; // idem per type per gene
  vector<vector<vector<double> > > subgroup_num_genes_; // idem per type and subgroup per gene
  vector<vector<double> > subgroup_denom_genes_; // idem per type per gene


  // data structure for control variables in squarem
  int squaremK_=1;//current implementation
  int method_=3;//1,2,3 indicates the types of step length to be used in squarem1,squarem2, 4,5 for "rre" and "mpe" in cyclem1 and cyclem2,  standing for reduced-rank ("rre") or minimal-polynomial ("mpe") extrapolation.
  double mstep_=4;
  //int maxiter=1500;max_nb_iters_
  bool square=true;
  //bool trace=true;verbose_
  double stepmin0_=1;
  double stepmax0_=1;
  double kr_=1;
  double objfninc_=1;//0 to enforce monotonicity, Inf for non-monotonic scheme, 1 for monotonicity far from solution and allows for non-monotonicity closer to solution
  
  Controller();
  Controller(const size_t & nb_subgroups, const string & model,
             const size_t & grid_size,
             const size_t & dim, const double & thresh,
             const size_t & max_nb_iters,
             const bool & squarem,
             const double & fixed_pi0, const int & nb_threads,
             const int & verbose);
  
  void load_data(const string & file_pattern,
                 const vector<string> & configs_tokeep,
                 const bool & keep_gen_abfs);
  
  void init_params(const size_t & seed);
  void init_params(const string & init_file);
  
  void randomize_parameter_sp(const size_t & seed);
  
  // params as args allows to calc lik with mix of old and new params
  double compute_log10_obs_lik(const double & pi0,
                               const vector<double> & grid_wts,
                               const vector<double> & config_prior,
                               const vector<double> & type_prior,
                               const vector<vector<double> > & subgroup_prior,
                               const bool & keep);
  
  void em_update_pi0();
  void em_update_config();
  void em_update_type();
  void em_update_subgroup();
  void em_update_grid();
  void em_update_snp_prior();

  void run_EM_fixedpoint();
  void update_params(int squaremstep=0);
  
  void show_state_EM(const size_t & iter);
  
  void run_EM_classic();
  void run_EM_square();
  void run_EM();
  
  void compute_log10_ICL();
  
  void compute_posterior();
  
  void estimate_profile_ci_pi0(const double & tick,
                               const double & max_l10_obs_lik);
  void estimate_profile_ci_configs(const double & tick,
                                   const double & max_l10_obs_lik);
  void estimate_profile_ci_types(const double & tick,
                                 const double & max_l10_obs_lik);
  void estimate_profile_ci_subgroups(const double & tick,
                                     const double & max_l10_obs_lik);
  void estimate_profile_ci_grids(const double & tick,
                                 const double & max_l10_obs_lik);
  void estimate_profile_ci();
  
  void save_result(const string & out_file, const bool & skip_bf);
};


Controller::Controller()
{
}

Controller::Controller(const size_t & nb_subgroups,
                       const string & model,
                       const size_t & gsize,
                       const size_t & dim,
                       const double & thresh,
                       const size_t & max_nb_iters,
                       const bool & squarem,
                       const double & fixed_pi0,
                       const int & nb_threads,
                       const int & verbose)
{
  grid_size_ = gsize;
  model_ = model;
  nb_subgroups_ = nb_subgroups;
  
  left_pi0_ = NaN;
  right_pi0_ = NaN;
  
  grid_wts_.assign(grid_size_, NaN);
  new_grid_wts_.assign(grid_size_, NaN);
  grid_wts_ones_.assign(grid_size_, 1.0);
  left_grids_.assign(grid_size_, NaN);
  right_grids_.assign(grid_size_, NaN);
  
  dim_ = dim;
  if(model_ == "configs"){
    config_prior_.assign(dim_, NaN);
    new_config_prior_.assign(dim_, NaN);
    config_wts_ones_.assign(dim_, 1.0);
    left_configs_.assign(dim_, NaN);
    right_configs_.assign(dim_, NaN);
  }
  else if(model_ == "types"){
    type_prior_.assign(dim_, NaN);
    new_type_prior_.assign(dim_, NaN);
    type_wts_ones_.assign(dim_, 1.0);
    left_types_.assign(dim_, NaN);
    right_types_.assign(dim_, NaN);
    subgroup_prior_.assign(dim_, vector<double>(nb_subgroups_, NaN));
    new_subgroup_prior_.assign(dim_, vector<double>(nb_subgroups_, NaN));
    left_subgroups_.assign(dim_, vector<double>(nb_subgroups_, NaN));
    right_subgroups_.assign(dim_, vector<double>(nb_subgroups_, NaN));
  }
  
  param2fixed_.insert(make_pair("pi0", false));
  if(! isNan(fixed_pi0)){
    pi0_ = fixed_pi0;
    param2fixed_["pi0"] = true;
  }
  param2fixed_.insert(make_pair("grid-points", false));
  param2fixed_.insert(make_pair("configs", false));
  param2fixed_.insert(make_pair("types", false));
  param2fixed_.insert(make_pair("subgroups-per-type", false));
  
  thresh_ = thresh;
  max_nb_iters_ = max_nb_iters;
  squarem_ = squarem;
  nb_threads_ = nb_threads;
  verbose_ = verbose;
}

void Controller::load_data_one_file(
  const string & file,
  const vector<string> & configs_tokeep,
  const bool & keep_gen_abfs,
  vector<string> & lines)
{
  // load the whole file at once (supposedly quicker)
  readFile(file, lines);
  if(verbose_ > 1)
    cout << lines.size() << " lines" << endl;
  
  // check the header
  vector<string> tokens;
  split(lines[0], " \t,", tokens);
  if(tokens[0].compare("gene") != 0 
     || tokens[1].compare("snp") != 0
     || tokens[2].compare("config") != 0){
    cerr << "ERROR: file " << file << " has wrong header line" << endl;
    exit(EXIT_FAILURE);
  }
  
  // parse the lines
  char * pch;
  string gene_id, snp_id, config, curr_gene, curr_snp;
  gene_eQTL geq;
  vector<vector<double> > raw_log10_bfs; // dim1 is config, dim2 is grid
  for(size_t i = 1; i < lines.size(); ++i){
    
    // skip line based on config if needed
    pch = strtok((char *) lines[i].c_str(), " \t,");
    gene_id = string(pch);
    pch = strtok(NULL, " \t,");
    snp_id = string(pch);
    pch = strtok(NULL, " \t,");
    config = string(pch);
    if((! keep_gen_abfs & config.find("gen") != string::npos)
       || (! configs_tokeep.empty()
           && find(configs_tokeep.begin(), configs_tokeep.end(), config)
           == configs_tokeep.end()))
      continue;
    if(keep_gen_abfs & (config == "gen-fix" || config == "gen-maxh"))
      continue;
    
    // record config names once (for output)
    if(config_names_.size() < config_prior_.size())
      config_names_.push_back(config);
    
    if(gene_id.compare(curr_gene) != 0){ // if new gene
      if(! geq.name_.empty()){
        geq.snps_.push_back(snp_eQTL(curr_snp, raw_log10_bfs,
                                     nb_subgroups_,
                                     (model_ == "types" ? dim_ : 0)));
        genes_.push_back(geq);
      }
      curr_gene = gene_id;
      geq = gene_eQTL(curr_gene, nb_subgroups_, dim_, grid_size_, &pi0_);
      raw_log10_bfs.clear();
    }
    
    if(snp_id.compare(curr_snp) != 0){ // if same gene but new snp
      if(! raw_log10_bfs.empty())
        geq.snps_.push_back(snp_eQTL(curr_snp, raw_log10_bfs,
                                     nb_subgroups_,
                                     (model_ == "types" ? dim_ : 0)));
      curr_snp = snp_id;
      raw_log10_bfs.clear();
    }
    
    // if same gene and snp but new config
    size_t config_idx = raw_log10_bfs.size();
    raw_log10_bfs.push_back(vector<double> (grid_wts_.size(), NaN));
    size_t j = 0;
    pch = strtok(NULL, " \t,");
    while(j < grid_wts_.size() && pch != NULL){
      raw_log10_bfs[config_idx][j] = atof(pch);
      ++j;
      pch = strtok(NULL, " \t,");
    }
  }
  
  geq.snps_.push_back(snp_eQTL(curr_snp, raw_log10_bfs,
                               nb_subgroups_,
                               (model_ == "types" ? dim_ : 0))); // last snp
  genes_.push_back(geq); // last gene
}

void Controller::load_data(
  const string & file_pattern,
  const vector<string> & configs_tokeep,
  const bool & keep_gen_abfs)
{
  if(verbose_ > 0)
    fprintf(stderr, "load data ...\n");
  clock_t startTime = clock();
  
  if(model_ == "configs" && ! configs_tokeep.empty()){
    fprintf(stderr, "configurations to keep: %s", configs_tokeep[0].c_str());
    for(size_t i = 1; i < configs_tokeep.size(); ++i)
      fprintf(stderr, " %s", configs_tokeep[i].c_str());
    fprintf(stderr, "\n");
  }
  
  vector<string> files = glob(file_pattern);
  if(files.size() == 0){
    cerr << "ERROR: no input file was found from pattern " << file_pattern << endl;
    exit(EXIT_FAILURE);
  }
  else if(verbose_ > 0)
    cout << "nb of input files: " << files.size() << endl << flush;
  
  vector<string> lines;
  for(size_t i = 0; i < files.size(); ++i){
    if(verbose_ == 1)
      progressBar("", i+1, files.size());
    else if(verbose_ > 1)
      cout << "file " << (i+1) << " " << files[i] << endl;
    load_data_one_file(files[i], configs_tokeep, keep_gen_abfs, lines);
    lines.clear();
    // vector<string>().swap(lines);
  }
  if(verbose_ == 1)
    cout << endl << flush;
  
  // initialization of useful data structures for the EM
  gene_wts_ones_.assign(genes_.size(), 1.0);
  grid_genes_ =
    vector<vector<double> >(grid_size_,
                            vector<double>(genes_.size(), NaN));
  if(model_ == "configs"){
    config_genes_ = 
      vector<vector<double> >(dim_,
                              vector<double>(genes_.size(), NaN));
  }
  else if(model_ == "types"){
    type_genes_ = 
      vector<vector<double> >(dim_,
                              vector<double>(genes_.size(), NaN));
    subgroup_num_genes_ =
      vector<vector<vector<double> > >(dim_,
                                       vector<vector<double> >(nb_subgroups_,
                                                               vector<double>(genes_.size(), NaN)));
    subgroup_denom_genes_ =
      vector<vector<double> >(dim_,
                              vector<double>(genes_.size(), NaN));;
  }
  
  if(verbose_ > 0){
    size_t nb_pairs = 0;
#pragma omp parallel for num_threads(nb_threads_) reduction(+:nb_pairs)
    for(int g = 0; g < (int)genes_.size(); ++g)
      nb_pairs += genes_[g].snps_.size();
    fprintf(stderr, "finish loading %zu genes and %zu gene-snp pairs (%f sec, %s)\n",
            genes_.size(),
            nb_pairs,
            getElapsedTime(startTime),
            getMaxMemUsedByProcess2Str().c_str());
  }
}

void Controller::init_params(const string & init_file)
{
#pragma omp parallel for num_threads(nb_threads_)
  for(int i = 0; i < (int) genes_.size(); ++i)
    genes_[i].set_snp_prior();
  
  if(verbose_)
    cout << "loading initialization file " << init_file << " ..." << endl;
  vector<string> lines;
  readFile(init_file, lines);
  
  vector<string> tokens;
  size_t idx_grid = 0, idx_config = 0, idx_type = 0,
    idx_subgroup1 = 0, idx_subgroup2 = 0;
  for(size_t i = 0; i < lines.size(); ++i){
    if(lines[i].find("#") == 0) // skip commented line
      continue;
    split(lines[i], "\t, ", tokens);
    if(tokens[0] == "param" && tokens[1] == "value") // skip header line
      continue;
    if(tokens[0].find("pi0") != string::npos){
      pi0_ = atof(tokens[1].c_str());
      if(tokens.size() == 3 && (tokens[2] == "TRUE" || tokens[2] == "true"))
        param2fixed_["pi0"] = true;
    }
    else if(tokens[0].find("grid") != string::npos){
      grid_wts_[idx_grid] = atof(tokens[1].c_str());
      ++idx_grid;
      if(tokens.size() == 3 && (tokens[2] == "TRUE" || tokens[2] == "true"))
        param2fixed_["grid-points"] = true;
    }
    else if(tokens[0].find("config") != string::npos){
      config_prior_[idx_config] = atof(tokens[1].c_str());
      ++idx_config;
      if(tokens.size() == 3 && (tokens[2] == "TRUE" || tokens[2] == "true"))
        param2fixed_["configs"] = true;
    }
    else if(tokens[0].find("type") != string::npos){
      type_prior_[idx_type] = atof(tokens[1].c_str());
      ++idx_type;
      if(tokens.size() == 3 && (tokens[2] == "TRUE" || tokens[2] == "true"))
        param2fixed_["types"] = true;
    }
    else if(tokens[0].find("subgroup") != string::npos){
      if(idx_subgroup2 == nb_subgroups_){
        idx_subgroup2 = 0;
        ++idx_subgroup1;
      }
      subgroup_prior_[idx_subgroup1][idx_subgroup2] = atof(tokens[1].c_str());
      ++idx_subgroup2;
      if(tokens.size() == 3 && (tokens[2] == "TRUE" || tokens[2] == "true"))
        param2fixed_["subgroups-per-type"] = true;
    }
  }
  
  if(verbose_ > 0){
    bool update_all = true;
    for(map<string,bool>::const_iterator it = param2fixed_.begin();
        it != param2fixed_.end(); ++it){
      if(it->second){
        update_all = false;
        break;
      }
    }
    if(update_all)
      cout << "update all parameters" << endl;
    else{
      cout << "parameters to update:";
      for(map<string,bool>::const_iterator it = param2fixed_.begin();
          it != param2fixed_.end(); ++it){
        if(model_ == "configs" && (it->first == "types"
                                   || it->first == "subgroups-per-type"))
          continue;
        else if(model_ == "types" && it->first == "configs")
          continue;
        if(! it->second)
          cout << " " << it->first;
      }
      cout << endl;
    }
  }
}

void Controller::randomize_parameter_sp(const size_t & seed)
{
  gsl_rng_env_setup();
  const gsl_rng_type * T = gsl_rng_default;
  gsl_rng * r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed); 
  
  // double x1 = gsl_ran_exponential(r,1.0);
  // double x2 = gsl_ran_exponential(r,1.0);
  // pi0_ = .90 + .1*x1/(x1 + x2);
  // pi0_ = x1/(x1 + x2);
  if(! param2fixed_["pi0"])
    pi0_ = gsl_rng_uniform(r);
  
  double sum = 0.0;
  for(size_t l = 0; l < grid_size_; ++l){
    grid_wts_[l] = gsl_ran_exponential(r, 1.0);
    sum += grid_wts_[l];
  }
  for(size_t l = 0; l < grid_wts_.size(); ++l)
    grid_wts_[l] /= sum;
  
  if(model_ == "configs"){
    if(config_prior_.size() == 1){
      new_config_prior_[0] = config_prior_[0] = 1.0;
    }
    else{
      double sum = 0.0;
      for(size_t k = 0; k < dim_; ++k){
        config_prior_[k] = gsl_ran_exponential(r, 1.0);
        sum+= config_prior_[k];
      }
      for(size_t k = 0; k < dim_; ++k)
        config_prior_[k] /= sum;
    }
  }
  else if(model_ == "types"){
    double sum = 0.0;
    for(size_t k = 0; k < dim_; ++k){
      type_prior_[k] = gsl_ran_exponential(r, 1.0);
      sum+= type_prior_[k];
    }
    for(size_t k = 0; k < dim_; ++k){
      type_prior_[k] /= sum;
      for(size_t s = 0; s < nb_subgroups_; ++s)
        subgroup_prior_[k][s] = gsl_rng_uniform(r);
    } 
  }
  
  gsl_rng_free(r);
}

void Controller::init_params(const size_t & seed)
{
#pragma omp parallel for num_threads(nb_threads_)
  for(int i = 0; i < (int)genes_.size(); ++i)
    genes_[i].set_snp_prior();
  
  if(seed != string::npos){
    randomize_parameter_sp(seed);
  }
  else{
    if(! param2fixed_["pi0"])
      pi0_ = 0.5;
    
    for(size_t i = 0; i < grid_wts_.size(); ++i)
      grid_wts_[i] = 1.0 / (double) grid_wts_.size();
    
    if(model_ == "configs"){
      if(dim_ == 1){
        config_prior_[0] = 1.0;
        new_config_prior_[0] = 1.0;
      }
      else
        for(size_t i = 0; i < dim_; ++i)
          config_prior_[i] = 1.0 / (double) dim_;
    }
    else if(model_ == "types"){
      for(size_t k = 0; k < dim_; ++k){
        type_prior_[k] = 1.0 / (double) dim_;
        for(size_t s = 0; s < nb_subgroups_; ++s)
          subgroup_prior_[k][s] = 0.5;
      }
    }
  } // end of "non-random initialization"
}

// log10(obslik) = sum_g log10(p(Y_g | X_g, Theta))
double Controller::compute_log10_obs_lik(
  const double & pi0,
  const vector<double> & grid_wts,
  const vector<double> & config_prior,
  const vector<double> & type_prior,
  const vector<vector<double> > & subgroup_prior,
  const bool & keep)
{
  double l10_sum = 0.0;
  
  if(model_ == "configs"){
#pragma omp parallel for num_threads(nb_threads_) reduction(+:l10_sum)
    for(int g = 0; g < (int)genes_.size(); ++g)
      l10_sum += genes_[g].compute_log10_obs_lik(pi0, grid_wts,
                                                 config_prior, keep);
  }
  else if(model_ == "types"){
#pragma omp parallel for num_threads(nb_threads_) reduction(+:l10_sum)
    for(int g = 0; g < (int)genes_.size(); ++g)
      l10_sum += genes_[g].compute_log10_obs_lik(pi0, grid_wts, type_prior,
                                                 subgroup_prior, keep);
  }
  
  if(isNan(l10_sum)){
    cerr << "ERROR: log10(obslik) is NaN" << endl;
    exit(EXIT_FAILURE);
  }
  if(isInfinite(l10_sum)){
    cerr << "ERROR: log10(obslik) is +-Inf" << endl;
    exit(EXIT_FAILURE);
  }
  
  return l10_sum;
}

void Controller::em_update_snp_prior()
{
#pragma omp parallel for num_threads(nb_threads_)
  for (int i = 0; i < (int)genes_.size(); ++i)
    genes_[i].em_update_snp_prior();
}

void Controller::em_update_pi0()
{
  if(! param2fixed_["pi0"]){
    double sum = 0;
#pragma omp parallel for num_threads(nb_threads_) reduction(+:sum)
    for(int g = 0; g < (int)genes_.size(); ++g){
      // genes_[g].update_snp_prior();
      sum += genes_[g].em_update_pi0(pi0_);
    }
    new_pi0_ = sum / genes_.size();
  }
  else
    new_pi0_ = pi0_;
#ifdef DEBUG
  fprintf(stdout, "pi_0 old=%.4e new=%.4e\n", pi0_, new_pi0_);
  fflush(stdout);
#endif
}

void Controller::em_update_config()
{
  if(config_prior_.size() == 1)
    return;
  
  if(! param2fixed_["configs"]){
    // compute the contribution of each gene to each config weight
#pragma omp parallel for num_threads(nb_threads_)
    for(int g = 0; g < (int)genes_.size(); ++g){
      vector<double> new_config_prior_tmp(dim_, NaN);
      genes_[g].em_update_config(grid_wts_, new_config_prior_tmp);
      for(size_t k = 0; k < dim_; ++k)
        config_genes_[k][g] = new_config_prior_tmp[k];
    }
    
    // average the contribution of all genes per config
    for(size_t k = 0; k < dim_; ++k)
      new_config_prior_[k] = log10_weighted_sum(&(config_genes_[k][0]),
                                                &(gene_wts_ones_[0]),
                                                genes_.size())
        + log10(config_prior_[k]);
    
    // compute the normalization constant (Lagrange multiplier)
    double l10_denom = log10_weighted_sum(&(new_config_prior_[0]),
                                          &(config_wts_ones_[0]),
                                          dim_);
    
    // compute each new config weight
    for(size_t k = 0; k < dim_; ++k)
      new_config_prior_[k] = pow(10, (new_config_prior_[k] - l10_denom));
  } // end of "if configs not fixed"
  else
    new_config_prior_ = config_prior_;
  
#ifdef DEBUG
  for(size_t k = 0; k < dim_; ++k)
    fprintf(stdout, "eta_%zu old=%.4e new=%.4e\n", k+1, config_prior_[k], new_config_prior_[k]);
  fflush(stdout);
#endif
}

void Controller::em_update_type()
{
  if(! param2fixed_["types"]){
    // compute the contribution of each gene to each type weight
#pragma omp parallel for num_threads(nb_threads_)
    for(int g = 0; g < (int)genes_.size(); ++g){
      vector<double> new_type_prior_tmp(dim_, NaN);
      genes_[g].em_update_type(grid_wts_, subgroup_prior_, new_type_prior_tmp);
      for(size_t k = 0; k < dim_; ++k)
        type_genes_[k][g] = new_type_prior_tmp[k];
    }
#ifdef DEBUG
    if(verbose_ > 1){
      for(size_t k = 0; k < dim_; ++k)
        for(int g = 0; g < (int)genes_.size(); ++g)
          cout << "type " << k+1 << " gene" << g+1 << " " << type_genes_[k][g] << endl;
    }
#endif
    
    // average the contribution of all genes per type
    for(size_t k = 0; k < dim_; ++k)
      new_type_prior_[k] = log10_weighted_sum(&(type_genes_[k][0]),
                                              &(gene_wts_ones_[0]),
                                              genes_.size())
        + log10(type_prior_[k]);
    
    // compute the normalization constant (Lagrange multiplier)
    double l10_denom = log10_weighted_sum(&(new_type_prior_[0]),
                                          &(type_wts_ones_[0]),
                                          dim_);
    
    // compute each new type weight
    for(size_t k = 0; k < dim_; ++k)
      new_type_prior_[k] = pow(10, (new_type_prior_[k] - l10_denom));
  } // end of "if types not fixed"
  else
    new_type_prior_ = type_prior_;
  
#ifdef DEBUG
  for(size_t k = 0; k < dim_; ++k)
    fprintf(stdout, "pi_%zu old=%.4e new=%.4e\n", k+1, type_prior_[k], new_type_prior_[k]);
  fflush(stdout);
#endif
}

void Controller::em_update_subgroup()
{
  if(! param2fixed_["subgroups-per-type"]){
    // compute the contribution of each gene to each "subgroup per type" weight
#pragma omp parallel for num_threads(nb_threads_)
    for(int g = 0; g < (int)genes_.size(); ++g){
      vector<vector<double> > exp_gpkls_gene(dim_, vector<double>(nb_subgroups_, NaN));
      vector<double> exp_gpkl_gene(dim_, NaN);
      genes_[g].em_update_subgroup(grid_wts_, subgroup_prior_, exp_gpkls_gene, exp_gpkl_gene);
      for(size_t k = 0; k < dim_; ++k){
        subgroup_denom_genes_[k][g] = exp_gpkl_gene[k];
        for(size_t s = 0; s < nb_subgroups_; ++s)
          subgroup_num_genes_[k][s][g] = exp_gpkls_gene[k][s];
      }
    }
    
    // average the contribution of all genes per "subgroup per type"
    for(size_t k = 0; k < dim_; ++k)
      for(size_t s = 0; s < nb_subgroups_; ++s)
        new_subgroup_prior_[k][s] = log10_weighted_sum(&(subgroup_num_genes_[k][s][0]),
                                                       &(gene_wts_ones_[0]),
                                                       genes_.size())
          + log10(subgroup_prior_[k][s]);
    
    // compute the normalization constant
    vector<double> l10_denom(dim_, NaN);
    for(size_t k = 0; k < dim_; ++k)
      l10_denom[k] = log10_weighted_sum(&(subgroup_denom_genes_[k][0]),
                                        &(gene_wts_ones_[0]),
                                        genes_.size());
    
    // compute each new weight per "subgroup per type"
    for(size_t k = 0; k < dim_; ++k)
      for(size_t s = 0; s < nb_subgroups_; ++s)
        new_subgroup_prior_[k][s] = pow(10, (new_subgroup_prior_[k][s] - l10_denom[k]));
  } // end of "if subgroups per type not fixed"
  else
    new_subgroup_prior_ = subgroup_prior_;
  
#ifdef DEBUG
  for(size_t k = 0; k < dim_; ++k)
    for(size_t s = 0; s < nb_subgroups_; ++s)
      fprintf(stdout, "q_%zu%zu old=%.4e new=%.4e\n", k+1, s+1, subgroup_prior_[k][s], new_subgroup_prior_[k][s]);
  fflush(stdout);
#endif
}

void Controller::em_update_grid()
{
  if(! param2fixed_["grid-points"]){
    // compute the contribution of each gene to each grid weight
    if(model_ == "configs"){
#pragma omp parallel for num_threads(nb_threads_)
      for(int g = 0; g < (int)genes_.size(); ++g){
        vector<double> new_grid_wts_tmp(grid_size_, NaN);
        genes_[g].em_update_grid(config_prior_, new_grid_wts_tmp);
        for(size_t l = 0; l < grid_size_; ++l)
          grid_genes_[l][g] = new_grid_wts_tmp[l];
      }
    }
    else if(model_ == "types"){
#pragma omp parallel for num_threads(nb_threads_)
      for(int g = 0; g < (int)genes_.size(); ++g){
        vector<double> new_grid_wts_tmp(grid_size_, NaN);
        genes_[g].em_update_grid(type_prior_, subgroup_prior_, new_grid_wts_tmp);
        for(size_t l = 0; l < grid_size_; ++l)
          grid_genes_[l][g] = new_grid_wts_tmp[l];
      }
    }
    
    // average the contribution of all genes per grid weight
    for(size_t l = 0; l < grid_size_; ++l)
      new_grid_wts_[l] = log10_weighted_sum(&(grid_genes_[l][0]),
                                            &(gene_wts_ones_[0]),
                                            genes_.size())
        + log10(grid_wts_[l]);
    
    // compute the normalization constant (Lagrange multiplier)
    double l10_denom = log10_weighted_sum(&(new_grid_wts_[0]),
                                          &(grid_wts_ones_[0]),
                                          grid_size_);
    
    // compute each new grid weight
    for(size_t l = 0; l < grid_size_; ++l)
      new_grid_wts_[l] = pow(10, (new_grid_wts_[l] - l10_denom));
  } // end of "if fixed grids"
  else
    new_grid_wts_ = grid_wts_;
}

void Controller::update_params(int squaremstep=0)
{
  pi0_ = new_pi0_;
  if(model_ == "configs"){
    if(dim_ > 1)
      for(size_t k = 0; k < dim_; ++k)
        config_prior_[k] = new_config_prior_[k];
  }
  else if(model_ == "types"){
    for(size_t k = 0; k < dim_; ++k){
      type_prior_[k] = new_type_prior_[k];
      for(size_t s = 0; s < nb_subgroups_; ++s)
        subgroup_prior_[k][s] = new_subgroup_prior_[k][s];
    }
  }
  for(size_t l = 0; l < grid_size_; ++l)
    grid_wts_[l] = new_grid_wts_[l];
  
  // for(size_t i = 0; i < genes_.size(); ++i)
  //  genes_[i].update_snp_prior();
    if (squaremstep!=0){
        switch(squaremstep){
            case 1:{
                log10_obs_lik_1=new_log10_obs_lik_;
                pi0_1 = new_pi0_;
                if(model_ == "configs"){
                    if(dim_ > 1)
                        for(size_t k = 0; k < dim_; ++k)
                            config_prior_1[k] = new_config_prior_[k];
                }
                else if(model_ == "types"){
                    for(size_t k = 0; k < dim_; ++k){
                        type_prior_1[k] = new_type_prior_[k];
                        for(size_t s = 0; s < nb_subgroups_; ++s)
                            subgroup_prior_1[k][s] = new_subgroup_prior_[k][s];
                    }
                }
                
                for(size_t l = 0; l < grid_size_; ++l)
                    grid_wts_1[l] = new_grid_wts_[l];
                
            }
            case 2:{
                log10_obs_lik_2=new_log10_obs_lik_;
                pi0_2 = new_pi0_;
                if(model_ == "configs"){
                    if(dim_ > 1)
                        for(size_t k = 0; k < dim_; ++k)
                            config_prior_2[k] = new_config_prior_[k];
                }
                else if(model_ == "types"){
                    for(size_t k = 0; k < dim_; ++k){
                        type_prior_2[k] = new_type_prior_[k];
                        for(size_t s = 0; s < nb_subgroups_; ++s)
                            subgroup_prior_2[k][s] = new_subgroup_prior_[k][s];
                    }
                }
                
                for(size_t l = 0; l < grid_size_; ++l)
                    grid_wts_2[l] = new_grid_wts_[l];
            }
        }
    }
}


void Controller::show_state_EM(const size_t & iter)
{
  fprintf(stdout, "iter %4zu", iter);
  
  if(iter == 0)
    fprintf(stdout, "  loglik %f", log10_obs_lik_);
  else
    fprintf(stdout, "  loglik %f", new_log10_obs_lik_);
  
  fprintf(stdout, "  pi0 %7.4e", pi0_);
  
  if(model_ == "configs"){
    fprintf(stdout, "  configs");
    if(iter == 0)
      for(size_t i = 0; i < dim_; ++i)
        fprintf(stdout, " %7.4e", config_prior_[i]);
    else
      for(size_t i = 0; i < dim_; ++i)
        fprintf(stdout, " %7.4e", new_config_prior_[i]);
  } // end of "if model == configs"
  else if(model_ == "types"){
    fprintf(stdout, "  types");
    if(iter == 0){
      for(size_t i = 0; i < dim_; ++i)
        fprintf(stdout, " %7.4e", type_prior_[i]);
      fprintf(stdout, "  subgroups-per-type");
      for(size_t i = 0; i < dim_; ++i)
        for(size_t j = 0; j < nb_subgroups_; ++j)
          fprintf(stdout, " %7.4e", subgroup_prior_[i][j]);
    }
    else{ // iter > 0
      for(size_t i = 0; i < dim_; ++i)
        fprintf(stdout, " %7.4e", new_type_prior_[i]);
      fprintf(stdout, "  subgroups-per-type");
      for(size_t i = 0; i < dim_; ++i)
        for(size_t j = 0; j < nb_subgroups_; ++j)
          fprintf(stdout, " %7.4e", new_subgroup_prior_[i][j]);
    }
  } // end of "if model == types"
  
  fprintf(stdout, "  grid-points");
  if(iter == 0)
    for(size_t i = 0; i < grid_size_; ++i)
      fprintf(stdout, " %7.4e", grid_wts_[i]);
  else
    for(size_t i = 0; i < grid_size_; ++i)
      fprintf(stdout, " %7.4e", new_grid_wts_[i]);
  
  fprintf(stdout, "\n");
  fflush(stdout);
}

void Controller::run_EM_classic()
{
  size_t iter = 0;
  log10_obs_lik_ = compute_log10_obs_lik(pi0_, grid_wts_, config_prior_,
                                         type_prior_, subgroup_prior_, true);
  if(verbose_)
    show_state_EM(iter);
  
  while(true){
    ++iter;
    
    em_update_pi0();
#ifdef DEBUG
    double tmp_log10_obs_lik = compute_log10_obs_lik(new_pi0_, grid_wts_,
                                                     config_prior_,
                                                     type_prior_,
                                                     subgroup_prior_, false);
    fprintf(stdout, "log obslik after pi0: %f\n", tmp_log10_obs_lik);
    fflush(stdout);
    if(tmp_log10_obs_lik < log10_obs_lik_){
      fprintf(stderr, "ERROR: observed log-likelihood is decreasing (%f < %f)\n" ,
              tmp_log10_obs_lik, log10_obs_lik_);
      exit(EXIT_FAILURE);
    }
#endif
    
    if(model_ == "configs"){
      em_update_config();
#ifdef DEBUG
      tmp_log10_obs_lik = compute_log10_obs_lik(new_pi0_, grid_wts_,
                                                new_config_prior_,
                                                type_prior_,
                                                subgroup_prior_, false);
      fprintf(stdout, "log obslik after configs: %f\n", tmp_log10_obs_lik);
      fflush(stdout);
      if(tmp_log10_obs_lik < log10_obs_lik_){
        fprintf(stderr, "ERROR: observed log-likelihood is decreasing (%f < %f)\n" ,
                tmp_log10_obs_lik, log10_obs_lik_);
        exit(EXIT_FAILURE);
      }
#endif
    }
    else if(model_ == "types"){
      em_update_type();
#ifdef DEBUG
      tmp_log10_obs_lik = compute_log10_obs_lik(new_pi0_, grid_wts_,
                                                config_prior_,
                                                new_type_prior_,
                                                subgroup_prior_, false);
      fprintf(stdout, "log obslik after types: %f\n", tmp_log10_obs_lik);
      fflush(stdout);
      if(tmp_log10_obs_lik < log10_obs_lik_){
        fprintf(stderr, "ERROR: observed log-likelihood is decreasing (%f < %f)\n" ,
                tmp_log10_obs_lik, log10_obs_lik_);
        exit(EXIT_FAILURE);
      }
#endif
      em_update_subgroup();
#ifdef DEBUG
      tmp_log10_obs_lik = compute_log10_obs_lik(new_pi0_, grid_wts_,
                                                config_prior_,
                                                new_type_prior_,
                                                new_subgroup_prior_, false);
      fprintf(stdout, "log obslik after subgroups-per-type: %f\n", tmp_log10_obs_lik);
      fflush(stdout);
      if(tmp_log10_obs_lik < log10_obs_lik_){
        fprintf(stderr, "ERROR: observed log-likelihood is decreasing (%f < %f)\n" ,
                tmp_log10_obs_lik, log10_obs_lik_);
        exit(EXIT_FAILURE);
      }
#endif
    }
    
    em_update_grid();
    em_update_snp_prior();
    update_params();
    new_log10_obs_lik_ = compute_log10_obs_lik(new_pi0_, new_grid_wts_,
                                               new_config_prior_,
                                               new_type_prior_,
                                               new_subgroup_prior_, true);
    if(verbose_)
      show_state_EM(iter);
    
    if(new_log10_obs_lik_ < log10_obs_lik_){
      fprintf(stderr, "ERROR: observed log-likelihood is decreasing (%f < %f)\n" ,
              new_log10_obs_lik_, log10_obs_lik_);
      exit(EXIT_FAILURE);
    }
    if(fabs(new_log10_obs_lik_ - log10_obs_lik_) < thresh_ ||
       (max_nb_iters_ != string::npos && iter == max_nb_iters_ - 1))
      break;
    
    log10_obs_lik_ = new_log10_obs_lik_;
  }
  
  // do a last update to get same results as William's initial implementation
  ++iter;
  em_update_pi0();
  if(model_ == "configs")
    em_update_config();
  else if(model_ == "types"){
    em_update_type();
    em_update_subgroup();
  }
  em_update_grid();
  em_update_snp_prior();
  update_params();
  new_log10_obs_lik_ = compute_log10_obs_lik(new_pi0_, new_grid_wts_,
                                             new_config_prior_,
                                             new_type_prior_,
                                             new_subgroup_prior_, true);
  if(verbose_)
    show_state_EM(iter);
}


//A wrapper function serving as the fixed point in SQUAREM
void Controller::run_EM_fixedpoint(){
    em_update_pi0();
#ifdef DEBUG
    double tmp_log10_obs_lik = compute_log10_obs_lik(new_pi0_, grid_wts_,
                                                     config_prior_,
                                                     type_prior_,
                                                     subgroup_prior_, false);
    fprintf(stdout, "log obslik after pi0: %f\n", tmp_log10_obs_lik);
    fflush(stdout);
    if(tmp_log10_obs_lik < log10_obs_lik_){
        fprintf(stderr, "ERROR: observed log-likelihood is decreasing (%f < %f)\n" ,
                tmp_log10_obs_lik, log10_obs_lik_);
        exit(EXIT_FAILURE);
    }
#endif
    
    if(model_ == "configs"){
        em_update_config();
#ifdef DEBUG
        tmp_log10_obs_lik = compute_log10_obs_lik(new_pi0_, grid_wts_,
                                                  new_config_prior_,
                                                  type_prior_,
                                                  subgroup_prior_, false);
        fprintf(stdout, "log obslik after configs: %f\n", tmp_log10_obs_lik);
        fflush(stdout);
        if(tmp_log10_obs_lik < log10_obs_lik_){
            fprintf(stderr, "ERROR: observed log-likelihood is decreasing (%f < %f)\n" ,
                    tmp_log10_obs_lik, log10_obs_lik_);
            exit(EXIT_FAILURE);
        }
#endif
    }
    else if(model_ == "types"){
        em_update_type();
#ifdef DEBUG
        tmp_log10_obs_lik = compute_log10_obs_lik(new_pi0_, grid_wts_,
                                                  config_prior_,
                                                  new_type_prior_,
                                                  subgroup_prior_, false);
        fprintf(stdout, "log obslik after types: %f\n", tmp_log10_obs_lik);
        fflush(stdout);
        if(tmp_log10_obs_lik < log10_obs_lik_){
            fprintf(stderr, "ERROR: observed log-likelihood is decreasing (%f < %f)\n" ,
                    tmp_log10_obs_lik, log10_obs_lik_);
            exit(EXIT_FAILURE);
        }
#endif
        em_update_subgroup();
#ifdef DEBUG
        tmp_log10_obs_lik = compute_log10_obs_lik(new_pi0_, grid_wts_,
                                                  config_prior_,
                                                  new_type_prior_,
                                                  new_subgroup_prior_, false);
        fprintf(stdout, "log obslik after subgroups-per-type: %f\n", tmp_log10_obs_lik);
        fflush(stdout);
        if(tmp_log10_obs_lik < log10_obs_lik_){
            fprintf(stderr, "ERROR: observed log-likelihood is decreasing (%f < %f)\n" ,
                    tmp_log10_obs_lik, log10_obs_lik_);
            exit(EXIT_FAILURE);
        }
#endif
    }
    
    em_update_grid();
    em_update_snp_prior();
    //update_params();
    new_log10_obs_lik_ = compute_log10_obs_lik(new_pi0_, new_grid_wts_,
                                               new_config_prior_,
                                               new_type_prior_,
                                               new_subgroup_prior_, true);
}


void Controller::run_EM_square()
{
  // TODO
    size_t iter = 0;
    double sr2_scalar,sv2_scalar,srv_scalar,squaremalpha_,stepmin,stepmax;
    bool extrap;
    stepmin=stepmin0_;
    stepmax=stepmax0_;
    
    log10_obs_lik_ = compute_log10_obs_lik(pi0_, grid_wts_, config_prior_,
                                           type_prior_, subgroup_prior_, true);
    if(verbose_)
        show_state_EM(iter);

    //main squarem loop
    while(true){
        ++iter;
        sr2_scalar=0;
        sv2_scalar=0;
        srv_scalar=0;
        extrap = true;
        //storing fixed point vector
        log10_obs_lik_0=log10_obs_lik_;
        pi0_0=pi0_;
        grid_wts_0=grid_wts_;
        config_prior_0=config_prior_;
        type_prior_0=type_prior_;
        subgroup_prior_0=subgroup_prior_;
        
        run_EM_fixedpoint();
        update_params(1);
        if(fabs(log10_obs_lik_1 - log10_obs_lik_) < thresh_ ||
           (max_nb_iters_ != string::npos && iter == max_nb_iters_ - 1))
            break;
        
        ++iter;
        run_EM_fixedpoint();
        update_params(2);
        if(fabs(log10_obs_lik_2 - log10_obs_lik_1) < thresh_ ||
           (max_nb_iters_ != string::npos && iter == max_nb_iters_ - 1))
            break;
        
        //Calculating alpha
        sr2_scalar+=pow(pi0_1-pi0_0,2);
        sv2_scalar+=pow(pi0_2-2*pi0_1+pi0_0,2);
        srv_scalar+=(pi0_2-2*pi0_1+pi0_0)*(pi0_1-pi0_0);

        if(model_ == "configs"){
            if(dim_ > 1)
                for(size_t k = 0; k < dim_; ++k){
                    sr2_scalar+=pow(config_prior_1[k]-config_prior_0[k],2);
                    sv2_scalar+=pow(config_prior_2[k]-2*config_prior_1[k]+config_prior_0[k],2);
                    srv_scalar+=(config_prior_2[k]-2*config_prior_1[k]+config_prior_0[k])*(config_prior_1[k]-config_prior_0[k]);
                }
        }
        else if(model_ == "types"){
            for(size_t k = 0; k < dim_; ++k){
                sr2_scalar+=pow(type_prior_1[k]-type_prior_0[k],2);
                sv2_scalar+=pow(type_prior_2[k]-2*type_prior_1[k]+type_prior_0[k],2);
                srv_scalar+=(type_prior_2[k]-2*type_prior_1[k]+type_prior_0[k])*(type_prior_1[k]-type_prior_0[k]);
                
                for(size_t s = 0; s < nb_subgroups_; ++s){
                    sr2_scalar+=pow(subgroup_prior_1[k][s]-subgroup_prior_0[k][s],2);
                    sv2_scalar+=pow(subgroup_prior_2[k][s]-2*subgroup_prior_1[k][s]+subgroup_prior_0[k][s],2);
                    srv_scalar+=(subgroup_prior_2[k][s]-2*subgroup_prior_1[k][s]+subgroup_prior_0[k][s])*(subgroup_prior_1[k][s]-subgroup_prior_0[k][s]);
                }
            }
        }
        
        for(size_t l = 0; l < grid_size_; ++l){
            sr2_scalar+=pow(grid_wts_1[l]-grid_wts_0[l],2);
            sv2_scalar+=pow(grid_wts_2[l]-2*grid_wts_1[l]+grid_wts_0[l],2);
            srv_scalar+=(grid_wts_2[l]-2*grid_wts_1[l]+grid_wts_0[l])*(grid_wts_1[l]-grid_wts_0[l]);
        }
        
        
        switch (method_){
            case 1:squaremalpha_=-srv_scalar/sv2_scalar;
            case 2:squaremalpha_=-sr2_scalar/srv_scalar;
            case 3:squaremalpha_=sqrt(sr2_scalar/sv2_scalar);
        }
        squaremalpha_=std::max(stepmin,std::min(stepmax,squaremalpha_));
        
        //Assigning new value
        pi0_=pi0_0+2.0*squaremalpha_*(pi0_1-pi0_0)+pow(squaremalpha_,2)*(pi0_2-2*pi0_1+pi0_0);
        if(model_ == "configs"){
            if(dim_ > 1)
                for(size_t k = 0; k < dim_; ++k){
                    config_prior_[k]=config_prior_0[k]+2.0*squaremalpha_*(config_prior_1[k]-config_prior_0[k])
                    +pow(squaremalpha_,2)*(config_prior_2[k]-2*config_prior_1[k]+config_prior_0[k]);
                }
        }
        else if(model_ == "types"){
            for(size_t k = 0; k < dim_; ++k){
                type_prior_[k]=type_prior_0[k]+2.0*squaremalpha_*(type_prior_1[k]-type_prior_0[k])
                +pow(squaremalpha_,2)*(type_prior_2[k]-2*type_prior_1[k]+type_prior_0[k]);
                for(size_t s = 0; s < nb_subgroups_; ++s){
                    subgroup_prior_[k][s]=subgroup_prior_0[k][s]+2.0*squaremalpha_*(subgroup_prior_1[k][s]-subgroup_prior_0[k][s])
                    +pow(squaremalpha_,2)*(subgroup_prior_2[k][s]-2*subgroup_prior_1[k][s]+subgroup_prior_0[k][s]);
                }
            }
        }
        
        for(size_t l = 0; l < grid_size_; ++l){
            grid_wts_[l]=grid_wts_0[l]+2*squaremalpha_*(grid_wts_1[l]-grid_wts_0[l])
            +pow(squaremalpha_,2)*(grid_wts_2[l]-2*grid_wts_1[l]+grid_wts_0[l]);
        }
        
        
        //Stalization
        if(std::abs(squaremalpha_-1)>0.01){
            try{run_EM_fixedpoint(); ++iter;}
            catch(...){
                log10_obs_lik_=log10_obs_lik_2;
                pi0_=pi0_2;
                grid_wts_=grid_wts_2;
                config_prior_=config_prior_2;
                type_prior_=type_prior_2;
                subgroup_prior_=subgroup_prior_2;
                if(squaremalpha_ == stepmax){
                    stepmax=std::max(stepmax0_,stepmax/mstep_);
                }
                squaremalpha_=1;
                extrap = false;
                if(squaremalpha_ == stepmax){stepmax=mstep_*stepmax;}
                if(stepmin<0 && squaremalpha_==stepmin){stepmin=mstep_*stepmin;}
                if(verbose_){std::cout<<"Objective fn: "<<log10_obs_lik_<<"  Extrapolation: "<<extrap<<"  Steplength: "<<squaremalpha_<<std::endl;}
                continue;
            }
            update_params(0);//assign _new to _
        }
        
        if (new_log10_obs_lik_<log10_obs_lik_0 - objfninc_){//maximize the objective function,less stringent criteria when objfninc_=1
            log10_obs_lik_=log10_obs_lik_2;
            pi0_=pi0_2;
            grid_wts_=grid_wts_2;
            config_prior_=config_prior_2;
            type_prior_=type_prior_2;
            subgroup_prior_=subgroup_prior_2;
            if(squaremalpha_ == stepmax){
                stepmax=std::max(stepmax0_,stepmax/mstep_);
            }
            squaremalpha_=1;
            extrap = false;
        }
        
        if(squaremalpha_==stepmax){stepmax=mstep_*stepmax;}
        if(stepmin<0 && squaremalpha_==stepmin){stepmin=mstep_*stepmin;}
        if(verbose_){std::cout<<"Objective fn: "<<log10_obs_lik_<<"  Extrapolation: "<<extrap<<"  Steplength: "<<squaremalpha_<<std::endl;}
    
    }
    //end of squarem loop

    // do a last update to get same results as William's initial implementation
    ++iter;
    em_update_pi0();
    if(model_ == "configs")
        em_update_config();
    else if(model_ == "types"){
        em_update_type();
        em_update_subgroup();
    }
    em_update_grid();
    em_update_snp_prior();
    update_params();
    new_log10_obs_lik_ = compute_log10_obs_lik(new_pi0_, new_grid_wts_,
                                               new_config_prior_,
                                               new_type_prior_,
                                               new_subgroup_prior_, true);
    if(verbose_)
        show_state_EM(iter);
}

void Controller::run_EM()
{
  if(verbose_ > 0){
    cout << "run EM algorithm (";
    if(squarem_)
      cout << "square";
    else
      cout << "classic";
    cout << ") ..." << endl << flush;
  }
  time_t startRawTime, endRawTime;
  time (&startRawTime);
  
  if(squarem_)
    run_EM_square();
  else
    run_EM_classic();
  
  time (&endRawTime);
  cout << "EM ran for " << getElapsedTime(startRawTime, endRawTime) << endl << flush;
}

void Controller::estimate_profile_ci_pi0(const double & tick,
                                         const double & max_l10_obs_lik)
{
  double pi0_mle = pi0_, curr_log10_obs_lik;
  left_pi0_ = pi0_;
  right_pi0_ = pi0_;
  
  while(left_pi0_ >= 0){
    left_pi0_ -= tick;
    if(left_pi0_ < 0){
      left_pi0_ = 0;
      break;
    }
    pi0_ = left_pi0_;
    curr_log10_obs_lik = compute_log10_obs_lik(pi0_, new_grid_wts_,
                                               new_config_prior_,
                                               new_type_prior_,
                                               new_subgroup_prior_, true);
    if(curr_log10_obs_lik / log10(exp(1))
       < max_l10_obs_lik / log10(exp(1)) - 2.0){  
      left_pi0_ += tick;
      break;
    }    
  }
  while(right_pi0_ <= 1){
    right_pi0_ += tick;
    if(right_pi0_ > 1){
      right_pi0_ = 1;
      break;
    }
    pi0_ = right_pi0_;  
    curr_log10_obs_lik = compute_log10_obs_lik(pi0_, new_grid_wts_,
                                               new_config_prior_,
                                               new_type_prior_,
                                               new_subgroup_prior_, true);
    if(curr_log10_obs_lik / log10(exp(1))
       < max_l10_obs_lik / log10(exp(1)) - 2.0){
      right_pi0_ -= tick;
      break;
    }
  }
  
  pi0_ = pi0_mle;
}

void Controller::estimate_profile_ci_configs(const double & tick,
                                             const double & max_l10_obs_lik)
{
  double curr_log10_obs_lik;
  vector<double> config_mle(config_prior_);
  
  for(size_t i = 0 ; i < config_prior_.size(); ++i){
    right_configs_[i] = config_mle[i];
    left_configs_[i] = config_mle[i];
    double cp_mle = config_mle[i];
    double st = 1 - cp_mle;
    while(left_configs_[i] >= 0){
      left_configs_[i] -= tick;
      if(left_configs_[i] < 0){
        left_configs_[i] = 0;
        break;
      }
      double diff = cp_mle - left_configs_[i];
      for(size_t j = 0; j < config_prior_.size(); ++j){
        if(j == i)
          continue;
        config_prior_[j] = config_mle[j] + diff*config_mle[j]/st;
      }
      config_prior_[i] = left_configs_[i];
      curr_log10_obs_lik = compute_log10_obs_lik(new_pi0_, new_grid_wts_,
                                                 config_prior_,
                                                 new_type_prior_,
                                                 new_subgroup_prior_, true);
      if(curr_log10_obs_lik / log10(exp(1))
         < max_l10_obs_lik / log10(exp(1)) - 2.0){
        left_configs_[i] += tick;
        break;
      }
    }
    while(right_configs_[i] <= 1){
      right_configs_[i] += tick;
      if(right_configs_[i] > 1){
        right_configs_[i] = 1;
        break;
      }
      double diff = cp_mle - right_configs_[i];
      for(size_t j = 0; j < config_prior_.size(); ++j){
        if(j == i)
          continue;
        config_prior_[j] = config_mle[j] + diff*config_mle[j]/st;
      }
      config_prior_[i] = right_configs_[i];
      curr_log10_obs_lik = compute_log10_obs_lik(new_pi0_, new_grid_wts_,
                                                 config_prior_,
                                                 new_type_prior_,
                                                 new_subgroup_prior_, true);
      if(curr_log10_obs_lik / log10(exp(1))
         < max_l10_obs_lik / log10(exp(1)) - 2.0){
        right_configs_[i] -= tick;
        break;
      }
    }
  }
  
  config_prior_ = config_mle;
}

void Controller::estimate_profile_ci_types(const double & tick,
                                           const double & max_l10_obs_lik)
{
  vector<double> type_mle(type_prior_);
  
  for(size_t k = 0 ; k < dim_; ++k){
    right_types_[k] = type_mle[k];
    left_types_[k] = type_mle[k];
  }
  
  type_prior_ = type_mle;
}

void Controller::estimate_profile_ci_subgroups(const double & tick,
                                               const double & max_l10_obs_lik)
{
  vector<vector<double> > subgroup_mle(dim_, vector<double>(nb_subgroups_, NaN));
  for(size_t k = 0; k < dim_; ++k)
    subgroup_mle[k] = subgroup_prior_[k];
  
  for(size_t k = 0 ; k < dim_; ++k){
    for(size_t s = 0; s < nb_subgroups_; ++s){
      right_subgroups_[k][s] = subgroup_mle[k][s];
      left_subgroups_[k][s] = subgroup_mle[k][s];
    }
  }
  
  subgroup_prior_ = subgroup_mle;
}

void Controller::estimate_profile_ci_grids(const double & tick,
                                           const double & max_l10_obs_lik)
{
  double curr_log10_obs_lik;
  vector<double> grid_mle(grid_wts_);
  
  for(size_t i = 0; i < grid_wts_.size(); ++i){
    right_grids_[i] = grid_mle[i];
    left_grids_[i] = grid_mle[i];
    double cp_mle = grid_mle[i];
    double st = 1 - cp_mle;
    while(left_grids_[i] >= 0){
      left_grids_[i] -= tick;
      if(left_grids_[i] < 0){
        left_grids_[i] = 0;
        break;
      }
      double diff = cp_mle - left_grids_[i];
      for(size_t j = 0; j < grid_wts_.size(); ++j){
        if(j == i)
          continue;
        grid_wts_[j] = grid_mle[j] + diff*grid_mle[j]/st;
      }
      grid_wts_[i] = left_grids_[i];
      curr_log10_obs_lik = compute_log10_obs_lik(new_pi0_, grid_wts_,
                                                 new_config_prior_,
                                                 new_type_prior_,
                                                 new_subgroup_prior_, true);
      if(curr_log10_obs_lik / log10(exp(1)) 
         < max_l10_obs_lik / log10(exp(1)) - 2.0){
        left_grids_[i] += tick;
        break;
      }
    }
    while(right_grids_[i] <= 1){
      right_grids_[i] += tick;
      if(right_grids_[i] > 1){
        right_grids_[i] = 1;
        break;
      }
      double diff = cp_mle - right_grids_[i];
      for(size_t j = 0; j < grid_wts_.size(); ++j){
        if(j == i)
          continue;
        grid_wts_[j] = grid_mle[j] + diff*grid_mle[j]/st;
      }
      grid_wts_[i] = right_grids_[i];
      curr_log10_obs_lik = compute_log10_obs_lik(new_pi0_, grid_wts_,
                                                 new_config_prior_,
                                                 new_type_prior_,
                                                 new_subgroup_prior_, true);
      if(curr_log10_obs_lik / log10(exp(1))
         < max_l10_obs_lik / log10(exp(1)) - 2.0){
        right_grids_[i] -= tick;
        break;
      }
    }
  }
  
  grid_wts_ = grid_mle;
}

void Controller::estimate_profile_ci()
{
  if(verbose_ > 0)
    cout << "compute profile-likelihood confidence intervals ..." << flush;
  clock_t startTime = clock();
  
  double tick = 0.001,
    max_l10_obs_lik = compute_log10_obs_lik(new_pi0_, new_grid_wts_,
                                            new_config_prior_,
                                            new_type_prior_,
                                            new_subgroup_prior_, true);
  
  estimate_profile_ci_pi0(tick, max_l10_obs_lik);
  
  if(model_ == "configs")
    estimate_profile_ci_configs(tick, max_l10_obs_lik);
  else if(model_ == "types"){
    estimate_profile_ci_types(tick, max_l10_obs_lik);
    estimate_profile_ci_subgroups(tick, max_l10_obs_lik);
  }
  
  estimate_profile_ci_grids(tick, max_l10_obs_lik);
  
  if(verbose_ > 0)
    cout << " (" << getElapsedTime(startTime) << " sec)" << endl;
}

// choose K via the ICL (based on the augmented log-lik minus a BIC-like penalty)
void Controller::compute_log10_ICL()
{
  if(model_ == "types"){
    if(verbose_ > 0)
      cout << "compute the ICL ..." << flush;
    clock_t startTime = clock();
    
    double log10_icl = - (double)(1 + grid_size_ + dim_ * (nb_subgroups_ + 1)) / 2.0
      * log10(genes_.size());
    
#pragma omp parallel for num_threads(nb_threads_) reduction(+:log10_icl)
    for(int g = 0; g < (int)genes_.size(); ++g)
      log10_icl += genes_[g].compute_log10_aug_lik(pi0_, grid_wts_, 
                                                   type_prior_,
                                                   subgroup_prior_);
    
    if(verbose_ > 0)
      cout << " log10(ICL)=" << setprecision(4) << scientific << log10_icl << fixed
           << " (" << getElapsedTime(startTime) << " sec)"
           << endl;
  }
}

void Controller::compute_posterior()
{
  if(verbose_ > 0)
    cout << "compute posteriors ..." << flush;
  clock_t startTime = clock();
  
#pragma omp parallel for num_threads(nb_threads_)
  for(int i = 0; i < (int)genes_.size(); ++i)
    genes_[i].compute_posterior(pi0_, grid_wts_, config_prior_);
  
  if(verbose_ > 0)
    cout << " (" << getElapsedTime(startTime) << " sec)" << endl;
}

void Controller::save_result(const string & out_file,
                             const bool & skip_bf)
{
  if(verbose_ > 0)
    cout << "save the results in " << out_file << " ..." << flush;
  clock_t startTime = clock();
  
  gzFile stream;
  openFile(out_file, stream, "wb");
  
  stringstream txt;
  size_t nb_lines = 0;
  
  // first lines are comments (#) with MLEs and CI of all parameters
  txt << "#param\tmle\tleft.ci\tright.ci\tfixed" << endl;
  ++nb_lines;
  gzwriteLine(stream, txt.str(), out_file, nb_lines);
  txt.str("");
  txt << setprecision(4);
  txt << "#pi0\t" << scientific << pi0_ << "\t" << left_pi0_ << "\t" << right_pi0_
      << "\t" << boolalpha << param2fixed_["pi0"] << endl;
  ++nb_lines;
  gzwriteLine(stream, txt.str(), out_file, nb_lines);
  if(model_ == "configs"){
    for(size_t k = 0; k < dim_; ++k){
      txt.str("");
      txt << "#config." << config_names_[k].c_str()
          << "\t" << config_prior_[k]
          << "\t" << left_configs_[k]
          << "\t" << right_configs_[k]
          << "\t" << param2fixed_["configs"]
          << endl;
      ++nb_lines;
      gzwriteLine(stream, txt.str(), out_file, nb_lines);
    }
  }
  else if(model_ == "types"){
    for(size_t k = 0; k < dim_; ++k){
      txt.str("");
      txt << "#type." << k+1
          << "\t" << type_prior_[k]
          << "\t" << left_types_[k]
          << "\t" << right_types_[k]
          << "\t" << param2fixed_["types"]
          << endl;
      ++nb_lines;
      gzwriteLine(stream, txt.str(), out_file, nb_lines);
    }
    for(size_t k = 0; k < dim_; ++k){
      for(size_t s = 0; s < nb_subgroups_; ++s){
        txt.str("");
        txt << "#subgroup." << k+1 << "-" << s+1
            << "\t" << subgroup_prior_[k][s]
            << "\t" << left_subgroups_[k][s]
            << "\t" << right_subgroups_[k][s]
            << "\t" << param2fixed_["subgroups-per-type"]
            << endl;
        ++nb_lines;
        gzwriteLine(stream, txt.str(), out_file, nb_lines);
      }
    }
  }
  for(size_t l = 0; l < grid_wts_.size(); ++l){
    txt.str("");
    txt << "#grid." << l+1
        << "\t" << grid_wts_[l]
        << "\t" << left_grids_[l]
        << "\t" << right_grids_[l]
        << "\t" << param2fixed_["grid-points"]
        << endl;
    ++nb_lines;
    gzwriteLine(stream, txt.str(), out_file, nb_lines);
  }
  
  // the rest of the file contains BFs and posteriors
  if(! skip_bf){
    txt.str("");
    txt << "gene\tgene.posterior.prob\tgene.log10.bf\tsnp\tsnp.log10.bf";
    if(model_ == "configs")
      for(size_t k = 0; k < dim_; ++k)
        txt << "\tlog10.bf." << config_names_[k];
    else if(model_ == "types")
      for(size_t k = 0; k < dim_; ++k)
        txt << "\tlog10.bf." << k+1;
    txt << "\n";
    ++nb_lines;
    gzwriteLine(stream, txt.str(), out_file, nb_lines);
    
    for(size_t g = 0; g < genes_.size(); ++g){ // loop over gene
      for(size_t p = 0; p < genes_[g].snps_.size(); ++p){ // loop over snp
        txt.str("");
        txt << setprecision(4)
            << genes_[g].name_
            << "\t"
            << genes_[g].post_prob_gene_;
	
        if(model_ == "configs")
          txt << "\t"
              << genes_[g].compute_log10_BF(grid_wts_, config_prior_, true);
        else if(model_ == "types")
          txt << "\t"
              << genes_[g].compute_log10_BF(grid_wts_, type_prior_, subgroup_prior_, true);
	
        txt << "\t"
            << genes_[g].snps_[p].name_;
	
        if(model_ == "configs"){
          txt << "\t"
              << genes_[g].snps_[p].compute_log10_BF(grid_wts_, config_prior_, true);
          for(size_t k = 0; k < dim_; ++k) // loop over configs
            txt << "\t"
                << genes_[g].snps_[p].compute_log10_config_BF(k, grid_wts_);
        }
        else if(model_ == "types"){
          txt << "\t"
              << genes_[g].snps_[p].compute_log10_BF(grid_wts_, type_prior_, subgroup_prior_, true);
          for(size_t k = 0; k < dim_; ++k) // loop over types
            txt << "\t"
                << genes_[g].snps_[p].compute_log10_type_BF(k, grid_wts_, subgroup_prior_);
        }
	
        txt << "\n";
        ++nb_lines;
        gzwriteLine(stream, txt.str(), out_file, nb_lines);
      }
    }
  }
  
  closeFile(out_file, stream);
  
  if(verbose_ > 0)
    cout << " (" << getElapsedTime(startTime) << " sec)" << endl;
}

//-----------------------------------------------------------------------------

/** \brief Display the help on stdout.
 */
void help(char ** argv)
{
  cout << "`" << argv[0] << "'"
       << " fits the hierarchical model of eQtlBma with an EM algorithm." << endl
       << endl
       << "Usage: " << argv[0] << " [OPTIONS] ..." << endl
       << endl
       << "Options:" << endl
       << "  -h, --help\tdisplay the help and exit" << endl
       << "  -V, --version\toutput version information and exit" << endl
       << "  -v, --verbose\tverbosity level (0/default=1/2/3)" << endl
       << "      --data\tinput data (usually output files from eqtlbma_bf)" << endl
       << "      --nsubgrp\tnumber of subgroups" << endl
       << "      --model\twhich model to fit (default=configs/types)" << endl
       << "      --dim\tdimension of the model (nb of active configs or types)" << endl
       << "      --ngrid\tnumber of grid points" << endl
       << "      --out\toutput file (gzipped)" << endl
       << "      --init\tfile for initialization" << endl
       << "\t\t3 columns: param<tab>value<tab>fixed (TRUE or FALSE)" << endl
       << "      --rand\trandom initialization" << endl
       << "      --seed\tseed used with --rand, otherwise use time" << endl
       << "      --thresh\tthreshold to stop the EM (default=0.05)" << endl
       << "      --maxit\tmaximum number of iterations (optional)" << endl
       << "\t\tuseful if wall-time limit (see also --init)" << endl
       << "      --sq\tuse square EM for speed-up" << endl
       << "      --thread\tnumber of threads (default=1)" << endl
       << "      --configs\tsubset of configurations to keep (e.g. \"1|3|1-3\")" << endl
       << "      --keepgen\tkeep 'general' ABFs (useful for BMAlite)" << endl
       << "      --getci\tcompute the confidence intervals (single thread, thus slow)" << endl
       << "      --getbf\tcompute the Bayes Factors using the estimated weights" << endl
       << "\t\tcan take some time, otherwise only the estimated weights are reported" << endl
       << "      --pi0\tfixed value for pi0 (pi0 hence won't be updated in the EM)" << endl
       << "      --ci\tfile with estimates of hyperparameters to only compute confidence intervals" << endl
       << endl;
}

/** \brief Display version and license information on stdout.
 */
void version(char ** argv)
{
  cout << argv[0] << " " << VERSION << endl
       << endl
       << "Copyright (C) 2012-2014 Xiaoquan Wen and Timothee Flutre." << endl
       << "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>" << endl
       << "This is free software; see the source for copying conditions.  There is NO" << endl
       << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl
       << endl
       << "Written by Xiaoquan Wen and Timothee Flutre." << endl;
}

/** \brief Parse the command-line arguments and check the values of the 
 *  compulsory ones.
 */
void parseCmdLine(
  int argc,
  char ** argv,
  string & file_pattern,
  size_t & nb_subgroups,
  string & model,
  size_t & dim,
  size_t & nb_grid_points,
  string & out_file,
  string & file_init,
  bool & rand_init,
  size_t & seed,
  double & thresh,
  size_t & max_nb_iters,
  bool & squarem,
  int & nb_threads,
  string & file_ci,
  vector<string> & configs_tokeep,
  bool & keep_gen_abfs,
  bool & skip_ci,
  bool & skip_bf,
  double & fixed_pi0,
  int & verbose)
{
  int c = 0;
  while(true){
    static struct option long_options[] = {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'V'},
      {"verbose", required_argument, 0, 'v'},
      {"data", required_argument, 0, 0},
      {"nsubgrp", required_argument, 0, 0},
      {"model", required_argument, 0, 0},
      {"dim", required_argument, 0, 0},
      {"ngrid", required_argument, 0, 0},
      {"out", required_argument, 0, 0},
      {"init", required_argument, 0, 0},
      {"rand", no_argument, 0, 0},
      {"seed", required_argument, 0, 0},
      {"thresh", required_argument, 0, 0},
      {"maxit", required_argument, 0, 0},
      {"squarem", no_argument, 0, 0},
      {"thread", required_argument, 0, 0},
      {"ci", required_argument, 0, 0},
      {"configs", required_argument, 0, 0},
      {"keepgen", no_argument, 0, 0},
      {"getci", no_argument, 0, 0},
      {"getbf", no_argument, 0, 0},
      {"pi0", required_argument, 0, 0},
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
      if(strcmp(long_options[option_index].name, "data") == 0){
        file_pattern = optarg;
        break;
      }
      if(strcmp(long_options[option_index].name, "nsubgrp") == 0){
        nb_subgroups = atol(optarg);
        break;
      }
      if(strcmp(long_options[option_index].name, "model") == 0){
        model = optarg;
        break;
      }
      if(strcmp(long_options[option_index].name, "dim") == 0){
        dim = atoi(optarg);
        break;
      }
      if(strcmp(long_options[option_index].name, "ngrid") == 0){
        nb_grid_points = atoi(optarg);
        break;
      }
      if(strcmp(long_options[option_index].name, "out") == 0){
        out_file = optarg;
        break;
      }
      if(strcmp(long_options[option_index].name, "init") == 0){
        file_init = optarg;
        break;
      }
      if(strcmp(long_options[option_index].name, "rand") == 0){
        rand_init = true;
        break;
      }
      if(strcmp(long_options[option_index].name, "seed") == 0){
        seed = atol(optarg);
        break;
      }
      if(strcmp(long_options[option_index].name, "thresh") == 0){
        thresh = atof(optarg);
        break;
      }
      if(strcmp(long_options[option_index].name, "maxit") == 0){
        max_nb_iters = atol(optarg);
        break;
      }
      if(strcmp(long_options[option_index].name, "sq") == 0){
        squarem = true;
        break;
      }
      if(strcmp(long_options[option_index].name, "thread") == 0){
        nb_threads = atoi(optarg);
        break;
      }
      if(strcmp(long_options[option_index].name, "ci") == 0){
        file_ci = optarg;
        break;
      }
      if(strcmp(long_options[option_index].name, "configs") == 0){
        split(optarg, "|", configs_tokeep);
        break;
      }
      if(strcmp(long_options[option_index].name, "keepgen") == 0){
        keep_gen_abfs = true;
        break;
      }
      if(strcmp(long_options[option_index].name, "getci") == 0){
        skip_ci = false;
        break;
      }
      if(strcmp(long_options[option_index].name, "getbf") == 0){
        skip_bf = false;
        break;
      }
      if(strcmp(long_options[option_index].name, "pi0") == 0){
        fixed_pi0 = atof(optarg);
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
         << "ERROR: missing compulsory option --data" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(model.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
         << "ERROR: missing compulsory option --model" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(model != "configs" && model != "types"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
         << "ERROR: --model " << model << " is invalid" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(nb_subgroups == string::npos){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
         << "ERROR: missing compulsory option --nsubgrp" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(dim == string::npos){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
         << "ERROR: missing compulsory option --dim" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(model == "configs" && dim >= (size_t) pow(2, (double)nb_subgroups)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
         << "ERROR: --model configs --nsubgrp " << nb_subgroups << " --dim "
         << dim << " is invalid" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(nb_grid_points == string::npos){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
         << "ERROR: missing compulsory option --ngrid" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(out_file.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
         << "ERROR: missing compulsory option --out" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(! file_init.empty() && ! doesFileExist(file_init)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
         << "ERROR: can't find " << file_init << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(rand_init && seed == string::npos)
    seed = getSeed();
  if(thresh <= 0.0){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
         << "ERROR: --thresh " << thresh << " is invalid" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(nb_threads <= 0)
    nb_threads = 1;
  if(! file_ci.empty() && ! doesFileExist(file_ci)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
         << "ERROR: can't find " << file_ci << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(! isNan(fixed_pi0) && (fixed_pi0 <= 0.0 || fixed_pi0 >= 1.0)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
         << "ERROR: --pi0 " << fixed_pi0 << " is invalid" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(! file_ci.empty() && ! doesFileExist(file_ci)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
         << "ERROR: can't find " << file_ci << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(! file_ci.empty() && skip_ci){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
         << "ERROR: --getci should be given with --ci" << endl << endl;
    help(argv);
    exit(EXIT_FAILURE);
  }
  if(model == "types"){
    configs_tokeep.clear(); // ignore --configs if model is "type"
    for(size_t s = 1; s <= nb_subgroups; ++s)
      configs_tokeep.push_back(toString(s));
  }
}

void run(
  const string & file_pattern,
  const size_t & nb_subgroups,
  const string & model,
  const size_t & dim,
  const size_t & nb_grid_points,
  const string & out_file,
  const string & file_init,
  const size_t & seed,
  const double & thresh,
  const size_t & max_nb_iters,
  const bool & squarem,
  const int & nb_threads,
  const string & file_ci,
  const vector<string> & configs_tokeep,
  const bool & keep_gen_abfs,
  const bool & skip_ci,
  const bool & skip_bf,
  const double & fixed_pi0,
  const int & verbose)
{
  Controller controller(nb_subgroups, model, nb_grid_points, dim, thresh,
                        max_nb_iters, squarem, fixed_pi0, nb_threads,
                        verbose);
  
  controller.load_data(file_pattern, configs_tokeep, keep_gen_abfs);
  
  if(! file_ci.empty()){
    controller.init_params(file_ci);
    controller.compute_posterior();
    controller.estimate_profile_ci();
    controller.save_result(out_file, true);
  }
  else{
    if(! file_init.empty())
      controller.init_params(file_init);
    else
      controller.init_params(seed);
    
    controller.run_EM();
    
    controller.compute_log10_ICL();
    
    if(! skip_bf)
      controller.compute_posterior();
    
    if( ! skip_ci)
      controller.estimate_profile_ci();
    
    controller.save_result(out_file, skip_bf);
  }
}

int main(int argc, char **argv)
{
#ifdef DEBUG
  cout << "DEBUG" << endl;
#endif
  int verbose = 1, nb_threads = 1;
  string file_pattern, model = "configs", out_file, file_init, file_ci;
  size_t nb_subgroups = string::npos, dim = string::npos,
    nb_grid_points = string::npos, seed = string::npos,
    max_nb_iters = string::npos;
  double thresh = 0.05, fixed_pi0 = NaN;
  vector<string> configs_tokeep;
  bool rand_init = false, squarem = false, keep_gen_abfs = false,
    skip_ci = true, skip_bf = true;
  
  parseCmdLine(argc, argv, file_pattern, nb_subgroups, model, dim,
               nb_grid_points, out_file, file_init, rand_init, seed, thresh,
               max_nb_iters, squarem, nb_threads, file_ci, configs_tokeep,
               keep_gen_abfs, skip_ci, skip_bf, fixed_pi0, verbose);
  
  time_t time_start, time_end;
  if(verbose > 0){
    time (&time_start);
    cout << "START " << basename(argv[0])
         << " " << getDateTime(time_start) << endl
         << "version " << VERSION << " compiled " << __DATE__
         << " " << __TIME__ << endl
         << "cmd-line: " << getCmdLine(argc, argv) << endl
         << "cwd: " << getCurrentDirectory() << endl << flush;
  }
  
  run(file_pattern, nb_subgroups, model, dim, nb_grid_points,
      out_file, file_init, seed, thresh, max_nb_iters, squarem, nb_threads,
      file_ci, configs_tokeep, keep_gen_abfs, skip_ci, skip_bf, fixed_pi0,
      verbose);
  
  if(verbose > 0){
    time (&time_end);
    cout << "END " << basename(argv[0])
         << " " << getDateTime(time_end) << endl
         << "elapsed -> " << getElapsedTime(time_start, time_end) << endl
         << "max.mem -> " << getMaxMemUsedByProcess2Str() << endl;
  }
  
  return EXIT_SUCCESS;
}
