/** \file eqtlbma_hm.cpp
 *
 *  `eqtlbma_hm' implements the EM algorithm to fit the hierarchical model from eQtlBma.
 *  Copyright (C) 2012-2013 Xiaoquan Wen, Timothee Flutre
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
			  vector<string> & lines);
  
public:
  vector<gene_eQTL> geneVec_;
  size_t config_size_;
  vector<string> config_names_;
  double thresh_; // on log-lik, to stop EM
  double fixed_pi0_;
  int nb_threads_; 
  int verbose_;
  
  // parameters to estimate
  double pi0_;
  double new_pi0_;
  vector<double> grid_wts_;
  vector<double> new_grid_wts_;
  vector<double> config_prior_;
  vector<double> new_config_prior_;
  
  // confidence intervals
  double left_pi0_, right_pi0_;
  vector<double> left_configs_, right_configs_;
  vector<double> left_grids_, right_grids_;
  
  // useful data structures for EM updates
  vector<double> gene_wts_;
  vector<double> grid_wts_ones_;
  vector<double> config_wts_ones_;
  double ** grid_mat_;
  double ** config_mat_;
  
  Controller();
  Controller(const size_t & grid_size, const size_t & config_size,
	     const double & thresh, const double & fixed_pi0,
	     const int & nb_threads, const int & verbose);
  ~Controller();
  
  void load_data(const string & file_pattern,
		 const vector<string> & configs_tokeep);
  
  void init_params(const size_t & seed);
  void init_params(const string & init_file);
  
  void randomize_parameter_sp(const size_t & seed);
  
  double compute_log10_lik();
  
  void em_update_pi0();
  void em_update_config();
  void em_update_grid();
  void em_update_snp_prior();
  
  void update_params();
  
  void estimate_profile_ci(const bool & skip_ci);
  
  void run_EM();
  void compute_posterior();
  
  void save_result(const string & out_file, const bool & skip_bf);
};

Controller::Controller()
{
}

Controller::Controller(const size_t & gsize, const size_t & csize,
		       const double & thresh, const double & fixed_pi0,
		       const int & nb_threads, const int & verbose)
{
  grid_wts_.assign(gsize, NaN);
  new_grid_wts_.assign(grid_wts_.size(), NaN);
  grid_wts_ones_.assign(grid_wts_.size(), 1.0);
  if(csize > 0){
    config_prior_.assign(csize, NaN);
    new_config_prior_.assign(config_prior_.size(), NaN);
    config_wts_ones_.assign(config_prior_.size(), 1.0);
  }
  
  thresh_ = thresh;
  fixed_pi0_ = fixed_pi0;
  nb_threads_ = nb_threads;
  verbose_ = verbose;
}

Controller::~Controller()
{
  for(size_t i = 0; i < grid_wts_.size(); ++i)
    delete[] grid_mat_[i];
  delete[] grid_mat_;
  if(config_prior_.size() > 1){
    for(size_t i = 0; i < config_prior_.size(); ++i)
      delete[] config_mat_[i];
    delete[] config_mat_;
  }
}

void Controller::load_data_one_file(
  const string & file,
  const vector<string> & configs_tokeep,
  vector<string> & lines)
{
  // load the whole file at once (supposedly quicker)
  readFile(file, lines);
  
  // check the header
  vector<string> tokens;
  split(lines[0], " \t,", tokens);
  if(tokens[0].compare("gene") != 0 
     || tokens[1].compare("snp") != 0
     || tokens[2].compare("config") != 0){
    cerr << "ERROR: file " << file << " has wrong header line" << endl;
    exit(1);
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
    if(config.find("gen") != string::npos
       || (! configs_tokeep.empty()
	   && find(configs_tokeep.begin(), configs_tokeep.end(), config)
	   == configs_tokeep.end()))
      continue;
    
    if(config_names_.size() < config_prior_.size()) // record config names once (for output)
      config_names_.push_back(config);
    
    if(gene_id.compare(curr_gene) != 0){ // if new gene
      if(! geq.name_.empty()){
	geq.snpVec.push_back(snp_eQTL(curr_snp, raw_log10_bfs));
	geq.init();
	geneVec_.push_back(geq);
      }
      curr_gene = gene_id;
      geq = gene_eQTL(curr_gene, &pi0_);
      raw_log10_bfs.clear();
    }
    
    if(snp_id.compare(curr_snp) != 0){ // if same gene but new snp
      if(! raw_log10_bfs.empty())
	geq.snpVec.push_back(snp_eQTL(curr_snp, raw_log10_bfs));
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
  
  geq.snpVec.push_back(snp_eQTL(curr_snp, raw_log10_bfs)); // last snp
  geq.init();
  geneVec_.push_back(geq); // last gene
}

void Controller::load_data(
  const string & file_pattern,
  const vector<string> & configs_tokeep)
{
  if(verbose_ > 0)
    fprintf(stderr, "load data ...\n");
  clock_t startTime = clock();
  
  if(! configs_tokeep.empty()){
    fprintf(stderr, "configurations to keep: %s", configs_tokeep[0].c_str());
    for(size_t i = 1; i < configs_tokeep.size(); ++i)
      fprintf(stderr, " %s", configs_tokeep[i].c_str());
    fprintf(stderr, "\n");
  }
  
  vector<string> files = glob(file_pattern);
  if(verbose_ > 0)
    cerr << "nb of input files: " << files.size() << endl << flush;
  
  vector<string> lines;
  for(size_t i = 0; i < files.size(); ++i){
    if(verbose_ > 0)
      progressBar("", i+1, files.size());
    load_data_one_file(files[i], configs_tokeep, lines);
    lines.clear();
  }
  if(verbose_ > 0)
    cerr << endl << flush;
  
  // initialization of useful data structures for the EM
  gene_wts_.assign(geneVec_.size(), 1.0);
  grid_mat_ = new double *[grid_wts_.size()];
  for(size_t i = 0; i < grid_wts_.size(); ++i){
    grid_mat_[i] = new double[geneVec_.size()];
    memset(grid_mat_[i], 0, sizeof(double)*geneVec_.size());
  }
  config_mat_ = new double *[config_prior_.size()];
  for(size_t i = 0; i < config_prior_.size(); ++i){
    config_mat_[i] = new double[geneVec_.size()];
    memset(config_mat_[i], 0, sizeof(double)*geneVec_.size());
  }
  
  if(verbose_ > 0){
    size_t nb_pairs = 0;
#pragma omp parallel for num_threads(nb_threads_) reduction(+:nb_pairs)
    for(int g = 0; g < (int)geneVec_.size(); ++g)
      nb_pairs += geneVec_[g].snpVec.size();
    fprintf(stderr, "finish loading %zu genes and %zu gene-snp pairs (%f sec)\n",
	    geneVec_.size(),
	    nb_pairs,
	    getElapsedTime(startTime));
  }
}

void Controller::estimate_profile_ci(const bool & skip_ci)
{
  if(verbose_ > 0)
    cerr << "compute profile-likelihood confidence intervals ..." << flush;
  clock_t startTime = clock();

  fprintf(stderr, "pi0=%7.4e config[0]=%7.4e grid[0]=%7.4e\n", pi0_, config_prior_[0], grid_wts_[0]);
  
  double tick = 0.001;
  double max = compute_log10_lik();
  
  // estimate CI for pi0
  double pi0_mle = pi0_;
  left_pi0_ = pi0_;
  right_pi0_ = pi0_;
  if (! skip_ci) {
    while(left_pi0_ >= 0){
      left_pi0_ -= tick;
      if(left_pi0_ < 0){
	left_pi0_ = 0;
	break;
      }
      pi0_ = left_pi0_;
      if(compute_log10_lik() / log10(exp(1)) < max / log10(exp(1)) - 2.0){  
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
      if(compute_log10_lik() / log10(exp(1)) < max / log10(exp(1)) - 2.0){
	right_pi0_ -= tick;
	break;
      }
    }
  }
  
  pi0_ = pi0_mle;
  
  // estimate CI for config
  vector<double> config_mle(config_prior_);
  left_configs_.assign(config_prior_.size(), NaN);
  right_configs_.assign(config_prior_.size(), NaN);
  for (size_t i = 0 ; i < config_prior_.size(); ++i){
    right_configs_[i] = config_mle[i];
    left_configs_[i] = config_mle[i];
    double cp_mle = config_mle[i];
    double st = 1 - cp_mle;
    if (! skip_ci) {
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
	if(compute_log10_lik() / log10(exp(1)) < max / log10(exp(1)) - 2.0){
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
	if(compute_log10_lik() / log10(exp(1)) < max / log10(exp(1)) - 2.0){
	  right_configs_[i] -= tick;
	  break;
	}
      }
    }
  }
  config_prior_ = config_mle;
  
  // estimate CI for grid
  vector<double> grid_mle(grid_wts_);
  left_grids_.assign(grid_wts_.size(), NaN);
  right_grids_.assign(grid_wts_.size(), NaN);
  for (size_t i = 0; i < grid_wts_.size(); ++i){
    right_grids_[i] = grid_mle[i];
    left_grids_[i] = grid_mle[i];
    double cp_mle = grid_mle[i];
    double st = 1 - cp_mle;
    if (! skip_ci) {
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
	
	if(compute_log10_lik() / log10(exp(1)) < max / log10(exp(1)) - 2.0){
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
	if(compute_log10_lik() / log10(exp(1)) < max / log10(exp(1)) - 2.0){
	  right_grids_[i] -= tick;
	  break;
	  
	}
      }
    }
  }
  grid_wts_ = grid_mle;
  
  if(verbose_ > 0)
    cerr << " (" << getElapsedTime(startTime) << " sec)" << endl;
}

void Controller::init_params(const string & init_file)
{
#pragma omp parallel for num_threads(nb_threads_)
  for(int i = 0; i < (int) geneVec_.size(); ++i)
    geneVec_[i].set_snp_prior();
  
  vector<string> lines;
  readFile(init_file, lines);
  
  vector<string> tokens;
  size_t idx_config = 0, idx_grid = 0;
  for(size_t i = 0; i < lines.size(); ++i){
    split(lines[i], "\t, ", tokens);
    if(tokens[0].find("pi0") != string::npos)
      pi0_ = atof(tokens[1].c_str());
    else if(tokens[0].find("config") != string::npos){
      config_prior_[idx_config] = atof(tokens[1].c_str());
      ++idx_config;
    }
    else if(tokens[0].find("grid") != string::npos){
      grid_wts_[idx_grid] = atof(tokens[1].c_str());
      ++idx_grid;
    }
  }
}

void Controller::randomize_parameter_sp(const size_t & seed)
{
  const gsl_rng_type * T;
  gsl_rng * r;
       
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  gsl_rng_set(r, seed); 
  
  //double x1 = gsl_ran_exponential(r,1.0);
  //double x2 = gsl_ran_exponential(r,1.0);
  //pi0_ = .90 + .1*x1/(x1 + x2);
  //pi0_ = x1/(x1 + x2);
  pi0_ = gsl_rng_uniform (r);
  
  double sum = 0;
  for(size_t i = 0; i < grid_wts_.size(); ++i){
    grid_wts_[i] = gsl_ran_exponential(r,1.0);
    sum += grid_wts_[i];
  }

  for(size_t i = 0; i < grid_wts_.size(); ++i){
    grid_wts_[i] /= sum;
  }

  sum = 0;
  if(config_prior_.size() == 1){
    new_config_prior_[0] = config_prior_[0] = 1.0;
  }else{
    for(size_t i = 0; i < config_prior_.size(); ++i){
      config_prior_[i] = gsl_ran_exponential(r,1.0);
      sum+= config_prior_[i];
    }

    for(size_t i = 0; i < config_prior_.size(); ++i){
      config_prior_[i] /= sum;
    }
  }
  
  gsl_rng_free(r);
}

void Controller::init_params(const size_t & seed)
{
#pragma omp parallel for num_threads(nb_threads_)
  for(int i = 0; i < (int)geneVec_.size(); ++i)
    geneVec_[i].set_snp_prior();
  
  if(seed != string::npos){
    randomize_parameter_sp(seed);
    return;
  }
  
  pi0_ = 0.5;
  
  if(config_prior_.size() == 1){
    config_prior_[0] = 1.0;
    new_config_prior_[0] = 1.0;
  }
  else
    for(size_t i = 0; i < config_prior_.size(); ++i)
      config_prior_[i] = 1.0 / (double) config_prior_.size();
  
  for(size_t i = 0; i < grid_wts_.size(); ++i)
    grid_wts_[i] = 1.0 / (double) grid_wts_.size();
}

double Controller::compute_log10_lik()
{
  double sum = 0;
  
#pragma omp parallel for num_threads(nb_threads_) reduction(+:sum)
  for(int i = 0; i < (int)geneVec_.size(); ++i)
    sum += geneVec_[i].compute_log10_lik(pi0_, grid_wts_, config_prior_);
  
  return sum;
}

void Controller::em_update_snp_prior()
{
#pragma omp parallel for num_threads(nb_threads_)
  for (int i = 0; i < (int)geneVec_.size(); ++i)
     geneVec_[i].em_update_snp_prior();
}

void Controller::em_update_pi0()
{
  if (isNan(fixed_pi0_)){
    double sum = 0;
#pragma omp parallel for num_threads(nb_threads_) reduction(+:sum)
    for(int i = 0; i < (int)geneVec_.size(); ++i){
      // geneVec_[i].update_snp_prior();
      sum += geneVec_[i].em_update_pi0(pi0_);
    }
    new_pi0_ = sum / geneVec_.size();
  }
  else
    new_pi0_ = fixed_pi0_;
}

void Controller::em_update_grid()
{
#pragma omp parallel for num_threads(nb_threads_)
  for(int i = 0; i < (int)geneVec_.size(); ++i){
    vector<double> new_grid_wts_tmp(grid_wts_.size(), NaN);
    geneVec_[i].em_update_grid(config_prior_, new_grid_wts_tmp);
    for(size_t j = 0; j < grid_wts_.size(); ++j)
      grid_mat_[j][i] = new_grid_wts_tmp[j];
  }
  
  for(size_t i = 0; i < grid_wts_.size(); ++i)
    new_grid_wts_[i] = log10_weighted_sum(grid_mat_[i], &(gene_wts_[0]), geneVec_.size())
      + log10(grid_wts_[i]);
  
  double sum = log10_weighted_sum(&(new_grid_wts_[0]), &(grid_wts_ones_[0]), grid_wts_.size());
  
  for(size_t i = 0; i < grid_wts_.size(); ++i)
    new_grid_wts_[i] = pow(10, (new_grid_wts_[i] - sum));
}

void Controller::em_update_config()
{
  if(config_prior_.size() == 1)
    return;
  
#pragma omp parallel for num_threads(nb_threads_)
  for(int i = 0; i < (int)geneVec_.size(); ++i){
    vector<double> new_config_prior_tmp(config_prior_.size(), NaN);
    geneVec_[i].em_update_config(grid_wts_, new_config_prior_tmp);
    for(size_t j = 0; j < config_prior_.size(); ++j)
      config_mat_[j][i] = new_config_prior_tmp[j];
  }
  
  for(size_t i = 0; i < config_prior_.size(); ++i)
    new_config_prior_[i] = log10_weighted_sum(config_mat_[i], &(gene_wts_[0]), geneVec_.size()) + log10(config_prior_[i]);
  
  double sum = log10_weighted_sum(&(new_config_prior_[0]), &(config_wts_ones_[0]), config_prior_.size());
  
  for(size_t i = 0; i < config_prior_.size(); ++i)
    new_config_prior_[i] = pow(10, (new_config_prior_[i] - sum));
}

void Controller::update_params()
{
  pi0_ = new_pi0_;
  if(config_prior_.size() > 1){
    for(size_t i = 0; i < config_prior_.size(); ++i)
      config_prior_[i] = new_config_prior_[i];
  }
  
  for(size_t i = 0; i < grid_wts_.size(); ++i)
    grid_wts_[i] = new_grid_wts_[i];
  
  // for(size_t i = 0; i < geneVec_.size(); ++i){
  //  geneVec_[i].update_snp_prior();
  // }
}

void Controller::run_EM()
{
  if(verbose_ > 0)
    cerr << "run EM algorithm ..." << endl << flush;
  time_t startRawTime, endRawTime;
  time (&startRawTime);
  
  size_t iter = 0;
  double prev_log10_lik = compute_log10_lik(), curr_log10_lik;
  
  if(verbose_){
    fprintf(stderr, "iter %4zu  loglik %.4f", iter, prev_log10_lik);
    fprintf(stderr, "  pi0 %7.4e  config", pi0_);
    for(size_t i = 0; i < config_prior_.size(); ++i)
      fprintf(stderr, " %7.4e", config_prior_[i]);
    fprintf(stderr, "  grid");
    for(size_t i = 0; i < grid_wts_.size(); ++i)
      fprintf(stderr, " %7.4e", grid_wts_[i]);
    fprintf(stderr, "\n");
  }
  
  while(true){
    ++iter;
    em_update_pi0();    
    em_update_config();
    em_update_grid();
    em_update_snp_prior();
    update_params();
    curr_log10_lik = compute_log10_lik();
    if(verbose_){
      fprintf(stderr, "iter %4zu  loglik %.4f", iter, curr_log10_lik);
      fprintf(stderr, "  pi0 %7.4e  config", pi0_);
      for(size_t i = 0; i < config_prior_.size(); ++i)
	fprintf(stderr, " %7.4e", new_config_prior_[i]);
      fprintf(stderr, "  grid");
      for(size_t i = 0; i < grid_wts_.size(); ++i)
	fprintf(stderr, " %7.4e", new_grid_wts_[i]);
      fprintf(stderr, "\n");
    }
    
    if(fabs(curr_log10_lik - prev_log10_lik) < thresh_)
      break;
    
    prev_log10_lik = curr_log10_lik;
  }
  
  // do a last update to get same results as William's initial implementation
  ++iter;
  em_update_pi0();    
  em_update_config();
  em_update_grid();
  em_update_snp_prior();
  update_params();
  curr_log10_lik = compute_log10_lik();
  if(verbose_){
    fprintf(stderr, "iter %4zu  loglik %.4f", iter, curr_log10_lik);
    fprintf(stderr, "  pi0 %7.4e  config", pi0_);
    for(size_t i = 0; i < config_prior_.size(); ++i)
      fprintf(stderr, " %7.4e", new_config_prior_[i]);
    fprintf(stderr, "  grid");
    for(size_t i = 0; i < grid_wts_.size(); ++i)
      fprintf(stderr, " %7.4e", new_grid_wts_[i]);
    fprintf(stderr, "\n");
  }
  
  time (&endRawTime);
  cerr << "EM ran for " << getElapsedTime(startRawTime, endRawTime) << endl << flush;
}

void Controller::compute_posterior()
{
  if(verbose_ > 0)
    cerr << "compute posteriors ..." << flush;
  clock_t startTime = clock();
  
#pragma omp parallel for num_threads(nb_threads_)
  for(int i = 0; i < (int)geneVec_.size(); ++i)
    geneVec_[i].compute_posterior(pi0_, grid_wts_, config_prior_);
  
  if(verbose_ > 0)
    cerr << " (" << getElapsedTime(startTime) << " sec)" << endl;
}

void Controller::save_result(const string & out_file,
			     const bool & skip_bf)
{
  if(verbose_ > 0)
    cerr << "save the results in " << out_file << " ..." << flush;
  clock_t startTime = clock();
  
  gzFile stream;
  openFile(out_file, stream, "wb");
  
  stringstream txt;
  size_t nb_lines = 0;
  
  txt << "#param\tmle\tleft.ci\tright.ci" << endl;
  ++nb_lines;
  gzwriteLine(stream, txt.str(), out_file, nb_lines);
  txt.str("");
  txt << setprecision(4);
  txt << "#pi0\t" << scientific << pi0_ << "\t" << left_pi0_ << "\t" << right_pi0_ << endl;
  ++nb_lines;
  gzwriteLine(stream, txt.str(), out_file, nb_lines);
  for(size_t i = 0; i < config_names_.size(); ++i){
    txt.str("");
    txt << "#config." << config_names_[i].c_str()
	<< "\t" << config_prior_[i]
	<< "\t" << left_configs_[i]
	<< "\t" << right_configs_[i]
	<< endl;
    ++nb_lines;
    gzwriteLine(stream, txt.str(), out_file, nb_lines);
  }
  for(size_t i = 0; i < grid_wts_.size(); ++i){
    txt.str("");
    txt << "#grid." << i+1
	<< "\t" << grid_wts_[i]
	<< "\t" << left_grids_[i]
	<< "\t" << right_grids_[i]
	<< endl;
    ++nb_lines;
    gzwriteLine(stream, txt.str(), out_file, nb_lines);
  }

  if(! skip_bf){
    txt.str("");
    txt << "gene\tgene.posterior.prob\tgene.log10.bf\tsnp\tsnp.log10.bf\t";
    for(size_t i = 0; i < config_names_.size(); ++i)
      txt << "log10.bf." << config_names_[i].c_str() << "\t";
    txt << "\n";
    ++nb_lines;
    gzwriteLine(stream, txt.str(), out_file, nb_lines);
    
    for(size_t i = 0; i < geneVec_.size(); ++i){ // loop over gene
      for(size_t j = 0; j < geneVec_[i].snpVec.size(); ++j){ // loop over snp
	txt.str("");
	txt << setprecision(4)
	    << geneVec_[i].name_
	    << "\t"
	    << geneVec_[i].post_prob_gene
	    << "\t"
	    << geneVec_[i].compute_log10_BF(grid_wts_, config_prior_)
	    << "\t"
	    << geneVec_[i].snpVec[j].name_
	    << "\t"
	    << geneVec_[i].snpVec[j].compute_log10_BF(grid_wts_, config_prior_);
	for(size_t k = 0; k < geneVec_[i].snpVec[j].config_size_; ++k) // loop over config
	  txt << "\t"
	      << log10_weighted_sum(&(geneVec_[i].snpVec[j].raw_log10_bfs_[k][0]),
				    &(grid_wts_[0]), grid_wts_.size());
	txt << "\n";
	++nb_lines;
	gzwriteLine(stream, txt.str(), out_file, nb_lines);
      }
    }
  }
  
  closeFile(out_file, stream);
  
  if(verbose_ > 0)
    cerr << " (" << getElapsedTime(startTime) << " sec)" << endl;
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
       << "      --nconf\tnumber of configurations" << endl
       << "      --ngrid\tnumber of grid points" << endl
       << "      --out\toutput file (gzipped)" << endl
       << "      --init\tfile for initialization" << endl
       << "      --rand\trandom initialization" << endl
       << "      --seed\tseed used with --rand, otherwise use time" << endl
       << "      --thresh\tthreshold to stop the EM (default=0.05)" << endl
       << "      --thread\tnumber of threads (default=1)" << endl
       << "      --configs\tsubset of configurations to keep (e.g. \"1|3|1-3\")" << endl
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
       << "Copyright (C) 2012-2013 Xiaoquan Wen and Timothee Flutre." << endl
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
  size_t & nb_active_configs,
  size_t & nb_grid_points,
  string & out_file,
  string & file_init,
  bool & rand_init,
  size_t & seed,
  double & thresh,
  int & nb_threads,
  string & file_ci,
  vector<string> & configs_tokeep,
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
      {"nconf", required_argument, 0, 0},
      {"ngrid", required_argument, 0, 0},
      {"out", required_argument, 0, 0},
      {"init", required_argument, 0, 0},
      {"rand", required_argument, 0, 0},
      {"seed", required_argument, 0, 0},
      {"thresh", required_argument, 0, 0},
      {"thread", required_argument, 0, 0},
      {"ci", required_argument, 0, 0},
      {"configs", required_argument, 0, 0},
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
      if(strcmp(long_options[option_index].name, "nconf") == 0){
	nb_active_configs = atoi(optarg);
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
    exit(1);
  }
  if(nb_active_configs == string::npos){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --nconf" << endl << endl;
    help(argv);
    exit(1);
  }
  if(nb_grid_points == string::npos){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --ngrid" << endl << endl;
    help(argv);
    exit(1);
  }
  if(out_file.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --out" << endl << endl;
    help(argv);
    exit(1);
  }
  if(! file_init.empty() && ! doesFileExist(file_init)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find " << file_init << endl << endl;
    help(argv);
    exit(1);
  }
  if(rand_init && seed == string::npos)
    seed = getSeed();
  if(thresh <= 0.0){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --thresh " << thresh << " is invalid" << endl << endl;
    help(argv);
    exit(1);
  }
  if(nb_threads <= 0)
    nb_threads = 1;
  if(! file_ci.empty() && ! doesFileExist(file_ci)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find " << file_ci << endl << endl;
    help(argv);
    exit(1);
  }
  if(! isNan(fixed_pi0) && (fixed_pi0 <= 0.0 || fixed_pi0 >= 1.0)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --pi0 " << fixed_pi0 << " is invalid" << endl << endl;
    help(argv);
    exit(1);
  }
  if(! file_ci.empty() && ! doesFileExist(file_ci)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find " << file_ci << endl << endl;
    help(argv);
    exit(1);
  }
  if(! file_ci.empty() && skip_ci){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --getci should be given with --ci" << endl << endl;
    help(argv);
    exit(1);
  }
}

void run(
  const string & file_pattern,
  const size_t & nb_active_configs,
  const size_t & nb_grid_points,
  const string & out_file,
  const string & file_init,
  const size_t & seed,
  const double & thresh,
  const int & nb_threads,
  const string & file_ci,
  const vector<string> & configs_tokeep,
  const bool & skip_ci,
  const bool & skip_bf,
  const double & fixed_pi0,
  const int & verbose)
{
  Controller controller(nb_grid_points, nb_active_configs,
			thresh, fixed_pi0, nb_threads, verbose);
  
  controller.load_data(file_pattern, configs_tokeep);
  
  if(! file_ci.empty()){
    controller.init_params(file_ci);
    controller.estimate_profile_ci(skip_ci);
    controller.save_result(out_file, true);
  }
  else{
    if(! file_init.empty())
      controller.init_params(file_init);
    else
      controller.init_params(seed);
    
    controller.run_EM();
    
    if(! skip_bf)
      controller.compute_posterior();
    
    controller.estimate_profile_ci(skip_ci);
    
    controller.save_result(out_file, skip_bf);
  }
}

int main(int argc, char **argv)
{
  int verbose = 1, nb_threads = 1;
  string file_pattern, out_file, file_init, file_ci;
  size_t nb_active_configs = string::npos, nb_grid_points = string::npos,
    seed = string::npos;
  double thresh = 0.05, fixed_pi0 = NaN;
  vector<string> configs_tokeep;
  bool rand_init = false, skip_ci = true, skip_bf = true;
  
  parseCmdLine(argc, argv, file_pattern, nb_active_configs, nb_grid_points,
	       out_file, file_init, rand_init, seed, thresh, nb_threads,
	       file_ci, configs_tokeep, skip_ci, skip_bf, fixed_pi0, verbose);
  
  time_t time_start, time_end;
  if(verbose > 0){
    time (&time_start);
    cerr << "START " << basename(argv[0])
	 << " " << getDateTime(time_start) << endl
	 << "version " << VERSION << " compiled " << __DATE__
	 << " " << __TIME__ << endl
	 << "cmd-line: " << getCmdLine(argc, argv) << endl
	 << "cwd: " << getCurrentDirectory() << endl << flush;
  }
  
  run(file_pattern, nb_active_configs, nb_grid_points,
      out_file, file_init, seed, thresh, nb_threads,
      file_ci, configs_tokeep, skip_ci, skip_bf, fixed_pi0, verbose);
  
  if(verbose > 0){
    time (&time_end);
    cerr << "END " << basename(argv[0])
	 << " " << getDateTime(time_end) << endl
	 << "elapsed -> " << getElapsedTime(time_start, time_end) << endl
	 << "max.mem -> " << getMaxMemUsedByProcess2Str() << endl;
  }
  
  return EXIT_SUCCESS;
}
