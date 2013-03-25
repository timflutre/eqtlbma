/** \file hm_classes.h
 *
 *  `hm_classes' contains the classes required by the `hm' program.
 *  Copyright (C) 2012-2013 Xiaoquan Wen
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

#ifndef __CLASSDEF_H_
#define __CLASSDEF_H_

#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


using namespace std;


class config_prob {
 public:
  size_t id;
  double prob;
  config_prob(size_t cid, double post_prob) { id = cid; prob = post_prob; };
};


class snp_eQTL {

 public:
  string snp;
  // possible snp information 
  
  // priors
  // snp specific prior (using information e.g. location to tss/tes)
  double snp_prior;
  double new_snp_prior;
  
  // prior for each tissue specific configuration -- global (same for all gene/snps) 
  size_t config_size;
  double* global_config_prior;
  
  // prior for snp-specific tissue specificity
  double* snp_config_prior; // not used yet

  // shared quantities - grid specification -- global (same same for all gene/snps)
  double* grid_wts; 
 
  double log10_BF;
  
  double post_prob_snp;
  double post_prob_config;
  
    
  // config x wts grid matrix
  double **gm;

  size_t grid_size;
  
  snp_eQTL(string name, const vector<vector<double> > & gm_vec,double *gcp, double *gw);
  

  double compute_log10_BF();
  double compute_log10_config_BF(size_t config);
  
  double *em_update_config();
  double *em_update_grid();
  
  void print_result(void);
  
 
};



class gene_eQTL {
  
 public:
  
  // possible gene information/annotation
  string gene;
  
  vector<snp_eQTL> snpVec;
  double* Pi0;
  double log10_lik; // log10 of likelihood
  
  double log10_BF;
  

  double post_prob_gene;
  vector<config_prob> post_prob_config;


  
  gene_eQTL(string name, double *p0){gene = name; Pi0 = p0;};
  gene_eQTL(){};
  
  double *em_update_config(size_t config_size);
  double *em_update_grid(size_t grid_size);
  double  em_update_pi0();
  void    em_update_snp_prior();
  
 
  void set_snp_prior();
  
  double compute_log10_BF();
  double compute_log10_lik();
  
  void update_snp_prior();
  
  void compute_posterior(double *config_prior, size_t config_size);
  
  void print_result(void);

};


class eQTL_controller {

 public:
  // storage
  vector<gene_eQTL> geqVec;
  
  double fixed_pi0_;
      
  // parameters need to be estimated
  double * pi0; // non-eqtl prob
  double new_pi0;

  double *grid_wts;
  double *new_grid_wts;

  double *config_prior;
  double *new_config_prior; 

  double *param_est;
  size_t     param_size;

  size_t types;
  size_t grid_size;
  size_t config_size;
 
  int output_option;

  vector<string> type_vec;

  void load_data(char *filename, size_t csize, size_t gsize);
  
  
  void init_params(int option=0);
  void init_params(char * init_file);
  void fix_pi0(const double & fixed_pi0) { fixed_pi0_ = fixed_pi0; };
  
  void randomize_parameter_sp();
  
  double compute_log10_lik();

  void print_estimate();
  void print_estimate(size_t index);

  void em_update_pi0();
  void em_update_config();
  void em_update_grid();
  void em_update_snp_prior();

  void update_params();
  
  void update_param_est(size_t index, double val);
  
  
  void estimate_profile_ci(const bool & skip_ci);
  
  void run_EM(double thresh);
  void compute_posterior();
  void print_result();
};


bool operator< (const config_prob& a, const config_prob& b);


double log10_weighted_sum(double* vec, double *wts, size_t size);

  
#endif
