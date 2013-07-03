/** \file hm_classes.hpp
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

class config_prob {
public:
  size_t id;
  double prob;
  config_prob(size_t cid, double post_prob) { id = cid; prob = post_prob; };
};

bool operator< (const config_prob & a, const config_prob & b);

class snp_eQTL {
public:
  std::string name_;
  std::vector<std::vector<double> > raw_log10_bfs_; // config x wts grid matrix
  size_t config_size_;
  size_t grid_size_;
  
  double log10_BF;
  double post_prob_snp;
  double post_prob_config;
  
  // snp specific prior (using information e.g. location to tss/tes)
  double snp_prior;
  double new_snp_prior;
  
  // prior for snp-specific tissue specificity
  double * snp_config_prior; // not used yet
  
  // useful data structures allocated only once
  std::vector<double> config_tmp_; // used for any computation requiring a vector of length "config_size_"
  
  snp_eQTL(const std::string & name,
	   const std::vector<std::vector<double> > & raw_log10_bfs);
  
  double compute_log10_BF(const std::vector<double> & grid_wts,
			  const std::vector<double> & global_config_prior);
  double compute_log10_config_BF(size_t config,
				 const std::vector<double> & grid_wts);
  
  void em_update_grid(const std::vector<double> & global_config_prior,
		      std::vector<double> & new_grid_wts);
  void em_update_config(const std::vector<double> & grid_wts,
			    std::vector<double> & new_global_config_prior);
  
  void print_result(const std::vector<double> & grid_wts);
};

class gene_eQTL {
public:
  std::string name_;
  size_t config_size_;
  size_t grid_size_;
  std::vector<snp_eQTL> snpVec;
  double * pi0_;
  double log10_BF;
  double log10_lik;
  double post_prob_gene;
  std::vector<config_prob> post_prob_config;
  
  // useful data structures allocated only once
  std::vector<double> snp_wts_; // snpVec.size
  std::vector<double> snp_log10_bfs_; // snpVec.size
  std::vector<std::vector<double> > grid_mat_; // grid_size x snpVec.size
  std::vector<std::vector<double> > config_mat_; // config_size x snpVec.size
  
  gene_eQTL(){};
  gene_eQTL(std::string name, double * pi0){name_ = name; pi0_ = pi0;};
  void init(void);
  
  void set_snp_prior();
  void update_snp_prior();
  
  void em_update_grid(const std::vector<double> & global_config_prior,
		      vector<double> & new_grid_wts);
  void em_update_config(const std::vector<double> & grid_wts,
			vector<double> & new_global_config_prior);
  double em_update_pi0(const double & pi0);
  void em_update_snp_prior();
  
  double compute_log10_BF(const std::vector<double> & grid_wts,
			  const std::vector<double> & global_config_prior);
  double compute_log10_lik(const double & pi0,
			   const std::vector<double> & grid_wts,
			   const std::vector<double> & global_config_prior);
  void compute_posterior(const double & pi0,
			 const std::vector<double> & grid_wts,
			 const std::vector<double> & config_prior);
  
  void print_result(const std::vector<double> & grid_wts);
};

#endif
