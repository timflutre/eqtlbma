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
  std::vector<std::vector<double> > raw_log10_bfs_; // dim1 is config/subgroup; dim2 is grid
  size_t nb_subgroups_;
  size_t grid_size_;
  size_t dim_;
  
  double log10_BF_;
  double post_prob_snp_;
  double post_prob_config_;
  
  // snp specific prior (using information e.g. location to tss/tes)
  double snp_prior;
  double new_snp_prior;
  
  // prior for snp-specific tissue specificity
  double * snp_config_prior; // not used yet
  
  // useful data structures allocated only once
  std::vector<double> grid_tmp_; // used for any computation requiring a vector of length "grid_size_"
  std::vector<double> dim_tmp_; // length "dim_"
  std::vector<double> subgroup_tmp_; // length "nb_subgroups_"
  
  snp_eQTL(const std::string & name,
	   const std::vector<std::vector<double> > & raw_log10_bfs,
	   const size_t & nb_subgroups,
	   const size_t & nb_types);
  
  double compute_log10_BF(const std::vector<double> & grid_wts,
			  const std::vector<double> & config_prior,
			  const bool & keep);
  double compute_log10_BF(const std::vector<double> & grid_wts,
			  const std::vector<double> & type_prior,
			  const std::vector<std::vector<double> > & subgroup_prior,
			  const bool & keep);
  double compute_log10_config_BF(const size_t & config_idx,
				 const std::vector<double> & grid_wts);
  double compute_log10_type_BF(const size_t & type_idx,
			       const std::vector<double> & grid_wts,
			       const std::vector<std::vector<double> > & subgroup_prior);
  
  void em_update_config(const std::vector<double> & grid_wts,
			std::vector<double> & new_global_config_prior);
  void em_update_type(const std::vector<double> & grid_wts,
		      const std::vector<std::vector<double> > & subgroup_prior,
		      std::vector<double> & new_type_prior);
  void em_update_subgroup(const std::vector<double> & grid_wts,
			  const std::vector<std::vector<double> > & subgroup_prior,
			  std::vector<std::vector<double> > & exp_gpkls_snp,
			  std::vector<double> & exp_gpkl_snp);
  void em_update_grid(const std::vector<double> & config_prior,
		      std::vector<double> & new_grid_wts);
  void em_update_grid(const std::vector<double> & type_prior,
		      const std::vector<std::vector<double> > & subgroup_prior,
		      std::vector<double> & new_grid_wts);
  
  void print_result(const std::vector<double> & grid_wts);
};

class gene_eQTL {
public:
  std::string name_;
  size_t nb_subgroups_;
  size_t dim_;
  size_t grid_size_;
  std::vector<snp_eQTL> snpVec;
  double * pi0_;
  double log10_BF_;
  double log10_obs_lik_;
  double post_prob_gene_;
  std::vector<config_prob> post_prob_config_;
  
  // useful data structures allocated only once
  std::vector<double> snp_wts_; // snpVec.size
  std::vector<double> snp_log10_bfs_; // snpVec.size
  std::vector<std::vector<double> > grid_snps_; // grids x snps
  std::vector<std::vector<double> > dim_snps_; // configs or types x snps
  std::vector<std::vector<std::vector<double> > > subgroup_num_snps_; // types x subgroups x snps
  std::vector<std::vector<double> > subgroup_denom_snps_; // types x snps
  
  gene_eQTL(){};
  gene_eQTL(const std::string name, const size_t & nb_subgroups, double * pi0){name_ = name; nb_subgroups_ = nb_subgroups; pi0_ = pi0;};
  void init(void);
  
  double compute_log10_BF(const std::vector<double> & grid_wts,
			  const std::vector<double> & config_prior,
			  const bool & keep);
  double compute_log10_BF(const std::vector<double> & grid_wts,
			  const std::vector<double> & type_prior,
			  const std::vector<std::vector<double> > & subgroup_prior,
			  const bool & keep);
  
  double compute_log10_obs_lik(const double & pi0,
			       const std::vector<double> & grid_wts,
			       const std::vector<double> & config_prior,
			       const bool & keep);
  double compute_log10_obs_lik(const double & pi0,
			       const std::vector<double> & grid_wts,
			       const std::vector<double> & type_prior,
			       const std::vector<std::vector<double> > & subgroup_prior,
			       const bool & keep);
  
  void set_snp_prior();
  void update_snp_prior();
  
  double em_update_pi0(const double & pi0);
  void em_update_config(const std::vector<double> & grid_wts,
			std::vector<double> & new_global_config_prior);
  void em_update_type(const std::vector<double> & grid_wts,
		      const std::vector<std::vector<double> > & subgroup_prior,
		      std::vector<double> & new_type_prior);
  void em_update_subgroup(const std::vector<double> & grid_wts,
			  const std::vector<std::vector<double> > & subgroup_prior,
			  std::vector<std::vector<double> > & exp_gpkls_gene,
			  std::vector<double> & exp_gpkl_gene);
  void em_update_grid(const std::vector<double> & config_prior,
		      std::vector<double> & new_grid_wts);
  void em_update_grid(const std::vector<double> & type_prior,
		      const std::vector<std::vector<double> > & subgroup_prior,
		      std::vector<double> & new_grid_wts);
  void em_update_snp_prior();
  
  void compute_posterior(const double & pi0,
			 const std::vector<double> & grid_wts,
			 const std::vector<double> & config_prior);
  
  void print_result(const std::vector<double> & grid_wts);
};

#endif
