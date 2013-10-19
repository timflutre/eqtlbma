/** \file hm_methods.cpp
 *
 *  `hm_methods' contains the methods of the classes defined in `hm_classes.hpp'.
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

#include <cmath>

using namespace std;

#include "utils/utils_math.hpp"
using namespace utils;

#include "hm_classes.hpp"

bool operator< (const config_prob & a, const config_prob & b)
{
  return a.prob < b.prob;
};

snp_eQTL::snp_eQTL(const string & name,
		   const vector<vector<double> > & raw_log10_bfs,
		   const size_t & nb_subgroups,
		   const size_t & nb_types)
{
  name_ = name;
  raw_log10_bfs_ = raw_log10_bfs;
  
#ifdef DEBUG
  for(size_t k = 0; k < raw_log10_bfs_.size(); ++k){
    for(size_t l = 0; l < raw_log10_bfs_[k].size(); ++l){
      if(isNan(raw_log10_bfs_[k][l])){
	cerr << "ERROR: raw log10(BF)[" << (k+1) << "," << (l+1) << "] of snp " << name_ << " is NaN" << endl;
	exit(1);
      }
      if(isInfinite(raw_log10_bfs_[k][l])){
	cerr << "ERROR: raw log10(BF)[" << (k+1) << "," << (l+1) << "] of snp " << name_ << " is +-Inf" << endl;
	exit(1);
      }
    }
  }
#endif
  
  grid_size_ = raw_log10_bfs_[0].size();
  
  if(nb_types == 0){ // i.e. model == "config"
    dim_ = raw_log10_bfs_.size();
    nb_subgroups_ = (size_t) log2(dim_ + 1);
  }
  else{
    nb_subgroups_ = nb_subgroups;
    if(raw_log10_bfs_.size() != nb_subgroups_){
      cerr << "ERROR: snp " << name_ << " has BFs for "
	   << raw_log10_bfs_.size() << " subgroups" << endl;
      exit(1);
    }
    dim_ = nb_types;
  }
}

double snp_eQTL::compute_log10_BF(
  const vector<double> & grid_wts,
  const vector<double> & config_prior,
  const bool & keep)
{
  double log10_BF;
  
  dim_tmp_.assign(dim_, NaN);
  
  // for each config, average BFs over grid
  for(size_t j = 0; j < dim_; ++j){
    dim_tmp_[j] = log10_weighted_sum(&(raw_log10_bfs_[j][0]),
				     &(grid_wts[0]),
				     grid_size_);
#ifdef DEBUG
    if(isNan(dim_tmp_[j])){
      cerr << "ERROR: log10(BF)[" << (j+1) << "] avg over grid of snp "
	   << name_ << " is NaN" << endl;
      exit(1);
    }
    if(isInfinite(dim_tmp_[j])){
      cerr << "ERROR: log10(BF)[" << (j+1) << "] avg over grid of snp "
	   << name_ << " is +-Inf" << endl;
      exit(1);
    }
#endif
  }
  
  // average over configs
  log10_BF = log10_weighted_sum(&(dim_tmp_[0]),
  				&(config_prior[0]),
  				dim_);
  
#ifdef DEBUG
  if(isNan(log10_BF)){
    cerr << "ERROR: log10(BF) avg over grid and configs"
	 << " of snp " << name_ << " is NaN" << endl;
    exit(1);
  }
  if(isInfinite(log10_BF)){
    cerr << "ERROR: log10(BF) avg over grid and configs"
	 << " of snp " << name_ << " is +-Inf" << endl;
    exit(1);
  }
#endif
  
  if(keep)
    log10_BF_ = log10_BF;
  
  vector<double>().swap(dim_tmp_); // free memory
  
  return log10_BF;
}

double snp_eQTL::compute_log10_BF(
  const vector<double> & grid_wts,
  const vector<double> & type_prior,
  const vector<vector<double> > & subgroup_prior,
  const bool & keep)
{
  double log10_BF, tmp;
  
  grid_tmp_.assign(grid_size_, NaN);
  dim_tmp_.assign(dim_, NaN);
  
  // for each type, average BFs over grid
  for(size_t k = 0; k < dim_; ++k){
    
    for(size_t l = 0; l < grid_size_; ++l){
      grid_tmp_[l] = 0.0;
      
      for(size_t s = 0; s < nb_subgroups_; ++s){
	tmp = log10(pow(10, log10(subgroup_prior[k][s])
			+ raw_log10_bfs_[s][l])
		    + 1 - subgroup_prior[k][s]);
#ifdef DEBUG
	if(isNan(tmp)){
	  cerr << "ERROR: log10(BF) for subgroup " << (s+1)
	       << " in grid " << (l+1) << " and type " << (k+1)
	       << " of snp " << name_ << " is NaN" << endl;
	  exit(1);
	}
	if(isInfinite(tmp)){
	  cerr << "ERROR: log10(BF) for subgroup " << (s+1)
	       << " in grid " << (l+1) << " and type " << (k+1)
	       << " of snp " << name_ << " is +-Inf" << endl;
	  exit(1);
	}
#endif
	
	grid_tmp_[l] += tmp;
      }
    }
    
    dim_tmp_[k] = log10_weighted_sum(&(grid_tmp_[0]),
				     &(grid_wts[0]),
				     grid_size_);
  }
  
  // average over types
  log10_BF = log10_weighted_sum(&(dim_tmp_[0]),
				&(type_prior[0]),
				dim_);
  
#ifdef DEBUG
  if(isNan(log10_BF)){
    cerr << "ERROR: log10(BF) avg over grid and types"
	 << " of snp " << name_ << " is NaN" << endl;
    exit(1);
  }
  if(isInfinite(log10_BF)){
    cerr << "ERROR: log10(BF) avg over grid and types"
	 << " of snp " << name_ << " is +-Inf" << endl;
    exit(1);
  }
#endif
  
  if(keep)
    log10_BF_ = log10_BF;
  
  vector<double>().swap(grid_tmp_);
  vector<double>().swap(dim_tmp_);
  
  return log10_BF;
}

void snp_eQTL::em_update_config(const vector<double> & grid_wts,
				vector<double> & new_config_prior)
{
  for(size_t k = 0; k < dim_; ++k)
    new_config_prior[k] = log10_weighted_sum(&(raw_log10_bfs_[k][0]),
					     &(grid_wts[0]),
					     grid_size_);
}

void snp_eQTL::em_update_type(const vector<double> & grid_wts,
			      const vector<vector<double> > & subgroup_prior,
			      vector<double> & new_type_prior)
{
  grid_tmp_.assign(grid_size_, NaN);
  
  for(size_t k = 0; k < dim_; ++k){
    for(size_t l = 0; l < grid_size_; ++l){
      grid_tmp_[l] = 0.0;
      for(size_t s = 0; s < nb_subgroups_; ++s)
	grid_tmp_[l] += log10(pow(10, log10(subgroup_prior[k][s])
				  + raw_log10_bfs_[s][l])
			      + 1 - subgroup_prior[k][s]);
    }
    new_type_prior[k] = log10_weighted_sum(&(grid_tmp_[0]),
					   &(grid_wts[0]),
					   grid_size_);
  }
  
  vector<double>().swap(grid_tmp_);
}

void snp_eQTL::em_update_subgroup(const vector<double> & grid_wts,
				  const vector<vector<double> > & subgroup_prior,
				  vector<vector<double> > & exp_gpkls_snp,
				  vector<double> & exp_gpkl_snp)
{
  grid_tmp_.assign(grid_size_, NaN);
  
  for(size_t k = 0; k < dim_; ++k){
    
    // compute the contribution of the SNP to the numerator of q_ks's update
    for(size_t s = 0; s < nb_subgroups_; ++s){
      for(size_t l = 0; l < grid_size_; ++l){
	grid_tmp_[l] = raw_log10_bfs_[s][l];
	for(size_t s2 = 0; s2 < nb_subgroups_; ++s2)
	  if(s2 != s)
	    grid_tmp_[l] += log10(pow(10, log10(subgroup_prior[k][s2])
				      + raw_log10_bfs_[s2][l])
				  + 1 - subgroup_prior[k][s2]);
      }
      exp_gpkls_snp[k][s] = log10_weighted_sum(&(grid_tmp_[0]),
					       &(grid_wts[0]),
					       grid_size_);
    }
    
    // compute the contribution of the SNP to the denominator of q_ks's update
    for(size_t l = 0; l < grid_size_; ++l){
      grid_tmp_[l] = 0.0;
      for(size_t s = 0; s < nb_subgroups_; ++s)
	grid_tmp_[l] += log10(pow(10, log10(subgroup_prior[k][s])
				  + raw_log10_bfs_[s][l])
			      + 1 - subgroup_prior[k][s]);
    }
    exp_gpkl_snp[k] = log10_weighted_sum(&(grid_tmp_[0]),
					 &(grid_wts[0]),
					 grid_size_);
  }
  
  vector<double>().swap(grid_tmp_);
}

void snp_eQTL::em_update_grid(const vector<double> & config_prior,
			      vector<double> & new_grid_wts)
{
  dim_tmp_.assign(dim_, NaN);
  
  for(size_t l = 0; l < grid_size_; ++l){
    for(size_t k = 0; k < dim_; ++k)
      dim_tmp_[k] = raw_log10_bfs_[k][l];
    new_grid_wts[l] = log10_weighted_sum(&(dim_tmp_[0]),
					 &(config_prior[0]),
					 dim_);
  }
  
  vector<double>().swap(dim_tmp_);
}

void snp_eQTL::em_update_grid(const vector<double> & type_prior,
			      const vector<vector<double> > & subgroup_prior,
			      vector<double> & new_grid_wts)
{
  dim_tmp_.assign(dim_, NaN);
  
  for(size_t l = 0; l < grid_size_; ++l){
    for(size_t k = 0; k < dim_; ++k){
      dim_tmp_[k] = 0.0;
      for(size_t s = 0; s < nb_subgroups_; ++s)
	dim_tmp_[k] += log10(pow(10, log10(subgroup_prior[k][s])
				 + raw_log10_bfs_[s][l])
			     + 1 - subgroup_prior[k][s]);
    }
    new_grid_wts[l] = log10_weighted_sum(&(dim_tmp_[0]),
					 &(type_prior[0]),
					 dim_);
  }
  
  vector<double>().swap(dim_tmp_);
}

double snp_eQTL::compute_log10_aug_BF(const double & pi0,
				      const vector<double> & grid_wts,
				      const vector<double> & type_prior,
				      const vector<vector<double> > & subgroup_prior,
				      const double & log10_obs_lik_gene)
{
  double log10_aug_BF = 0.0;
  
  for(size_t l = 0; l < grid_size_; ++l)
    log10_aug_BF += log10(grid_wts[l]);
  
  for(size_t k = 0; k < dim_; ++k)
    log10_aug_BF += log10(type_prior[k]);
  
  for(size_t l = 0; l < grid_size_; ++l)
    for(size_t k = 0; k < dim_; ++k)
      for(size_t s = 0; s < nb_subgroups_; ++s){
	
	double log10_p_gpks = log10(subgroup_prior[k][s])
	  + raw_log10_bfs_[s][l]
	  - log10(pow(10, log10(subgroup_prior[k][s])
		      + raw_log10_bfs_[s][l])
		  + (1-subgroup_prior[k][s]));
	
	if(pow(10, log10_p_gpks) >= 0.5)
	  log10_aug_BF += log10(subgroup_prior[k][s])
	    + raw_log10_bfs_[s][l];
	else
	  log10_aug_BF += log10(1 - subgroup_prior[k][s]);
      }
  
  return log10_aug_BF;
}

double snp_eQTL::compute_log10_config_BF(const size_t & config_idx,
					 const vector<double> & grid_wts)
{
  return log10_weighted_sum(&(raw_log10_bfs_[config_idx][0]),
			    &(grid_wts[0]),
			    grid_size_);
}

double snp_eQTL::compute_log10_type_BF(const size_t & type_idx,
				       const vector<double> & grid_wts,
				       const vector<vector<double> > & subgroup_prior)
{
  double log10_BF;
  
  grid_tmp_.assign(grid_size_, NaN);
  
  for(size_t l = 0; l < grid_size_; ++l){
    grid_tmp_[l] = 0.0;
    for(size_t s = 0; s < nb_subgroups_; ++s)
      grid_tmp_[l] += log10(pow(10, log10(subgroup_prior[type_idx][s])
				+ raw_log10_bfs_[s][l])
			    + 1 - subgroup_prior[type_idx][s]);
  }
  
  log10_BF = log10_weighted_sum(&(grid_tmp_[0]),
				&(grid_wts[0]),
				grid_size_);
  
  vector<double>().swap(grid_tmp_);
  
  return log10_BF;
}

gene_eQTL::gene_eQTL(const std::string name,
		     const size_t & nb_subgroups,
		     const size_t & dim,
		     const size_t & grid_size,
		     double * pi0)
{
  name_ = name;
  nb_subgroups_ = nb_subgroups;
  dim_ = dim;
  grid_size_ = grid_size;
  pi0_ = pi0;
}

void gene_eQTL::set_snp_prior()
{
  for(size_t p = 0; p < snps_.size(); ++p)
    snps_[p].snp_prior = 1.0 / snps_.size(); // for now
}

double gene_eQTL::compute_log10_BF(const vector<double> & grid_wts,
				   const vector<double> & config_prior,
				   const bool & keep)
{
  double log10_BF;
  
  snp_wts_.assign(snps_.size(), NaN);
  snp_log10_bfs_.assign(snps_.size(), NaN);
  
  for(size_t p = 0; p < snps_.size(); ++p){
    snp_wts_[p] = snps_[p].snp_prior;
    snp_log10_bfs_[p] = snps_[p].compute_log10_BF(grid_wts,
						   config_prior,
						   keep);
  }
  
  log10_BF = log10_weighted_sum(&(snp_log10_bfs_[0]),
				&(snp_wts_[0]),
				snps_.size());
#ifdef DEBUG
  if(isNan(log10_BF)){
    cerr << "ERROR: log10(BF) avg over SNPs of gene " << name_ << " is NaN" << endl;
    exit(1);
  }
  if(isInfinite(log10_BF)){
    cerr << "ERROR: log10(BF) avg over SNPs of gene " << name_ << " is +-Inf" << endl;
    exit(1);
  }
#endif
  
  if(keep)
    log10_BF_ = log10_BF;
  
  vector<double>().swap(snp_wts_);
  vector<double>().swap(snp_log10_bfs_);
  
  return log10_BF;
}

double gene_eQTL::compute_log10_BF(
  const vector<double> & grid_wts,
  const vector<double> & type_prior,
  const vector<vector<double> > & subgroup_prior,
  const bool & keep)
{
  double log10_BF;
  
  snp_wts_.assign(snps_.size(), NaN);
  snp_log10_bfs_.assign(snps_.size(), NaN);
  
  for(size_t p = 0; p < snps_.size(); ++p){
    snp_wts_[p] = snps_[p].snp_prior;
    snp_log10_bfs_[p] = snps_[p].compute_log10_BF(grid_wts,
						   type_prior,
						   subgroup_prior,
						   keep);
  }
  
  log10_BF = log10_weighted_sum(&(snp_log10_bfs_[0]),
				&(snp_wts_[0]),
				snps_.size());
#ifdef DEBUG
  if(isNan(log10_BF)){
    cerr << "ERROR: log10(BF) avg over SNPs of gene " << name_ << " is NaN" << endl;
    exit(1);
  }
  if(isInfinite(log10_BF)){
    cerr << "ERROR: log10(BF) avg over SNPs of gene " << name_ << " is +-Inf" << endl;
    exit(1);
  }
#endif
  
  if(keep)
    log10_BF_ = log10_BF;
  
  vector<double>().swap(snp_wts_);
  vector<double>().swap(snp_log10_bfs_);
  
  return log10_BF;
}

// log10(obslik) = log10(p(Y | X, Theta))
//               = log10(pi0 x p(Y|X,Theta,z=0) + (1-pi0) x p(Y|X,Theta,z=1))
//               = log10(pi0 x 10^0 + (1-pi0) x 10^log10(BF))
double gene_eQTL::compute_log10_obs_lik(
  const double & pi0,
  const vector<double> & grid_wts,
  const vector<double> & config_prior,
  const bool & keep)
{
  double log10_obs_lik;
  
  double vec[2];
  double wts[2];
  wts[0] = pi0;
  wts[1] = 1 - pi0;
  vec[0] = 0;
  vec[1] = compute_log10_BF(grid_wts, config_prior, keep);
  
  log10_obs_lik = log10_weighted_sum(vec, wts, 2);
  
  if(keep)
    log10_obs_lik_ = log10_obs_lik;
  
  return log10_obs_lik;
}

double gene_eQTL::compute_log10_obs_lik(
  const double & pi0,
  const vector<double> & grid_wts,
  const vector<double> & type_prior,
  const vector<vector<double> > & subgroup_prior,
  const bool & keep)
{
  double log10_obs_lik;
  
  double vec[2];
  double wts[2];
  wts[0] = pi0;
  wts[1] = 1 - pi0;
  vec[0] = 0;
  vec[1] = compute_log10_BF(grid_wts, type_prior, subgroup_prior, keep);
  
  log10_obs_lik = log10_weighted_sum(vec, wts, 2);
  
  if(keep)
    log10_obs_lik_ = log10_obs_lik;
  
  return log10_obs_lik;
}

double gene_eQTL::em_update_pi0(const double & pi0)
{
  return pow(10, log10(pi0) - log10_obs_lik_);
}

// last to do
void gene_eQTL::em_update_snp_prior()
{
  snp_wts_.assign(snps_.size(), NaN);
  snp_log10_bfs_.assign(snps_.size(), NaN);
  
  for(size_t p = 0; p < snps_.size(); ++p){
    snp_wts_[p] = snps_[p].snp_prior;
    snp_log10_bfs_[p] = snps_[p].log10_BF_;
  }
  
  double factor = log10_weighted_sum(&(snp_log10_bfs_[0]),
				     &(snp_wts_[0]),
				     snps_.size());         
  
  for(size_t p = 0; p < snps_.size(); ++p)
    snps_[p].new_snp_prior = pow(10, snps_[p].log10_BF_ - factor) * snps_[p].snp_prior;
  
  vector<double>().swap(snp_wts_);
  vector<double>().swap(snp_log10_bfs_);
}

void gene_eQTL::update_snp_prior()
{
  double sum = 0;
  for(size_t p = 0; p < snps_.size(); ++p){
    if(snps_[p].new_snp_prior < 0.001)
      snps_[p].new_snp_prior = 0.001;
    sum += snps_[p].new_snp_prior;
  }
  
 for(size_t p = 0; p < snps_.size(); ++p)
    snps_[p].snp_prior = snps_[p].new_snp_prior / sum;
}

void gene_eQTL::em_update_config(const vector<double> & grid_wts,
				 vector<double> & new_global_config_prior)
{
  vector<double> new_global_config_prior_tmp(dim_, NaN);
  
  snp_wts_.assign(snps_.size(), NaN);
  dim_snps_.assign(dim_,
		   vector<double>(snps_.size(), NaN));
  
  for(size_t p = 0; p < snps_.size(); ++p){
    snp_wts_[p] = snps_[p].snp_prior;
    snps_[p].em_update_config(grid_wts, new_global_config_prior_tmp);
    for(size_t k = 0; k < dim_; ++k)
      dim_snps_[k][p] = new_global_config_prior_tmp[k];
  }
  
  for(size_t k = 0; k < dim_; ++k)
    new_global_config_prior[k] = log10_weighted_sum(&(dim_snps_[k][0]),
						    &(snp_wts_[0]),
						    snps_.size())
      - log10_obs_lik_;
  
  vector<double>().swap(snp_wts_);
  vector<vector<double> >().swap(dim_snps_);
}

void gene_eQTL::em_update_type(const vector<double> & grid_wts,
			       const vector<vector<double> > & subgroup_prior,
			       vector<double> & new_type_prior)
{
  vector<double> new_type_prior_tmp(dim_, NaN);
  
  snp_wts_.assign(snps_.size(), NaN);
  dim_snps_.assign(dim_,
		   vector<double>(snps_.size(), NaN));
  
  for(size_t p = 0; p < snps_.size(); ++p){
    snp_wts_[p] = snps_[p].snp_prior;
    snps_[p].em_update_type(grid_wts, subgroup_prior, new_type_prior_tmp);
    for(size_t k = 0; k < dim_; ++k)
      dim_snps_[k][p] = new_type_prior_tmp[k];
  }
  
  for(size_t k = 0; k < dim_; ++k)
    new_type_prior[k] = log10_weighted_sum(&(dim_snps_[k][0]),
					   &(snp_wts_[0]),
					   snps_.size())
      - log10_obs_lik_;
  
  vector<double>().swap(snp_wts_); // free memory
  vector<vector<double> >().swap(dim_snps_);
}

void gene_eQTL::em_update_subgroup(const vector<double> & grid_wts,
				   const vector<vector<double> > & subgroup_prior,
				   vector<vector<double> > & exp_gpkls_gene,
				   vector<double> & exp_gpkl_gene)
{
  snp_wts_.assign(snps_.size(), NaN);
  subgroup_num_snps_.assign(dim_,
			    vector<vector<double> >(nb_subgroups_,
						    vector<double>(snps_.size(), NaN)));
  subgroup_denom_snps_.assign(dim_,
			      vector<double>(snps_.size(), NaN));
  
  // compute the contribution of each SNP to each "subgroup per type" weight
  for(size_t p = 0; p < snps_.size(); ++p){
    vector<vector<double> > exp_gpkls_snp(dim_, vector<double>(nb_subgroups_, NaN));
    vector<double> exp_gpkl_snp(dim_, NaN);
    snp_wts_[p] = snps_[p].snp_prior;
    snps_[p].em_update_subgroup(grid_wts, subgroup_prior, exp_gpkls_snp, exp_gpkl_snp);
    for(size_t k = 0; k < dim_; ++k){
      subgroup_denom_snps_[k][p] = exp_gpkl_snp[k];
      for(size_t s = 0; s < nb_subgroups_; ++s)
	subgroup_num_snps_[k][s][p] = exp_gpkls_snp[k][s];
    }
  }
  
  // average the contribution of all SNPs per "subgroup per type"
  for(size_t k = 0; k < dim_; ++k)
    for(size_t s = 0; s < nb_subgroups_; ++s)
      exp_gpkls_gene[k][s] = log10_weighted_sum(&(subgroup_num_snps_[k][s][0]),
						&(snp_wts_[0]),
						snps_.size())
	- log10_obs_lik_;
  
  // compute the normalization constant per type (same for all subgroups within a type)
  for(size_t k = 0; k < dim_; ++k)
    exp_gpkl_gene[k] = log10_weighted_sum(&(subgroup_denom_snps_[k][0]),
					  &(snp_wts_[0]),
					  snps_.size())
      - log10_obs_lik_;
  
  vector<double>().swap(snp_wts_);
  vector<vector<vector<double> > >().swap(subgroup_num_snps_);
  vector<vector<double> >().swap(subgroup_denom_snps_);
}

void gene_eQTL::em_update_grid(const vector<double> & config_prior,
			       vector<double> & new_grid_wts)
{
  snp_wts_.assign(snps_.size(), NaN);
  grid_snps_.assign(grid_size_,
   		    vector<double>(snps_.size(), NaN));
  
  // compute the contribution of each SNP to each grid weight
  vector<double> new_grid_wts_tmp(grid_size_, NaN);
  for(size_t p = 0; p < snps_.size(); ++p){
    snp_wts_[p] = snps_[p].snp_prior;
    snps_[p].em_update_grid(config_prior, new_grid_wts_tmp);
    for(size_t l = 0; l < grid_size_; ++l)
      grid_snps_[l][p] = new_grid_wts_tmp[l];
  }
  
  // compute the normalization constant (per grid weight)
  for(size_t l = 0; l < grid_size_; ++l)
    new_grid_wts[l] = log10_weighted_sum(&(grid_snps_[l][0]),
					 &(snp_wts_[0]),
					 snps_.size())
      - log10_obs_lik_;
  
  vector<double>().swap(snp_wts_);
  vector<vector<double> >().swap(grid_snps_);
}

void gene_eQTL::em_update_grid(const vector<double> & type_prior,
			       const vector<vector<double> > & subgroup_prior,
			       vector<double> & new_grid_wts)
{
  snp_wts_.assign(snps_.size(), NaN);
  grid_snps_.assign(grid_size_,
   		    vector<double>(snps_.size(), NaN));
  
  // compute the contribution of each SNP to each grid weight
  vector<double> new_grid_wts_tmp(grid_size_, NaN);
  for(size_t i = 0; i < snps_.size(); ++i){
    snp_wts_[i] = snps_[i].snp_prior;
    snps_[i].em_update_grid(type_prior, subgroup_prior, new_grid_wts_tmp);
    for(size_t j = 0; j < grid_size_; ++j)
      grid_snps_[j][i] = new_grid_wts_tmp[j];
  }
  
  // compute the normalization constant (per grid weight)
  for(size_t i = 0; i < grid_size_; ++i)
    new_grid_wts[i] = log10_weighted_sum(&(grid_snps_[i][0]),
					 &(snp_wts_[0]),
					 snps_.size())
      - log10_obs_lik_;
  
  vector<double>().swap(snp_wts_);
  vector<vector<double> >().swap(grid_snps_);
}

double gene_eQTL::compute_log10_aug_BF(const double & pi0,
				       const vector<double> & grid_wts,
				       const vector<double> & type_prior,
				       const vector<vector<double> > & subgroup_prior)
{
  double log10_aug_BF = 0.0;
  
  for(size_t p = 0; p < snps_.size(); ++p)
    log10_aug_BF += log10(snps_[p].snp_prior)
      + snps_[p].compute_log10_aug_BF(pi0, grid_wts, type_prior,
				       subgroup_prior, log10_obs_lik_);
  
  return log10_aug_BF;
}

// ignore p(Y_g | X_g, z_g=0)
double gene_eQTL::compute_log10_aug_lik(
  const double & pi0,
  const vector<double> & grid_wts,
  const vector<double> & type_prior,
  const vector<vector<double> > & subgroup_prior)
{
  double log10_aug_lik;
  
  if(pow(10, log10(1 - pi0) + log10_BF_ - log10_obs_lik_) >= 0.5)
    log10_aug_lik = log10(1 - pi0)
      + compute_log10_aug_BF(pi0, grid_wts, type_prior, subgroup_prior);
  else
    log10_aug_lik = log10(pi0);
  
  return log10_aug_lik;
}

void gene_eQTL::compute_posterior(const double & pi0,
				  const vector<double> & grid_wts,
				  const vector<double> & global_config_prior)
{
  // gene level posterior
  post_prob_gene_ = pow(10, log10(1.0-pi0) + log10_BF_ - log10_obs_lik_);
  if(post_prob_gene_ > 1.0)
    post_prob_gene_ = 1.0;
  
  // posterior for snp eQTL P(S = i, Z =1 | Y)
  for(size_t i = 0; i < snps_.size(); ++i){
    snps_[i].post_prob_snp_ = pow(10, log10(1-pi0) + log10(snps_[i].snp_prior) + snps_[i].log10_BF_ - log10_obs_lik_);
    if(snps_[i].post_prob_snp_ > 1.0)
      snps_[i].post_prob_snp_ = 1.0;
  }
  
  // posterior for configuration
  for(size_t i = 0; i < global_config_prior.size(); ++i){
    double prob = 0;
    for(size_t j = 0; j < snps_.size(); ++j){
      double bfc = snps_[j].compute_log10_config_BF(i, grid_wts);
      double sprior = snps_[j].snp_prior * (1-pi0) * global_config_prior[i];
      prob += sprior * pow(10, bfc - log10_obs_lik_);
    }
    
    if(prob > 1.0)
      prob = 1.0;
    post_prob_config_.push_back(config_prob(i+1, prob));
  }
}
