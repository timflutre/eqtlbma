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
		     const vector<vector<double> > & raw_log10_bfs)
{
  name_ = name;
  raw_log10_bfs_ = raw_log10_bfs;
  
  config_size_ = raw_log10_bfs_.size();
  grid_size_ = raw_log10_bfs_[0].size();
  
  config_tmp_.assign(config_size_, NaN);
}

void snp_eQTL::em_update_config(const vector<double> & grid_wts,
				  vector<double> & new_global_config_prior)
{
  for(size_t i = 0; i < config_size_; ++i)
    new_global_config_prior[i] = log10_weighted_sum(&(raw_log10_bfs_[i][0]),
						    &(grid_wts[0]),
						      grid_size_);
}

void snp_eQTL::em_update_grid(const vector<double> & global_config_prior,
				std::vector<double> & new_grid_wts)
{
  for(size_t i = 0; i < grid_size_; ++i){
    for(size_t j = 0; j < config_size_; ++j)
      config_tmp_[j] = raw_log10_bfs_[j][i];
     new_grid_wts[i] = log10_weighted_sum(&(config_tmp_[0]),
					    &(global_config_prior[0]),
					    config_size_);
  }
}

double snp_eQTL::compute_log10_BF(const vector<double> & grid_wts,
				    const vector<double> & global_config_prior)
{
  for(size_t i = 0; i < config_size_; ++i)
    config_tmp_[i] = log10_weighted_sum(&(raw_log10_bfs_[i][0]), &(grid_wts[0]), grid_size_);
  
  double wsum = log10_weighted_sum(&(config_tmp_[0]), &(global_config_prior[0]), config_size_);
  
  log10_BF = wsum; // save for EM update
  
  return wsum;
}

double snp_eQTL::compute_log10_config_BF(size_t config, const vector<double> & grid_wts)
{
return log10_weighted_sum(&(raw_log10_bfs_[config][0]), &(grid_wts[0]), grid_size_);
}

void snp_eQTL::print_result(const vector<double> & grid_wts)
{
  printf("%12s\t%9.4f\t\t", name_.c_str(), log10_BF);

  for(size_t i = 0; i < config_size_; ++i)
    config_tmp_[i] = log10_weighted_sum(&(raw_log10_bfs_[i][0]), &(grid_wts[0]), grid_size_); 
  
  //double *cprob = new double[config_size_];
  for(size_t i = 0; i < config_size_; ++i){
    //cprob[i] = pow(10, (log10(global_config_prior[i])+rst[i]-log10_BF)); 
    //if(cprob[i]>1)
    //  cprob[i] = 1.0;   
    printf("%7.4f\t", config_tmp_[i]);
  }
}

void gene_eQTL::init()
{
  grid_size_ = snpVec[0].grid_size_;
  config_size_ = snpVec[0].config_size_;
  snp_wts_.assign(snpVec.size(), NaN);
  snp_log10_bfs_.assign(snpVec.size(), NaN);
  grid_mat_.assign(grid_size_,
		     vector<double>(snpVec.size(), NaN));
  config_mat_.assign(config_size_,
		     vector<double>(snpVec.size(), NaN));
}

void gene_eQTL::set_snp_prior()
{
  for(size_t i = 0; i < snpVec.size(); ++i)
    snpVec[i].snp_prior = 1.0 / snpVec.size(); // for now
}

double gene_eQTL::compute_log10_BF(const vector<double> & grid_wts,
				     const vector<double> & global_config_prior)
{
  for(size_t i = 0; i < snpVec.size(); ++i){
    snp_wts_[i] = snpVec[i].snp_prior;
    snp_log10_bfs_[i] = snpVec[i].compute_log10_BF(grid_wts, global_config_prior);
  }
  
  double rst = log10_weighted_sum(&(snp_log10_bfs_[0]), &(snp_wts_[0]), snpVec.size());
  
  log10_BF = rst;

  return rst;
}

double gene_eQTL::compute_log10_lik(const double & pi0,
				      const vector<double> & grid_wts,
				      const vector<double> & global_config_prior)
{
  double vec[2];
  double wts[2];
  wts[0] = pi0;
  wts[1] = 1 - pi0;
  vec[0] = 0;
  vec[1] = compute_log10_BF(grid_wts, global_config_prior);

  log10_lik = log10_weighted_sum(vec, wts, 2);
  
  return log10_lik;
}

double gene_eQTL::em_update_pi0(const double & pi0)
{
  return pow(10, log10(pi0) - log10_lik);
}

// last to do
void gene_eQTL::em_update_snp_prior()
{
  for(size_t i = 0; i < snpVec.size(); ++i){ 
    snp_wts_[i] = snpVec[i].snp_prior;
    snp_log10_bfs_[i] = snpVec[i].log10_BF;
  }
  
  double factor = log10_weighted_sum(&(snp_log10_bfs_[0]), &(snp_wts_[0]), snpVec.size());         
  
  for(size_t i=0;i<snpVec.size();i++)
    snpVec[i].new_snp_prior = pow(10, snpVec[i].log10_BF - factor) * snpVec[i].snp_prior;
}

void gene_eQTL::update_snp_prior()
{
  double sum = 0;
  for(size_t i = 0; i < snpVec.size(); ++i){
    if(snpVec[i].new_snp_prior<0.001)
      snpVec[i].new_snp_prior = 0.001;
    sum += snpVec[i].new_snp_prior;
  }

 for(size_t i = 0; i < snpVec.size(); ++i)
    snpVec[i].snp_prior = snpVec[i].new_snp_prior / sum;
}

void gene_eQTL::em_update_grid(const vector<double> & global_config_prior,
			       vector<double> & new_grid_wts)
{
  vector<double> new_grid_wts_tmp(grid_size_, NaN);
  for(size_t i = 0; i < snpVec.size(); ++i){
    snp_wts_[i] = snpVec[i].snp_prior;
    snpVec[i].em_update_grid(global_config_prior, new_grid_wts_tmp);
    for(size_t j = 0; j < grid_size_; ++j)
      grid_mat_[j][i] = new_grid_wts_tmp[j];
  }
  
  for(size_t i = 0; i < grid_size_; ++i)
    new_grid_wts[i] = log10_weighted_sum(&(grid_mat_[i][0]),
					 &(snp_wts_[0]),
					 snpVec.size())
      - log10_lik;
}

void gene_eQTL::em_update_config(const vector<double> & grid_wts,
				 vector<double> & new_global_config_prior)
{
  vector<double> new_global_config_prior_tmp(config_size_, NaN);
  for(size_t i = 0; i < snpVec.size(); ++i){
    snp_wts_[i] = snpVec[i].snp_prior;
    snpVec[i].em_update_config(grid_wts, new_global_config_prior_tmp);
    for(size_t j = 0; j < config_size_; ++j)
      config_mat_[j][i] = new_global_config_prior_tmp[j];
  }
  
  for(size_t i = 0; i < config_size_; ++i)
    new_global_config_prior[i] = log10_weighted_sum(&(config_mat_[i][0]),
						    &(snp_wts_[0]),
						    snpVec.size())
      - log10_lik;
}

void gene_eQTL::compute_posterior(const double & pi0,
				  const vector<double> & grid_wts,
				  const vector<double> & global_config_prior)
{
  // gene level posterior
  post_prob_gene = pow(10, log10(1.0-pi0) + log10_BF - log10_lik);
  if(post_prob_gene > 1.0)
    post_prob_gene = 1.0;
  
  // posterior for snp eQTL P(S = i, Z =1 | Y)
  for(size_t i = 0; i < snpVec.size(); ++i){
    snpVec[i].post_prob_snp = pow(10, log10(1-pi0) + log10(snpVec[i].snp_prior) + snpVec[i].log10_BF - log10_lik);
    if(snpVec[i].post_prob_snp > 1.0)
      snpVec[i].post_prob_snp = 1.0;
  }

  // posterior for configuration
  for(size_t i = 0; i < global_config_prior.size(); ++i){
    double prob = 0;
    for(size_t j = 0; j < snpVec.size(); ++j){
      double bfc = snpVec[j].compute_log10_config_BF(i, grid_wts);
      double sprior = snpVec[j].snp_prior * (1-pi0) * global_config_prior[i];
      prob += sprior * pow(10, bfc - log10_lik);
    }

    if(prob > 1.0)
      prob = 1.0;
    post_prob_config.push_back(config_prob(i+1, prob));
  }
}

void gene_eQTL::print_result(const vector<double> & grid_wts)
{
  for(size_t i = 0; i < snpVec.size(); ++i){
    printf("%12s\t%7.3f\t%9.3f\t\t", name_.c_str(), post_prob_gene, log10_BF);
    snpVec[i].print_result(grid_wts);
    printf("\n");
  }
}
