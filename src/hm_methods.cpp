/** \file hm_methods.cpp
 *
 *  `hm_methods' contains the methods of the classes defined in `hm_classes.h'.
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

#include "hm_classes.h"
#include <math.h>

snp_eQTL::snp_eQTL(string name, const vector<vector<double> > & gm_vec,double *gcp, double *gw){
  
  snp = name;
  global_config_prior = gcp;
  grid_wts = gw;  
  config_size = gm_vec.size();
  grid_size = gm_vec[0].size();
  
  gm = new double *[config_size];
  
  for(size_t i=0;i<gm_vec.size();i++){
    gm[i] = new double[grid_size];
    for(size_t j=0;j<grid_size;j++)
      gm[i][j] = gm_vec[i][j];
  }
  

}

double *snp_eQTL::em_update_config(){  

  double *rst = new double[config_size];
  for(size_t i=0;i<config_size;i++)
    rst[i] = log10_weighted_sum(gm[i],grid_wts,grid_size);  
  return rst;
}


double *snp_eQTL::em_update_grid(){
  double *rst = new double[grid_size];
  double *config = new double[config_size];
  for(size_t i=0;i<grid_size;i++){
    // for i-th weight function   
    for(size_t j=0;j<config_size;j++)
      config[j] = gm[j][i];
    rst[i] = log10_weighted_sum(config,global_config_prior,config_size);   
  }
  delete[] config;
  return rst;
}



double snp_eQTL::compute_log10_BF(){
  
  // integrate out phi, omega according to current wts
  double *rst = new double[config_size];
  for(size_t i=0;i<config_size;i++)
    rst[i] = log10_weighted_sum(gm[i],grid_wts,grid_size);
  
  double wsum = log10_weighted_sum(rst,global_config_prior,config_size);
  delete[] rst;
  // save for em update
  log10_BF = wsum;
  return wsum;
  
}



double snp_eQTL::compute_log10_config_BF(size_t config){
  return log10_weighted_sum(gm[config],grid_wts,grid_size);  
}


void snp_eQTL::print_result(void){
  
  printf("%12s\t%9.4f\t\t", snp.c_str(), log10_BF);

  double *rst = new double[config_size];
  for(size_t i=0;i<config_size;i++)
    rst[i] = log10_weighted_sum(gm[i],grid_wts,grid_size); 
  
  //double *cprob = new double[config_size];
  for(size_t i=0;i<config_size;i++){
    //cprob[i] = pow(10, (log10(global_config_prior[i])+rst[i]-log10_BF)); 
    //if(cprob[i]>1)
    //  cprob[i] = 1.0;   
    printf("%7.4f\t", rst[i]);
  }
  delete[] rst;
  
}


void gene_eQTL::set_snp_prior(){

  for(size_t i=0;i<snpVec.size();i++){
    // for now
    snpVec[i].snp_prior = 1.0/snpVec.size();    
  }

}


double gene_eQTL::compute_log10_BF(){
   
  double *snp_wts = new double[snpVec.size()];
  double *snp_logBF = new double[snpVec.size()];
  
  for(size_t i=0;i<snpVec.size();i++){
    snp_wts[i] = snpVec[i].snp_prior;
    snp_logBF[i] = snpVec[i].compute_log10_BF();
  }
  
  double rst = log10_weighted_sum(snp_logBF,snp_wts,snpVec.size());
  
  log10_BF = rst;

  delete[] snp_wts;
  delete[] snp_logBF;
  
  return rst;
}
  

double gene_eQTL::compute_log10_lik(){

  double vec[2];
  double wts[2];
  wts[0] = *Pi0;
  wts[1] = 1-(*Pi0);
  vec[0] = 0;
  vec[1] = compute_log10_BF();

  log10_lik = log10_weighted_sum(vec, wts, 2);
  
  return log10_lik;
}

double gene_eQTL::em_update_pi0(){  
  
  return pow(10,(log10(*Pi0)-log10_lik));
}


// last to do
void gene_eQTL::em_update_snp_prior(){
  
  double *snp_wts = new double[snpVec.size()];
  double *snp_logBF = new double[snpVec.size()];
  
  for(size_t i=0;i<snpVec.size();i++){ 
    snp_wts[i] = snpVec[i].snp_prior;
    snp_logBF[i] = snpVec[i].log10_BF;                                                                                                                                  
  }                                                                               
  
  double factor = log10_weighted_sum(snp_logBF,snp_wts,snpVec.size());         
  
  delete[] snp_wts;
  delete[] snp_logBF;

  for(size_t i=0;i<snpVec.size();i++)
    snpVec[i].new_snp_prior = pow(10,(snpVec[i].log10_BF-factor))*snpVec[i].snp_prior;
  
}

void gene_eQTL::update_snp_prior(){
  double sum = 0;
  for(size_t i=0;i<snpVec.size();i++){
    if(snpVec[i].new_snp_prior<0.001)
      snpVec[i].new_snp_prior = 0.001;
    sum += snpVec[i].new_snp_prior;
  }

 for(size_t i=0;i<snpVec.size();i++){
    snpVec[i].snp_prior = snpVec[i].new_snp_prior/sum;
 } 

}



double* gene_eQTL::em_update_grid(size_t grid_size){

  double **matrix = new double *[grid_size];
  double *sprior = new double[snpVec.size()];
  for(size_t i=0;i<grid_size;i++){
    matrix[i] = new double[snpVec.size()];
    memset(matrix[i],0,sizeof(double)*snpVec.size());
  }

  double *new_grid = new double[grid_size];
  memset(new_grid,0,sizeof(double)*grid_size);

  for(size_t i=0;i<snpVec.size();i++){
    sprior[i] = snpVec[i].snp_prior;
    double *rst = snpVec[i].em_update_grid();
    for(size_t j=0;j<grid_size;j++)
      matrix[j][i] = rst[j];
    delete[] rst;
  }

  for(size_t i=0;i<grid_size;i++){
    new_grid[i] = log10_weighted_sum(matrix[i],sprior,snpVec.size())-log10_lik;
    delete[] matrix[i];
  }
  

  delete[] matrix;
  delete[] sprior;

  return new_grid;
}


double* gene_eQTL::em_update_config(size_t config_size){
  
  double **matrix = new double *[config_size];
  double *sprior = new double[snpVec.size()];
  for(size_t i=0;i<config_size;i++){
    matrix[i] = new double[snpVec.size()];
    memset(matrix[i],0,sizeof(double)*snpVec.size());
  }
  
  double *new_config = new double[config_size];
  memset(new_config,0,sizeof(double)*config_size);
  
  for(size_t i=0;i<snpVec.size();i++){
    sprior[i] = snpVec[i].snp_prior;
    double *rst = snpVec[i].em_update_config();
    for(size_t j=0;j<config_size;j++)
      matrix[j][i] = rst[j];
    delete[] rst;
  }
  
  for(size_t i=0;i<config_size;i++){
    new_config[i] = log10_weighted_sum(matrix[i],sprior,snpVec.size())-log10_lik;
    delete[] matrix[i];
  }

  delete[] matrix;
  delete[] sprior;
  
  return new_config;
}


void gene_eQTL::compute_posterior(double *config_prior, size_t config_size){

  // gene level posterior
  post_prob_gene = pow(10,(log10(1.0-(*Pi0))+log10_BF-log10_lik));
  if(post_prob_gene > 1.0)
    post_prob_gene = 1.0;
  
  // posterior for snp eQTL P(S = i, Z =1 | Y)
  for(size_t i=0;i<snpVec.size();i++){
    snpVec[i].post_prob_snp = pow(10,(log10(1-(*Pi0))+log10(snpVec[i].snp_prior)+snpVec[i].log10_BF-log10_lik));
    if(snpVec[i].post_prob_snp>1.0)
      snpVec[i].post_prob_snp = 1.0;
  }
  


  // posterior for configuration
  for(size_t i=0;i<config_size;i++){
    double prob = 0;
    for(size_t j=0;j<snpVec.size();j++){
      double bfc = snpVec[j].compute_log10_config_BF(i);
      double sprior = snpVec[j].snp_prior*(1-(*Pi0))*config_prior[i];
      prob += sprior*pow(10,bfc-log10_lik);
    }

    if(prob > 1.0)
      prob = 1.0;
    post_prob_config.push_back(config_prob(i+1,prob));
  }
  
  
}


void gene_eQTL::print_result(void){
  
  for(size_t i=0;i<snpVec.size();i++){
    printf("%12s\t%7.3f\t%9.3f\t\t",gene.c_str(),post_prob_gene,log10_BF);
    snpVec[i].print_result();
    printf("\n");
  }
}


bool operator< (const config_prob& a, const config_prob& b){  
  return a.prob < b.prob;
};



// This computes log_{10}(\sum_i w_i 10^vec_i)
double log10_weighted_sum(double* vec, double *wts, size_t size){

  double max = vec[0];
  for(size_t i=1;i<size;i++){
    if(vec[i]>max)
      max = vec[i];
  }
  double sum = 0;
  for(size_t i=0;i<size;i++){
    sum += wts[i]*pow(10, (vec[i]-max));
  }
  
  return (max+log10(sum));
}


  
