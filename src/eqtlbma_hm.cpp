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
 */

#include <cmath>
#include <cstring>

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

#include "hm_classes.h"

//-----------------------------------------------------------------------------

class eQTL_controller
{
  void load_data_init_format(char *filename, size_t csize, size_t gsize);
  void load_data_new_format_one_file(const string & file,
				     size_t csize, size_t gsize,
				     vector<string> & lines,
				     vector<string> & configs_tokeep);
  void load_data_new_format(const string & file_pattern, 
			    size_t csize, size_t gsize,
			    vector<string> & configs_tokeep);
  
public:
  // storage
  vector<gene_eQTL> geqVec;
  
  double fixed_pi0_;
  
  // parameters need to be estimated
  double *pi0; // non-eqtl prob
  double new_pi0;
  
  double *grid_wts;
  double *new_grid_wts;
  
  double *config_prior;
  double *new_config_prior; 
  
  double *param_est;
  size_t param_size;
  
  // confidence intervals
  double left_pi0, right_pi0;
  vector<double> left_configs, right_configs;
  vector<double> left_grids, right_grids;
  
  int nthread; 
  
  size_t types;
  size_t grid_size;
  size_t config_size;
  
  int output_option;
  
  int verbose;
  
  vector<string> type_vec;
  
  void load_data(char *filename, size_t csize, size_t gsize,
		 vector<string> & configs_tokeep);
  ~eQTL_controller();
  
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
  void save_result(const string & out_file, const bool & skip_bf);
};

eQTL_controller::~eQTL_controller()
{
  delete[] pi0;
  delete[] grid_wts;
  delete[] new_grid_wts;
  if(config_prior != NULL){
    delete[] config_prior;
    delete[] new_config_prior;
  }
}

void eQTL_controller::load_data_init_format(
  char *filename, 
  size_t csize,
  size_t gsize)
{
  gene_eQTL geq;
  
  string stype;
  
  string curr_sid;   // current gene-snp 
  string curr_gene;  // current gene
  string curr_snp;
  
  ifstream infile(filename);
  string line;
  char delim[] = " ";
  char *gene_id, *snp_id;
  char sid[128];
  
  vector<vector<double> > sbf_vec;
  
  while(getline(infile,line)){
     
    const char *cline = line.c_str();
    char *c_line = new char[strlen(cline)+1];
    memset(c_line,0,strlen(cline)+1);
    strcpy(c_line,cline);

    char *token;
  
    
    // parse the line
    gene_id = strtok(c_line, delim);
    snp_id = strtok(NULL, delim);
    memset(sid, 0, 128);
    snprintf(sid, 128, "%s_|_%s", snp_id, gene_id);
    
    if(strcmp(sid, curr_sid.c_str())!=0){
      
      // if gene changes
      if(strcmp(gene_id,curr_gene.c_str())!=0){
	
	if(sbf_vec.size()!=config_size && sbf_vec.size()!=0){
	  fprintf(stderr, "Error: %s contains information for %d types of model BFs, expecting %zu \n", curr_sid.c_str(),int(sbf_vec.size()),config_size);
	  exit(1);
	}

	if(sbf_vec.size()!=0){
	  // process current recorded gene
	  geq.snpVec.push_back(snp_eQTL(curr_snp,sbf_vec,config_prior,grid_wts));
	  geqVec.push_back(geq);
	}
	
	curr_gene = string(gene_id);
	curr_snp = string(snp_id);
	geq = gene_eQTL(curr_gene,pi0);
	
	// create new gene_eqtl structure
      }else{      
	// if only snp changes
	if(sbf_vec.size()!=0)
	  // create new snp_eqtl structure
	  geq.snpVec.push_back(snp_eQTL(curr_snp,sbf_vec,config_prior,grid_wts));
	
	curr_snp = string(snp_id);
      }
      
      curr_sid = string(sid);
      
      // clean up
      for(size_t i=0;i<sbf_vec.size();i++)
	sbf_vec[i].clear();
      sbf_vec.clear();
      
    }

    char *mtype = strtok(NULL, delim);  
    if(type_vec.size() < csize){
      
      string ts(mtype);
      type_vec.push_back(ts);
      //printf("%s\n",ts.c_str());
    }
    
    vector<double> bfv;
    
    while((token=strtok(0, delim))!=0){
      bfv.push_back(atof(token));
    }
    
    if(bfv.size()!= grid_size){
      fprintf(stderr, "Error: Model %s of %s has %d BFs, expecting %zu\n",mtype,sid,int(bfv.size()),grid_size);
      exit(1);
    }

    sbf_vec.push_back(bfv);
    
    delete[] c_line; 
  }

  // last gene
  if(sbf_vec.size()!=0){  
    
    if(sbf_vec.size()!=config_size){
      fprintf(stderr, "Error: %s contains information for %d types of model BFs, expecting %zu \n", curr_sid.c_str(),int(sbf_vec.size()),config_size);
      exit(1);
    }   
    geq.snpVec.push_back(snp_eQTL(curr_snp,sbf_vec,config_prior,grid_wts));
    geqVec.push_back(geq);
  
  }
}

void eQTL_controller::load_data_new_format_one_file(
  const string & file,
  size_t csize,
  size_t gsize,
  vector<string> & lines,
  vector<string> & configs_tokeep)
{
  // load the whole file at once
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
  vector<vector<double> > sbf_vec; // dim1 is config, dim2 is grid
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
    
    if(type_vec.size() < csize) // record config names once (for output)
      type_vec.push_back(config);
    
    if(gene_id.compare(curr_gene) != 0){ // if new gene
      if(! geq.gene.empty()){
	geq.snpVec.push_back(snp_eQTL(curr_snp, sbf_vec, config_prior, grid_wts));
	geqVec.push_back(geq);
      }
      curr_gene = gene_id;
      geq = gene_eQTL(curr_gene, pi0);
      sbf_vec.clear();
    }
    
    if(snp_id.compare(curr_snp) != 0){ // if same gene but new snp
      if(! sbf_vec.empty())
	geq.snpVec.push_back(snp_eQTL(curr_snp, sbf_vec, config_prior, grid_wts));
      curr_snp = snp_id;
      sbf_vec.clear();
    }
    
    // if same gene and snp but new config
    size_t config_idx = sbf_vec.size();
    sbf_vec.push_back(vector<double> (gsize, NaN));
    size_t j = 0;
    pch = strtok(NULL, " \t,");
    while(j < gsize && pch != NULL){
      sbf_vec[config_idx][j] = atof(pch);
      ++j;
      pch = strtok(NULL, " \t,");
    }
  }
  
  geq.snpVec.push_back(snp_eQTL(curr_snp, sbf_vec, config_prior, grid_wts)); // last snp
  geqVec.push_back(geq); // last gene
}

void eQTL_controller::load_data_new_format(
  const string & file_pattern, 
  size_t csize,
  size_t gsize,
  vector<string> & configs_tokeep)
{
  vector<string> files = glob(file_pattern);
  if(verbose > 0)
    cerr << "nb of input files: " << files.size() << endl << flush;
  
  vector<string> lines;
  for(size_t i = 0; i < files.size(); ++i){
    if(verbose > 0)
      progressBar("", i+1, files.size());
    load_data_new_format_one_file(files[i], csize, gsize, lines, configs_tokeep);
    lines.clear();
  }
  if(verbose > 0)
    cerr << endl << flush;
}

void eQTL_controller::load_data(
  char *filename,
  size_t csize,
  size_t gsize,
  vector<string> & configs_tokeep)
{
  if(verbose > 0)
    fprintf(stderr, "load data ...\n");
  clock_t startTime = clock();
  
  // initialization of class parameters
  pi0 = new double[1];
  grid_size = gsize;
  grid_wts = new double[gsize]; // grid prior starts w/ equal weights
  new_grid_wts = new double[gsize];
  config_size = csize;
  if(config_size > 0){
    config_prior = new double[config_size];
    new_config_prior = new double[config_size];
  }
  
  if(configs_tokeep.empty() && false)
    load_data_init_format(filename, csize, gsize);
  else{
    if(! configs_tokeep.empty()){
      fprintf(stderr, "configurations to keep: %s", configs_tokeep[0].c_str());
      for(size_t i = 1; i < configs_tokeep.size(); ++i)
	fprintf(stderr, " %s", configs_tokeep[i].c_str());
      fprintf(stderr, "\n");
    }
    string file = string(filename);
    load_data_new_format(file, csize, gsize, configs_tokeep);
  }
  
  size_t nb_pairs = 0;
  for(size_t g = 0; g < geqVec.size(); ++g){
    nb_pairs += geqVec[g].snpVec.size();
    // cerr << g+1 << "/" << geqVec.size() << " " << geqVec[g].gene << ": " << geqVec[g].snpVec.size() << " snps" << endl;
    // for(size_t v = 0; v < geqVec[g].snpVec.size(); ++v)
    //   cerr << v+1 << "/" << geqVec[g].snpVec.size() << " " << geqVec[g].snpVec[v].snp << endl;
  }
  
  if(verbose > 0)
    fprintf(stderr, "finish loading %zu genes and %zu gene-snp pairs (%f sec)\n",
	    geqVec.size(),
	    nb_pairs,
	    getElapsedTime(startTime));
}

void eQTL_controller::estimate_profile_ci(const bool & skip_ci)
{
  if(verbose > 0)
    cerr << "compute profile-likelihood confidence intervals ..." << flush;
  clock_t startTime = clock();
  
  double tick = 0.001;
  double max = compute_log10_lik();
  
  // estimate CI for pi0
  double pi0_mle = pi0[0];
  left_pi0 = pi0[0];
  right_pi0 = pi0[0];
  if (! skip_ci) {
    while(left_pi0>=0){
      left_pi0 -= tick;
      if(left_pi0<0){
	left_pi0 = 0;
	break;
      }
      pi0[0] = left_pi0;  
      if(compute_log10_lik()/log10(exp(1))<max/log10(exp(1))-2.0){  
	left_pi0 += tick; 
	break;
      }    
    }
    while(right_pi0<=1){
      right_pi0 += tick;
      if(right_pi0 >1){
	right_pi0 = 1;
	break;
      }
      pi0[0] = right_pi0;  
      if(compute_log10_lik()/log10(exp(1))< max/log10(exp(1))-2.0){
	right_pi0 -= tick;
	break;
      }
    }
  }
  // fprintf(stderr,"\nProfile-likelihood Confidence Intervals\n");
  // fprintf(stderr,"pi0: %7.3e [%7.3e, %7.3e]\n",pi0_mle,left_pi0,right_pi0);
  // fprintf(stderr,"config: ");
  
  pi0[0] = pi0_mle;
  
  // estimate CI for config
  double *config_mle = new double[config_size];
  memcpy(config_mle, config_prior, config_size*sizeof(double));
  left_configs.resize(config_size);
  right_configs.resize(config_size);
  for (size_t i=0;i<config_size;i++){
    right_configs[i] = config_mle[i];
    left_configs[i] = config_mle[i];
    double cp_mle = config_mle[i];
    double st = 1 - cp_mle;
    if (! skip_ci) {
      while(left_configs[i]>=0){
	left_configs[i] -= tick;
	if(left_configs[i]<0){
	  left_configs[i] = 0;
	  break;
	}
	double diff = cp_mle - left_configs[i];
	for(size_t j=0;j<config_size;j++){
	  if(j==i)
	    continue;
	  config_prior[j] = config_mle[j] + diff*config_mle[j]/st;
	}
	config_prior[i] = left_configs[i];
	if(compute_log10_lik()/log10(exp(1))<max/log10(exp(1))-2.0){
	  left_configs[i] += tick;
	  break;
	}
      }
      while(right_configs[i]<=1){
	right_configs[i] += tick;
	if(right_configs[i]>1){
	  right_configs[i] = 1;
	  break;
	}
	double diff = cp_mle - right_configs[i];
	for(size_t j=0;j<config_size;j++){
	  if(j==i)
	    continue;
	  config_prior[j] = config_mle[j]+ diff*config_mle[j]/st;
	}
	config_prior[i] = right_configs[i];
	if(compute_log10_lik()/log10(exp(1))<max/log10(exp(1))-2.0){
	  right_configs[i] -= tick;
	  break;
	}
      }
    }
    // fprintf(stderr,"%s  %.3f [%.3f, %.3f]  ",type_vec[i].c_str(), config_mle[i],left_configs[i], right_configs[i]);
    // fflush(stdout);
  }
  memcpy(config_prior, config_mle,config_size*sizeof(double));
  // fprintf(stderr,"\n");
  delete[] config_mle;
  
  // estimate CI for grid
  // fprintf(stderr,"grid:  ");
  double *grid_mle = new double[grid_size];
  memcpy(grid_mle, grid_wts, grid_size*sizeof(double));
  left_grids.resize(grid_size);
  right_grids.resize(grid_size);
  for (size_t i=0;i<grid_size;i++){
    right_grids[i] = grid_mle[i];
    left_grids[i] = grid_mle[i];
    double cp_mle = grid_mle[i];
    double st = 1 - cp_mle;
    if (! skip_ci) {
      while(left_grids[i]>=0){
	left_grids[i] -= tick;
	if(left_grids[i]<0){
	  left_grids[i] = 0;
	  break;
	}
	double diff = cp_mle - left_grids[i];
	for(size_t j=0;j<grid_size;j++){
	  if(j==i)
	    continue;
	  grid_wts[j] = grid_mle[j] + diff*grid_mle[j]/st;
	}
	grid_wts[i] = left_grids[i];
	
	if(compute_log10_lik()/log10(exp(1))<max/log10(exp(1))-2.0){
	  left_grids[i] += tick;
	  break;
	}
      }
      while(right_grids[i]<=1){
	right_grids[i] += tick;
	if(right_grids[i]>1){
	  right_grids[i] = 1;
	  break;
	}
	double diff = cp_mle - right_grids[i];
	for(size_t j=0;j<grid_size;j++){
	  if(j==i)
	    continue;
	  grid_wts[j] = grid_mle[j]+ diff*grid_mle[j]/st;
	}
	grid_wts[i] = right_grids[i];
	if(compute_log10_lik()/log10(exp(1))<max/log10(exp(1))-2.0){
	  right_grids[i] -= tick;
	  break;
	  
	}
      }
    }
    // fprintf(stderr,"%.3f [%.3f, %.3f]  ",grid_mle[i],left_grids[i], right_grids[i]);
    // fflush(stdout);
  }
  memcpy(grid_wts, grid_mle, grid_size*sizeof(double));
  // fprintf(stderr,"\n");
  delete[] grid_mle;
  
  if(verbose > 0)
    cerr << " (" << getElapsedTime(startTime) << " sec)" << endl;
}

void eQTL_controller::init_params(char *init_file)
{
  ifstream ifile(init_file);
  string line;
  istringstream ins;
  
  while(getline(ifile,line)){

    ins.clear();
    ins.str(line);
    double value;
    string type;
    ins>>type;
    
    if(strcmp(type.c_str(),"pi0")==0||strcmp(type.c_str(),"Pi0")==0){
      ins>>value;
      pi0[0] = value;
      continue;
    }
    
    if(strcmp(type.c_str(),"config")==0){
      for(size_t i=0;i<config_size;i++){
	ins>>value;
	config_prior[i] = value;
      }
      continue;
    }
    
    if(strcmp(type.c_str(),"grid")==0){
      for(size_t i=0;i<grid_size;i++){
        ins>>value;
	grid_wts[i] = value;
      }
      continue;
    }
  }

#pragma omp parallel for num_threads(nthread)
  for(int i=0;i<geqVec.size();i++)
    geqVec[i].set_snp_prior();

  //printf("%f %f %f\n",pi0[0],config_prior[0],grid_wts[0]);
  
  ifile.close();
}

void eQTL_controller::init_params(int option)
{
#pragma omp parallel for num_threads(nthread)
  for(int i=0;i<geqVec.size();i++)
    geqVec[i].set_snp_prior();
  
  if(option==1){
    randomize_parameter_sp();
    return;
  }
  
  pi0[0] = 0.5;
  
  if(config_size == 1){
    config_prior[0] = 1.0;
    new_config_prior[0] = 1.0;
  }else{
    for(size_t i=0;i<config_size;i++){
      config_prior[i] = 1.0/(config_size);
    }
    
  }
  
  for(size_t i=0;i<grid_size;i++){
    grid_wts[i] = 1.0/grid_size;
  }
}

void eQTL_controller::randomize_parameter_sp()
{
  const gsl_rng_type * T;
  gsl_rng * r;
       
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  long seed = time (NULL);    
  gsl_rng_set (r, seed); 
  
  //double x1 = gsl_ran_exponential(r,1.0);
  //double x2 = gsl_ran_exponential(r,1.0);
  //pi0[0] = .90 + .1*x1/(x1 + x2);
  //pi0[0] = x1/(x1 + x2);
  pi0[0] = gsl_rng_uniform (r);
  
  double sum = 0;
  for(size_t i=0;i<grid_size;i++){
    grid_wts[i] = gsl_ran_exponential(r,1.0);
    sum += grid_wts[i];
  }

  for(size_t i=0;i<grid_size;i++){
    grid_wts[i] /= sum;
  }

  sum = 0;
  if(config_size==1){
    new_config_prior[0]=config_prior[0] = 1.0;
  }else{
    for(size_t i=0;i<config_size;i++){
      config_prior[i] = gsl_ran_exponential(r,1.0);
      sum+= config_prior[i];
    }

    for(size_t i=0;i<config_size;i++){
      config_prior[i] /= sum;
    }
  }
  
  gsl_rng_free (r);

  return;
}

double eQTL_controller::compute_log10_lik()
{
  double sum = 0;
  
#pragma omp parallel for num_threads(nthread) reduction(+:sum)
  for (int i=0;i<geqVec.size();i++){
    sum += geqVec[i].compute_log10_lik();
  }
  
  return sum;
}

void eQTL_controller::em_update_snp_prior()
{
#pragma omp parallel for num_threads(nthread)
  for (int i=0;i<geqVec.size();i++){
     geqVec[i].em_update_snp_prior();
  }
}

void eQTL_controller::em_update_pi0()
{
  if (fixed_pi0_ == -1.0){
    double sum = 0;
#pragma omp parallel for num_threads(nthread) reduction(+:sum)
    for(int i=0;i<geqVec.size();i++){
      //geqVec[i].update_snp_prior();
      sum += geqVec[i].em_update_pi0();
    }
    new_pi0 = sum/geqVec.size();
  }
  else
    new_pi0 = fixed_pi0_;
}

void eQTL_controller::em_update_grid()
{
  double **matrix = new double *[grid_size];
  double *gene_wts = new double[geqVec.size()];
  for(size_t i=0;i<grid_size;i++){
    matrix[i] = new double[geqVec.size()];
    memset(matrix[i],0,sizeof(double)*geqVec.size());
  }

  double *new_grid = new double[grid_size];
  memset(new_grid,0,sizeof(double)*grid_size);

#pragma omp parallel for num_threads(nthread)
  for(int i=0;i<geqVec.size();i++){
    gene_wts[i] = 1.0;
    double *rst = geqVec[i].em_update_grid(grid_size);
    for(size_t j=0;j<grid_size;j++)
      matrix[j][i] = rst[j];
    delete[] rst;
  }

  for(size_t i=0;i<grid_size;i++){
    new_grid[i] = log10_weighted_sum(matrix[i],gene_wts,geqVec.size())+log10(grid_wts[i]);
    delete[] matrix[i];
  }

  delete[] matrix;
  delete[] gene_wts;

  double *wts = new double[grid_size];
  for(size_t i=0;i<grid_size;i++){
    wts[i] = 1.0;
  }

  double sum = log10_weighted_sum(new_grid,wts,grid_size);
  delete[] wts;

  for(size_t i=0;i<grid_size;i++)
    new_grid_wts[i] = pow(10, (new_grid[i] - sum));

  delete[] new_grid;
}

void eQTL_controller::em_update_config()
{
  if(config_size==1)
    return;
  
  double **matrix = new double *[config_size];
  double *gene_wts = new double[geqVec.size()];
  for(size_t i=0;i<config_size;i++){
    matrix[i] = new double[geqVec.size()];
    memset(matrix[i],0,sizeof(double)*geqVec.size());
  }
  
  double *new_config = new double[config_size];
  memset(new_config,0,sizeof(double)*config_size);
  
#pragma omp parallel for num_threads(nthread)
  for(int i=0;i<geqVec.size();i++){
    gene_wts[i] = 1.0;   
    double *rst = geqVec[i].em_update_config(config_size);
    for(size_t j=0;j<config_size;j++)
      matrix[j][i] = rst[j];
    delete[] rst;
  }
  
  for(size_t i=0;i<config_size;i++){
    new_config[i] = log10_weighted_sum(matrix[i],gene_wts,geqVec.size())+log10(config_prior[i]);
    delete[] matrix[i];
  }
  
  delete[] matrix;
  delete[] gene_wts;
  
  double *wts = new double[config_size];
  for(size_t i=0;i<config_size;i++){
    wts[i] = 1.0;
  }

  double sum = log10_weighted_sum(new_config,wts,config_size);
  delete[] wts;
    
  for(size_t i=0;i<config_size;i++)
    new_config_prior[i] = pow(10, (new_config[i] - sum));
  
  delete[] new_config;
}

void eQTL_controller::update_params()
{
  *pi0 = new_pi0;
  if(config_size>1){
    for(size_t i=0;i<config_size;i++)
      config_prior[i] = new_config_prior[i];
  }
  
  for(size_t i=0;i<grid_size;i++)
    grid_wts[i] = new_grid_wts[i];
  
  //for (size_t i=0;i<geqVec.size();i++){
  //  geqVec[i].update_snp_prior();
  //}
}

void eQTL_controller::run_EM(double thresh)
{
  if(verbose > 0)
    cerr << "run EM algorithm ..." << endl << flush;

  // start iteration
  time_t startRawTime, endRawTime;
  time (&startRawTime);
  
  size_t count = 1;
  double last_log10_lik = -9999999;
  //new_pi0 = *pi0 = 0.50;
  while(true){
    
    double curr_log10_lik = compute_log10_lik();
    em_update_pi0();    
    em_update_config();
    em_update_grid();
   
    //em_update_snp_prior();

    update_params();
    
    if(output_option==1){
      // output 
      fprintf(stderr,"iter  %4zu  loglik = %.4f   ",count++,curr_log10_lik);
      fprintf(stderr,"pi0  %7.3e     config ", pi0[0]);
	      
      for(size_t i=0;i<config_size;i++)
	fprintf(stderr,"%7.3e  ",new_config_prior[i]);

      fprintf(stderr,"  grid ");

      for(size_t i=0;i<grid_size;i++)
      	fprintf(stderr,"%7.3e  ",new_grid_wts[i]);

      fprintf(stderr,"\n");
    }
    
    if(fabs(curr_log10_lik-last_log10_lik)<thresh)
      break;
    
    last_log10_lik = curr_log10_lik;
  }

  time (&endRawTime);
  cerr << "EM ran for " << getElapsedTime(startRawTime, endRawTime) << endl << flush;
}

void eQTL_controller::compute_posterior()
{
  if(verbose > 0)
    cerr << "compute posteriors ..." << flush;
  clock_t startTime = clock();
  
#pragma omp parallel for num_threads(nthread)
  for(int i=0;i<geqVec.size();i++)
    geqVec[i].compute_posterior(config_prior, config_size);
  
  if(verbose > 0)
    cerr << " (" << getElapsedTime(startTime) << " sec)" << endl;
}

void eQTL_controller::print_result()
{
  printf("gene\tgene.posterior.prob\tgene.log10.bf\t\tsnp\tsnp.log10.bf\t\t");
  for(size_t i=0;i<type_vec.size();i++){
    printf("log10.bf.%s\t",type_vec[i].c_str());
  }
  printf("\n");
  
  for(size_t i=0;i<geqVec.size();i++)
    geqVec[i].print_result();
}

void eQTL_controller::save_result(const string & out_file,
				  const bool & skip_bf)
{
  if(verbose > 0)
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
  txt << "#pi0\t" << scientific << pi0[0] << "\t" << left_pi0 << "\t" << right_pi0 << endl;
  ++nb_lines;
  gzwriteLine(stream, txt.str(), out_file, nb_lines);
  for(size_t i = 0; i < type_vec.size(); ++i){
    txt.str("");
    txt << "#config." << type_vec[i].c_str()
	<< "\t" << config_prior[i]
	<< "\t" << left_configs[i]
	<< "\t" << right_configs[i]
	<< endl;
    ++nb_lines;
    gzwriteLine(stream, txt.str(), out_file, nb_lines);
  }
  for(size_t i = 0; i < grid_size; ++i){
    txt.str("");
    txt << "#grid." << i+1
	<< "\t" << grid_wts[i]
	<< "\t" << left_grids[i]
	<< "\t" << right_grids[i]
	<< endl;
    ++nb_lines;
    gzwriteLine(stream, txt.str(), out_file, nb_lines);
  }

  if(! skip_bf){
    txt.str("");
    txt << "gene\tgene.posterior.prob\tgene.log10.bf\tsnp\tsnp.log10.bf\t";
    for(size_t i = 0; i < type_vec.size(); ++i)
      txt << "log10.bf." << type_vec[i].c_str() << "\t";
    txt << "\n";
    ++nb_lines;
    gzwriteLine(stream, txt.str(), out_file, nb_lines);
    
    for(size_t i = 0; i < geqVec.size(); ++i){ // loop over gene
      for(size_t j = 0; j < geqVec[i].snpVec.size(); ++j){ // loop over snp
	txt.str("");
	txt << setprecision(4)
	    << geqVec[i].gene
	    << "\t"
	    << geqVec[i].post_prob_gene
	    << "\t"
	    << geqVec[i].compute_log10_BF()
	    << "\t"
	    << geqVec[i].snpVec[j].snp
	    << "\t"
	    << geqVec[i].snpVec[j].compute_log10_BF();
	for(size_t k = 0; k < geqVec[i].snpVec[j].config_size; ++k) // loop over config
	  txt << "\t"
	      << log10_weighted_sum(geqVec[i].snpVec[j].gm[k],
				    geqVec[i].snpVec[j].grid_wts,
				    geqVec[i].snpVec[j].grid_size);
	txt << "\n";
	++nb_lines;
	gzwriteLine(stream, txt.str(), out_file, nb_lines);
      }
    }
  }
  
  closeFile(out_file, stream);
  
  if(verbose > 0)
    cerr << " (" << getElapsedTime(startTime) << " sec)" << endl;
}

//-----------------------------------------------------------------------------

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
       << "  -d, --data\tinput data (usually output files from eqtlbma_bf)" << endl
       << "  -s\t\tnumber of subgroup configurations" << endl
       << "  -g, --grid\tnumber of grid points" << endl
       << "  -o, --out\toutput file (gzipped)" << endl
       << "  -i, --init\tfile for initialization" << endl
       << "  -r, --ran\trandom initialization" << endl
       << "  -t, --thresh\tthreshold to stop the EM (default=0.05)" << endl
       << "      --thread\tnumber of threads (default=1)" << endl
       << "  -c, --ci\t" << endl
       << "      --configs\tsubset of configurations to keep (e.g. \"1|3|1-3\")" << endl
       << "      --getci\tcompute the confidence intervals (single thread, thus slow)" << endl
       << "      --getbf\tcompute the Bayes Factors using the estimated weights" << endl
       << "\t\tcan take some time, otherwise only the estimated weights are reported" << endl
       << "      --pi0\tfixed value for pi0" << endl
       << endl;
}

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

int main(int argc, char **argv)
{
  int verbose = 1;
  
  char data_file[128];
  memset(data_file,0,128);
  size_t csize = string::npos;
  size_t gsize = string::npos;
  string out_file;
  
  int nthread = 1;
  
  char init_file[128];
  int option = 0;
  
  char ci_file[128];
  memset(ci_file,0,128); 
  memset(data_file,0,128); 
  memset(init_file,0,128);
  
  double thresh = 0.05;
  
  vector<string> configs_tokeep;
  
  bool skip_ci = true;
  bool skip_bf = true;
  
  double fixed_pi0 = -1.0;
  
  for(int i=1;i<argc;i++){
    if(strcmp(argv[i], "-h")==0 || strcmp(argv[i], "--help")==0){
      help(argv);
      exit(0);
    }
    if(strcmp(argv[i], "-V")==0 || strcmp(argv[i], "--version")==0){
      version(argv);
      continue;
    }
    if(strcmp(argv[i], "-v")==0 || strcmp(argv[i], "--verbose")==0){
      verbose = atoi(argv[++i]);
      continue;
    }
    if(strcmp(argv[i], "-d")==0 || strcmp(argv[i], "--data")==0){
      strcpy(data_file, argv[++i]);
      continue;
    }
    if(strcmp(argv[i], "-s")==0){
      csize = atoi(argv[++i]);
      continue;	
    }
    if(strcmp(argv[i], "-g")==0 || strcmp(argv[i], "--grid")==0){
      gsize = atoi(argv[++i]);
      continue;
    }
    if(strcmp(argv[i], "-o")==0 || strcmp(argv[i], "--out")==0){
      out_file = string(argv[++i]);
      continue;
    }
    
    ////////////////////// optional ///////////////////////////
    if(strcmp(argv[i],"-i")==0 || strcmp(argv[i],"--init")==0){
      strcpy(init_file, argv[++i]);
      continue;
    }
    if(strcmp(argv[i], "-r")==0 || strcmp(argv[i],"--ran")==0){
      option = 1;
      continue;
    }
    if(strcmp(argv[i], "-t")==0 || strcmp(argv[i], "--thresh")==0){
      thresh = atof(argv[++i]);
      continue;
    }
    if(strcmp(argv[i], "--thread") == 0){
      nthread = atoi(argv[++i]);
      continue;
    }
    if(strcmp(argv[i], "-c")==0 || strcmp(argv[i],"--ci")==0){
      strcpy(ci_file, argv[++i]);
      continue;
    }
    if(strcmp(argv[i], "--configs")==0){
      split(argv[++i], "|", configs_tokeep);
      continue;
    }
    if(strcmp(argv[i],"--getci")==0){
      skip_ci = false;
      continue;
    }
    if(strcmp(argv[i],"--getbf")==0){
      skip_bf = false;
      continue;
    }
    if(strcmp(argv[i], "--pi0")==0){
      fixed_pi0 = atof(argv[++i]);
      continue;
    }
  }
  
  // checking mandatory arguments
  if(strlen(data_file)==0){
    fprintf(stderr,"Error: data file unspecified\n");
    exit(1);
  }
  if(csize == string::npos){
    fprintf(stderr,"Error: number of model types unspecified\n");
    exit(1);
  }
  if(gsize == string::npos){
    fprintf(stderr,"Error: number of model grids unspecified\n");
    exit(1);
  }
  if(out_file.empty()){
    fprintf(stderr,"Error: output file unspecified\n");
    exit(1);
  }
  
  // fit the hierarchical model via EM
  time_t startRawTime, endRawTime;
  time (&startRawTime);
  if(verbose > 0)
    cerr << "START " << basename(argv[0])
	 << " " << getDateTime(startRawTime) << endl
	 << "version " << VERSION << " compiled " << __DATE__
	 << " " << __TIME__ << endl
	 << "cmd-line: " << getCmdLine(argc, argv) << endl
	 << "cwd: " << getCurrentDirectory() << endl << flush;
  
  eQTL_controller controller;
  controller.nthread = nthread;
  controller.output_option = 1;
  controller.verbose = verbose;
  
  controller.load_data(data_file,csize,gsize,configs_tokeep);
  
  if(strlen(ci_file)>0){
    controller.init_params(ci_file);
    controller.estimate_profile_ci(skip_ci);
    return 0;
  }
  
  if(strlen(init_file)>0)
    controller.init_params(init_file);
  else
    controller.init_params(option);
  
  controller.fix_pi0(fixed_pi0);
  
  controller.run_EM(thresh);
  
  if(! skip_bf)
    controller.compute_posterior();
  
  controller.estimate_profile_ci(skip_ci);
  
  controller.save_result(out_file, skip_bf);
  
  time (&endRawTime);
  cerr << "END " << basename(argv[0]) << " " << getDateTime(endRawTime) << endl
       << "elapsed -> " << getElapsedTime(startRawTime, endRawTime) << endl
       << "max.mem -> " << getMaxMemUsedByProcess2Str () << endl;
}
