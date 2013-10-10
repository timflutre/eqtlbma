/** \file gene_snp_pair.cpp
 *
 *  `GeneSnpPair' is a class 
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

#include <numeric>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_machine.h> // GSL_DBL_EPSILON
#include <gsl/gsl_statistics_double.h>

#include "quantgen/gene_snp_pair.hpp"

using namespace std;

using namespace utils;

namespace quantgen {

  GeneSnpPair::GeneSnpPair(const string & gene_name, const string & snp_name)
  {
    gene_name_ = gene_name;
    snp_name_ = snp_name;
  }

  bool GeneSnpPair::HasResults(const string & subgroup) const
  {
    bool res = false;
    if(subgroup2samplesize_.find(subgroup) != subgroup2samplesize_.end()
       && subgroup2samplesize_.find(subgroup)->second > 0)
      res = true;
    return res;
  }

/** \brief
 *  \param Y it contains one vector of expression levels per subgroup
 *  \param Xg if there is a single subgroup or if there are several subgroups 
 *   and same_individuals==true, it contains one vector of genotypes; if there 
 *   are several subgroups and same_individuals==false, it contains one vector
 *   of genotypes per subgroup (always at the same SNP)
 *  \param Xc it contains one vector of values per covariate per subgroup; 
 *   whether there is any covariate or not, the intercept is not added
 */
  void GeneSnpPair::FillStlContainers(
    const Samples & samples,
    const Gene & gene,
    const Snp & snp,
    const Covariates & covariates,
    const vector<string> & subgroups,
    const bool & same_individuals,
    const bool & need_qnorm,
    const gsl_permutation * perm,
    vector<vector<double> > & Y,
    vector<vector<double> > & Xg,
    vector<vector<vector<double> > > & Xc,
    vector<string> & subgroups_with_data)
  {
    string subgroup;
    for(vector<string>::const_iterator it_sbgrp = subgroups.begin();
	it_sbgrp != subgroups.end(); ++it_sbgrp){
      subgroup = *it_sbgrp;
    
      vector<double> Y_tmp, Xg_tmp;
      vector<vector<double> > Xc_tmp;
      size_t idx_all, idx_explevel, idx_genotype, idx_covar;
      vector<size_t> indices_all; // used for covariates
      for(vector<string>::const_iterator it_all = samples.begin();
	  it_all != samples.end(); ++it_all){
	if(perm == NULL)
	  idx_all = it_all - samples.begin();
	else
	  idx_all = gsl_permutation_get(perm, it_all - samples.begin());
	idx_explevel = samples.GetIndexExplevel(idx_all, subgroup);
	idx_genotype = samples.GetIndexGenotype(it_all - samples.begin(), subgroup);
	if(idx_explevel != string::npos && idx_genotype != string::npos
	   && ! isNan(gene.GetExplevel(subgroup, idx_explevel))){
	  Y_tmp.push_back(gene.GetExplevel(subgroup, idx_explevel));
	  Xg_tmp.push_back(snp.GetGenotype(subgroup, idx_genotype));
	  indices_all.push_back(it_all - samples.begin());
	}
      }
    
      if(Y_tmp.empty())
	continue;
    
      subgroups_with_data.push_back(subgroup);
    
      if(need_qnorm)
	qqnorm(&Y_tmp[0], Y_tmp.size());
      Y.push_back(Y_tmp);
    
      if((! same_individuals) || (same_individuals && Xg.size() == 0)){
	Xg.push_back(Xg_tmp);
      
	for(map<string,vector<double> >::const_iterator it_covars
	      = covariates.begin(subgroup);
	    it_covars != covariates.end(subgroup); ++it_covars){
	  vector<double> covar_values;
	  for(vector<size_t>::const_iterator it_idx = indices_all.begin();
	      it_idx != indices_all.end(); ++it_idx){
	    idx_covar = samples.GetIndexCovariate(*it_idx, subgroup);
	    if(idx_covar == string::npos){
	      cerr << "ERROR: missing covariate " << it_covars->first
		   << " of gene " << gene_name_ << " and snp " << snp_name_ 
		   << " for sample " << samples.GetSample(*it_idx)
		   << " from subgroup " << subgroup << endl;
	      exit(1);
	    }
	    covar_values.push_back(
	      covariates.GetCovariate(subgroup, it_covars->first, idx_covar));
	  }
	  Xc_tmp.push_back(covar_values);
	}
	Xc.push_back(Xc_tmp);
      }
    
      subgroup2samplesize_.insert(make_pair(subgroup, Y_tmp.size()));
      subgroup2nbcovariates_.insert(make_pair(subgroup, Xc_tmp.size()));
      subgroup2pve_.insert(make_pair(subgroup, NaN));
      subgroup2sigmahat_.insert(make_pair(subgroup, NaN));
      subgroup2sstats_.insert(make_pair(subgroup, vector<double>(3,NaN)));
    }
  }

/** \brief
 *  \note both gene and snp are supposed to have data from the given subgroup
 */
  void GeneSnpPair::CalcSstatsOneSbgrp(const Samples & samples,
				       const Gene & gene,
				       const Snp & snp,
				       const Covariates & covariates,
				       const string & subgroup,
				       const string & likelihood,
				       const bool & need_qnorm,
				       const gsl_permutation * perm)
  {
    vector<vector<double> > Y, Xg;
    vector<vector<vector<double> > > Xc;
    vector<string> subgroups_with_data;
    FillStlContainers(samples, gene, snp, covariates, vector<string>(1, subgroup),
		      false, need_qnorm, perm, Y, Xg, Xc, subgroups_with_data);
  
    if(likelihood.compare("normal") == 0){
      size_t N = Xg[0].size(); // nb of individuals
      gsl_matrix * X = gsl_matrix_alloc(N, 1 + 1 + Xc[0].size());
      gsl_vector * y = gsl_vector_alloc(N);
      for(size_t i = 0; i < N; ++i){
	gsl_vector_set(y, i, Y[0][i]);
	gsl_matrix_set(X, i, 0, 1.0); // intercept
	gsl_matrix_set(X, i, 1, Xg[0][i]);
	for(size_t j = 0; j < Xc[0].size(); ++j)
	  gsl_matrix_set(X, i, j+2, Xc[0][j][i]);
      }
      FitSingleGeneWithSingleSnp(X, y, subgroup2pve_[subgroup],
				 subgroup2sigmahat_[subgroup],
				 subgroup2sstats_[subgroup][0],
				 subgroup2sstats_[subgroup][1],
				 subgroup2sstats_[subgroup][2]);
      gsl_matrix_free(X);
      gsl_vector_free(y);
    }
    else if(likelihood.find("poisson") != string::npos){
      vector<vector<double> > X(1 + Xc[0].size(), vector<double>(Xg[0].size(), NaN)); // P x N
      for(size_t i = 0; i < X[0].size(); ++i) // fill genotypes
	X[0][i] = Xg[0][i];
      for(size_t j = 1; j < X.size(); ++j) // fill other covariates
	for(size_t i = 0; i < X[0].size(); ++i)
	  X[j][i] = Xc[0][j-1][i];
      IRLS irls("log-link");
      if(likelihood.compare("quasipoisson") == 0){
	irls.link->quasi = true;
      }
      else{
	irls.link->quasi = false;
      }
      vector<double> offv(Y[0].size(), 0.0);
      irls.load_data(Y[0], X, offv);
      irls.fit_model();
      vector<double> coef = irls.get_coef(),
	se_coef = irls.get_stderr();
      subgroup2sigmahat_[subgroup] = sqrt(irls.get_dispersion());
      subgroup2sstats_[subgroup][0] = coef[1];
      subgroup2sstats_[subgroup][1] = se_coef[1];
      if(likelihood.compare("quasipoisson") == 0){
	subgroup2sstats_[subgroup][2] = 2 * gsl_cdf_tdist_P(
	  -fabs(coef[1] / se_coef[1]), Xg[0].size() - irls.get_rank_X());
      } else{
	subgroup2sstats_[subgroup][2] = 2 * gsl_cdf_gaussian_P(
	  -fabs(coef[1] / se_coef[1]), 1.0);
      }
    }
  }

  void GeneSnpPair::StandardizeSstatsAndCorrectSmallSampleSize(
    map<string,vector<double> > & subgroup2stdsstats)
  {
    string subgroup;
  
    for(map<string,vector<double> >::const_iterator it
	  = subgroup2sstats_.begin(); it != subgroup2sstats_.end(); ++it) {
      subgroup = it->first;
      subgroup2stdsstats.insert(make_pair(subgroup, vector<double>(3,NaN)));
    
      double N = subgroup2samplesize_[subgroup],
	bhat = subgroup2sstats_[subgroup][0] / subgroup2sigmahat_[subgroup],
	sebhat = subgroup2sstats_[subgroup][1] / subgroup2sigmahat_[subgroup],
	t = bhat / sebhat;
    
      // apply quantile-quantile transformation (Wen and Stephens, arXiv 2011)
      double nu = N - 2 - subgroup2nbcovariates_[subgroup]; // degrees of freedom
      t = gsl_cdf_gaussian_Pinv(gsl_cdf_tdist_P(-fabs(bhat/sebhat), nu), 1.0);
    
      if(fabs(t) > 1e-8) {
	double sigmahat = fabs(subgroup2sstats_[subgroup][0])
	  / (fabs(t) * sebhat);
	bhat = subgroup2sstats_[subgroup][0] / sigmahat,
	  sebhat = bhat / t;
      }
      else {
	bhat = 0;
	sebhat = numeric_limits<double>::infinity();
      }
      subgroup2stdsstats[subgroup][0] = bhat;
      subgroup2stdsstats[subgroup][1] = sebhat;
      subgroup2stdsstats[subgroup][2] = t;
    }
  }

/** \brief Return the log10 of the approximate Bayes Factor 
 *  from Wen and Stephens (arXiv 2011)
 *  \note this is the univariate version of the ABF
 *  \note gamma[s]=1 means the eQTL is present in subgroup s
 */
  double CalcLog10AbfUvlr(
    const vector<int> & gamma,
    const vector<vector<double> > & stdsstats,
    const double phi2,
    const double oma2)
  {
    double l10AbfAll = 0.0, bhat = 0.0, varbhat = 0.0, t = 0.0,
      bbarhat_num = 0.0, bbarhat_denom = 0.0, varbbarhat = 0.0;
    vector<double> l10AbfsSingleSbgrp;
  
    for(size_t s = 0; s < gamma.size(); ++s) {
      if(gamma[s] == 0)
	continue; // skip inactive eQTL or unexpressed gene
      bhat = stdsstats[s][0];
      varbhat = pow(stdsstats[s][1], 2);
      t = stdsstats[s][2];
      double lABF_single;
      if(fabs(t) < 1e-8) {
	lABF_single = 0;
      }
      else {
	bbarhat_num += bhat / (varbhat + phi2);
	bbarhat_denom += 1 / (varbhat + phi2);
	varbbarhat += 1 / (varbhat + phi2);
	lABF_single = 0.5 * log10(varbhat)
	  - 0.5 * log10(varbhat + phi2)
	  + (0.5 * pow(t,2) * phi2 / (varbhat + phi2)) / log(10);
      }
      l10AbfsSingleSbgrp.push_back(lABF_single);
#ifdef DEBUG
      fprintf(stderr, "l10AbfsSingleSbgrp[%zu]=%e\n",
	      s+1, l10AbfsSingleSbgrp[s]);
#endif
    }
  
    double bbarhat = (bbarhat_denom != 0.0) ?
      bbarhat_num / bbarhat_denom
      : 0.0;
    varbbarhat = (varbbarhat != 0.0) ?
      1 / varbbarhat
      : numeric_limits<double>::infinity();
    if(bbarhat != 0.0 && varbbarhat < numeric_limits<double>::infinity()) {
      double T2 = pow(bbarhat, 2.0) / varbbarhat;
      double lABF_bbar = (T2 != 0) ?
	0.5 * log10(varbbarhat) - 0.5 * log10(varbbarhat + oma2)
	+ (0.5 * T2 * oma2 / (varbbarhat + oma2)) / log(10)
	: 0;
#ifdef DEBUG
      fprintf(stderr, "bbarhat=%e varbbarhat=%e T2=%e lABF_bbar=%e\n",
	      bbarhat, varbbarhat, T2, lABF_bbar);
#endif
      l10AbfAll = lABF_bbar;
      for(size_t i = 0; i < l10AbfsSingleSbgrp.size(); ++i)
	l10AbfAll += l10AbfsSingleSbgrp[i];
    }
    else // minor allele freq very close, or equal, to 0
      l10AbfAll = 0.0;
  
    return l10AbfAll;
  }

/** \brief Calculate the ABF for the consistent configuration,
 *  averaged over the large grid for the general model, the fixed-effect
 *  model and the maximum-heterogeneity model
 *  \note the configuration is represented by 'gen', 'gen-fix'
 *  and 'gen-maxh'
 */
  void GeneSnpPair::CalcAbfsUvlrForConsistentConfiguration(
    const Grid & grid,
    const map<string,vector<double> > & subgroup2stdsstats,
    const vector<string> & subgroups)
  {
    vector<int> gamma(subgroups.size(), string::npos);
    vector<vector<double> > stdsstats;
    for(vector<string>::const_iterator it = subgroups.begin();
	it != subgroups.end(); ++it) {
      if(HasResults(*it)) {
	gamma[it - subgroups.begin()] = 1;
	stdsstats.push_back(subgroup2stdsstats.find(*it)->second);
      }
      else { // e.g. when gene is absent from some subgroups
	gamma[it - subgroups.begin()] = 0;
	stdsstats.push_back(vector<double>(3, 0));
      }
    }
  
    vector<double> l10_abfs(grid.size(), NaN),
      l10_abfs_fix(grid.size(), NaN),
      l10_abfs_maxh(grid.size(), NaN);
  
    for(size_t grid_id = 0; grid_id < grid.size(); ++grid_id) {
      l10_abfs[grid_id] = CalcLog10AbfUvlr(gamma,
					   stdsstats,
					   grid.phi2s[grid_id],
					   grid.oma2s[grid_id]);
      l10_abfs_fix[grid_id] = CalcLog10AbfUvlr(gamma,
					       stdsstats,
					       0.0,
					       grid.phi2s[grid_id]
					       + grid.oma2s[grid_id]);
      l10_abfs_maxh[grid_id] = CalcLog10AbfUvlr(gamma,
						stdsstats,
						grid.phi2s[grid_id]
						+ grid.oma2s[grid_id],
						0.0);
    }
  
    unweighted_abfs_.insert(make_pair("gen", l10_abfs));
    unweighted_abfs_.insert(make_pair("gen-fix", l10_abfs_fix));
    unweighted_abfs_.insert(make_pair("gen-maxh", l10_abfs_maxh));
  
    weighted_abfs_.insert(
      make_pair("gen", log10_weighted_sum(&(l10_abfs[0]), l10_abfs.size())));
    weighted_abfs_.insert(
      make_pair("gen-fix",
		log10_weighted_sum(&(l10_abfs_fix[0]), l10_abfs_fix.size())));
    weighted_abfs_.insert(
      make_pair("gen-maxh",
		log10_weighted_sum(&(l10_abfs_maxh[0]), l10_abfs_maxh.size())));
  }

/** \brief Calculate the ABF for each singleton, averaged over the small grid
 *  \note the configuration is represented by the index of the given subgroup
 *  (eg. '1' if in the first subgroup)
 */
  void GeneSnpPair::CalcAbfsUvlrForSingletons(
    const Grid & grid,
    const map<string,vector<double> > & subgroup2stdsstats,
    const vector<string> & subgroups)
  {
    stringstream config_name;
    vector<int> gamma;
    vector<vector<double> > stdsstats;
    vector<double> l10_abfs;
  
    for(vector<string>::const_iterator it = subgroups.begin();
	it != subgroups.end(); ++it) {
      config_name.str("");
      config_name << it - subgroups.begin() + 1;
      l10_abfs.assign(grid.size(), 0.0);
    
      if(HasResults(*it)) {
	gamma.assign(subgroups.size(), 0);
	stdsstats.clear();
	for(vector<string>::const_iterator it2 = subgroups.begin();
	    it2 != subgroups.end(); ++it2) {
	  if(*it2 == *it) {
	    gamma[it2 - subgroups.begin()] = 1;
	    stdsstats.push_back(subgroup2stdsstats.find(*it)->second);
	  }
	  else {
	    gamma[it2 - subgroups.begin()] = 0;
	    stdsstats.push_back(vector<double>(3, 0.0));
	  }
	}
	for(size_t grid_id = 0; grid_id < grid.size(); ++grid_id)
	  l10_abfs[grid_id] = CalcLog10AbfUvlr(gamma,
					       stdsstats,
					       grid.phi2s[grid_id],
					       grid.oma2s[grid_id]);
      }
      
      unweighted_abfs_.insert(make_pair(config_name.str(), l10_abfs));
      weighted_abfs_.insert(
	make_pair(config_name.str(),
		  log10_weighted_sum(&(l10_abfs[0]), l10_abfs.size())));
    }
  }

/** \brief The configuration is represented only by the index/indices
 *  of the subgroup(s) in which the eQTL is active (eg. '2' or '6-7-11')
 */
  static void prepare_config(const gsl_combination * comb,
			     stringstream & config_name,
			     vector<int> & gamma)
  {
    config_name.str("");
    gamma.assign(comb->n, 0);
  
    config_name << gsl_combination_get(comb, 0) + 1;
    gamma[gsl_combination_get(comb, 0)] = 1;
  
    if(comb->k > 1) {
      for(size_t i = 1; i < comb->k; ++i) {
	config_name << "-" << gsl_combination_get(comb, i) + 1;
	gamma[gsl_combination_get(comb, i)] = 1;
      }
    }
  }
  static void prepare_config(const gsl_combination * comb,
			     stringstream & config_name,
			     gsl_vector * gamma)
  {
    config_name.str("");
    gsl_vector_set_all(gamma, 0.0);
  
    config_name << gsl_combination_get(comb, 0) + 1;
    gsl_vector_set(gamma, gsl_combination_get(comb, 0), 1);
  
    if(comb->k > 1){
      for(size_t i = 1; i < comb->k; ++i){
	config_name << "-" << gsl_combination_get(comb, i) + 1;
	gsl_vector_set(gamma, gsl_combination_get(comb, i), 1);
      }
    }
  }

  void GeneSnpPair::CalcAbfsUvlrForEachConfiguration(
    const Grid & grid,
    const map<string,vector<double> > & subgroup2stdsstats,
    const vector<string> & subgroups)
  {
    gsl_combination * comb;
    stringstream config_name;
    vector<int> gamma(subgroups.size(), 0); // 1,1,0 if S=3 and config="1-2"
    vector<vector<double> > stdsstats;
    vector<double> l10_abfs(grid.size(), NaN);
  
    for(vector<string>::const_iterator it = subgroups.begin();
	it != subgroups.end(); ++it) {
      comb = gsl_combination_calloc(subgroups.size(),
				    it - subgroups.begin() + 1);
      if(comb == NULL) {
	cerr << "ERROR: can't allocate memory for the combination" << endl;
	exit(1);
      }
      while(true) {
	prepare_config(comb, config_name, gamma);
	stdsstats.clear();
	for(vector<string>::const_iterator it2 = subgroups.begin();
	    it2 != subgroups.end(); ++it2) {
	  if(gamma[it2 - subgroups.begin()] == 1 && HasResults(*it2))
	    stdsstats.push_back(subgroup2stdsstats.find(*it2)->second);
	  else {
	    gamma[it2 - subgroups.begin()] = 0;
	    stdsstats.push_back(vector<double>(3, 0.0));
	  }
	}
	l10_abfs.assign(grid.size(), 0.0);
	for(size_t grid_id = 0; grid_id < grid.size(); ++grid_id)
	  l10_abfs[grid_id] = CalcLog10AbfUvlr(gamma,
					       stdsstats,
					       grid.phi2s[grid_id],
					       grid.oma2s[grid_id]);
	unweighted_abfs_.insert(make_pair(config_name.str(), l10_abfs));
	weighted_abfs_.insert(
	  make_pair(config_name.str(), log10_weighted_sum(&(l10_abfs[0]),
							  l10_abfs.size())));
	if(gsl_combination_next(comb) != GSL_SUCCESS)
	  break;
      }
      gsl_combination_free(comb);
    }
  }

  void GeneSnpPair::CalcBMAlite(const vector<string> & subgroups)
  {
    stringstream config_name;
    vector<double> l10_abfs, weights;
  
    for(vector<string>::const_iterator it = subgroups.begin();
	it != subgroups.end(); ++it) {
      config_name.str("");
      config_name << it - subgroups.begin() + 1;
      l10_abfs.push_back(weighted_abfs_[config_name.str()]);
      weights.push_back((1.0 / 2.0) * (1.0 / subgroups.size()));
    }
  
    l10_abfs.push_back(weighted_abfs_["gen"]);
    weights.push_back(1.0 / 2.0);
    weighted_abfs_.insert(
      make_pair("gen-sin", log10_weighted_sum(&(l10_abfs[0]), &(weights[0]),
					      l10_abfs.size())));
  }

  void GeneSnpPair::CalcBMA(const vector<string> & subgroups)
  {
    gsl_combination * comb;
    stringstream config_name;
    vector<int> gamma(subgroups.size(), false); // 1,1,0 if S=3 and config="1-2"
    vector<double> l10_abfs, weights;
  
    for(vector<string>::const_iterator it = subgroups.begin();
	it != subgroups.end(); ++it) {
      comb = gsl_combination_calloc(subgroups.size(),
				    it - subgroups.begin() + 1);
      if(comb == NULL) {
	cerr << "ERROR: can't allocate memory for the combination" << endl;
	exit(1);
      }
      while(true) {
	prepare_config(comb, config_name, gamma);
	l10_abfs.push_back(weighted_abfs_[config_name.str()]);
	weights.push_back((1.0 / (double) subgroups.size())
			  * (1.0 / gsl_sf_choose(subgroups.size(),
						 it - subgroups.begin() + 1)));
	if(gsl_combination_next(comb) != GSL_SUCCESS)
	  break;
      }
      gsl_combination_free(comb);
    }
  
    weighted_abfs_.insert(
      make_pair("all", log10_weighted_sum(&(l10_abfs[0]), &(weights[0]),
					  l10_abfs.size())));
  }

  void GeneSnpPair::CalcAbfsUvlr(const vector<string> & subgroups,
				 const string & whichBfs,
				 const Grid & iGridL,
				 const Grid & iGridS)
  {
    map<string,vector<double> > subgroup2stdsstats;
    StandardizeSstatsAndCorrectSmallSampleSize(subgroup2stdsstats);
    CalcAbfsUvlrForConsistentConfiguration(iGridL, subgroup2stdsstats,
					   subgroups);
    if(whichBfs.find("sin") != string::npos) {
      CalcAbfsUvlrForSingletons(iGridS, subgroup2stdsstats, subgroups);
      CalcBMAlite(subgroups);
    }
    else if(whichBfs.compare("all") == 0) {
      CalcAbfsUvlrForEachConfiguration(iGridS, subgroup2stdsstats, subgroups);
      CalcBMAlite(subgroups);
      CalcBMA(subgroups);
    }
  }

  void GeneSnpPair::CalcAbfsMvlrForConsistentConfiguration(
    const Grid & grid,
    const double & propFitSigma,
    vector<vector<double> > & Y,
    vector<vector<double> > & Xg,
    vector<vector<double> > & Xc)
  {
    MVLR iMvlr;
    iMvlr.set_sigma_option(propFitSigma);
    iMvlr.init(Y, Xg, Xc);
    vector<vector<int> > vvGamma(1, vector<int>(Y.size(), 1));
  
    iMvlr.set_effect_vec(grid.phi2s, grid.oma2s);
    vector<double> l10_abfs = iMvlr.compute_log10_ABF_vec(vvGamma);
    iMvlr.set_effect_vec(grid.phi2s_fix, grid.oma2s_fix);
    vector<double> l10_abfs_fix = iMvlr.compute_log10_ABF_vec(vvGamma);
    iMvlr.set_effect_vec(grid.phi2s_maxh, grid.oma2s_maxh);
    vector<double> l10_abfs_maxh = iMvlr.compute_log10_ABF_vec(vvGamma);
  
    unweighted_abfs_.insert(make_pair("gen", l10_abfs));
    unweighted_abfs_.insert(make_pair("gen-fix", l10_abfs_fix));
    unweighted_abfs_.insert(make_pair("gen-maxh", l10_abfs_maxh));
  
    weighted_abfs_.insert(
      make_pair("gen", log10_weighted_sum(&(l10_abfs[0]),
					  l10_abfs.size())));
    weighted_abfs_.insert(
      make_pair("gen-fix", log10_weighted_sum(&(l10_abfs_fix[0]),
					      l10_abfs_fix.size())));
    weighted_abfs_.insert(
      make_pair("gen-maxh", log10_weighted_sum(&(l10_abfs_maxh[0]),
					       l10_abfs_maxh.size())));
  }

  void GeneSnpPair::CalcAbfsMvlrForSingletons(
    const Grid & grid,
    const double & propFitSigma,
    vector<vector<double> > & Y,
    vector<vector<double> > & Xg,
    vector<vector<double> > & Xc)
  {
    stringstream config_name;
    vector<vector<int> > vvGamma(1, vector<int>());
    vector<double> l10_abfs;
  
    for(size_t s = 0; s < Y.size(); ++s){ // for each subgroup
      config_name.str("");
      config_name << s + 1;
      l10_abfs.assign(grid.size(), 0.0);
    
      vvGamma[0].assign(Y.size(), 0);
      vvGamma[0][s] = 1;
      MVLR iMvlr;
      iMvlr.set_sigma_option(propFitSigma);
      iMvlr.init(Y, Xg, Xc);
      iMvlr.set_effect_vec(grid.phi2s, grid.oma2s);
      l10_abfs = iMvlr.compute_log10_ABF_vec(vvGamma);
    
      unweighted_abfs_.insert(make_pair(config_name.str(), l10_abfs));
      weighted_abfs_.insert(
	make_pair(config_name.str(), log10_weighted_sum(&(l10_abfs[0]),
							l10_abfs.size())));
    }
  }

  void GeneSnpPair::CalcAbfsMvlrForEachConfiguration(
    const Grid & grid,
    const double & propFitSigma,
    vector<vector<double> > & Y,
    vector<vector<double> > & Xg,
    vector<vector<double> > & Xc)
  {
    gsl_combination * comb;
    stringstream config_name;
    vector<int> gamma;
    vector<vector<int> > vvGamma (1, vector<int>());
    vector<double> l10_abfs;
  
    for(size_t k = 1; k <= Y.size(); ++k) {
      comb = gsl_combination_calloc(Y.size(), k);
      if(comb == NULL) {
	cerr << "ERROR: can't allocate memory for the combination" << endl;
	exit(1);
      }
    
      while(true) {
	prepare_config(comb, config_name, gamma);
	vvGamma[0] = gamma;
	MVLR iMvlr;
	iMvlr.set_sigma_option(propFitSigma);
	iMvlr.init(Y, Xg, Xc);
      
	iMvlr.set_effect_vec(grid.phi2s, grid.oma2s);
	l10_abfs = iMvlr.compute_log10_ABF_vec(vvGamma);
	unweighted_abfs_.insert(
	  make_pair(config_name.str(), l10_abfs));
	weighted_abfs_.insert(
	  make_pair(config_name.str(), log10_weighted_sum(&(l10_abfs[0]),
							  l10_abfs.size())));
	if(gsl_combination_next(comb) != GSL_SUCCESS)
	  break;
      }
      gsl_combination_free(comb);
    }
  }

  void GeneSnpPair::CalcAbfsMvlr(const vector<string> & subgroups,
				 const Samples & samples,
				 const Gene & gene,
				 const Snp & snp,
				 const Covariates & covariates,
				 const bool & need_qnorm,
				 const string & whichBfs,
				 const Grid & iGridL,
				 const Grid & iGridS,
				 const double & propFitSigma,
				 const gsl_permutation * perm)
  {
    vector<vector<double> > Y, Xg;
    vector<vector<vector<double> > > Xc;
    vector<string> subgroups_with_data;
    FillStlContainers(samples, gene, snp, covariates, subgroups, true,
		      need_qnorm, perm, Y, Xg, Xc, subgroups_with_data);
  
    CalcAbfsMvlrForConsistentConfiguration(iGridL, propFitSigma, Y, Xg, Xc[0]);
    if(whichBfs.find("sin") != string::npos){ // can also be 'gen-sin' (permutations)
      CalcAbfsMvlrForSingletons(iGridS, propFitSigma, Y, Xg, Xc[0]);
      CalcBMAlite(subgroups);
    }
    else if(whichBfs.compare("all") == 0){
      CalcAbfsMvlrForEachConfiguration(iGridS, propFitSigma, Y, Xg, Xc[0]);
      CalcBMAlite(subgroups);
      CalcBMA(subgroups);
    }
  }

  void GeneSnpPair::CalcBetahatsAndDiagsPerSubgroup(
    const vector<vector<double> > & Y,
    const vector<vector<double> > & Xg,
    const vector<vector<vector<double> > > & Xc,
    const vector<string> & subgroups_with_data,
    const double propFitSigma,
    gsl_matrix * & betas_g_hat,
    gsl_vector * & Sigma_hat_diag,
    gsl_vector * & Vg_diag)
  {
    size_t S = Y.size(); // nb of subgroups
  
    for(size_t s = 0; s < S; ++s){
      size_t N = Y[s].size(); // nb of individuals
      size_t Q = Xc[s].size(), Q1 = Q + 1, Q2 = Q + 2;
      gsl_vector * y = gsl_vector_alloc(N);
      gsl_matrix * X = gsl_matrix_alloc(N, Q2);
      gsl_matrix * Xc_gsl = gsl_matrix_alloc(N, Q1);
      for(size_t i = 0; i < N; ++i) {
	gsl_vector_set(y, i, Y[s][i]);
	gsl_matrix_set(X, i, 0, 1.0); // intercept
	gsl_matrix_set(Xc_gsl, i, 0, 1.0); // intercept
	gsl_matrix_set(X, i, 1, Xg[s][i]);
	for(size_t j = 0; j < Q; ++j){
	  gsl_matrix_set(X, i, j+2, Xc[s][j][i]);
	  gsl_matrix_set(Xc_gsl, i, j+1, Xc[s][j][i]);
	}
      }
    
      // full model
      gsl_matrix * X_ps = gsl_matrix_alloc(Q2, N); // pseudo-inverse of X
      gsl_matrix * U = gsl_matrix_alloc(N, Q2), * V = gsl_matrix_alloc(Q2, Q2);
      gsl_vector * D_diag = gsl_vector_alloc(Q2), * work = gsl_vector_alloc(Q2);
      gsl_matrix_memcpy(U, X);
      gsl_linalg_SV_decomp(U, V, D_diag, work);
      size_t rank_X_full = 0;
      for(size_t i = 0; i < D_diag->size; ++i)
	if(gsl_vector_get(D_diag, i) > GSL_DBL_EPSILON)
	  rank_X_full += 1;
      mygsl_vector_pow(D_diag, -1.0);
      gsl_matrix * D_inv = mygsl_matrix_diagalloc(D_diag, 0.0);
      gsl_matrix * VDinv = gsl_matrix_alloc(Q2, Q2);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, D_inv, 0.0, VDinv);
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, VDinv, U, 0.0, X_ps);
    
      mygsl_vector_pow(D_diag, 2.0);
      gsl_matrix * D_inv2 = mygsl_matrix_diagalloc(D_diag, 0.0);
      gsl_matrix * V_Dinv2 = gsl_matrix_alloc(Q2, Q2);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, D_inv2, 0.0, V_Dinv2);
      gsl_matrix * V_Dinv2_tV = gsl_matrix_alloc(Q2, Q2);
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, V_Dinv2, V, 0.0, V_Dinv2_tV);
    
      gsl_vector * B_hat_full = gsl_vector_alloc(Q2);
      gsl_blas_dgemv(CblasNoTrans, 1.0, X_ps, y, 0.0, B_hat_full);
      gsl_vector * XB_hat_full = gsl_vector_alloc(N);
      gsl_blas_dgemv(CblasNoTrans, 1.0, X, B_hat_full, 0.0, XB_hat_full);
      gsl_vector * E_hat_full = mygsl_vector_alloc(y);
      gsl_vector_sub(E_hat_full, XB_hat_full);
      double rss_full, sigma2_hat_full;
      gsl_blas_ddot(E_hat_full, E_hat_full, &rss_full);
      sigma2_hat_full = rss_full / (double) N;
    
      // record summary stats of full model
      subgroup2pve_[subgroups_with_data[s]] =
	1 - rss_full / gsl_stats_tss(y->data, y->stride, y->size);
      subgroup2sigmahat_[subgroups_with_data[s]] =
	sqrt(rss_full / (double)(N - rank_X_full));
      subgroup2sstats_[subgroups_with_data[s]][0] =
	gsl_vector_get(B_hat_full, 1);
      subgroup2sstats_[subgroups_with_data[s]][1] =
	subgroup2sigmahat_[subgroups_with_data[s]] * sqrt(gsl_matrix_get(V_Dinv2_tV, 1, 1));
      subgroup2sstats_[subgroups_with_data[s]][2] =
	2 * gsl_cdf_tdist_Q(fabs(subgroup2sstats_[subgroups_with_data[s]][0] /
				 subgroup2sstats_[subgroups_with_data[s]][1]),
			    N - rank_X_full);
    
      // null model
      gsl_matrix * Xc_ps = gsl_matrix_alloc(Q1, N);
      mygsl_linalg_pseudoinverse(Xc_gsl, Xc_ps);
      gsl_vector * B_hat_null = gsl_vector_alloc(Q1);
      gsl_blas_dgemv(CblasNoTrans, 1.0, Xc_ps, y, 0.0, B_hat_null);
      gsl_vector * XcB_hat_null = gsl_vector_alloc(N);
      gsl_blas_dgemv(CblasNoTrans, 1.0, Xc_gsl, B_hat_null, 0.0, XcB_hat_null);
      gsl_vector * E_hat_null = mygsl_vector_alloc(y);
      gsl_vector_sub(E_hat_null, XcB_hat_null);
      double sigma2_hat_null;
      gsl_blas_ddot(E_hat_null, E_hat_null, &sigma2_hat_null);
      sigma2_hat_null /= (double) N;
    
      // final estimates
      gsl_matrix_set(betas_g_hat, s, 0, gsl_vector_get(B_hat_full, 1));
      gsl_vector_set(Sigma_hat_diag, s, propFitSigma * sigma2_hat_full
		     + (1-propFitSigma) * sigma2_hat_null);
      gsl_vector_set(Vg_diag, s, gsl_vector_get(Sigma_hat_diag, s)
		     * gsl_matrix_get(V_Dinv2_tV, 1, 1));
    
      gsl_vector_free(y);
      gsl_matrix_free(X);
      gsl_matrix_free(Xc_gsl);
      gsl_matrix_free(X_ps);
      gsl_matrix_free(U);
      gsl_matrix_free(V);
      gsl_vector_free(D_diag);
      gsl_vector_free(work);
      gsl_matrix_free(D_inv);
      gsl_matrix_free(VDinv);
      gsl_matrix_free(D_inv2);
      gsl_matrix_free(V_Dinv2);
      gsl_matrix_free(V_Dinv2_tV);
      gsl_vector_free(B_hat_full);
      gsl_vector_free(XB_hat_full);
      gsl_vector_free(E_hat_full);
      gsl_matrix_free(Xc_ps);
      gsl_vector_free(B_hat_null);
      gsl_vector_free(XcB_hat_null);
      gsl_vector_free(E_hat_null);
    } // end for loop over subgroups
  }

  void GeneSnpPair::FillGslStructuresForPairOfSubgroup(
    const string & subgroup1,
    const string & subgroup2,
    const Samples & samples,
    const Gene & gene,
    const Snp & snp,
    const Covariates & covariates,
    gsl_matrix * & Y_s1s2,
    gsl_matrix * & X_s1s2,
    gsl_matrix * & X_s1,
    gsl_matrix * & X_s2)
  {
    vector<size_t> inds_s1s2, inds_s1, inds_s2;
    inds_s1s2.clear(); inds_s1.clear(); inds_s2.clear();
    samples.GetCommonAndUniqueIndividualsBetweenPairOfSubgroups(
      subgroup1, subgroup2, gene, inds_s1s2, inds_s1, inds_s2);
    if(inds_s1s2.empty()){
      cerr << "ERROR: subgroup " << subgroup1 << " and subgroup "
	   << subgroup2 << " have no individuals in common" << endl;
      exit(1);
    }
  
    size_t Q = covariates.GetNbCovariates(subgroup1);
  
    // fill matrices with individuals common to both subgroup
    Y_s1s2 = gsl_matrix_alloc(inds_s1s2.size(), 2);
    X_s1s2 = gsl_matrix_alloc(inds_s1s2.size(), 1+1+Q);
    for(size_t i = 0; i < inds_s1s2.size(); ++i){
      gsl_matrix_set(Y_s1s2, i, 0,
		     gene.GetExplevel(subgroup1,
				      samples.GetIndexExplevel(inds_s1s2[i],
							       subgroup1)));
      gsl_matrix_set(Y_s1s2, i, 1,
		     gene.GetExplevel(subgroup2,
				      samples.GetIndexExplevel(inds_s1s2[i],
							       subgroup2)));
      gsl_matrix_set(X_s1s2, i, 0, 1.0);
      gsl_matrix_set(X_s1s2, i, 1,
		     snp.GetGenotype(subgroup1,
				     samples.GetIndexGenotype(inds_s1s2[i],
							      subgroup1)));
      size_t j = 0;
      for(map<string,vector<double> >::const_iterator it_covars
	    = covariates.begin(subgroup1);
	  it_covars != covariates.end(subgroup1); ++it_covars){
	gsl_matrix_set(X_s1s2, i, 2 + j,
		       covariates.GetCovariate(subgroup1,
					       it_covars->first, inds_s1s2[i]));
	++j;
      }
    }
  
    // fill design matrix with individuals unique to subgroup s1
    if(inds_s1.size() > 0){
      X_s1 = gsl_matrix_alloc(inds_s1.size(), 1+1+Q);
      for(size_t i = 0; i < inds_s1.size(); ++i){
	gsl_matrix_set(X_s1, i, 0, 1.0);
	gsl_matrix_set(X_s1, i, 1,
		       snp.GetGenotype(subgroup1,
				       samples.GetIndexGenotype(inds_s1[i],
								subgroup1)));
	size_t j = 0;
	for(map<string,vector<double> >::const_iterator it_covars
	      = covariates.begin(subgroup1);
	    it_covars != covariates.end(subgroup1); ++it_covars){
	  gsl_matrix_set(X_s1, i, 2 + j,
			 covariates.GetCovariate(subgroup1,
						 it_covars->first, inds_s1[i]));
	  ++j;
	}
      }
    }
    else
      X_s1 = NULL;
  
    // fill design matrix with individuals unique to subgroup s2
    if(inds_s2.size() > 0){
      X_s2 = gsl_matrix_alloc(inds_s2.size(), 1+1+Q);
      for(size_t i = 0; i < inds_s2.size(); ++i){
	gsl_matrix_set(X_s2, i, 0, 1.0);
	gsl_matrix_set(X_s2, i, 1,
		       snp.GetGenotype(subgroup1,
				       samples.GetIndexGenotype(inds_s2[i],
								subgroup1)));
	size_t j = 0;
	for(map<string,vector<double> >::const_iterator it_covars
	      = covariates.begin(subgroup1);
	    it_covars != covariates.end(subgroup1); ++it_covars){
	  gsl_matrix_set(X_s2, i, 2 + j,
			 covariates.GetCovariate(subgroup1,
						 it_covars->first, inds_s2[i]));
	  ++j;
	}
      }
    }
    else
      X_s2 = NULL;
  }

  void GeneSnpPair::GetMatrixA(
    const gsl_matrix * X_s1s2,
    const gsl_matrix * X_us,
    gsl_matrix * & tXs1s2Xs1s2,
    gsl_matrix * & A_s)
  {
    gsl_matrix * tmp_inv = gsl_matrix_alloc(X_s1s2->size2, X_s1s2->size2);
    gsl_permutation * perm = gsl_permutation_alloc(X_s1s2->size2);
  
    if(X_us == NULL){ // if subgroup has no unique individuals
      mygsl_linalg_pseudoinverse(tXs1s2Xs1s2, tmp_inv);
    }
    else{
      gsl_matrix * tXusXus = gsl_matrix_alloc(X_us->size2, X_us->size2);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X_us, X_us, 0.0, tXusXus);
      gsl_matrix * tmp = gsl_matrix_alloc(tXs1s2Xs1s2->size1, X_us->size2);
      gsl_matrix_memcpy(tmp, tXs1s2Xs1s2);
      gsl_matrix_add(tmp, tXusXus);
      mygsl_linalg_pseudoinverse(tmp, tmp_inv);
      gsl_matrix_free(tXusXus);
      gsl_matrix_free(tmp);
    }
  
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tmp_inv, X_s1s2, 0.0, A_s);
  
    gsl_matrix_free(tmp_inv);
    gsl_permutation_free(perm);
  }

  void GeneSnpPair::GetMatricesA(
    const gsl_matrix * X_s1s2,
    const gsl_matrix * X_s1,
    const gsl_matrix * X_s2,
    gsl_matrix * & tXs1s2Xs1s2,
    gsl_matrix * & A_s1,
    gsl_matrix * & A_s2)
  {
    tXs1s2Xs1s2 = gsl_matrix_alloc(X_s1s2->size2, X_s1s2->size2);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X_s1s2, X_s1s2, 0.0, tXs1s2Xs1s2);
  
    A_s1 = gsl_matrix_alloc(X_s1s2->size2, X_s1s2->size1);
    A_s2 = gsl_matrix_alloc(X_s1s2->size2, X_s1s2->size1);
  
    GetMatrixA(X_s1s2, X_s1, tXs1s2Xs1s2, A_s1);
    GetMatrixA(X_s1s2, X_s2, tXs1s2Xs1s2, A_s2);
  }

/** \brief Estimate the error covariance matrix in a bivariate linear 
 *  regression, under the full and null model
 *  \note Y_s1s2 is N x 2, X_s1s2 is N x Q2 (intercept, genotype and
 *        Q other covariates), tXs1s2Xs1s2 is Q2 x Q2
 */
  void GeneSnpPair::GetErrCovSigmaBtwPairSubgroups(
    const gsl_matrix * Y_s1s2,
    const gsl_matrix * X_s1s2,
    gsl_matrix * & tXs1s2Xs1s2,
    gsl_matrix * & Sigma_s1s2_hat_full,
    gsl_matrix * & Sigma_s1s2_hat_null)
  {
    // estimate Sigma under the full model
    CalcMleErrorCovariance(Y_s1s2, X_s1s2, tXs1s2Xs1s2, Sigma_s1s2_hat_full);
  
    // estimate Sigma under the null model (no genotype)
    gsl_matrix * Xc_s1s2 = gsl_matrix_alloc(X_s1s2->size1, X_s1s2->size2 - 1);
    for(size_t i = 0; i < X_s1s2->size1; ++i){
      gsl_matrix_set(Xc_s1s2, i, 0, 1.0);
      for(size_t j = 2; j < X_s1s2->size2; ++j)
	if(j != 1)
	  gsl_matrix_set(Xc_s1s2, i, j-1, gsl_matrix_get(X_s1s2, i, j));
    }
    CalcMleErrorCovariance(Y_s1s2, Xc_s1s2, NULL, Sigma_s1s2_hat_null);
    gsl_matrix_free(Xc_s1s2);
  }

  void GeneSnpPair::CalcOffDiagCovarsFromPairsOfSubgroups(
    const vector<string> & subgroups,
    const Samples & samples,
    const Gene & gene,
    const Snp & snp,
    const Covariates & covariates,
    const gsl_vector * Sigma_hat_diag,
    const gsl_vector * Vg_diag,
    const double propFitSigma,
    gsl_matrix * & Sigma_hat,
    gsl_matrix * & Vg)
  {
    gsl_vector_view diag = gsl_matrix_diagonal(Sigma_hat);
    gsl_vector_memcpy(&diag.vector, Sigma_hat_diag);
    diag = gsl_matrix_diagonal(Vg);
    gsl_vector_memcpy(&diag.vector, Vg_diag);
  
    gsl_matrix * Sigma_s1s2_hat_full = gsl_matrix_calloc(2, 2),
      * Sigma_s1s2_hat_null = gsl_matrix_calloc(2, 2);
  
    size_t S = subgroups.size();
    gsl_matrix * Y_s1s2, * X_s1s2, * X_s1, * X_s2, * tXs1s2Xs1s2, * A_s1, * A_s2,
      * cov_betahats_s1s2;
    for(size_t s1 = 0; s1 < S-1; ++s1) {
      for(size_t s2 = s1+1; s2 < S; ++s2) {
	// cerr << "s1=" << s1 << " " << "s2=" << s2 << endl;
	FillGslStructuresForPairOfSubgroup(subgroups[s1], subgroups[s2],
					   samples, gene, snp, covariates,
					   Y_s1s2, X_s1s2, X_s1, X_s2);
      
	GetMatricesA(X_s1s2, X_s1, X_s2, tXs1s2Xs1s2, A_s1, A_s2);
	// cerr << "A_s1" << endl;
	// print_matrix(A_s1, 4, 4);
	// cerr << "A_s2" << endl;
	// print_matrix(A_s2, 4, 4);
      
	GetErrCovSigmaBtwPairSubgroups(Y_s1s2, X_s1s2, tXs1s2Xs1s2, Sigma_s1s2_hat_full, Sigma_s1s2_hat_null);
	// cerr << "Sigma_s1s2_hat_full" << endl;
	// print_matrix(Sigma_s1s2_hat_full, 4, 4);
	// cerr << "Sigma_s1s2_hat_null" << endl;
	// print_matrix(Sigma_s1s2_hat_null, 4, 4);
      
	gsl_matrix_set(
	  Sigma_hat, s1, s2,
	  propFitSigma * gsl_matrix_get(Sigma_s1s2_hat_full, 0, 1)
	  + (1-propFitSigma)*gsl_matrix_get(Sigma_s1s2_hat_null, 0, 1));
	gsl_matrix_set(Sigma_hat, s2, s1, gsl_matrix_get(Sigma_hat, s1, s2));
      
	cov_betahats_s1s2 = gsl_matrix_alloc(A_s1->size1, A_s2->size1);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, gsl_matrix_get(Sigma_hat, s1, s2), A_s1, A_s2, 0.0, cov_betahats_s1s2);
	// cerr << "cov_betahats_s1s2" << endl;
	// print_matrix(cov_betahats_s1s2, 4, 4);
	gsl_matrix_set(Vg, s1, s2, gsl_matrix_get(cov_betahats_s1s2, 1, 1));
	gsl_matrix_set(Vg, s2, s1, gsl_matrix_get(Vg, s1, s2));
      
	// cerr << "Sigma_hat" << endl;
	// print_matrix(Sigma_hat, 4, 4);
	// cerr << "Vg" << endl;
	// print_matrix(Vg, 4, 4);
      
	gsl_matrix_free(Y_s1s2);
	gsl_matrix_free(X_s1s2);
	gsl_matrix_free(X_s1);
	gsl_matrix_free(X_s2);
	gsl_matrix_free(tXs1s2Xs1s2);
	gsl_matrix_free(A_s1);
	gsl_matrix_free(A_s2);
	gsl_matrix_free(cov_betahats_s1s2);
      }
    }
  
    gsl_matrix_free(Sigma_s1s2_hat_full);
    gsl_matrix_free(Sigma_s1s2_hat_null);
  }

  void GeneSnpPair::CalcSstatsHybrid(
    const vector<string> & subgroups,
    const Samples & samples,
    const Gene & gene,
    const Snp & snp,
    const Covariates & covariates,
    const bool & need_qnorm,
    const double & propFitSigma,
    const gsl_permutation * perm,
    gsl_matrix * & betas_g_hat,
    gsl_matrix * & Sigma_hat,
    gsl_matrix * & Vg)
  {
    vector<vector<double> > Y, Xg;
    vector<vector<vector<double> > > Xc;
    vector<string> subgroups_with_data;
    FillStlContainers(samples, gene, snp, covariates, subgroups, false,
		      need_qnorm, perm, Y, Xg, Xc, subgroups_with_data);
  
    gsl_vector * Sigma_hat_diag = gsl_vector_alloc(Y.size()),
      * Vg_diag = gsl_vector_alloc(Y.size());
    CalcBetahatsAndDiagsPerSubgroup(Y, Xg, Xc, subgroups_with_data, propFitSigma,
				    betas_g_hat, Sigma_hat_diag, Vg_diag);
  
    CalcOffDiagCovarsFromPairsOfSubgroups(subgroups, samples, gene, snp,
					  covariates, Sigma_hat_diag, Vg_diag,
					  propFitSigma, Sigma_hat, Vg);
  
    gsl_vector_free(Sigma_hat_diag);
    gsl_vector_free(Vg_diag);
  }

  double CalcLog10AbfMvlr(
    const gsl_vector * gamma,
    const gsl_matrix * betas_g_hat,
    const gsl_matrix * Sigma_hat,
    const gsl_matrix * Vg,
    const double phi2,
    const double oma2,
    const bool debug/*=false*/)
  {
    double log10_abf = NaN;
#ifdef DEBUG
    if(debug){fprintf(stderr, "phi2=%.6e oma2=%.6e\n", phi2, oma2);}
#endif
  
    size_t S = gamma->size; // nb of subgroups
    gsl_matrix * Wg = gsl_matrix_alloc(S, S);
    for(size_t i = 0; i < S; ++i){
      for(size_t j = 0; j < S; ++j){
	if(i == j)
	  gsl_matrix_set(Wg, i, j, phi2 + oma2);
	else
	  gsl_matrix_set(Wg, i, j, oma2);
      }
    }
    gsl_matrix * gamma2 = gsl_matrix_alloc(S, S);
    mygsl_linalg_outer(gamma, gamma, gamma2);
    gsl_matrix_mul_elements(Wg, gamma2);
  
    gsl_matrix * Vg_inv = gsl_matrix_alloc(S, S);
    mygsl_linalg_invert(Vg, Vg_inv);
  
    gsl_matrix * bVg = gsl_matrix_alloc(1, S);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, betas_g_hat, Vg_inv, 0.0, bVg);
  
    // update Wg for the ES model
    gsl_matrix * tmp1 = mygsl_matrix_diagalloc(Sigma_hat, 0.0);
    mygsl_matrix_pow(tmp1, 0.5);
    gsl_matrix * tmp2 = gsl_matrix_alloc(S, S);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmp1, Wg, 0.0, tmp2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmp2, tmp1, 0.0, Wg);
#ifdef DEBUG
    if(debug){cerr<<"Wg="<<endl;print_matrix(Wg, 4, 4);}
#endif
  
    gsl_matrix * tmp3 = gsl_matrix_alloc(S, S);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Vg_inv, Wg, 0.0, tmp3);
    gsl_matrix * ivw = gsl_matrix_alloc(S, S);
    gsl_matrix_set_identity(ivw);
    gsl_matrix_add(ivw, tmp3);
    gsl_matrix * ivw_lu = mygsl_matrix_alloc(ivw);
    gsl_permutation * perm = gsl_permutation_alloc(S);
    int signum;
    gsl_linalg_LU_decomp(ivw_lu, perm, &signum);
  
    log10_abf = - 0.5 * gsl_linalg_LU_lndet(ivw_lu);
  
    gsl_matrix * tmp4 = gsl_matrix_alloc(1, S);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, bVg, Wg, 0.0, tmp4);
    gsl_matrix * ivw_inv = gsl_matrix_alloc(S, S);
    mygsl_linalg_invert(ivw, ivw_inv);
    gsl_matrix * tmp5 = gsl_matrix_alloc(1, S);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmp4, ivw_inv, 0.0, tmp5);
    gsl_matrix * tmp6 = gsl_matrix_alloc(1, 1);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tmp5, bVg, 0.0, tmp6);
  
    log10_abf += 0.5 * gsl_matrix_get(tmp6, 0, 0);
  
    log10_abf /= log(10);
#ifdef DEBUG
    if(debug){cerr<<"log10_abf="<<log10_abf<<endl;}
#endif
  
    gsl_matrix_free(Wg);
    gsl_matrix_free(gamma2);
    gsl_matrix_free(Vg_inv);
    gsl_matrix_free(bVg);
    gsl_matrix_free(tmp1);
    gsl_matrix_free(tmp2);
    gsl_matrix_free(tmp3);
    gsl_matrix_free(ivw);
    gsl_matrix_free(ivw_lu);
    gsl_permutation_free(perm);
    gsl_matrix_free(tmp4);
    gsl_matrix_free(ivw_inv);
    gsl_matrix_free(tmp5);
    gsl_matrix_free(tmp6);
  
    return log10_abf;
  }

  void GeneSnpPair::CalcAbfsHybridForConsistentConfiguration(
    const Grid & grid,
    const gsl_matrix * betas_g_hat,
    const gsl_matrix * Sigma_hat,
    const gsl_matrix * Vg)
  {
    size_t S = betas_g_hat->size1; // nb of subgroups
    gsl_vector * gamma = gsl_vector_alloc(S);
    gsl_vector_set_all(gamma, 1);
    vector<double> l10_abfs(grid.size(), NaN),
      l10_abfs_fix(grid.size(), NaN),
      l10_abfs_maxh(grid.size(), NaN);
  
    for(size_t grid_id = 0; grid_id < grid.size(); ++grid_id){
      l10_abfs[grid_id] = CalcLog10AbfMvlr(gamma,
					   betas_g_hat,
					   Sigma_hat,
					   Vg,
					   grid.phi2s[grid_id],
					   grid.oma2s[grid_id]);
      l10_abfs_fix[grid_id] = CalcLog10AbfMvlr(gamma,
					       betas_g_hat,
					       Sigma_hat,
					       Vg,
					       0.0,
					       grid.phi2s[grid_id]
					       + grid.oma2s[grid_id]);
      l10_abfs_maxh[grid_id] = CalcLog10AbfMvlr(gamma,
						betas_g_hat,
						Sigma_hat,
						Vg,
						grid.phi2s[grid_id]
						+ grid.oma2s[grid_id],
						0.0);
    }
  
    unweighted_abfs_.insert(make_pair("gen", l10_abfs));
    unweighted_abfs_.insert(make_pair("gen-fix", l10_abfs_fix));
    unweighted_abfs_.insert(make_pair("gen-maxh", l10_abfs_maxh));
  
    weighted_abfs_.insert(
      make_pair("gen", log10_weighted_sum(&(l10_abfs[0]), l10_abfs.size())));
    weighted_abfs_.insert(
      make_pair("gen-fix",
		log10_weighted_sum(&(l10_abfs_fix[0]), l10_abfs_fix.size())));
    weighted_abfs_.insert(
      make_pair("gen-maxh",
		log10_weighted_sum(&(l10_abfs_maxh[0]), l10_abfs_maxh.size())));
  
    gsl_vector_free(gamma);
  }

  void GeneSnpPair::CalcAbfsHybridForSingletons(
    const Grid & grid,
    const vector<string> & subgroups,
    const gsl_matrix * betas_g_hat,
    const gsl_matrix * Sigma_hat,
    const gsl_matrix * Vg)
  {
    stringstream config_name;
    gsl_vector * gamma = gsl_vector_alloc(subgroups.size());
    vector<double> l10_abfs;
  
    for(vector<string>::const_iterator it_sbgrp = subgroups.begin();
	it_sbgrp != subgroups.end(); ++it_sbgrp){
      config_name.str("");
      config_name << it_sbgrp - subgroups.begin() + 1;
      l10_abfs.assign(grid.size(), 0.0);
      
      gsl_vector_set_all(gamma, 0.0);
      gsl_vector_set(gamma, it_sbgrp - subgroups.begin(), 1.0);
      for(size_t grid_id = 0; grid_id < grid.size(); ++grid_id)
	l10_abfs[grid_id] = CalcLog10AbfMvlr(gamma,
					     betas_g_hat,
					     Sigma_hat,
					     Vg,
					     grid.phi2s[grid_id],
					     grid.oma2s[grid_id]);
      
      unweighted_abfs_.insert(make_pair(config_name.str(), l10_abfs));
      weighted_abfs_.insert(
	make_pair(config_name.str(),
		  log10_weighted_sum(&(l10_abfs[0]), l10_abfs.size())));
    }
  
    gsl_vector_free(gamma);
  }

  void GeneSnpPair::CalcAbfsHybridForEachConfiguration(
    const Grid & grid,
    const vector<string> & subgroups,
    const gsl_matrix * betas_g_hat,
    const gsl_matrix * Sigma_hat,
    const gsl_matrix * Vg)
  {
    gsl_combination * comb;
    stringstream config_name;
    gsl_vector * gamma = gsl_vector_alloc(subgroups.size());
    vector<double> l10_abfs(grid.size(), 0.0);
  
    for(vector<string>::const_iterator it_sbgrp = subgroups.begin();
	it_sbgrp != subgroups.end(); ++it_sbgrp) {
      comb = gsl_combination_calloc(subgroups.size(),
				    it_sbgrp - subgroups.begin() + 1);
      if(comb == NULL) {
	cerr << "ERROR: can't allocate memory for the combination" << endl;
	exit(1);
      }
      while(true) {
	prepare_config(comb, config_name, gamma);
	l10_abfs.assign(grid.size(), 0.0);
	for(size_t grid_id = 0; grid_id < grid.size(); ++grid_id)
	  l10_abfs[grid_id] = CalcLog10AbfMvlr(gamma,
					       betas_g_hat,
					       Sigma_hat,
					       Vg,
					       grid.phi2s[grid_id],
					       grid.oma2s[grid_id]);
	unweighted_abfs_.insert(make_pair(config_name.str(), l10_abfs));
	weighted_abfs_.insert(
	  make_pair(config_name.str(), log10_weighted_sum(&(l10_abfs[0]),
							  l10_abfs.size())));
	if(gsl_combination_next(comb) != GSL_SUCCESS)
	  break;
      }
      gsl_combination_free(comb);
    }
  
    gsl_vector_free(gamma);
  }

  void GeneSnpPair::CalcAbfsHybrid(const vector<string> & subgroups,
				   const Samples & samples,
				   const Gene & gene,
				   const Snp & snp,
				   const Covariates & covariates,
				   const bool & need_qnorm,
				   const string & whichBfs,
				   const Grid & iGridL,
				   const Grid & iGridS,
				   const double & propFitSigma,
				   const gsl_permutation * perm)
  {
    size_t S = subgroups.size();
    gsl_matrix * betas_g_hat = gsl_matrix_alloc(S, 1);
    gsl_matrix * Sigma_hat = gsl_matrix_calloc(S, S),
      * Vg = gsl_matrix_calloc(S, S);
    CalcSstatsHybrid(subgroups, samples, gene, snp, covariates, need_qnorm,
		     propFitSigma, perm, betas_g_hat, Sigma_hat, Vg);
  
    CalcAbfsHybridForConsistentConfiguration(iGridL, betas_g_hat, Sigma_hat, Vg);
    if(whichBfs.find("sin") != string::npos){ // can also be 'gen-sin' (permutations)
      CalcAbfsHybridForSingletons(iGridS, subgroups, betas_g_hat, Sigma_hat, Vg);
      CalcBMAlite(subgroups);
    }
    else if(whichBfs.compare("all") == 0){
      CalcAbfsHybridForEachConfiguration(iGridS, subgroups, betas_g_hat, Sigma_hat,
					 Vg);
      CalcBMAlite(subgroups);
      CalcBMA(subgroups);
    }
    
    gsl_matrix_free(betas_g_hat);
    gsl_matrix_free(Sigma_hat);
    gsl_matrix_free(Vg);
  }

  size_t GeneSnpPair::GetSampleSize(const string & subgroup) const
  {
    return subgroup2samplesize_.find(subgroup)->second;
  }

  double GeneSnpPair::GetPve(const string & subgroup) const
  {
    return subgroup2pve_.find(subgroup)->second;
  }

  double GeneSnpPair::GetSigmahat(const string & subgroup) const
  {
    return subgroup2sigmahat_.find(subgroup)->second;
  }

  double GeneSnpPair::GetBetahatGeno(const string & subgroup) const
  {
    return subgroup2sstats_.find(subgroup)->second[0];
  }

  double GeneSnpPair::GetSebetahatGeno(const string & subgroup) const
  {
    return subgroup2sstats_.find(subgroup)->second[1];
  }

  double GeneSnpPair::GetBetapvalGeno(const string & subgroup) const
  {
    return subgroup2sstats_.find(subgroup)->second[2];
  }

  vector<double>::const_iterator GeneSnpPair::BeginUnweightedAbf(
    const string & abf_type) const
  {
    return unweighted_abfs_.find(abf_type)->second.begin();
  }

  vector<double>::const_iterator GeneSnpPair::EndUnweightedAbf(
    const string & abf_type) const
  {
    return unweighted_abfs_.find(abf_type)->second.end();
  }

  double GeneSnpPair::GetWeightedAbf(const string & abf_type) const
  {
    return weighted_abfs_.find(abf_type)->second;
  }

} //namespace quantgen
