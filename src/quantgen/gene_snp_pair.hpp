/** \file gene_snp_pair.hpp
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

#ifndef QUANTGEN_GENE_SNP_PAIR_HPP
#define QUANTGEN_GENE_SNP_PAIR_HPP

#include <cstdio>

#include <vector>
#include <map>
#include <string>
#include <numeric>

#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "quantgen/utils_math.hpp"
#include "quantgen/gene.hpp"
#include "quantgen/snp.hpp"
#include "quantgen/samples.hpp"
#include "quantgen/grid.hpp"
#include "quantgen/MVLR.hpp"

namespace quantgen {
  
  // forward declaration, e.g. in GetSstatsOneSbgrp()
  class Gene;
  class Snp;
  class Covariates;
  
  class GeneSnpPair {
  private:
    string gene_name_;
    string snp_name_;
    string analysis_type_; // uvlr or mvlr
    
    map<string,size_t> subgroup2samplesize_;
    map<string,size_t> subgroup2nbcovariates_;
    map<string,double> subgroup2pve_;
    map<string,double> subgroup2sigmahat_;
    map<string,vector<double> > subgroup2sstats_; // 0:betahat ; 1:sebetahat ; 2:betapval
  
    // raw ABFs
    // keys: 'gen', 'gen-fix', 'gen-maxh', '1-2-3', '1', '2', etc
    map<string,vector<double> > unweighted_abfs_;
  
    // averaged ABFs (over a grid for all, over configurations for some)
    // keys: same as above
    map<string,double> weighted_abfs_;
  
    void FillStlContainers(
      const Samples & samples,
      const Gene & gene,
      const Snp & snp,
      const Covariates & covariates,
      const vector<string> & subgroups,
      const bool & same_individuals,
      const bool & needQnorm,
      const gsl_permutation * perm,
      vector<vector<double> > & Y,
      vector<vector<double> > & Xg,
      vector<vector<vector<double> > > & Xc,
      vector<string> & subgroups_with_data);
    void FillGslStructuresForPairOfSubgroup(
      const string & subgroup1,
      const string & subgroup2,
      const Samples & samples,
      const Gene & gene,
      const Snp & snp,
      const Covariates & covariates,
      gsl_matrix * & Y_s1s2,
      gsl_matrix * & X_s1s2,
      gsl_matrix * & X_s1,
      gsl_matrix * & X_s2);
    void CalcBetahatsAndDiagsPerSubgroup(
      const vector<vector<double> > & Y,
      const vector<vector<double> > & Xg,
      const vector<vector<vector<double> > > & Xc,
      const vector<string> & subgroups_with_data,
      const double propFitSigma,
      gsl_matrix * & betas_g_hat,
      gsl_vector * & Sigma_hat_diag,
      gsl_vector * & Vg_diag);
    void CalcOffDiagCovarsFromPairsOfSubgroups(
      const vector<string> & subgroups,
      const Samples & samples,
      const Gene & gene,
      const Snp & snp,
      const Covariates & covariates,
      const gsl_vector * Sigma_hat_diag,
      const gsl_vector * Vg_diag,
      const double propFitSigma,
      gsl_matrix * & Sigma_hat,
      gsl_matrix * & Vg);
    void GetMatrixA(
      const gsl_matrix * X_s1s2,
      const gsl_matrix * X_us,
      gsl_matrix * & tXs1s2Xs1s2,
      gsl_matrix * & A_s);
    void GetMatricesA(
      const gsl_matrix * X_s1s2,
      const gsl_matrix * X_s1,
      const gsl_matrix * X_s2,
      gsl_matrix * & tXs1s2Xs1s2,
      gsl_matrix * & A_s1,
      gsl_matrix * & A_s2);
    void GetErrCovSigmaBtwPairSubgroups(
      const gsl_matrix * Y_s1s2,
      const gsl_matrix * X_s1s2,
      gsl_matrix * & tXs1s2Xs1s2,
      gsl_matrix * & Sigma_s1s2_hat_full,
      gsl_matrix * & Sigma_s1s2_hat_null);
  
  public:
    GeneSnpPair();
    GeneSnpPair(const string & gene_name, const string & snp_name);
    string GetGeneName(void) const { return gene_name_; };
    string GetSnpName(void) const { return snp_name_; };
    bool HasResults(const string & subgroup) const;
    void CalcSstatsOneSbgrp(
      const Samples & samples,
      const Gene & gene,
      const Snp & snp,
      const Covariates & covariates,
      const string & subgroup,
      const bool & needQnorm,
      const gsl_permutation * perm);
    void StandardizeSstatsAndCorrectSmallSampleSize(
      map<string,vector<double> > & subgroup2stdsstats);
    void CalcAbfsUvlrForConsistentConfiguration(
      const Grid & grid,
      const map<string,vector<double> > & subgroup2stdsstats,
      const vector<string> & subgroups);
    void CalcAbfsUvlrForSingletons(
      const Grid & grid,
      const map<string,vector<double> > & subgroup2stdsstats,
      const vector<string> & subgroups);
    void CalcAbfsUvlrForEachConfiguration(
      const Grid & grid,
      const map<string,vector<double> > & subgroup2stdsstats,
      const vector<string> & subgroups);
    void CalcBMAlite(const vector<string> & subgroups);
    void CalcBMA(const vector<string> & subgroups);
    void CalcAbfsUvlr(
      const vector<string> & subgroups,
      const string & whichBfs,
      const Grid & iGridL,
      const Grid & iGridS);
    void CalcAbfsMvlrForConsistentConfiguration(
      const Grid & iGridL,
      const double & propFitSigma,
      vector<vector<double> > & Y,
      vector<vector<double> > & Xg,
      vector<vector<double> > & Xc);
    void CalcAbfsMvlrForSingletons(
      const Grid & grid,
      const double & propFitSigma,
      vector<vector<double> > & Y,
      vector<vector<double> > & Xg,
      vector<vector<double> > & Xc);
    void CalcAbfsMvlrForEachConfiguration(
      const Grid & grid,
      const double & propFitSigma,
      vector<vector<double> > & Y,
      vector<vector<double> > & Xg,
      vector<vector<double> > & Xc);
    void CalcAbfsMvlr(
      const vector<string> & subgroups,
      const Samples & samples,
      const Gene & gene,
      const Snp & snp,
      const Covariates & covariates,
      const bool & needQnorm,
      const string & whichBfs,
      const Grid & iGridL,
      const Grid & iGridS,
      const double & propFitSigma,
      const gsl_permutation * perm);
    void CalcSstatsHybrid(
      const vector<string> & subgroups,
      const Samples & samples,
      const Gene & gene,
      const Snp & snp,
      const Covariates & covariates,
      const bool & needQnorm,
      const double & propFitSigma,
      const gsl_permutation * perm,
      gsl_matrix * & betas_g_hat,
      gsl_matrix * & Sigma_hat,
      gsl_matrix * & Vg);
    void CalcAbfsHybridForConsistentConfiguration(
      const Grid & grid,
      const gsl_matrix * betas_g_hat,
      const gsl_matrix * Sigma_hat,
      const gsl_matrix * Vg);
    void CalcAbfsHybridForSingletons(
      const Grid & grid,
      const vector<string> & subgroups,
      const gsl_matrix * betas_g_hat,
      const gsl_matrix * Sigma_hat,
      const gsl_matrix * Vg);
    void CalcAbfsHybridForEachConfiguration(
      const Grid & grid,
      const vector<string> & subgroups,
      const gsl_matrix * betas_g_hat,
      const gsl_matrix * Sigma_hat,
      const gsl_matrix * Vg);
    void CalcAbfsHybrid(
      const vector<string> & subgroups,
      const Samples & samples,
      const Gene & gene,
      const Snp & snp,
      const Covariates & covariates,
      const bool & needQnorm,
      const string & whichBfs,
      const Grid & iGridL,
      const Grid & iGridS,
      const double & propFitSigma,
      const gsl_permutation * perm);
    size_t GetNbSubgroups(void) const { return subgroup2samplesize_.size(); };
    size_t GetSampleSize(const string & subgroup) const;
    double GetPve(const string & subgroup) const;
    double GetSigmahat(const string & subgroup) const;
    double GetBetahatGeno(const string & subgroup) const;
    double GetSebetahatGeno(const string & subgroup) const;
    double GetBetapvalGeno(const string & subgroup) const;
    vector<double>::const_iterator BeginUnweightedAbf(
      const string & abf_type) const;
    vector<double>::const_iterator EndUnweightedAbf(
      const string & abf_type) const;
    double GetWeightedAbf(const string & abf_type) const;
  };

  double CalcLog10AbfUvlr(
    const vector<int> & gamma,
    const vector<vector<double> > & stdsstats,
    const double phi2,
    const double oma2);

  double CalcLog10AbfMvlr(
    const gsl_vector * gamma,
    const gsl_matrix * betas_g_hat,
    const gsl_matrix * Sigma_hat,
    const gsl_matrix * Vg,
    const double phi2,
    const double oma2,
    const bool debug=false);

} // namespace quantgen

#endif // QUANTGEN_GENE_SNP_PAIR_HPP
