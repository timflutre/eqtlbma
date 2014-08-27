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

#include <cstdlib>

#include <vector>
#include <map>
#include <string>
#include <sstream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "utils/utils_math.hpp"

#include "quantgen/gene.hpp"
#include "quantgen/snp.hpp"
#include "quantgen/samples.hpp"
#include "quantgen/grid.hpp"
#include "quantgen/MVLR.hpp"

#include "irls/IRLS.h"

namespace quantgen {
  
  // forward declaration, e.g. in GetSstatsOneSbgrp()
  class Gene;
  class Snp;
  class Covariates;
  
  class GeneSnpPair {
  private:
    std::string gene_name_;
    std::string snp_name_;
    std::string error_model_; // uvlr or mvlr or hybrid (optional)
    
    std::map<std::string,size_t> subgroup2samplesize_;
    std::map<std::string,size_t> subgroup2nbcovariates_;
    std::map<std::string,double> subgroup2pve_;
    std::map<std::string,double> subgroup2sigmahat_;
    std::map<std::string,std::vector<double> > subgroup2sstats_; // 0:betahat ; 1:sebetahat ; 2:betapval
  
    // raw ABFs
    // keys: 'gen', 'gen-fix', 'gen-maxh', '1-2-3', '1', '2', etc
    std::map<std::string,std::vector<double> > unweighted_abfs_;
  
    // averaged ABFs (over a grid for all, over configurations for some)
    // keys: same as above
    std::map<std::string,double> weighted_abfs_;
  
    void FillStlContainers(
      const Samples & samples,
      const Gene & gene,
      const Snp & snp,
      const Covariates & covariates,
      const std::vector<std::string> & subgroups,
      const bool & same_individuals,
      const bool & needQnorm,
      const gsl_permutation * perm,
      std::vector<std::vector<double> > & Y,
      std::vector<std::vector<double> > & Xg,
      std::vector<std::vector<std::vector<double> > > & Xc,
      std::vector<std::string> & subgroups_with_data);
    void FillGslStructuresForPairOfSubgroup(
      const std::string & subgroup1,
      const std::string & subgroup2,
      const Samples & samples,
      const Gene & gene,
      const Snp & snp,
      const Covariates & covariates,
      gsl_matrix * & Y_s1s2,
      gsl_matrix * & X_s1s2,
      gsl_matrix * & X_s1,
      gsl_matrix * & X_s2);
    void CalcBetahatsAndDiagsPerSubgroup(
      const std::vector<std::vector<double> > & Y,
      const std::vector<std::vector<double> > & Xg,
      const std::vector<std::vector<std::vector<double> > > & Xc,
      const std::vector<std::string> & subgroups_with_data,
      const double propFitSigma,
      gsl_matrix * & betas_g_hat,
      gsl_vector * & Sigma_hat_diag,
      gsl_vector * & Vg_diag);
    void CalcOffDiagCovarsFromPairsOfSubgroups(
      const std::vector<std::string> & subgroups,
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
    GeneSnpPair(void);
    GeneSnpPair(const std::string & gene_name, const std::string & snp_name);
    GeneSnpPair(const std::string & gene_name, const std::string & snp_name,
		const std::string & error_model);
    void SetGeneName(const std::string & gene_name) { gene_name_ = gene_name; };
    void SetSnpName(const std::string & snp_name) { snp_name_ = snp_name; };
    void SetErrorModel(const std::string & error_model) { error_model_ = error_model; };
    std::string GetGeneName(void) const { return gene_name_; };
    std::string GetSnpName(void) const { return snp_name_; };
    std::string GetErrorModel(void) { return error_model_; };
    bool HasResults(const std::string & subgroup) const;
    void CalcSstatsOneSbgrp(
      const Samples & samples,
      const Gene & gene,
      const Snp & snp,
      const Covariates & covariates,
      const std::string & subgroup,
      const std::string & likelihood,
      const bool & needQnorm,
      const gsl_permutation * perm);
    void SetSstats(const std::string & subgroup,
		   const size_t & n,
		   const double & sigmahat,
		   const double & betahatgeno,
		   const double & sebetahatgeno);
    void StandardizeSstatsAndCorrectSmallSampleSize(
      std::map<std::string,std::vector<double> > & subgroup2stdsstats);
    void CalcAbfsUvlrForConsistentConfiguration(
      const Grid & grid,
      const std::map<std::string,std::vector<double> > & subgroup2stdsstats,
      const std::vector<std::string> & subgroups);
    void CalcAbfsUvlrForSingletons(
      const Grid & grid,
      const std::map<std::string,std::vector<double> > & subgroup2stdsstats,
      const std::vector<std::string> & subgroups);
    void CalcAbfsUvlrForEachConfiguration(
      const Grid & grid,
      const std::map<std::string,std::vector<double> > & subgroup2stdsstats,
      const std::vector<std::string> & subgroups);
    void CalcBMAlite(const std::vector<std::string> & subgroups);
    void CalcBMA(const std::vector<std::string> & subgroups);
    void CalcAbfsUvlr(
      const std::vector<std::string> & subgroups,
      const std::string & whichBfs,
      const Grid & iGridL,
      const Grid & iGridS);
    void CalcAbfsMvlrForConsistentConfiguration(
      const Grid & iGridL,
      const double & propFitSigma,
      std::vector<std::vector<double> > & Y,
      std::vector<std::vector<double> > & Xg,
      std::vector<std::vector<double> > & Xc);
    void CalcAbfsMvlrForSingletons(
      const Grid & grid,
      const double & propFitSigma,
      std::vector<std::vector<double> > & Y,
      std::vector<std::vector<double> > & Xg,
      std::vector<std::vector<double> > & Xc);
    void CalcAbfsMvlrForEachConfiguration(
      const Grid & grid,
      const double & propFitSigma,
      std::vector<std::vector<double> > & Y,
      std::vector<std::vector<double> > & Xg,
      std::vector<std::vector<double> > & Xc);
    void CalcAbfsMvlr(
      const std::vector<std::string> & subgroups,
      const Samples & samples,
      const Gene & gene,
      const Snp & snp,
      const Covariates & covariates,
      const bool & needQnorm,
      const std::string & whichBfs,
      const Grid & iGridL,
      const Grid & iGridS,
      const double & propFitSigma,
      const gsl_permutation * perm);
    void CalcSstatsHybrid(
      const std::vector<std::string> & subgroups,
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
      const std::vector<std::string> & subgroups,
      const gsl_matrix * betas_g_hat,
      const gsl_matrix * Sigma_hat,
      const gsl_matrix * Vg);
    void CalcAbfsHybridForEachConfiguration(
      const Grid & grid,
      const std::vector<std::string> & subgroups,
      const gsl_matrix * betas_g_hat,
      const gsl_matrix * Sigma_hat,
      const gsl_matrix * Vg);
    void CalcAbfsHybrid(
      const std::vector<std::string> & subgroups,
      const Samples & samples,
      const Gene & gene,
      const Snp & snp,
      const Covariates & covariates,
      const bool & needQnorm,
      const std::string & whichBfs,
      const Grid & iGridL,
      const Grid & iGridS,
      const double & propFitSigma,
      const gsl_permutation * perm);
    size_t GetNbSubgroups(void) const { return subgroup2samplesize_.size(); };
    size_t GetSampleSize(const std::string & subgroup) const;
    double GetPve(const std::string & subgroup) const;
    double GetSigmahat(const std::string & subgroup) const;
    double GetBetahatGeno(const std::string & subgroup) const;
    double GetSebetahatGeno(const std::string & subgroup) const;
    double GetBetapvalGeno(const std::string & subgroup) const;
    std::vector<double>::const_iterator BeginUnweightedAbf(
      const std::string & abf_type) const;
    std::vector<double>::const_iterator EndUnweightedAbf(
      const std::string & abf_type) const;
    double GetWeightedAbf(const std::string & abf_type) const;
  };

  double CalcLog10AbfUvlr(
    const std::vector<int> & gamma,
    const std::vector<std::vector<double> > & stdsstats,
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

  bool operator==(const GeneSnpPair& lhs, const GeneSnpPair& rhs);
  bool operator!=(const GeneSnpPair& lhs, const GeneSnpPair& rhs);

} // namespace quantgen

#endif // QUANTGEN_GENE_SNP_PAIR_HPP
