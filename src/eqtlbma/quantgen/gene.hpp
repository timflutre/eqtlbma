/** \file gene.hpp
 *
 *  `Gene' is a class 
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

#ifndef QUANTGEN_GENE_HPP
#define QUANTGEN_GENE_HPP

#include <vector>
#include <map>
#include <string>
#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <omp.h>

#include "utils/utils_io.hpp"
#include "utils/utils_math.hpp"

#include "quantgen/samples.hpp"
#include "quantgen/snp.hpp"
#include "quantgen/gene_snp_pair.hpp"
#include "quantgen/covariates.hpp"
#include "quantgen/grid.hpp"

namespace quantgen {
  
  class GeneSnpPair; // forward declaration
  
  class Gene {
  private:
    std::string name_;
    std::string chromosome_;
    size_t start_; // 1-based coordinate
    size_t end_; // idem
    
    std::map<std::string,std::vector<double> > subgroup2explevels_;
    
    std::vector<Snp*> snps_; // cis snps
    std::vector<GeneSnpPair> gene_snp_pairs_; // contain results
    
    // test statistic: min P-value over SNPs in each subgroup
    std::map<std::string,double> subgroup2permpval_;
    std::map<std::string,size_t> subgroup2nbperms_;
    std::map<std::string,double> subgroup2trueminpval_;
    
    // test statistic: min P-value over SNPs across subgroups
    size_t nbpermutations_sep_allsbgrps_;
    double pval_perm_sep_allsbgrps_;
    double pval_true_min_allsbgrps_;
    
    // test statistic: ABF (max or avg over SNPs)
    size_t nbpermutations_join_;
    double pval_perm_join_;
    double l10_abf_true_max_;
    double l10_abf_true_avg_;
    double l10_abf_perm_med_;
    
    void FindMinTruePvaluePerSubgroup(const std::string & subgroup);
    void FindMinTruePvalueAllSubgroups(void);
    void FindMaxTrueL10Abf(const std::string & whichPermBf);
    void AvgTrueL10Abfs(const std::string & whichPermBf);
    double CalcPermutationPvalue(const size_t & nbperms_total,
				 const size_t & nbperms_sofar,
				 const double & nbperms_more_extreme,
				 const size_t & trick_cutoff,
				 const gsl_rng * rng) const;
    
  public:
    Gene(void);
    Gene(const std::string & name, const std::string & chromosome,
	 const std::string & start, const std::string & end);
    std::string GetName(void) const { return name_; };
    std::string GetChromosome(void) const { return chromosome_; };
    size_t GetStart(void) const { return start_; };
    size_t GetEnd(void) const { return end_; };
    std::string GetRegionInTabixFormat(const std::string & anchor,
				       const size_t & radius) const;
    size_t GetNbSubgroups(void) const { return subgroup2explevels_.size(); };
    size_t GetNbGeneSnpPairs(void) const { return gene_snp_pairs_.size(); };
    bool HasExplevelsInAtLeastOneSubgroup(void) const;
    bool HasExplevels(const std::string & subgroup) const;
    bool HasExplevelsInAllSubgroups(const std::vector<std::string> & subgroups) const;
    bool HasAtLeastOneCisSnpInAtLeastOneSubgroup(void) const;
    bool HasAtLeastOneCisSnp(const std::string & subgroup) const;
    void AddSubgroup(const std::string & subgroup,
		     const std::vector<std::string>::const_iterator & begin,
		     const std::vector<std::string>::const_iterator & end);
    void Show(std::ostream & os);
    void SetCisSnps(const std::map<std::string, std::vector<Snp*> > & mChr2VecPtSnps,
		    const std::string & anchor, const size_t & radius);
    double GetExplevel(const std::string & subgroup, const size_t & idx) const;
    void TestForAssociations(const std::vector<std::string> & subgroups,
			     const Samples & samples,
			     const std::string & likelihood,
			     const std::string & type_analysis,
			     const bool & need_qnorm,
			     const Covariates & covariates,
			     const Grid & iGridL,
			     const Grid & iGridS,
			     const std::string & whichBfs,
			     const std::string & covErrors,
			     const float & propFitSigma,
			     const int & verbose);
    std::vector<GeneSnpPair>::const_iterator BeginPair(void) const;
    std::vector<GeneSnpPair>::const_iterator EndPair(void) const;
    void MakePermutationsSepPerSubgroup(
      const std::string & subgroup,
      const Samples & samples,
      const std::string & likelihood,
      const bool & need_qnorm,
      const Covariates & covariates,
      const size_t & nb_permutations,
      const int & trick,
      const size_t & trick_cutoff,
      const gsl_rng * rngPerm,
      const gsl_rng * rngTrick);
    size_t GetNbGeneSnpPairs(const std::string & subgroup) const;
    double GetPermutationPvalueSep(const std::string & subgroup) const;
    size_t GetNbPermutationsSep(const std::string & subgroup) const;
    double GetTrueMinPval(const std::string & subgroup) const;
    void MakePermutationsSepAllSubgroups(
      const std::vector<std::string> & subgroups,
      const Samples & samples,
      const std::string & likelihood,
      const bool & need_qnorm,
      const Covariates & covariates,
      const size_t & nb_permutations,
      const int & trick,
      const size_t & trick_cutoff,
      const gsl_rng * rngPerm,
      const gsl_rng * rngTrick);
    double GetPermutationPvalueSep(void) const { return pval_perm_sep_allsbgrps_; };
    size_t GetNbPermutationsSep(void) const { return nbpermutations_sep_allsbgrps_; };
    double GetTrueMinPval(void) const { return pval_true_min_allsbgrps_; };
    void MakePermutationsJoin(const std::vector<std::string> & subgroups,
			      const Samples & samples,
			      const std::string & likelihood,
			      const bool & need_qnorm,
			      const Covariates & covariates,
			      const Grid & iGridL,
			      const Grid & iGridS,
			      const std::string & covErrors,
			      const float & propFitSigma,
			      const size_t & nbPerms,
			      const int & trick,
			      const size_t & trick_cutoff,
			      const std::string & whichPermBf,
			      const bool & useMaxBfOverSnps,
			      const gsl_rng * rngPerm,
			      const gsl_rng * rngTrick);
    double GetPermutationPvalueJoin(void) const { return pval_perm_join_; };
    size_t GetNbPermutationsJoin(void) const { return nbpermutations_join_; };
    double GetTrueL10Abf(const bool & use_max_bf) const;
    double GetMedianPermL10Abf(void) const { return l10_abf_perm_med_; };
  };
  
  bool operator==(const Gene& lhs, const Gene& rhs);
  bool operator!=(const Gene& lhs, const Gene& rhs);
  bool operator< (const Gene& lhs, const Gene& rhs);
  bool operator> (const Gene& lhs, const Gene& rhs);
  bool operator<=(const Gene& lhs, const Gene& rhs);
  bool operator>=(const Gene& lhs, const Gene& rhs);
  
  // 'lt' means '<', not '<=' (that would be 'le')
  bool pt_gene_lt_pt_gene(const Gene* pt_lhs, const Gene* pt_rhs);
  
} //namespace quantgen

#endif // QUANTGEN_GENE_HPP
