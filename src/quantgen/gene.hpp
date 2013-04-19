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

#include "quantgen/samples.hpp"
#include "quantgen/snp.hpp"
#include "quantgen/gene_snp_pair.hpp"
#include "quantgen/covariates.hpp"
#include "quantgen/grid.hpp"
#include "quantgen/utils_io.hpp"
#include "quantgen/utils_math.hpp"

namespace quantgen {
  
  class GeneSnpPair; // forward declaration
  
  class Gene {
  private:
    string name_;
    string chromosome_;
    size_t start_; // 1-based coordinate
    size_t end_; // idem
    
    map<string,vector<double> > subgroup2explevels_;
    
    vector<Snp*> snps_; // cis snps
    vector<GeneSnpPair> gene_snp_pairs_; // contain results
    
    // test statistic: min P-value over SNPs in each subgroup
    map<string,double> subgroup2permpval_;
    map<string,size_t> subgroup2nbperms_;
    map<string,double> subgroup2trueminpval_;
    
    // test statistic: min P-value over SNPs across subgroups
    size_t nbpermutations_sep_allsbgrps_;
    double pval_perm_sep_allsbgrps_;
    double pval_true_min_allsbgrps_;
    
    // test statistic: ABF (max or avg over SNPs)
    size_t nbpermutations_join_;
    double pval_perm_join_;
    double l10_abf_true_max_;
    double l10_abf_true_avg_;
    
    void FindMinTruePvaluePerSubgroup(const string & subgroup);
    void FindMinTruePvalueAllSubgroups(void);
    void FindMaxTrueL10Abf(const string & whichPermBf);
    void AvgTrueL10Abfs(const string & whichPermBf);
    double CalcPermutationPvalue(const size_t & nbperms_total,
				 const size_t & nbperms_sofar,
				 const double & nbperms_more_extreme,
				 const size_t & trick_cutoff,
				 const gsl_rng * rng) const;
    
  public:
    Gene(void);
    Gene(const string & name, const string & chromosome,
	 const string & start, const string & end);
    string GetName(void) const { return name_; };
    string GetChromosome(void) const { return chromosome_; };
    size_t GetStart(void) const { return start_; };
    size_t GetEnd(void) const { return end_; };
    size_t GetNbSubgroups(void) const { return subgroup2explevels_.size(); };
    size_t GetNbGeneSnpPairs(void) const { return gene_snp_pairs_.size(); };
    bool HasExplevelsInAtLeastOneSubgroup(void) const;
    bool HasExplevels(const string & subgroup) const;
    bool HasExplevelsInAllSubgroups(const vector<string> & subgroups) const;
    bool HasAtLeastOneCisSnpInAtLeastOneSubgroup(void) const;
    bool HasAtLeastOneCisSnp(const string & subgroup) const;
    void AddSubgroup(const string & subgroup,
		     const vector<string>::const_iterator & begin,
		     const vector<string>::const_iterator & end);
    void Show(ostream & os);
    void SetCisSnps(const map<string, vector<Snp*> > & mChr2VecPtSnps,
		    const string & anchor, const size_t & radius);
    double GetExplevel(const string & subgroup, const size_t & idx) const;
    void TestForAssociations(const vector<string> & subgroups,
			     const Samples & samples,
			     const string & type_analysis,
			     const bool & need_qnorm,
			     const Covariates & covariates,
			     const Grid & iGridL,
			     const Grid & iGridS,
			     const string & whichBfs,
			     const string & covErrors,
			     const float & propFitSigma,
			     const int & verbose);
    vector<GeneSnpPair>::const_iterator BeginPair(void) const;
    vector<GeneSnpPair>::const_iterator EndPair(void) const;
    void MakePermutationsSepPerSubgroup(
      const string & subgroup,
      const Samples & samples,
      const bool & need_qnorm,
      const Covariates & covariates,
      const size_t & nb_permutations,
      const int & trick,
      const size_t & trick_cutoff,
      const int & nb_threads,
      const gsl_rng * rngPerm,
      const gsl_rng * rngTrick);
    size_t GetNbGeneSnpPairs(const string & subgroup) const;
    double GetPermutationPvalueSep(const string & subgroup) const;
    size_t GetNbPermutationsSep(const string & subgroup) const;
    double GetTrueMinPval(const string & subgroup) const;
    void MakePermutationsSepAllSubgroups(
      const vector<string> & subgroups,
      const Samples & samples,
      const bool & need_qnorm,
      const Covariates & covariates,
      const size_t & nb_permutations,
      const int & trick,
      const size_t & trick_cutoff,
      const int & nb_threads,
      const gsl_rng * rngPerm,
      const gsl_rng * rngTrick);
    double GetPermutationPvalueSep(void) const { return pval_perm_sep_allsbgrps_; };
    size_t GetNbPermutationsSep(void) const { return nbpermutations_sep_allsbgrps_; };
    double GetTrueMinPval(void) const { return pval_true_min_allsbgrps_; };
    void MakePermutationsJoin(const vector<string> & subgroups,
			      const Samples & samples,
			      const bool & need_qnorm,
			      const Covariates & covariates,
			      const Grid & iGridL,
			      const Grid & iGridS,
			      const string & covErrors,
			      const float & propFitSigma,
			      const size_t & nbPerms,
			      const int & trick,
			      const size_t & trick_cutoff,
			      const string & whichPermBf,
			      const bool & useMaxBfOverSnps,
			      const int & nb_threads,
			      const gsl_rng * rngPerm,
			      const gsl_rng * rngTrick);
    double GetPermutationPvalueJoin(void) const { return pval_perm_join_; };
    size_t GetNbPermutationsJoin(void) const { return nbpermutations_join_; };
    double GetTrueL10Abf(const bool & use_max_bf) const;
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
