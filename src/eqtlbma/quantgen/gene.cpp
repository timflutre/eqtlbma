/** \file gene.cpp
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

#include <algorithm>

#include "quantgen/gene.hpp"

using namespace std;

using namespace utils;

namespace quantgen {
  
  Gene::Gene(void)
  {
  }
  
  Gene::Gene(const string & name, const string & chromosome,
	     const string & start, const string & end)
  {
    name_ = name;
    chromosome_ = chromosome;
    start_ = atol(start.c_str());
    end_ = atol(end.c_str());
  }
  
  bool operator==(const Gene& lhs, const Gene& rhs)
  {
    return(lhs.GetName().compare(rhs.GetName()) == 0
	   && lhs.GetChromosome().compare(rhs.GetChromosome()) == 0
	   && lhs.GetStart() == rhs.GetStart()
	   && lhs.GetEnd() == rhs.GetEnd()
	   && lhs.GetNbSubgroups() == rhs.GetNbSubgroups());
  }
  
  bool operator!=(const Gene& lhs, const Gene& rhs)
  {
    return !operator==(lhs,rhs);
  }

  bool operator<(const Gene& lhs, const Gene& rhs)
  {
    if(lhs.GetChromosome().compare(rhs.GetChromosome()) != 0){
      fprintf(stderr, "ERROR: %s and %s are on different chromosomes, thus they can't be sorted\n",
	      lhs.GetName().c_str(), rhs.GetName().c_str());
      exit(1);
    }
    return((lhs.GetStart() < rhs.GetStart()) ||
	   (lhs.GetStart() == rhs.GetStart() && 
	    lhs.GetEnd() < rhs.GetEnd() ));
  }

  bool operator>(const Gene& lhs, const Gene& rhs)
  {
    return operator<(rhs,lhs);
  }

  bool operator<=(const Gene& lhs, const Gene& rhs)
  {
    return !operator>(lhs,rhs);
  }

  bool operator>=(const Gene& lhs, const Gene& rhs)
  {
    return !operator<(lhs,rhs);
  }

  bool pt_gene_lt_pt_gene(const Gene* pt_lhs, const Gene* pt_rhs)
  {
    return *pt_lhs < *pt_rhs;
  }

  string Gene::GetRegionInTabixFormat(const string & anchor,
				      const size_t & radius) const
  {
    stringstream region;
    region << chromosome_ << ":";
    if(anchor.compare("TSS+TES") == 0)
      region << start_ - radius << "-" << end_ + radius;
    else if(anchor.compare("TSS") == 0)
      region << start_ - radius << "-" << start_ + radius;
    return region.str();
  }

  void Gene::AddSubgroup(const string & subgroup,
			 const vector<string>::const_iterator & begin,
			 const vector<string>::const_iterator & end)
  {
    vector<double> explevels(end - begin, NaN);
    size_t i = 0;
    bool any_missing_value = false;
    for(vector<string>::const_iterator it = begin; it != end; ++it){
      if(it->compare("NA") == 0 || it->compare("na") == 0
	 || it->compare("NaN") == 0 || it->compare("nan") == 0){
	cerr << "WARNING: skip gene " << name_ << " in subgroup " << subgroup
	     << " because of missing value" << endl;
	any_missing_value = true;
	break;
      }
      explevels[i] = atof(it->c_str());
      ++i;
    }
    if(! any_missing_value)
      subgroup2explevels_.insert(make_pair(subgroup, explevels));
  }

  void Gene::Show(ostream & os)
  {
    os << name_ << " " << chromosome_ << " " << start_ << " " << end_ << endl
       << GetNbSubgroups() << " subgroups" << endl;
    for(map<string,vector<double> >::const_iterator it =
	  subgroup2explevels_.begin(); it != subgroup2explevels_.end(); ++it)
      os << it->first << ": " << it->second.size() << " samples" << endl;
  }

  void Gene::SetCisSnps(
    const map<string, vector<Snp*> > & mChr2VecPtSnps,
    const string & anchor, const size_t & radius)
  {
    map<string, vector<Snp*> >::const_iterator itVecPtSnps =
      mChr2VecPtSnps.find(chromosome_);
    if(itVecPtSnps != mChr2VecPtSnps.end()){
      for(size_t snp_id = 0; snp_id < itVecPtSnps->second.size(); ++snp_id){
	Snp * pt_snp = (itVecPtSnps->second)[snp_id];
	int is_in_cis = pt_snp->IsInCis(start_, end_, anchor, radius);
	if(is_in_cis == 1) // snp.coord > gene.end + radius
	  break;
	else if(is_in_cis == -1) // snp.coord < gene.end + radius
	  continue;
	snps_.push_back((itVecPtSnps->second)[snp_id]);
      }
    }
  }

  bool Gene::HasExplevelsInAtLeastOneSubgroup(void) const
  {
    bool res = false;
    for(map<string,vector<double> >::const_iterator it =
	  subgroup2explevels_.begin(); it != subgroup2explevels_.end(); ++it)
      if(it->second.size() > 0){
	res = true;
	break;
      }
    return res;
  }

  bool Gene::HasExplevels(const string & subgroup) const
  {
    return subgroup2explevels_.find(subgroup) != subgroup2explevels_.end();
  }

  bool Gene::HasExplevelsInAllSubgroups(const vector<string> & subgroups) const
  {
    bool res = true;
    for(vector<string>::const_iterator it = subgroups.begin();
	it != subgroups.end(); ++it)
      if(! HasExplevels(*it)){
	res = false;
	break;
      }
    return res;
  }

  bool Gene::HasAtLeastOneCisSnpInAtLeastOneSubgroup(void) const
  {
    bool res = false;
    Snp * pt_snp;
    for(size_t snp_id = 0; snp_id < snps_.size(); ++snp_id){
      pt_snp = snps_[snp_id];
      for(map<string,vector<double> >::const_iterator it =
	    subgroup2explevels_.begin(); it != subgroup2explevels_.end(); ++it)
	if(pt_snp->HasGenotypes(it->first)){
	  res = true;
	  break;
	}
    }
    return res;
  }

  bool Gene::HasAtLeastOneCisSnp(const string & subgroup) const
  {
    bool res = false;
    if(HasExplevels(subgroup)
       && HasAtLeastOneCisSnpInAtLeastOneSubgroup()){
      Snp * pt_snp;
      for(size_t snp_id = 0; snp_id < snps_.size(); ++snp_id){
	pt_snp = snps_[snp_id];
	if(pt_snp->HasGenotypes(subgroup)){
	  res = true;
	  break;
	}
      }
    }
    return res;
  }

  double Gene::GetExplevel(const string & subgroup, const size_t & idx) const
  {
    return subgroup2explevels_.find(subgroup)->second[idx];
  }

  void Gene::TestForAssociations(
    const vector<string> & subgroups,
    const Samples & samples,
    const string & likelihood,
    const string & type_analysis,
    const bool & need_qnorm,
    const Covariates & covariates,
    const Grid & iGridL,
    const Grid & iGridS,
    const string & whichBfs,
    const string & type_errors,
    const float & prop_cov_errors,
    const int & verbose)
  {
#pragma omp parallel for
    for(int idx_snp = 0; idx_snp < snps_.size(); ++idx_snp){
      
      Snp * pt_snp = snps_[idx_snp];
      if(verbose > 0)
	cout << name_ << " (" << GetNbSubgroups() << " subgroups) versus "
	     << pt_snp->name_ << " (" << pt_snp->GetNbSubgroups()
	     << " subgroups)" << endl;
    
      GeneSnpPair gene_snp_pair(name_, pt_snp->name_);
      gsl_permutation * perm = NULL;
    
      if(type_analysis.compare("sep") == 0 ||
	 (type_analysis.compare("join") == 0 && type_errors.compare("uvlr") == 0)){
	for(vector<string>::const_iterator it = subgroups.begin();
	    it != subgroups.end(); ++it){
	  if(HasExplevels(*it) && pt_snp->HasGenotypes(*it))
	    gene_snp_pair.CalcSstatsOneSbgrp(samples, *this, *pt_snp,
					     covariates, *it,
					     likelihood, need_qnorm,
					     perm);
	}
	if(type_analysis.compare("join") == 0)
	  gene_snp_pair.CalcAbfsUvlr(subgroups, whichBfs, iGridL, iGridS);
      }
      else{ // if mvlr or hybrid
	if(! pt_snp->HasGenotypesInAllSubgroups(subgroups)){
	  if(verbose > 0)
	    cerr << "WARNING: skip pair " << name_ << "-" << pt_snp->name_
		 << " because option --error mvlr/hybrid"
		 << " requires genotypes in all subgroups" << endl;
	  continue;
	}
	if(type_errors.compare("mvlr") == 0)
	  gene_snp_pair.CalcAbfsMvlr(subgroups, samples, *this, *pt_snp,
				     covariates, need_qnorm, whichBfs, iGridL,
				     iGridS, prop_cov_errors, perm);
	else if(type_errors.compare("hybrid") == 0)
	  gene_snp_pair.CalcAbfsHybrid(subgroups, samples, *this, *pt_snp,
				       covariates, need_qnorm, whichBfs, iGridL,
				       iGridS, prop_cov_errors, perm);
      }
      gene_snp_pairs_.push_back(gene_snp_pair);
    }
  }

  vector<GeneSnpPair>::const_iterator Gene::BeginPair(void) const
  {
    return gene_snp_pairs_.begin();
  }

  vector<GeneSnpPair>::const_iterator Gene::EndPair(void) const
  {
    return gene_snp_pairs_.end();
  }

/** \brief 
 *  \note nbperms_total is usually 10000, nbperms_sofar may be different 
 *  because of the trick
 */
  double Gene::CalcPermutationPvalue(const size_t & nbperms_total,
				     const size_t & nbperms_sofar,
				     const double & nbperms_more_extreme,
				     const size_t & trick_cutoff,
				     const gsl_rng * rng) const
  {
    double res;
    if(nbperms_sofar == nbperms_total)
      res = nbperms_more_extreme / (nbperms_total + 1);
    else
      res = gsl_ran_flat(rng,
			 ((1 + trick_cutoff) /
			  ((double) (nbperms_sofar + 2))),
			 ((1 + trick_cutoff) /
			  ((double) (nbperms_sofar + 1))));
    return res;
  }

/** \brief Retrieve the lowest genotype p-value over SNPs of the given gene
 *  in the given subgroup
 */
  void Gene::FindMinTruePvaluePerSubgroup(const string & subgroup)
  {
    double pval_true_min_persbgrp_ = 1.0;
    for(vector<GeneSnpPair>::const_iterator it_pair = gene_snp_pairs_.begin();
	it_pair != gene_snp_pairs_.end(); ++it_pair)
      if(it_pair->HasResults(subgroup) &&
	 it_pair->GetBetapvalGeno(subgroup) < pval_true_min_persbgrp_)
	pval_true_min_persbgrp_ = it_pair->GetBetapvalGeno(subgroup);
    subgroup2trueminpval_.insert(make_pair(subgroup, pval_true_min_persbgrp_));
  }

  void Gene::MakePermutationsSepPerSubgroup(
    const string & subgroup,
    const Samples & samples,
    const string & likelihood,
    const bool & need_qnorm,
    const Covariates & covariates,
    const size_t & nb_permutations,
    const int & trick,
    const size_t & trick_cutoff,
    const gsl_rng * rngPerm,
    const gsl_rng * rngTrick)
  {
    gsl_permutation * perm = NULL;
  
    perm = gsl_permutation_calloc(samples.GetTotalNbSamples());
    if(perm == NULL){
      cerr << "ERROR: can't allocate memory for the permutation" << endl;
      exit(1);
    }
  
    subgroup2nbperms_.insert(make_pair(subgroup, 0));
    subgroup2permpval_.insert(make_pair(subgroup, 1));
  
    vector<double> pvals_perm; // per SNP
    double pval_perm_min;
    bool shuffle_only = false;
  
    FindMinTruePvaluePerSubgroup(subgroup);
  
    for(size_t perm_id = 0; perm_id < nb_permutations; ++perm_id){
      gsl_ran_shuffle(rngPerm, perm->data, perm->size, sizeof(size_t));
      if(shuffle_only)
	continue;
    
      ++subgroup2nbperms_[subgroup];
      pvals_perm.assign(snps_.size(), 1.0);
      pval_perm_min = 1;

#pragma omp parallel for shared(pvals_perm)
      for(int idx_snp = 0; idx_snp < snps_.size(); ++idx_snp){
	Snp * pt_snp = snps_[idx_snp];
	GeneSnpPair gene_snp_pair(name_, pt_snp->name_);
	if(HasExplevels(subgroup) && pt_snp->HasGenotypes(subgroup)){
	  gene_snp_pair.CalcSstatsOneSbgrp(samples, *this, *pt_snp, covariates,
					   subgroup, likelihood, need_qnorm,
					   perm);
	  pvals_perm[idx_snp] = gene_snp_pair.GetBetapvalGeno(subgroup);
	}
      }
    
      pval_perm_min = *min_element(pvals_perm.begin(), pvals_perm.end());
    
      if(pval_perm_min <= subgroup2trueminpval_[subgroup])
	++subgroup2permpval_[subgroup];
      if(trick != 0 && subgroup2permpval_[subgroup] == 1 + trick_cutoff){
	if(trick == 1)
	  break;
	else if(trick == 2)
	  shuffle_only = true;
      }
    }
  
    subgroup2permpval_[subgroup] = CalcPermutationPvalue(
      nb_permutations, subgroup2nbperms_[subgroup], subgroup2permpval_[subgroup],
      trick_cutoff, rngTrick);
  
    gsl_permutation_free(perm);
  }

  size_t Gene::GetNbGeneSnpPairs(const string & subgroup) const
  {
    size_t res = 0;
    for(vector<GeneSnpPair>::const_iterator it = gene_snp_pairs_.begin();
	it != gene_snp_pairs_.end(); ++it)
      if(it->HasResults(subgroup))
	++res;
    return res;
  }

  double Gene::GetPermutationPvalueSep(const string & subgroup) const
  {
    return subgroup2permpval_.find(subgroup)->second;
  }

  size_t Gene::GetNbPermutationsSep(const string & subgroup) const
  {
    return subgroup2nbperms_.find(subgroup)->second;
  }

  double Gene::GetTrueMinPval(const string & subgroup) const
  {
    return subgroup2trueminpval_.find(subgroup)->second;
  }

/** \brief Retrieve the lowest genotype p-value over SNPs and subgroups 
 *  of the given gene
 */
  void Gene::FindMinTruePvalueAllSubgroups(void)
  {
    pval_true_min_allsbgrps_ = 1.0;
    for(vector<GeneSnpPair>::const_iterator it_pair = gene_snp_pairs_.begin();
	it_pair != gene_snp_pairs_.end(); ++it_pair)
      for(map<string,vector<double> >::const_iterator it_sbgrp =
	    subgroup2explevels_.begin(); it_sbgrp != subgroup2explevels_.end();
	  ++it_sbgrp)
	if(it_pair->HasResults(it_sbgrp->first) &&
	   it_pair->GetBetapvalGeno(it_sbgrp->first) < pval_true_min_allsbgrps_)
	  pval_true_min_allsbgrps_ = it_pair->GetBetapvalGeno(it_sbgrp->first);
  }

  void Gene::MakePermutationsSepAllSubgroups(
    const vector<string> & subgroups,
    const Samples & samples,
    const string & likelihood,
    const bool & need_qnorm,
    const Covariates & covariates,
    const size_t & nb_permutations,
    const int & trick,
    const size_t & trick_cutoff,
    const gsl_rng * rngPerm,
    const gsl_rng * rngTrick)
  {
    gsl_permutation * perm = NULL;
  
    perm = gsl_permutation_calloc(samples.GetTotalNbSamples());
    if(perm == NULL){
      cerr << "ERROR: can't allocate memory for the permutation" << endl;
      exit(1);
    }
  
    nbpermutations_sep_allsbgrps_ = 0;
    pval_perm_sep_allsbgrps_ = 1;
  
    vector<double> pvals_perm; // per SNP over subgroups
    double pval_perm_min;
    bool shuffle_only = false;
  
    FindMinTruePvalueAllSubgroups();
  
    for(size_t perm_id = 0; perm_id < nb_permutations; ++perm_id){
      gsl_ran_shuffle(rngPerm, perm->data, perm->size, sizeof(size_t));
      if(shuffle_only)
	continue;
    
      ++nbpermutations_sep_allsbgrps_;
      pvals_perm.assign(snps_.size(), 1.0);
      pval_perm_min = 1;
    
#pragma omp parallel for shared(pvals_perm)
      for(int idx_snp = 0; idx_snp < snps_.size(); ++idx_snp){
	Snp * pt_snp = snps_[idx_snp];
	GeneSnpPair gene_snp_pair(name_, pt_snp->name_);
	double pval_perm_sbgrp, pval_perm_min_sbgrp = 1;
	for(vector<string>::const_iterator it_sbgrp = subgroups.begin();
	    it_sbgrp != subgroups.end(); ++it_sbgrp){
	  if(HasExplevels(*it_sbgrp) && pt_snp->HasGenotypes(*it_sbgrp)){
	    gene_snp_pair.CalcSstatsOneSbgrp(samples, *this, *pt_snp, covariates,
					     *it_sbgrp, likelihood, need_qnorm,
					     perm);
	    pval_perm_sbgrp = gene_snp_pair.GetBetapvalGeno(*it_sbgrp);
	    if(pval_perm_sbgrp < pval_perm_min_sbgrp)
	      pval_perm_min_sbgrp = pval_perm_sbgrp;
	  }
	  pvals_perm[idx_snp] = pval_perm_min_sbgrp;
	}
      }
    
      pval_perm_min = *min_element(pvals_perm.begin(), pvals_perm.end());
    
      if(pval_perm_min <= pval_true_min_allsbgrps_)
	++pval_perm_sep_allsbgrps_;
      if(trick != 0 && pval_perm_sep_allsbgrps_ == 1 + trick_cutoff){
	if(trick == 1)
	  break;
	else if(trick == 2)
	  shuffle_only = true;
      }
    }
  
    pval_perm_sep_allsbgrps_ = CalcPermutationPvalue(
      nb_permutations, nbpermutations_sep_allsbgrps_, pval_perm_sep_allsbgrps_,
      trick_cutoff, rngTrick);
  
    gsl_permutation_free(perm);
  }

/** \brief Retrieve the highest log10(ABF) over SNPs of the given gene
 *  \note whichPermBf is 'gen', 'sin', 'gen-sin' or 'all'
 */
  void Gene::FindMaxTrueL10Abf(const string & whichPermBf)
  {
    l10_abf_true_max_ = - numeric_limits<double>::infinity();
    for(vector<GeneSnpPair>::const_iterator it_pair = gene_snp_pairs_.begin();
	it_pair != gene_snp_pairs_.end(); ++it_pair)
      if(it_pair->GetWeightedAbf(whichPermBf) > l10_abf_true_max_)
	l10_abf_true_max_ = it_pair->GetWeightedAbf(whichPermBf);
  }

/** \brief Average the log10(ABF) over SNPs of the given gene
 *  \note whichPermBf is 'gen', 'sin', 'gen-sin' or 'all'
 */
  void Gene::AvgTrueL10Abfs(const string & whichPermBf)
  {
    vector<double> l10_abfs;
    for(vector<GeneSnpPair>::const_iterator it_pair = gene_snp_pairs_.begin();
	it_pair != gene_snp_pairs_.end(); ++it_pair){
      if(! isNan(it_pair->GetWeightedAbf(whichPermBf)))
	l10_abfs.push_back(it_pair->GetWeightedAbf(whichPermBf));
    }
    l10_abf_true_avg_= log10_weighted_sum(&(l10_abfs[0]), l10_abfs.size());
  }

  void Gene::MakePermutationsJoin(const vector<string> & subgroups,
				  const Samples & samples,
				  const string & likelihood,
				  const bool & need_qnorm,
				  const Covariates & covariates,
				  const Grid & iGridL,
				  const Grid & iGridS,
				  const string & type_errors,
				  const float & prop_cov_errors,
				  const size_t & nb_permutations,
				  const int & trick,
				  const size_t & trick_cutoff,
				  const string & whichPermBf,
				  const bool & useMaxBfOverSnps,
				  const gsl_rng * rngPerm,
				  const gsl_rng * rngTrick)
  {
    gsl_permutation * perm = NULL;
  
    perm = gsl_permutation_calloc(samples.GetTotalNbSamples());
    if(perm == NULL){
      cerr << "ERROR: can't allocate memory for the permutation" << endl;
      exit(1);
    }
  
    nbpermutations_join_ = 0;
    pval_perm_join_ = 1;
    l10_abf_perm_med_ = NaN;
  
    vector<double> l10_abfs_perm_snps, // per SNP
      l10_abfs_perms(nb_permutations, NaN); // for the median
    double l10_abf_perm_max, l10_abf_perm_avg;
    bool shuffle_only = false;
  
    if(useMaxBfOverSnps)
      FindMaxTrueL10Abf(whichPermBf);
    else
      AvgTrueL10Abfs(whichPermBf);
    
    for(size_t perm_id = 0; perm_id < nb_permutations; ++perm_id){
      gsl_ran_shuffle(rngPerm, perm->data, perm->size, sizeof(size_t));
      if(shuffle_only)
	continue;
    
      ++nbpermutations_join_;
      l10_abfs_perm_snps.assign(snps_.size(), 0.0);
      if(useMaxBfOverSnps)
	l10_abf_perm_max = - numeric_limits<double>::infinity();
      else
	l10_abf_perm_avg = NaN;
    
#pragma omp parallel for shared(l10_abfs_perm_snps)
      for(int idx_snp = 0; idx_snp < snps_.size(); ++idx_snp){
	Snp * pt_snp = snps_[idx_snp];
	GeneSnpPair gene_snp_pair(name_, pt_snp->name_);
	if(type_errors.compare("uvlr") == 0){
	  for(vector<string>::const_iterator it_sbgrp = subgroups.begin();
	      it_sbgrp != subgroups.end(); ++it_sbgrp)
	    if(HasExplevels(*it_sbgrp) && pt_snp->HasGenotypes(*it_sbgrp))
	      gene_snp_pair.CalcSstatsOneSbgrp(samples, *this, *pt_snp,
					       covariates, *it_sbgrp,
					       likelihood, need_qnorm, perm);
	  gene_snp_pair.CalcAbfsUvlr(subgroups, whichPermBf, iGridL, iGridS);
	}
	else{ // if mvlr or hybrid
	  if(! pt_snp->HasGenotypesInAllSubgroups(subgroups))
	    continue;
	  if(type_errors.compare("mvlr") == 0)
	    gene_snp_pair.CalcAbfsMvlr(subgroups, samples, *this, *pt_snp,
				       covariates, need_qnorm, whichPermBf,
				       iGridL, iGridS, prop_cov_errors, perm);
	  else if(type_errors.compare("hybrid") == 0)
	    gene_snp_pair.CalcAbfsHybrid(subgroups, samples, *this, *pt_snp,
					 covariates, need_qnorm, whichPermBf,
					 iGridL, iGridS, prop_cov_errors, perm);
	}
	l10_abfs_perm_snps[idx_snp] = gene_snp_pair.GetWeightedAbf(whichPermBf);
      }
    
      if(useMaxBfOverSnps){
	l10_abf_perm_max = *max_element(l10_abfs_perm_snps.begin(),
					l10_abfs_perm_snps.end());
	if(l10_abf_perm_max >= l10_abf_true_max_)
	  ++pval_perm_join_;
	l10_abfs_perms[perm_id] = l10_abf_perm_max;
      }
      else{ // if avg over SNPs
	l10_abf_perm_avg = log10_weighted_sum(&(l10_abfs_perm_snps[0]),
					      l10_abfs_perm_snps.size());
	if(l10_abf_perm_avg >= l10_abf_true_avg_)
	  ++pval_perm_join_;
	l10_abfs_perms[perm_id] = l10_abf_perm_avg;
      }
      if(trick != 0 && pval_perm_join_ == 1 + trick_cutoff){
	if(trick == 1)
	  break;
	else if(trick == 2)
	  shuffle_only = true;
      }
    }
  
    pval_perm_join_ = CalcPermutationPvalue(
      nb_permutations, nbpermutations_join_, pval_perm_join_, trick_cutoff,
      rngTrick);
  
    l10_abf_perm_med_ = median(l10_abfs_perms.begin(),
			       l10_abfs_perms.begin() + nbpermutations_join_ + 1);
  
    gsl_permutation_free(perm);
  }

  double Gene::GetTrueL10Abf(const bool & use_max_bf) const
  {
    if(use_max_bf)
      return l10_abf_true_max_;
    else
      return l10_abf_true_avg_;
  }

} // namespace quantgen
