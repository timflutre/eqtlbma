/** \file gene.cpp
 *
 *  `Gene' is a class 
 *  Copyright (C) 2013-2015 Timothée Flutre
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
  
  Gene::Gene(const string & name)
  {
    name_ = name;
  }
  
  Gene::Gene(const string & name, const string & chromosome,
	     const string & start, const string & end)
  {
    name_ = name;
    chromosome_ = chromosome;
    start_ = atol(start.c_str()) + 1; // in BED format, start starts at 0
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
      exit(EXIT_FAILURE);
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
	 || it->compare("NaN") == 0 || it->compare("nan") == 0)
	++i;
      else{
	explevels[i] = atof(it->c_str());
	++i;
      }
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

  void Gene::AddCisSnp(const Snp * pt_snp)
  {
    snps_.push_back(pt_snp);
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
	AddCisSnp((itVecPtSnps->second)[snp_id]);
      }
    }
  }

  bool Gene::HasCisSnp(const Snp * pt_snp)
  {
    return find(snps_.begin(), snps_.end(), pt_snp) != snps_.end();
  }

  size_t Gene::FindIdxSnp(const Snp * pt_snp)
  {
    size_t idx_snp = string::npos;
    vector<const Snp*>::iterator it = find(snps_.begin(), snps_.end(), pt_snp);
    if(it != snps_.end())
      idx_snp = it - snps_.begin();
    return idx_snp;
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
    const Snp * pt_snp;
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
      const Snp * pt_snp;
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

  void Gene::GetSubgroupsWithExpLevels(vector<string> & subgroups_with_exp)
  {
    subgroups_with_exp.clear();
    for(map<string,vector<double> >::const_iterator it =
          subgroup2explevels_.begin(); it != subgroup2explevels_.end(); ++it)
      if(HasExplevels(it->first))
        subgroups_with_exp.push_back(it->first);
  }

  double Gene::GetExplevel(const string & subgroup, const size_t & idx) const
  {
    return subgroup2explevels_.find(subgroup)->second[idx];
  }

  vector<GeneSnpPair>::iterator Gene::AddGeneSnpPair(
    const string & snp_name,
    const string & error_model)
  {
    gene_snp_pairs_.push_back(GeneSnpPair(name_, snp_name, error_model));
    return gene_snp_pairs_.end() - 1;
  }

  vector<GeneSnpPair>::iterator Gene::FindGeneSnpPair(
    const size_t & idx_snp)
  {
    return gene_snp_pairs_.begin() + idx_snp;
  }

  void Gene::TestForAssociations(
    const bool & hasDataNotSstats,
    const vector<string> & subgroups,
    const Samples & samples,
    const string & likelihood,
    const string & analysis,
    const bool & need_qnorm,
    const Covariates & covariates,
    const Grid & iGridL,
    const Grid & iGridS,
    const string & whichBfs,
    const string & error_model,
    const float & prop_cov_errors,
    const int & verbose)
  {
    const Snp * pt_snp = NULL;
    vector<GeneSnpPair>::iterator it_gsp;
    gsl_permutation * perm = NULL;
    
    // To use OpenMP here, need to have a temporary vector shared between all
    // threads, replacing gene_snp_pairs_ inside the loop.
    // As this requires more memory, I prefer not to use OpenMP here.
    for(size_t idx_snp = 0; idx_snp < snps_.size(); ++idx_snp){
    
      pt_snp = snps_[idx_snp];
      if(verbose > 0)
	cout << name_ << " (" << GetNbSubgroups() << " subgroups) versus "
	     << pt_snp->name_ << " (" << pt_snp->GetNbSubgroups()
	     << " subgroups)" << endl;
    
      if(hasDataNotSstats)
	it_gsp = AddGeneSnpPair(pt_snp->name_, analysis);
      else
	it_gsp = FindGeneSnpPair(idx_snp);
      perm = NULL;
    
      if(analysis.compare("sep") == 0 ||
	 (analysis.compare("join") == 0 && error_model.compare("uvlr") == 0)){
	if(hasDataNotSstats){
	  for(vector<string>::const_iterator it = subgroups.begin();
	      it != subgroups.end(); ++it){
	    if(HasExplevels(*it) && pt_snp->HasGenotypes(*it))
	      it_gsp->CalcSstatsOneSbgrp(samples, *this, *pt_snp,
					 covariates, *it,
					 likelihood, need_qnorm,
					 perm);
	  }
	}
	if(analysis.compare("join") == 0)
	  it_gsp->CalcAbfsUvlr(subgroups, whichBfs, iGridL, iGridS);
      }
      else{ // if mvlr or hybrid
	if(! pt_snp->HasGenotypesInAllSubgroups(subgroups)){
	  if(verbose > 0)
	    cerr << "WARNING: skip pair " << name_ << "-" << pt_snp->name_
		 << " because option --error mvlr/hybrid"
		 << " requires genotypes in all subgroups" << endl;
	  continue;
	}
	if(error_model.compare("mvlr") == 0)
	  it_gsp->CalcAbfsMvlr(subgroups, samples, *this, *pt_snp,
			       covariates, need_qnorm, whichBfs, iGridL,
			       iGridS, prop_cov_errors, perm);
	else if(error_model.compare("hybrid") == 0)
	  it_gsp->CalcAbfsHybrid(subgroups, samples, *this, *pt_snp,
				 covariates, need_qnorm, whichBfs, iGridL,
				 iGridS, prop_cov_errors, perm);
      }
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
    size_t & nb_permutations,
    const int & trick,
    const size_t & trick_cutoff,
    const gsl_rng * rngPerm,
    const gsl_rng * rngTrick)
  {
    gsl_permutation * perm = NULL;
  
    perm = gsl_permutation_calloc(samples.GetTotalNbSamples());
    if(perm == NULL){
      cerr << "ERROR: can't allocate memory for the permutation" << endl;
      exit(EXIT_FAILURE);
    }
  
    subgroup2nbperms_.insert(make_pair(subgroup, 0));
    subgroup2permpval_.insert(make_pair(subgroup, 1));
  
    vector<double> pvals_perm; // per SNP
    double pval_perm_min;
    bool shuffle_only = false;
  
    FindMinTruePvaluePerSubgroup(subgroup);
    size_t nb_all_permutations_ = nb_permutations;
    for(size_t perm_id = 0; perm_id < nb_all_permutations_; ++perm_id){
      gsl_ran_shuffle(rngPerm, perm->data, perm->size, sizeof(size_t));
      if(shuffle_only)
	continue;
    
      pvals_perm.assign(snps_.size(), 1.0);
      pval_perm_min = 1;

#pragma omp parallel for shared(pvals_perm)
      for(int idx_snp = 0; idx_snp < (int) snps_.size(); ++idx_snp){
	const Snp * pt_snp = snps_[idx_snp];
	GeneSnpPair gene_snp_pair(name_, pt_snp->name_);
	if(HasExplevels(subgroup) && pt_snp->HasGenotypes(subgroup)){
	  gene_snp_pair.CalcSstatsOneSbgrp(samples, *this, *pt_snp, covariates,
					   subgroup, likelihood, need_qnorm,
					   perm);
	  pvals_perm[idx_snp] = gene_snp_pair.GetBetapvalGeno(subgroup);
	}
      }
    
      pval_perm_min = *min_element(pvals_perm.begin(), pvals_perm.end());
      if(isNan(pval_perm_min)) {
        nb_permutations--;
        continue;
      }
      ++subgroup2nbperms_[subgroup];
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
    size_t & nb_permutations,
    const int & trick,
    const size_t & trick_cutoff,
    const gsl_rng * rngPerm,
    const gsl_rng * rngTrick)
  {
    gsl_permutation * perm = NULL;
  
    perm = gsl_permutation_calloc(samples.GetTotalNbSamples());
    if(perm == NULL){
      cerr << "ERROR: can't allocate memory for the permutation" << endl;
      exit(EXIT_FAILURE);
    }
  
    nbpermutations_sep_allsbgrps_ = 0;
    pval_perm_sep_allsbgrps_ = 1;
  
    vector<double> pvals_perm; // per SNP over subgroups
    double pval_perm_min;
    bool shuffle_only = false;
  
    FindMinTruePvalueAllSubgroups();
    size_t nb_all_permutations_ = nb_permutations;
    for(size_t perm_id = 0; perm_id < nb_all_permutations_; ++perm_id){
      gsl_ran_shuffle(rngPerm, perm->data, perm->size, sizeof(size_t));
      if(shuffle_only)
	continue;
    
      pvals_perm.assign(snps_.size(), 1.0);
      pval_perm_min = 1;
    
#pragma omp parallel for shared(pvals_perm)
      for(int idx_snp = 0; idx_snp < (int) snps_.size(); ++idx_snp){
	const Snp * pt_snp = snps_[idx_snp];
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
      if (isNan(pval_perm_min)) {
        nb_permutations--;
        continue;
      }
      ++nbpermutations_sep_allsbgrps_;
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
				  const string & error_model,
				  const float & prop_cov_errors,
				  size_t & nb_permutations,
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
      exit(EXIT_FAILURE);
    }
  
    nbpermutations_join_ = 0;
    pval_perm_join_ = 1;
    l10_abf_perm_med_ = NaN;
  
    vector<double> l10_abfs_perm_snps, // per SNP
      l10_abfs_perms; // for the median
    double l10_abf_perm_max, l10_abf_perm_avg;
    bool shuffle_only = false;
  
    if(useMaxBfOverSnps)
      FindMaxTrueL10Abf(whichPermBf);
    else
      AvgTrueL10Abfs(whichPermBf);

    size_t nb_all_permutations_ = nb_permutations;
    for(size_t perm_id = 0; perm_id < nb_all_permutations_; ++perm_id){
      gsl_ran_shuffle(rngPerm, perm->data, perm->size, sizeof(size_t));
      if(shuffle_only)
	continue;
    
      l10_abfs_perm_snps.assign(snps_.size(), 0.0);
      if(useMaxBfOverSnps)
	l10_abf_perm_max = - numeric_limits<double>::infinity();
      else
	l10_abf_perm_avg = NaN;
    
#pragma omp parallel for shared(l10_abfs_perm_snps)
      for(int idx_snp = 0; idx_snp < (int) snps_.size(); ++idx_snp){
	const Snp * pt_snp = snps_[idx_snp];
	GeneSnpPair gene_snp_pair(name_, pt_snp->name_);
	if(error_model.compare("uvlr") == 0){
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
	  if(error_model.compare("mvlr") == 0)
	    gene_snp_pair.CalcAbfsMvlr(subgroups, samples, *this, *pt_snp,
				       covariates, need_qnorm, whichPermBf,
				       iGridL, iGridS, prop_cov_errors, perm);
	  else if(error_model.compare("hybrid") == 0)
	    gene_snp_pair.CalcAbfsHybrid(subgroups, samples, *this, *pt_snp,
					 covariates, need_qnorm, whichPermBf,
					 iGridL, iGridS, prop_cov_errors, perm);
	}
	l10_abfs_perm_snps[idx_snp] = gene_snp_pair.GetWeightedAbf(whichPermBf);
      }
    
      if(useMaxBfOverSnps){
	l10_abf_perm_max = *max_element(l10_abfs_perm_snps.begin(),
					l10_abfs_perm_snps.end());
    if (isNan(l10_abf_perm_max)) {
      nb_permutations--;
      continue;
    }
    ++nbpermutations_join_;
	if(l10_abf_perm_max >= l10_abf_true_max_)
	  ++pval_perm_join_;
	l10_abfs_perms.push_back(l10_abf_perm_max);
      }
      else{ // if avg over SNPs
	l10_abf_perm_avg = log10_weighted_sum(&(l10_abfs_perm_snps[0]),
					      l10_abfs_perm_snps.size());
    if (isNan(l10_abf_perm_avg)) {
      nb_permutations--;
      continue;
    }
    ++nbpermutations_join_;
	if(l10_abf_perm_avg >= l10_abf_true_avg_)
	  ++pval_perm_join_;
	l10_abfs_perms.push_back(l10_abf_perm_avg);
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
