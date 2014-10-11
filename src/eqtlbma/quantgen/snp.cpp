/** \file snp.cpp
 *
 *  `Snp' is a class 
 *  Copyright (C) 2013 Timothee Flutre
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Snpral Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Snpral Public License for more details.
 *
 *  You should have received a copy of the GNU Snpral Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "quantgen/snp.hpp"

using namespace std;

using namespace utils;

namespace quantgen {

  Snp::Snp(void)
  {
  }

  Snp::Snp(const string & name)
  {
    name_ = name;
  }

  Snp::Snp(const string & name, const string & chromosome,
	   const string & pos)
  {
    name_ = name;
    chromosome_ = chromosome;
    pos_ = atol(pos.c_str());
  }

  bool operator==(const Snp& lhs, const Snp& rhs)
  {
    return(lhs.GetName().compare(rhs.GetName()) == 0
	   && lhs.GetChromosome().compare(rhs.GetChromosome()) == 0
	   && lhs.GetPosition() == rhs.GetPosition()
	   && lhs.GetNbSubgroups() == rhs.GetNbSubgroups());
  }

  bool operator!=(const Snp& lhs, const Snp& rhs)
  {
    return !operator==(lhs,rhs);
  }

  bool operator< (const Snp& lhs, const Snp& rhs)
  {
    if (lhs.GetChromosome().compare(rhs.GetChromosome()) != 0) {
      cerr << "ERROR: " << lhs.GetName() << " and " << rhs.GetName()
	   << " are on different chromosomes, thus they can't be sorted" << endl;
      exit(EXIT_FAILURE);
    }
    return(lhs.GetPosition() < rhs.GetPosition());
  }

  bool operator> (const Snp& lhs, const Snp& rhs)
  {
    return  operator< (rhs,lhs);
  }

  bool operator<=(const Snp& lhs, const Snp& rhs)
  {
    return !operator> (lhs,rhs);
  }

  bool operator>=(const Snp& lhs, const Snp& rhs)
  {
    return !operator< (lhs,rhs);
  }

  bool pt_snp_lt_pt_snp(const Snp* pt_lhs, const Snp* pt_rhs)
  {
    return *pt_lhs < *pt_rhs;
  }

  void Snp::AddSubgroupFromImputeLine(
    const vector<string>::const_iterator & begin,
    const vector<string>::const_iterator & end,
    vector<double> & genotypes,
    double & minor_allele_freq)
  {
    if((end - begin) % 3 != 0){
      cerr << "ERROR: SNP " << name_
	   << " from IMPUTE file has not the right number of columns" << endl;
      exit(EXIT_FAILURE);
    }
    size_t nb_samples = static_cast<size_t>((end - begin) / 3);
    genotypes.assign(nb_samples, NaN);
    minor_allele_freq = 0.0;
    double AA, BB, AB;
    for(size_t i = 0; i < nb_samples; ++i){
      AA = atof((begin+3*i)->c_str());
      AB = atof((begin+3*i+1)->c_str());
      BB = atof((begin+3*i+2)->c_str());
      if(AA == 0 && AB == 0 && BB == 0){
	minor_allele_freq = NaN;
      }
      else{
	genotypes[i] = 0 * AA + 1 * AB + 2 * BB;
	minor_allele_freq += genotypes[i];
      }
    }
    if(! isNan(minor_allele_freq)){
      minor_allele_freq /= (2 * nb_samples);
      minor_allele_freq = (minor_allele_freq <= 0.5 ?
			   minor_allele_freq
			   : 1 - minor_allele_freq);
    }
  }

  void Snp::AddSubgroupFromVcfLine(
    const vector<string>::const_iterator & begin,
    const vector<string>::const_iterator & end,
    const size_t & idx_gt,
    vector<double> & genotypes,
    double & minor_allele_freq)
  {
    size_t nb_samples = static_cast<size_t>(end - begin);
    genotypes.assign(nb_samples, NaN);
    minor_allele_freq = 0.0;
    vector<string> tokens2, tokens3;
    for(size_t i = 0; i < nb_samples; ++i){
      split(*(begin+i), ":", tokens2);
      if(tokens2[idx_gt].find(".") != string::npos) {
	minor_allele_freq = NaN;
      }
      else{
	split(tokens2[idx_gt], "|/", tokens3);
	genotypes[i] = 0;
	if(tokens3[0].compare("1") == 0)
	  genotypes[i] += 1;
	if(tokens3[1].compare("1") == 0)
	  genotypes[i] += 1;
	minor_allele_freq += genotypes[i];
      }
    }
    if(! isNan(minor_allele_freq)){
      minor_allele_freq /= (2 * nb_samples);
      minor_allele_freq = (minor_allele_freq <= 0.5 ?
			   minor_allele_freq
			   : 1 - minor_allele_freq);
    }
  }

  void Snp::AddSubgroupFromDoseLine(
    const vector<string>::const_iterator & begin,
    const vector<string>::const_iterator & end,
    vector<double> & genotypes,
    double & minor_allele_freq)
  {
    size_t nb_samples = static_cast<size_t>(end - begin);
    genotypes.assign(nb_samples, NaN);
    minor_allele_freq = 0.0;
    for(size_t i = 0; i < nb_samples; ++i){
      if((begin+i)->compare("NA") == 0 || (begin+i)->compare("na") == 0
	  || (begin+i)->compare("NaN") == 0 || (begin+i)->compare("nan") == 0){
	minor_allele_freq = NaN;
      }
      else{
	genotypes[i] = atof((begin+i)->c_str());
	minor_allele_freq += genotypes[i];
      }
    }
    if(! isNan(minor_allele_freq)){
      minor_allele_freq /= (2 * nb_samples);
      minor_allele_freq = (minor_allele_freq <= 0.5 ?
			   minor_allele_freq
			   : 1 - minor_allele_freq);
    }
  }

  void Snp::AddSubgroup(const string & subgroup,
			const vector<string>::const_iterator & begin,
			const vector<string>::const_iterator & end,
			const string & format,
			const size_t & idx_gt)
  {
    vector<double> genotypes;
    double minor_allele_freq;
    
    if(format.compare("impute") == 0)
      AddSubgroupFromImputeLine(begin, end, genotypes, minor_allele_freq);
    else if(format.compare("vcf") == 0)
      AddSubgroupFromVcfLine(begin, end, idx_gt, genotypes, minor_allele_freq);
    else if(format.compare("dose") == 0)
      AddSubgroupFromDoseLine(begin, end, genotypes, minor_allele_freq);
    else{
      cerr << "ERROR: genotype format '" << format << "' is not recognized" << endl;
      exit(EXIT_FAILURE);
    }
    
    Genotypes genos;
    genos.values = genotypes;
    genos.maf = minor_allele_freq;
    subgroup2genotypes_.insert(make_pair(subgroup, genos));
  }

  size_t Snp::GetNbSamples(const string & subgroup) const
  {
    return subgroup2genotypes_.find(subgroup)->second.values.size();
  }

  double Snp::GetMinorAlleleFreq(const string & subgroup) const
  {
    return subgroup2genotypes_.find(subgroup)->second.maf;
  }

  void Snp::Show(ostream & os)
  {
    os << name_ << " " << chromosome_ << " " << pos_ << endl
       << GetNbSubgroups() << " subgroups" << endl;
    for (map<string,Genotypes>::const_iterator it =
	   subgroup2genotypes_.begin(); it != subgroup2genotypes_.end(); ++it)
      os << it->first << ": " << GetNbSamples(it->first) << " samples "
	 << " (maf=" << GetMinorAlleleFreq(it->first) << ")" << endl;
  }

  void Snp::EraseIfMissingValuesPerSubgroup()
  {
    map<string,Genotypes>::iterator it = subgroup2genotypes_.begin();
    while(it != subgroup2genotypes_.end()){
      if(isNan(GetMinorAlleleFreq(it->first)))
	subgroup2genotypes_.erase(it++);
      else
	++it;
    }
  }

  void Snp::EraseIfLowMafPerSubgroup(const double & min_maf)
  {
    map<string,Genotypes>::iterator it = subgroup2genotypes_.begin();
    while(it != subgroup2genotypes_.end()){
      if(GetMinorAlleleFreq(it->first) < min_maf)
	subgroup2genotypes_.erase(it++);
      else
	++it;
    }
  }

  void Snp::DuplicateGenotypesFromFirstSubgroup(const string & subgroup_old,
						const string & subgroup_new)
  {
    subgroup2genotypes_.insert(make_pair(subgroup_new,
					 subgroup2genotypes_[subgroup_old]));
  }

  bool Snp::HasGenotypesInAtLeastOneSubgroup(void) const
  {
    bool res = false;
    for (map<string,Genotypes>::const_iterator it =
	   subgroup2genotypes_.begin(); it != subgroup2genotypes_.end(); ++it)
      if (GetNbSamples(it->first) > 0){
	res = true;
	break;
      }
    return res;
  }

  int Snp::IsInCis(const size_t & start, const size_t & end,
		   const string & anchor, const size_t & radius) const
  {
    int res = -1; // snp.pos < anchor - radius
    if (anchor.compare("TSS+TES") == 0) {
      if (((start >= radius &&
	    pos_ >= start - radius) ||
	   (start < radius)) &&
	  pos_ <= end + radius)
	res = 0; // snp.pos in [anchor-radius,anchor+radius]
      else if (pos_ > end + radius)
	res = 1; // snp.pos > anchor - radius
    }
    else if (anchor.compare("TSS") == 0) {
      if (((start >= radius &&
	    pos_ >= start - radius) ||
	   (start < radius)) &&
	  pos_ <= start + radius)
	res = 0;
      else if (pos_ > start + radius)
	res = 1;
    }
    return res;
  }

  bool Snp::HasGenotypes(const string & subgroup) const
  {
    return (subgroup2genotypes_.find(subgroup) != subgroup2genotypes_.end() &&
	    GetNbSamples(subgroup) > 0);
  }

  bool Snp::HasGenotypesInAllSubgroups(const vector<string> & subgroups) const
  {
    bool res = true;
    for(vector<string>::const_iterator it = subgroups.begin();
	it != subgroups.end(); ++it)
      if(! HasGenotypes(*it)){
	res = false;
	break;
      }
    return res;
  }

  double Snp::GetGenotype(const string & subgroup, const size_t & idx) const
  {
    return subgroup2genotypes_.find(subgroup)->second.values[idx];
  }

} // namespace quantgen
