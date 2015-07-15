/** \file data_loader.cpp
 *
 *  `data_loader' gathers functions useful to load gen-etics/-omics data
 *  Copyright (C) 2011-2015 Timothée Flutre
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

#include <sys/stat.h>

#include <algorithm>

#include "quantgen/data_loader.hpp"
#include "quantgen/samples.hpp"

#include "utils/utils_io.hpp"

#include "tabix/bgzf.h"
#include "tabix/tabix.h"

using namespace std;

using namespace utils;

namespace quantgen {
  
  map<string,string> loadTwoColumnFile(const string & file,
                                       const int & verbose)
  {
    map<string,string> mItems;
    if(file.empty())
      return mItems;
    
    string line;
    gzFile stream;
    vector<string> tokens;
    size_t nb_lines = 0;
    
    openFile(file, stream, "rb");
    if(verbose > 0)
      cout <<"load file " << file << " ..." << endl;
    
    while(getline(stream, line)){
      nb_lines++;
      split(line, " \t,", tokens);
      if(tokens.size() != 2){
        cerr << "ERROR: file " << file << " should have only two columns"
             << " at line " << nb_lines << endl;
        exit(EXIT_FAILURE);
      }
      if(tokens[0][0] == '#')
        continue;
      if(mItems.find(tokens[0]) == mItems.end())
        mItems.insert(make_pair(tokens[0], tokens[1]));
    }
    
    if(! gzeof(stream)){
      cerr << "ERROR: can't read successfully file "
           << file << " up to the end" << endl;
      exit (1);
    }
    closeFile(file, stream);
    
    if(verbose > 0)
      cout << "items loaded: " << mItems.size() << endl;
    
    return mItems;
  }
  
  void loadListsGenoExplevelAndCovarFiles(
    const string & file_genopaths,
    const string & file_exppaths,
    const string & file_covarpaths,
    const vector<string> & subgroups_tokeep,
    const string & error_model,
    const int & verbose,
    map<string,string> & subgroup2genofile,
    map<string,string> & subgroup2explevelfile,
    map<string,string> & subgroup2covarfile,
    vector<string> & subgroups)
  {
    map<string,string>::iterator it;
    
    // load paths to exp levels files
    subgroup2explevelfile = loadTwoColumnFile(file_exppaths, verbose);
    
    // erase subgroups with exp levels but not to keep
    it = subgroup2explevelfile.begin();
    while(it != subgroup2explevelfile.end()){
      if(! subgroups_tokeep.empty()
         && find(subgroups_tokeep.begin(), subgroups_tokeep.end(), it->first)
         == subgroups_tokeep.end())
        subgroup2explevelfile.erase(it++);
      else
        ++it;
    }
    
    // load paths to geno files
    subgroup2genofile = loadTwoColumnFile(file_genopaths, verbose);
    
    // erase subgroups with geno but no exp levels
    it = subgroup2genofile.begin();
    while(it != subgroup2genofile.end()){
      if(subgroup2explevelfile.find(it->first) == subgroup2explevelfile.end())
        subgroup2genofile.erase(it++);
      else
        ++it;
    }
    
    // erase subgroups with exp levels but no geno
    it = subgroup2explevelfile.begin();
    while(it != subgroup2explevelfile.end()){
      if(subgroup2genofile.find(it->first) == subgroup2genofile.end())
        subgroup2explevelfile.erase(it++);
      else
        ++it;
    }
    
    // for correlated errors, check that, at least, all subgroups 
    // have genotypes in the same file
    if(error_model != "uvlr")
      for(it = subgroup2genofile.begin(); it != subgroup2genofile.end(); ++it)
        if(it->second.compare(subgroup2genofile.begin()->second) != 0){
          cerr << "ERROR: --error mvlr/hybrid requires the same genotypes in a single file for all subgroups"
               << endl;
          exit(EXIT_FAILURE);
        }
    
    // identify the names of all remaining subgroups to be considered
    for(it = subgroup2explevelfile.begin(); it != subgroup2explevelfile.end();
        ++it)
      subgroups.push_back(it->first);
    
    // load paths to covariate files
    subgroup2covarfile = loadTwoColumnFile(file_covarpaths, verbose);
    
    // erase subgroups with covariates but neither geno nor exp levels
    it = subgroup2covarfile.begin();
    while(it != subgroup2covarfile.end()){
      if(find(subgroups.begin(), subgroups.end(), it->first) == subgroups.end())
        subgroup2covarfile.erase(it++);
      else++it;
    }
    
    if(verbose > 0)
      cout << "analyze " << subgroups.size() << " subgroup"
           << (subgroups.size() > 1 ? "s" : "")
           <<" (identifier):" << endl;
    for(vector<string>::const_iterator it = subgroups.begin();
        it != subgroups.end(); ++it)
      cout << *it << " (" << it - subgroups.begin() + 1 << ")" << endl;
  }

  void loadSamplesFromMatrixEqtl(
    const map<string, string> & subgroup2file,
    const string & data_type,
    const int & verbose,
    vector<string> & samples,
    map<string,vector<string> > & subgroup2samples)
  {
    gzFile fileStream;
    string line;
    vector<string> tokens;
    for(map<string,string>::const_iterator it = subgroup2file.begin();
        it != subgroup2file.end(); ++it){
      openFile(it->second, fileStream, "rb");
      if(! getline(fileStream, line)){
        cerr << "ERROR: problem with the header of file " << it->second << endl;
        exit(EXIT_FAILURE);
      }
      if(line.empty()){
        cerr << "ERROR: file " << it->second << " is empty" << endl;
        exit(EXIT_FAILURE);
      }
      closeFile(it->second, fileStream);
      if(it == subgroup2file.begin()){
        split(line, " \t", samples);
        if(samples[0] == "Id" || samples[0] == "id" || samples[0] == "ID")
          samples.erase(samples.begin());
        if(! is_unique(samples)){
          cerr << "ERROR: file " << it->second << " has redundant samples"
               << " in its header";
          exit(EXIT_FAILURE);
        }
        subgroup2samples.insert(make_pair(it->first, samples));
      }
      else{
        split(line, " \t", tokens);
        if(tokens[0] == "Id" || tokens[0] == "id" || tokens[0] == "ID")
          tokens.erase(tokens.begin());
        if(! is_unique(tokens)){
          cerr << "ERROR: file " << it->second << " has redundant samples"
               << " in its header";
          exit(EXIT_FAILURE);
        }
        subgroup2samples.insert(make_pair(it->first, tokens));
        for(vector<string>::const_iterator itT = tokens.begin();
            itT != tokens.end(); ++itT)
          if(find(samples.begin(), samples.end(), *itT) == samples.end())
            samples.push_back(*itT);
      }
    }
    if(verbose > 0){
      cout << "nb of samples (" << data_type << "): " << samples.size() << endl
           << flush;
      for(map<string,vector<string> >::const_iterator it =
            subgroup2samples.begin(); it != subgroup2samples.end(); ++it){
        cout << it->first << ": " << it->second.size() << " samples" << endl;
        if(verbose > 1){
          size_t count = 0;
          for(vector<string>::const_iterator itS = it->second.begin();
              itS != it->second.end(); ++itS){
            ++count;
            cout << count << "/" << it->second.size() << " " << *itS << endl;
          }
        }
      }
    }
  }

  void loadSamplesFromGenotypes(
    const map<string, string> & subgroup2genofile,
    const int & verbose,
    vector<string> & samples,
    map<string,vector<string> > & subgroup2samples_genotypes)
  {
    gzFile fileStream;
    string line, sample;
    vector<string> tokens, tokens2;
    size_t i;
    for(map<string,string>::const_iterator it = subgroup2genofile.begin();
        it != subgroup2genofile.end(); ++it){
      openFile(it->second, fileStream, "rb");
      if(! getline(fileStream, line)){
        cerr << "ERROR: problem with the header of file " << it->second << endl;
        exit(EXIT_FAILURE);
      }
      if(line.empty()){
        cerr << "ERROR: file " << it->second << " is empty" << endl;
        exit(EXIT_FAILURE);
      }
      
      // if file in VCF format
      if(line.find("##fileformat=VCF") != string::npos){
        while(getline(fileStream, line)){
          if(line.find("#CHROM") == string::npos)
            continue;
          closeFile(it->second, fileStream);
          split(line, " \t", tokens);
          if(! is_unique(tokens)){
            cerr << "ERROR: file " << it->second << " has redundant samples"
                 << " in its header";
            exit(EXIT_FAILURE);
          }
          subgroup2samples_genotypes.insert(
            make_pair(it->first, vector<string>(tokens.begin()+9, tokens.end())));
          break;
        }
      }
      else{ // not VCF
        split(line, " \t", tokens);
        closeFile(it->second, fileStream);
	
        // if file in IMPUTE format
        if(tokens[0] == "chr"
           && (tokens[1] == "name" || tokens[1] == "id")
           && tokens[2] == "coord"
           && tokens[3] == "a1"
           && tokens[4] == "a2"){
          if((tokens.size() - 5) % 3 != 0){
            cerr << "ERROR: the header of IMPUTE file " << it->second
                 << " is badly formatted" << endl;
            exit(EXIT_FAILURE);
          }
          tokens2.clear();
          i = 5;
          while(i < tokens.size()){
            sample = split(tokens[i], "_a", 0); // sampleX_a1a1, sampleX_a1a2 or sampleX_a2a2
            tokens2.push_back(sample);
            i = i + 3;
          }
          if(! is_unique(tokens2)){
            cerr << "ERROR: file " << it->second << " has redundant samples"
                 << " in its header";
            exit(EXIT_FAILURE);
          }
          subgroup2samples_genotypes.insert(make_pair(it->first, tokens2));
          for(vector<string>::const_iterator itT = tokens2.begin();
              itT != tokens2.end(); ++itT)
            if(find(samples.begin(), samples.end(), *itT) == samples.end())
              samples.push_back(*itT);
        }
        else{ // if file in MatrixEQTL format (allele dosage)
          if(tokens[0] == "Id" || tokens[0] == "id" || tokens[0] == "ID")
            tokens.erase(tokens.begin());
          if(! is_unique(tokens)){
            cerr << "ERROR: file " << it->second << " has redundant samples"
                 << " in its header";
            exit(EXIT_FAILURE);
          }
          subgroup2samples_genotypes.insert(make_pair(it->first, tokens));
        }
      } // if file is not VCF
      
      for(vector<string>::const_iterator itS =
            subgroup2samples_genotypes[it->first].begin();
          itS != subgroup2samples_genotypes[it->first].end(); ++itS){
        if(find(samples.begin(), samples.end(), *itS) == samples.end())
          samples.push_back(*itS);
      }
    } // end for loop over subgroups
    
    if(verbose > 0){
      cout << "nb of samples (genotypes): " << samples.size() << endl << flush;
      for(map<string,vector<string> >::const_iterator it =
            subgroup2samples_genotypes.begin();
          it != subgroup2samples_genotypes.end(); ++it){
        cout << it->first << ": " << it->second.size() << " samples" << endl;
        // if(verbose > 1){
        // 	for(vector<string>::const_iterator itS = it->second.begin();
        // 	     itS != it->second.end(); ++itS)
        // 	  cout << *itS << endl;
        // }
      }
    }
  }
  
  void loadSamples(const map<string,string> & subgroup2genofile,
                   const map<string,string> & subgroup2explevelfile,
                   const map<string,string> & subgroup2covarfile,
                   const string & error_model,
                   const int & verbose,
                   Samples & samples)
  {
    if(verbose > 0)
      cout << "load samples ..." << endl << flush;
    
    vector<string> samples_explevels;
    map<string,vector<string> > subgroup2samples_explevels;
    loadSamplesFromMatrixEqtl(subgroup2explevelfile, "explevels", verbose,
                              samples_explevels, subgroup2samples_explevels);
    
    vector<string> samples_genotypes;
    map<string,vector<string> > subgroup2samples_genotypes;
    loadSamplesFromGenotypes(subgroup2genofile, verbose, samples_genotypes,
                             subgroup2samples_genotypes);
    
    // fill samples by merging samples_explevels and samples_genotypes
    samples.AddSamplesIfNew(samples_explevels);
    samples.AddSamplesIfNew(samples_genotypes);
    if(verbose > 0)
      cout << "total nb of samples: " << samples.GetTotalNbSamples()
           << endl << flush;
    
    vector<string> samples_covariates;
    map<string,vector<string> > subgroup2samples_covariates;
    loadSamplesFromMatrixEqtl(subgroup2covarfile, "covariates", verbose,
                              samples_covariates, subgroup2samples_covariates);
    
    // check that each covariate sample has both geno and pheno
    for(vector<string>::const_iterator it = samples_covariates.begin();
        it != samples_covariates.end(); ++it)
      if(samples.IsAbsent(*it)){
        cerr << "ERROR: sample " << *it << " has covariates but neither expression levels nor genotypes" << endl;
        exit(EXIT_FAILURE);
      }
    
    // do the mapping between the vector of all samples and each subgroup
    samples.AddSamplesFromData(subgroup2samples_genotypes, "genotype");
    samples.AddSamplesFromData(subgroup2samples_explevels, "explevel");
    samples.AddSamplesFromData(subgroup2samples_covariates, "covariate");
    
    if(error_model == "hybrid"){
      if(verbose > 0){
        cout << "nb of samples with genotype and exp level (pairs of subgroups):" << endl;
        samples.ShowPairs(cout);
      }
      if(verbose > 1)
        samples.ShowAllMappings(cerr);
    }
  }
  
  /** \brief Parse the BED file
   */
  void loadGeneInfo(const string & file_genecoords, const int & verbose,
                    map<string,Gene> & gene2object,
                    map<string,vector<Gene*> > & mChr2VecPtGenes)
  {
    if(verbose > 0)
      cout << "load gene coordinates ..." << endl << flush;
    
    vector<string> lines;
    readFile(file_genecoords, lines);
    
    vector<string> tokens;
    for(size_t i = 0; i < lines.size(); ++i){
      split(lines[i], " \t", tokens);
      if(gene2object.find(tokens[3]) != gene2object.end())
        continue; // in case of redundancy
      if(tokens[1] == tokens[2]){
        cerr << "ERROR: start and end coordinates of " << tokens[3]
             << " should be different (at least 1 bp)" << endl;
        exit(1);
      }
      Gene gene(tokens[3], tokens[0], tokens[1], tokens[2]);
      if(tokens.size() >= 6)
        gene.SetStrand(tokens[5]);
      gene2object.insert(make_pair(gene.GetName(), gene));
      if(mChr2VecPtGenes.find(gene.GetChromosome()) == mChr2VecPtGenes.end())
        mChr2VecPtGenes.insert(make_pair(gene.GetChromosome(),
                                         vector<Gene*>()));
      mChr2VecPtGenes[gene.GetChromosome()].push_back(
        &(gene2object[gene.GetName()]));
    }
    
    // sort the genes per chr
    for(map<string,vector<Gene*> >::iterator it = mChr2VecPtGenes.begin();
        it != mChr2VecPtGenes.end(); ++it)
      sort(it->second.begin(), it->second.end(), pt_gene_lt_pt_gene);
    
    if(verbose > 0)
      cout << "total nb of genes with coordinates: " << gene2object.size()
           << endl;
  }
  
  void loadExplevels(const map<string,string> & subgroup2explevelfile,
                     const int & verbose,
                     map<string,Gene> & gene2object)
  {
    if(gene2object.empty())
      return;
    
    if(verbose > 0)
      cout << "load gene expression levels ..." << endl << flush;
    
    gzFile explevelstream;
    string subgroup, explevelfile, line;
    vector<string> tokens;
    size_t nb_samples, nb_lines, nb_genes_tokeep_per_subgroup;
    
    for(map<string,string>::const_iterator it = subgroup2explevelfile.begin();
        it != subgroup2explevelfile.end(); ++it){
      nb_genes_tokeep_per_subgroup = 0;
      nb_lines = 0;
      subgroup = it->first;
      explevelfile = it->second;
      openFile(explevelfile, explevelstream, "rb");
      if(! getline(explevelstream, line)){
        cerr << "ERROR: problem with the header of file " << explevelfile
             << endl;
        exit(EXIT_FAILURE);
      }
      ++nb_lines;
      split(line, " \t", tokens);
      if(tokens[0] == "Id" || tokens[0] == "id" || tokens[0] == "ID")
        nb_samples = tokens.size() - 1;
      else
        nb_samples = tokens.size();
      
      while(getline(explevelstream, line)){
        ++nb_lines;
        split(line, " \t", tokens);
        if(tokens.size() != nb_samples + 1){
          cerr << "ERROR: not enough columns on line " << nb_lines << " of file "
               << explevelfile << " (" << tokens.size() << " != "
               << nb_samples + 1 << ")" << endl;
          exit(EXIT_FAILURE);
        }
        if(gene2object.find(tokens[0]) == gene2object.end())
          continue; // skip if no coordinate
        gene2object[tokens[0]].AddSubgroup(subgroup, tokens.begin()+1,
                                           tokens.end());
        ++nb_genes_tokeep_per_subgroup;
      }
      if(! gzeof(explevelstream)){
        cerr << "ERROR: can't read successfully file " << explevelfile
             << " up to the end" << endl;
        exit(EXIT_FAILURE);
      }
      
      closeFile(explevelfile, explevelstream);
      if(verbose > 0)
        cout << subgroup << " (" << explevelfile << "): "
             << (nb_lines-1) << " genes (to keep: "
             << nb_genes_tokeep_per_subgroup << ")"
             << endl << flush;
    }
    
    map<string,Gene>::iterator it = gene2object.begin();
    while(it != gene2object.end()){
      if(! it->second.HasExplevelsInAtLeastOneSubgroup()){
        if(verbose > 1)
          cerr << "WARNING: skip gene " << it->second.GetName()
               << " because it has no expression level in any subgroup" << endl;
        gene2object.erase(it++);
      }
      else
        ++it;
    }
    
    if(verbose > 0)
      cout << "total nb of genes to analyze: " << gene2object.size() << endl;
    if(verbose > 1){
      vector<string> tmp;
      size_t count = 0;
      for(map<string,Gene>::iterator it = gene2object.begin();
          it != gene2object.end(); ++it){
        ++count;
        cout << count << "/" << gene2object.size() << " " << it->first << ":";
        it->second.GetSubgroupsWithExpLevels(tmp);
        for(size_t i = 0; i < tmp.size(); ++i)
          cout << " " << tmp[i];
        cout << endl;
      }
    }
  }
  
  void loadSnpsToKeep(const string & file_snpstokeep, const int & verbose,
                      set<string> & sSnpsToKeep)
  {
    if(! file_snpstokeep.empty())
    {  
      string line;
      gzFile stream;
      vector<string> tokens;
      size_t nb_lines = 0;
      
      openFile(file_snpstokeep, stream, "rb");
      if(verbose > 0)
        cout <<"load file " << file_snpstokeep << " ..." << endl;
      
      while(getline(stream, line)){
        nb_lines++;
        split(line, " \t,", tokens);
        if(tokens.size() != 1){
          cerr << "ERROR: file " << file_snpstokeep
               << " should have only one column"
               << " at line " << nb_lines << endl;
          exit(EXIT_FAILURE);
        }
        if(tokens[0][0] == '#')
          continue;
        if(sSnpsToKeep.find(tokens[0]) == sSnpsToKeep.end())
          sSnpsToKeep.insert(tokens[0]);
      }
      
      if(! gzeof(stream)){
        cerr << "ERROR: can't read successfully file "
             << file_snpstokeep << " up to the end" << endl;
        exit(EXIT_FAILURE);
      }
      closeFile(file_snpstokeep, stream);
      
      if(verbose > 0)
        cout << "nb of SNPs to keep: " << sSnpsToKeep.size() << endl;
    }
  }
  
  void loadGenosAndSnpInfoFromImpute(
    const map<string,string>::const_iterator & it_subgroup2genofile,
    const set<string> & sSnpsToKeep,
    const map<string, vector<Gene*> > & mChr2VecPtGenes,
    gzFile & genoStream,
    string & line,
    size_t & tot_nb_snps,
    size_t & nb_snps_tokeep_per_subgroup,
    map<string, Snp> & snp2object,
    map<string, vector<Snp*> > & mChr2VecPtSnps)
  {
    string subgroup = it_subgroup2genofile->first,
      genofile = it_subgroup2genofile->second;
    vector<string> tokens;
    size_t nb_lines = 1; // header line already read
    split(line, " \t", tokens); // header line
    if((tokens.size() - 5) % 3 != 0){
      cerr << "ERROR: wrong number of columns on line " << nb_lines
           << " of file " << genofile << endl;
      exit(EXIT_FAILURE);
    }
    size_t nb_samples = static_cast<size_t>((tokens.size() - 5) / 3);
    
    while(getline(genoStream, line)){
      ++nb_lines;
      ++tot_nb_snps;
      split(line, " \t", tokens);
      if(tokens.size() != 3 * nb_samples + 5){
        cerr << "ERROR: not enough columns on line " << nb_lines << " of file "
             << genofile << " (" << tokens.size() << " != "
             << (3 * nb_samples + 5) << ")" << endl;
        exit(EXIT_FAILURE);
      }
      if(mChr2VecPtGenes.find(tokens[0]) == mChr2VecPtGenes.end())
        continue; // no gene on this chromosome
      if(! sSnpsToKeep.empty() && sSnpsToKeep.find(tokens[1])
         == sSnpsToKeep.end())
        continue;
      if(snp2object.find(tokens[1]) == snp2object.end()){
        Snp snp(tokens[1], tokens[0], tokens[2]);
        snp2object.insert(make_pair(snp.GetName(), snp));
      }
      if(snp2object[tokens[1]].HasGenotypes(subgroup)){
        cerr << "ERROR: SNP " << tokens[1] << " is duplicated in file "
             << genofile << endl;
        exit(EXIT_FAILURE);
      }
      snp2object[tokens[1]].AddSubgroup(subgroup, tokens.begin()+5,
                                        tokens.end(), "impute");
      ++nb_snps_tokeep_per_subgroup;
    }
  }
  
  void
  loadGenosAndSnpInfoFromVcf (
    const map<string,string>::const_iterator & it_subgroup2genofile,
    const set<string> & sSnpsToKeep,
    const map<string, vector<Gene*> > & mChr2VecPtGenes,
    gzFile & genoStream,
    string & line,
    size_t & tot_nb_snps,
    size_t & nb_snps_tokeep_per_subgroup,
    map<string, Snp> & snp2object,
    map<string, vector<Snp*> > & mChr2VecPtSnps)
  {
    string subgroup = it_subgroup2genofile->first,
      genofile = it_subgroup2genofile->second;
    vector<string> tokens, tokens2;
    size_t nb_lines = 1; // header line already read
    
    // skip the first lines of meta-data
    while(getline(genoStream, line)){
      ++nb_lines;
      if(line.find("#CHROM") != string::npos)
        break;
    }
    split(line, " \t", tokens);
    size_t nb_samples = tokens.size() - 9;
    
    while(getline(genoStream, line)){
      ++nb_lines;
      ++tot_nb_snps;
      split(line, " \t", tokens);
      if(tokens.size() != nb_samples + 9){
        cerr << "ERROR: not enough columns on line " << nb_lines
             << " of file " << genofile << " (" << tokens.size() << " != "
             << nb_samples + 9 << ")" << endl;
        exit(EXIT_FAILURE);
      }
      if(tokens[8].find("GT") == string::npos){
        cerr << "ERROR: missing GT in 9-th field on line " << nb_lines
             << " of file " << genofile << endl;
        exit(EXIT_FAILURE);
      }
      if(mChr2VecPtGenes.find(tokens[0]) == mChr2VecPtGenes.end())
        continue; // no gene on this chromosome
      if(! sSnpsToKeep.empty() && sSnpsToKeep.find(tokens[2])
         == sSnpsToKeep.end())
        continue;
      if(snp2object.find(tokens[2]) == snp2object.end()){
        Snp snp(tokens[2], tokens[0], tokens[1]);
        snp2object.insert(make_pair(snp.GetName(), snp));
      }
      split(tokens[8], ":", tokens2); // if several fields in column FORMAT
      size_t idx_gt = 0;
      while(idx_gt < tokens2.size()){
        if(tokens2[idx_gt] == "GT")
          break;
        ++idx_gt;
      }
      if(snp2object[tokens[2]].HasGenotypes(subgroup)){
        cerr << "ERROR: SNP " << tokens[2] << " is duplicated in file "
             << genofile << endl;
        exit(EXIT_FAILURE);
      }
      snp2object[tokens[2]].AddSubgroup(subgroup, tokens.begin()+9,
                                        tokens.end(), "vcf", idx_gt);
      ++nb_snps_tokeep_per_subgroup;
    }
  }
  
  void duplicateGenosPerSnpInAllSubgroups(
    const map<string,string> & subgroup2genofile,
    const int & verbose,
    map<string, Snp> & snp2object)
  {
    if(verbose > 0)
      cout << "duplicate genotypes for other subgroups (same file) ..."
           << flush << endl;
    map<string,string>::const_iterator it = subgroup2genofile.begin();
    string first_subgroup = it->first;
    ++it;
    while(it != subgroup2genofile.end()){
      clock_t startTime = clock();
      for(map<string,Snp>::iterator it_snp = snp2object.begin();
          it_snp != snp2object.end(); ++it_snp)
        it_snp->second.DuplicateGenotypesFromFirstSubgroup(first_subgroup,
                                                           it->first);
      if(verbose > 0)
        cout << it->first << " (" << it->second << "): "
             << snp2object.size() << " SNPs duplicated in "
             << fixed << setprecision(2) << getElapsedTime(startTime) << " sec"
             << endl << flush;
      ++it;
    }
  }
  
  void loadGenosAndSnpInfo(
    const map<string, string> & subgroup2genofile,
    const float & min_maf,
    const set<string> & sSnpsToKeep,
    const map<string, vector<Gene*> > & mChr2VecPtGenes,
    const int & verbose,
    map<string, Snp> & snp2object,
    map<string, vector<Snp*> > & mChr2VecPtSnps)
  {
    if(verbose > 0)
      cout << "load genotypes and SNP coordinates ..." << endl << flush;
    
    gzFile genoStream;
    string line;
    size_t tot_nb_snps, nb_snps_tokeep_per_subgroup;
    
    bool same_files = false;
    for(map<string,string>::const_iterator it = subgroup2genofile.begin();
        it != subgroup2genofile.end(); ++it){
      if(it != subgroup2genofile.begin() && it->second.compare(
           subgroup2genofile.begin()->second) == 0){
        same_files = true;
        break; // avoid loading same file several times
      }
      
      clock_t startTime = clock();
      tot_nb_snps = 0;
      nb_snps_tokeep_per_subgroup = 0;
      openFile(it->second, genoStream, "rb");
      if(! getline(genoStream, line)){
        cerr << "ERROR: problem with the header of file " << it->second << endl;
        exit(EXIT_FAILURE);
      }
      
      if(line.find("##fileformat=VCF") != string::npos) // VCF format
        loadGenosAndSnpInfoFromVcf(it, sSnpsToKeep, mChr2VecPtGenes,
                                   genoStream, line, tot_nb_snps,
                                   nb_snps_tokeep_per_subgroup, snp2object,
                                   mChr2VecPtSnps);
      else if(line.find("chr") != string::npos
              && (line.find("name") != string::npos
                  || line.find("id") != string::npos)
              && line.find("coord") != string::npos
              && line.find("a1") != string::npos
              && line.find("a2") != string::npos) // IMPUTE format
        loadGenosAndSnpInfoFromImpute(it, sSnpsToKeep, mChr2VecPtGenes,
                                      genoStream, line, tot_nb_snps,
                                      nb_snps_tokeep_per_subgroup,
                                      snp2object, mChr2VecPtSnps);
      else if(line.compare(0, 2, "id") == 0 || line.compare(0, 2, "Id") == 0
              || line.compare(0, 2, "ID") == 0){
        cerr << "ERROR: file " << it->second
             << " seems to be in the custom format but --scoord is missing"
             << endl;
        exit(EXIT_FAILURE);
      }
      else{
        cerr << "ERROR: can't recognize the format of file " << it->second
             << endl;
        exit(EXIT_FAILURE);
      }
      
      if(! gzeof(genoStream)){
        cerr << "ERROR: can't read successfully file " << it->second
             << " up to the end" << endl;
        exit(EXIT_FAILURE);
      }
      
      closeFile(it->second , genoStream);
      if(verbose > 0)
        cout << it->first << " (" << it->second << "): " << tot_nb_snps
             << " SNPs (to keep: " << nb_snps_tokeep_per_subgroup << ", loaded in "
             << fixed << setprecision(2) << getElapsedTime(startTime) << " sec)"
             << endl << flush;
    }
    
    if(verbose > 0)
      cout << "discard SNPs with missing values ..." << endl << flush;
    map<string,Snp>::iterator it = snp2object.begin();
    while(it != snp2object.end()){
      it->second.EraseIfMissingValuesPerSubgroup();
      if(! it->second.HasGenotypesInAtLeastOneSubgroup()){
        if(verbose > 1)
          cerr << "WARNING: skip SNP " << it->second.GetName()
               << " because it has missing values in each subgroup" << endl;
        snp2object.erase(it++);
      } else
        ++it;
    }
    
    if(min_maf > 0){
      if(verbose > 0)
        cout << "filter SNPs with MAF < " << min_maf << " ..." << endl << flush;
      map<string,Snp>::iterator it = snp2object.begin();
      while(it != snp2object.end()){
        it->second.EraseIfLowMafPerSubgroup(min_maf);
        if(! it->second.HasGenotypesInAtLeastOneSubgroup()){
          if(verbose > 1)
            cerr << "WARNING: skip SNP " << it->second.GetName()
                 << " because it has a low MAF in each subgroup" << endl;
          snp2object.erase(it++);
        } else
          ++it;
      }
    }
    
    if(same_files)
      duplicateGenosPerSnpInAllSubgroups(subgroup2genofile, verbose,
                                         snp2object);
    
    for(map<string, Snp>::const_iterator it = snp2object.begin();
        it != snp2object.end(); ++it){
      if(mChr2VecPtSnps.find(it->second.GetChromosome()) == mChr2VecPtSnps.end())
        mChr2VecPtSnps.insert(make_pair(it->second.GetChromosome(), vector<Snp*>()));
      mChr2VecPtSnps[it->second.GetChromosome()].push_back(&(snp2object[it->second.GetName()]));
    }
    
    // sort the SNPs per chr
    for(map<string, vector<Snp*> >::iterator it = mChr2VecPtSnps.begin();
        it != mChr2VecPtSnps.end(); ++it)
      sort(it->second.begin(), it->second.end(), pt_snp_lt_pt_snp);
    
    if(verbose > 0)
      cout << "nb of SNPs: " << snp2object.size() << endl;
  }
  
  /** \brief Parse the unindexed BED file
   */
  void loadSnpInfo(const string & snpCoordsFile,
                   const set<string> & sSnpsToKeep,
                   const int & verbose,
                   map<string, Snp> & snp2object)
  {
    gzFile snpCoordsStream;
    openFile(snpCoordsFile, snpCoordsStream, "rb");
    string line;
    vector<string> tokens;
    while(getline(snpCoordsStream, line)){
      split (line, " \t", tokens);
      if(! sSnpsToKeep.empty() && sSnpsToKeep.find(tokens[3])
         == sSnpsToKeep.end())
        continue;
      if(snp2object.find(tokens[3]) != snp2object.end())
        continue; // in case of redundancy
      if(tokens[1] == tokens[2]){
        cerr << "ERROR: start and end coordinates of " << tokens[3]
             << " should be different (at least 1 bp)" << endl;
        exit(1);
      }
      Snp snp(tokens[3], tokens[0], tokens[2]);
      snp2object.insert(make_pair(snp.GetName(), snp));
    }
    if(! gzeof(snpCoordsStream))
    {
      cerr << "ERROR: can't read successfully file " << snpCoordsFile
           << " up to the end" << endl;
      exit(EXIT_FAILURE);
    }
    closeFile(snpCoordsFile, snpCoordsStream);
  }
  
  void loadGenos(const map<string, string> & subgroup2genofile,
                 const float & min_maf, const int & verbose,
                 map<string, Snp> & snp2object,
                 map<string, vector<Snp*> > & mChr2VecPtSnps)
  {
    if(snp2object.empty())
      return;
    
    if(verbose > 0)
      cout << "load genotypes (custom format) ..." << endl << flush;
    
    gzFile genoStream;
    string line;
    vector<string> tokens;
    size_t nb_samples, nb_lines, nb_snps_tokeep_per_subgroup;
    
    bool same_files = false;
    for(map<string,string>::const_iterator it = subgroup2genofile.begin();
        it != subgroup2genofile.end(); ++it){
      if(it != subgroup2genofile.begin() && it->second.compare(
           subgroup2genofile.begin()->second) == 0){
        same_files = true;
        break; // avoid loading same file several times
      }
      
      clock_t startTime = clock();
      nb_snps_tokeep_per_subgroup = 0;
      nb_lines = 0;
      openFile(it->second, genoStream, "rb");
      if(! getline(genoStream, line)){
        cerr << "ERROR: problem with the header of file " << it->second << endl;
        exit(EXIT_FAILURE);
      }
      ++nb_lines;
      split(line, " \t", tokens);
      if(tokens[0] == "chr"){
        cerr << "ERROR: don't use --scoord if genotypes in IMPUTE format" << endl;
        exit(1);
      }
      if(tokens[0].find("##fileformat=VCF") != string::npos){
        cerr << "ERROR: don't use --scoord if genotypes in VCF format" << endl;
        exit(1);
      }
      if(tokens[0] == "Id" || tokens[0] == "id" || tokens[0] == "ID")
        nb_samples = tokens.size() - 1;
      else
        nb_samples = tokens.size();
      
      while(getline(genoStream, line)){
        ++nb_lines;
        split(line, " \t", tokens);
        if(tokens.size() != nb_samples + 1){
          cerr << "ERROR: not enough columns on line " << nb_lines << " of file "
               << it->second << " (" << tokens.size() << " != "
               << nb_samples + 1 << ")" << endl;
          exit(EXIT_FAILURE);
        }
        if(snp2object.find(tokens[0]) == snp2object.end())
          continue; // skip if no coordinate
        if(find(tokens.begin(), tokens.end(), "NA") != tokens.end())
          continue; // skip if missing genotype(s)
        if(snp2object[tokens[0]].HasGenotypes(it->first)){
          cerr << "ERROR: SNP " << tokens[0] << " is duplicated in file "
               << it->second << endl;
          exit(EXIT_FAILURE);
        }
        snp2object[tokens[0]].AddSubgroup(it->first, tokens.begin()+1,
                                          tokens.end(), "dose");
        ++nb_snps_tokeep_per_subgroup;
      }
      if(! gzeof(genoStream)){
        cerr << "ERROR: can't read successfully file " << it->second
             << " up to the end" << endl;
        exit(EXIT_FAILURE);
      }
      
      closeFile(it->second, genoStream);
      if(verbose > 0)
        cout << it->first << " (" << it->second << "): " << (nb_lines-1)
             << " SNPs (to keep: " << nb_snps_tokeep_per_subgroup << ", loaded in "
             << fixed << setprecision(2) << getElapsedTime(startTime) << " sec)"
             << endl << flush;
    }
    
    if(verbose > 0)
      cout << "discard SNPs with missing values ..." << endl << flush;
    map<string,Snp>::iterator it = snp2object.begin();
    while(it != snp2object.end()){
      it->second.EraseIfMissingValuesPerSubgroup();
      if(! it->second.HasGenotypesInAtLeastOneSubgroup()){
        if(verbose > 1)
          cerr << "WARNING: skip SNP " << it->second.GetName()
               << " because it has missing values in each subgroup" << endl;
        snp2object.erase(it++);
      } else
        ++it;
    }
    
    if(min_maf > 0){
      if(verbose > 0)
        cout << "filter SNPs with MAF < " << min_maf << " ..." << endl << flush;
      map<string,Snp>::iterator it = snp2object.begin();
      while(it != snp2object.end()){
        it->second.EraseIfLowMafPerSubgroup(min_maf);
        if(! it->second.HasGenotypesInAtLeastOneSubgroup()){
          if(verbose > 1)
            cerr << "WARNING: skip SNP " << it->second.GetName()
                 << " because it has a low MAF in each subgroup" << endl;
          snp2object.erase(it++);
        } else
          ++it;
      }
    }
    
    if(same_files)
      duplicateGenosPerSnpInAllSubgroups(subgroup2genofile, verbose,
                                         snp2object);
    
    for(map<string, Snp>::const_iterator it = snp2object.begin();
        it != snp2object.end(); ++it){
      if(mChr2VecPtSnps.find(it->second.GetChromosome()) == mChr2VecPtSnps.end())
        mChr2VecPtSnps.insert(make_pair(it->second.GetChromosome(), vector<Snp*>()));
      mChr2VecPtSnps[it->second.GetChromosome()].push_back(&(snp2object[it->second.GetName()]));
    }
    
    // sort the SNPs per chr
    for(map<string,vector<Snp*> >::iterator it = mChr2VecPtSnps.begin();
        it != mChr2VecPtSnps.end(); ++it)
      sort(it->second.begin(), it->second.end(), pt_snp_lt_pt_snp);
    
    if(verbose > 0)
      cout << "total nb of SNPs to analyze: " << snp2object.size() << endl;
  }
  
  void loadListCovarFiles(
    const string & file_covarpaths,
    const string & sbgrpToKeep,
    const vector<string> & subgroups,
    map<string, string> & mCovarPaths,
    const int & verbose)
  {
    mCovarPaths = loadTwoColumnFile (file_covarpaths, verbose);
    
    map<string, string>::iterator it = mCovarPaths.begin();
    while(it != mCovarPaths.end())
    {
      if(! sbgrpToKeep.empty() && sbgrpToKeep.compare(it->first) != 0)
      {
        mCovarPaths.erase (it++);
        continue;
      }
      if(find (subgroups.begin(), subgroups.end(), it->first)
         == subgroups.end())
      {
        cerr << "WARNING: skip covariates of subgroup " << it->first
             << " as there is no corresponding genotype nor phenotype files"
             << endl;
        mCovarPaths.erase (it++);
      }
      else
        ++it;
    }
  }
  
  /** \brief Allow to load several covariate files per subgroup.
   */
  void loadListCovarFiles(
    const string & file_covarpaths,
    const string & sbgrpToKeep,
    const vector<string> & subgroups,
    map<string, vector<string> > & mCovarPaths,
    const int & verbose)
  {
    string line;
    gzFile stream;
    vector<string> tokens;
    size_t nb_lines = 0, nbLoadedCovarFiles = 0;
    
    openFile(file_covarpaths, stream, "rb");
    if(verbose > 0)
      cout <<"load file " << file_covarpaths << " ..." << endl;
    
    while(getline(stream, line)){
      nb_lines;
      split(line, " \t,", tokens);
      if(tokens.size() != 2){
        cerr << "ERROR: file " << file_covarpaths
             << " should have only two columns at line " << nb_lines << endl;
        exit(EXIT_FAILURE);
      }
      if(tokens[0][0] == '#')
        continue;
      if(! sbgrpToKeep.empty() && sbgrpToKeep.compare(tokens[0]) != 0)
        continue;
      if(find(subgroups.begin(), subgroups.end(), tokens[0])
         == subgroups.end()){
        cerr << "WARNING: skip covariates of subgroup " << tokens[0]
             << " as there is no corresponding genotype / phenotype files"
             << endl;
        continue;
      }
      if(mCovarPaths.find(tokens[0]) == mCovarPaths.end())
        mCovarPaths.insert(make_pair(tokens[0], vector<string> ()));
      mCovarPaths[tokens[0]].push_back(tokens[1]);
      ++nbLoadedCovarFiles;
    }
    
    if(! gzeof(stream)){
      cerr << "ERROR: can't read successfully file "
           << file_covarpaths << " up to the end" << endl;
      exit(EXIT_FAILURE);
    }
    closeFile(file_covarpaths, stream);
    
    if(verbose > 0)
      cout << "items loaded: " << nbLoadedCovarFiles << " files for "
           << mCovarPaths.size() << " subgroups" << endl;
  }
  
  void loadCovariates(
    const map<string,string> subgroup2covarfile,
    const int & verbose,
    Covariates & covariates)
  {
    if(subgroup2covarfile.empty())
      return;
    
    if(verbose > 0)
      cout << "load covariates ..." << endl << flush;
    
    gzFile covarstream;
    vector<string> tokens;
    string subgroup, covarfile, line;
    size_t nb_samples, nb_lines;
    
    for(map<string,string>::const_iterator it = subgroup2covarfile.begin();
        it != subgroup2covarfile.end(); ++it){
      map<string,vector<string> > covariate2values;
      nb_lines = 0;
      subgroup = it->first;
      covarfile = it->second;
      openFile(covarfile, covarstream, "rb");
      if(! getline(covarstream, line)){
        cerr << "ERROR: problem with the header of file " << covarfile << endl;
        exit(EXIT_FAILURE);
      }
      ++nb_lines;
      split(line, " \t", tokens);
      if(tokens[0] == "Id" || tokens[0] == "id" || tokens[0] == "ID")
        nb_samples = tokens.size() - 1;
      else
        nb_samples = tokens.size();
      
      while(getline(covarstream, line)){
        ++nb_lines;
        split(line, " \t", tokens);
        if(tokens.size() != nb_samples + 1){
          cerr << "ERROR: not enough columns on line " << nb_lines << " of file "
               << covarfile << " (" << tokens.size() << " != "
               << nb_samples + 1 << ")" << endl;
          exit(EXIT_FAILURE);
        }
        covariate2values.insert(
          make_pair(tokens[0], vector<string>(tokens.begin()+1, tokens.end())));
      }
      if(! gzeof(covarstream)){
        cerr << "ERROR: can't read successfully file " << covarfile
             << " up to the end" << endl;
        exit(EXIT_FAILURE);
      }
      
      closeFile(covarfile, covarstream);
      covariates.AddSubgroup(subgroup, covariate2values);
    }
    
    if(verbose > 0)
      for(map<string,string>::const_iterator it = subgroup2covarfile.begin();
          it != subgroup2covarfile.end(); ++it)
        cout << it->first << " (" << it->second << "): "
             << covariates.GetNbCovariates(it->first)
             << " covariates" << endl;
  }
  
  void loadRawInputData(
    const string & file_genopaths,
    const string & file_snpcoords,
    const string & file_exppaths,
    const string & file_genecoords,
    const string & anchor,
    const size_t & radius,
    const float & min_maf,
    const string & file_covarpaths,
    const string & error_model,
    const vector<string> & subgroups_tokeep,
    const set<string> & sSnpsToKeep,
    const int & verbose,
    vector<string> & subgroups,
    Samples & samples,
    map<string,Snp> & snp2object,
    map<string,vector<Snp*> > & mChr2VecPtSnps,
    Covariates & covariates,
    map<string,Gene> & gene2object)
  {
    map<string,string> subgroup2genofile, subgroup2explevelfile,
      subgroup2covarfile;
    loadListsGenoExplevelAndCovarFiles(file_genopaths, file_exppaths,
                                       file_covarpaths, subgroups_tokeep,
                                       error_model, verbose, subgroup2genofile,
                                       subgroup2explevelfile,
                                       subgroup2covarfile, subgroups);
    
    loadSamples(subgroup2genofile, subgroup2explevelfile, subgroup2covarfile,
                error_model, verbose, samples);
    
    loadCovariates(subgroup2covarfile, verbose, covariates);
    
    map<string,vector<Gene*> > mChr2VecPtGenes;
    loadGeneInfo(file_genecoords, verbose, gene2object, mChr2VecPtGenes);
    loadExplevels(subgroup2explevelfile, verbose, gene2object);
    if(gene2object.empty())
      return;
    
    if(file_snpcoords.empty())
      loadGenosAndSnpInfo(subgroup2genofile, min_maf, sSnpsToKeep,
                          mChr2VecPtGenes, verbose, snp2object, mChr2VecPtSnps);
    else{
      loadSnpInfo(file_snpcoords, sSnpsToKeep, gene2object, anchor, radius,
                  verbose, snp2object);
      loadGenos(subgroup2genofile, min_maf, verbose, snp2object, mChr2VecPtSnps);
    }
    if(snp2object.empty())
      return;
  }
  
  void loadListSstatsFile(
    const string & file_sstats,
    const int & verbose,
    map<string,string> & subgroup2sstatsfile)
  {
    string line;
    gzFile stream;
    vector<string> tokens;
    size_t nb_lines = 0;
    
    openFile(file_sstats, stream, "rb");
    if(verbose > 0)
      cout <<"load file " << file_sstats << " ..." << endl;
    
    while(getline(stream, line)){
      ++nb_lines;
      split(line, " \t,", tokens);
      if(tokens.size() != 2){
        cerr << "ERROR: file " << file_sstats 
             << " should have only two columns at line " << nb_lines << endl;
        exit(EXIT_FAILURE);
      }
      if(tokens[0][0] == '#')
        continue;
      if(subgroup2sstatsfile.find(tokens[0]) == subgroup2sstatsfile.end())
        subgroup2sstatsfile.insert(make_pair(tokens[0], tokens[1]));
    }
    
    if(! gzeof(stream)){
      cerr << "ERROR: can't read successfully file "
           << file_sstats << " up to the end" << endl;
      exit(EXIT_FAILURE);
    }
    closeFile(file_sstats, stream);
    
    if(verbose > 0)
      cout << "items loaded: " << subgroup2sstatsfile.size() << " files" << endl;
  }
  
  void fillGeneSnpPairsWithSstats(
    const map<string,string> & subgroup2sstatsfile,
    const int & verbose,
    map<string,Gene> & gene2object,
    map<string,Snp> & snp2object)
  {
    vector<string> lines, tokens;
    map<string,size_t> col2idx;
    col2idx["gene"] = string::npos;
    col2idx["snp"] = string::npos;
    col2idx["n"] = string::npos;
    col2idx["sigmahat"] = string::npos;
    col2idx["betahat.geno"] = string::npos;
    col2idx["sebetahat.geno"] = string::npos;
    Gene * pt_gene = NULL;
    Snp * pt_snp = NULL;
    vector<GeneSnpPair>::iterator it_gsp;
    size_t idx_snp = string::npos;
    
    // loop over subgroups
    for(map<string,string>::const_iterator it_sf = subgroup2sstatsfile.begin();
        it_sf != subgroup2sstatsfile.end(); ++it_sf){
      if(verbose > 0)
        cout << "load summary statistics for subgroup "
             << it_sf->first << " ..." << endl;
      readFile(it_sf->second, lines);
      
      // parse file header
      for(map<string,size_t>::iterator it_c = col2idx.begin();
          it_c != col2idx.end(); ++it_c)
        it_c->second = string::npos;
      split(lines[0], "\t", tokens);
      for(size_t col_id = 0; col_id < tokens.size(); ++col_id)
        if(col2idx.find(tokens[col_id]) != col2idx.end())
          col2idx[tokens[col_id]] = col_id;
      for(map<string,size_t>::const_iterator it_c = col2idx.begin();
          it_c != col2idx.end(); ++it_c)
        if(it_c->second == string::npos){
          cerr << "ERROR: missing " << it_c->second << " in header of "
               << it_sf->second << endl;
          exit(EXIT_FAILURE);
        }
      
      // parse file content
      for(size_t line_id = 1; line_id < lines.size(); ++line_id){
        split(lines[line_id], "\t", tokens);
	
        // get gene (create it if necessary)
        if(gene2object.find(tokens[col2idx["gene"]]) == gene2object.end())
          gene2object.insert(make_pair(tokens[col2idx["gene"]],
                                       Gene(tokens[col2idx["gene"]])));
        pt_gene = &(gene2object[tokens[col2idx["gene"]]]);
	
        // get snp (create it if necessary)
        if(snp2object.find(tokens[col2idx["snp"]]) == snp2object.end())
          snp2object.insert(make_pair(tokens[col2idx["snp"]],
                                      Snp(tokens[col2idx["snp"]])));
        pt_snp = &(snp2object[tokens[col2idx["snp"]]]);
	
        // get gene-snp pair (create it if necessary)
        idx_snp = pt_gene->FindIdxSnp(pt_snp);
        if(idx_snp == string::npos){
          pt_gene->AddCisSnp(pt_snp);
          it_gsp = pt_gene->AddGeneSnpPair(pt_snp->GetName(), "uvlr");
        } else
          it_gsp = pt_gene->FindGeneSnpPair(idx_snp);
	
        it_gsp->SetSstats(it_sf->first,
                          atol(tokens[col2idx["n"]].c_str()),
                          atof(tokens[col2idx["sigmahat"]].c_str()),
                          atof(tokens[col2idx["betahat.geno"]].c_str()),
                          atof(tokens[col2idx["sebetahat.geno"]].c_str()));
	
      } // end of loop over lines
      
    } // end of loop over subgroups
  }
  
  void loadSummaryStats(
    const string & file_sstats,
    const int & verbose,
    vector<string> & subgroups,
    map<string,Gene> & gene2object,
    map<string,Snp> & snp2object)
  {
    map<string,string> subgroup2sstatsfile;
    loadListSstatsFile(file_sstats, verbose, subgroup2sstatsfile);
    
    keys2vec(subgroup2sstatsfile, subgroups);
    
    fillGeneSnpPairsWithSstats(subgroup2sstatsfile, verbose, gene2object,
                               snp2object);
  }
  
  /** \brief Parse the tabix-indexed BED file
   */
  void loadSnpInfo(const string & file_snpcoords,
                   const string & file_snpcoords_idx,
                   const set<string> & sSnpsToKeep,
                   const map<string, Gene> gene2object,
                   const string & anchor,
                   const size_t & radius,
                   const int & verbose,
                   map<string, Snp> & snp2object)
  {
    struct stat stat_bed, stat_tbi;
    stat(file_snpcoords.c_str(), &stat_bed);
    stat(file_snpcoords_idx.c_str(), &stat_tbi);
    if(stat_bed.st_mtime > stat_tbi.st_mtime){
      cerr << "ERROR: index file (" << file_snpcoords_idx
           << ") is older than data file (" << file_snpcoords << ")" << endl;
      exit(EXIT_FAILURE);
    }
    
    tabix_t * t;
    if((t = ti_open(file_snpcoords.c_str(), 0)) == 0){
      cerr << "ERROR: fail to open the data file (tabix)" << endl;
      exit(EXIT_FAILURE);
    }
    if(ti_lazy_index_load(t) < 0){
      cerr << "ERROR: failed to load the index file (tabix)" << endl;
      exit(EXIT_FAILURE);
    }
    
    const char *s;
    ti_iter_t iter;
    vector<string> tokens;
    for(map<string,Gene>::const_iterator it = gene2object.begin();
        it != gene2object.end(); ++it){
      int t_id, t_beg, t_end, len;
      if(ti_parse_region(t->idx, it->second.GetRegionInTabixFormat(anchor, radius).c_str(),
                         &t_id, &t_beg, &t_end) == 0){
        iter = ti_queryi(t, t_id, t_beg, t_end);
        while((s = ti_read(t, iter, &len)) != 0){
	  
          split(string(s), "\t", tokens);
          if(! sSnpsToKeep.empty() && sSnpsToKeep.find(tokens[3])
             == sSnpsToKeep.end())
            continue;
          if(snp2object.find(tokens[3]) != snp2object.end())
            continue; // in case of redundancy
          Snp snp(tokens[3], tokens[0], tokens[2]);
          snp2object.insert(make_pair(snp.GetName(), snp));
	  
        }
        ti_iter_destroy(iter);
      }
    }
    ti_close(t);
  }
  
  /** \brief Parse the BED file (indexed or not)
   */
  void loadSnpInfo(const string & file_snpcoords,
                   const set<string> & sSnpsToKeep,
                   const map<string,Gene> gene2object,
                   const string & anchor,
                   const size_t & radius,
                   const int & verbose,
                   map<string, Snp> & snp2object)
  {
    if(verbose > 0)
      cout << "load SNP coordinates";
    clock_t startTime = clock();
    
    stringstream file_snpcoords_idx;
    file_snpcoords_idx << file_snpcoords << ".tbi";
    if(doesFileExist(file_snpcoords_idx.str())){
      cout << " (tabix-indexed BED file) ..." << endl << flush;
      loadSnpInfo(file_snpcoords, file_snpcoords_idx.str(), sSnpsToKeep,
                  gene2object, anchor, radius, verbose, snp2object);
    }
    else{
      cout << " (unindexed BED file) ..." << endl << flush;
      loadSnpInfo(file_snpcoords, sSnpsToKeep, verbose, snp2object);
    }
    
    if(verbose > 0)
      cout << "total nb of SNPs with coordinates: " << snp2object.size()
           << " (loaded in " << fixed << setprecision(2)
           << getElapsedTime(startTime) << " sec)" << endl;
  }
  
} // namespace quantgen
