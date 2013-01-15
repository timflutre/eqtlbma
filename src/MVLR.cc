/** \file MVLR.cc
 *
 *  `MVLR' is a class implementing the multivariate linear regression
 *  Copyright (C) 2012-2013 Xioaquan Wen
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

using namespace std;

#include "MVLR.h"
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>


void MVLR::init(vector<vector<double> > & Y_in, vector<vector<double> > & Xg_in, vector<vector<double> > & Xc_in){
  
  s = Y_in.size();
  n = Y_in[0].size();
  p = Xg_in.size();
  q = Xc_in.size();
  
  // for the intercept 
  q++;

  Y  = gsl_matrix_calloc(n,s);
  Xg = gsl_matrix_calloc(n,p);
  Xc = gsl_matrix_calloc(n,q);
  
  for(int i=0;i<s;i++){
    for(int j=0;j<n;j++){
      gsl_matrix_set(Y, j,i,Y_in[i][j]);
    }
  }    
  

  for(int i=0;i<p;i++){
    for(int j=0;j<n;j++){
      gsl_matrix_set(Xg,j,i,Xg_in[i][j]);
    }
  }

  if(q>1){
    for(int i=1;i<q;i++){
      for(int j=0;j<n;j++){
	gsl_matrix_set(Xc,j,i,Xc_in[i-1][j]);
      }
    }
  }
  
  for(int j=0;j<n;j++){
    gsl_matrix_set(Xc,j,0,1.0);
  }

  // default value for IW prior on Sigma H = diag(1e-4), m = s-1
  m = s-1;
  H = gsl_matrix_calloc(s,s);
  for(int i=0;i<s;i++){
    gsl_matrix_set(H,i,i,1e-4);
  }

  
  T = eVb = eVg_inv = Gamma = Wg = Sigma = Sigma_inv = Sigma0 = Sigma0_inv=0;  
  compute_common(); 
  compute_Sigma_null();

}


// ================ Setting parmaeters/options ======================== //


void MVLR::set_effect_vec(const vector<double> &phi2_vec_in , const vector<double> &omg2_vec_in){ 

  phi2_vec = phi2_vec_in;
  omg2_vec = omg2_vec_in;
  
}


void MVLR::set_IW_prior(gsl_matrix *H_in, int m_in){
  if(H!=0)
    gsl_matrix_free(H);
  H = H_in;
  m = m_in;
}


// =============== Core Computation ================================ //


void MVLR::compute_common(){
 
  gsl_matrix *XctXc = gsl_matrix_calloc(q,q);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Xc,Xc,0,XctXc);
  gsl_matrix *XctXc_inv = gsl_matrix_calloc(q,q);
 
  if(q==1)
    gsl_matrix_set(XctXc_inv,0,0,1.0/gsl_matrix_get(XctXc,0,0));
  else{
    gsl_matrix *t = gsl_matrix_calloc(q,q);
    gsl_matrix_memcpy(t,XctXc);
    int ss;
    gsl_permutation * pp = gsl_permutation_alloc(q);
    gsl_linalg_LU_decomp (t, pp, &ss);
    gsl_linalg_LU_invert (t, pp, XctXc_inv);
    gsl_permutation_free(pp);
    gsl_matrix_free(t);
  }
  
  
     
  gsl_matrix *t1 = gsl_matrix_calloc(q,n);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,XctXc_inv,Xc,0,t1);
  
  gsl_matrix *t2 = gsl_matrix_calloc(n,n);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Xc,t1,0,t2);
  
  T = gsl_matrix_calloc(n,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      double v = -gsl_matrix_get(t2,i,j);
      if(i==j)
	v += 1;
      gsl_matrix_set(T,i,j,v);
    }
  }
  
  gsl_matrix_free(XctXc);
  gsl_matrix_free(XctXc_inv);
  gsl_matrix_free(t1);
  gsl_matrix_free(t2);
  
}



void MVLR::compute_Sigma(vector<vector<int> >& config){

  // use the null Sigma estimated initially
  
  if(sigma_option < 1e-6){
    Sigma = Sigma0;
    Sigma_inv = Sigma0_inv;
    return;
  }
  
  // else

  if(Sigma!=0){
    gsl_matrix_free(Sigma);
    Sigma = 0;
  }
  
  if(Sigma_inv !=0){
    gsl_matrix_free(Sigma_inv);
    Sigma_inv = 0;
  }

  compute_Sigma_mle(config);
  return;
  
}


void MVLR::compute_Sigma_null(){
  
  // compute S0, the under the null
  gsl_matrix *t1 = gsl_matrix_calloc(n,s);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,T,Y,0,t1);
  
  Sigma0 = gsl_matrix_calloc(s,s);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Y,t1,0,Sigma0);
  
  gsl_matrix_scale(Sigma0,1.0/(n+m-q-s-1));
  Sigma = Sigma0;

  invert_Sigma();
  
  Sigma0_inv = Sigma_inv;
  Sigma = Sigma_inv = 0;
  
  gsl_matrix_free(t1);
  
}





void MVLR::compute_Sigma_mle(vector<vector<int> >& config){
  

  gsl_matrix *E = gsl_matrix_calloc(n,s);
  gsl_vector *cpv = gsl_vector_calloc(n);
  vector<double> fac_vec;
  
  for(int i=0;i<s;i++){
    
    vector<int> index_vec;
    for(int j=0;j<p;j++){
      if(config[j][i] == 1){
	index_vec.push_back(j);
      }
    }
   

    int msize = q+index_vec.size();

    
    gsl_matrix *yv = gsl_matrix_calloc(n,1);
    gsl_matrix_get_col(cpv,Y,i);
    gsl_matrix_set_col(yv,0,cpv);
    
    gsl_matrix *Xv = gsl_matrix_calloc(n,msize);
    for(int j=0;j<q;j++){
      gsl_matrix_get_col(cpv,Xc,j);
      gsl_matrix_set_col(Xv,j,cpv);
    }
    
    for(size_t j=0;j<index_vec.size();j++){
      gsl_matrix_get_col(cpv,Xg,index_vec[j]);
      gsl_matrix_set_col(Xv,q+j,cpv);
    }

       
    double factor = 1;
    gsl_matrix *resv = compute_residual(yv,Xv,msize,factor);
    fac_vec.push_back(sqrt(factor));
    

    gsl_matrix_get_col(cpv,resv,0);
    gsl_matrix_set_col(E,i,cpv);
    
    gsl_matrix_free(yv);
    gsl_matrix_free(Xv);
    gsl_matrix_free(resv);
    
  }
  
  gsl_matrix *S =gsl_matrix_calloc(s,s);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,E,E,0,S);
  gsl_matrix_scale(S,double(1.0)/double(m+n));

  
  Sigma = gsl_matrix_calloc(s,s);
  gsl_matrix_memcpy(Sigma,H);
  gsl_matrix_scale(Sigma,double(m)/double(m+n));
  
  gsl_matrix_add(Sigma,S);
  
      
  gsl_matrix *C = gsl_matrix_calloc(s,s);
  for(int i=0;i<s;i++){
    gsl_matrix_set(C,i,i,fac_vec[i]);
  }

  gsl_matrix *t1 = gsl_matrix_calloc(s,s);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,C,Sigma,0,t1);

  gsl_matrix *t2 = gsl_matrix_calloc(s,s);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Sigma,C,0,t2);
  gsl_matrix_memcpy(Sigma,t2);
  

  gsl_matrix *t3 = gsl_matrix_calloc(s,s);
  gsl_matrix_memcpy(t3,Sigma0);
  gsl_matrix_scale(t3,1.0-sigma_option);
  gsl_matrix_scale(Sigma,sigma_option);
  gsl_matrix_add(Sigma,t3);



  

  gsl_matrix_free(C);
  gsl_matrix_free(t1);
  gsl_matrix_free(t2);
  
  gsl_matrix_free(S);
  
  gsl_matrix_free(E);
  gsl_vector_free(cpv);
  gsl_matrix_free(t3);

  invert_Sigma();

}
    

void MVLR::invert_Sigma(){
  
  if(Sigma_inv !=0){
    gsl_matrix_free(Sigma_inv);
  }
  
  int ss;
  gsl_permutation * pp = gsl_permutation_alloc(s);
  gsl_matrix *t1 = gsl_matrix_calloc(s,s);
  
  gsl_matrix_memcpy(t1,Sigma);
  gsl_linalg_LU_decomp (t1, pp, &ss);
  Sigma_inv = gsl_matrix_calloc(s,s);
  gsl_linalg_LU_invert (t1, pp, Sigma_inv);
  
  gsl_permutation_free(pp);
  gsl_matrix_free(t1);
}




gsl_matrix *MVLR::compute_residual(gsl_matrix *y, gsl_matrix *X, int size, double &factor){
  
  gsl_matrix *XtX = gsl_matrix_calloc(size, size);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,X,X,0,XtX);
    
  // compute inverse of XtX (generalized inverse version)
  gsl_matrix *V = gsl_matrix_calloc(size,size);
  gsl_vector *S = gsl_vector_calloc(size);
  gsl_vector *work = gsl_vector_calloc(size);
  gsl_linalg_SV_decomp (XtX, V, S,work);
  
  gsl_matrix *t1 = gsl_matrix_calloc(size,size);
  for(int i=0;i<size;i++){
    double v = gsl_vector_get(S,i);
    if(v>1e-8){
      gsl_matrix_set(t1,i,i,1/v);
    }
  }
  gsl_matrix *t2 = gsl_matrix_calloc(size,size);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,t1,0,t2);
  
  gsl_matrix *XtX_inv = gsl_matrix_calloc(size,size);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,t2,V,0,XtX_inv);
  
  
  // (X'X)^{-1)X'
  gsl_matrix *t3 = gsl_matrix_calloc(size,n);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,XtX_inv,X,0,t3);
  
  gsl_matrix *hB = gsl_matrix_calloc(size,1);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t3,y,0,hB);

  gsl_matrix *fy = gsl_matrix_calloc(n,1);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,X,hB,0,fy);
  
  gsl_matrix *res= gsl_matrix_calloc(n,1);
  gsl_matrix_memcpy(res,y);
  gsl_matrix_sub(res,fy);
  
    
  
  //comput correction factor
  if(size>q){
    
    gsl_matrix *t4 = gsl_matrix_calloc(1,1);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,res,res,0,t4);
    double sigma1 = gsl_matrix_get(t4,0,0)/(n-size);
      
    gsl_matrix *t5 = gsl_matrix_calloc(n,size-q);
    for(int i=0;i<size-q;i++){
      for(int j=0;j<n;j++){
	gsl_matrix_set(t5,j,i,gsl_matrix_get(X,j,i+q));
      }
    }
    
    gsl_matrix *t6 = gsl_matrix_calloc(size-q,n);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,t5,T,0,t6);
    
    gsl_matrix *t7 = gsl_matrix_calloc(size-q,size-q);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t6,t5,0,t7);
    
    gsl_matrix *t8 = gsl_matrix_calloc(size-q,1);
    for(int i=0;i<size-q;i++){
      gsl_matrix_set(t8,i,0,gsl_matrix_get(hB,q+i,0));
    }
    
    gsl_matrix *t9 = gsl_matrix_calloc(1,size-q);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,t8,t7,0,t9);
    
    gsl_matrix *t10 = gsl_matrix_calloc(1,1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t9,t8,0,t10);
    

    // hotelling's T2
    double T2 = gsl_matrix_get(t10,0,0)/pow(sigma1,2);
    

    // convert to F
    double v1 = size-q;
    double v2 = n-size;
    double F = (v2-v1+1)*T2/(v1*v2);
    double qv = gsl_cdf_fdist_Q(F, v1, v2-v1+1);
    double new_F = gsl_cdf_chisq_Qinv (qv,v1)/v1;
    if(F<1e-8)
      factor = 1;
    else
      factor = F/new_F;


    gsl_matrix_free(t4);
    gsl_matrix_free(t5);
    gsl_matrix_free(t6);
    gsl_matrix_free(t7);
    gsl_matrix_free(t8);
    gsl_matrix_free(t9);
    gsl_matrix_free(t10);

  }
    
  gsl_matrix_free(t1);
  gsl_matrix_free(t2);
  gsl_matrix_free(t3);
  
  gsl_matrix_free(XtX);
  gsl_matrix_free(XtX_inv);
  
  gsl_matrix_free(hB);
  gsl_matrix_free(fy);
  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_vector_free(work);

  return res;
  
}


// construct priors according to skeleton information

void MVLR::construct_Gamma(vector<vector<int> >& indicator, vector<int> &noz_vec){
    
  if(Gamma != 0){
    gsl_matrix_free(Gamma);
    Gamma = 0;
  }
  
  Gamma = gsl_matrix_calloc(ep*s,ep*s);
  
  if(prior_option == 1)
    construct_meta_Gamma(indicator,noz_vec);
  if(prior_option == 2)
    construct_diag_Gamma(indicator,noz_vec);

}



void MVLR::construct_diag_Gamma(vector<vector<int> >& indicator, vector<int> &noz_vec){

  for(int i=0;i<ep;i++){
    int index = noz_vec[i];
    for(int k=0;k<s;k++){
      gsl_matrix_set(Gamma,i*s+k,i*s+k,indicator[index][k]*gsl_matrix_get(Sigma,k,k));
    }
  }
  return;

}


void MVLR::construct_meta_Gamma(vector<vector<int> >& indicator, vector<int> &noz_vec){
    
  for(int i=0;i<ep;i++){ 
    int index = noz_vec[i];
    gsl_matrix *t1 = gsl_matrix_calloc(s,1);
    gsl_matrix *t2 = gsl_matrix_calloc(s,s);
    
    for(int j=0;j<s;j++){
      gsl_matrix_set(t1,j,0,indicator[index][j]*sqrt(gsl_matrix_get(Sigma,j,j)));    
    }
    
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,t1,t1,0,t2);
    for(int j=0;j<s;j++){
      for(int k=0;k<s;k++){
	gsl_matrix_set(Gamma,i*s+j,i*s+k,gsl_matrix_get(t2,j,k));
      }
    }
    gsl_matrix_free(t1);
    gsl_matrix_free(t2);
  }
  return;
}


void MVLR::set_Wg(double phi2, double omg2){
  
  if(Wg != 0){
    gsl_matrix_free(Wg);
    Wg = 0;
  }
  
  Wg = gsl_matrix_calloc(ep*s,ep*s);  
  gsl_matrix_memcpy(Wg, Gamma);
  gsl_matrix_scale(Wg, omg2);
  
  for(int j=0;j<ep*s;j++){
    gsl_matrix_set(Wg,j,j,(omg2+phi2)*gsl_matrix_get(Gamma,j,j));
  }
  
}




void MVLR::compute_stats(vector<int>& noz_vec){ 
  
  gsl_matrix *eXg = gsl_matrix_calloc(n,ep);
  gsl_vector *mv = gsl_vector_calloc(n);
  for(int i=0;i<ep;i++){
    int index = noz_vec[i];
    gsl_matrix_get_col(mv,Xg,index);
    gsl_matrix_set_col(eXg,i,mv);
  }

  
  gsl_matrix *G = gsl_matrix_calloc(n,ep);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,T,eXg,0,G);


  gsl_matrix *eKg = gsl_matrix_calloc(ep,ep);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,G,G,0,eKg);
  
  
  eVg_inv = kron(eKg,Sigma_inv, ep,s);
  
  gsl_matrix *t1 = gsl_matrix_calloc(s,n);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,Sigma_inv,Y,0,t1);
  
  gsl_matrix *t2 = gsl_matrix_calloc(s,ep);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t1,G,0,t2);
  
  eVb = vec(t2,s,ep);
    
  
  gsl_matrix_free(eXg);
  gsl_vector_free(mv);
  gsl_matrix_free(G);
  gsl_matrix_free(eKg);
  gsl_matrix_free(t1);
  gsl_matrix_free(t2);
 
}




double MVLR::compute_log10_ABF(gsl_matrix *Wg){
   
  // 1. I+Vg^-1Wg
  gsl_matrix *t1 = gsl_matrix_calloc(ep*s,ep*s);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,eVg_inv,Wg,0,t1);
  
  for(int i=0;i<ep*s;i++){
    gsl_matrix_set(t1,i,i,gsl_matrix_get(t1,i,i)+1);
  }
 
  
  // 2.  (I+Vg^-1Wg)^-1  and determinant
  int ss;
  gsl_permutation *pp = gsl_permutation_alloc(ep*s);
  gsl_linalg_LU_decomp (t1, pp, &ss);
  double log_detVal = gsl_linalg_LU_lndet(t1);  
  gsl_matrix *t2 = gsl_matrix_calloc(ep*s,ep*s);
  gsl_linalg_LU_invert (t1, pp, t2);
  
  
  //3. Wg(I+Vg^-1Wg)^-1
  gsl_matrix *t3 = gsl_matrix_calloc(ep*s,ep*s);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Wg,t2,0,t3);
  
  //4. quadratic form
  gsl_matrix *t4 = gsl_matrix_calloc(1,ep*s);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,eVb,t3,0,t4);
  gsl_matrix *t5 = gsl_matrix_calloc(1,1);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t4,eVb,0,t5);
  
  double rst = .5*gsl_matrix_get(t5,0,0);


  rst += -0.5*log_detVal;
  //printf("log_detVal = %f\n",log_detVal*(-0.5));

  gsl_permutation_free(pp);

  gsl_matrix_free(t1);
  gsl_matrix_free(t2);
  gsl_matrix_free(t3);
  gsl_matrix_free(t4);
  gsl_matrix_free(t5);

  return rst/log(10.0);

}


  




double MVLR::compute_log10_ABF(vector<vector<int> >& indicator){
  
    
  // construct the skeleton from the indicator
  vector<int> noz_vec; // non_zero prior SNP index
  
  for(int i=0;i<p;i++){
    int flag = 0;
    for(int j=0;j<s;j++){
      if(indicator[i][j]!=0){
	flag=1;
	break;
      }
    } 
    if(flag==1)
      noz_vec.push_back(i);
  }
  
  ep = noz_vec.size();
  if(ep == 0)
    return 0.0;


  // 1. first compute Sigma
  compute_Sigma(indicator);
  // 2. pre-compute necessary pieces common to all (phi2,omg2)
  compute_stats(noz_vec);

  // 2. scaled skeleton cov
  construct_Gamma(indicator,noz_vec); 
   
  vector<double> rst_vec;
  vector<double> wts_vec;
  
  int size = omg2_vec.size();
  
  for(int i=0;i<size;i++){
    set_Wg(phi2_vec[i], omg2_vec[i]);
    rst_vec.push_back(compute_log10_ABF(Wg));
    wts_vec.push_back(1.0/double(size));
    gsl_matrix_free(Wg);
    Wg = 0;
  }
  
  gsl_matrix_free(eVb);
  gsl_matrix_free(eVg_inv);
  gsl_matrix_free(Gamma);
  eVb=eVg_inv=Gamma=0;

  return log10_weighted_sum(rst_vec, wts_vec);

} 
		      




double MVLR::compute_log10_ABF(vector<vector<int> >& indicator, double phi2, double omg2){
  
    
  // construct the skeleton from the indicator
  vector<int> noz_vec; // non_zero prior SNP index
  
  for(int i=0;i<p;i++){
    int flag = 0;
    for(int j=0;j<s;j++){
      if(indicator[i][j]!=0){
	flag=1;
	break;
      }
    } 
    if(flag==1)
      noz_vec.push_back(i);
  }
  
  ep = noz_vec.size();
  if(ep == 0)
    return 0.0;


  // 1. first compute Sigma
  compute_Sigma(indicator);
  // 2. pre-compute necessary pieces common to all (phi2,omg2)
  compute_stats(noz_vec);

  // 2. scaled skeleton cov
  construct_Gamma(indicator,noz_vec); 
   
  set_Wg(phi2, omg2);
  double rst = compute_log10_ABF(Wg);
  
  gsl_matrix_free(Wg);
  gsl_matrix_free(eVb);
  gsl_matrix_free(eVg_inv);
  gsl_matrix_free(Gamma);
  Wg=eVb=eVg_inv=Gamma=0;

  return rst;

} 


vector<double> MVLR::compute_log10_ABF_vec(vector<vector<int> >& indicator){
  
  vector<double> rst_vec;
  
    
  // construct the skeleton from the indicator
  vector<int> noz_vec; // non_zero prior SNP index
  
  for(int i=0;i<p;i++){
    int flag = 0;
    for(int j=0;j<s;j++){
      if(indicator[i][j]!=0){
	flag=1;
	break;
      }
    } 
    if(flag==1)
      noz_vec.push_back(i);
  }
  
  ep = noz_vec.size();
  if(ep == 0)
    return rst_vec;


  // 1. first compute Sigma
  compute_Sigma(indicator);
  // 2. pre-compute necessary pieces common to all (phi2,omg2)
  compute_stats(noz_vec);

  // 2. scaled skeleton cov
  construct_Gamma(indicator,noz_vec); 
   
  int size = omg2_vec.size();
  
  for(int i=0;i<size;i++){
    set_Wg(phi2_vec[i], omg2_vec[i]);
    rst_vec.push_back(compute_log10_ABF(Wg));
    gsl_matrix_free(Wg);
    Wg = 0;
  }
  
  gsl_matrix_free(eVb);
  gsl_matrix_free(eVg_inv);
  gsl_matrix_free(Gamma);
  eVb=eVg_inv=Gamma=0;

  return rst_vec;

} 
		      




// utility


double MVLR::log10_weighted_sum(vector<double> &vec, vector<double> &wts){

    double max = vec[0];
    for(size_t i=0;i<vec.size();i++){
      if(vec[i]>max)
	max = vec[i];
    }
    double sum = 0;
    for(size_t i=0;i<vec.size();i++){
      sum += wts[i]*pow(10, (vec[i]-max));
    }

    return (max+log10(sum));
}




gsl_matrix *MVLR::vec (gsl_matrix *M, int a, int b){
  
  gsl_matrix *v = gsl_matrix_calloc(a*b, 1);
  int k = 0;
  for(int j=0;j<b;j++){
    for(int i=0;i<a;i++){
      gsl_matrix_set(v,k++,0,gsl_matrix_get(M,i,j));
    }
  }
  
  return v;
}


gsl_matrix *MVLR::kron (gsl_matrix *M, gsl_matrix *L, int a, int b){
  
  gsl_matrix *R = gsl_matrix_calloc (a*b,a*b);
  for (int i = 0; i < a; i++){
    for (int j = 0; j < a; j++){ 
	  for (int k = 0; k < b; k++){
	      for (int l = 0; l < b; l++){
		  gsl_matrix_set (R, i*b+k,j*b+l, gsl_matrix_get (M, i, j)*gsl_matrix_get (L, k, l));
	      }
	  }
    }
  }
  return R;
}

gsl_matrix *MVLR::kron2(gsl_matrix *A, int nrows, int ncols,gsl_matrix *B, int mrows, int mcols){
  
  gsl_matrix *C = gsl_matrix_calloc(nrows*mrows, ncols*mcols);
  
  int cr = 0;
  int cc = 0;
  for (int i = 0; i<nrows;i++){
    for(int j=0;j<ncols;j++){
      cr = i*mrows;
      cc = j*mcols;
      double a = gsl_matrix_get(A,i,j);
      
      for(int k=0;k<mrows;k++){
	for(int l=0;l<mcols;l++){
	  gsl_matrix_set(C,cr+k,cc+l,a*gsl_matrix_get(B,k,l));
	}
      }
      
    }
    
  }
  
  return C;
}




void MVLR::print_matrix(gsl_matrix *M, int a, int b){
  
  for(int i=0;i<a;i++){
    for(int j=0;j<b;j++){
      //printf("%f  ",gsl_matrix_get(M,i,j));
      printf("%e  ",gsl_matrix_get(M,i,j));
    }
    printf("\n");
  }
  printf("\n");
}



MVLR::~MVLR(){
  
  
  gsl_matrix_free(Y);
  gsl_matrix_free(Xg);
  gsl_matrix_free(Xc);
  gsl_matrix_free(H);
  
  if(T!=0)
    gsl_matrix_free(T);
  
  if(Wg!=0){
    gsl_matrix_free(Wg);
  }
    
  
  if(Sigma!=0)
    gsl_matrix_free(Sigma);
  
  if(Sigma_inv != 0)
    gsl_matrix_free(Sigma_inv);
  

  if(eVg_inv!=0)
    gsl_matrix_free(eVg_inv);  
  
  
}




