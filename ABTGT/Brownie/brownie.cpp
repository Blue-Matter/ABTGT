#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () 
{
  
  DATA_VECTOR(N_rel);
  DATA_VECTOR(N_NR);
  DATA_MATRIX(N_recap);
  DATA_INTEGER(n_group);
  DATA_INTEGER(n_y);
  DATA_INTEGER(latency);
  
  PARAMETER(log_M);
  PARAMETER_VECTOR(log_F);
  PARAMETER(log_retain); // Proportion of tags retained from year to year
  PARAMETER(log_report_rate);
  PARAMETER_VECTOR(log_latent_report);
  
  matrix<Type> N_pred(n_group, n_y);
  matrix<Type> C_pred(n_group, n_y);
  C_pred.setZero();
  N_pred.setZero();
  
  vector<Type> F(n_y);
  vector<Type> Z(n_y);
  Type M = exp(log_M);
  for(int j=0;j<n_y;j++) {
    F(j) = exp(log_F(j));
    Z(j) = F(j) + M;
  }
  Type retain = exp(log_retain);
  Type report_rate = exp(log_report_rate);
  vector<Type> latent_report(log_latent_report.size());
  for(int ii=0;ii<log_latent_report.size();ii++) latent_report(ii) = exp(log_latent_report(ii));
  
  matrix<Type> surv(n_group, n_y);
  matrix<Type> surv_cum(n_group, n_y);
  
  for(int i=0;i<n_group;i++) { // i indexes release year
    for(int j=0;j<n_y;j++) { // j indexes recapture year
      if(i == j) {
        surv(i,j) = 1;
        surv_cum(i,j) = 1;
      } else if(i < j) {
        int ii = i + latency;
        if(ii > j) {
          surv(i,j) = exp(-latent_report(j-i) * F(j-1) - M) * retain;
        } else {
          surv(i,j) = exp(-F(j-1) - M) * retain;
        }
        surv_cum(i,j) = 1;
      }
      
      if(i < j) {
        for(int jj=i;jj<j;jj++) surv_cum(i,j) *= surv(i,jj); //Double checked in R
      }
      
      
      if(i<=j) {
        N_pred(i,j) = N_rel(i) * surv_cum(i,j);
        int ii = i + latency;
        if(ii > j) {
          Type ZZ = latent_report(j-i) * F(j) + M;
          C_pred(i,j) = latent_report(j-i) * F(j)/ZZ * (1 - exp(-ZZ)) * N_pred(i,j) * report_rate;
        } else {
          C_pred(i,j) = F(j)/Z(j) * (1 - exp(-Z(j))) * N_pred(i,j) * report_rate;
        }
      }
    }
  }
  
  vector<Type> C_pred_i(n_group);
  C_pred_i.setZero();
  matrix<Type> P_pred(n_group, n_y);
  
  vector<Type> NLL(n_group);
  NLL.setZero();
  
  for(int i=0;i<n_group;i++) {
    for(int j=0;j<n_y;j++) {
      if(i <= j) C_pred_i(i) += C_pred(i,j);
    }
    for(int j=0;j<n_y;j++) {
      if(i <= j) {
        P_pred(i,j) = C_pred(i,j)/N_rel(i);
        NLL(i) -= N_recap(i,j) * log(P_pred(i,j));
      }
    }
    NLL(i) -= N_NR(i) * log(1 - C_pred_i(i)/N_rel(i));
  }
  
  REPORT(N_pred);
  REPORT(P_pred);
  REPORT(surv);
  REPORT(surv_cum);
  REPORT(C_pred);
  REPORT(C_pred_i);
  REPORT(NLL);
  REPORT(P_pred);
  ADREPORT(M);
  ADREPORT(F);
  ADREPORT(retain);
  ADREPORT(latent_report);
  REPORT(M);
  REPORT(F);
  REPORT(latent_report);
  
  Type fn = NLL.sum();
  return fn;
}
