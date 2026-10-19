// Audit only. Enumerate all supports after the practical search, retaining
// regional sufficient summaries rather than a million-row support matrix.
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <chrono>
using namespace Rcpp;
double logadd(double a,double b){if(!std::isfinite(a))return b;double m=std::max(a,b);return m+std::log1p(std::exp(-std::abs(a-b)));}
// [[Rcpp::export]]
List pilot_exact_reference(NumericMatrix gram,NumericVector xy,double yty,int n,
 double sigma,double tau2,double pi,IntegerVector group,double seconds=180){
  int p=gram.nrow();if(p>20||p<1||gram.ncol()!=p||xy.size()!=p||group.size()!=p)stop("Invalid audit dimension");
  int h=max(group)+1;std::vector<int> power(h,1);int leaves=1;
  for(int j=0;j<h;++j){power[j]=leaves;leaves*=3;}
  int total=1<<p;NumericVector logs(total),leaflogs(leaves,R_NegInf);
  IntegerVector codes(total),counts(leaves);NumericMatrix moments(leaves,p);
  std::vector<int> active(p),gc(h);std::vector<double>L(p*p),v(p);
  double constant=-.5*n*std::log(2*M_PI*sigma*sigma)-yty/(2*sigma*sigma);
  double lp=std::log(pi),lq=std::log1p(-pi),logZ=R_NegInf;
  auto began=std::chrono::steady_clock::now();
  for(int mask=0;mask<total;++mask){
    if(mask%1024==0){Rcpp::checkUserInterrupt();
      if(std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()>seconds)stop("Audit time limit");}
    int k=0;std::fill(gc.begin(),gc.end(),0);
    for(int j=0;j<p;++j)if(mask&(1<<j)){active[k++]=j;gc[group[j]]++;}
    double ld=0,gain=0;
    for(int i=0;i<k;++i){
      for(int j=0;j<=i;++j){double a=(i==j?1.0:0.0)+tau2*gram(active[i],active[j]);
        for(int z=0;z<j;++z)a-=L[i*p+z]*L[j*p+z];
        if(i==j){if(!(a>0))stop("Non-positive audit Cholesky pivot");L[i*p+j]=std::sqrt(a);ld+=2*std::log(L[i*p+j]);}
        else L[i*p+j]=a/L[j*p+j];
      }
      double a=xy[active[i]];for(int j=0;j<i;++j)a-=L[i*p+j]*v[j];
      v[i]=a/L[i*p+i];gain+=v[i]*v[i];
    }
    double value=constant-.5*ld+tau2*gain/(2*sigma*sigma)+k*lp+(p-k)*lq;
    logs[mask]=value;logZ=logadd(logZ,value);
    int code=0;for(int j=0;j<h;++j)code+=std::min(2,gc[j])*power[j];
    codes[mask]=code;counts[code]++;leaflogs[code]=logadd(leaflogs[code],value);
  }
  for(int mask=0;mask<total;++mask){double w=std::exp(logs[mask]-logZ);int code=codes[mask];
    for(int j=0;j<p;++j)if(mask&(1<<j))moments(code,j)+=w;
  }
  return List::create(_["logZ"]=logZ,_["leaf_log_mass"]=leaflogs,
    _["leaf_model_count"]=counts,_["leaf_pip_contribution"]=moments,
    _["support_log_scores"]=logs,_["enumerated_models"]=total,
    _["seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count());
}
