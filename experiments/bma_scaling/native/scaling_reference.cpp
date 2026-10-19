// Offline audit only. Streaming extension of the independently checked audit
// enumerator: retain regional masses and selected score checks, not every model.
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <chrono>
#include <cstdint>
using namespace Rcpp;
double scaling_logadd(double a,double b) {
  if(!std::isfinite(a)) return b;
  double m=std::max(a,b);return m+std::log1p(std::exp(-std::abs(a-b)));
}
// [[Rcpp::export]]
List scaling_prefix_log_masses(NumericVector leaflogs,int h) {
  if(h<1||h>12)stop("Invalid number of groups");
  int expected=1;for(int j=0;j<h;++j)expected*=3;
  if(leaflogs.size()!=expected)stop("Invalid leaf table");
  List levels(h+1);levels[h]=leaflogs;
  NumericVector previous=leaflogs;
  for(int depth=h-1;depth>=0;--depth) {
    int size=previous.size()/3;NumericVector values(size);
    for(int j=0;j<size;++j)
      values[j]=scaling_logadd(scaling_logadd(previous[j],previous[j+size]),previous[j+2*size]);
    levels[depth]=values;previous=values;
  }
  return levels;
}
// [[Rcpp::export]]
List scaling_exact_reference(NumericMatrix gram,NumericVector xy,double yty,int n,
 double sigma,double tau2,double pi,IntegerVector group,IntegerVector check_masks,double seconds=900) {
  int p=gram.nrow();
  if(p<1||p>24||gram.ncol()!=p||xy.size()!=p||group.size()!=p||
     sigma<=0||tau2<=0||pi<=0||pi>=1) stop("Invalid audit input");
  int h=max(group)+1;
  if(h<1||h>12||min(group)<0)stop("Audit supports at most 12 groups");
  std::vector<int> powers(h,1);int leaves=1;
  for(int j=0;j<h;++j){powers[j]=leaves;leaves*=3;}
  int total=1<<p;
  for(int i=0;i<check_masks.size();++i)
    if(check_masks[i]<0||check_masks[i]>=total||(i&&check_masks[i]<=check_masks[i-1]))
      stop("Check masks must be sorted, unique, and in range");
  NumericVector leaflogs(leaves,R_NegInf),checked(check_masks.size());
  IntegerVector counts(leaves);
  std::vector<int> active(p),gc(h);std::vector<double>L(p*p),v(p);
  double constant=-.5*n*std::log(2*M_PI*sigma*sigma)-yty/(2*sigma*sigma);
  double lp=std::log(pi),lq=std::log1p(-pi),logZ=R_NegInf;
  int check_position=0;
  auto began=std::chrono::steady_clock::now();
  for(int mask=0;mask<total;++mask) {
    if(mask%1024==0) {
      Rcpp::checkUserInterrupt();
      if(std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()>seconds)
        stop("Audit time limit");
    }
    int k=0;std::fill(gc.begin(),gc.end(),0);
    for(int j=0;j<p;++j)if(mask&(1<<j)){active[k++]=j;gc[group[j]]++;}
    double ld=0,gain=0;
    for(int i=0;i<k;++i) {
      for(int j=0;j<=i;++j) {
        double a=(i==j?1.0:0.0)+tau2*gram(active[i],active[j]);
        for(int z=0;z<j;++z)a-=L[i*p+z]*L[j*p+z];
        if(i==j) {
          if(!(a>0))stop("Non-positive audit Cholesky pivot");
          L[i*p+j]=std::sqrt(a);ld+=2*std::log(L[i*p+j]);
        } else L[i*p+j]=a/L[j*p+j];
      }
      double a=xy[active[i]];for(int j=0;j<i;++j)a-=L[i*p+j]*v[j];
      v[i]=a/L[i*p+i];gain+=v[i]*v[i];
    }
    double value=constant-.5*ld+tau2*gain/(2*sigma*sigma)+k*lp+(p-k)*lq;
    logZ=scaling_logadd(logZ,value);
    int code=0;for(int j=0;j<h;++j)code+=std::min(2,gc[j])*powers[j];
    counts[code]++;leaflogs[code]=scaling_logadd(leaflogs[code],value);
    if(check_position<check_masks.size()&&mask==check_masks[check_position])checked[check_position++]=value;
  }
  // Avoid millions of sequential log-additions determining normalization.
  // Regional aggregation has far fewer terms; accumulate them in long double.
  double shift=max(leaflogs);long double sum=0;
  for(int j=0;j<leaves;++j)sum+=std::exp((long double)leaflogs[j]-shift);
  double stable_logZ=shift+(double)std::log(sum);
  return List::create(_["logZ"]=stable_logZ,_["sequential_logZ"]=logZ,_["leaf_log_mass"]=leaflogs,
    _["leaf_model_count"]=counts,_["check_masks"]=check_masks,_["check_log_scores"]=checked,
    _["enumerated_models"]=total,
    _["seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count());
}
