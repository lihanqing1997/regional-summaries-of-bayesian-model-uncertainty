// Conservative verification of an externally ranked prefix; no optimization.
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <mpfr.h>
#include <vector>
#include <cmath>
using namespace Rcpp;
struct Number {
  mpfr_t x;
  Number(){mpfr_init2(x,256);mpfr_set_zero(x,0);}
  Number(const Number& o){mpfr_init2(x,256);mpfr_set(x,o.x,MPFR_RNDN);}
  ~Number(){mpfr_clear(x);}
};
// [[Rcpp::export]]
IntegerVector verified_ranked_prefix(NumericVector lo, NumericVector hi,
                                    IntegerVector ranking, double epsilon) {
  int n=lo.size(),m=ranking.size();
  if(hi.size()!=n || epsilon<=0 || epsilon>=1)stop("Invalid inputs");
  if(m==0)return IntegerVector();
  double shift=R_NegInf;
  std::vector<int> rank(n,m+1);
  for(int k=0;k<m;++k){int i=ranking[k]-1;
    if(i<0||i>=n||rank[i]!=m+1)stop("Invalid ranking");rank[i]=k;}
  for(int i=0;i<n;++i){
    if(!std::isfinite(hi[i]) || std::isnan(lo[i]) || lo[i]>hi[i])stop("Invalid bounds");
    shift=std::max(shift,hi[i]);
  }
  std::vector<Number> lower(n),upper(n);
  Number ep,q,tmp;
  mpfr_set_d(ep.x,epsilon,MPFR_RNDN);
  mpfr_ui_sub(q.x,1,ep.x,MPFR_RNDU);
  for(int i=0;i<n;++i){
    if(std::isfinite(lo[i])){
      mpfr_set_d(lower[i].x,lo[i],MPFR_RNDN);
      mpfr_sub_d(lower[i].x,lower[i].x,shift,MPFR_RNDD);
      mpfr_exp(lower[i].x,lower[i].x,MPFR_RNDD);
    }
    mpfr_set_d(upper[i].x,hi[i],MPFR_RNDN);
    mpfr_sub_d(upper[i].x,upper[i].x,shift,MPFR_RNDU);
    mpfr_exp(upper[i].x,upper[i].x,MPFR_RNDU);
  }
  auto passes=[&](int k){
    Number retained,omitted,a,b;
    for(int i=0;i<n;++i){
      if(rank[i]<k)mpfr_add(retained.x,retained.x,lower[i].x,MPFR_RNDD);
      else mpfr_add(omitted.x,omitted.x,upper[i].x,MPFR_RNDU);
    }
    mpfr_mul(a.x,ep.x,retained.x,MPFR_RNDD);
    mpfr_mul(b.x,q.x,omitted.x,MPFR_RNDU);
    return mpfr_cmp(a.x,b.x)>=0;
  };
  if(!passes(m))return IntegerVector();
  int left=1,right=m;
  while(left<right){int middle=left+(right-left)/2;
    if(passes(middle))right=middle;else left=middle+1;}
  return ranking[Range(0,left-1)];
}
