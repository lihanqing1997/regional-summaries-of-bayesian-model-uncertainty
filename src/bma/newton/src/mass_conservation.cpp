// Deterministic propagation of supplied log-mass endpoints through a tree.
// MPFR encloses the arithmetic; validity of the supplied endpoints is external.
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <mpfr.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
namespace {
struct Bound {
  mpfr_t x;
  Bound(){mpfr_init2(x,160);mpfr_set_zero(x,0);}
  Bound(const Bound& a){mpfr_init2(x,160);mpfr_set(x,a.x,MPFR_RNDN);}
  ~Bound(){mpfr_clear(x);}
};
void from_log(Bound& x,double v,mpfr_rnd_t rounding){
  if(v==R_NegInf){mpfr_set_zero(x.x,0);return;}
  if(v==R_PosInf){mpfr_set_inf(x.x,1);return;}
  mpfr_set_d(x.x,v,MPFR_RNDN);mpfr_exp(x.x,x.x,rounding);
}
double to_log(Bound& x,mpfr_rnd_t rounding){
  if(mpfr_zero_p(x.x))return R_NegInf;
  if(mpfr_inf_p(x.x))return R_PosInf;
  Bound y;mpfr_log(y.x,x.x,rounding);return mpfr_get_d(y.x,rounding);
}
}

// Each children entry is a COMPLETE disjoint partition of the named parent.
// Stopping before a fixed point preserves every enclosure already obtained.
// [[Rcpp::export]]
List ars_propagate_log_mass(NumericVector log_lower,NumericVector log_upper,
                           IntegerVector parents,List children,int max_passes=4){
  const int n=log_lower.size(),m=parents.size();
  if(log_upper.size()!=n || children.size()!=m || max_passes<1)
    stop("Invalid mass-conservation dimensions");
  std::vector<Bound> lo(n),hi(n);
  std::vector<std::vector<int>> kids(m);
  for(int i=0;i<n;++i){
    if(std::isnan(log_lower[i]) || std::isnan(log_upper[i]) ||
       log_lower[i]>log_upper[i] || log_lower[i]==R_PosInf)
      stop("Invalid mass-conservation endpoints");
    from_log(lo[i],log_lower[i],MPFR_RNDD);from_log(hi[i],log_upper[i],MPFR_RNDU);
  }
  for(int k=0;k<m;++k){
    if(parents[k]<1 || parents[k]>n)stop("Invalid parent index");
    IntegerVector ids=children[k];
    if(ids.size()==0)stop("Empty child partition");
    for(int i:ids){if(i<1 || i>n || i==parents[k])stop("Invalid child index");kids[k].push_back(i-1);}
  }
  bool changed=false;int passes=0;long updates=0;
  auto tighten=[&](int i,Bound& a,Bound& b){
    if(mpfr_cmp(a.x,lo[i].x)>0){mpfr_set(lo[i].x,a.x,MPFR_RNDN);changed=true;++updates;}
    if(mpfr_cmp(b.x,hi[i].x)<0){mpfr_set(hi[i].x,b.x,MPFR_RNDN);changed=true;++updates;}
    if(mpfr_cmp(lo[i].x,hi[i].x)>0)stop("Mass-conservation interval conflict");
  };
  for(int pass=0;pass<max_passes;++pass){
    Rcpp::checkUserInterrupt();changed=false;++passes;
    // A bottom-up pass first transmits child evidence to ancestors.
    for(int k=m-1;k>=0;--k){
      Bound a,b;
      for(int i:kids[k]){mpfr_add(a.x,a.x,lo[i].x,MPFR_RNDD);mpfr_add(b.x,b.x,hi[i].x,MPFR_RNDU);}
      tighten(parents[k]-1,a,b);
    }
    // A top-down pass transmits ancestor evidence and sibling restrictions.
    for(int k=0;k<m;++k){
      int p=parents[k]-1;
      for(int i:kids[k]){
        Bound a,b,sl,su;
        for(int j:kids[k])if(i!=j){mpfr_add(sl.x,sl.x,lo[j].x,MPFR_RNDD);mpfr_add(su.x,su.x,hi[j].x,MPFR_RNDU);}
        if(!mpfr_inf_p(su.x)){
          mpfr_sub(a.x,lo[p].x,su.x,MPFR_RNDD);
          if(mpfr_sgn(a.x)<0)mpfr_set_zero(a.x,0);
        }
        if(mpfr_inf_p(hi[p].x))mpfr_set_inf(b.x,1);
        else mpfr_sub(b.x,hi[p].x,sl.x,MPFR_RNDU);
        tighten(i,a,b);
      }
    }
    if(!changed)break;
  }
  NumericVector lower(n),upper(n);int changed_nodes=0;
  for(int i=0;i<n;++i){
    lower[i]=std::max(log_lower[i],to_log(lo[i],MPFR_RNDD));
    upper[i]=std::min(log_upper[i],to_log(hi[i],MPFR_RNDU));
    if(lower[i]>upper[i])stop("Mass-conservation interval conflict after conversion");
    if(lower[i]>log_lower[i] || upper[i]<log_upper[i])++changed_nodes;
  }
  return List::create(_["log_lower"]=lower,_["log_upper"]=upper,
    _["passes"]=passes,_["updates"]=static_cast<double>(updates),
    _["changed_nodes"]=changed_nodes,_["fixed_point"]=!changed);
}
