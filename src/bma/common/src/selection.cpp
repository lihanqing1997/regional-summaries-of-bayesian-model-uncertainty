// Directed-rounding interval dynamic program for the original mass constraint.
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <mpfr.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
using namespace Rcpp;
static thread_local mpfr_prec_t PREC=128;
struct M {
  mpfr_t x;
  M(){mpfr_init2(x,PREC);mpfr_set_zero(x,0);}
  M(const M& o){mpfr_init2(x,mpfr_get_prec(o.x));mpfr_set(x,o.x,MPFR_RNDN);}
  M(M&& o) noexcept {mpfr_init2(x,mpfr_get_prec(o.x));mpfr_swap(x,o.x);}
  M& operator=(const M& o){if(this!=&o){mpfr_set_prec(x,mpfr_get_prec(o.x));mpfr_set(x,o.x,MPFR_RNDN);}return *this;}
  ~M(){mpfr_clear(x);}
};
void endpoint(M& out,double logv,double shift,mpfr_rnd_t rnd){
  if(!std::isfinite(logv)){mpfr_set_zero(out.x,0);return;}
  mpfr_set_d(out.x,logv,MPFR_RNDN);
  mpfr_sub_d(out.x,out.x,shift,rnd);mpfr_exp(out.x,out.x,rnd);
}
IntegerVector onebased(const std::vector<int>& s){IntegerVector x(s.size());for(size_t j=0;j<s.size();++j)x[j]=s[j]+1;return x;}

// The interval sums enclose the mathematical exponentials of the supplied
// floating-point log endpoints. Input-bound derivations remain the caller's job.
// [[Rcpp::export]]
List ars_interval_knapsack(NumericVector log_lower,NumericVector log_upper,
                          IntegerVector costs,LogicalVector eligible,double epsilon,
                          bool conservative=true,int maximum_bits=16384,
                          double memory_mb=192){
  const int n=costs.size();std::vector<int> ids;
  double largest=R_NegInf,smallest=R_PosInf;
  for(int i=0;i<n;++i){
    if(eligible[i])ids.push_back(i);
    largest=std::max(largest,log_upper[i]);
    if(std::isfinite(log_lower[i]))smallest=std::min(smallest,log_lower[i]);
  }
  if(ids.empty())return List::create(_["feasible"]=false,_["lower_cost"]=R_PosInf,
      _["upper_cost"]=R_PosInf,_["selected"]=IntegerVector(),_["status"]="infeasible",_["bits"]=0);
  int cheapest=costs[ids[0]];for(int i:ids)cheapest=std::min(cheapest,costs[i]);
  double span=std::isfinite(smallest)?largest-smallest:0;
  int bits=static_cast<int>(std::min(static_cast<double>(maximum_bits),
                std::max(128.0,std::ceil(span/std::log(2.0))+96.0)));
  for(;;){
    Rcpp::checkUserInterrupt();PREC=bits;
    std::vector<M> ll(n),lh(n),ul(n),uh(n),bl(n),bh(n);
    M ep,q,tl,th,tmp,tmp2;
    mpfr_set_d(ep.x,epsilon,MPFR_RNDN);mpfr_ui_sub(q.x,1,ep.x,MPFR_RNDN);
    for(int i=0;i<n;++i){
      endpoint(ll[i],log_lower[i],largest,MPFR_RNDD);
      endpoint(lh[i],log_lower[i],largest,MPFR_RNDU);
      endpoint(ul[i],log_upper[i],largest,MPFR_RNDD);
      endpoint(uh[i],log_upper[i],largest,MPFR_RNDU);
      // Conservative benefit q*U + ep*L; optimistic benefit q*L + ep*U.
      M &al=conservative?ul[i]:ll[i], &ah=conservative?uh[i]:lh[i];
      M &cl=conservative?ll[i]:ul[i], &ch=conservative?lh[i]:uh[i];
      mpfr_mul(bl[i].x,q.x,al.x,MPFR_RNDD);mpfr_mul(tmp.x,ep.x,cl.x,MPFR_RNDD);
      mpfr_add(bl[i].x,bl[i].x,tmp.x,MPFR_RNDD);
      mpfr_mul(bh[i].x,q.x,ah.x,MPFR_RNDU);mpfr_mul(tmp.x,ep.x,ch.x,MPFR_RNDU);
      mpfr_add(bh[i].x,bh[i].x,tmp.x,MPFR_RNDU);
      mpfr_add(tl.x,tl.x,al.x,MPFR_RNDD);mpfr_add(th.x,th.x,ah.x,MPFR_RNDU);
    }
    mpfr_mul(tl.x,tl.x,q.x,MPFR_RNDD);mpfr_mul(th.x,th.x,q.x,MPFR_RNDU);
    // Direct test avoids cancellation between the large selected upper masses
    // on both sides of the transformed covering inequality.
    auto direct=[&](const std::vector<int>& chosen){
      std::vector<bool> yes(n,false);for(int i:chosen)yes[i]=true;
      M retained_lo,retained_hi,omitted_lo,omitted_hi,a,b;
      for(int i=0;i<n;++i){
        if(yes[i]){
          mpfr_add(retained_lo.x,retained_lo.x,(conservative?ll[i]:ul[i]).x,MPFR_RNDD);
          mpfr_add(retained_hi.x,retained_hi.x,(conservative?lh[i]:uh[i]).x,MPFR_RNDU);
        }else{
          mpfr_add(omitted_lo.x,omitted_lo.x,(conservative?ul[i]:ll[i]).x,MPFR_RNDD);
          mpfr_add(omitted_hi.x,omitted_hi.x,(conservative?uh[i]:lh[i]).x,MPFR_RNDU);
        }
      }
      mpfr_mul(a.x,ep.x,retained_lo.x,MPFR_RNDD);mpfr_mul(b.x,q.x,omitted_hi.x,MPFR_RNDU);
      if(mpfr_cmp(a.x,b.x)>=0)return 1;
      mpfr_mul(a.x,ep.x,retained_hi.x,MPFR_RNDU);mpfr_mul(b.x,q.x,omitted_lo.x,MPFR_RNDD);
      if(mpfr_cmp(a.x,b.x)<0)return -1;
      return 0;
    };
    int all=direct(ids);
    if(all<0)return List::create(_["feasible"]=false,_["lower_cost"]=R_PosInf,
      _["upper_cost"]=R_PosInf,_["selected"]=IntegerVector(),_["status"]="infeasible",_["bits"]=bits);
    if(all==0){
      if(bits<maximum_bits){bits=std::min(maximum_bits,bits*2);continue;}
      return List::create(_["feasible"]=false,_["lower_cost"]=cheapest,
        _["upper_cost"]=R_PosInf,_["selected"]=IntegerVector(),_["status"]="precision-unresolved",_["bits"]=bits);
    }
    if(mpfr_zero_p(th.x)) {
      int j=ids[0];for(int i:ids)if(costs[i]<costs[j])j=i;
      std::vector<int> singleton{j};
      if(direct(singleton)==1)return List::create(_["feasible"]=true,
        _["lower_cost"]=costs[j],_["upper_cost"]=costs[j],
        _["selected"]=onebased(singleton),_["status"]="optimal",_["bits"]=bits);
    }
    // Obtain a directly verified feasible cost cap; either greedy ordering may
    // be approximate, but the cap is used only after exact directed verification.
    std::vector<int> chosen=ids;long long cap=0;for(int i:chosen)cap+=costs[i];
    for(int rule=0;rule<2;++rule){
      std::vector<int> order=ids;
      std::sort(order.begin(),order.end(),[&](int a,int b){
        M aa=bl[a],bb=bl[b];
        if(rule){mpfr_div_ui(aa.x,aa.x,costs[a],MPFR_RNDN);mpfr_div_ui(bb.x,bb.x,costs[b],MPFR_RNDN);}
        int cmp=mpfr_cmp(aa.x,bb.x);return cmp==0?a<b:cmp>0;
      });
      M cumulative;std::vector<int> prefix;long long cost=0;
      for(int i:order){prefix.push_back(i);cost+=costs[i];mpfr_add(cumulative.x,cumulative.x,bl[i].x,MPFR_RNDD);
        if(mpfr_cmp(cumulative.x,th.x)>=0 && direct(prefix)==1){
          if(cost<cap){cap=cost;chosen=prefix;}break;
        }
      }
    }
    if(cap>INT_MAX-1)stop("RESOURCE_LIMIT: integer cost dimension too large");
    size_t width=static_cast<size_t>(cap)+1;
    double bytes=ids.size()*static_cast<double>(width)+2.0*width*(sizeof(M)+bits/8.0+32.0);
    if(bytes>memory_mb*1024*1024)return List::create(_["feasible"]=true,
      _["lower_cost"]=cheapest,_["upper_cost"]=static_cast<double>(cap),
      _["selected"]=onebased(chosen),_["status"]="resource-bounded",_["bits"]=bits);
    std::vector<M> bestlo(width),besthi(width);
    std::vector<unsigned char> take(ids.size()*width,0);
    for(size_t k=0;k<ids.size();++k){
      if(k%16==0)Rcpp::checkUserInterrupt();int i=ids[k],c=costs[i];
      for(int b=static_cast<int>(cap);b>=c;--b){
        mpfr_add(tmp.x,bestlo[b-c].x,bl[i].x,MPFR_RNDD);
        if(mpfr_cmp(tmp.x,bestlo[b].x)>0){mpfr_set(bestlo[b].x,tmp.x,MPFR_RNDN);take[k*width+b]=1;}
        mpfr_add(tmp.x,besthi[b-c].x,bh[i].x,MPFR_RNDU);
        if(mpfr_cmp(tmp.x,besthi[b].x)>0)mpfr_set(besthi[b].x,tmp.x,MPFR_RNDN);
      }
    }
    int lb=-1,ub=-1;
    for(int b=cheapest;b<=cap;++b){
      if(lb<0 && mpfr_cmp(besthi[b].x,tl.x)>=0)lb=b;
      if(ub<0 && mpfr_cmp(bestlo[b].x,th.x)>=0)ub=b;
    }
    if(lb<0)lb=cheapest; // Numerical uncertainty can weaken, never raise, this bound.
    std::vector<int> selected=chosen;int uppercost=static_cast<int>(cap);
    if(ub>=0){
      std::vector<int> trial;int b=ub;
      for(int k=static_cast<int>(ids.size())-1;k>=0;--k)if(take[k*width+b]){
        trial.push_back(ids[k]);b-=costs[ids[k]];
      }
      std::reverse(trial.begin(),trial.end());
      if(!trial.empty() && direct(trial)==1){selected=trial;uppercost=0;for(int i:selected)uppercost+=costs[i];}
    }
    if(lb==uppercost || bits>=maximum_bits)return List::create(_["feasible"]=true,
      _["lower_cost"]=lb,_["upper_cost"]=uppercost,_["selected"]=onebased(selected),
      _["status"]=lb==uppercost?"optimal":"precision-bounded",_["bits"]=bits);
    bits=std::min(maximum_bits,bits*2);
  }
}
