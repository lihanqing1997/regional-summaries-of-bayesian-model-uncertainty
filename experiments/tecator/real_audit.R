# The native audit has no Cholesky ridge; the inherited R score adds 1e-10.
# Use the same scale-relative criterion as node audits for the real-data score
# comparison, recording the absolute discrepancy in every reference object.
real_reference <- function(d) {
  f<-make_production_reference
  replace_check<-function(x){
    if(identical(x,quote(discrepancy>1e-8)))
      return(quote(discrepancy>1e-9*max(1,abs(ans$support_log_scores[masks+1L]))))
    if(is.call(x))for(i in seq_along(x)){
      y<-replace_check(x[[i]]);if(!is.null(y))x[[i]]<-y
    }
    x
  }
  body(f)<-replace_check(body(f))
  f(d)
}
