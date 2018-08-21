#include <vic_R.h>

// A auxiraly function for the routing of Lohmann model
//
// [[Rcpp::export]]
NumericVector aux_Lohmann_conv(NumericMatrix tmpm) {
  int len = tmpm.rows();
  int cols = tmpm.cols();
  int i, j, clen;
  NumericVector Q(len, 0.);

  for(i = 0; i < len; i++) {
    clen = i + 1;
    if (clen > cols) clen = cols;
    for(j = 0; j < clen; j ++) {
      Q[i] += tmpm(i-j, j);
    }
  }
  return Q;
}
