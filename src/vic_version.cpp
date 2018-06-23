#include <vic_R.h>

//' Display version of VIC model and this package.
//'
//'
//' @export
// [[Rcpp::export]]
void vic_version() {
  Rprintf("VIC Driver      : %s\n", VIC_DRIVER);
  Rprintf("Package version : 0.1.0\n");
  Rprintf("VIC Version     : %s\n", VERSION);
  print_license();
}

