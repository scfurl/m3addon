// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp17")]]


// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//


namespace Internal {

// struct IndexValue{
//   int index;
//   int value;
//   IndexValue(int ind, int val){
//     index = ind;
//     value = val;
//   }
// };

arma::sp_mat sp_submatrix(const arma::sp_mat *A, std::vector<std::size_t> *rows, std::vector<std::size_t> *cols) {
  
  std::size_t n_rows = rows->size();
  std::size_t n_cols = cols->size();
  
  bool found = false;
  std::size_t n = 0;
  std::size_t p = 0;
  std::size_t found_idx = 0;
  
  arma::vec new_val(A->n_nonzero);
  arma::uvec new_row_ind(A->n_nonzero);
  arma::uvec new_col_ptr(n_cols + 1);
  
  new_col_ptr(p) = 0;
  for (auto const& j: *cols) { // for every column in the cols vector
    for (std::size_t k = A->col_ptrs[j]; k < A->col_ptrs[j + 1]; k++) {  // k is the index of the "values" and "row_indices" that corresponds to the column j
      // search row_ind[k] in rows
      found = false;
      found_idx = 0;
      while (!found && found_idx < n_rows) {
        if (A->row_indices[k] == rows->at(found_idx))
          found = true;
        found_idx++;
      }
      // store the values if the row was found in rows
      if (found) { // if the row index is in the designated rows...
        new_val(n) = A->values[k]; // store the value
        new_row_ind(n) = found_idx - 1;  // store the index where the original index was found inside "rows"
        n++;
      }
    }
    p++;
    new_col_ptr(p) = n;
  }
  new_col_ptr(p) = n ;
  // reshape the vectors to the actual number of elements
  new_val.reshape(n, 1);
  new_row_ind.reshape(n, 1);
  return arma::sp_mat(new_row_ind, new_col_ptr, new_val, n_rows, n_cols);
}

double vst(double q, const double eP, const double aD){
  double r = std::log( (1 + eP + 2 * aD * q + 2 * std::sqrt(aD * q * (1 + eP + aD * q)))/(4 * aD))/std::log(2);
  return r;
}


}


// [[Rcpp::export]]
Rcpp::NumericMatrix sp_show_storage(const arma::sp_mat& x) {
  Rcpp::NumericMatrix ans(x.n_nonzero, 3u);
  int i = 0;
  for(arma::sp_mat::const_iterator it = x.begin(); it != x.end(); ++it) {
    ans(i, 0) = it.row(); // Row position
    ans(i, 1) = it.col(); // Col position
    ans(i++, 2) = *it;    // Value
  }
  // Adding colnames
  colnames(ans) = Rcpp::CharacterVector::create("row", "col", "val");
  return ans;
}


// [[Rcpp::export]]
arma::sp_mat subsetSM(const arma::sp_mat& a, Rcpp::NumericVector rind, Rcpp::NumericVector cind) {
  std::vector<size_t> rind_out = Rcpp::as<std::vector<size_t>>(rind);
  std::vector<size_t> cind_out = Rcpp::as<std::vector<size_t>>(cind);
  return Internal::sp_submatrix(&a, &rind_out, &cind_out);
}

// [[Rcpp::export]]
arma::vec calculateScoresCPP(const arma::sp_mat& a, const Rcpp::NumericVector Rgenevec) {
  arma::sp_mat A = a;
  Rcpp::NumericVector Cgenevec = Rgenevec - 1;
  A.for_each( [](arma::sp_mat::elem_type& val) { val = std::log10(val + 1); });
  std::vector<size_t> cind_out = Rcpp::as<std::vector<size_t> >(Cgenevec);
  std::vector<size_t> rvec(A.n_rows);
  std::iota(rvec.begin(), rvec.end(), 0);
  arma::sp_mat outmat = Internal::sp_submatrix(&A, &rvec, &cind_out) ;
  arma::vec svec( sum( outmat, 1) );
  arma::vec outvec = svec/A.n_rows;
  //arma::sp_mat submat = A.cols(cind);
  return outvec;
}

// [[Rcpp::export]]
arma::sp_mat logNormCPP(const arma::sp_mat& a) {
  arma::sp_mat A = a;
  A.for_each( [](arma::sp_mat::elem_type& val) { val = std::log10(val + 1); });
  return A;
}

// [[Rcpp::export]]
arma::mat colStdDev(arma::sp_mat& a) {
  arma::mat D(a);
  arma::mat V = stddev(D, 0, 0);
  return V;
}


// [[Rcpp::export]]
arma::mat rowStdDev(arma::sp_mat& a) {
  arma::mat D(a.t());
  arma::mat V = stddev(D, 0, 0);
  return V;
}


// [[Rcpp::export]]
arma::mat vstExprsCPP(const arma::sp_mat& a, const double eP, const double aD, bool return_sparse = true, bool calc_zeros = false) {
  arma::sp_mat A = a;
  A.for_each( [eP, aD](arma::sp_mat::elem_type& val) { val = Internal::vst(val, eP, aD); });
  arma::mat D(A);
  if (calc_zeros) {
    double r = Internal::vst(0, eP, aD);
    D.replace(0, r);
  }
  return D;
}


// [[Rcpp::export]]
void iterateSparseMatrixCPP(const arma::sp_mat& a) {
  //mat<-as.matrix(log10((Matrix::t(Matrix::t(Biobase:::exprs(cds))/pData(cds)$Size_Factor)) + 1))
  arma::sp_mat::const_iterator it     = a.begin();
  arma::sp_mat::const_iterator it_end = a.end();
  for(; it != it_end; ++it)
  {
    arma::cout << "val: " << (*it)    << arma::endl;
    arma::cout << "row: " << it.row() << arma::endl;
    arma::cout << "col: " << it.col() << arma::endl;
  }
  ///arma::mat Y(mana);
  //Rcpp::Rcout << "DenseMat res:\n" << Y.print("Y:") << std::endl;
  //Rcpp::Rcout << "DenseMat res:\n" << std::endl;
  return;
}


#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::export]]
Rcpp::NumericVector computeSparseRowVariances(IntegerVector j, NumericVector val, NumericVector rm, int n) {
  const int nv = j.size();
  const int nm = rm.size();
  Rcpp::NumericVector rv(nm);
  Rcpp::NumericVector rit(nm);
  int current;
  // Calculate RowVars Initial
  for (int i = 0; i < nv; ++i) {
    current = j(i) - 1;
    rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
    rit(current) = rit(current) + 1;
  }
  // Calculate Remainder Variance
  for (int i = 0; i < nm; ++i) {
    rv(i) = rv(i) + (n - rit(i))*rm(i)*rm(i);
  }
  rv = rv / (n - 1);
  return(rv);
}

