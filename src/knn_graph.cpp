#include <Rcpp.h>
#include <queue>

using namespace Rcpp;
using namespace std;

// Euclidean Distance
double euc (NumericVector x, NumericVector y){
  int n = y.length();
  double total = 0;
  for (int i = 0; i < n ; ++i) {
    total += pow(x(i)-y(i),2.0);
  }
  total = sqrt(total);
  return total;
}

// Covariance of two vectors
double covariance(NumericVector v1,NumericVector v2){
  double mean1 = mean(v1), mean2 = mean(v2);
  double sum = (double(v1[0]) - mean1) * (double(v2[0]) - mean2);
  for (unsigned int i=1; i < v1.size(); i++){
    sum += (double(v1[i]) - mean1) * (double(v2[i]) - mean2);
  }
  return double(sum) / double(v1.size()-1);
}

// Standard Deviation
double std_dev(NumericVector v1){
  return sqrt(covariance(v1, v1));
}

// Pearson Distance
double pearson_dist(NumericVector v1, NumericVector v2){
  if (std_dev(v1) * std_dev(v2) == 0){
    return 0;
  }
  return 1 - (covariance(v1,v2) / ( std_dev(v1) * std_dev(v2)));
}

// Thanks for Martin Morgan et al.
// http://gallery.rcpp.org/articles/top-elements-from-vectors-using-priority-queue/
template <int RTYPE>
class IndexComparator {
public:
  typedef typename Rcpp::traits::storage_type<RTYPE>::type STORAGE;
  
  IndexComparator(const Vector<RTYPE>& data_) : data(data_.begin()) {}
  
  inline bool operator()(int i, int j) const {
    return data[i] > data[j] || (data[i] == data[j] && j > i);
  }
  
private:
  const STORAGE* data;
};

template <>
class IndexComparator<STRSXP> {
public:
  IndexComparator( const CharacterVector& data_ ) : data(data_.begin()) {}
  
  inline bool operator()(int i, int j) const {
    return (String)data[i] > (String)data[j] || (data[i] == data[j] && j > i );
  }
  
private:
  const Vector<STRSXP>::const_iterator data;
};

template <int RTYPE>
class IndexQueue {
public:
  typedef std::priority_queue<int, std::vector<int>, IndexComparator<RTYPE> > Queue;
  
  IndexQueue(const Vector<RTYPE>& data_) : comparator(data_), q(comparator), data(data_) {}
  
  inline operator IntegerVector() {
    int n = q.size();
    IntegerVector res(n);
    for( int i=0; i<n; i++) {
      // +1 for 1-based R indexing
      res[i] = q.top() + 1;
      q.pop();
    }
    return res;
  }
  inline void input( int i) {
    // if( data[ q.top() ] < data[i] ) {
    if (comparator(i, q.top())) {
      q.pop();
      q.push(i);
    }
  }
  inline void pop() { q.pop(); }
  inline void push(int i) { q.push(i); }
  
private:
  IndexComparator<RTYPE> comparator;
  Queue q;
  const Vector<RTYPE>& data;
};


template <int RTYPE>
IntegerVector top_index(Vector<RTYPE> v, int n) {
  int size = v.size();
  
  // not interesting case. Less data than n
  if( size < n){
    return seq( 0, n-1 );
  }
  
  IndexQueue<RTYPE> q( v );
  for( int i=0; i<n; i++) q.push(i);
  for( int i=n; i<size; i++) q.input(i);
  return q;
}

IntegerVector top_index(SEXP x, int n) {
  switch (TYPEOF(x)) {
  case INTSXP: return top_index<INTSXP>(x, n);
  case REALSXP: return top_index<REALSXP>(x, n);
  case STRSXP: return top_index<STRSXP>(x, n);
  default: stop("type not handled");
  }
  return IntegerVector() ; // not used
}

//' Calculate k Nearest Neighbors from Pearson distance metric
//'
//' Each distance metric has its own function for speed/efficiency
//' This takes a sample X feature matrix and returns
//' a matrix of k nearest neighbors. This is the one for Pearson.
//'
//' @param x An m x n numeric matrix
//' @param k The number of nearest neighbors to return
//' per sample
//' @return An m x k matrix of indicies 1...k of the
//' nearest neighbors for the specified row based on
//' Pearson distance.
//' @export
// [[Rcpp::export]]
IntegerMatrix calcKNNgraph_pearson (NumericMatrix x, int k = 1){
  int outrows = x.nrow();
  int outcols = k;
  IntegerMatrix out(outrows,outcols);
  
  for (int i = 0 ; i < outrows; i++){
    NumericVector xx(outrows);
    for (int j = 0 ; j < outrows; j ++) {
      NumericVector v1 = x.row(i);
      NumericVector v2 = x.row(j);
      double d = pearson_dist(v1, v2);
      xx(j) = (-1) * d; // To get biggest values
    }
    IntegerVector ti = top_index(xx, k);
    for (int p = 0 ; p < k; p ++){
      out(i,p) = ti(k-p-1);
    }
  }
  return out;
}


//' Calculate k Nearest Neighbors from Euclidean distance metric
//'
//' Each distance metric has its own function for speed/efficiency
//' This takes a sample X feature matrix and returns
//' a matrix of k nearest neighbors. This is the one for Euclidean.
//'
//' @param x An m x n numeric matrix
//' @param k The number of nearest neighbors to return
//' per sample
//' @return An m x k matrix of indicies 1...k of the
//' nearest neighbors for the specified row based on
//' Euclidean distance.
//' @export
// [[Rcpp::export]]
IntegerMatrix calcKNNgraph_euclidean (NumericMatrix x, int k = 1){
  int outrows = x.nrow();
  int outcols = k;
  IntegerMatrix out(outrows,outcols);
  
  for (int i = 0 ; i < outrows; i++){
    NumericVector xx(outrows);
    for (int j = 0 ; j < outrows; j ++) {
      NumericVector v1 = x.row(i);
      NumericVector v2 = x.row(j);
      double d = euc(v1, v2);
      xx(j) = (-1) * d; // To get biggest values
    }
    IntegerVector ti = top_index(xx, k);
    for (int p = 0 ; p < k; p ++){
      out(i,p) = ti(k-p-1);
    }
  }
  return out;
}