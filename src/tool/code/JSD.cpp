// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash.hpp>

using namespace Rcpp;

template <class T>
inline void hash_combine(std::size_t& seed, T const& v)
{
  seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{  
  template<typename T, size_t N>
  struct hash<array<T, N> >
  {
    typedef array<T, N> argument_type;
    typedef size_t result_type;
    result_type operator()(const argument_type& a) const
    {
      hash<T> hasher;
      result_type h = 0;
      for (result_type i = 0; i < N; ++i)
	{
	  h = h * 31 + hasher(a[i]);
	}
      return h;
    }
  };
  
  // How to use std::vector? //
  template<typename T>
  struct hash<vector<T> >
  {
    typedef vector<T> argument_type;
    typedef size_t result_type;
    result_type operator()(const argument_type& a) const
    {
      size_t size = a.size();
      size_t seed = 0;
      for (size_t i = 0; i < size; i++)
	//Combine the hash of the current vector with the hashes of the previous ones
	hash_combine(seed, a[i]);
      return seed;
    }    
  };
}

// [[Rcpp::export]]
NumericVector naomitInt(NumericVector x){
  std::vector<int> r(x.size());
  int k=0;
  for(int i = 0; i < x.size(); ++i) {
    if(x[i]==x[i]) {
      r[k] = x[i];
      k++;
    }
  }
  r.resize(k);
  return Rcpp::wrap(r);    
}

//[[Rcpp::export]]
NumericVector computeJSDforEachCell(NumericVector& target_cells_expr, NumericVector& distr_vector, NumericVector& ideal_vector){
  int size_per = target_cells_expr.size();
  NumericVector target_cells_jsd(size_per);
  NumericVector distr_vector_copy(distr_vector.size());
  NumericVector mean_distr_vector_copy(distr_vector.size());
  NumericVector kl2_vector(distr_vector.size());
  for(int i=0; i < size_per; i++){
    for(int k=0; k < distr_vector_copy.size(); k++){ // recreate distr_vector_copy
      distr_vector_copy[k] = distr_vector[k];
    }      
    distr_vector_copy[0] = target_cells_expr[i];
    double distr_sum = sum(distr_vector_copy);
    if(distr_sum==0){
      target_cells_jsd[i] = 1;
      continue;
    }
    for(int k=0; k < distr_vector_copy.size(); k++){
      distr_vector_copy[k] = distr_vector_copy[k] / distr_sum;
    }    
    for(int k=0; k < distr_vector_copy.size(); k++){
      mean_distr_vector_copy[k] = 0.5 * (ideal_vector[k] + distr_vector_copy[k]);
      if(mean_distr_vector_copy[k]!=0 && distr_vector_copy[k]!=0){
	kl2_vector[k] = distr_vector_copy[k] * log2(distr_vector_copy[k] / mean_distr_vector_copy[k]);
      }else{
	kl2_vector[k] = 0;
      }
    }
    double KL1 = log2(1/mean_distr_vector_copy[0]);
    double KL2 = sum(kl2_vector);
    target_cells_jsd[i] = 0.5*(KL1+KL2);
  }
  return target_cells_jsd;
}

//[[Rcpp::export]]
NumericMatrix computeJSDforEachTF(NumericMatrix& target_subpopulation_data, NumericMatrix& JSD_backround_vectors, NumericVector& ideal_vector){
  NumericMatrix JSD_values(target_subpopulation_data.ncol(), target_subpopulation_data.nrow());
  for(int i=0; i < target_subpopulation_data.nrow(); i++){
    NumericVector temp1 = target_subpopulation_data(i,_);
    NumericVector temp2(JSD_backround_vectors.nrow()+1);
    temp2(0) = 1;
    for(int j=0; j < JSD_backround_vectors.nrow(); j++){
      temp2(j+1)= JSD_backround_vectors(j,i);
    }
    NumericVector target_cells_jsd = computeJSDforEachCell(temp1, temp2, ideal_vector);
    JSD_values(_,i) = target_cells_jsd;
  }  
  return(JSD_values);
}

