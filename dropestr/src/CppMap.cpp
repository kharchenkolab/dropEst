#include "dropestr.h"
#include "CppMap.h"

CppMap::CppMap()
{}

CppMap::CppMap(s_vec_t names, std::vector<double> values)
{
  if (names.size() != values.size())
    Rcpp::stop("Names and values must have the same size");

  for (int i = 0; i < values.size(); ++i) {
    this->_data.emplace(std::string(names[i]), values[i]);
  }
}

CppMap::CppMap(SEXP sexp){
  Rcpp::XPtr<CppMap> xp = UnwrapRobject<CppMap>(sexp);
  this->_data = xp->_data;
}

void CppMap::set(const std::string &key, double value) {
  try {
    this->_data[key] = value;
  }
  catch(std::exception err) {
    Rcpp::stop("Error: " + std::string(err.what()));
  }
}

Rcpp::NumericVector CppMap::at(const s_vec_t &keys) const {
  Rcpp::NumericVector res(keys.size());
  for (size_t i = 0; i < keys.size(); ++i) {
    auto iter = this->_data.find(keys[i]);
    if (iter == this->_data.end())
      Rcpp::stop("Can't find key: " + keys[i]);

    res[i] = iter->second;
  }

  res.attr("names") = Rcpp::wrap(keys);
  return res;
}

size_t CppMap::size() const {
  return this->_data.size();
}

void CppMap::clear() {
  this->_data.clear();
}

const CppMap::data_t& CppMap::data() const {
  return this->_data;
}

RCPP_MODULE(CppMapModule) {
  using namespace Rcpp;

  class_<CppMap>("CppMap")
    .default_constructor("Default constructor")
    .constructor<SEXP>()
    .constructor<s_vec_t, std::vector<double>>()
    .method("at", &CppMap::at)
    .method("set", &CppMap::set)
    .method("size", &CppMap::size)
    .method("clear", &CppMap::clear)
  ;
}
