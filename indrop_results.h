#ifndef INDROP_RESULTS_H
#define INDROP_RESULTS_H 1

#include <vector>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/assume_abstract.hpp>

using namespace std;
class count_matrix {
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int /* file_version */){
    ar & cell_names & gene_names & counts;
  }
public:
  vector<string> cell_names;
  vector<string> gene_names;
  vector<int> counts;
  count_matrix() {};
  count_matrix(vector<string>& cn,vector<string>& gn,vector<int>& c): cell_names(cn), gene_names(gn), counts(c) {}
};

class indrop_results {
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int /* file_version */){
    ar & cm & non_exon_counts & non_exon_count_names & reads_per_umi & umig_covered & exon_counts & exon_count_names & merge_n;
  }
public:
  count_matrix cm;
  vector<int> non_exon_counts;
  vector<string> non_exon_count_names;
  vector<int> exon_counts;
  vector<string> exon_count_names;
  vector<double> reads_per_umi;
  vector<int> umig_covered;
  vector<int> merge_n;
  indrop_results() {};
 indrop_results(count_matrix& _cm,vector<int>& nec,vector<string>& necm,vector<double> rpu,vector<int> u, vector<int>& ec, vector<string>& ecn,vector<int>& mn): cm(_cm), non_exon_counts(nec), non_exon_count_names(necm), reads_per_umi(rpu), umig_covered(u), exon_counts(ec), exon_count_names(ecn), merge_n(mn) {}
};


#endif
