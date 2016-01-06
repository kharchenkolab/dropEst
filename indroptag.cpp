#include <iostream>
#include <iomanip>
#include <ostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>
#include <climits>
#include <list>
#include <functional>
#include <algorithm>
#include <sstream>
#include <numeric>
#include <utility>
#include <set>
#include <getopt.h>
#include <ext/hash_set>
#include <ext/hash_map>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "edit_distance.h"

using namespace std;
using namespace __gnu_cxx; 

#undef DEBUG
//#define DEBUG 1

static void usage() {
  cerr << "\tindroptag -- generate tagged indrop fastq files for alignment"<<endl;
  cerr << "SYNOPSIS\n";
  cerr << "\tindroptag [-m|--max-reads 10] [-n|--name baseName] [-f|--tag-filter] [-v|--verbose] read_1.fastq read_2.fastq"<<endl;
  cerr << "OPTIONS:\n";
  cerr << "\t-f, --tag-filter perform tag reduction/filtering step"<<endl;
  cerr << "\t-m, --max-reads n split fastq files to contain no more than n million reads"<<endl;
  cerr << "\t-n, --name BASE_NAME specify alternative output base name"<<endl;
  cerr << "\t-l, --min-length n minimum read length to output after trimming"<<endl;
  cerr << "\t-r, --rclength n specify reverse complement R1 length to lookup in R2"<<endl;
}

// reverse complement
string rc(string& s) {
  char rcs[s.length()];
  
  for(int i=0;i<s.length();i++) {
    switch(s.at(s.length()-i-1)) {
    case 'A':
      rcs[i]='T';
      break;
    case 'T':
      rcs[i]='A';
      break;
    case 'C':
      rcs[i]='G';
      break;
    case 'G':
      rcs[i]='C';
      break;
    default:
      rcs[i]='N';
      break;
    }
  }
  string rcss=string(rcs,s.length()); 
  return(rcss);
}


int main(int argc,char **argv) {
  bool verbose=false;
  bool filter=false;
  int maxreads=INT_MAX;
  string basename;
  bool basenameset=false;

  // inDrop read structure: [C8-11][S][C8][M6], S=GAGTGATTGCTTGTGACGCCTT

  string spacer="GAGTGATTGCTTGTGACGCCTT";
  string spacer2="GAGTGATTGCTTGTGCCGCCTT";
  int spacerminpos=8;
  int spacermaxpos=11;
  int firstlen=8;
  int umilen=6;


  /*
// examples of Huidan's reads
TGAGGTCT GTGATTGCTTGTGACGCCTT GTTTGTTT ATAGCG AGCC
TAACCCGT GTGATTGCTTGTGACGCCTC TCCTTATT TCATCG TCCTTTTTTTTTTTTTTTTTTTTG
GCATGGGT GTGATTGCTTGTGACGCCTT TCGACGGT CCCCTA AATATTTTTTTTTTTTTTTTTTTT
TGAGGTCT GTGATTGCTTGTGACGCCTT GTTTGTTT CCGACA TGAGTTTTTTTTTTTTTTTTTCACG
CCGCTGTT GTGATTGCTTGTGACGCCTT TGCTCCGT ACCCGA CAGTTTTTTTTTTTTTTTTTTTTTT
CCAGCAGT GTGATTGCTTGTGACGCCTT CATTTGTT CCCGCC TCATTAACCTATGGATTCAGTTAAT
CTTCAGGT GTGATTGCTTGTGACGCCTT AGATGGCT TCCAAT CCGATTTTTTTTTTTTTTTTTTTTT
  */
  // Huidan read structure: [C8][S][C8][M6], S=GTGATTGCTTGTGACGCCTT
  /*
  string spacer="GTGATTGCTTGTGACGCCTT";
  string spacer2="GTGATTGCTTGTGACGCCTC";
  int spacerminpos=7;
  int spacermaxpos=9;
  int firstlen=8;
  int umilen=6;
  */

  int rclen=8; // RC of R1 to lookup in R2 for trimming
  string polyA("AAAAAAAA");
  int minalignlen=10; // minimum length of R2 after trimming to output for the alignment
  int max_spacer_ed=3; // maximum number of spacer mismatches
  int prefix_len=5; // length of the spacer prefix/suffix to use as a seed
  
  int option_index=0;
  int c;
  static struct option long_options[] = {
    {"verbose", no_argument, 0, 'v'},
    {"name", required_argument, 0, 'n'},
    {"max-reads", required_argument, 0, 'm'},
    {"min-length", required_argument, 0, 'l'},
    {"spacer", required_argument, 0, 's'},
    {"rclength", required_argument, 0, 'r'},
    {"filter", no_argument, 0, 'f'},
    { 0, 0, 0, 0 }
  };
  while (( c=getopt_long(argc, argv, "vfn:m:s:l:r:",long_options, &option_index )) != -1) {
    switch(c) {
    case 'v' :
      verbose=true;
      break;
    case 'f' :
      filter=true;
      break;
    case 'm' :
      maxreads = atoi( optarg )*1000000;
      break;
    case 'l' :
      minalignlen = atoi( optarg );
      break;
    case 'r' :
      rclen = atoi( optarg );
      break;
    case 'n' :
      basename = string( optarg );
      basenameset=true;
      break;
    case 's' :
      spacer = string( optarg );
      break;
    default:
      cerr<<"indroptag: unknown arguments passed"<<endl;
      usage();
      return 1;
    }
  }

  int minlength=spacerminpos+firstlen+umilen+spacer.length();

  if(optind != argc - 2) {
    cerr << "indroptag: two read files must be provided"<<endl;
    usage();
    return 1;
  }

  string r1_file = string(argv[optind++]);
  string r2_file = string(argv[optind++]);
  if(!basenameset) {
    basename=r2_file+".tagged";
  }
  
  if(verbose) cout<<"reading reads "<<flush;

  
  // open input files
  //ifstream r1fh(r1_file.c_str());
  //ifstream r2fh(r2_file.c_str());
  ifstream r1f(r1_file.c_str(),std::ios_base::in | std::ios_base::binary);
  if(!r1f) {
    cerr<<"can't open R1 file \""<<r1_file<<"\""<<endl;
    return 1;
  }
  ifstream r2f(r2_file.c_str(),std::ios_base::in | std::ios_base::binary);  
  if(!r2f) {
    cerr<<"can't open R2 file \""<<r2_file<<"\""<<endl;
    return 1;
  }

  boost::iostreams::filtering_istream r1fh;
  if(boost::ends_with(r1_file,".gz") || boost::ends_with(r1_file,".gzip")) {
    r1fh.push(boost::iostreams::gzip_decompressor());
  }
  r1fh.push(r1f);
  boost::iostreams::filtering_istream r2fh;
  if(boost::ends_with(r2_file,".gz") || boost::ends_with(r2_file,".gzip")) {
    r2fh.push(boost::iostreams::gzip_decompressor());
  }
  r2fh.push(r2f);

  int ofi=1;
  stringstream ss; 
  ss<<basename<<"."<<ofi<<".fastq.gz";
  string ofname=ss.str();

  ofstream ofs(ofname.c_str(), ios_base::out | ios_base::binary);
  //boost::iostreams::filtering_streambuf<boost::iostreams::output> ofb;
  boost::iostreams::filtering_ostream ofile;
  ofile.push(boost::iostreams::gzip_compressor());
  ofile.push(ofs);
  //ostream ofile(&ofb);
  
  //ofstream ofile;
  //ofile.open(ofname.c_str());

  long int readn=0; // reads read
  long int readc=0; // reads written out in to the current output file
  
  
  int outcomes[5]; for(int i=0;i<5;i++) { outcomes[i]=0; };
  int trimst[4]; for(int i=0;i<4;i++) { trimst[i]=0; };
  
  string r1l1,r1l2,r1l3,r1l4,r2l1,r2l2,r2l3,r2l4;
  while(getline(r1fh,r1l1)) {
    readn++;
    if(!getline(r1fh,r1l2)) { cerr<<"read "<<readn<<": R1 fastq ended prematurely!"<<endl; break;}
    if(!getline(r1fh,r1l3)) { cerr<<"read "<<readn<<": R1 fastq ended prematurely!"<<endl; break;}
    if(!getline(r1fh,r1l4)) { cerr<<"read "<<readn<<": R1 fastq ended prematurely!"<<endl; break;}
    if(!getline(r2fh,r2l1)) { cerr<<"read "<<readn<<": R2 fastq ended prematurely!"<<endl; break;}
    if(!getline(r2fh,r2l2)) { cerr<<"read "<<readn<<": R2 fastq ended prematurely!"<<endl; break;}
    if(!getline(r2fh,r2l3)) { cerr<<"read "<<readn<<": R2 fastq ended prematurely!"<<endl; break;}
    if(!getline(r2fh,r2l4)) { cerr<<"read "<<readn<<": R2 fastq ended prematurely!"<<endl; break;}
    if(r1l1.at(0)!='@') { cerr<<"read "<<readn<<": R1 fastq malformed!"<<endl; break;}
    if(r2l1.at(0)!='@') { cerr<<"read "<<readn<<": R2 fastq malformed!"<<endl; break;}
    
    if(verbose && (readn % 1000000 == 0)) { 
      cout<<"."<<flush; 
    }
    
    if(readc>maxreads) {
      ofile.pop();
      ofs.close();
      ss.str(string());
      ofi++;
      ss<<basename<<"."<<ofi<<".fastq.gz";
      ofname=ss.str();
      ofs.open(ofname.c_str(), ios_base::out | ios_base::binary);
      ofile.push(ofs);
      
      //ofile.open(ofname.c_str());
      readc=0;
      if(verbose) cout<<"|"<<flush;
    }
    
    // parse out R1
#ifdef DEBUG
    cout<<r1l2<<":"<<endl;
#endif
    
    if(r1l2.length()<minlength) { 
      outcomes[2]++;
#ifdef DEBUG
      cout<<"-- read is too short"<<endl;
#endif
      continue; 
    }
    
    int spacerpos=r1l2.find(spacer);
    if(spacerpos==string::npos) {
      // try looking for the secondary spacer
      spacerpos=r1l2.find(spacer2);
      if(spacerpos==string::npos) {
#ifdef DEBUG
	cout<<"-- spacer not found"<<endl;
#endif
	// try partial spacer matches
	
	// try suffix lookup
	int prepos;
	int postpos=r1l2.rfind(spacer.substr(spacer.length()-prefix_len,prefix_len));
	if(postpos!=string::npos) { 
	  prepos=postpos-spacer.length()+prefix_len; 
	
	  if(prepos>=spacerminpos && prepos<(spacermaxpos-spacer.length()+prefix_len)) {
	  // check edit distance
	    int ed=edit_distance(spacer.c_str(),r1l2.substr(prepos,spacer.length()).c_str());
#ifdef DEBUG
	    cout<<"-- postfix match at "<<postpos<<" ed="<<ed<<endl;
	    cout<<"-- given:"<<spacer<<endl;;
	    cout<<"-- match:"<<r1l2.substr(prepos,spacer.length())<<endl;
#endif
	    if(ed<=max_spacer_ed) {
	      spacerpos=prepos;
	    }
	  } 
	}
	  
	if(spacerpos==string::npos) { 
	  prepos=r1l2.rfind(spacer.substr(0,prefix_len));
	  if(prepos>=spacerminpos && prepos<(spacermaxpos-spacer.length()+prefix_len)) {
	    // check edit distance
	    int ed=edit_distance(spacer.c_str(),r1l2.substr(prepos,spacer.length()).c_str());
#ifdef DEBUG
	    cout<<"-- prefix match at "<<prepos<<" ed="<<ed<<endl;
#endif
	    if(ed<=max_spacer_ed) {
	      spacerpos=prepos;
	    }
	  } 
	}
	if(spacerpos==string::npos) {
	  outcomes[1]++;
	  continue; 
	}
      } else {
	outcomes[4]++;
#ifdef DEBUG
	cout<<"-- secondary spacer"<<endl;
#endif
      }
    }
    
    if(spacerpos<spacerminpos || spacerpos>spacermaxpos) { 
      outcomes[3]++;
#ifdef DEBUG
      cout<<"-- invalid spacer position"<<endl;
#endif
      continue; 
    }

    if(r1l2.length() < spacerpos+spacer.length()+firstlen+umilen+1) {
      outcomes[2]++;
#ifdef DEBUG
      cout<<"-- read is too short"<<endl;
#endif
      continue; 

    }
    
    string cellbarcode=r1l2.substr(0,spacerpos) + r1l2.substr(spacerpos+spacer.length(),firstlen);
#ifdef DEBUG
    cout<<"-- cell barcode: "<<cellbarcode<<" ("<<cellbarcode.length()<<"nt)"<<endl;
#endif

    string umibarcode=r1l2.substr(spacerpos+spacer.length()+firstlen,umilen);
#ifdef DEBUG
    cout<<"-- umi barcode: "<<umibarcode<<endl;
    cout<<"R2: "<<r2l2<<endl;
#endif
    outcomes[0]++;
    
    // clean up R2
    
    int r2trim=r2l2.length();
    // attempt 1: check for reverse complement of the UMI+second barcode, remove trailing As
    // RC of UMI+second barcode (up to a length rclen)
    string rcb=r1l2.substr(spacerpos+spacer.length()+firstlen+umilen-rclen,rclen);
    rcb=rc(rcb);
    
#ifdef DEBUG
    cout<<"-- barcode RC: "<<rcb<<endl;
#endif    
    int rcpos=r2l2.find(rcb);
    if(rcpos!=-1) {
      r2trim=rcpos;
      trimst[1]++;
#ifdef DEBUG
      cout<<"-- found barcode RC at "<<rcpos<<endl;
#endif    
    } 

    else {
      // attempt 2: find polyA block
      rcpos=r2l2.find(polyA);
      if(rcpos!=-1) {
	r2trim=rcpos;
	trimst[2]++;
#ifdef DEBUG
	cout<<"-- found polyA at "<<rcpos<<endl;
#endif
      }
    } 

    // attempt 3: trim trailing As
    bool at=false;
    while(r2trim>0 && (r2l2.at(r2trim-1)=='A' || r2l2.at(r2trim-1)=='N')) {
      r2trim--;
      at=true;
#ifdef DEBUG
      cout<<"-"<<flush;
#endif
    }
    if(at) { trimst[3]++;}
#ifdef DEBUG
    cout<<"   trimming "<<(r2l2.length()-r2trim)<<endl;
#endif
    
    int tm=r2l2.length()-r2trim;
    if(tm>0) {
      r2l2=r2l2.substr(0,r2trim);
      r2l4=r2l4.substr(0,r2trim);
    } else {
      trimst[0]++;
    }
    
#ifdef DEBUG
    cout<<" trimmed:"<<r2l2<<endl;
#endif
    // output
    if(r2trim>minalignlen) {
      //ofile<<r2l1<<"\tCB:"<<cellbarcode<<"\tUB:"<<umibarcode<<endl;
      ofile<<'@'<<readn<<'!'<<cellbarcode<<'#'<<umibarcode<<endl;
      ofile<<r2l2<<endl<<r2l3<<endl<<r2l4<<endl;
      readc++;
    }
  }

  //ofile.close();
  ofile.pop();
  ofs.close();


  if(verbose) cout<<" ("<<readn<<" reads)"<<endl<<flush;
  
  if(verbose) { 
    cout<<" outcomes:[ (OK) (no spacer) (short) (spacer misplaced) (spacer2)]\n";
    cout<<" outcomes:["; for(int i=0;i<5;i++) { cout<<outcomes[i]<<" ";} cout<<"]"<<endl; 
    cout<<" outcomes:["<<setprecision(1); for(int i=0;i<5;i++) { cout<<(((double)outcomes[i])/((double)readn)*100.0)<<" ";} cout<<"] %"<<endl; 
  }


  if(verbose) { 
    cout<<" trimst:[(no trim) (RC) (polyA) (-A)]\n";
    int ttream=0;
    cout<<" trimst:["; for(int i=0;i<4;i++) { cout<<trimst[i]<<" "; ttream+=trimst[i]; } cout<<"]"<<endl; 
    cout<<" trimst:["<<setprecision(1); for(int i=0;i<4;i++) { cout<<(((double)trimst[i])/((double)ttream)*100.0)<<" ";} cout<<"] %"<<endl; 
  }
}
