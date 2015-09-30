#g++ -O3 -m64 -Wno-deprecated -I /usr/include/R -o indroptag indroptag.cpp  -lm -lboost_iostreams
g++ -O3 -m64 -std=c++11 -I /home/pkharchenko/samtools/trunk -I /usr/include/R -o indropest indropest.cpp  -lz -lm -lbam -lboost_iostreams
