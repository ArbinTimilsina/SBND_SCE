rm -rf MakeMapSCE
find ../../HistoDirectory/ -name "*.pdf" | xargs rm
rm -rf ../../OutputFiles/Result_*
g++ -Wall -std=c++11 MakeMapSCE.C ../../MakeMapSource/MapSCE.C `root-config --cflags --glibs` -o MakeMapSCE;
./MakeMapSCE
rm -rf ../../OutputFiles/Result_*
rm -rf MakeMapSCE