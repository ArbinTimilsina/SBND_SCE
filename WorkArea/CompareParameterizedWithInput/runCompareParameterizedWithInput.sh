rm -rf CompareParameterizedWithInput
g++ -Wall -std=c++11 CompareParameterizedWithInput.C ../../ForLarSoft/SpaceCharge/SpaceChargeSBND.cxx `root-config --cflags --glibs` -o CompareParameterizedWithInput
./CompareParameterizedWithInput
rm -rf CompareParameterizedWithInput
