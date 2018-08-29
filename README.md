# Parameterization of Space Charge Effect (SCE) for SBND

### Space Charge Effect
The space charge effect is the build-up of slow-moving positive ions in a detector due to, for instance, ionization from cosmic rays, leading to a distortion of the electric field within the detector. This effect leads to a displacement in the reconstructed position of signal ionization electrons in LArTPC (Liquid Argon Time Projection Chamber) detectors.

## Brief instructions 
MapSCE::PerformTransformation performs the transformation. Field to transform can be Spatial or E-Field and dimension to transform can be X, Y, or Z. Produced histograms are saved in HistoDirectory and ROOT output files in OutputFiles directory.
All the parameters are set in MakeMapSCE.C and 
```
./WorkArea/MakeMap/runMakeMapSCE.sh 
```
runs the the code to produce the SCE offsets. 

CompareParameterizedWithInput/CompareParameterizedWithInput.C accesses the parameterized SCE and compares with the input distribution. To run it, do
```
./WorkArea/CompareParameterizedWithInput/runCompareParameterizedWithInput.sh
```

