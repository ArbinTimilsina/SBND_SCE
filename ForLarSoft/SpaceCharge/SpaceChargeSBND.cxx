////////////////////////////////////////////////////////////////////////////////////
// SpaceChargeSBND.cxx                                                            //
// Brief implementation of class for storing/accessing space charge distortions  //
////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <TFile.h>
#include <TCanvas.h>

#include "SpaceChargeSBND.h"

SpaceChargeSBND::SpaceChargeSBND(std::string filename,
                                 const int initialSpatialPolN[3], const int intermediateSpatialPolN[3], const int initialEFieldPolN[3], const int intermediateEFieldPolN[3],
                                 const double drift)
    : DriftField(drift)
{
    for(int i = 0; i < 3; i++)
        {
            initialSpatialFitPolN[i] = initialSpatialPolN[i];
            intermediateSpatialFitPolN[i] = intermediateSpatialPolN[i];

            initialEFieldFitPolN[i] = initialEFieldPolN[i];
            intermediateEFieldFitPolN[i] = intermediateEFieldPolN[i];
        }
    Configure(filename);
}

SpaceChargeSBND::~SpaceChargeSBND()
{
}

bool SpaceChargeSBND::Configure(std::string filename)
{

    TFile inputfile(filename.c_str(), "READ");

    for(int i = 0; i < initialSpatialFitPolN[0] + 1; i++)
        {
            for(int j = 0; j < intermediateSpatialFitPolN[0] + 1; j++)
                {
                    gSpatialGraphX[i][j] = (TGraph*)inputfile.Get(Form("deltaX/g%i_%i", i, j));

                }
            intermediateSpatialFitFunctionX[i] = new TF1(Form("intermediateSpatialFitFunctionX_%i", i), Form("pol%i", intermediateSpatialFitPolN[0]));
        }
    for(int i = 0; i < initialSpatialFitPolN[1] + 1; i++)
        {
            for(int j = 0; j < intermediateSpatialFitPolN[1] + 1; j++)
                {
                    gSpatialGraphY[i][j] = (TGraph*)inputfile.Get(Form("deltaY/g%i_%i", i, j));

                }
            intermediateSpatialFitFunctionY[i] = new TF1(Form("intermediateSpatialFitFunctionY_%i", i), Form("pol%i", intermediateSpatialFitPolN[1]));
        }
    for(int i = 0; i < initialSpatialFitPolN[2] + 1; i++)
        {
            for(int j = 0; j < intermediateSpatialFitPolN[2] + 1; j++)
                {
                    gSpatialGraphZ[i][j] = (TGraph*)inputfile.Get(Form("deltaZ/g%i_%i", i, j));

                }
            intermediateSpatialFitFunctionZ[i] = new TF1(Form("intermediateSpatialFitFunctionZ_%i", i), Form("pol%i", intermediateSpatialFitPolN[2]));
        }

    initialSpatialFitFunctionX =  new TF1("initialSpatialFitFunctionX", Form("pol%i", initialSpatialFitPolN[0]));
    initialSpatialFitFunctionY =  new TF1("initialSpatialFitFunctionY", Form("pol%i", initialSpatialFitPolN[1]));
    initialSpatialFitFunctionZ =  new TF1("initialSpatialFitFunctionZ", Form("pol%i", initialSpatialFitPolN[2]));

    for(int i = 0; i < initialEFieldFitPolN[0] + 1; i++)
        {
            for(int j = 0; j < intermediateEFieldFitPolN[0] + 1; j++)
                {
                    gEFieldGraphX[i][j] = (TGraph*)inputfile.Get(Form("deltaEx/g%i_%i", i, j));

                }
            intermediateEFieldFitFunctionX[i] = new TF1(Form("intermediateEFieldFitFunctionX_%i", i), Form("pol%i", intermediateEFieldFitPolN[0]));
        }
    for(int i = 0; i < initialEFieldFitPolN[1] + 1; i++)
        {
            for(int j = 0; j < intermediateEFieldFitPolN[1] + 1; j++)
                {
                    gEFieldGraphY[i][j] = (TGraph*)inputfile.Get(Form("deltaEy/g%i_%i", i, j));

                }
            intermediateEFieldFitFunctionY[i] = new TF1(Form("intermediateEFieldFitFunctionY_%i", i), Form("pol%i", intermediateEFieldFitPolN[1]));
        }
    for(int i = 0; i < initialEFieldFitPolN[2] + 1; i++)
        {
            for(int j = 0; j < intermediateEFieldFitPolN[2] + 1; j++)
                {
                    gEFieldGraphZ[i][j] = (TGraph*)inputfile.Get(Form("deltaEz/g%i_%i", i, j));

                }
            intermediateEFieldFitFunctionZ[i] = new TF1(Form("intermediateEFieldFitFunctionZ_%i", i), Form("pol%i", intermediateEFieldFitPolN[2]));
        }

    initialEFieldFitFunctionX =  new TF1("initialEFieldFitFunctionX", Form("pol%i", initialEFieldFitPolN[0]));
    initialEFieldFitFunctionY =  new TF1("initialEFieldFitFunctionY", Form("pol%i", initialEFieldFitPolN[1]));
    initialEFieldFitFunctionZ =  new TF1("initialEFieldFitFunctionZ", Form("pol%i", initialEFieldFitPolN[2]));

    inputfile.Close();

    return true;
}

// Primary working method of service that provides position offsets to be
// used in ionization electron drift
std::vector<double> SpaceChargeSBND::GetPosOffsets(double xVal, double yVal, double zVal) const
{
    std::vector<double> thePosOffsets;

    if(IsInsideBoundaries(xVal, yVal, zVal) == false)
        {
            thePosOffsets.resize(3, 0.0);
        }
    else
        {
            thePosOffsets = GetPosOffsetsParametric(xVal, yVal, zVal);
        }

    // GetPosOffsetsParametric returns m
    // The PosOffsets are returned as cm
    thePosOffsets.at(0) = 100.0 * thePosOffsets.at(0);
    thePosOffsets.at(1) = 100.0 * thePosOffsets.at(1);
    thePosOffsets.at(2) = 100.0 * thePosOffsets.at(2);

    return thePosOffsets;
}

// Provides position offsets using a parametric representation
std::vector<double> SpaceChargeSBND::GetPosOffsetsParametric(double xVal, double yVal, double zVal) const
{
    std::vector<double> thePosOffsetsParametric;

    double xValNew = TransformX(xVal);
    double yValNew = TransformY(yVal);
    double zValNew = TransformZ(zVal);

    thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew, yValNew, zValNew, "X"));
    thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew, yValNew, zValNew, "Y"));
    thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew, yValNew, zValNew, "Z"));

    return thePosOffsetsParametric;
}


// Provides one position offset using a parametric representation, for a given axis
double SpaceChargeSBND::GetOnePosOffsetParametric(double xValNew, double yValNew, double zValNew, std::string axis) const
{
    double parA[99][99];
    double parB[99];

    for(int i = 0; i < 99; i++)
        {
            for(int j = 0; j < 99; j++)
                {
                    parA[i][j] = 0.0;
                }
            parB[i] = 0.0;
        }
    if(axis == "X")
        {
            for(int i = 0; i < initialSpatialFitPolN[0] + 1; i++)
                {
                    for(int j = 0; j < intermediateSpatialFitPolN[0] + 1; j++)
                        {
                            parA[i][j] = gSpatialGraphX[i][j]->Eval(zValNew);
                        }
                    intermediateSpatialFitFunctionX[i]->SetParameters(parA[i]);
                }
        }
    else if(axis == "Y")
        {
            for(int i = 0; i < initialSpatialFitPolN[1] + 1; i++)
                {
                    for(int j = 0; j < intermediateSpatialFitPolN[1] + 1; j++)
                        {
                            parA[i][j] = gSpatialGraphY[i][j]->Eval(zValNew);
                        }
                    intermediateSpatialFitFunctionY[i]->SetParameters(parA[i]);
                }
        }
    else if(axis == "Z")
        {
            for(int i = 0; i < initialSpatialFitPolN[2] + 1; i++)
                {
                    for(int j = 0; j < intermediateSpatialFitPolN[2] + 1; j++)
                        {
                            parA[i][j] = gSpatialGraphZ[i][j]->Eval(zValNew);
                        }
                    intermediateSpatialFitFunctionZ[i]->SetParameters(parA[i]);
                }
        }

    double aValNew;
    double bValNew;

    if(axis == "Y")
        {
            aValNew = xValNew;
            bValNew = yValNew;
        }
    else
        {
            aValNew = yValNew;
            bValNew = xValNew;
        }
    double offsetValNew = 0.0;
    if(axis == "X")
        {
            for(int i = 0; i < initialSpatialFitPolN[0] + 1; i++)
                {
                    parB[i] = intermediateSpatialFitFunctionX[i]->Eval(aValNew);
                }
            initialSpatialFitFunctionX->SetParameters(parB);
            offsetValNew = initialSpatialFitFunctionX->Eval(bValNew);
        }
    else if(axis == "Y")
        {
            for(int i = 0; i < initialSpatialFitPolN[1] + 1; i++)
                {
                    parB[i] = intermediateSpatialFitFunctionY[i]->Eval(aValNew);
                }
            initialSpatialFitFunctionY->SetParameters(parB);
            offsetValNew = initialSpatialFitFunctionY->Eval(bValNew);
        }
    else if(axis == "Z")
        {
            for(int i = 0; i < initialSpatialFitPolN[2] + 1; i++)
                {
                    parB[i] = intermediateSpatialFitFunctionZ[i]->Eval(aValNew);
                }
            initialSpatialFitFunctionZ->SetParameters(parB);
            offsetValNew = initialSpatialFitFunctionZ->Eval(bValNew);
        }

    return offsetValNew;
}

// Primary working method of service that provides E field offsets to be
// used in charge/light yield calculation (e.g.)
std::vector<double> SpaceChargeSBND::GetEfieldOffsets(double xVal, double yVal, double zVal) const
{
    std::vector<double> theEfieldOffsets;

    if(IsInsideBoundaries(xVal, yVal, zVal) == false)
        {
            theEfieldOffsets.resize(3, 0.0);
        }
    else
        {
            theEfieldOffsets = GetEfieldOffsetsParametric(xVal, yVal, zVal);
        }

    // GetOneEfieldOffsetParametric returns V/m
    // The E-field offsets are returned as -dEx/|E_nominal|, -dEy/|E_nominal|, and -dEz/|E_nominal| where |E_nominal| is DriftField
    theEfieldOffsets.at(0) = -1.0 * theEfieldOffsets.at(0) / (100.0 * DriftField);
    theEfieldOffsets.at(1) = -1.0 * theEfieldOffsets.at(1) / (100.0 * DriftField);
    theEfieldOffsets.at(2) = -1.0 * theEfieldOffsets.at(2) / (100.0 * DriftField);

    return theEfieldOffsets;
}

// Provides E field offsets using a parametric representation
std::vector<double> SpaceChargeSBND::GetEfieldOffsetsParametric(double xVal, double yVal, double zVal) const
{
    std::vector<double> theEfieldOffsetsParametric;

    double xValNew = TransformX(xVal);
    double yValNew = TransformY(yVal);
    double zValNew = TransformZ(zVal);

    theEfieldOffsetsParametric.push_back(GetOneEfieldOffsetParametric(xValNew, yValNew, zValNew, "X"));
    theEfieldOffsetsParametric.push_back(GetOneEfieldOffsetParametric(xValNew, yValNew, zValNew, "Y"));
    theEfieldOffsetsParametric.push_back(GetOneEfieldOffsetParametric(xValNew, yValNew, zValNew, "Z"));

    return theEfieldOffsetsParametric;
}

// Provides one E field offset using a parametric representation, for a given axis
double SpaceChargeSBND::GetOneEfieldOffsetParametric(double xValNew, double yValNew, double zValNew, std::string axis) const
{
    double parA[99][99];
    double parB[99];

    for(int i = 0; i < 99; i++)
        {
            for(int j = 0; j < 99; j++)
                {
                    parA[i][j] = 0.0;
                }
            parB[i] = 0.0;
        }
    if(axis == "X")
        {
            for(int i = 0; i < initialEFieldFitPolN[0] + 1; i++)
                {
                    for(int j = 0; j < intermediateEFieldFitPolN[0] + 1; j++)
                        {
                            parA[i][j] = gEFieldGraphX[i][j]->Eval(zValNew);
                        }
                    intermediateEFieldFitFunctionX[i]->SetParameters(parA[i]);
                }
        }
    else if(axis == "Y")
        {
            for(int i = 0; i < initialEFieldFitPolN[1] + 1; i++)
                {
                    for(int j = 0; j < intermediateEFieldFitPolN[1] + 1; j++)
                        {
                            parA[i][j] = gEFieldGraphY[i][j]->Eval(zValNew);
                        }
                    intermediateEFieldFitFunctionY[i]->SetParameters(parA[i]);
                }
        }
    else if(axis == "Z")
        {
            for(int i = 0; i < initialEFieldFitPolN[2] + 1; i++)
                {
                    for(int j = 0; j < intermediateEFieldFitPolN[2] + 1; j++)
                        {
                            parA[i][j] = gEFieldGraphZ[i][j]->Eval(zValNew);
                        }
                    intermediateEFieldFitFunctionZ[i]->SetParameters(parA[i]);
                }
        }

    double aValNew;
    double bValNew;

    if(axis == "Y")
        {
            aValNew = xValNew;
            bValNew = yValNew;
        }
    else
        {
            aValNew = yValNew;
            bValNew = xValNew;
        }

    double offsetValNew = 0.0;
    if(axis == "X")
        {
            for(int i = 0; i < initialEFieldFitPolN[0] + 1; i++)
                {
                    parB[i] = intermediateEFieldFitFunctionX[i]->Eval(aValNew);
                }
            initialEFieldFitFunctionX->SetParameters(parB);
            offsetValNew = initialEFieldFitFunctionX->Eval(bValNew);
        }
    else if(axis == "Y")
        {
            for(int i = 0; i < initialEFieldFitPolN[1] + 1; i++)
                {
                    parB[i] = intermediateEFieldFitFunctionY[i]->Eval(aValNew);
                }
            initialEFieldFitFunctionY->SetParameters(parB);
            offsetValNew = initialEFieldFitFunctionY->Eval(bValNew);
        }
    else if(axis == "Z")
        {
            for(int i = 0; i < initialEFieldFitPolN[2] + 1; i++)
                {
                    parB[i] = intermediateEFieldFitFunctionZ[i]->Eval(aValNew);
                }
            initialEFieldFitFunctionZ->SetParameters(parB);
            offsetValNew = initialEFieldFitFunctionZ->Eval(bValNew);
        }

    return offsetValNew;
}

// Transform LarSoft-X (cm) to SCE-X (m) coordinate
// [-196.5, 196.5] to [0, 2.0]
double SpaceChargeSBND::TransformX(double xVal) const
{
    xVal = xVal / 100.0;

    // We use the same map twice; 0 is cathod
    // Map [0, 196.5] to [0, 2.0]
    double xValNew = (2.0 / 1.965) * fabs(xVal);

    return xValNew;
}

// Transform LarSoft-Y (cm) to SCE-Y (m) coordinate
// [-200.0, 200.0] to [0, 4.0]
double SpaceChargeSBND::TransformY(double yVal) const
{
    yVal = yVal / 100.0;

    return (yVal + 2.0);
}

// Transform LarSoft-Z (cm) to SCE-Z (m) coordinate
// [0, 500.0] to [0, 5.0]
double SpaceChargeSBND::TransformZ(double zVal) const
{
    return (zVal / 100.0);
}

// Check to see if point is inside boundaries of map
// x = [-196.5, 196.5], y = [-200, 200], and z = [0, 500]  is the active volume boundary
// For now, don't allow to go (slightly) out of range; will change if necessary
bool SpaceChargeSBND::IsInsideBoundaries(double xVal, double yVal, double zVal) const
{
    bool isInside = true;

    if((xVal < -196.5) || (xVal > 196.5) || (yVal < -200.0) || (yVal > 200.0) || (zVal < 0) || (zVal > 500.0))
        {
            isInside = false;
        }

    return isInside;
}
