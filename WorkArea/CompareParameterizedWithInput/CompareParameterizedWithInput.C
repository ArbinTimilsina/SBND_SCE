// C++ stuff
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

// Root stuff
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TRatioPlot.h>

#include "../../ForLarSoft/SpaceCharge/SpaceChargeSBND.h"

using namespace std;

void setHisto(TH1 *histo, TString xAxis, TString yAxis)
{
    histo->GetXaxis()->SetTitle(xAxis);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetXaxis()->SetTitleOffset(0.75);
    histo->GetXaxis()->SetTitleFont(62);
    histo->GetXaxis()->SetTickLength(0.04);
    histo->GetXaxis()->SetLabelFont(62);
    histo->GetXaxis()->SetLabelSize(0.045);

    histo->GetYaxis()->SetTitle(yAxis);
    histo->GetYaxis()->SetTitleSize(0.06);
    histo->GetYaxis()->SetTitleOffset(0.85);
    histo->GetYaxis()->SetTitleFont(62);
    histo->GetYaxis()->SetTickLength(0.04);
    histo->GetYaxis()->SetLabelFont(62);
    histo->GetYaxis()->SetLabelSize(0.045);

    histo->SetMinimum(1.0);
    histo->SetTitle("");
}

int main()
{
    bool setLog = true;

    TFile *fileInput = new TFile("../../InputFiles/dispOutput_SBND_E500.root");

    int initialSpatialN[3] = {3, 4, 3};
    int intermediateSpatialN[3] = {4, 4, 4};

    int initialEFieldN[3] = {3, 3, 3};
    int intermediateEFieldN[3] = {6, 4, 4};

    SpaceChargeSBND *mapSBND = new SpaceChargeSBND("../../OutputFiles/SCEoffsets_SBND_E500.root",
						   initialSpatialN, intermediateSpatialN, initialEFieldN, intermediateEFieldN,
						   500.0);
    int xMin = -206;
    int xMax = 206;
    int yMin = -210;
    int yMax = 210;
    int zMin = -10;
    int zMax = 510;

    const int dBins = 150;
    const double minD = -6.0;
    const double maxD = 6.0;
    const int eBins = 100;
    const double minE = -0.04;
    const double maxE = 0.04;

    TH1D *iDx = new TH1D("iDx", "", dBins, -0.01, 0.006);
    TH1D *iDy = new TH1D("iDy", "", dBins, minD / 100.0, maxD / 100.0);
    TH1D *iDz = new TH1D("iDz", "", dBins, minD / 100.0, maxD / 100.0);
    TH1D *iEx = new TH1D("iEx", "", eBins, -3000.0, 3000.0);
    TH1D *iEy = new TH1D("iEy", "", eBins, minE * 100 * 500, maxE * 100 * 500);
    TH1D *iEz = new TH1D("iEz", "", eBins, minE * 100 * 500, maxE * 100 * 500);

    TH1D *mDx = new TH1D("mDx", "", dBins, -1.0, 0.6);
    TH1D *mDy = new TH1D("mDy", "", dBins, minD, maxD);
    TH1D *mDz = new TH1D("mDz", "", dBins, minD, maxD);
    TH1D *mEx = new TH1D("mEx", "", eBins, -0.06, 0.06);
    TH1D *mEy = new TH1D("mEy", "", eBins, minE, maxE);
    TH1D *mEz = new TH1D("mEz", "", eBins, minE, maxE);

    TTree *TreeD = (TTree*)fileInput->Get("SpaCEtree_fwdDisp");
    TreeD->Draw("Dx>>iDx");
    TreeD->Draw("Dy>>iDy");
    TreeD->Draw("Dz>>iDz");
    TTree *TreeE = (TTree*)fileInput->Get("SpaCEtree");
    TreeE->Draw("Ex>>iEx");
    TreeE->Draw("Ey>>iEy");
    TreeE->Draw("Ez>>iEz");

    int nSkip = 5;
    for(int iX = xMin; iX <= xMax - nSkip; iX++)
        {
            iX = iX + nSkip;
            for(int iY = yMin; iY <= yMax - nSkip; iY++)
                {
                    iY = iY + nSkip;
                    for(int iZ = zMin; iZ <= zMax - nSkip; iZ++)
                        {
                            iZ = iZ + nSkip;
                            cout << iX << ", " << iY << ", " << iZ << endl;

                            vector<double> mSpatialOffsets = mapSBND->GetPosOffsets(iX, iY, iZ);
                            if(!((mSpatialOffsets.at(0) == mSpatialOffsets.at(1)) &&
				 (mSpatialOffsets.at(1) == mSpatialOffsets.at(2)) &&
				 (mSpatialOffsets.at(2) == 0.0)))
                                {
                                    mDx->Fill(mSpatialOffsets.at(0));
                                    mDy->Fill(mSpatialOffsets.at(1));
                                    mDz->Fill(mSpatialOffsets.at(2));
                                }

                            vector<double> mEfieldOffsets = mapSBND->GetEfieldOffsets(iX, iY, iZ);
                            if(!((mEfieldOffsets.at(0) == mEfieldOffsets.at(1)) &&
				 (mEfieldOffsets.at(1) == mEfieldOffsets.at(2)) &&
				 (mEfieldOffsets.at(2) == 0.0)))
                                {
                                    mEx->Fill(mEfieldOffsets.at(0));
                                    mEy->Fill(mEfieldOffsets.at(1));
                                    mEz->Fill(mEfieldOffsets.at(2));
                                }
                        }
                }
        }


    TCanvas *cComparisionD = new TCanvas("cComparisionD", "", 1300, 1500);
    cComparisionD->Divide(2, 3);

    setHisto(iDx, "dDx (m)", "Count");
    iDx->SetTitle("Input");
    setHisto(mDx, "dDx (cm)", "");
    mDx->SetTitle("Parameterized");
    setHisto(iDy, "dDy (m)", "Count");
    setHisto(mDy, "dDy (cm)", "");
    setHisto(iDz, "dDz (m)", "Count");
    setHisto(mDz, "dDz (cm)", "");

    cComparisionD->cd(1);
    if(setLog)
        {
            gPad->SetLogy();
        }
    iDx->Draw();
    cComparisionD->cd(2);
    if(setLog)
        {
            gPad->SetLogy();
        }
    mDx->Draw();
    cComparisionD->cd(3);
    if(setLog)
        {
            gPad->SetLogy();
        }
    iDy->Draw();
    cComparisionD->cd(4);
    if(setLog)
        {
            gPad->SetLogy();
        }
    mDy->Draw();
    cComparisionD->cd(5);
    if(setLog)
        {
            gPad->SetLogy();
        }
    iDz->Draw();
    cComparisionD->cd(6);
    if(setLog)
        {
            gPad->SetLogy();
        }
    mDz->Draw();

    TCanvas *cComparisionE = new TCanvas("cComparisionE", "", 1300, 1500);
    cComparisionE->Divide(2, 3);

    setHisto(iEx, "dEx (V/m)", "Count");
    iEx->SetTitle("Input");
    setHisto(mEx, "-dEx/E", "");
    mEx->SetTitle("Parameterized");
    setHisto(iEy, "dEy (V/m)", "Count");
    setHisto(mEy, "-dEy/E", "");
    setHisto(iEz, "dEz (V/m)", "Count");
    setHisto(mEz, "-dEz/E", "");


    cComparisionE->cd(1);
    if(setLog)
        {
            gPad->SetLogy();
        }
    iEx->Draw();
    cComparisionE->cd(2);
    if(setLog)
        {
            gPad->SetLogy();
        }
    mEx->Draw();
    cComparisionE->cd(3);
    if(setLog)
        {
            gPad->SetLogy();
        }
    iEy->Draw();
    cComparisionE->cd(4);
    if(setLog)
        {
            gPad->SetLogy();
        }
    mEy->Draw();
    cComparisionE->cd(5);
    if(setLog)
        {
            gPad->SetLogy();
        }
    iEz->Draw();
    cComparisionE->cd(6);
    if(setLog)
        {
            gPad->SetLogy();
        }
    mEz->Draw();

    cComparisionD->SaveAs("SBND_Spatial.pdf");
    cComparisionE->SaveAs("SBND_EField.pdf");
}
