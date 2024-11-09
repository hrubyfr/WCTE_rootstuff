#include "Math/Vector2D.h"
#include "Math/Point2D.h"
#include <TMath.h>
#include <iostream>
#include <vector>
#include "TH2D.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TMarker.h"

std::vector<ROOT::Math::Polar2DVector> Generate_SiPMTVector(double radius){
    std::vector<ROOT::Math::Polar2DVector> SiPMTVector;
    for (int i = 0; i < 16; i++){
        int pm_number = i;
        double angle = i * TMath::Pi()/8;
        ROOT::Math::Polar2DVector helpVector(radius, angle);
        SiPMTVector.push_back(helpVector);
        // std::cout << "______________SiPMT " << pmt_number << "______________" << endl;
        // std::cout << "SiPMTVector check: " << SiPMTVector[pmt_number] << endl;
        // std::cout << "SiPMTVector coordinates:\t x:" << SiPMTVector[pmt_number].X() << "\t y:" << SiPMTVector[pmt_number].Y() << endl;
    }
    return SiPMTVector;
}

void fillBeamPosition(TH2D* histogram, std::vector<double> SiPM_charge, std::vector<ROOT::Math::Polar2DVector> SiPMVector){
    ROOT::Math::XYPoint BeamPosition;
    double Total_charge = 0;
    for(int i = 0; i < SiPMVector.size(); i++){
        BeamPosition += SiPM_charge[i] * SiPMVector[i];
        Total_charge += SiPM_charge[i];
    }
    BeamPosition /= Total_charge;
    histogram->Fill(BeamPosition.X(), BeamPosition.Y());
}

void DrawBeamHistogram(std::vector<ROOT::Math::Polar2DVector> SiPMVector, std::vector<std::vector<double>> SiPM_charges, double Scint_radius){
    TString can_title = "Beam histogram no attenuation";
    auto Beam_canvas = new TCanvas(can_title, can_title, 800, 800);

    TString Hist_title = "Histogram of approximated beam position from SiPM integrated charge";
    TString Hist_xAxis = "x [cm]";
    TString Hist_yAxis = "y [cm]";

    auto Beam_histogram = new TH2D("", "", 100, -Scint_radius, Scint_radius, 100, -Scint_radius, Scint_radius);
    Beam_histogram->SetTitle(Hist_title);
    Beam_histogram->GetXaxis()->SetTitle(Hist_xAxis);
    Beam_histogram->GetYaxis()->SetTitle(Hist_yAxis);

    

    for (int i = 0; i < SiPM_charges.size(); i++){
        fillBeamPosition(Beam_histogram, SiPM_charges[i], SiPMVector);
    }
    Beam_histogram->Draw("colz");

    for(int i = 0; i < SiPMVector.size(); i++){
        auto SiPM_position = new TMarker(SiPMVector[i].X(), SiPMVector[i].Y(), kFullCircle);
        SiPM_position->SetMarkerColor(kBlue);
        SiPM_position->Draw("same");
    }

    auto ellipse = new TEllipse(0, 0, Scint_radius, Scint_radius);
    ellipse->SetFillStyle(0);
    ellipse->SetLineColor(kBlack);
    ellipse->Draw("same");
}

void BeamPositionMonitor() {

    //### CONSTANTS ###

    double radius = 3.0;
    int point_number = 1000;

    //### GRID PART ###

    const char* filename = "BeamMonitorGrid.pdf";

    std::vector<std::vector<double>> PMT_signal_NoAtt;

    auto SiPMTVector = Generate_SiPMTVector(radius);
    std::vector<std::vector<double>> SiPM_charges;
    for (int i = 0; i < 10; i++){
        std::vector<double> SiPM_charge = {3, 2.5, 2, 1.5, 1, 0, 1, 1.5, 3, 2.5, 2, 1.5, 1, 0, 1, 1.5};
        SiPM_charges.push_back(SiPM_charge);
    }

    //### BEAM POSITION HISTOGRAM FROM INTEGRATED CHARGE ###
    DrawBeamHistogram(SiPMTVector, SiPM_charges, radius);
}