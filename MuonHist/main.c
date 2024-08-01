#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <TH2.h>
#include <TCanvas.h>
#include "TMath.h"
#include <THStack.h>
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TF1.h"


struct DataPoint {
    double x;
    double y;
    double y_error;
};


void CreateHistogramFromData(const char* filename) {
    // Open the input file
    std::ifstream input(filename);

    // Check if the file was successfully opened
    if (!input.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    // Declare variables for storing data from each line of the file
    std::string line;
    std::vector<double> binCos, binPhi, flux, Emean;

    // Read lines from the file and process them
    while (getline(input, line)) {
        std::istringstream iss(line);
        double cos, phi, flx, mean;
        std::string token;

    // Tokenize the line
        if (!(iss >> cos >> phi)) { 
            // If tokenization of first and second columns fails, skip this line
            continue; 
    }

    // Read and ignore next six columns
        for (int i = 0; i < 6; ++i) {
            if (!(iss >> token)) {
                // If tokenization fails, skip this line
                continue;
        }
    }

    // Read flux and mean from ninth and tenth columns
        if (!(iss >> flx >> mean)) { 
            // If tokenization of ninth and tenth columns fails, skip this line
            continue; 
    }
        // Store tokenized values
        binCos.push_back(cos);
        binPhi.push_back(phi);
        flux.push_back(flx);
        Emean.push_back(mean);
    }

    // Create histograms from the tokenized data
    TH2D *hFluxCosmics = new TH2D("hFluxCosmics", "HK Flux", 180, 0, 360, 100, 0, 1);
    TH2D *hEmeanCosmics = new TH2D("hEmeanCosmics", "HK Emean", 180, 0, 360, 100, 0, 1);

    hFluxCosmics->GetXaxis()->SetTitle("#phi (deg)");
    hFluxCosmics->GetYaxis()->SetTitle("cos #theta");

    hEmeanCosmics->GetXaxis()->SetTitle("#phi (deg)");
    hEmeanCosmics->GetYaxis()->SetTitle("cos #theta");


    for (size_t i = 0; i < binCos.size(); ++i) {
        //std::cout << binPhi[i] << " " << binCos[i] << " " << flux[i] << " " << Emean[i] << endl;
        hFluxCosmics->SetBinContent(binPhi[i]+1, binCos[i]+1, flux[i]);
        hEmeanCosmics->SetBinContent(binPhi[i]+1, binCos[i]+1, Emean[i]);
    }

    gStyle->SetOptStat(0);


    // Draw histograms
    TCanvas *canvas = new TCanvas("canvas1", "Cosmic Ray Histograms", 1200, 600);
    canvas->Divide(2, 1);
    canvas->cd(1)->SetRightMargin(1);
    canvas->cd(2)->SetRightMargin(1);
    canvas->cd(1);
    hFluxCosmics->Draw("colz");
    canvas->cd(2);
    hEmeanCosmics->Draw("colz");
    canvas->Draw();
    canvas->SaveAs("HyperKMuonHist.jpg");

    //############# PROJECTIONS ##################

    TCanvas* ProjectionCanvas = new TCanvas("ProjectionCanvas", "Muon Histogram Projections");
    ProjectionCanvas -> Divide(2, 2);
    
    ProjectionCanvas -> cd(1);

    TH1D *projY = hFluxCosmics->ProjectionY("Muon Flux cos theta projection");
    projY->GetYaxis()->SetTitle("Flux [s^{-1}m^{-2}]");
    projY->SetMinimum(0);
    projY -> Draw("");
        
    ProjectionCanvas -> cd(2);

    TH1D* projX = hFluxCosmics->ProjectionX("Muon Flux Phi projection");
    projX->GetYaxis()->SetTitle("Flux [s^{-1}m^{-2}]");
    projX->SetMinimum(0);
    projX -> Draw("");

    ProjectionCanvas-> cd(3);

    TH1D* EmeanprojY = hEmeanCosmics -> ProjectionY("Muon mean energy on cos theta");
    EmeanprojY->GetYaxis()->SetTitle("Emean [MeV]");
    EmeanprojY->SetMinimum(0);
    EmeanprojY->Draw("");

    ProjectionCanvas->cd(4);

    TH1D* EmeanprojX = hEmeanCosmics->ProjectionX("Muon mean energy on phi");
    EmeanprojX->GetYaxis()->SetTitle("Emean [MeV]");
    EmeanprojX->SetMinimum(0);
    EmeanprojX->Draw("");

    ProjectionCanvas->SaveAs("Projections.pdf");
    ProjectionCanvas->SaveAs("Projections.jpg");


}


void SeaLevel2DHist(){
    int Can_number = 10;

    TCanvas *muoncanvas[Can_number];

    muoncanvas[0] = new TCanvas("Muon canvas", "Muon Canvas", 800, 600);
    muoncanvas[0]->cd();


    TH2D *MuonHist[2];
    MuonHist[0] = new TH2D("SeaLevel2DHist", "MuonHist", 100, 0, 1, 1000, 1, 10);
    for (int i = 1; i <= (MuonHist[0]->GetNbinsX()); ++i) {
        double cos_theta = MuonHist[0]->GetXaxis()->GetBinCenter(i);

        for (int j = 1; j <= MuonHist[0]->GetNbinsY(); ++j) {
            double E_mu = MuonHist[0]->GetYaxis()->GetBinCenter(j);
            double dN_dE_dOmega = 0.14 * pow(E_mu, -2.7) * (1 / (1 + (1.1 * E_mu * cos_theta)/115) +  0.054 / (1 + (1.1 * E_mu * cos_theta / 850)));
            //std::cout << "E_mu is " << E_mu << "\t" << "cos theta is " << cos_theta << "\t" << "Flux is " << dN_dE_dOmega << endl;
            MuonHist[0]->SetBinContent(i, j, dN_dE_dOmega);
        }
    }

    MuonHist[0]->GetYaxis()->SetTitle("E_{#mu} [GeV]");
    MuonHist[0]->GetXaxis()->SetTitle("cos(#theta)");

    MuonHist[0]->Draw("colz");

    muoncanvas[0]->Draw();


    // muoncanvas[1] = new TCanvas("cos theta phi canvas", "cos theta phi canvas", 80000000);
    // muoncanvas[1]->Divide(2);
    // muoncanvas[1]->cd(1);

    // MuonHist[1] = new TH2D("Sea level muon flux", "Sea level muon flux", 180, 0, 360, 100, 0, 1);
    // MuonHist[2] = new TH2D("Sea level muon energy", "Sea level muon energy", 180, 0, 360, 100, 0, 1);


    // for (int i = 0; i<180; i++){
    //     int Phi = 2*i;
    //     for (int j = 0; j < 100; ++j){
    //         double cos_theta = j;
    //         double dN_dE_dOmega = dN_dE_dOmega_value()
    //     }
    //  }    

    }
    
double dN_dE_dOmega_value(double cos_theta, double E_mu){
    double dN_dE_dOmega = 0.14 * pow(E_mu, -2.7) * (1 / (1 + (1.1 * E_mu * cos_theta)/115) +  0.054 / (1 + (1.1 * E_mu * cos_theta / 850)));
    return dN_dE_dOmega;
}
    
double dN_dE_dOmega_value2(double cos_theta, double E_mu){
    double dN_dE_dOmega = 1400 * (1 / (1 + (1.1 * E_mu * cos_theta)/115) +  0.054 / (1 + (1.1 * E_mu * cos_theta / 850)));
    return dN_dE_dOmega;
}    

TH1D* Histogram_from_theta(double cos_theta){
    TH1D *histogram = new TH1D("Histogram_from_theta", ";Energy [GeV];#frac{dN}{dEd#Omega} [m^{-2}s^{-1}sr^{-1}GeV^{1.7}]", 10000, 1, 1000);
    for (int i = 0; i < histogram->GetNbinsX(); ++i){
        double dN_dE_dOmega  = dN_dE_dOmega_value2(cos_theta, i);
        histogram->SetBinContent(i, dN_dE_dOmega);
    }
    return histogram;
     

}





void SeaLevel1DHist(){

    std::vector<TH1D*> MuonHists;

    TString CanvasName = "name";
    TString CanvasTitle = "title";
    TCanvas *MuonCanvas = new TCanvas(CanvasName, CanvasTitle, 800, 600);
    THStack *MuonHS = new THStack("enery on cos theta", ";Energy[GeV];Flux");

    MuonCanvas->cd();

    double help = 0.1;
    while (help <= 1){
        TString HistName = "Cos_theta =" + std::to_string(help); 
        TString HistTitle = "Title " + std::to_string(help); 
        MuonHists.push_back(new TH1D(HistName, HistTitle + ";Energy[GeV];Flux" , 1000, 1, 10));
        help += 0.05;
        
    } 

    double cos_theta = 0.10;

    for (int i = 0; i < MuonHists.size(); ++i){

        for (int j = 1; j <= MuonHists[i]->GetNbinsX(); ++j){
            double E_mu = MuonHists[i]->GetXaxis()->GetBinCenter(j);
            
            auto dN_dE_dOmega = dN_dE_dOmega_value(cos_theta, E_mu);

            MuonHists[i]->SetBinContent(j, dN_dE_dOmega);
        }
        cos_theta += 0.05;

    }
    for (int i = 0; i < MuonHists.size(); ++i){
        MuonHS->Add(MuonHists[i]);
    }

    MuonHS->Draw("NOSTACk");

    MuonCanvas->Update();
    MuonCanvas->Draw();

    
    // ################# PDG FUNCTION ########################

    TCanvas *pdgcanvas = new TCanvas("pdgcanvas", "pdgcanvas", 800, 600);
    pdgcanvas->cd();
    pdgcanvas->SetLogx();
    pdgcanvas->SetLogy();

    double theta1 = 0;
    double theta2deg = 75;
    double theta2rad = theta2deg * M_PI / 180;
    double theta2 = std::cos(theta2rad);
    TH1D *pdghist  = Histogram_from_theta(theta1);
    TH1D *pdghist2 = Histogram_from_theta(theta2);
    
    //pdghist2->Scale(2);
    pdghist2->SetLineColor(kRed);
    pdgcanvas->SetGrid();
    pdghist->Draw("hist W");
    pdghist2->Draw("same");

    TLegend *legenda = new TLegend(0.1, 0.1, 0.4, 0.3);
    legenda -> AddEntry(pdghist, "cos_theta = 0", "l");
    legenda -> AddEntry (pdghist2, "cos_theta = 0.75", "l");   

    legenda->Draw();

    pdgcanvas->Draw();

    pdgcanvas->SaveAs("pdgpicture.pdf");
    pdgcanvas->SaveAs("pdgpicture.jpg");

    
}

void flux_on_theta(){
    TCanvas *flux_on_theta_canvas = new TCanvas("lol", "flux_on_theta_canvas");
    TH1D *flux_on_theta_hist = new TH1D("flux_on_theta_hist", "flux_on_theta_hist;cos#theta;dN/dEd#Omega", 100, 0, 1);

    double E_mu = 4;
    double cos_theta_max = 1;
    double cos_theta = 0;
    int i = 1;
    while(cos_theta < cos_theta_max){
        double dN_dE_dOmega = dN_dE_dOmega_value(cos_theta, E_mu);
        flux_on_theta_hist->SetBinContent(i, dN_dE_dOmega);
        i++;
        cos_theta +=0.01;
    }
    flux_on_theta_hist->Draw("hist");
    flux_on_theta_canvas->Draw();
    flux_on_theta_canvas->SaveAs("ThetaFlux.pdf");
    flux_on_theta_canvas->SaveAs("ThetaFlux.jpg");


}

void MuonDataFit(const char* file_1, const char* file_2){
    std::ifstream file1(file_1);
    std::ifstream file2(file_2);

    if (!file1.is_open() && !file2.is_open()){
        std:cerr << "Error: unable to open a file" << endl;
        return;
    }

    std::vector<DataPoint> data_mu_plus;
    std::vector<DataPoint> data_mu_minus;
    double x, y, y_err;

    while (file1 >> x >> y >> y_err) {
        DataPoint point;
        point.x = x;
        point.y = y;
        point.y_error = y_err;
        data_mu_plus.push_back(point);
    }
    while (file2 >> x >> y >> y_err) {
        DataPoint point;
        point.x = x;
        point.y = y;
        point.y_error = y_err;
        data_mu_minus.push_back(point);
    }

    TGraphErrors *mu_plus = new TGraphErrors((data_mu_plus.size()));
    for (size_t i = 0; i < data_mu_plus.size(); ++i){
        mu_plus->SetPoint(i, data_mu_plus[i].x, data_mu_plus[i].y /** pow(data_mu_plus[i].x, 2.7)*/);
        mu_plus->SetPointError(i, 0, data_mu_plus[i].y_error);
    }
    TGraphErrors *mu_minus = new TGraphErrors((data_mu_minus.size()));
    for (size_t i = 0; i < data_mu_minus.size(); ++i){
        mu_minus->SetPoint(i, data_mu_minus[i].x, data_mu_minus[i].y * pow(data_mu_minus[i].x, 2.7));
        mu_minus->SetPointError(i, 0, data_mu_minus[i].y_error);
    }
    TCanvas *mucan = new TCanvas("canvas", "canvas");
    //mucan -> SetLogx();
    //mucan -> SetLogy();
    mu_plus -> SetTitle("");
    mu_plus -> GetXaxis()->SetTitle("p[GeV/c]");
    mu_plus -> GetYaxis()->SetTitle("#mu+ Flux x [m^{-2}sr^{-1}sec^{-1}(GeV/c)^{1.7}]");
    
    TString fit_fun = "[0] * pow(x, -([1] * pow(log10(x), 3) + [2] * pow(log10(x), 2) + [3] * pow(log10(x), 1) + [4]))";
    TF1 *fit = new TF1("fit", fit_fun, 0, 500);
    fit -> SetParameters(0.0029, 0.0252, -0.263, 1.2743, 0.3061);
    fit -> SetNpx(10000);

    mu_plus->Fit(fit, "R");

    mu_plus->Draw("AP*");
    //mu_minus->Draw("P same");
    fit -> Draw("same L");
    TLegend* mu_leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    mu_leg->AddEntry(mu_plus,"data", "AP*");
    mu_leg->AddEntry(fit, "fit", "l");

    mu_leg -> Draw();
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);

    mucan->Draw();

    file1.close();
    file2.close();

}



int main() {
    // Specify the input file containing the tokenized data
    const char* filename = "MuonFlux-HyperK-ThetaPhi.dat";
    const char* muon_plus_file = "../../Data_files/mu_plus.txt";
    const char* muon_minus_file = "/../../Data_files/mu_minus.txt";


    // Call the function to create histograms from the data
    CreateHistogramFromData(filename);

    SeaLevel2DHist();
    SeaLevel1DHist();
    flux_on_theta();

    MuonDataFit(muon_plus_file, muon_minus_file);

    return 0;
}
