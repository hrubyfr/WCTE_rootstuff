#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include <iostream>
#include "TH2D.h"
#include "TAttPad.h"
#include "TH1D.h"
#include "TLegend.h"
#include <fstream>



class Reyna_params{
public: 
	double c1 = 0.00253;
	double c2 = 0.2455;
	double c3 = 1.288;
	double c4 = -0.2555;
	double c5 = 0.0209;
	double c_0 = 3e8;
	double m_mu = 0.105659;
};



Double_t Reyna_phi(Double_t *x, Double_t *par){
	Double_t E_mu = x[0];
	Double_t theta = par[0];
	Reyna_params params;
	Double_t p_mu = 1 / params.c_0 * sqrt(pow(E_mu, 2) - (pow(params.m_mu * pow(params.c_0, 2), 2)));

	Double_t zeta = p_mu * cos(theta);
	Double_t I_V = params.c1 * pow(zeta, -1 * (params.c2 + params.c3 * log10(zeta) + params.c4 * pow(log10(zeta), 2) + params. c5 * pow(log10(zeta), 3)));
	Double_t I = pow(cos(theta), 3) * I_V;

	//#### debugging ####
	std::cout << "P_mu is : " << p_mu << endl;
	std::cout << "Zeta is : " << zeta << endl;
	std::cout << "Vertical Flux is : " << I_V << endl;
	std::cout << "Final flux is : " << I << endl;

	return I;

}

Double_t x_Reyna_phi(Double_t *x, Double_t *par){

	Double_t E_mu = x[0];
	Double_t theta = par[0];
	Reyna_params params;
	Double_t p_mu = 1 / params.c_0 * sqrt(pow(E_mu, 2) - (pow(params.m_mu * pow(params.c_0, 2), 2)));

	Double_t zeta = p_mu * cos(theta);
	Double_t I_V = params.c1 * pow(zeta, -1 * (params.c2 + params.c3 * log10(zeta) + params.c4 * pow(log10(zeta), 2) + params. c5 * pow(log10(zeta), 3)));
	Double_t I = E_mu * pow(cos(theta), 3) * I_V;

	//#### debugging ####

	// std::cout << "Zeta is : " << zeta << endl;
	// std::cout << "Vertical Flux is : " << I_V << endl;
	// std::cout << "Final flux is : " << I << endl;

	return I;

}




void th2d_Reyna_flux(){
	Reyna_params params;
	auto canvas = new TCanvas();
	canvas -> cd();

	Double_t cos_theta = 1;
	Double_t theta = acos(cos_theta);
	Double_t E0 = 1000;


	


	
	TString flux_title = "cos ^{2} (#theta) Muon Flux";
	TString flux_xaxis_title = ";#phi [deg]";
	TString flux_yaxis_title = ";cos (#theta)";
	TString flux_zaxis_title = ";normallized flux [GeV^{-1} cm^{-2} s^{-1} sr^{-1}]";
	auto TH2hist_flux = new TH2D("", flux_title + flux_xaxis_title + flux_yaxis_title + flux_zaxis_title, 180, 0, 360, 100, 0, 1);
	for (int i = 0; i < TH2hist_flux->GetNbinsX(); i++){
		for (int j = 0; j < TH2hist_flux->GetNbinsY(); j++){
			TH2hist_flux->SetBinContent(i+1, j+1, j*j);
		}
	}
	TH2hist_flux->Scale(1.0/TH2hist_flux->Integral());
	canvas->SetRightMargin(0.15);

	TH2hist_flux->Draw("colz");
	canvas->SaveAs("TH2Flux.jpg");
	canvas->SaveAs("TH2Flux.pdf");




	canvas -> Update();

	TString title = "Muon Energy";
	TString xaxis_title = ";#phi [deg]";
	TString yaxis_title = ";cos (#theta)";

	auto energy_canvas = new TCanvas();
	energy_canvas->cd();
	energy_canvas->SetLogz();

	auto TH2hist = new TH2D("", title + xaxis_title + yaxis_title,180, 0, 360, 100, 0, 1);
	for (double j = 0; j < TH2hist->GetNbinsY(); j++){
		std::cout << "cos theta is : " << j / 100 << endl;
		double Theta = acos(j / 100);
		std::cout << "Theta is : " << Theta * TMath::RadToDeg() << endl;
		auto fun = new TF1("Reyna model muon energy spectrum", Reyna_phi, params.m_mu, 1000, 1);
		fun->SetParameter(0, Theta);
		fun->SetNpx(10000);

		auto weighted_fun = new TF1("", x_Reyna_phi, params.m_mu, 1000, 1);
		weighted_fun->SetParameter(0, Theta);
		weighted_fun->SetNpx(10000);

		auto fun_int = fun->Integral(0, 1000);
		auto weighted_fun_int = weighted_fun ->Integral(0, 1000);

		Double_t mean = 0;
		if (fun_int !=0){mean = weighted_fun_int / fun_int;}
		std::cout << "Mean is : " << mean << endl;

		for (double i = 0; i < TH2hist->GetNbinsX(); i++){
			

			TH2hist->SetBinContent(i+1, j+1, mean);
		}
	}

	canvas->SetRightMargin(0.1);

	TH2hist->GetZaxis()->SetTitle("p mean [GeV/c]");
	TH2hist->Draw("colz");

	canvas ->SaveAs("Energy_hist.pdf");
	canvas->SaveAs("Energy_hist.jpg");

	// OUTPUT HISTOGRAM TO DATA FILE

	std::ofstream output("E_mean_hist.dat");

	double cos_mean, cos_min, cos_max, phi_mean, phi_min, phi_max, flux, Emean;

	for (int i = 0; i < TH2hist->GetNbinsY(); ++i){
		cos_mean = TH2hist -> GetYaxis() -> GetBinCenter(i+1);
		cos_min = cos_mean - 0.005;
		cos_max	= cos_mean + 0.005;
		for (int j = 0; j < TH2hist->GetNbinsX(); ++j){
			phi_mean = TH2hist -> GetXaxis()->GetBinCenter(j+1);
			phi_min = phi_mean - 1;
			phi_max = phi_mean + 1;
			flux = TH2hist_flux -> GetBinContent(j+1, i+1);
			Emean = TH2hist -> GetBinContent(j+1, i+1);

			output << i << "\t" << j << "\t" << cos_mean << "\t" << cos_min << "\t" << cos_max << "\t" << phi_mean << "\t" << phi_min << "\t" << phi_max << "\t" << flux << "\t" << Emean << endl;
		}
	}
	output.close();
}



int cosmics(){
	
	// th2denergy();
	th2d_Reyna_flux();

	return 0;
}