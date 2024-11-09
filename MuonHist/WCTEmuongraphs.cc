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

class Parameters{
public:
	double A = 0.002382;
	double r = 0.76;
	double a = 2.5;
	double y0 = 1000;
	double gamma = 8.0/3.0;
	double b_mu = 0.800;
	double m_mu = 105.659;
	double t_mu0 = 2.2e-6;
	double rho_0 = 0.00123;
	double c = 3e10;
	double lambda_pi = 120;
	double b = 0.771;
	double tau_0 = 2.61e-8;
	double m_pi = 139.580;
	double j_pi = 148160;
};

class Reyna_params{
public: 
	double c1 = 0.00253;
	double c2 = 0.2455;
	double c3 = 1.288;
	double c4 = -0.2555;
	double c5 = 0.0209;
};



Double_t phi_s(/*Double_t *x, Double_t *par*/Double_t E0, Double_t theta){

	// Double_t E0 = x[0];
	// Double_t theta = par[0];
	Parameters params;

	Double_t B_mu = params.b_mu * params.m_mu * params.y0 / (params.c * params.t_mu0 * params.rho_0);
	std::cout << "B_mu is : " << B_mu << endl;
	Double_t E_pi = (1 / params.r) * (E0 + params.a * params.y0 * (1 / cos(theta) - 0.1));
	std::cout << "E_pi is : " << E_pi << endl;
	Double_t P_mu = pow(0.1 * cos(theta) * (1 - (params.a *(params.y0 * (1/cos(theta)) - 100) / (params.r * E_pi) )), B_mu / ((params.r * E_pi + 100 * params.a) * cos(theta)));
	std::cout << "P_mu is : " << P_mu << endl;
	Double_t phis = params.A * pow(E_pi, -params.gamma) * P_mu * params.lambda_pi * params.b * params.j_pi / (E_pi * cos(theta) + params.b * params.j_pi);
	std::cout << "Phi_s is : " << phis << endl;
	return phis;
}	

Double_t Reyna_phi(Double_t *x, Double_t *par){
	Double_t p_mu = x[0];
	Double_t theta = par[0];
	Reyna_params params;

	Double_t zeta = p_mu * cos(theta);
	Double_t I_V = params.c1 * pow(zeta, -1 * (params.c2 + params.c3 * log10(zeta) + params.c4 * pow(log10(zeta), 2) + params. c5 * pow(log10(zeta), 3)));
	Double_t I = pow(cos(theta), 3) * I_V;

	//#### debugging ####

	// std::cout << "Zeta is : " << zeta << endl;
	// std::cout << "Vertical Flux is : " << I_V << endl;
	// std::cout << "Final flux is : " << I << endl;

	return I;

}

Double_t x_Reyna_phi(Double_t *x, Double_t *par){
	Double_t p_mu = x[0];
	Double_t theta = par[0];
	Reyna_params params;

	Double_t zeta = p_mu * cos(theta);
	Double_t I_V = params.c1 * pow(zeta, -1 * (params.c2 + params.c3 * log10(zeta) + params.c4 * pow(log10(zeta), 2) + params. c5 * pow(log10(zeta), 3)));
	Double_t I = p_mu * pow(cos(theta), 3) * I_V;

	//#### debugging ####

	// std::cout << "Zeta is : " << zeta << endl;
	// std::cout << "Vertical Flux is : " << I_V << endl;
	// std::cout << "Final flux is : " << I << endl;

	return I;

}


void th2dflux(){
	gStyle->SetOptStat(0);
	auto canvas = new TCanvas();
	canvas->cd();

	TString title = "cos ^{2} (#theta) Muon Flux";
	TString xaxis_title = ";#phi [deg]";
	TString yaxis_title = ";cos (#theta)";
	TString zaxis_title = ";normallized flux [GeV^{-1} cm^{-2} s^{-1} sr^{-1}]";
	auto TH2hist = new TH2D("", title + xaxis_title + yaxis_title + zaxis_title, 180, 0, 360, 100, 0, 1);
	for (int i = 0; i < TH2hist->GetNbinsX(); i++){
		for (int j = 0; j < TH2hist->GetNbinsY(); j++){
			TH2hist->SetBinContent(i+1, j+1, j*j);
		}
	}
	TH2hist->Scale(1.0/TH2hist->Integral());
	canvas->SetRightMargin(0.15);

	TH2hist->Draw("colz");
	canvas->SaveAs("TH2Flux.jpg");
	canvas->SaveAs("TH2Flux.pdf");

}

void th2denergy(){
	auto canvas = new TCanvas();
	canvas->cd();
	Double_t cos_theta = 1;
	Double_t theta = acos(cos_theta);
	Double_t E0 = 1;

	std::cout << phi_s(E0, theta) << endl;

	canvas->SetLogx();
	canvas->SetLogy();



	// auto fun = new TF1("fun", phi_s, 0, 1000, 1);
	// fun->SetParameter(0, theta);
	// fun->SetNpx(10000);

	// fun->GetYaxis()->SetRangeUser(1e-11, 0.01);
	// TString title = "cos ^{2} (#theta) Muon Energy";
	// TString xaxis_title = ";#phi [deg]";
	// TString yaxis_title = ";cos (#theta)";

	// fun->Draw("L");
// 	auto TH2hist = new TH2D("", title + xaxis_title + yaxis_title, 180, 0, 360, 100, 0, 1);
// 	for (int i = 0; i < TH2hist->GetNbinsX(); i++){
// 		for (int j = 0; j < TH2hist->GetNbinsY(); j++){

// 			TH2hist->SetBinContent(i+1, j+1, j*j);
// 		}
// 	}
// 	TH2hist->Scale(1.0/TH2hist->Integral());
// 	canvas->SetRightMargin(0.15);

// 	TH2hist->Draw("colz");
// 	canvas->SaveAs("TH2Energy.jpg");
// 	canvas->SaveAs("TH2Energy.pdf");

}


void th2d_Reyna_flux(){
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




	Double_t par[1] = {0 * TMath::DegToRad()};
	auto hist = new TH1D("", "", 1000, 0.1, 1000);
	for (int i = 0; i < hist->GetNbinsX(); i++){
		Double_t x[1] = {hist -> GetBinCenter(i+1)};
		Double_t value = Reyna_phi(x, par);

		hist -> SetBinContent(i+1, value);
	}

	
	double hist_mean = hist ->GetMean();
	std::cout << "hist mean is : " << hist_mean << endl;
	
	auto test_canvas = new TCanvas();
	test_canvas->cd();
	
	auto test_fun1 = new TF1("", Reyna_phi, 0, 1000, 1);
	auto test_theta1 = 0;
	test_fun1->SetParameter(0, test_theta1);

	auto test_fun2 = new TF1("", Reyna_phi, 0, 1000, 1);
	auto test_theta2 = 60 * TMath::DegToRad();	
	test_fun2->SetParameter(0, test_theta2);
	// auto legend = new TLegend(0.6, 0.7, 0.9, 0.9);
	// legend -> AddEntry(fun, "#theta = 0^{o}", "l");
	// legend -> AddEntry(fun2, "#theta = 60^{o}", "l");
	// legend -> Draw("same");

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
		auto fun = new TF1("Reyna model muon energy spectrum", Reyna_phi, 0, 1000, 1);
		fun->SetParameter(0, Theta);
		fun->SetNpx(10000);

		auto weighted_fun = new TF1("", x_Reyna_phi, 0, 1000, 1);
		weighted_fun->SetParameter(0, Theta);
		weighted_fun->SetNpx(10000);

		auto fun_int = fun->Integral(0, 1000);
		auto weighted_fun_int = weighted_fun ->Integral(0, 1000);

		Double_t mean = 0;
		if (fun_int !=0){mean = weighted_fun_int / fun_int;}
		Double_t mean_energy = sqrt(pow(mean, 2) + pow(0.105659, 2));  
		std::cout << "Mean energy is : " << mean << endl;

		for (double i = 0; i < TH2hist->GetNbinsX(); i++){
			

			TH2hist->SetBinContent(i+1, j+1, mean);
		}
	}

	canvas->SetRightMargin(0.1);

	TH2hist->GetZaxis()->SetTitle("p mean [GeV/c]");
	TH2hist->Draw("colz");

	canvas ->SaveAs("Energy_hist2.pdf");
	canvas->SaveAs("Energy_hist2.jpg");

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

			output << i << " " << j << " " << cos_mean << " " << cos_min << " " << cos_max << " " << phi_mean << " " << phi_min << " " << phi_max << " " << flux << " " << Emean << endl;
		}
	}
	output.close();
}



int WCTEmuongraphs(){
	th2dflux();
	// th2denergy();
	th2d_Reyna_flux();

	return 0;
}