#include "ostream"
#include "fstream"
#include "vector"

#include "math.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <Rtypes.h>
#include <TFile.h>
#include <TPad.h>
#include <TString.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>

#include "TGraph.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

using std::vector;
using std::cout;
using std::endl;

double chi_squared(const double *vars, const vector<double>& data, const vector<double>& data_errors, const vector<double>& length){
	const double v_eff = vars[0];
	const double sigma_sipm = 0.267;
	double chi2 = 0;
	for (int i = 0; i < data.size(); i++){
		//if (i%8 == 1 || i%8 == 2 || i%8 == 5 || i%8 == 6) continue;
		double sigma_sq_meas = data[i];
		double dim = length[i];

		double sigma_sq_pred = sqrt(dim * dim / (3 * v_eff * v_eff) + vars[i%8+1] * sigma_sipm * sigma_sipm);
		double value = pow((sigma_sq_pred - sigma_sq_meas) / (data_errors[i]), 2) ;
		chi2 += value;
	}
	return chi2;
}

double chi_squared_nolambda(const double *vars, const vector<double>& data, const vector<double>& data_errors, const vector<double>& length){
	const double v_eff = vars[0];
	const double sigma_sipm = 0.267;
	double chi2 = 0;
	for (int i = 0; i < data.size(); i++){
		//if (i%8 == 1 || i%8 == 2 || i%8 == 5 || i%8 == 6) continue;
		double sigma_sq_meas = data[i];
		double dim = length[i];

		double sigma_sq_pred = sqrt(dim * dim / (3 * v_eff * v_eff) + sigma_sipm * sigma_sipm);
		double value = pow((sigma_sq_pred - sigma_sq_meas) / (data_errors[i]), 2) ;
		chi2 += value;
	}
	return chi2;
}

double chi_nuisance(const double *vars, const vector<double>& lambda_errors, const vector<double>& lambdas, const vector<double>& data, const vector<double>& data_errors, const vector<double>& length){
	const double v_eff = vars[0];
	const double sigma_sipm = 0.267;
	
	double chi2 = 0;
	for (int i = 0; i < data.size(); i++){
		//if (i%8 == 1 || i%8 == 2 || i%8 == 5 || i%8 == 6) continue;
		double sigma_sq_meas = data[i];
		double dim = length[i];

		double sigma_sq_pred = sqrt(dim * dim / (3 * v_eff * v_eff) + vars[i%8+1] * sigma_sipm * sigma_sipm);
		double sum = 0;
		for (int j = 0; j < 8; j++){
			double lambda = lambdas[j];
			double error = lambda_errors[j];
			sum += log((1/(error * sqrt(2 * M_PI))) * exp(-pow(vars[j+1] - lambda, 2) / (2 * error * error) )); 
		}
		double value = pow((sigma_sq_pred - sigma_sq_meas) / (data_errors[i]), 2) * sum;
		chi2 += value;
	}
	return chi2;
}
double chi_scan(const double v_eff, const vector<double>& lambdas, const vector<double>& data, const vector<double>& data_errors, const vector<double>& length){
	const double sigma_sipm = 0.267;
	double chi2 = 0;
	for (int i = 0; i < data.size(); i++){
		double sigma_sq_meas = data[i];
		double dim = length[i];

		double sigma_sq_pred = sqrt(dim * dim / (3 * v_eff * v_eff) + lambdas[i%8] * sigma_sipm * sigma_sipm);
		double value = pow((sigma_sq_pred - sigma_sq_meas) / (data_errors[i]), 2) ;
		chi2 += value;
	}
	return chi2;
}

void data_minimizer(const char* filename){
	
	//CREATE HISTOGRAMS
	TGraph2D* g_scan2D = new TGraph2D();
	g_scan2D->SetTitle("2D Scan of effective speed and #lambda_{0} values of minimization;v_eff;#lambda_{0}");
	
	TString fname = filename;
	fname = fname(0, fname.Length() - 4);

	// #############################################################
	TString plots_folder = "analysis_plots/data_minimizer_plots_" + fname;
	if(gSystem->AccessPathName(plots_folder))gSystem->Exec("mkdir -p " + plots_folder);

	std::ifstream file(filename);
	int run_number, energy;
	vector<int> run_numbers, energies;
	double sigma_tot, sigma_error, length;
	vector<double> sigmas, errors, lengths;
	unsigned int NData = 0;
	while (file >> run_number >> energy >> sigma_tot >> sigma_error >> length){
	
		if(NData%8 == 0){
		run_numbers.push_back(run_number);
		energies.push_back(energy);
		}

		sigmas.push_back(sigma_tot);
		errors.push_back(sigma_error);
		lengths.push_back(length);
		NData++;
	}

	gSystem->cd(plots_folder);

	auto chi_lambda = [&](const double *x){
		return chi_squared(x, sigmas, errors, lengths);
	};
	
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

	min->SetTolerance(10e-6);
	min->SetPrintLevel(1);
	min->SetStrategy(2);

	ROOT::Math::Functor func(chi_lambda, 9);
	min->SetFunction(func);

	min->SetVariable(0, "v_eff", 172, 0.1);
	//min->SetVariableLimits(0, 0, 200);

	
	for (int i = 0; i < 8; i++){
		if (i == 2) min->SetVariable(i+1, Form("lambda_%i", i), 2.0, 0.001);
		else min->SetVariable(i+1, Form("lambda_%i", i), 1.0, 0.001);
		min->SetVariableLimits(i+1, 0.0, 100.0);
	}
	
	min->SetMaxFunctionCalls(1000000);
	min->SetMaxIterations(100000);
	min->SetPrintLevel(1);

	min->Minimize();
	if (min->Status() != 0) {
		std::cerr << "Warning: Minimization did NOT converge! Status = " 
			<< min->Status() << endl;
	}


	const double* results = min->X();
	double fit_errors[9];
	for(int i = 0; i < 9; i++){
		fit_errors[i] = min->Errors()[i];
	}
	cout << "Best v_eff: " << results[0] << "\t smallest sigma_sipm: " << sqrt(results[5] * 0.267 * 0.267) << endl;    
	
	cout << "Size of pointer is: " << sizeof(results) << endl;

	vector<double> lambdas;
	vector<double> lambda_errors;

	for (int i = 1; i <= sizeof(results); i++){
		lambdas.push_back(results[i]);
		lambda_errors.push_back(fit_errors[i]);
	}

	for (int i = 0; i < 200; i++){
		for (int j = 0; j < 200; j++){
			double v_eff = 170.0 + 2.0/5.0 * i/2.0;
			double lambda_0 = 0.4 + 16.0/10.0 * j/200.0;
			lambdas[0] = lambda_0;
			double chi_value = chi_scan(v_eff, lambdas, sigmas, errors, lengths);
			g_scan2D->AddPoint(v_eff, lambda_0, chi_value);
		}
	}
	TCanvas* c_2Dscan = new TCanvas("c_2Dscan", "c_2Dscan", 1800, 900);
	g_scan2D->Draw("con1z lego2");
	c_2Dscan->Print("c_2D_scan.pdf");

	lambdas[0] = results[1];

	unsigned int nsteps = 100;
	double v[nsteps];
	double chi[nsteps];

	double MinValue = min->MinValue();

	cout << min->MinValue() << endl;
	min->Scan(0, nsteps, v, chi, 160, 200);

	TGraph* h_vscan = new TGraph(nsteps, v, chi);
	TCanvas* c_vscan = new TCanvas("c_vscan", "c_vscan", 1800, 900);
	cout << " SCAN RESULTS: " << endl;
	h_vscan->SetTitle("Scan of #chi^{2} values for different effective speeds;v_{eff};#chi^{2}");
	h_vscan->Draw("AC*");
	c_vscan->Print("v_speed_scan.pdf");

	cout << "Chi^2/NDF = " << MinValue/(NData-10) << endl;

	gStyle->SetOptStat(0);
	TCanvas* can = new TCanvas("can", "can", 1800, 900);
	int n_hists = sigmas.size() / 8;
	TH1D* hist = new TH1D("h_pred_data", "h_pred_data", 8, 0, 8);
	for (int i = 0; i < 8; i++){
		//cout << "lambda results: " << results[i+1] << " Lambda errors: " << fit_errors[i+1] << endl;
		double sigma_sq_pred = sqrt(lengths[i] * lengths[i] / (3 * results[0] * results[0]) + results[i + 1] * 0.267 * 0.267);
		hist->SetBinContent(i+1, sigma_sq_pred);
		double error = 1/sigma_sq_pred  * sqrt(pow(0.267, 4) / 4 * pow(fit_errors[i + 1], 2) + pow(lengths[i], 4) / (9 * pow(results[0], 6)) * pow(fit_errors[0], 2));
		cout << "Predicted sigma for scintillator: " << i << " Sigma_tot  = " << sigma_sq_pred << " +- " << error << endl;
		hist->SetBinError(i+1, error);
	}


	TH1D* hist_data[n_hists];
	for(int i = 0; i < n_hists; i++){
		hist_data[i] = new TH1D(Form("h_data_%i", i), Form("h_data_%i", i), n_hists * 8, 0, 8);
	}
	hist_data[0]->GetYaxis()->SetRangeUser(0, 0.7);
	hist_data[0]->SetTitle("Comparison of predicted and measured #sigma_{tot} values");
	hist_data[0]->GetXaxis()->SetTitle("Scintillator number");
	hist_data[0]->GetYaxis()->SetTitle("#sigma_{tot} [ns]");
	for (int histn = 0; histn < n_hists; histn++){
		for(int i = 0; i < 8; i++){
			hist_data[histn]->SetBinContent(histn + (n_hists*i+1), sigmas[8*histn + i]);
			hist_data[histn]->SetBinError(histn + (n_hists*i+1), errors[8*histn + i]);
		}
		hist_data[histn]->SetLineColor(kRed);
		if (histn == 0) hist_data[histn]->Draw("e1x0");
		else hist_data[histn]->Draw("same");
	}
	hist->SetLineWidth(2);
	hist->Draw("same histe");
	TLegend* leg  =new TLegend(0.7, 0.7, 0.9, 0.9);
	leg->AddEntry(hist_data[0], "data entries", "l");
	leg->AddEntry(hist, "Predicted value", "l");
	leg->Draw("same");

	can->Print("Data_vs_Fit.pdf");

	TCanvas *c_data = new TCanvas("c_data", "c_data", 1800, 900);
	vector<TH1D*> h_data(sigmas.size()/8);
	for(int i = 0; i < sigmas.size()/8; i++){
		h_data[i] = new TH1D(Form("h_data_only_%i", i), "histogram of data values with errors;scintillator;#sigma_{meas}", n_hists * 8, 0, 8);
		for(int j = 0; j < 8; j++){
			h_data[i]->SetBinContent(i + (n_hists*j+1), sigmas[8*i + j]);
			h_data[i]->SetBinError(i + (n_hists*j+1), errors[8*i + j]);
		}
		h_data[i]->SetLineColor(kPink + i);
		h_data[i]->SetLineWidth(2);
		if (i == 0) {
			h_data[i]->GetYaxis()->SetRangeUser(0, 0.8);
			h_data[i]->Draw("e1x0");

		}
		else h_data[i]->Draw("e1x0 same");
	}
	TLegend *l_data = new TLegend(0.7, 0.7, 0.9, 0.9);
	for(int i = 0; i < n_hists; i++){
		l_data->AddEntry(h_data[i], Form("Run %i, energy %i", run_numbers[i], energies[i]), "le");
	}
	l_data->Draw("same");

	c_data->Print("h_data.png");
	c_data->Print("h_data.pdf");

	ROOT::Math::Minimizer* min_nuisance = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

	auto chi_nuis = [&](const double *x){
		return chi_nuisance(x ,lambda_errors ,lambdas ,sigmas, errors, lengths);
	};

	/*
	   min_nuisance->SetTolerance(10e-6);
	   min_nuisance->SetPrintLevel(1);
	   min_nuisance->SetStrategy(2);

	   ROOT::Math::Functor func_nuisance(chi_nuis, 9);
	   min_nuisance->SetFunction(func_nuisance);

	   min_nuisance->SetVariable(0, "v_eff", 172, 0.1);
	//min_nuisance->SetVariableLimits(0, 0, 200);


	for (int i = 0; i < 8; i++){
	min_nuisance->SetVariable(i+1, Form("lambda_%i", i), lambdas[i], 0.001);
	min_nuisance->SetVariableLimits(i+1, 0.0, 100.0);
	}

	min_nuisance->SetMaxFunctionCalls(1000000);
	min_nuisance->SetMaxIterations(100000);

	min_nuisance->Minimize();
	if (min_nuisance->Status() != 0) {
	std::cerr << "Warning: Minimization did NOT converge! Status = " 
	<< min_nuisance->Status() << endl;
	}

*/


} // end of code
