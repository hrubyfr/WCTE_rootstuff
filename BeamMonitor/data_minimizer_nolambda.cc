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

double chi_squared_nolambda(const double *vars, const vector<double>& data, const vector<double>& data_errors, const vector<double>& length){
	const double v_eff = vars[0];
	const double sigma_sipm = vars[1];
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

void data_minimizer_nolambda(const char* filename){
	
	//CREATE HISTOGRAMS
	TString fname = filename;
	fname = fname(0, fname.Length() - 4);

	// #############################################################
	TString plots_folder = "analysis_plots/data_minimizer_plots_nolambda_" + fname;
	if(gSystem->AccessPathName(plots_folder))gSystem->Exec("mkdir -p " + plots_folder);

	std::ifstream file(filename);
	double sigma_tot, sigma_error, length;
	vector<double> sigmas, errors, lengths;
	unsigned int NData = 0;
	for (int i = 0; i < 8; i++){
		file >> sigma_tot >> sigma_error >> length;
		sigmas.push_back(sigma_tot);
		errors.push_back(sigma_error);
		lengths.push_back(length);
		NData++;
	}

	gSystem->cd(plots_folder);

	auto chi_lambda = [&](const double *x){
		return chi_squared_nolambda(x, sigmas, errors, lengths);
	};
	
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

	min->SetTolerance(10e-6);
	min->SetPrintLevel(1);
	min->SetStrategy(2);

	ROOT::Math::Functor func(chi_lambda, 2);
	min->SetFunction(func);

	min->SetVariable(0, "v_eff", 172, 0.1);
	min->SetVariable(1, "sigma_sipm", 0.2, 0.001);
	//min->SetVariableLimits(0, 0, 200);

	
	
	min->SetMaxFunctionCalls(1000000);
	min->SetMaxIterations(100000);

	min->Minimize();
	if (min->Status() != 0) {
		std::cerr << "Warning: Minimization did NOT converge! Status = " 
			<< min->Status() << endl;
	}


	const double* results = min->X();
	double fit_errors[2];
	for(int i = 0; i < 2; i++){
		fit_errors[i] = min->Errors()[i];
	}
	cout << "Best v_eff: " << results[0] << "\t best sigma_sipm: " << results[1] << endl;
	
	cout << "Chi square of the fit is " << min->MinValue() << endl;

	cout << "Chi^2/NDF of the fit is " << min->MinValue() / (NData - 2) << endl;



} // end of code
