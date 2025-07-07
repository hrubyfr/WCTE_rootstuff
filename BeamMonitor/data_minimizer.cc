#include "ostream"
#include "fstream"
#include "vector"

#include "math.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <TFile.h>
#include <fstream>
#include <vector>

using std::vector;
using std::cout;
using std::endl;

double chi_squared(const double *vars, const vector<double>& data, const vector<double>& data_errors, const vector<double>& length){
	const double v_eff = vars[0];
	const double sigma_sipm = vars[1];
	double chi2 = 0;
	for (int i = 0; i < data.size(); i++){
		if (i%8 == 1 || i%8 == 2 || i%8 == 5 || i%8 == 6) continue;
		double sigma_sq_meas = data[i] * data[i];
		double dim = length[i];

		double sigma_sq_pred = dim * dim / (3 * v_eff * v_eff) + /*vars[i%8+2] * */ sigma_sipm * sigma_sipm;
		double value = pow((sigma_sq_pred - sigma_sq_meas) / (data_errors[i] * data_errors[i]), 2) ;
		chi2 += value;
	}
	return chi2;
}

void data_minimizer(const char* filename){
	std::ifstream file(filename);
	double sigma_tot, sigma_error, length;
	vector<double> sigmas, errors, lengths;
	unsigned int NData = 0;
	while (file >> sigma_tot >> sigma_error >> length){
		sigmas.push_back(sigma_tot);
		errors.push_back(sigma_error);
		lengths.push_back(length);
		NData++;
	}

	auto chi_lambda = [&](const double *x){
		return chi_squared(x, sigmas, errors, lengths);
	};
	
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

	min->SetTolerance(0.01);
	min->SetPrintLevel(1);
	min->SetStrategy(2);

	ROOT::Math::Functor func(chi_lambda, 2);
	min->SetFunction(func);

	min->SetVariable(0, "v_eff", 172, 0.1);
	//min->SetVariableLimits(0, 0, 200);

	min->SetVariable(1, "sigma_sipm", 0.267, 0.0001);
	/*
	for (int i = 0; i < 8; i++){
		if (i == 2) min->SetVariable(i+2, Form("lambda_%i", i), 2.0, 0.001);
		else min->SetVariable(i+2, Form("lambda_%i", i), 1.0, 0.001);
		min->SetVariableLimits(i+2, 1.0, 10.0);
	}
	*/
	min->SetMaxFunctionCalls(1000000);
	min->SetMaxIterations(100000);

	min->Minimize();
	if (min->Status() != 0) {
		std::cerr << "Warning: Minimization did NOT converge! Status = " 
			<< min->Status() << endl;
	}
	min->Hesse();


	const double* results = min->X();
	cout << "Best v_eff: " << results[0] << "\t best sigma_sipm: " << results[1] << endl;    

	unsigned int nsteps = 100;
	double v[nsteps];
	double chi[nsteps];

	double MinValue = min->MinValue();

	cout << min->MinValue() << endl;
	min->Scan(0, nsteps, v, chi, 160, 200);

	cout << "Chi^2/NDF = " << MinValue/(NData-10) << endl;
}
