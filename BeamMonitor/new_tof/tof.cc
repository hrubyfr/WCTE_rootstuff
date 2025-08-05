#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "math.h"

#include <TCanvas.h>
#include "TFile.h"
#include "TH2D.h"
#include <TGraph.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TArrow.h>

#include <cstdio>
#include <vector>
#include <array>

using std::cout;
using std::endl;
using std::array;
using std::vector;

array<double, 2> construct_PMT(double x, double y){
	array<double, 2> arr = {x, y};
	return arr;
}

double distance(double x, double y){
	return sqrt(x * x + y * y);
}

//create a grid of points 
vector<array<double, 2>> get_grid(double size, int max_enum){
	vector<array<double, 2>> helpvector;
	for(int i = 0; i <= max_enum; i++){
		double y = size/2.0 - (size / max_enum) * i;
		for (int j = 0; j <= max_enum; j++){
			double x = -size/2 + j * (size / max_enum);
				array<double, 2> helppoint{x, y};
				helpvector.push_back(helppoint);
			}
		}	
	return helpvector;

}


double Det_response(array<double, 2> PMT_position, array<double, 2> position){
	array<double, 2> rel_position = {PMT_position[0] - position[0], PMT_position[1] - position[1]};
	// std::cout << "Relative position of point to detector is: " << rel_position << endl;
	// std::cout << "Distance of point from detector is " << distance(rel_position.X(), rel_position.Y()) << " cm" << endl;
	double r = distance(rel_position[0], rel_position[1]);
	auto response = 1 / (r * r);
	return response;
}

//################################################ MAIN FUNCTION ###########################################

void tof(){
	int NDet = 4;
	double size = 7.0;
	vector<array<double, 2>> PMTs;
	gSystem->Exec("mkdir -p ./T1_plots/");

	vector<array<double, 2>> PMT_positions = {{-size/2.0, size/4.0}, {size/2.0, size/4.0}, {-size/2.0, -size/4.0}, {size/2.0, -size/4.0}};
	
	TRandom3 rand;
	
	vector<array<double, 2>> new_coordinates;  //data arrays definition
	std::vector<double> chi2_values;
	


	int n_enum = 25;
	auto points = get_grid(size, n_enum);	//Generating grid


	int n_crash = 0;
	double blur = 1;			//defining parameters used in minimization and signal blurring
	double sigma = 0.001;

	for(int i = 0; i < points.size(); i++){
		auto point = points[i];
		std::vector<double> signals;
		for (int j = 0; j < PMT_positions.size(); j++){
			auto signal = Det_response(PMT_positions[j], point);
			//signal = signal + sigma * rand.Gaus(0, blur); // blurring the signal
			signals.push_back(signal);
		}
		auto chi2function = [&](const double* params) -> double{	//inline function declaration of chi2 function, that gets minimized 
			double chi2 = 0;
			double x = params[0];
			double y = params[1];
			for (int i = 0; i < NDet; i++){
				double dx = x - PMT_positions[i][0];
				double dy = y - PMT_positions[i][1];			//defining the chi2 function	
				double dist = sqrt(dx * dx + dy * dy);

				double theory_sig = 1.0 / (dist * dist);

				double residual = (signals[i] - theory_sig);
				chi2 += residual * residual;

			}
			// cout << "FUMILI was used" << endl;
			return chi2; 
		};

		ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

		minimizer->SetTolerance(0.001);
		minimizer->SetPrintLevel(1);

		ROOT::Math::Functor func(chi2function, 2);

		minimizer->SetFunction(func);

		minimizer->SetVariable(0, "x", 0.0, 0.001);
		minimizer->SetVariable(1, "y", 0.0, 0.001);
		minimizer->SetVariableLimits(0, -size/2.0, size/2.0);
		minimizer->SetVariableLimits(1, -size/2.0, size/2.0);
				
		minimizer->SetMaxFunctionCalls(10000);		

		minimizer->Minimize();			//Minimizing the chi2
		if (minimizer->Status() != 0) {
			std::cerr << "Warning: Minimization did NOT converge! Status = " 
				<< minimizer->Status() << std::endl;
			n_crash++;
		}
		double chi2value = minimizer->MinValue();				// Reading out the minimized values
		chi2_values.push_back(chi2value);		


		const double* results = minimizer->X();
		//cout << "Primary point X: " << points[i][0] << " Y: " << points[i][1] << endl;
		//cout << "Best x: " << results[0] << "\t best y: " << results[1] << endl;    
		
		std::array<double, 2> help_array = {{results[0], results[1]}};
		new_coordinates.push_back(help_array);

		
	}	//end of event loop

//##################################################### drawing plots #####################################################

	cout << "Minimization did not converge " << n_crash << " times out of " << new_coordinates.size() << "points"  << endl;
	
	
	TGraph* Primary_grid = new TGraph(points.size());
	Primary_grid->SetName("Primary_grid");
	Primary_grid->SetTitle("Primary grid;x [cm];y [cm]");
	for(int i = 0; i < points.size(); i++){
		Primary_grid->SetPoint(i, points[i][0], points[i][1]);
	}
	
	TGraph* Secondary_grid = new TGraph(new_coordinates.size());

	Secondary_grid->SetName("Secondary_grid");
	Secondary_grid->SetTitle("Grid after reconstruction through minimization;x [cm]; y[cm]");
	for(int i = 0; i < new_coordinates.size(); i++){
		Secondary_grid->SetPoint(i, new_coordinates[i][0], new_coordinates[i][1]);
	}
	
	TCanvas* c_grid = new TCanvas("", "", 1800, 900);
	c_grid->Divide(2);
	c_grid->cd(1);
	Primary_grid->Draw("AP*");	
	c_grid->cd(2);
	Secondary_grid->Draw("AP*");
	
	c_grid->Print("T1_plots/g_grid.png");
	delete c_grid;
	
	// ################### ARROW CANVAS #########################
	TCanvas* c_arrow = new TCanvas("", "", 1800, 900);
	
	
	TGraph* g_arrow = new TGraph(new_coordinates.size());
	
	g_arrow->SetName("g_arrow");
	g_arrow->SetTitle("point movement after reconstruction;x [cm]; y[cm]");

	for(int i = 0; i < new_coordinates.size();i++){
		g_arrow->SetPoint (i, new_coordinates[i][0], new_coordinates[i][1]);
	}
	g_arrow->Draw("ap*");
	for (int i = 0 ; i < points.size(); i++){
		TArrow* arrow = new TArrow(points[i][0], points[i][1], new_coordinates[i][0], new_coordinates[i][1], 0.005, "|>");
		arrow->SetLineColor(kBlue);
		arrow->Draw();
		}
	c_arrow->Print("T1_plots/g_arrow.png");
	delete c_arrow;




} //end of code
