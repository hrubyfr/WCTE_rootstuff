#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "math.h"

#include "TFile.h"
#include "TBox.h"
#include <Math/Vector2Dfwd.h>
#include <Rtypes.h>
#include <TGraph2D.h>
#include "TFile.h"
#include "TH2D.h"
#include "Math/Vector2D.h"
#include "Math/Point2D.h"
#include <TGraph.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TArrow.h>

#include <TVirtualPad.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <vector>
#include <array>

using std::array;
using std::vector;
using std::cout;
using std::endl;
using ROOT::Math::XYVector;

XYVector construct_SiPM(double x, double y){
	return XYVector(x, y);
}

ROOT::Math::XYPoint get_point(){
	return ROOT::Math::XYPoint(2.0, -0.5);
}
double distance(double x, double y){
	return sqrt(pow(x, 2) + pow(y, 2));
}

//create a grid of points 
vector<XYVector> get_grid(double radius, int max_enum){
	vector<XYVector> helpvector;
	for(int i = 0; i <= max_enum; i++){
		double y = radius - i * (1.0 / max_enum) * 2 * radius;
		for (int j = 0; j <= max_enum; j++){
			double x = -radius + j * (1.0 / max_enum) * 2 * radius;
			if (distance(x, y) < radius){
				XYVector helppoint(x, y);
				helpvector.push_back(helppoint);
			}
		}
	
	}
	return helpvector;

}


vector<vector<XYVector>> sort_points(vector<XYVector> points, vector<array<double, 2>> scint_dimensions, vector<double> start_y){
	vector<vector<XYVector>> sorted_points(9);
	for (int i = 0; i < points.size(); i++){
		if (abs(points[i].X()) < abs(scint_dimensions[0][0]/2) && points[i].Y() < start_y[0] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[0] - scint_dimensions[0][1]/2) sorted_points[0].push_back(points[i]);
		else if (abs(points[i].X()) < abs(scint_dimensions[1][0]/2) && points[i].Y() < start_y[1] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[1] - scint_dimensions[0][1]/2) sorted_points[1].push_back(points[i]);
		else if (abs(points[i].X()) < abs(scint_dimensions[2][0]/2) && points[i].Y() < start_y[2] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[2] - scint_dimensions[0][1]/2) sorted_points[2].push_back(points[i]);
		else if (abs(points[i].X()) < abs(scint_dimensions[3][0]/2) && points[i].Y() < start_y[3] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[3] - scint_dimensions[0][1]/2) sorted_points[3].push_back(points[i]);
		else if (abs(points[i].X()) < abs(scint_dimensions[4][0]/2) && points[i].Y() < start_y[4] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[4] - scint_dimensions[0][1]/2) sorted_points[4].push_back(points[i]);
		else if (abs(points[i].X()) < abs(scint_dimensions[5][0]/2) && points[i].Y() < start_y[5] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[5] - scint_dimensions[0][1]/2) sorted_points[5].push_back(points[i]);
		else if (abs(points[i].X()) < abs(scint_dimensions[6][0]/2) && points[i].Y() < start_y[6] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[6] - scint_dimensions[0][1]/2) sorted_points[6].push_back(points[i]);
		else if (abs(points[i].X()) < abs(scint_dimensions[7][0]/2) && points[i].Y() < start_y[7] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[7] - scint_dimensions[0][1]/2) sorted_points[7].push_back(points[i]);
		else sorted_points[8].push_back(points[i]);
	}//end for loop over points
	return sorted_points;
}//end sort point function



double Det_response(XYVector SiPM_position, XYVector position){
	XYVector Det_position(SiPM_position.X(), SiPM_position.Y());
	auto rel_position = Det_position.X() - position.X();
	// std::cout << "Relative position of point to detector is: " << rel_position << endl;
	// std::cout << "Distance of point from detector is " << distance(rel_position.X(), rel_position.Y()) << " cm" << endl;
	double r = rel_position;
	auto response = 1 / pow(r, 2);
	return response;
}

void draw_boxes(vector<TBox*> boxes){
	for(int i = 0; i < boxes.size(); i++){
		boxes[i]->SetFillStyle(0);
		boxes[i]->SetLineColor(kBlack);
		boxes[i]->Draw("same");
	}
}

double chi2plot(double x, XYVector point, XYVector SiPM1, XYVector SiPM2){
	double chi2 = 0;
	
	double dx1 = x - SiPM1.X();
	double dx2 = x - SiPM2.X();
	double theory_sig1 = 1.0 / (dx1 * dx1);
	double theory_sig2 = 1.0 / (dx2 * dx2);
	double signal1 = Det_response(SiPM1, point);
	double signal2 = Det_response(SiPM2, point);


	double residual1, residual2;
	residual1 = signal1 - theory_sig1;
	residual2 = signal2 - theory_sig2;

	chi2 += residual1 * residual1 + residual2 * residual2;

	return chi2; 

}

double GetMinimum_chi2(const double tolerance, const double signals[2], const XYVector SiPM1, const XYVector SiPM2){
	double x = (SiPM1.X() * signals[0] + SiPM2.X() * signals[1])/(signals[0] + signals[1]); //Start point guess from weighted charge		
	double step_size[5] = {10, 10e0, 10e-1, 10e-2, 10e-3};

	double chi2 = 0;
	double chi2_history = 0;
	auto evaluate = [x, SiPM1, SiPM2, signals, &chi2](){
		double dx1 = x - SiPM1.X();
		double dx2 = x - SiPM2.X();

		double theory_sig1 = 1.0 / (dx1 * dx1);
		double theory_sig2 = 1.0 / (dx2 * dx2);

		double residual1, residual2;					//Get chi2 at start point
		residual1 = signals[0] - theory_sig1;
		residual2 = signals[1] - theory_sig2;

		chi2 = residual1 * residual1 + residual2 * residual2;
	};
	evaluate();
	int i_step = 0;
	int direction;
	if (signals[0] > signals[1]) direction = -1;
	else direction = 1;
	while (chi2 > tolerance){
		x += step_size[i_step] * direction;
		chi2_history = chi2;
		evaluate();


		if (chi2 > chi2_history) {
			direction *=-1;
			i_step++;
		}
		


	}
	return x;
}

//################################################ MAIN FUNCTION ###########################################

void new_tof(double sigma = 0){

	int NDet = 16;
	vector<array<double, 2>> scint_dimensions = {{51, 16.25}, 
		{94, 16.25},
		{112, 16.25},
		{123, 16.25},
		{123, 16.25},
		{112, 16.25},
		{94, 16.25},
		{51, 16.25}
	};
	vector<double> start_y = {-58.875, -42.125, -25.375, -8.625, 8.625, 25.375, 42.125, 58.875};
	vector<array<double, 2>> det_positions;
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < NDet/2; j++){
			array<double, 2> position;
			if(i == 0) position = {scint_dimensions[j][0]/2, start_y[j]};
			if(i == 1) position = {-scint_dimensions[j][0]/2, start_y[j]};

			det_positions.push_back(position);
		}
	}//end for loop

	double radius = 75;
	std::vector<XYVector> SiPMs;
	for (int i = 0; i  < NDet; i++){
		SiPMs.push_back(construct_SiPM(det_positions[i][0], det_positions[i][1]));
	}

	TRandom3 rand;



	vector<vector<array<double, 2>>> rec_coordinates(8);  //data arrays definition
	vector<double> chi2_values;



	int n_enum = 50;
	auto points = sort_points(get_grid(radius, n_enum), scint_dimensions, start_y);		//Generating grid



	long npoints = 0;
	int n_crash = 0;
	double blur = 1;			//defining parameters used in minimization and signal blurring

	bool use_minuit = false;
	std::vector<std::vector<double>> signal_memory;
	for(int iscint = 0; iscint < points.size() - 1; iscint++){
		for (int ievent = 0; ievent < points[iscint].size(); ievent++){
			npoints++;
			auto point = points[iscint][ievent];
			std::vector<double> signals;


			signals.push_back(Det_response(SiPMs[iscint], point));
			signals.push_back(Det_response(SiPMs[iscint+8], point));
			signal_memory.push_back(signals);

			if (use_minuit){
			auto chi2function = [&](const double* params) -> double{	//inline function declaration of chi2 function, that gets minimized 
				double chi2 = 0;
				double x = params[0];

				if (abs(x) > abs(SiPMs[iscint].X()))  chi2+=10e6;
				//	double dx1 = x - SiPMs[iscint].X();
				//	double dy1 = y - SiPMs[iscint].Y();
				//	double dx2 = x - SiPMs[iscint+8].X();
				//	double dy2 = y - SiPMs[iscint+8].Y();
				//	double theory_sig1 = 1.0 / (dx1 * dx1 + dy1 * dy1);
				//	double theory_sig2 = 1.0 / (dx2 * dx2 + dy2 * dy2);
				double dx1 = x - SiPMs[iscint].X();
				double dx2 = x - SiPMs[iscint+8].X();
				double theory_sig1 = 1.0 / (dx1 * dx1);
				double theory_sig2 = 1.0 / (dx2 * dx2);

				double residual1, residual2;
				residual1 = 10e6 * (signals[0] - theory_sig1);
				residual2 = 10e6* (signals[1] - theory_sig2);

				chi2 += (residual1 * residual1 + residual2 * residual2)/2.0;


				// cout << "FUMILI was used" << endl;
				return chi2; 
			};

			ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

			minimizer->SetTolerance(10e-15);
			minimizer->SetPrintLevel(1);

			ROOT::Math::Functor func(chi2function, 1);

			minimizer->SetFunction(func);

			minimizer->SetVariable(0, "x", 0.0, 0.00001);
			//minimizer->SetVariable(1, "y", start_y[iscint], 0.001);

			minimizer->SetVariableLimits(0, SiPMs[iscint + 8].X(), SiPMs[iscint].X());
			//minimizer->SetVariableLimits(1, start_y[iscint] - 1, start_y[iscint] + 1);

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

			std::array<double, 2> help_array = {{results[0], start_y[iscint]}};
			rec_coordinates[iscint].push_back(help_array);

			cout << "Primary point: (" << point.X() << ", " << point.Y() << ")" << endl;
			cout << "reconstructed point: (" << results[0] << ", " << start_y[iscint] << ")" << endl;
			}
			else{
				double xyz = GetMinimum_chi2(1, signals, SiPMs[iscint], SiPMs[iscint+8]);
				double x_rec = GetMinimum_chi2(tolerance, signals, SiPMs[iscint], SiPMs[iscint + 8]);

			}
		}
	}	//end of event loop

	//##################################################### drawing plots #####################################################

	std::cout << "Minimization did not converge " << n_crash << " times " << std::endl;
	cout << points[8].size() << " Points did not hit the scintillator" << endl;
	TGraph* g_start = new TGraph(npoints);
	g_start->SetName("g_start");
	g_start->SetTitle("Primary position of points; x[mm]; y[mm]");
	TGraph* g_end = new TGraph(npoints);
	g_end->SetName("g_end");
	g_end->SetTitle("reconstructed position of points; x[mm]; y[mm]");
	int iter  = 0;
	cout << "There are " << npoints << " insiide the scintilators" << endl;
	for (int j = 0; j < points.size() - 1; j++){

		for(int i = 0; i < points[j].size(); i++){
			g_start->SetPoint(iter, points[j][i].X(), points[j][i].Y());
			g_end->SetPoint(iter, rec_coordinates[j][i][0], rec_coordinates[j][i][1]);
			iter++;
		}}//end for loop over first scintillator


	vector<TBox*> scints;
	for(int i = 0; i < scint_dimensions.size(); i++){
		TBox* help_box = new TBox(-scint_dimensions[i][0]/2, start_y[i] - scint_dimensions[i][1]/2, scint_dimensions[i][0]/2, start_y[i] + scint_dimensions[i][1]/2) ;
		scints.push_back(help_box);
	}//end create scintillator tiles
	TCanvas* can = new TCanvas("", "", 1800, 900);
	can->Divide(2);
	can->cd(1);
	g_start->GetXaxis()->SetRangeUser(-75, 75);
	g_start->GetYaxis()->SetRangeUser(-75, 75);

	g_start->Draw("AP*");
	draw_boxes(scints);
	can->cd(2);

	g_end->GetXaxis()->SetRangeUser(-75, 75);
	g_end->GetYaxis()->SetRangeUser(-75, 75);

	g_end->Draw("AP*");
	draw_boxes(scints);

	// ############ weighted charge #######################################

	//	TGraph* g_weight = new TGraph(npoints);
	//	iter = 0;
	//	for (int i = 0; i < points.size() - 1; i++){
	//		for(int j = 0; j < points[i].size(); j++) {
	//			double charge[2] = {Det_response(SiPMs[i], points[i][j]), Det_response(SiPMs[i+8], points[i][j])};
	//			XYVector pos[2] = {SiPMs[i] - points[i][j], SiPMs[i + 8] - points[i][j]};
	//			g_weight->SetPoint(iter, (charge[0] * pos[0].X() + charge[1] * pos[1].X())/(charge[0] + charge[1]), start_y[i]);
	//			iter++;
	//		}
	//	}//endl over points
	//
	//	TCanvas* c_weight = new TCanvas("", "", 1800, 900);
	//	c_weight->Divide(2);
	//	c_weight->cd(1);
	//	g_start->GetXaxis()->SetRangeUser(-75, 75);
	//	g_start->GetYaxis()->SetRangeUser(-75, 75);
	//
	//	g_start->Draw("AP*");
	//	draw_boxes(scints);
	//	c_weight->cd(2);
	//
	//	g_weight->GetXaxis()->SetRangeUser(-75, 75);
	//	g_weight->GetYaxis()->SetRangeUser(-75, 75);
	//
	//	g_weight->Draw("AP*");
	//	draw_boxes(scints);
	//########################################################################
} //end of code
