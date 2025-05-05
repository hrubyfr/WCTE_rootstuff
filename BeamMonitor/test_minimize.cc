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
#include <TH2.h>
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
#include <ios>
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
	double lambda = 1200.0;
	double atten = exp(-abs(rel_position)/lambda);
	auto response = atten / pow(r, 2);
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

double GetMinimum_chi2(const int niter, const double tolerance, const double signals[2], const XYVector SiPM1, const XYVector SiPM2){
	int iter = 0;
	double x = (SiPM1.X() * signals[0] + SiPM2.X() * signals[1])/(signals[0] + signals[1]); //Start point guess from weighted charge		
	double step_size = 10.0;
	double chi2 = 0;
	double chi2_history = 0;
	auto evaluate = [&x, SiPM1, SiPM2, signals, &chi2](){
		double dx1 = x - SiPM1.X();
		double dx2 = x - SiPM2.X();

		double theory_sig1 = 1.0 / (dx1 * dx1);
		double theory_sig2 = 1.0 / (dx2 * dx2);

		double residual1, residual2;					//Get chi2 at start point
		residual1 = signals[0] - theory_sig1;
		residual2 = signals[1] - theory_sig2;

		chi2 = residual1 * residual1 + residual2 * residual2;
		return chi2;
	};
	evaluate();
	int i_step = 0;
	int direction;
	if (signals[0] > signals[1]) direction = -1;
	else direction = 1;
	while (chi2 > tolerance){
		x += step_size * direction;
		chi2_history = chi2;
		chi2 = evaluate();


		if (chi2 > chi2_history) {
			direction *=-1;
			step_size /= 10.0;
		}
		if (abs(x) > abs(SiPM1.X())) {
			std::cerr << "ERROR: X OUT OF BOUNDS" << endl;
			return -1;
		}
		if(iter == niter) break;
		iter++;
	}
	return x;
}

//################################################ MAIN FUNCTION ###########################################

void test_minimize(double sigma = 0){

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
	double blur = 0.001;			//defining parameters used in minimization and signal blurring
	
	std::vector<std::vector<double>> signal_memory;
	for(int iscint = 0; iscint < points.size() - 1; iscint++){
		for (int ievent = 0; ievent < points[iscint].size(); ievent++){
			npoints++;
			auto point = points[iscint][ievent];
			std::vector<double> signals;
			double signal1 = Det_response(SiPMs[iscint], point);
			double signal2 = Det_response(SiPMs[iscint+8], point);
		//	signal1 += blur * rand.Gaus(0, 1); 
		//	signal1 += blur * rand.Gaus(0, 1); 
			signals.push_back(signal1);
			signals.push_back(signal2);
			signal_memory.push_back(signals);

			cout << endl;
			cout << "#########################################" << endl;
			cout << "Start point: (" << points[iscint][ievent].X() << ", " << points[iscint][ievent].Y() << ")" << endl;
			double arr[2] = {signals[0], signals[1]};
			double tolerance = 10e-12;
			int max_iter = 1000;
			double x_rec = GetMinimum_chi2(max_iter, tolerance, arr, SiPMs[iscint], SiPMs[iscint+8]);
			
			array<double, 2> help_arr = {x_rec, start_y[iscint]};
			rec_coordinates[iscint].push_back(help_arr);

			cout << "End point: (" << x_rec << ", " << start_y[iscint] << ")" << endl; 
			npoints++;
		} //end loop over event in scintillator
	}	//end loop over scintillators
	//##############################################
	vector<TBox*> scints;
	for(int i = 0; i < scint_dimensions.size(); i++){
		TBox* help_box = new TBox(-scint_dimensions[i][0]/2, start_y[i] - scint_dimensions[i][1]/2, scint_dimensions[i][0]/2, start_y[i] + scint_dimensions[i][1]/2) ;
		scints.push_back(help_box);
	}//end create scintillator tiles
	 //
	 //#################################################
	
	TCanvas* can = new TCanvas("", "", 1800, 900);

	can->Divide(2);
	can->cd(1);
	TGraph* g_start = new TGraph(npoints);
	TGraph* g_rec = new TGraph(npoints);
	int ip = 0;
	for (int i = 0; i < SiPMs.size() / 2; i++){
		for(int j = 0; j < points[i].size(); j++){
		g_rec->SetPoint(ip, rec_coordinates[i][j][0], rec_coordinates[i][j][1]);
		g_start->SetPoint(ip, points[i][j].X(), points[i][j].Y());
		ip++;
		}
	}

	g_start->Draw("AP*");
	draw_boxes(scints);
	can->cd(2);
	g_rec->Draw("AP*");
	draw_boxes(scints);

	//############################################################################
	TH2D* h_hits_init = new TH2D("h_hist_init", "Histogram of hits of primary points; x[mm]; y[mm]", 41, -scint_dimensions[4][0]/2, scint_dimensions[4][0]/2,
				8, start_y[0] - scint_dimensions[0][1]/2, start_y[7] + scint_dimensions[0][1]/2);
	TH2D* h_hits_final = new TH2D("h_hist_final", "Histogram of hits of reconstructed points+ x[mm]; y[mm]", 41, -scint_dimensions[4][0]/2, scint_dimensions[4][0]/2,
				8, start_y[0] - scint_dimensions[0][1]/2, start_y[7] + scint_dimensions[0][1]/2);
	for (int i = 0; i < SiPMs.size() / 2; i++){
		for(int j = 0; j < points[i].size(); j++){
		h_hits_final->Fill(rec_coordinates[i][j][0], rec_coordinates[i][j][1]);
		h_hits_init->Fill(points[i][j].X(), points[i][j].Y());
		}
	}
	
	TCanvas* c_hist = new TCanvas("", "", 1800, 900);
	c_hist->Divide(2);
	c_hist->cd(1);
	h_hits_init->Draw("colz");
	draw_boxes(scints);
	c_hist->cd(2);
	h_hits_final->Draw("colz");
	draw_boxes(scints);
} //end of code
