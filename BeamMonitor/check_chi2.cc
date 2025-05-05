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
//################################################ MAIN FUNCTION ###########################################

void check_chi2(double sigma = 0){

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



	//TFile* file = new TFile("chi2.root", "RECREATE");
	TGraph* g_chi2 = new TGraph(201);
	for (int ip = 1; ip < 200; ip++){
		double point = SiPMs[4 + 8].X() + 2 * SiPMs[4].X() * ip / 200;
		double chi2 = chi2plot(point, points[4][30], SiPMs[12], SiPMs[4]);
		g_chi2->SetPoint(ip, point, chi2);
		cout << "Tested point has X = " << point << endl;
		cout << "chi squared for this point is " << chi2 << endl;
	}
	cout << "Primary point is: (" << points[4][30].X() << ", " << points[4][30].Y() << ")" << endl;
	TCanvas* c_check = new TCanvas("", "", 1800, 900);
	g_chi2->Draw("AL");
	long npoints = 0;
	int n_crash = 0;
	double blur = 1;			//defining parameters used in minimization and signal blurring

	std::vector<std::vector<double>> signal_memory;
	for(int iscint = 0; iscint < points.size() - 1; iscint++){
		for (int ievent = 0; ievent < points[iscint].size(); ievent++){
			npoints++;
			auto point = points[iscint][ievent];
			std::vector<double> signals;


			signals.push_back(Det_response(SiPMs[iscint], point));
			signals.push_back(Det_response(SiPMs[iscint+8], point));
			signal_memory.push_back(signals);


		}
	}	//end of event loop


} //end of code
