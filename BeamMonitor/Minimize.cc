#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "math.h"

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

ROOT::Math::Polar2DVector construct_SiPM(double r, double phi){
	return ROOT::Math::Polar2DVector(r, phi);
}

ROOT::Math::XYPoint get_point(){
	return ROOT::Math::XYPoint(2.0, -0.5);
}
double distance(double x, double y){
	return sqrt(pow(x, 2) + pow(y, 2));
}

//create a grid of points 
std::vector<ROOT::Math::XYPoint> get_grid(double radius, int max_enum){
	std::vector<ROOT::Math::XYPoint> helpvector;
	for(int i = 0; i <= max_enum; i++){
		double y = radius - i * (1.0 / max_enum) * 2 * radius;
		for (int j = 0; j <= max_enum; j++){
			double x = -radius + j * (1.0 / max_enum) * 2 * radius;
			if (distance(x, y) < radius){
				ROOT::Math::XYPoint helppoint(x, y);
				helpvector.push_back(helppoint);
			}
		}
	
	}
	return helpvector;

}


double Det_response(ROOT::Math::Polar2DVector SiPM_position, ROOT::Math::XYPoint position){
	ROOT::Math::XYPoint Det_position(SiPM_position.X(), SiPM_position.Y());
	auto rel_position = Det_position - position;
	// std::cout << "Relative position of point to detector is: " << rel_position << endl;
	// std::cout << "Distance of point from detector is " << distance(rel_position.X(), rel_position.Y()) << " cm" << endl;
	double r = distance(rel_position.X(), rel_position.Y());
	auto response = 1 / pow(r, 2);
	return response;
}

//################################################ MAIN FUNCTION ###########################################

void Minimize(double sigma = 0){
	int NDet = 16;
	std::vector<ROOT::Math::Polar2DVector> SiPMs;
	TString folder = Form("plots_blur_%.3f/", sigma);
	gSystem->Exec(Form("mkdir -p ./plots_blur_%.3f/", sigma));
	

	double radius = 3.0;
	for(int i = 0; i < NDet; i++){
		auto r = radius;
		auto phi = i * 2 * M_PI /NDet;
		SiPMs.push_back(construct_SiPM(r, phi));		//create SiPMs
	}					
	TRandom3 rand;

	std::vector<std::array<double, 2>> rec_coordinates;  //data arrays definition
	std::vector<double> chi2_values;
	


	int n_enum = 100;
	auto points = get_grid(radius, n_enum);		//Generating grid


	int n_crash = 0;
	double blur = 1;			//defining parameters used in minimization and signal blurring

	std::vector<std::vector<double>> signal_memory;
	for(int i = 0; i < points.size(); i++){
		auto point = points[i];
		std::vector<double> signals;
		for (int j = 0; j < SiPMs.size(); j++){
			auto signal = Det_response(SiPMs[j], point);
			signal = signal + sigma * rand.Gaus(0, blur); // blurring the signal
			signals.push_back(signal);
		}
		signal_memory.push_back(signals);
		auto chi2function = [&](const double* params) -> double{	//inline function declaration of chi2 function, that gets minimized 
			double chi2 = 0;
			double x = params[0];
			double y = params[1];
			if (sqrt(x * x + y * y) >= radius)  chi2+=10e6;
			for (int i = 0; i < NDet; i++){
				double dx = x - SiPMs[i].X();
				double dy = y - SiPMs[i].Y();			//defining the chi2 function	
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
		minimizer->SetVariableLimits(0, -3, 3);
		minimizer->SetVariableLimits(1, -3, 3);
				
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
		std::cout << "Primary point X: " << points[i].X() << " Y: " << points[i].Y() << std::endl;
		std::cout << "Best x: " << results[0] << "\t best y: " << results[1] << std::endl;    
		
		std::array<double, 2> help_array = {{results[0], results[1]}};
		rec_coordinates.push_back(help_array);

		
	}	//end of event loop

//##################################################### drawing plots #####################################################

	std::cout << "Minimization did not converge " << n_crash << " times out of " << rec_coordinates.size() << "points"  << std::endl;
	TGraph* Primary_grid = new TGraph(points.size());
	Primary_grid->SetName("Primary_grid");
	Primary_grid->SetTitle("Primary grid;x [cm];y [cm]");
	for(int i = 0; i < points.size(); i++){
		Primary_grid->SetPoint(i, points[i].X(), points[i].Y());
	}
	Primary_grid->SetMarkerStyle(20);
	
	TGraph* Secondary_grid = new TGraph(rec_coordinates.size());

	Secondary_grid->SetName("Secondary_grid");
	Secondary_grid->SetTitle("Grid after reconstruction through minimization;x [cm]; y[cm]");
	for(int i = 0; i < rec_coordinates.size(); i++){
		Secondary_grid->SetPoint(i, rec_coordinates[i][0], rec_coordinates[i][1]);
	}
	Secondary_grid->SetMarkerStyle(20);	
	
	TCanvas* c_grid = new TCanvas("", "", 1800, 900);
	c_grid->Divide(2);
	c_grid->cd(1);
	Primary_grid->Draw("AP");	
	c_grid->cd(2);
	Secondary_grid->Draw("AP");
	
	c_grid->Print(folder + "g_grid.png");
	delete c_grid;

	//################### arrow plots ################# 
	
	TCanvas* c_arrow = new TCanvas ("", "", 1800, 900);
	TGraph* arrow_graph = new TGraph(rec_coordinates.size());
	arrow_graph->GetXaxis()->SetRangeUser(-radius, radius); 
	arrow_graph->GetYaxis()->SetRangeUser(-radius, radius); 
	
	arrow_graph->SetName("arrow_graph");
	arrow_graph->SetTitle("point movement after reconstruction;x [cm]; y[cm]");
	
	arrow_graph->SetMarkerStyle(20);

	for(int i = 0; i < rec_coordinates.size();i++){
		arrow_graph->SetPoint (i, rec_coordinates[i][0], rec_coordinates[i][1]);
	}
	arrow_graph->Draw("ap");
	for (int i = 0 ; i < points.size(); i++){
		TArrow* arrow = new TArrow(points[i].X(), points[i].Y(), rec_coordinates[i][0], rec_coordinates[i][1], 0.005, "|>");
		arrow->SetLineColor(kBlue);
		arrow->Draw();
		}
	c_arrow->Print(folder + "g_arrow.png");
	delete c_arrow;
	
	//########################################### printing residuals ###################################################
	
	TGraph2D* residuals_x = new TGraph2D(rec_coordinates.size());
	residuals_x->SetName("Residuals_x");
	residuals_x->SetTitle("Residuals x;x[cm];y[cm];x_residual[cm]");
	TGraph2D* residuals_y = new TGraph2D(rec_coordinates.size());
	residuals_y->SetName("Residuals_y");
	residuals_y->SetTitle("Residuals y;x[cm];y[cm];y_residual[cm]");
	
	for(int i = 0; i < rec_coordinates.size(); i++){
		residuals_x->SetPoint(i, points[i].X(), points[i].Y(), points[i].X() - rec_coordinates[i][0]);
		residuals_y->SetPoint(i, points[i].X(), points[i].Y(), points[i].Y() - rec_coordinates[i][1]);
	}

	TCanvas* c_residual = new TCanvas("", "", 1800, 900);
	c_residual->Divide(2);
	c_residual->cd(1);
	residuals_x->Draw("TRI");
	c_residual->cd(2);
	residuals_y->Draw("TRI");
	
	c_residual->Print(folder + "residuals.png");
	delete c_residual;

	TH2D* h_x_residuals = new TH2D("h_x_residuals", "x Residuals;x [cm];y[cm];residual_x [cm]",
					n_enum + 1, -radius, radius, n_enum + 1, - radius, radius);
	TH2D* h_y_residuals = new TH2D("h_y_residuals", "y Residuals;x [cm];y[cm];residual_y [cm]",
					n_enum + 1, -radius, radius, n_enum + 1, - radius, radius);
	
	for(int i = 0; i < points.size(); i++){
		h_x_residuals->Fill(points[i].X(), points[i].Y(), rec_coordinates[i][0] - points[i].X()); 
		h_y_residuals->Fill(points[i].X(), points[i].Y(), rec_coordinates[i][1] - points[i].Y()); 
		}//end filling residuals histogram

	TCanvas* c_h_residuals = new TCanvas("", "", 1800, 900);
	c_h_residuals->Divide(2);
	c_h_residuals->cd(1);
	
	h_x_residuals->Draw("colz");
	
	c_h_residuals->cd(2);
	h_y_residuals->Draw("colz");
	
	c_h_residuals->Print(folder + "c_h_residuals.png");
	delete c_h_residuals;
	

	std::cout << "Start projections" << std::endl;
	TH1D* h_xx_projection = h_x_residuals->ProjectionX();
	h_xx_projection->SetName("h_xx_projection");
	h_xx_projection->SetTitle("x residual projection on x axis;x [cm]; x residual");	

	TH1D* h_yx_projection = h_x_residuals->ProjectionY();
	h_yx_projection->SetName("h_yx_projection");
	h_yx_projection->SetTitle("x residual projection on y axis;y [cm]; x residual");	

	TH1D* h_xy_projection = h_y_residuals->ProjectionX();
	h_xy_projection->SetName("h_xy_projection");
	h_xy_projection->SetTitle("y residual projection on x axis;x [cm]; y residual");	

	TH1D* h_yy_projection = h_y_residuals->ProjectionY();
	h_yy_projection->SetName("h_yy_projection");
	h_yy_projection->SetTitle("y residual projection on y axis;y [cm]; y residual");	

	std::cout << "Projections finished" << std::endl;
	TCanvas* c_h_projections = new TCanvas("", "", 1800, 900);
	
	c_h_projections->Divide(2, 2);
	c_h_projections->cd(1);
	
	h_xx_projection->Draw("hist");
	
	c_h_projections->cd(2);
	h_yx_projection->Draw("hist");

	c_h_projections->cd(3);
	
	h_xy_projection->Draw("hist");
	
	c_h_projections->cd(4);
	h_yy_projection->Draw("hist");
	
	c_h_projections->Print(folder + "c_h_projections.png");
	delete c_h_projections;
	TH1D* h_xx2_projection = new TH1D("h_xx2_projection", "x^2 residual projection on x axis;x [cm]; x^2 residual",
						n_enum + 1, -radius, radius);
	TH1D* h_yx2_projection = new TH1D("h_yx2_projection", "x^2 residual projection on y axis;y [cm]; x^2 residual",
						n_enum + 1, -radius, radius);
	TH1D* h_xy2_projection = new TH1D("h_xy2_projection", "y^2 residual projection on x axis;x [cm]; y^2 residual",
						n_enum + 1, -radius, radius);
	TH1D* h_yy2_projection = new TH1D("h_yy2_projection", "y^2 residual projection on y axis;y [cm]; y^2 residual",
						n_enum + 1, -radius, radius);
	for(int i = 0; i <= n_enum; i++){ //start for loop
		h_xx2_projection->SetBinContent(i+1, (h_xx_projection->GetBinContent(i)) * (h_xx_projection->GetBinContent(i)));
		h_yx2_projection->SetBinContent(i, h_yx_projection->GetBinContent(i) * h_yx_projection->GetBinContent(i));
		h_xy2_projection->SetBinContent(i, h_xy_projection->GetBinContent(i) * h_xy_projection->GetBinContent(i));
		h_yy2_projection->SetBinContent(i, h_yy_projection->GetBinContent(i) * h_yy_projection->GetBinContent(i));
	} //end for loop


	TCanvas* c_h2_projections = new TCanvas("", "", 1800, 900);
	
	c_h2_projections->Divide(2, 2);
	c_h2_projections->cd(1);
	
	h_xx2_projection->Draw("hist");
	
	c_h2_projections->cd(2);
	h_yx2_projection->Draw("hist");

	c_h2_projections->cd(3);
	
	h_xy2_projection->Draw("hist");
	
	c_h2_projections->cd(4);
	h_yy2_projection->Draw("hist");
	
	c_h2_projections->Print(folder + "c_h2_projections.png");
	delete c_h2_projections;

	//###

	TH2D* h_rx = new TH2D("h_rx", Form("occurence of x residuals on x with blurring = %f; y [cm]; residual x [cm]; counts", blur),
				n_enum + 1, - radius, radius, 5*n_enum + 1, - 0.3, 0.3);
	TH2D* h_ry= new TH2D("h_ry", Form("occurence of x residuals on y with blurring = %f; y [cm]; residual x [cm]; counts", blur),
				n_enum + 1, - radius, radius, 5*n_enum + 1, - 0.3, 0.3);
	for(int i = 0; i < points.size(); i++){
		double check = distance(points[i].X(), points[i].Y());
		if (check <= radius - 1){
		h_rx->Fill(points[i].X(), points[i].X() - rec_coordinates[i][0]);		
		h_ry->Fill(points[i].Y(), points[i].X() - rec_coordinates[i][0]);		
}} //end loop	
	TCanvas* c_res = new TCanvas("", "", 1800, 900);
	c_res->Divide(2);
	c_res->cd(1);
	
	h_rx->Draw("colz");
	
	c_res->cd(2);

	h_ry->Draw("colz");

	c_res->Print(folder + "h_res.png");	
	
	delete c_res;

	TCanvas* c_chi2 = new TCanvas("", "", 1800, 900);
	
	TH1D* h_chi2 = new TH1D("h_chi2", "Final #chi ^{2};event;#chi ^{2}", chi2_values.size(), 0, chi2_values.size()+1);
	for(int i = 0; i < chi2_values.size(); i++){
		h_chi2->SetBinContent(i+1, chi2_values[i]);	
	}//end loop
	h_chi2->Draw("hist");
	c_chi2->Print(folder + "h_chi2.png");
	delete c_chi2;


	TCanvas* c_dist = new TCanvas("", "", 1800, 900);
	TGraph2D* g_dist = new TGraph2D(points.size());
	for(int i = 0; i < points.size(); i++){
		
		double x_dif = points[i].X() - rec_coordinates[i][0];
		double y_dif = points[i].Y() - rec_coordinates[i][1];
		double deviation = sqrt(x_dif * x_dif + y_dif * y_dif);
		/*if (distance(points[i].X(), points[i].Y()) <= radius - 1) */g_dist->SetPoint(i, points[i].X(), points[i].Y(), deviation);
	}//end loop
	g_dist->SetName("g_dist");
	g_dist->SetTitle("distance of reconstructed point from its primary point;x[cm];y[cm];deviation[cm]"); 
	g_dist->Draw("TRI");
	c_dist->Print(folder + "g_dist.png");
	delete c_dist;

	TCanvas* c_xy_chi2 = new TCanvas("", "", 1800, 900);
	
	c_xy_chi2->Divide(2);

	c_xy_chi2->cd(1);
	
	g_dist->GetZaxis()->SetTitleOffset(2.0);
	g_dist->Draw("TRI");

	c_xy_chi2->cd(2);
	
	TH2D* h_xy_chi2 = new TH2D("h_xy_chi2", "chi2 values of points; x[cm];y[cm];#chi^{2}"
					, n_enum + 1, -radius, radius, n_enum+1, -radius, radius);
	for(int i = 0; i < points.size(); i++){
		
		if (distance(points[i].X(), points[i].Y()) <= radius - 1)h_xy_chi2->Fill(points[i].X(), points[i].Y(), chi2_values[i]);
		} //end loop
	gPad->SetLogz();
	gPad->SetRightMargin(0.15);
	gStyle->SetOptStat(0);
	h_xy_chi2->Draw("colz");
	c_xy_chi2->Print(folder + "c_xy_chi2.png");	
	delete c_xy_chi2;

	//################## COMPARISON OF WEIGHTED SUM ###########################

	TGraph* g_weighted = new TGraph(points.size());
	g_weighted->SetName("g_weighted");
	g_weighted->SetTitle("Position as sum of weighted charge;x[cm];y[cm]");
	for(int i = 0; i < points.size(); i++){
		std::array<double, 2> pos = {0, 0};
		double total_charge = 0.0;
		for(int j = 0; j < signal_memory[i].size(); j++){
			pos[0] += (signal_memory[i][j] * SiPMs[j].X());
			pos[1] += (signal_memory[i][j] * SiPMs[j].Y());
			total_charge += signal_memory[i][j];
		}
		pos[0] /= total_charge;
		pos[1] /= total_charge;
		g_weighted->SetPoint(i , pos[0], pos[1]);
		}//end for cycle
	TCanvas* c_comp = new TCanvas("", "", 1800, 900);
	c_comp->Divide(2);
	c_comp->cd(1);
	g_weighted->Draw("ap*");

	c_comp->cd(2);
	Secondary_grid->Draw("ap*");
	c_comp->Print(folder + "g_comp.png");
	delete c_comp;
	
	// ############## NEW MINIMIZATION? APPROXIMATE STARTING POSITION USING WEIGHTED CHARGE ###################
	
	std::vector<std::array<double, 2>> new_coordinates;
	for(int i = 0; i < points.size();i++){
		auto point = points[i];
		std::vector<double> signals = signal_memory[i];
		
		std::array<double, 2> start_position;
		double total_charge;
		for (int j = 0; j < signals.size(); j++) {
			start_position[0] += signals[j] * SiPMs[j].X();
			start_position[1] += signals[j] * SiPMs[j].Y();
			total_charge += signals[j];
		}
		start_position[0] /= total_charge;
		start_position[1] /= total_charge;

		auto chi2function = [&](const double* params) -> double{	//inline function declaration of chi2 function, that gets minimized 
			double chi2 = 0;
			double x = params[0];
			double y = params[1];
			if (sqrt(x * x + y * y) >= radius)  chi2+=10e6;
			for (int i = 0; i < NDet; i++){
				double dx = x - SiPMs[i].X();
				double dy = y - SiPMs[i].Y();			//defining the chi2 function	
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

		minimizer->SetVariable(0, "x", start_position[0], 0.001);
		minimizer->SetVariable(1, "y", start_position[1], 0.001);
		minimizer->SetVariableLimits(0, -3, 3);
		minimizer->SetVariableLimits(1, -3, 3);
				
		minimizer->SetMaxFunctionCalls(10000);		

		minimizer->Minimize();			//Minimizing the chi2
		if (minimizer->Status() != 0) {
			std::cerr << "Warning: Minimization did NOT converge! Status = " 
				<< minimizer->Status() << std::endl;
			n_crash++;
		}
		std::array<double, 2> help_arr = {minimizer->X()[0], minimizer->X()[1]};
		new_coordinates.push_back(help_arr);
	}//end loop over events
	
	TGraph2D* g_new = new TGraph2D(points.size());
	g_new->SetName("g_new");
	g_new->SetTitle("difference in end distances from primary point of reconstruction via primary approximation versus no primary approximation; x[cm]; y[cm]; difference[cm]");
	for(int i = 0; i < points.size(); i++){
		
		double x_dif = points[i].X() - rec_coordinates[i][0];
		double y_dif = points[i].Y() - rec_coordinates[i][1];
		double deviation1 = sqrt(x_dif * x_dif + y_dif * y_dif);
		x_dif = points[i].X() - new_coordinates[i][0];
		y_dif = points[i].Y() - new_coordinates[i][1];
		double deviation2 = sqrt(x_dif * x_dif + y_dif * y_dif);
		double final_deviation = deviation2 - deviation1;
		g_new->SetPoint(i, points[i].X(), points[i].Y(), final_deviation);
	}//end loop
	
	TCanvas* c_new = new TCanvas ("", "", 1800, 900);
	g_new->Draw("TRI");
	c_new->Print(folder + "c_new.png");
	delete c_new;
	
	TCanvas* c_dist2 = new TCanvas("", "", 1800, 900);
	TGraph2D* g_dist2 = new TGraph2D(points.size());
	for(int i = 0; i < points.size(); i++){
		double x_dif = points[i].X() - new_coordinates[i][0];
		double y_dif = points[i].Y() - new_coordinates[i][1];
		double deviation = sqrt(x_dif * x_dif + y_dif * y_dif);
		g_dist2->SetPoint(i, points[i].X(), points[i].Y(), deviation);
		}//end of loop	
	g_dist2->Draw("TRI");
	c_dist2->Print(folder + "c_dist_new.png");

	} //end of code
