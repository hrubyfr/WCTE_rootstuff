#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "math.h"

#include "TFile.h"
#include "TH2D.h"
#include "Math/Vector2D.h"
#include "Math/Point2D.h"
#include <TGraph.h>
#include <TSystem.h>

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

void Minimize(){
	int NDet = 16;
	std::vector<ROOT::Math::Polar2DVector> SiPMs;

	gSystem->Exec("mkdir -p ./plots/");

	double radius = 3.0;
	for(int i = 0; i < NDet; i++){
		auto r = radius;
		auto phi = i * 2 * M_PI /NDet;
		SiPMs.push_back(construct_SiPM(r, phi));		//create SiPMs
	}					


	int n_enum = 100;
	auto points = get_grid(radius, n_enum);

	std::vector<std::array<double, 2>> rec_coordinates;
	int n_crash = 0;
	for(int i = 0; i < points.size(); i++){
		auto point = points[i];
		std::vector<double> signals;
		for (int j = 0; j < SiPMs.size(); j++){
			auto signal = Det_response(SiPMs[j], point);
			signals.push_back(signal);
		}
		auto chi2function = [&](const double* params) -> double{
			double chi2 = 0;
			double x = params[0];
			double y = params[1];
			if (sqrt(x * x + y * y) > radius)  chi2+=1e6;
			for (int i = 0; i < NDet; i++){
				double dx = x - SiPMs[i].X();
				double dy = y - SiPMs[i].Y();
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

		minimizer->SetVariable(0, "x", 0.0, 0.01);
		minimizer->SetVariable(1, "y", 0.0, 0.01);
		minimizer->SetVariableLimits(0, -3, 3);
		minimizer->SetVariableLimits(1, -3, 3);
				
		minimizer->SetMaxFunctionCalls(5000);		

		minimizer->Minimize();
		if (minimizer->Status() != 0) {
			std::cerr << "Warning: Minimization did NOT converge! Status = " 
				<< minimizer->Status() << std::endl;
			n_crash++;
		}

		const double* results = minimizer->X();
		std::cout << "Primary point X: " << points[i].X() << " Y: " << points[i].Y() << endl;
		std::cout << "Best x: " << results[0] << "\t best y: " << results[1] << endl;    
		
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
	Secondary_grid->SetTitle("Grid after reconstruction;x [cm]; y[cm]");
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
	
	c_grid->Print("plots/g_grid.png");

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
	
	c_residual->Print("plots/residuals.png");
	

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
	
	c_h_residuals->Print("plots/c_h_residuals.png");
	

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
	
	c_h_projections->Print("plots/c_h_projections.png");
	
	TH1D* h_xx2_projection = new TH1D("h_xx2_projection", "x^2 residual projection on x axis;x [cm]; x^2 residual",
						n_enum + 1, -radius, radius);
	TH1D* h_yx2_projection = new TH1D("h_yx2_projection", "x^2 residual projection on y axis;y [cm]; x^2 residual",
						n_enum + 1, -radius, radius);
	TH1D* h_xy2_projection = new TH1D("h_xy2_projection", "y^2 residual projection on x axis;x [cm]; y^2 residual",
						n_enum + 1, -radius, radius);
	TH1D* h_yy2_projection = new TH1D("h_yy2_projection", "y^2 residual projection on y axis;y [cm]; y^2 residual",
						n_enum + 1, -radius, radius);
	for(int i = 0; i <= n_enum; i++){ //start for loop
		std::cout << "########### debug ########### bin_content = " << h_xx_projection->GetBinContent(i) * h_xx_projection->GetBinContent(i)<< std::endl;
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
	
	c_h2_projections->Print("plots/c_h2_projections.png");
	


} //end of code
