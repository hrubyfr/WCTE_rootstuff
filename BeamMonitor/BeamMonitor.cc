#include "Math/Vector2D.h"
#include "Math/Point2D.h"
#include <Math/Point2Dfwd.h>
#include <Math/Vector2Dfwd.h>
#include <TAttMarker.h>
#include <TMath.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "TRandom3.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TMarker.h"
#include <cstdio>
#include "TFile.h"
#include "TGraph.h"
#include "TArrow.h"
#include "TStyle.h"

#include "TSystem.h"
#include "TText.h"
#include "TF1.h"

//___________________________HELP______FUNCTIONS_______________________________________________________

double distance(double x, double y){
    return sqrt(x*x + y*y);
}

TString double_to_string(double number){
    char buffer[50];
    snprintf(buffer, sizeof(buffer), "%.2f", number);
    return TString(buffer);
}
//######################################################################################################

//____________________________MC____POINT____GENERATION_____________________________________________________

ROOT::Math::XYPoint gen_point(double radius){
	TRandom3 rand;
	ROOT::Math::XYPoint point;
	
	point.SetX(rand.Uniform(-radius, radius));
	point.SetY(rand.Uniform(-radius, radius));
	while (distance(point.X(), point.Y()) > radius/2){
	point.SetX(rand.Uniform(-radius, radius));
	point.SetY(rand.Uniform(-radius, radius));
	}
	return point;
}

std::vector<ROOT::Math::XYPoint> Generate_Points(double x, double y, int point_number, double radius){

	//	GENERATES A GAUSSIAN SHAPE (MAYBE) TO TEST FOR MEAN POSITION BEFORE AND AFTER RECONSTRUCTION

    std::vector<ROOT::Math::XYPoint> points;
    TRandom3 random;
   

    while(true){
        
        // std::cout << "check x = \t" << x << "\t y = " << y << endl;
        if(distance(x, y) < radius/2){break;}
        
    }
    
    for (int i = 0; i < point_number; ++i){
        ROOT::Math::XYPoint point;
        double random_x = random.Gaus(x, 0.05);
        double random_y = random.Gaus(y, 0.05);
        if (distance(random_x,random_y) >= radius){i--;}
        else{
            point.SetX(random_x);
            point.SetY(random_y);
            points.push_back(point);
        }
    }
    return points;
}

std::vector<ROOT::Math::XYPoint> Generate_Shape(int point_number, double radius){ 

	// GENERATES A SEMICIRCLE + A RECTANGLE SHAPE FOR RECONSTRUCTION TESTING

    std::vector<ROOT::Math::XYPoint> points;
    TRandom3 random;
    double x, y;

    while(true){
        x = random.Uniform(-0.25, 0.25);
        y = random.Uniform(-0.25, 0.25);
        // std::cout << "check x = \t" << x << "\t y = " << y << endl;
        if(distance(x, y) < radius/4){break;}
    }
    
    for (int i = 0; i < point_number/2; ++i){
        ROOT::Math::XYPoint point;
        double random_x = random.Uniform(x - 0.1, x +0.1);
        double random_y = random.Uniform(y - 0.05, y + 0.25);
        if (distance(random_x,random_y) >= 0.5 || distance(random_x, random_y) <= 0.2){i--;}
        else{
            point.SetX(random_x);
            point.SetY(random_y);
            points.push_back(point);
        }
    }

    x = x * 0.92388 - y * 0.38268;
    y = x * 0.38268 + y * 0.92388;

    for (int i = 0; i < point_number/2; ++i){
        ROOT::Math::XYPoint point;
        double random_x = random.Uniform(x - 0.1, x +0.1);
        double random_y = random.Uniform(y - 0.25, y + 0.25);
        if (distance(random_x,random_y) >= 0.25|| distance(random_x, random_y) <= 0.17){i--;}
        else{
            point.SetX(random_x);
            point.SetY(random_y);
            points.push_back(point);
        }
    }
    return points;
}

std::vector<ROOT::Math::XYPoint> Generate_Grid(double radius){

    std::vector<ROOT::Math::XYPoint> points;

    for (int i = 0; i <= 20; ++i)
    {
        ROOT::Math::XYPoint point;
        double y = i * 0.1 * radius - 1*radius;
        for (int j = 0; j <= 20; ++j)
        {
            double x  = j * 0.1 * radius - 1*radius;
            point.SetX(x);
            point.SetY(y);
            if (distance(x, y) <= radius){points.push_back(point);}
            // std::cout<<"point " << i << " " << j << point << endl;
        }
    }
    return points;
}

std::vector<ROOT::Math::XYPoint> Generate_Grid_sim(double radius){

    std::vector<ROOT::Math::XYPoint> points;

    for (int i = 0; i <= 100; ++i)
    {
        ROOT::Math::XYPoint point;
        double y = i * 0.02 * radius - 1*radius;
        for (int j = 0; j <= 100; ++j)
        {
            double x  = j * 0.02 * radius - 1*radius;
            point.SetX(x);
            point.SetY(y);
            if (distance(x, y) <= radius){points.push_back(point);
            // std::cout<<"point X: " << x << " Y: " << y << endl;
            }
            
        }
    }
    return points;
}
//##########################################################################################################

//__________________________________GENERATION_____OF____SIGNAL_____IN____SIPMS_____________________________

double generate_signal_str(ROOT::Math::XYPoint point, ROOT::Math::Polar2DVector det_position, double lambda){

    double max_sig_strength = 1.0;
    auto distance_vector = point - det_position;
    double Sig_strength = max_sig_strength * exp((-1/lambda) * distance(distance_vector.X(), distance_vector.Y()));
    // std::cout << Sig_strength << endl;
    return Sig_strength;
}
double generate_signal_str_NoAtt(ROOT::Math::XYPoint point, ROOT::Math::Polar2DVector det_position){
    double max_sig_strength = 1.0;


    auto distance_vector = point - det_position;
    double Sig_strength = max_sig_strength * exp(-1 * distance(distance_vector.X(), distance_vector.Y()));
    // std::cout << Sig_strength << endl;
    return Sig_strength;
}
double generate_signal_onDistance(ROOT::Math::XYPoint point, ROOT::Math::Polar2DVector det_position){

    double max_sig_strength = 1.0;
    auto distance_vector = point - det_position;
    double Sig_strength = max_sig_strength * (1/pow(2,distance(distance_vector.X(), distance_vector.Y())));
    // std::cout << Sig_strength << endl;
    return Sig_strength;
}
//#####################################################################################################

//_______________________________________DIFFERENT____HISTOGRAMS_______________________________________
void draw_OP(std::vector<ROOT::Math::XYPoint> points, double radius){
    auto canvas = new TCanvas("Original points", "original points", 800, 800);

    canvas->Clear();
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    auto points_histogram = new TH2D("", "", 1000, -radius, radius, 1000, -radius, radius);
    for (int i = 0; i < points.size(); i++){
        points_histogram->Fill(points[i].X(), points[i].Y());
    }
    points_histogram->SetTitle("Original MC points");
    auto circle = new TEllipse(0, 0, radius, radius);
    circle->SetFillStyle(0);
    points_histogram->Draw("colz");
    circle->Draw("same");
}

void draw_arrows(std::vector<ROOT::Math::XYPoint> points, std::vector<ROOT::Math::XYPoint> reconstructed_points, TCanvas* canvas, double lambda, TString filename, double radius){
    canvas->Clear();
    TString lambda_truncated = double_to_string(lambda);
    TString canvas_title = "Points movement after reconstruction " + lambda_truncated;
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);


    auto points_histogram = new TH2D("", "", 100, -radius, radius, 100, -radius, radius);
    for (int i = 0; i < points.size(); i++){
        points_histogram->Fill(points[i].X(), points[i].Y());
    }
    points_histogram->SetTitle("Original MC points");

    auto reconstructed_histogram = new TH2D("", "", 100, -radius, radius, 100, -radius, radius);
    for (int i = 0; i < points.size(); i++){
        reconstructed_histogram->Fill(reconstructed_points[i].X(), reconstructed_points[i].Y());
    }
    reconstructed_histogram->SetTitle("Reconstruction #lambda = " + lambda_truncated);

    auto circle = new TEllipse(0, 0, radius, radius);
    circle->SetFillStyle(0);

    
    
    canvas->SetTitle(canvas_title);
    canvas->Divide(2);
    canvas->cd(1);
    points_histogram->Draw("colz");
    circle->Draw("same");
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->Update();
    canvas->cd(2);
    

    reconstructed_histogram->Draw("colz");
    for(int i = 0; i < points.size(); i++){
        auto arrow = new TArrow(points[i].X(), points[i].Y(), reconstructed_points[i].X(), reconstructed_points[i].Y(), 0.005, "|>");
        arrow->Draw();
    }
    circle->Draw("same");
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->Update();
    canvas->Print(filename);
    TString pngname = "Beam Monitor Grid lambda " + lambda_truncated + ".png";
    canvas->SaveAs(pngname);
}


void draw_histogram(std::vector<ROOT::Math::XYPoint> points, std::vector<ROOT::Math::XYPoint> reconstructed_points, std::vector<ROOT::Math::Polar2DVector> SiPMTVector, TString Title, double lambda, double radius){
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);

    auto points_histogram = new TH2D("", "", 100, -radius, radius, 100, -radius, radius);

    for (int i = 0; i < points.size(); i++){
        points_histogram->Fill(points[i].X(), points[i].Y());
    }
    points_histogram->SetTitle("Original MC points");

    auto reconstructed_histogram = new TH2D("", "", 100, -1, 1, 100, -1, 1);

    for (int i = 0; i < points.size(); i++){
        reconstructed_histogram->Fill(reconstructed_points[i].X(), reconstructed_points[i].Y());
    }
    reconstructed_histogram->SetTitle("Reconstruction, " + Title);

    auto circle = new TEllipse(0, 0, radius, radius);
    circle->SetFillStyle(0);

    auto canvas = new TCanvas();
    canvas->SetTitle(Title);
    canvas->Divide(2);

    canvas->cd(1);
    points_histogram->Draw("colz");
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->Update();

    canvas->cd(2);
    reconstructed_histogram->Draw("colz");
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->Update();

    for(int i = 1; i < 3; ++i){
        canvas->cd(i);
        circle->Draw("same");
        for (int j = 0; j < SiPMTVector.size(); j++){
            auto pmt = new TMarker(SiPMTVector[j].X(), SiPMTVector[j].Y(), kFullCircle);
            pmt->SetMarkerColor(kRed);
            pmt->Draw("same");
    }
        auto center = new TMarker(0, 0, kFullCircle);
        center->SetMarkerColor(kBlack);
        center->Draw("same");

    }
    TString filename = "Att. length = " + std::to_string(lambda) +".pdf";
    TString filenamejpg = "Att. length = " + std::to_string(lambda) +".jpg";

    if (lambda == 0){canvas->SaveAs("NoAtt.pdf");}
    else {canvas->SaveAs(filename);}
    if (lambda == 0){canvas->SaveAs("NoAtt.jpg");}
    else {canvas->SaveAs(filenamejpg);}

    
}

// _______________ DRAW TH2D OF TOF0 AND TOF1 RESPONSE VIA GRID POINTS__________________________
void FillCharge(TH2D* histogram, double SiPMCharge, ROOT::Math::XYPoint point){
    histogram->Fill(point.X(), point.Y(), SiPMCharge);
}

void TOFsim_hist(double radius, std::vector<ROOT::Math::XYPoint> points, std::vector<double> SiPMCharge, TString title){
    auto hist = new TH2D(title, title, 101, -radius, radius, 101, -radius, radius);
    hist->GetXaxis()->SetTitle("x [cm]");
    hist->GetYaxis()->SetTitle("y [cm]");
    hist->GetZaxis()->SetTitle("arbitrary charge");
    for (int i = 0; i < points.size(); i++){
        FillCharge(hist, SiPMCharge[i], points[i]);
    }
    hist->Draw("colz");
    hist->Write();
}

//____________________GET___MEAN___RANGE___AFTER___RECONSTRUCTION________________________________

double Get_Mean_Range(std::vector<ROOT::Math::XYPoint> reconstructed_points, double radius){
    auto reconstructed_histogram = new TH2D("", "", 100, -radius, radius, 100, -radius, radius);

    for (int i = 0; i < reconstructed_points.size(); i++){
        reconstructed_histogram->Fill(reconstructed_points[i].X(), reconstructed_points[i].Y());
    }
    double mean_x = reconstructed_histogram->GetMean(1);
    double mean_y = reconstructed_histogram->GetMean(2);
    double mean_distance = distance(mean_x, mean_y);
    return mean_distance;
}


//#############################################################################################

//____________________________CREATE___TH1D___RANGE___ON___LAMBDA_____________________________
// void graph_RLambda(double mean_range, double lambda, TCanvas* canvas){
//     canvas->Clear();
//     auto graph = new TGraph(50, mean_range, lambda);
//     graph->Draw();
//     canvas->SaveAs("graph_RLambda.png");
//     canvas->SaveAs("graph_RLambda.pdf");
// }
//_________________________________DETECTOR____POSITIONS_____________________________________

std::vector<ROOT::Math::Polar2DVector> Generate_8SiPMVector(double radius){
    std::vector<ROOT::Math::Polar2DVector> SiPMTVector;
    for (int i = 0; i < 8; i++){
        int pm_number = i;
        double angle = i * TMath::Pi()/4;
        ROOT::Math::Polar2DVector helpVector(radius, angle);
        SiPMTVector.push_back(helpVector);
        // std::cout << "______________SiPMT " << pmt_number << "______________" << endl;
        // std::cout << "SiPMTVector check: " << SiPMTVector[pmt_number] << endl;
        // std::cout << "SiPMTVector coordinates:\t x:" << SiPMTVector[pmt_number].X() << "\t y:" << SiPMTVector[pmt_number].Y() << endl;
    }
    return SiPMTVector;
}
std::vector<ROOT::Math::Polar2DVector> Generate_SiPMTVector(double radius){
    std::vector<ROOT::Math::Polar2DVector> SiPMTVector;
    for (int i = 0; i < 16; i++){
        int pm_number = i;
        double angle = i * TMath::Pi()/8;
        ROOT::Math::Polar2DVector helpVector(radius, angle);
        SiPMTVector.push_back(helpVector);
        // std::cout << "______________SiPMT " << pmt_number << "______________" << endl;
        // std::cout << "SiPMTVector check: " << SiPMTVector[pmt_number] << endl;
        // std::cout << "SiPMTVector coordinates:\t x:" << SiPMTVector[pmt_number].X() << "\t y:" << SiPMTVector[pmt_number].Y() << endl;
    }
    return SiPMTVector;
}
//###########################################################################################

//________________________________________POINT____RECONSTRUCTION___________________________________

ROOT::Math::XYPoint Reconstruct_Point(std::vector<double> PMT_signal, std::vector<ROOT::Math::Polar2DVector> SiPMTVector){
    ROOT::Math::XYPoint reconstructed_point(0,0);
    
    double total_signal;
    
        for (int j = 0; j < SiPMTVector.size(); ++j)
        {
            reconstructed_point += PMT_signal[j] * SiPMTVector[j];
            total_signal += PMT_signal[j];
        }
        reconstructed_point/=total_signal;
        // std::cout << "reconstructed point position:" << reconstructed_points[i] << endl;
        

    
    return reconstructed_point;
}
std::vector<ROOT::Math::XYPoint> Reconstruct_Points(std::vector<std::vector<double>> PMT_signal, std::vector<ROOT::Math::Polar2DVector> SiPMTVector, std::vector<ROOT::Math::XYPoint> points){
    std::vector<ROOT::Math::XYPoint> reconstructed_points;
    ROOT::Math::XYPoint init_point;
    double total_signal;
    for(int i = 0; i < points.size(); i++){
        reconstructed_points.push_back(init_point);
        for (int j = 0; j < SiPMTVector.size(); ++j)
        {
            reconstructed_points[i] += PMT_signal[i][j] * SiPMTVector[j];
            /* total_signal += PMT_signal[i][j];
            reconstructed_points[i]/=total_signal; */
        }
        // std::cout << "reconstructed point position:" << reconstructed_points[i] << endl;
        

    }
    return reconstructed_points;
}

//################################################################################################

//_________________________FILL___BEAM___FROM___CHARGE___IN___SIPM_______________________________


void fillBeamPosition(TH2D* histogram, std::vector<double> SiPM_charge, std::vector<ROOT::Math::Polar2DVector> SiPMVector){
    ROOT::Math::XYPoint BeamPosition;
    double Total_charge = 0;
    for(int i = 0; i < SiPMVector.size(); i++){
        BeamPosition += SiPM_charge[i] * SiPMVector[i];
        Total_charge += SiPM_charge[i];
    }
    BeamPosition /= Total_charge;
    histogram->Fill(BeamPosition.X(), BeamPosition.Y());
}

void DrawBeamHistogram(std::vector<ROOT::Math::Polar2DVector> SiPMVector, std::vector<std::vector<double>> SiPM_charges, double Scint_radius){
    TString can_title = "Beam histogram no attenuation";
    auto Beam_canvas = new TCanvas(can_title, can_title, 800, 800);

    TString Hist_title = "Histogram of approximated beam position from SiPM integrated charge";
    TString Hist_xAxis = "x [cm]";
    TString Hist_yAxis = "y [cm]";

    auto Beam_histogram = new TH2D("", "", 100, -Scint_radius, Scint_radius, 100, -Scint_radius, Scint_radius);
    Beam_histogram->SetTitle(Hist_title);
    Beam_histogram->GetXaxis()->SetTitle(Hist_xAxis);
    Beam_histogram->GetYaxis()->SetTitle(Hist_yAxis);

    

    for (int i = 0; i < SiPM_charges.size(); i++){
        fillBeamPosition(Beam_histogram, SiPM_charges[i], SiPMVector);
    }
    Beam_histogram->Draw("colz");

    for(int i = 0; i < SiPMVector.size(); i++){
        auto SiPM_position = new TMarker(SiPMVector[i].X(), SiPMVector[i].Y(), kFullCircle);
        SiPM_position->SetMarkerColor(kBlue);
        SiPM_position->Draw("same");
    }

    auto ellipse = new TEllipse(0, 0, Scint_radius, Scint_radius);
    ellipse->SetFillStyle(0);
    ellipse->SetLineColor(kBlack);
    ellipse->Draw("same");
}

//____________________________________MAIN____FUNCTION__________________________________________

void BeamMonitor() {

    //### CONSTANTS ###
    //
	TString plots_folder = "beam_monitor_simulation";
	
	if (gSystem->AccessPathName(plots_folder)) gSystem->Exec("mkdir -p " + plots_folder);
	gSystem->cd(plots_folder);


	double radius = 1.0;
	int point_number = 1000;


	auto shape_points = Generate_Shape(point_number, radius);


	//### GRID PART ###

	const char* filename = "BeamMonitorGrid.pdf";

	std::vector<std::vector<double>> PMT_signal_NoAtt;

	auto SiPMTVector = Generate_SiPMTVector(radius);
	auto SiPMVector_8 = Generate_8SiPMVector(radius);

	// SHAPE RECONSTRUCTION WOHOOOOOOOOOOOOO
	//
	TCanvas* c_shape_reconstruction = new TCanvas("c_shape_reconstruction", "c_shape_reconstruction", 1800, 900);
	TH2D* h_shape_original = new TH2D("h_shape_original", "Original points", 100, - radius, radius, 100, -radius, radius);
	for (int i = 0; i < shape_points.size(); i++){
		h_shape_original->Fill(shape_points[i].X(), shape_points[i].Y());
	}
	c_shape_reconstruction->Divide(4, 2);
	auto ellipse_shape = new TEllipse(0, 0, radius, radius);
	ellipse_shape->SetFillStyle(0);
	ellipse_shape->SetLineColor(kBlack);
	for (int i = 0; i < 4; i++){
		int j = 2 * i;
		c_shape_reconstruction->cd(j + 1);
		h_shape_original->Draw("colz");
		ellipse_shape->Draw("same");
		for(int i = 0; i < SiPMTVector.size(); i++){
			auto SiPM_position = new TMarker(SiPMVector_8[i].X(), SiPMVector_8[i].Y(), kFullCircle);
			SiPM_position->SetMarkerColor(kBlue);
			SiPM_position->Draw("same");
		}
		c_shape_reconstruction->cd(j+2);
		double lambda;
		std::vector<std::vector<double>> SiPMs_signals;
		if (i == 0) lambda = 0.5;
		else lambda = i;

		for (int ipoint = 0; ipoint < shape_points.size(); ipoint++){
			std::vector<double> point_responses;
			for (int idet = 0; idet < SiPMVector_8.size(); idet++){
				double point_response = generate_signal_str(shape_points[ipoint], SiPMVector_8[idet], lambda);
				point_responses.push_back(point_response);
			}
			SiPMs_signals.push_back(point_responses);

		}
		auto reconstructed_points = Reconstruct_Points(SiPMs_signals, SiPMVector_8, shape_points);
		TH2D* h_shape_reconstructed = new TH2D(Form("h_shape_reconstructed_%i", i), Form("Reconstructed points with #lambda = %0.1f", lambda), 100, - radius, radius, 100, -radius, radius);

		for (int ipoint = 0; ipoint  < reconstructed_points.size(); ipoint++){
			h_shape_reconstructed->Fill(reconstructed_points[ipoint].X(), reconstructed_points[ipoint].Y());
		}
		h_shape_reconstructed->Draw("colz");
		ellipse_shape->Draw("same");
		for(int i = 0; i < SiPMVector_8.size(); i++){
			auto SiPM_position = new TMarker(SiPMVector_8[i].X(), SiPMVector_8[i].Y(), kFullCircle);
			SiPM_position->SetMarkerColor(kBlue);
			SiPM_position->Draw("same");
		}

	}
	c_shape_reconstruction->Print("h_shape_reconstruction.png");
	c_shape_reconstruction->Print("h_shape_reconstruction.pdf");
	delete c_shape_reconstruction;

	// #######################################################################
	//
	//
	// 	GRID RECONSTRUCTION BEFORE SIPM_NUMBER CHANGE
	
	auto grid_points = Generate_Grid(radius);
	
	TCanvas* c_grid_reconstruction = new TCanvas("c_grid_reconstruction", "c_grid_reconstruction", 1800, 900);
	TH2D* h_grid_original = new TH2D("h_grid_original", "Original points", 100, - radius, radius, 100, -radius, radius);
	double mean_distance;
	TGraph *g_distances = new TGraph();
	for (int i = 0; i < grid_points.size(); i++){
		h_grid_original->Fill(grid_points[i].X(), grid_points[i].Y());
	}
	c_grid_reconstruction->Divide(2);
	for (int i = 0; i < 10; i++){
		double lambda = 0.5 * (i+1);
		std::vector<double> distances;
		for (int j = 0; j < 2; j++){
			c_grid_reconstruction->cd(j+1);
			if (j == 0){
				h_grid_original->Draw("colz");
			}
			else{
				std::vector<std::vector<double>> SiPMs_signals;

				for (int ipoint = 0; ipoint < grid_points.size(); ipoint++){
					std::vector<double> point_responses;
					for (int idet = 0; idet < SiPMVector_8.size(); idet++){
						double point_response = generate_signal_str(grid_points[ipoint], SiPMVector_8[idet], lambda);
						point_responses.push_back(point_response);
					}
					SiPMs_signals.push_back(point_responses);

				}
				auto reconstructed_points = Reconstruct_Points(SiPMs_signals, SiPMVector_8, grid_points);
				TH2D* h_grid_reconstructed = new TH2D("h_grid_reconstructed", Form("Reconstructed points with #lambda = %0.1f", lambda), 100, - radius, radius, 100, -radius, radius);

				for (int ipoint = 0; ipoint  < reconstructed_points.size(); ipoint++){
					h_grid_reconstructed->Fill(reconstructed_points[ipoint].X(), reconstructed_points[ipoint].Y());
				}
				h_grid_reconstructed->Draw("colz");
				for (int iarrow = 0; iarrow < grid_points.size(); iarrow++){
					auto arrow = new TArrow(grid_points[iarrow].X(), grid_points[iarrow].Y(), reconstructed_points[iarrow].X(), reconstructed_points[iarrow].Y(), 0.005, "|>");
					arrow->Draw();
				distances.push_back(distance(grid_points[iarrow].X() - reconstructed_points[iarrow].X(), grid_points[iarrow].Y() - reconstructed_points[iarrow].Y()));
				}


			}

			ellipse_shape->Draw("same");
			for(int k = 0; k < SiPMVector_8.size(); k++){
				auto SiPM_position = new TMarker(SiPMVector_8[k].X(), SiPMVector_8[k].Y(), kFullCircle);
				SiPM_position->SetMarkerColor(kBlue);
				SiPM_position->Draw("same");
			}
			

		}
		mean_distance = 0;
		for (int ipoint = 0; ipoint < distances.size(); ipoint++){
			mean_distance+=distances[ipoint];
		}
		mean_distance/=distances.size();

		c_grid_reconstruction->Print(Form("h_grid_reconstruction_%0.1f.png", lambda));
		c_grid_reconstruction->Print(Form("h_grid_reconstruction_%0.1f.pdf", lambda));
		g_distances->AddPoint(lambda, mean_distance);
	}
	delete c_grid_reconstruction;

	TCanvas *c_mean_distance = new TCanvas("c_mean_distance", "c_mean_distance", 1800, 900);
	g_distances->SetTitle("Graph of mean distance of all reconstructed points from the original points;#lambda;mean distance");
	g_distances->Draw("AC*");
	c_mean_distance->Print("g_mean_distance.png");
	c_mean_distance->Print("g_mean_distance.pdf");

	// #####################################################################

	auto point = gen_point(radius);
	std::vector<double> responses(8);
	for (int i = 0; i < 8; i++){
		responses[i] = generate_signal_str_NoAtt(point, SiPMVector_8[i]);
	}
	auto rec_point = Reconstruct_Point(responses, SiPMVector_8);

	std::cout << "point position: " << point.X() << ", " << point.Y() << std::endl; 
	std::cout << "Reconstructed point position: " << rec_point.X() << ", " << rec_point.Y() << std::endl; 


	TCanvas* c_point = new TCanvas("c_point", "c_point", 1800, 900);
	c_point->Divide(2);
	double x[1] = {point.X()};
	double y[1] = {point.Y()};
	c_point->cd(1);
	gPad->SetLeftMargin(0.15);
	TGraph* g_point = new TGraph(1, x, y);

	g_point->SetTitle("Primary particle position");
	g_point->GetXaxis()->SetTitle("x[cm]");
	g_point->GetYaxis()->SetTitle("y[cm]");
	g_point->GetXaxis()->SetLimits(-radius, radius);
	g_point->SetMinimum(-radius);
	g_point->SetMaximum(radius);
	g_point->SetMarkerStyle(kFullCircle);
	g_point->SetMarkerColor(kRed);
	g_point->SetMarkerSize(1.5);

	g_point->Draw("AP");
	auto ellipse = new TEllipse(0, 0, radius, radius);
	ellipse->SetFillStyle(0);
	ellipse->SetLineColor(kBlack);
	ellipse->Draw("same");
	for(int i = 0; i < SiPMTVector.size(); i++){
		auto SiPM_position = new TMarker(SiPMVector_8[i].X(), SiPMVector_8[i].Y(), kFullCircle);
		SiPM_position->SetMarkerColor(kBlue);
		SiPM_position->Draw("same");
	}
	c_point->cd(2);

	gPad->SetLeftMargin(0.15);
	x[0] = rec_point.X();
	y[0] = rec_point.Y();
	TGraph* g_rec_point = new TGraph(1, x, y);
	g_rec_point->SetTitle("Reconstructed particle position");
	g_rec_point->GetXaxis()->SetTitle("x[cm]");
	g_rec_point->GetYaxis()->SetTitle("y[cm]");
	g_rec_point->GetXaxis()->SetLimits(-radius, radius);
	g_rec_point->SetMinimum(-radius);
	g_rec_point->SetMaximum(radius);
	g_rec_point->SetMarkerStyle(kFullCircle);
	g_rec_point->SetMarkerColor(kRed);
	g_rec_point->SetMarkerSize(1.5);

	g_rec_point->Draw("AP");
	ellipse->Draw("same");
	for(int i = 0; i < SiPMTVector.size(); i++){
		auto SiPM_position = new TMarker(SiPMVector_8[i].X(), SiPMVector_8[i].Y(), kFullCircle);
		SiPM_position->SetMarkerColor(kBlue);
		SiPM_position->Draw("same");
	}
	c_point->Print("plots/point_reconstruction.pdf");


	//#### GRID SIGNAL RESPONSE SIMULATION

	double TOF0_radius = 3.0;
	double TOF1_radius = 4.0;

	auto TOF0_SiPM = Generate_SiPMTVector(TOF0_radius);
	auto TOF1_SiPM = Generate_SiPMTVector(TOF1_radius);

	std::vector<ROOT::Math::XYPoint> sim_points_TOF0 = Generate_Grid_sim(TOF0_radius);
	std::vector<ROOT::Math::XYPoint> sim_points_TOF1 = Generate_Grid_sim(TOF1_radius);
	std::vector<std::vector<double>> TOF0_signal;
	std::vector<std::vector<double>> TOF1_signal;

	draw_OP(sim_points_TOF0, TOF0_radius);

	for(int i = 0; i < TOF0_SiPM.size(); i++){
		std::vector<double> v1;

		for (int j = 0; j < sim_points_TOF0.size(); j++)
		{
			v1.push_back(generate_signal_onDistance(sim_points_TOF0[j], TOF0_SiPM[i]));
		}
		TOF0_signal.push_back(v1);
	}
	for(int i = 0; i < TOF1_SiPM.size(); i++){
		std::vector<double> v1;

		for (int j = 0; j < sim_points_TOF1.size(); j++)
		{
			v1.push_back(generate_signal_onDistance(sim_points_TOF1[j], TOF1_SiPM[i]));
		}
		TOF1_signal.push_back(v1);
	}
	auto file = TFile::Open("output.root", "RECREATE");
	auto can1 = new TCanvas("TOF0_sim", "TOF0_sim", 1980, 800);
	can1->Divide(4, 4);
	for(int i = 0; i < TOF0_signal.size(); i++){
		TString hist_title = "SiPM_TOF0_" + std::to_string(i);
		can1->cd(i+1);
		file->cd();
		TOFsim_hist(TOF0_radius, sim_points_TOF0, TOF0_signal[i], hist_title);

	}
	can1->Update();
	auto can2 = new TCanvas("TOF1_sim", "TOF1_sim", 1980, 800);
	can2->Divide(4, 4);
	for(int i = 0; i < TOF1_signal.size(); i++){
		TString hist_title = "SiPM_TOF1_" + std::to_string(i);
		can2->cd(i+1);
		TOFsim_hist(TOF1_radius, sim_points_TOF1, TOF1_signal[i], hist_title);
	}
	can2->Update();

	can1->SaveAs("TOF0 SiPM response.png");
	can2->SaveAs("TOF1 SiPM response.png");
	can1->SaveAs("TOF0 SiPM response.pdf");
	can2->SaveAs("TOF1 SiPM response.pdf");

	file->Close();

	// auto can3 = new TCanvas();
	// can3->Divide(2);
	// can3->cd(1);
	// double NDet = 16;
	// TH1D* Yhist[NDet];
	// TH1D* Xhist[NDet];
	// Yhist[0] = new TH1D("", "", 101, -TOF0_radius, TOF0_radius);
	// Xhist[0] = new TH1D("", "", 101, -TOF0_radius, TOF0_radius);
	// auto coordx = SiPMVector[0]->GetX();
	// auto coordy = SiPMVector[0]->GetY();
	// for(int i = 0; i < sim_points_TOF0.size(); i++){
	//     if(sim_points_TOF0[i].X() == TOF0_SiPM[i])
	// }

	// for (int i = 0; i < Yhist[0].GetNBins(); i++){
	//     Yhist[0]->SetBinContent(i, TOFsim_hist)
	// }


	//#########################################################

	std::vector<std::vector<double>> SiPM_charges;
	for (int i = 0; i < 10; i++){
		std::vector<double> SiPM_charge = {3, 2.5, 2, 1.5, 1, 0, 1, 1.5, 3, 2.5, 2, 1.5, 1, 0, 1, 1.5};
		SiPM_charges.push_back(SiPM_charge);
	}



	//### BEAM POSITION HISTOGRAM FROM INTEGRATED CHARGE ###
	DrawBeamHistogram(SiPMTVector, SiPM_charges, radius);


	// std::vector<ROOT::Math::XYPoint> points = Generate_Shape(point_number, radius);


	std::vector<ROOT::Math::XYPoint> points = Generate_Grid(radius);

	auto canvas = new TCanvas("", "", 1920, 860);
	canvas->Print(Form("%s[", filename));
	for (double lambda = 0.5; lambda <= 5.0; lambda += 0.5){
		std::vector<std::vector<double>> PMT_signal;
		for (int i = 0; i < points.size(); i++){
			std::vector<double> v1;
			// std::cout << points[i] << endl;
			for (int j = 0; j < SiPMTVector.size(); j++){
				// std::cout << "DEBUG 1" << endl;
				v1.push_back(generate_signal_str(points[i], SiPMTVector[j], lambda));
				// std::cout << "DEBUG 2" << endl;
			}
			PMT_signal.push_back(v1);
		}

		auto reconstructed_points_Att = Reconstruct_Points(PMT_signal, SiPMTVector, points);
		draw_arrows(points, reconstructed_points_Att, canvas, lambda, filename, radius);

	}
	canvas->Print(Form("%s]", filename));


	// for (int i = 0; i < points.size(); i++){
	//     std::vector<double> v1;
	//     // std::cout << points[i] << endl;
	//     for (int j = 0; j < SiPMTVector.size(); j++){
	//         v1.push_back(generate_signal_str(points[i], SiPMTVector[j], lambda));
	// }
	//     PMT_signal.push_back(v1);
	// }

	//### MEAN DISTANCE AFTER RECONSTRUCTION WITH DIFFERENT LAMBDA PART ###

	points = Generate_Points(radius/5, radius/5, point_number, radius);
	double mean_range[499];
	double lambdas[499];
	int j = 0;
	for (double lambda = 0.1; lambda < 50.0; lambda += 0.1){

		std::vector<std::vector<double>> PMT_signal;
		for (int i = 0; i < points.size(); i++){
			std::vector<double> v1;
			// std::cout << points[i] << endl;
			for (int j = 0; j < SiPMVector_8.size(); j++){
				v1.push_back(generate_signal_str(points[i], SiPMVector_8[j], lambda));
			}
			PMT_signal.push_back(v1);
		}

		auto reconstructed_points_Att = Reconstruct_Points(PMT_signal, SiPMVector_8, points);
		double R = Get_Mean_Range(reconstructed_points_Att, radius);
		// std::cout<<"mean range with lambda = " << lambda << " is R = " << R << endl;
		mean_range[j] = R;
		// std::cout<<"mean range array has: " << mean_range[j] <<endl; 

		lambdas[j] = lambda;
		// std::cout<<"lambdas array has: " << lambdas[j] <<endl; 
		j++;
	}
	// for (int i = 0; i < 50; i++)
	// {
	//     // std::cout << "mean_range array: " << mean_range[i] << endl;
	//     // std::cout << "lambdas array: " << lambdas[i] << endl;
	// }

	auto graph = new TGraph();
	for (int i = 0; i < 499; i++){
		graph->AddPoint(lambdas[i], mean_range[i]);
		std::cout << "LAMBDA: " << lambdas[i] << " RANGE: " << mean_range[i] << std::endl;
	}

	canvas->Clear();

	double original_mean = distance(radius/5, radius/5);
	auto line = new TLine(0, original_mean, 10, original_mean);
	line->SetLineColor(kBlue);
	line->SetLineWidth(2);

	auto text = new TText(7, original_mean, "Original points mean");
	text->SetTextColor(kBlue);
	graph->GetXaxis()->SetTitle("#lambda");
	graph->GetYaxis()->SetTitle("<R>");
	graph->SetTitle("Mean range on attenuation length #lambda");
	graph->GetXaxis()->SetRangeUser(0, 10);
	graph->Draw("AC*");
	line->Draw("same");
	text->Draw();




	canvas->SaveAs("RangeonLambda.png");
	canvas->SaveAs("RangeonLambda.pdf");

	// graph_RLambda(mean_range, lambdas, canvas);
	// for (int i = 0; i < points.size(); i++){
	//     std::vector<double> v1;
	//     // std::cout << points[i] << endl;

	//     for (int j = 0; j < SiPMTVector.size(); j++){
	//         v1.push_back(generate_signal_str_NoAtt(points[i], SiPMTVector[j]));

	// }


	//     PMT_signal_NoAtt.push_back(v1);


	// }





	// auto reconstructed_points_Att = Reconstruct_Points(PMT_signal, SiPMTVector, points);

	// auto reconstructed_points_NoAtt = Reconstruct_Points(PMT_signal_NoAtt, SiPMTVector, points);

	// std::cout << "DEBUG: original point: " << points[0] << " reconstructed point with attenuation: " << reconstructed_points_Att[0] << " reconstructed point with no attenuation " << reconstructed_points_NoAtt[0];


	//     //CREATE HISTOGRAMS
	//     TString title_att = "With attenuation #lambda = " + std::to_string(lambda);

	//     draw_histogram(points, reconstructed_points_Att, SiPMTVector, title_att, lambda);

	//     lambda = 0;

	//     draw_histogram(points, reconstructed_points_NoAtt, SiPMTVector, "No Attenuation", lambda);

	//     draw_arrows(points, reconstructed_points_Att);
}
//###################################################################################################
