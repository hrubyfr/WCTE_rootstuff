#include "Math/Vector2D.h"
#include "Math/Point2D.h"
#include <TMath.h>
#include <iostream>
#include <vector>
#include "TRandom3.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TMarker.h"
#include <cstdio>
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


std::vector<ROOT::Math::XYPoint> Generate_Points(double x, double y, int point_number, double radius){

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
//#####################################################################################################

//_______________________________________DIFFERENT____HISTOGRAMS_______________________________________

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

std::vector<ROOT::Math::XYPoint> Reconstruct_Points(std::vector<std::vector<double>> PMT_signal, std::vector<ROOT::Math::Polar2DVector> SiPMTVector, std::vector<ROOT::Math::XYPoint> points){
    std::vector<ROOT::Math::XYPoint> reconstructed_points;
    ROOT::Math::XYPoint init_point;
    double total_signal;
    for(int i = 0; i < points.size(); i++){
        reconstructed_points.push_back(init_point);
        for (int j = 0; j < SiPMTVector.size(); ++j)
        {
            reconstructed_points[i] += PMT_signal[i][j] * SiPMTVector[j];
            total_signal += PMT_signal[i][j];
            reconstructed_points[i]/=total_signal;
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

    double radius = 3.0;
    int point_number = 1000;

    //### GRID PART ###

    const char* filename = "BeamMonitorGrid.pdf";

    std::vector<std::vector<double>> PMT_signal_NoAtt;

    auto SiPMTVector = Generate_SiPMTVector(radius);
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
    double mean_range[500];
    double lambdas[500];
    int j = 0;
    for (double lambda = 0.1; lambda <= 50.0; lambda += 0.1){
        
        std::vector<std::vector<double>> PMT_signal;
        for (int i = 0; i < points.size(); i++){
            std::vector<double> v1;
            // std::cout << points[i] << endl;
            for (int j = 0; j < SiPMTVector.size(); j++){
                v1.push_back(generate_signal_str(points[i], SiPMTVector[j], lambda));
    }
            PMT_signal.push_back(v1);
    }
        
        auto reconstructed_points_Att = Reconstruct_Points(PMT_signal, SiPMTVector, points);
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
    
    auto graph = new TGraph(500, lambdas, mean_range);
    
    canvas->Clear();

    double original_mean = distance(radius/5, radius/5);
    auto line = new TLine(0, original_mean, 5, original_mean);
    line->SetLineColor(kBlue);
    line->SetLineWidth(2);

    auto text = new TText(5, original_mean, "Original points mean");
    text->SetTextColor(kBlue);
    graph->GetXaxis()->SetTitle("#lambda");
    graph->GetYaxis()->SetTitle("<R>");
    graph->SetTitle("Mean range on attenuation length #lambda");
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