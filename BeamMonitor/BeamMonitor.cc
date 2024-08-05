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
double distance(double x, double y){
    return sqrt(x*x + y*y);
}




std::vector<ROOT::Math::XYPoint> Generate_Points(int point_number, double radius){
    std::vector<ROOT::Math::XYPoint> points;

    TRandom3 random;
    double x, y;
    while(true){
        x = random.Uniform(-1, 1);
        y = random.Uniform(-1, 1);
        std::cout << "check x = \t" << x << "\t y = " << y << endl;
        if(distance(x, y) < radius){break;}
        
    }
    
    for (int i = 0; i < point_number; ++i){
        ROOT::Math::XYPoint point;
        double random_x = random.Gaus(x, 0.1);
        double random_y = random.Gaus(y, 0.1);
        if (distance(random_x,random_y) >= 1){i--;}
        else{
            point.SetX(random.Gaus(x, 0.01));
            point.SetY(random.Gaus(y, 0.01));
            points.push_back(point);
        }
        
    }

    return points;

}

double generate_signal_str(ROOT::Math::XYPoint point, ROOT::Math::Polar2DVector det_position){
    double max_sig_strength = 1.0;


    auto distance_vector = point - det_position;
    double Sig_strength = max_sig_strength * exp((-1) * distance(distance_vector.X(), distance_vector.Y()));
    // std::cout << Sig_strength << endl;
    return Sig_strength;
}

void BeamMonitor() {

    double radius = 1.0;
    std::vector<ROOT::Math::Polar2DVector> SiPMTVector;
    std::vector<std::vector<double>> PMT_signal;

    int point_number = 10;
    std::vector<ROOT::Math::XYPoint> points = Generate_Points(point_number, radius);

    for (int i = -3; i < 5; i++){
        int pmt_number = i + 3;
        double angle = i * TMath::Pi()/4;
        ROOT::Math::Polar2DVector helpVector(radius, angle);
        SiPMTVector.push_back(helpVector);
        // std::cout << "______________SiPMT " << pmt_number << "______________" << endl;
        // std::cout << "SiPMTVector check: " << SiPMTVector[pmt_number] << endl;
        // std::cout << "SiPMTVector coordinates:\t x:" << SiPMTVector[pmt_number].X() << "\t y:" << SiPMTVector[pmt_number].Y() << endl;
    }
    for (int i = 0; i < points.size(); i++){
        std::vector<double> v1;
        std::cout << points[i] << endl;

        for (int j = 0; j < SiPMTVector.size(); j++){
            v1.push_back(generate_signal_str(points[i], SiPMTVector[j]));

    }


        PMT_signal.push_back(v1);


    }


    std::vector<ROOT::Math::XYPoint> reconstructed_points;
    ROOT::Math::XYPoint init_point;
    for(int i = 0; i < points.size(); i++){
        reconstructed_points.push_back(init_point);
        for (int j = 0; j < SiPMTVector.size(); ++j)
        {
            reconstructed_points[i] += PMT_signal[i][j] * SiPMTVector[j];
        }
        std::cout << "reconstructed point position:" << reconstructed_points[i] << endl;
        

    }
    
    //CREATE HISTOGRAMS

auto points_histogram = new TH2D("", "", 100, -1, 1, 100, -1, 1);
for (int i = 0; i < points.size(); i++){
    points_histogram->Fill(points[i].X(), points[i].Y());
}

auto reconstructed_histogram = new TH2D("", "", 100, -1, 1, 100, -1, 1);
for (int i = 0; i < points.size(); i++){
    reconstructed_histogram->Fill(reconstructed_points[i].X(), reconstructed_points[i].Y());
}
auto circle = new TEllipse(0, 0, 1, 1);
circle->SetFillStyle(0);

auto canvas = new TCanvas();
canvas->Divide(2);
canvas->cd(1);
points_histogram->Draw("colz");
circle->Draw("same");
for (int i = 0; i < SiPMTVector.size(); i++){
    auto pmt = new TMarker(SiPMTVector[i].X(), SiPMTVector[i].Y(), kFullCircle);
    pmt->SetMarkerColor(kRed);
    pmt->Draw("same");
}

canvas->cd(2);
reconstructed_histogram->Draw("colz");
circle->Draw("same");
for (int i = 0; i < SiPMTVector.size(); i++){
    auto pmt = new TMarker(SiPMTVector[i].X(), SiPMTVector[i].Y(), kFullCircle);
    pmt->SetMarkerColor(kRed);
    pmt->Draw("same");
}


    
}
