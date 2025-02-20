#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "math.h"
#include "TFile.h"
#include <cstdio>
#include <vector>
#include "TH2D.h"
#include "Math/Vector2D.h"
#include "Math/Point2D.h"

ROOT::Math::Polar2DVector construct_SiPM(double r, double phi){
    return ROOT::Math::Polar2DVector(r, phi);
}

ROOT::Math::XYPoint get_point(){
    return ROOT::Math::XYPoint(2.0, -0.5);
}
double distance(double x, double y){
    return sqrt(pow(x, 2) + pow(y, 2));
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
    std::vector<double> signals;
    auto point = get_point();
    double radius = 3.0;
    for(int i = 0; i < NDet; i++){
        auto r = radius;
        auto phi = i * 2 * M_PI /NDet;
        SiPMs.push_back(construct_SiPM(r, phi));
    }

    std::vector<double> sigma;
    for(int i = 0; i < NDet; i++){
        signals.push_back(Det_response(SiPMs[i], point));
        sigma.push_back(sqrt(signals[i]));
    }
    for (int i = 0; i < NDet; ++i){
        // std::cout << "Detector position " << i << " coordinates are: X: " << SiPMs[i].X() << " Y: " << SiPMs[i].Y() << endl;;
        // std::cout << "Signal in SiPM " << i << " is " << signals[i] << " AU" << endl;
    }

    

    auto chi2function = [&](const double* params) -> double{
        double chi2 = 0;
        double x = params[0];
        double y = params[1];

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

    minimizer->Minimize();
    if (minimizer->Status() != 0) {
    std::cerr << "Warning: Minimization did NOT converge! Status = " 
              << minimizer->Status() << std::endl;
}

    const double* results = minimizer->X();
    std::cout << "Best x: " << results[0] << "\t best y: " << results[1] << endl;    
}
