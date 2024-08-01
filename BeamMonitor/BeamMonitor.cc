#include "Math/Vector2D.h"
#include "Math/Point2D.h"
#include <TMath.h>
#include <iostream>
#include <vector>

// Example function to calculate beamline position
// TVector2 GetBeamlinePosition(const std::vector<TVector2>& detectorPositions, const std::vector<double>& signalStrengths) {
//     if (detectorPositions.size() != 8 || signalStrengths.size() != 8) {
//         std::cerr << "Error: Expected 8 detector positions and 8 signal strengths." << std::endl;
//         return TVector2();
//     }

//     double totalSignal = 0;
//     TVector2 weightedSum(0, 0);

//     for (size_t i = 0; i < detectorPositions.size(); ++i) {
//         weightedSum += detectorPositions[i] * signalStrengths[i];
//         totalSignal += signalStrengths[i];
//     }

//     if (totalSignal == 0) {
//         std::cerr << "Error: Total signal strength is zero." << std::endl;
//         return TVector2();
//     }

//     TVector2 beamlinePosition = weightedSum * (1.0 / totalSignal);
//     return beamlinePosition;
// }
double distance(double x, double y){
    return sqrt(x*x + y*y);
}


double generate_signal_str(ROOT::Math::XYPoint point, ROOT::Math::Polar2DVector det_position){
    double max_sig_strength = 1.0;
    auto distance_vector = point - det_position;
    double Sig_strength = max_sig_strength * exp((-1) * distance(distance_vector.X(), distance_vector.Y()));
    std::cout << Sig_strength << endl;
    return Sig_strength;
}




void BeamMonitor() {

    double radius = 1.0;
    std::vector<ROOT::Math::Polar2DVector> SiPMTVector;
    std::vector<double> PMT_signal;

    ROOT::Math::XYPoint point1(0.5, 0.5);

    for (int i = -3; i < 5; i++){
        int pmt_number = i + 3;
        double angle = i * TMath::Pi()/4;
        ROOT::Math::Polar2DVector helpVector(radius, angle);
        SiPMTVector.push_back(helpVector);
        std::cout << "______________SiPMT " << pmt_number << "______________" << endl;
        std::cout << "SiPMTVector check: " << SiPMTVector[pmt_number] << endl;
        std::cout << "SiPMTVector coordinates:\t x:" << SiPMTVector[pmt_number].X() << "\t y:" << SiPMTVector[pmt_number].Y() << endl;
    }
    for (int i = 0; i < SiPMTVector.size(); i++){
        PMT_signal.push_back(generate_signal_str(point1, SiPMTVector[i]));
    }
    ROOT::Math::XYPoint reconstructed_point;
    for(int i = 0; i < SiPMTVector.size(); i++){
        reconstructed_point += PMT_signal[i] * SiPMTVector[i];

    }
    std::cout << "reconstructed point position:" << reconstructed_point << endl;

    

    

    // Example detector positions arranged in a circle
    // double radius = 1.0;
    // std::vector<TVector2> detectorPositions = {
    //     TVector2(radius, 0),
    //     TVector2(radius * TMath::Cos(TMath::Pi() / 4), radius * TMath::Sin(TMath::Pi() / 4)),
    //     TVector2(0, radius),
    //     TVector2(-radius * TMath::Cos(TMath::Pi() / 4), radius * TMath::Sin(TMath::Pi() / 4)),
    //     TVector2(-radius, 0),
    //     TVector2(-radius * TMath::Cos(TMath::Pi() / 4), -radius * TMath::Sin(TMath::Pi() / 4)),
    //     TVector2(0, -radius),
    //     TVector2(radius * TMath::Cos(TMath::Pi() / 4), -radius * TMath::Sin(TMath::Pi() / 4))
    // };

    // // Example signal strengths from the detectors
    // std::vector<double> signalStrengths = {1.0, 2.0, 1.5, 3.0, 1.0, 0.5, 2.5, 2.0};

    // TVector2 beamlinePosition = GetBeamlinePosition(detectorPositions, signalStrengths);
    // std::cout << "Beamline Position: (" 
    //           << beamlinePosition.X() << ", " 
    //           << beamlinePosition.Y() << ")" << std::endl;

    
}
