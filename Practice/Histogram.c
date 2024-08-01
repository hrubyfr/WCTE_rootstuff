#include <iostream>
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"

// Define the user-defined function that creates and returns a histogram
TH1F* CreateHistogram() {
    // Create a TF1 object
    TF1 *myFunc = new TF1("myFunc", "gaus", -5, 5);
    
    // Set initial parameter values
    myFunc->SetParameters(1.0, 0.0, 1.0); // Amplitude, Mean, Sigma
    
    // Create a histogram using the TF1 object
    TH1F *histogram = new TH1F("histogram", "Histogram of my function", 100, -5, 5);
    
    // Fill the histogram with the TF1 function values
    histogram->FillRandom("myFunc", 10000); // Generate 10000 random values
    
    // Clean up
    delete myFunc;

    // Return the histogram pointer
    return histogram;
}

int Histogram() {
    // Call the function to create the histogram
    TH1F *histogram = CreateHistogram();
    
    // Plot the histogram
    TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 600);
    histogram->Draw();
    canvas->Draw();
    
    

    return 0;
}
