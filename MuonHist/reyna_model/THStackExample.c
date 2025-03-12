#include <TH1F.h>
#include <THStack.h>
#include <TCanvas.h>

void THStackExample() {
    // Create a canvas to draw on
    TCanvas *canvas = new TCanvas("canvas", "THStack Example", 800, 600);

    // Create a THStack
    THStack *stack = new THStack("stack", "Stacked Histograms");

    // Create some example histograms
    TH1F *hist1 = new TH1F("hist1", "Histogram 1", 100, 0, 10);
    TH1F *hist2 = new TH1F("hist2", "Histogram 2", 100, 0, 10);

    // Fill the histograms with some example data
    for (int i = 0; i < 10000; ++i) {
        hist1->Fill(gRandom->Gaus(5, 1));  // Gaussian distribution around 5
        hist2->Fill(gRandom->Exp(1));      // Exponential distribution
    }

    // Set colors for the histograms
    hist1->SetFillColor(kBlue);
    hist2->SetFillColor(kRed);

    // Add histograms to the stack
    stack->Add(hist1);
    stack->Add(hist2);

    // Draw the stack with options to show the legend
    stack->Draw("nostack");  // "nostack" option keeps histograms stacked
    stack->GetXaxis()->SetTitle("X-axis");
    stack->GetYaxis()->SetTitle("Y-axis");
    stack->SetTitle("Stacked Histograms");
    stack->Draw("hist");

    // Add a legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hist1, "Histogram 1", "f");
    legend->AddEntry(hist2, "Histogram 2", "f");
    legend->SetBorderSize(0);
    legend->Draw();

    // Update the canvas
    canvas->Update();
    canvas->Draw();
}
