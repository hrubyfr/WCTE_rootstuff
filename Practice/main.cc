#include <TCanvas.h>
#include <TH1F.h>
#include <string>

// Function to create a histogram and draw it on the provided canvas, adding it to the PDF
void addCanvasToPdf(TCanvas* c, int index, const char* pdfFilename) {
    // Clear the canvas before drawing a new histogram
    c->Clear();
    
    // Create a histogram with unique properties for each call
    TH1F *hist = new TH1F("hist", Form("Histogram %d", index + 1), 100, -4, 4);
    hist->FillRandom("gaus", 1000 + index * 100);  // Adjust filling to vary histograms
    
    // Customize histogram appearance (optional)
    hist->SetLineColor(kBlue + index);
    
    // Draw histogram on the canvas
    hist->Draw();
    c->Update();
    
    // Add this canvas as a new page in the PDF
    c->Print(pdfFilename);
    
    // Clean up histogram
    delete hist;
}

// Main function
void main() {
    // Define the number of loops and the output PDF filename
    int numLoops = 5;
    const char* pdfFilename = "histograms.pdf";
    
    // Create a single canvas and open the PDF for writing
    TCanvas *c = new TCanvas("c", "Canvas", 800, 600);
    c->Print(Form("%s[", pdfFilename));  // Open the PDF

    // Loop to generate multiple histograms and add each to the PDF
    for (int i = 0; i < numLoops; i++) {
        addCanvasToPdf(c, i, pdfFilename);
    }
    
    // Close the PDF after all pages have been added
    c->Print(Form("%s]", pdfFilename));
    
    // Clean up canvas
    delete c;
}