
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <cmath>

// Constants
const double K = 0.307075; // MeV cm^2 / mol
const double z = 1; // Charge of muon
const double Z_over_A = 0.5; // For example, for silicon
const double me = 0.511; // Electron mass in MeV
const double I = 173e-6; // Mean excitation energy in MeV (example value)
const double density = 1; // g/cm^3 water density
const double c = 3e8; // light speed  
const double mu = 105.66; // approximate muon mass

// Function to calculate dE/dx
double BetheBloch(double beta, double gamma) {
  double bloch = K * density  * 1 / (beta * beta) * Z_over_A * z * z * (0.5 * log(2 * me * c * c * beta * beta * gamma * gamma / I) - beta * beta);
  return bloch;
}

void Bragg_muons() {
	TGraph* graph = new TGraph();
	TGraph* graph2 = new TGraph(); 
	int nPoints = 1000;

	for (int i = 0; i < nPoints; ++i) {
		double beta = 0.01 + i * (0.99 / nPoints);
		double gamma = 1 / sqrt(1 - beta * beta);
		double dEdx = BetheBloch(beta, gamma);
		double Ekin = (gamma - 1) * mu * c * c;
		graph->SetPoint(i, beta * gamma, dEdx);
		graph2->SetPoint (i, Ekin, dEdx);
	}

	TCanvas* canvas = new TCanvas("canvas", "Bethe-Bloch Curve", 1800, 900);
	canvas->Divide(2);
	canvas->cd(1);
	graph->SetTitle("Bethe-Bloch Curve for Muon; #beta#gamma; -dE/dx (MeV / cm)");
	graph->Draw("AL");

	canvas->cd(2);
	graph2->SetTitle("Bethe-Bloch Curve for Muon; Ekin; -dE/dx (MeV / cm)");
	graph2->Draw("AL");

	TGraph* g_range = new TGraph();
	
	


}
