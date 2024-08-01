#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"


class Parameters{
public:
	Double_t par1 = 0.5;
};

Double_t myfun(Double_t *x, Double_t *par){
	Double_t E0 = x[0];
	Double_t theta = par[0];
	Parameters params;

	Double_t f = E0 * params.par1 + theta;
	return f;
}

int TFfun(){
	auto can = new TCanvas();
	auto fun = new TF1("", myfun, 0, 100, 1);
	fun->SetParameter(0, 1);

	fun->Draw("L");
	return 0;
}