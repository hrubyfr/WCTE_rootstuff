#include "TH3D.h"
#include <iostream>
#include "TCanvas.h"

void hist_3D(){
	TH3D* hist3D = new TH3D("hist3D", "3D histogram;x;y;z", 50, -250, 250, 50, -250, 250, 50, -250, 250);
	hist3D->Fill(1, 1, 1);

	TCanvas* can = new TCanvas("", "", 1800, 900);
	hist3D->Draw("BOX");
	

}//code end
