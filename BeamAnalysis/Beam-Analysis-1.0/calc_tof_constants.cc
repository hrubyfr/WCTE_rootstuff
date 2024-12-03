{

  TFile *fi  = new TFile("tof_hists.root");

	TH1D *hists[32];

	TVectorD  offsets(32);

	TF1 *gauss = new TF1("gauss","[0]*TMath::Exp(-0.5*pow((x-[1])/[2],2))",0,200);

	for(int i=0; i<32; i++){
	  hists[i] = (TH1D*)fi->Get(Form("t_PMT%d",i+28));
		int maxbin = hists[i]->GetMaximumBin();
		double minrange = hists[i]->GetXaxis()->GetBinLowEdge(maxbin-2); 
		double maxrange = hists[i]->GetXaxis()->GetBinUpEdge(maxbin+2); 
		hists[i]->GetXaxis()->SetRangeUser(minrange,maxrange);
		gauss->SetParameter(0,hists[i]->GetMaximum());
		gauss->SetParameter(1,hists[i]->GetXaxis()->GetBinCenter(maxbin));
		gauss->SetParameter(2,1.0);
		hists[i]->Fit("gauss");
		hists[i]->Draw();
		if(hists[i]->GetMean()>10) offsets(i) = gauss->GetParameter(1);
		else offsets(i) = 0.;
	}

	double tof0Mean = 0.;
	double tof1Mean = 0.;
	int tof0iter = 0;
	int tof1iter = 0;

	for(int i=0; i<16; i++){

	  if(offsets(i)>0.){
		  tof0Mean +=offsets(i);
			tof0iter++;
		}
	  if(offsets(i+16)>0.){
		  tof1Mean +=offsets(i+16);
			tof1iter++;
		}
	}

	tof0Mean = tof0Mean/(double)tof0iter;
	tof1Mean = tof1Mean/(double)tof1iter;

	for(int i=0; i<16; i++){

	  if(offsets(i)>0.) offsets(i) = offsets(i)-tof0Mean;
	  if(offsets(i+16)>0.) offsets(i+16) = offsets(i+16)-tof1Mean;

	}

	offsets.Print();

	TFile *fout = new TFile("tof_offsets.root","RECREATE");
	offsets.Write("tof_offsets");
	fout->Close();

}	
