#include "math.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "TLine.h"
#include "TF1.h"
#include "TStyle.h"
#include "TFile.h"
#include "TBox.h"
#include <GuiTypes.h>
#include <Math/Vector2Dfwd.h>
#include <Rtypes.h>
#include <RtypesCore.h>
#include <TGraph2D.h>
#include "TFile.h"
#include "TH2D.h"
#include "Math/Vector2D.h"
#include "Math/Point2D.h"
#include <TGraph.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TArrow.h>
#include "TLegend.h"
#include "TTree.h"
#include <TVirtualPad.h>
#include <TMarker.h>

#include <cstddef>
#include <fstream>
#include <iterator>
#include <map>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <ostream>
#include <pthread.h>
#include <vector>
#include <set>


using std::array;
using std::string;
using std::vector;
using std::cout;
using std::endl;
using ROOT::Math::XYVector;
using std::set;

void draw_boxes(vector<TBox*> boxes){
	for(int i = 0; i < boxes.size(); i++){
		boxes[i]->SetFillStyle(0);
		boxes[i]->SetLineColor(kBlack);
		boxes[i]->Draw("same");
	}
}
void draw_detector(double size, vector<array<double, 2>> PMTVector){
	TBox* scint = new TBox(-size/2, -size/2, size/2, size/2);
	for (int j = 0; j < PMTVector.size(); j++){
		TMarker *pmt = new TMarker(PMTVector[j][0], PMTVector[j][1], kFullCircle);
		pmt->SetMarkerColor(kBlue);
		pmt->Draw("same");
	}
	scint->SetFillStyle(0);
	scint->SetLineColor(kBlack);
	scint->Draw("same");
}

double distance(double x, double y){
	return sqrt(x * x + y * y);
}


vector<array<double, 2>> get_grid(double size, int max_enum){
	vector<array<double, 2>> helpvector;
	for(int i = 0; i <= max_enum; i++){
		double y = size/2.0 - (size / max_enum) * i;
		for (int j = 0; j <= max_enum; j++){
			double x = -size/2 + j * (size / max_enum);
				array<double, 2> helppoint{x, y};
				helpvector.push_back(helppoint);
			}
		}	
	return helpvector;

}
double Det_response(array<double, 2> PMT_position, array<double, 2> position){
	array<double, 2> rel_position = {PMT_position[0] - position[0], PMT_position[1] - position[1]};
	// std::cout << "Relative position of point to detector is: " << rel_position << endl;
	// std::cout << "Distance of point from detector is " << distance(rel_position.X(), rel_position.Y()) << " cm" << endl;
	double r = distance(rel_position[0], rel_position[1]);
	auto response = 1 / (r * r);
	return response;
}


double chi2function_noamp(const double* params,const vector<double>& signal, const vector<array<double, 2>>& PMT_positions){
	int NDet = 4;
	double chi2 = 0;
	double x = params[0];
	double y = params[1];
	for (int i = 0; i < NDet; i++){
		double dx = x - PMT_positions[i][0];
		double dy = y - PMT_positions[i][1];			//defining the chi2 function	
		double dist = sqrt(dx * dx + dy * dy);

		double theory_sig = 1.0 / (dist * dist);

		double residual = (signal[i] - theory_sig);
		chi2 += residual * residual;

	}
	// cout << "FUMILI was used" << endl;
	return chi2; 
};
double chi2function_amp(const double* params,const vector<double>& signal, const vector<array<double, 2>>& PMT_positions){
	int NDet = 4;
	double chi2 = 0;
	double x = params[0];
	double y = params[1];
	for (int i = 0; i < NDet; i++){
		double dx = x - PMT_positions[i][0];
		double dy = y - PMT_positions[i][1];			//defining the chi2 function	
		double dist = sqrt(dx * dx + dy * dy);

		double theory_sig = 15000 * 1.0 / (dist * dist)* exp(-dist / 120);

		double residual = (signal[i] - theory_sig);
		chi2 += residual * residual;

	}
	// cout << "FUMILI was used" << endl;
	return chi2; 
};
double chi2function_frac(const double* params,const vector<double>& signal, const vector<array<double, 2>>& PMT_positions){
	int NDet = 4;
	double chi2 = 0;
	double x = params[0];
	double y = params[1];
	cout << "Theoretical signals: (";
	double theory_sum = 0;
	vector<double> theory_sig;
	for (int i = 0; i < NDet; i++){
		double dx = x - PMT_positions[i][0];
		double dy = y - PMT_positions[i][1];			//defining the chi2 function	
		double dist = sqrt(dx * dx + dy * dy);
		theory_sig.push_back(1.0 / (dist * dist) * exp(-dist / 120));

		theory_sum += theory_sig[i];

	}
	for(int i = 0; i < NDet; i++){
		theory_sig[i] /= theory_sum;
		cout << theory_sig[i] << ", ";  
		double residual = (signal[i] - theory_sig[i]);
		chi2 += residual * residual;
	}


	cout << ")" << endl;
	cout << "PMT signals: (";
	for (int i = 0; i < 4; i++){
		cout << signal[i] << ", ";
	}
	cout << ")" << endl;
	// cout << "FUMILI was used" << endl;
	return chi2; 
};
//################################################ MAIN FUNCTION ###########################################


void reconstruct_T1(TString filename = "/media/frantisek/T7\ Shield/WCTE_data/WCTE_offline_R1375S0_VME1450.root"){

	int NDet = 4;
	double size = 7.0;
	vector<array<double, 2>> PMTs;

	vector<array<double, 2>> PMT_positions = {{-size/2.0, -size/4.0}, {-size/2.0, size/4.0}, {size/2.0, -size/4.0}, {size/2.0, size/4.0}};
	
	TRandom3 rand;
	
	vector<array<double, 2>> new_coordinates;  //data arrays definition
	vector<double> chi2_values;


	int n_enum = 25;
	auto points = get_grid(size, n_enum);	//Generating grid


	int n_crash = 0;
	double blur = 1;			//defining parameters used in minimization and signal blurring
	double sigma = 0.001;

	TGraph* g_sim_weighed = new TGraph();
	g_sim_weighed->SetTitle("A graph of positions reconstructed via weighed charge method;x[cm];y[cm]");
	for(int i = 0; i < points.size(); i++){
		double x = 0, y = 0;
		auto point = points[i];

		std::vector<double> signals;
		double total_signal = 0;
		for (int j = 0; j < PMT_positions.size(); j++){
			auto signal = Det_response(PMT_positions[j], point);
			x += signal * PMT_positions[j][0];
			y += signal * PMT_positions[j][1];
			total_signal += signal;
			//signal = signal + sigma * rand.Gaus(0, blur); // blurring the signal
			signals.push_back(signal);
		}
		x/=total_signal;
		y/=total_signal;
		g_sim_weighed->AddPoint(x, y);
		auto chi2function = [&](const double* params) -> double{	//inline function declaration of chi2 function, that gets minimized 
			return chi2function_noamp(params, signals, PMT_positions);
		};

		ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

		minimizer->SetTolerance(0.0001);
		minimizer->SetPrintLevel(1);

		ROOT::Math::Functor func(chi2function, 2);

		minimizer->SetFunction(func);

		minimizer->SetVariable(0, "x", 0.0, 0.001);
		minimizer->SetVariable(1, "y", 0.0, 0.001);
		minimizer->SetVariableLimits(0, -size/2.0, size/2.0);
		minimizer->SetVariableLimits(1, -size/2.0, size/2.0);

		minimizer->SetMaxFunctionCalls(1000000);
		minimizer->SetMaxIterations(100000);


		minimizer->Minimize();			//Minimizing the chi2
		if (minimizer->Status() != 0) {
			std::cerr << "Warning: Minimization did NOT converge! Status = " 
				<< minimizer->Status() << std::endl;
			n_crash++;
		}
		double chi2value = minimizer->MinValue();				// Reading out the minimized values
		chi2_values.push_back(chi2value);		


		const double* results = minimizer->X();
		//cout << "Primary point X: " << points[i][0] << " Y: " << points[i][1] << endl;
		//cout << "Best x: " << results[0] << "\t best y: " << results[1] << endl;    

		std::array<double, 2> help_array = {{results[0], results[1]}};
		new_coordinates.push_back(help_array);


	}	//end of event loop
	
	TGraph *g_sim_points = new TGraph();
	TGraph *g_sim_reconstruction = new TGraph();
	for (int i = 0; i < points.size(); i++){
		g_sim_points->AddPoint(points[i][0], points[i][1]);
		g_sim_reconstruction->AddPoint(new_coordinates[i][0], new_coordinates[i][1]);

	}


	TCanvas *c_sim_reconstruction = new TCanvas("c_sim_reconstruction", "c_sim_reconstruction", 1800, 900);
	c_sim_reconstruction->Divide(2);
	c_sim_reconstruction->cd(1);
	g_sim_points->SetTitle("Original points;x[cm];y[cm]");
	g_sim_points->Draw("AP*");
	draw_detector(size, PMT_positions);
	c_sim_reconstruction->cd(2);
	g_sim_reconstruction->SetTitle("Reconstructed points through Minimization;x[cm];y[cm]");
	g_sim_reconstruction->Draw("AP*");
	draw_detector(size, PMT_positions);
	for(int i = 0; i < new_coordinates.size(); i++){
		TArrow* arrow = new TArrow(points[i][0], points[i][1], new_coordinates[i][0], new_coordinates[i][1], 0.005, "|>");
		arrow->Draw();
	}
	TCanvas *c_sim_weighed = new TCanvas("c_sim_weighed", "c_sim_weighed", 1800, 900);
	c_sim_weighed->Divide(2);
	c_sim_weighed->cd(1);
	g_sim_points->Draw("AP*");
	draw_detector(size, PMT_positions);
	c_sim_weighed->cd(2);
	g_sim_weighed->SetMinimum(-4.2);
	g_sim_weighed->SetMaximum(4.2);
	g_sim_weighed->Draw("AP*");
	draw_detector(size, PMT_positions);


	// END SIMULATION, BEGIN ANALYSIS
	//
	bool apply_cuts = true;  // Set true if we want cuts during analysis
	int verbose = 0;

	TFile* file = new TFile(filename, "READ");	//open file


	// ########################### EXTRACT RUN NUMBER FOR NAMING ################
	size_t posR = filename.Index("_R");		
	size_t posS = filename.Index("S", posR);

	int run_number;
	if (posR != TString::kNPOS && posS != TString::kNPOS) {
		// Extract the substring between _R and S
		TString run_number_str = filename(posR + 2, posS - (posR + 2));
		run_number = run_number_str.Atoi(); // Convert to int using Atoi()

		cout << "Run number: " << run_number << endl;
	} else {
		cout << "WRONG FILE FORMAT!" << endl;
		return;
	}
	// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	// ################# LOAD TREES, EXTRACT DATA TO VECTORS ##############################
	TTree* tree = (TTree*) file->Get("WCTEReadoutWindows");
	tree->Print();
	TString plots_folder;

	plots_folder = Form("T1_rec_plots/file_%i", run_number);
	if(gSystem->AccessPathName(plots_folder)) gSystem->Exec("mkdir -p " + plots_folder);  // make a folder to save plots to
	gSystem->cd(plots_folder);

	long nentries = tree->GetEntries();


	cout << "There are " << nentries << " entries" << endl;

	vector<int>* beamline_pmt_qdc_ids = nullptr;
	tree->SetBranchAddress("beamline_pmt_qdc_ids", &beamline_pmt_qdc_ids);
	vector<float>* beamline_pmt_qdc_charge = nullptr;
	tree->SetBranchAddress("beamline_pmt_qdc_charges", &beamline_pmt_qdc_charge);
	vector<float>* beamline_pmt_tdc_times = nullptr;
	tree->SetBranchAddress("beamline_pmt_tdc_times", &beamline_pmt_tdc_times);
	vector<int>* beamline_pmt_tdc_ids = nullptr;
	tree->SetBranchAddress("beamline_pmt_tdc_ids", &beamline_pmt_tdc_ids);
	vector<int>* hit_mpmt_card_ids = nullptr;
	tree->SetBranchAddress("hit_mpmt_card_ids", &hit_mpmt_card_ids);
	vector<int>* hit_pmt_channel_ids = nullptr;
	tree->SetBranchAddress("hit_pmt_channel_ids", &hit_pmt_channel_ids);
	vector<float>* hit_pmt_charges = nullptr;
	tree->SetBranchAddress("hit_pmt_charges", &hit_pmt_charges);
	vector<float>* hit_pmt_times = nullptr;
	tree->SetBranchAddress("hit_pmt_times", &hit_pmt_times);
	vector<double>* trigger_times = nullptr;
	tree->SetBranchAddress("trigger_times", &trigger_times);
	vector<vector<double>>* pmt_waveforms = nullptr;
	tree->SetBranchAddress("pmt_waveforms", &pmt_waveforms);
	vector<int>* pmt_waveform_mpmt_card_ids = nullptr;
	tree->SetBranchAddress("pmt_waveform_mpmt_card_ids", &pmt_waveform_mpmt_card_ids);
	vector<int>* pmt_waveform_pmt_channel_ids = nullptr;
	tree->SetBranchAddress("pmt_waveform_pmt_channel_ids", &pmt_waveform_pmt_channel_ids);
	vector<int>* pmt_waveform_times = nullptr;
	tree->SetBranchAddress("pmt_waveform_times", &pmt_waveform_times);

	cout << "Loading entries" << endl;

	set<int> beam_cards = {130, 131, 132};
	set<int> tof_pmt_ids = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16};
	set<TString> tof_pmt_ids_names = {"TOF_0", "TOF_1", "TOF-2", "TOF-3", "TOF-4", "TOF-5", "TOF_6", "TOF-7", "TOF-8", "TOF-9", "TOF-10", "TOF-11", "TOF-12", "TOF-13", "TOF-14", "TOF-15"};
	set<int> tof_beamline_pmt_ids = {48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63};
	set<int> T0_beamline_pmt_ids = {0, 1, 2, 3};
	set<int> T1_beamline_pmt_ids = {4, 5, 6, 7};
	std::map<int, int> corresponding_pmt = {{0, 8}, {1, 10},{2, 11},{3, 12}, {4, 13},{5, 14},{6, 15},{7, 16},}; 		//Map of opposite SiPMs
	std::map<int, int> inverted_pmt = {{8, 0}, {10, 1},{11, 2},{12, 3}, {13, 4},{14, 5},{15, 6},{16, 7},}; 		//Map of opposite SiPMs
	set<int> T0_BRB_IDs = {12, 13, 14, 15};
	set<int> T1_BRB_IDs = {13, 14, 15, 16};

	// ################ Initialise histograms #########################
	TH2D *h_position_reconstruction = new TH2D("h_position_reconstruction", "Histogram of particle positions reconstructed from charge data by minimization;x[cm];y[cm];count", 50, -size/2, size/2, 50, -size/2, size/2);
	TH2D *h_weighed_reconstruction = new TH2D("h_weighed_reconstruction", "Histogram of particle positions reconstructed from charge data by weighed charge;x[cm];y[cm];count", 50, -size/2, size/2, 50, -size/2, size/2);
	// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	int verb = 1000;

	// ######################### Initialise the map for cuts ############################
	std::map<int, bool> Trigger_hits;

	for (int i = 0; i < 8; i++){
		Trigger_hits[i] = false;		//Fill map with detctor IDs we want to do cuts on (0, 1, 2, 3 - T0; 4, 5, 6, 7 - T1; 43, 44 - T4)
	}
	Trigger_hits[42] = false;
	Trigger_hits[43] = false;

	std::map<int, bool> HC_hits;
	HC_hits[9] = false;				//Fill map with HC IDs - we want to exclude all event when HCs were hit (9 - HC-0; 10 - HC-1)
	HC_hits[10] = false;

	// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	vector<vector<vector<float>>> BRB_pmt_charge(nentries, vector<vector<float>>(16));
	vector<vector<vector<float>>> BRB_pmt_time(nentries, vector<vector<float>>(16));
	vector<vector<vector<double>>> beamline_pmt_time(nentries, vector<vector<double>>(16));
	vector<vector<vector<double>>> beamline_T1_charges(nentries, vector<vector<double>>(4));
	vector<vector<vector<double>>> beamline_T0_time(nentries, vector<vector<double>>(4));
	vector<vector<vector<double>>> beamline_T1_time(nentries, vector<vector<double>>(4));
	vector<vector<double>> TDCT0_hits(nentries);
	vector<vector<double>> TDCT1_hits(nentries);


	for(long ievent = 0; ievent < 1000; ievent++){
		tree->GetEntry(ievent);

		// ################## console print out progress ########################

		if (ievent % verb == 0){
			printf("\rProcessed %ld / %ld events (%.1f%%)", ievent, nentries, 100.0 * ievent / nentries);		//output progress on console
			fflush(stdout);
		}

		// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

		// ######################### CUTS ###################################

		bool pass_cut = true;
		if (verbose){
			cout << endl;
			cout << "#####################################" << endl;
			cout << "Event " << ievent << " beamline data:" << endl;
			cout << "number of TDC event IDs: " << beamline_pmt_tdc_ids->size() << " Number of TDC events: " << beamline_pmt_tdc_times->size() << endl; 
			cout << "number of QDC event IDs: " << beamline_pmt_qdc_ids->size() << " Number of QDC events: " << beamline_pmt_qdc_charge->size() << endl; 
			cout << " TDC channels: ";
			for (int i = 0; i < beamline_pmt_tdc_ids->size(); i++){
				cout << beamline_pmt_tdc_ids->at(i) << " ";
			}
			cout << endl;
			cout << " QDC channels: ";
			for (int i = 0; i < beamline_pmt_qdc_ids->size(); i++){
				cout << beamline_pmt_qdc_ids->at(i) << " ";
			}
			cout << endl;
		}	
		if(apply_cuts){ //set if we want to apply cuts for this analysis
			for (const auto& [ID, was_hit] : Trigger_hits){
				Trigger_hits[ID] = false;
			}							//Initialise all detectors that we cut on to NOT HIT (false)
			for (const auto& [ID, was_hit] : HC_hits){
				HC_hits[ID] = false;
			}


			for (int i = 0; i < beamline_pmt_tdc_ids->size(); i++){
				int tdc_ID = beamline_pmt_tdc_ids->at(i);
				// ######################## T0 CUTS ######################################

				if (tdc_ID == 0) Trigger_hits[0] = true;
				else if (tdc_ID == 1) Trigger_hits[1] = true;
				else if (tdc_ID == 2) Trigger_hits[2] = true;
				else if (tdc_ID == 3) Trigger_hits[3] = true;
				// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

				// ######################## T1 CUTS ######################################

				else if (tdc_ID == 4) Trigger_hits[4] = true;
				else if (tdc_ID == 5) Trigger_hits[5] = true;
				else if (tdc_ID == 6) Trigger_hits[6] = true;
				else if (tdc_ID == 7) Trigger_hits[7] = true;
				// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

				// ######################## T4 CUTS ######################################
			}
			for (int i = 0; i < beamline_pmt_qdc_ids->size(); i++){
				int qdc_id = beamline_pmt_qdc_ids->at(i);

				if (qdc_id == 42){
					if (beamline_pmt_qdc_charge->at(i) > 150){
						Trigger_hits[42] = true;
					}
				} 	

				else if (qdc_id == 43){
					if (beamline_pmt_qdc_charge->at(i) > 100){
						Trigger_hits[43] = true;
					}
				} 	
				// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

				// ######################### HOLE COUNTER CUTS ###########################
				else if (qdc_id == 9){
					if (beamline_pmt_qdc_charge->at(i) > 150) HC_hits[qdc_id] = true;	
				}
				else if (qdc_id == 10){
					if (beamline_pmt_qdc_charge->at(i) > 100) HC_hits[qdc_id] = true;	
				}
				// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


			}

		}
		// xxxxxxxxxxxxxxxxxxxxxxx   END CUTS xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx




		if (!pass_cut) continue;
		// Set trigger times for the beam monitors cards

		// ####################### defining useful variables for analysis #######################

		double trigger_time[3];
		int card_ID;
		int pmt_ID;

		// ######################## VME data analysis ##############################

		double bm_trig_time_0;
		double bm_trig_time_1;
		double Ntrigs_0 = 0;
		double Ntrigs_1 = 0;
		for (int i = 0; i < beamline_pmt_tdc_ids->size(); i++){
			pmt_ID = beamline_pmt_tdc_ids->at(i);

			if (pmt_ID == 31) {
				bm_trig_time_0 = beamline_pmt_tdc_times->at(i); 
				Ntrigs_0++;
			}
			else if (pmt_ID == 46){
				bm_trig_time_1 = beamline_pmt_tdc_times->at(i);
				Ntrigs_1++;
			}
		}
		if (Ntrigs_1 > 1 || Ntrigs_0 > 1) continue; 

		TDCT0_hits[ievent].push_back(bm_trig_time_0);
		TDCT1_hits[ievent].push_back(bm_trig_time_1);

		for (int i = 0; i < beamline_pmt_qdc_charge->size(); i++){
			cout << beamline_pmt_qdc_charge->size() << endl;
			pmt_ID = beamline_pmt_qdc_ids->at(i);
			if (T1_beamline_pmt_ids.count(pmt_ID)){
				cout << "Charge in PMT " << pmt_ID << " is " << beamline_pmt_qdc_ids->at(i) << endl;

				beamline_T1_charges[ievent][pmt_ID - 4].push_back(beamline_pmt_qdc_charge->at(i));
			}
		}

		for (int i = 0; i < beamline_pmt_tdc_ids->size(); i++){
			pmt_ID = beamline_pmt_tdc_ids->at(i);
			double pmt_time = beamline_pmt_tdc_times->at(i);


		}
		// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	} //end loop over all hits

	gStyle->SetPalette(1);
	gStyle->SetOptFit(1);

	for (int ievent = 0; ievent < 1000; ievent++) {
		cout << "Start event " << ievent << endl;
		vector<double> T1_signals;
		bool T1_hit[4];
		for(int i = 0; i < 4; i++){
			T1_hit[i] = false;
		}
		bool T1_all_hit = true;
		if (beamline_T1_charges[ievent].size() < 1) continue;
		for (int pmt_i = 0; pmt_i < beamline_T1_charges[ievent].size(); pmt_i++) {
			for (int i = 0; i < beamline_T1_charges[ievent][pmt_i].size(); i++){

				cout << "Size of event " << ievent << " pmt " << pmt_i << " is " <<  beamline_T1_charges[ievent][pmt_i].size() << endl;
				if(beamline_T1_charges[ievent][pmt_i].size() == 1){
					double T1_charge = beamline_T1_charges[ievent][pmt_i][i];
					T1_signals.push_back(T1_charge);
					cout << "Charge of T1: " << T1_charge << endl;
					T1_hit[pmt_i] = true;
				}
			}
			if (!T1_hit[pmt_i]) T1_all_hit = false;
		}
		if (!T1_all_hit) continue;
		double x_rec = 0;
		double y_rec = 0;
		double signal_sum = 0;
		for (int i = 0; i < T1_signals.size(); i++){
			x_rec += T1_signals[i] * PMT_positions[i][0];
			y_rec += T1_signals[i] * PMT_positions[i][1];
			signal_sum+=T1_signals[i];
		}
		x_rec/= signal_sum;
		y_rec/= signal_sum;
		h_weighed_reconstruction->Fill(x_rec, y_rec);

		for(int i = 0; i < T1_signals.size(); i++){
			T1_signals[i] /= signal_sum;
		}

		auto chi_amp = [&](const double *x){
			return chi2function_frac(x, T1_signals, PMT_positions);
		};

		ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

		min->SetTolerance(10e-10);
		min->SetPrintLevel(1);
		min->SetStrategy(2);

		ROOT::Math::Functor func(chi_amp, 2);
		min->SetFunction(func);

		min->SetVariable(0, "x", 0, 0.1);
		min->SetVariable(1, "y", 0, 0.1);
		//min->SetVariableLimits(0, 0, 200);

		min->SetVariableLimits(0, -size/2, size/2);
		min->SetVariableLimits(1, -size/2, size/2);


		min->SetMaxFunctionCalls(1000000);
		min->SetMaxIterations(100000);
		min->SetPrintLevel(1);

		cout << "Start Minimizing" << endl;
		min->Minimize();
		cout << "Finished Minimizing" << endl;
		if (min->Status() != 0) {
			std::cerr << "Warning: Minimization did NOT converge! Status = " 
				<< min->Status() << endl;
		}
		cout << "Event " << ievent << " finished" << endl;
		const double* results = min->X();
		cout << "X position of particle: " << results[0] << " Y position of particle: " << results[1] << endl;
		double x = results[0];
		double y = results[1];
		h_position_reconstruction->Fill(x, y);

	}
	gStyle->SetOptStat(0);
	// Start drawing histograms
	TCanvas* c_reconstruction = new TCanvas("c_reconstruction", "c_reconstruction", 1800, 900); 
	c_reconstruction->Divide(2);
	c_reconstruction->cd(1);
	gPad->SetRightMargin(0.15);
	h_weighed_reconstruction->GetXaxis()->SetRangeUser(-4, 4);
	h_weighed_reconstruction->GetYaxis()->SetRangeUser(-4, 4);
	h_weighed_reconstruction->Draw("colz");
	draw_detector(size, PMT_positions);


	c_reconstruction->cd(2);

	gPad->SetRightMargin(0.15);
	h_position_reconstruction->GetXaxis()->SetRangeUser(-4, 4);
	h_position_reconstruction->GetYaxis()->SetRangeUser(-4, 4);
	h_position_reconstruction->Draw("colz");
	draw_detector(size, PMT_positions);

	TCanvas* c_partials_reconstruction = new TCanvas("c_partials_reconstruction", "c_partials_reconstruction", 1000, 900); 
	c_partials_reconstruction->SetRightMargin(0.15);
	h_position_reconstruction->Draw("colz");
	draw_detector(size, PMT_positions);

	c_sim_reconstruction->Print("sim_min_reconstruction.pdf");
	c_sim_weighed->Print("sim_weigh_reconstruction.pdf");
	c_reconstruction->Print("sim_reconstruction_charges.pdf");
	c_partials_reconstruction->Print("sim_reconstruction_partial_charges.pdf");
	//Finish drawing histograms



} //end of code
