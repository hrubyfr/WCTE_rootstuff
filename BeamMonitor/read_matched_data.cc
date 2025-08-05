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

double calculate_speed(double sigma_tot1, double sigma_tot2, double length1, double length2){
	double v_eff;
	v_eff = sqrt(((length1 * length1 - length2 * length2) / (3 * sigma_tot1 * sigma_tot1 - 3 * sigma_tot2 * sigma_tot2)));
	return v_eff;
}

double calculate_speed_uncertainty(double Sigma1, double Sigma2, double length1, double length2, double Sigma1_error, double Sigma2_error, double v_eff){
	double sigma_d = 0.289;
	double v_eff_sigma = 1/ (3* v_eff) * 1 / (Sigma1 * Sigma1 + Sigma2 * Sigma2) * sqrt(sigma_d * sigma_d * length1 * length1 + sigma_d * sigma_d * length2 * length2 + Sigma1_error * Sigma1_error * Sigma1 * Sigma1 * pow((length1 * length1 - length2 * length2)/(Sigma1 * Sigma1 - Sigma2 * Sigma2), 2) + Sigma2_error * Sigma2_error * Sigma2 * Sigma2 * pow((length1 * length1 - length2 * length2)/(Sigma1 * Sigma1 - Sigma2 * Sigma2), 2));
	return v_eff_sigma;
}

double chi_squared_nolambda(const double *vars, const vector<double>& data, const vector<double>& data_errors, const vector<double>& length){
	const double v_eff = vars[0];
	const double sigma_sipm = vars[1];
	double chi2 = 0;
	for (int i = 0; i < 8; i++){
		double sigma_sq_meas = data[i] * data[i];
		double dim = length[i];
		double sigma_sq_pred = dim * dim / (3 * v_eff * v_eff) + sigma_sipm * sigma_sipm;
		double value = pow((sigma_sq_pred - sigma_sq_meas) / (data_errors[i] * data_errors[i]), 2) ;
		chi2 += value;
	}
	return chi2;
}
double chi_squared(const double *vars, const vector<double>& data, const vector<double>& data_errors, const vector<double>& length){
	const double v_eff = vars[0];
	const double sigma_sipm = vars[1];
	double chi2 = 0;
	for (int i = 0; i < 8; i++){
		double sigma_sq_meas = data[i] * data[i];
		double dim = length[i];
		double sigma_sq_pred = dim * dim / (3 * v_eff * v_eff) + vars[i+2] * sigma_sipm * sigma_sipm;
		double value = pow((sigma_sq_pred - sigma_sq_meas) / (data_errors[i] * data_errors[i]), 2) ;
		chi2 += value;
	}
	return chi2;
}

//################################################ MAIN FUNCTION ###########################################


void read_matched_data(TString filename = "/media/frantisek/T7\ Shield/WCTE_data/WCTE_offline_R1362S0_VME1443.root", const char* out_file = "output.dat", const int energy = 0){

	vector<std::array<double, 2>> scint_dimensions = {{51, 16.25}, 
		{94, 16.25},
		{112, 16.25},
		{123, 16.25},
		{123, 16.25},
		{112, 16.25},
		{94, 16.25},
		{51, 16.25}
	};
	vector<double> start_y = {-58.875, -42.125, -25.375, -8.625, 8.625, 25.375, 42.125, 58.875};
	vector<double> scint_lengths;
	for (const auto pair : scint_dimensions){
		scint_lengths.push_back(pair[0]);
	}
	bool apply_cuts = true;  // Set true if we want cuts during analysis
	int verbose = 0;

	double ACT_g1_cut = 1200;
	double ACT_g2_cut = 1050;
	
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

	if (energy == 0) plots_folder = Form("analysis_plots/matching/file_%i", run_number);
	else plots_folder = Form("analysis_plots/matching/file_%i_E_%i", run_number, energy);
	if(gSystem->AccessPathName(plots_folder)) gSystem->Exec("mkdir -p " + plots_folder);  // make a folder to save plots to
	gSystem->cd(plots_folder);
	
	long nentries = tree->GetEntries();
	//long nentries = 50000;


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

	TH1D* h_beamline_trigger_time_0 = new TH1D("h_beamline_trigger_time_0", "Histogram of trigger times (VME); time[ns]; counts", 500, 0, 500);
	TH1D* h_beamline_trigger_time_1 = new TH1D("h_beamline_trigger_time_1", "Histogram of trigger times (VME); time[ns]; counts", 500, 0, 500);
	TH1D* h_beamline_T0_avgtime = new TH1D("h_beamline_T0_avgtime", "Histogram of T0 average times (VME); time[ns]; counts", 500, 0, 500);
	TH1D* h_beamline_T1_avgtime = new TH1D("h_beamline_T1_avgtime", "Histogram of T1 average times (VME); time[ns]; counts", 500, 0, 500);
	TH1D* h_beamline_T0_counts = new TH1D("h_beamline_T0_counts", "Histogram of T0 counts (VME); counts; ##", 10, 0, 10);
	TH1D* h_beamline_T0_T1_time = new TH1D("h_beamline_T0_T1_time", "Histogram of time of flight from T0 to T1 (VME); T1 - T0 time[ns]; counts", 100, 10, 20);
	vector<TH1D*> h_beamline_T0_times(4);
	for (int i = 0; i < 4; i++){
 		h_beamline_T0_times[i] = new TH1D(Form("h_beamline_T0_times_%i", i), "Histogram of T0 times (VME); time[ns]; counts", 500, 0, 500);
	}
	vector<TH1D*> h_beamline_T1_times(4);
	for (int i = 0; i < 4; i++){
 		h_beamline_T1_times[i] = new TH1D(Form("h_beamline_T1_times_%i", i), "Histogram of T1 times (VME); time[ns]; counts", 500, 0, 500);
	}
	TH1D* h_BRB_trigger_time = new TH1D("h_BRB_trigger_time", "Histogram of trigger times on board 132 (BRB); time[ns]; counts", 500, -500, 2000);
	vector<TH1D*> h_beamline_tof_times(16);
	for (int i = 0; i < 16; i++){
		h_beamline_tof_times[i] = new TH1D(Form("h_beamline_tof_time_%i", i), Form("Time of hit in TOF-%i (VME);time[ns]; counts", i), 500, 0, 500);
	}
	vector<TH1D*> h_BRB_tof_times(16);
	for (int i = 0; i < 16; i++){
		h_BRB_tof_times[i] = new TH1D(Form("h_BRB_tof_time_%i", i), Form("Time of hit in TOF-%i (BRB);time-trigger[ns];count", i), 500, -300, 0);
	}
	vector<TH1D*> h_beamline_diff_time(8);
	for (int i = 0; i < 8; i++){
		h_beamline_diff_time[i] = new TH1D(Form("h_beamline_diff_time_%i", i), 
				Form("Run %i Difference of time hits in TOF-%i and TOF-%i(VME);time TOF-%i - TOF-%i[ns]; counts",run_number, i, i+8, i, i+8), 
				500, -10, 10);
	}
	vector<TH1D*> h_BRB_diff_time(8);
	for (int i = 0; i < 8; i++){
		h_BRB_diff_time[i] = new TH1D(Form("h_BRB_diff_time_%i", i), 
				Form("Run %i Difference of time hits in TOF-%i and TOF-%i (BRB);time TOF-%i - TOF-%i[ns];counts", run_number, i, i+8, i, i+8),
				500, -10, 10);
	}
	vector<TH1D*> h_beamline_T0_tof_time(16);
	for (int i = 0; i < 16; i++){
		h_beamline_T0_tof_time[i] = new TH1D(Form("h_beamline_T0_tof_time_%i", i), 
				Form("Run %i TOF between T0 and TOF-%i (VME);time TOF - T0[ns]; counts",run_number, i), 
				200, 30, 45);
	}
	vector<TH1D*> h_beamline_T1_tof_time(16);
	for (int i = 0; i < 16; i++){
		h_beamline_T1_tof_time[i] = new TH1D(Form("h_beamline_T1_tof_time_%i", i), 
				Form("Run %i TOF between T1 and TOF-%i (VME);time TOF - T1[ns]; counts",run_number, i), 
				200, 0, 40);
	}
	vector<TH2D*> h_beamline_T0_tof_timediff(8);
	for (int i = 0; i < 8; i++){
		h_beamline_T0_tof_timediff[i] = new TH2D(Form("h_beamline_T0_tof_timediff_%i", i), 
				Form("Run %i Histogram of time hits in TOF-%i - T0 and TOF-%i - T0 (VME);time TOF-%i - T0;time TOF-%i - T0;counts", run_number, i, i+8, i, i+8),
				200, 30, 45, 200, 30, 45);
	}
	vector<TH2D*> h_beamline_T1_tof_timediff(8);
	for (int i = 0; i < 8; i++){
		h_beamline_T1_tof_timediff[i] = new TH2D(Form("h_beamline_T1_tof_timediff_%i", i), 
				Form("Run %i Histogram of time hits in TOF-%i - T1 and TOF-%i - T1 (VME);time TOF-%i - T1;time TOF-%i - T1;counts", run_number, i, i+8, i, i+8),
				200, 10, 35, 200, 10, 35);
	}
	vector<TH1D*> h_BRB_T0_tof_time(16);
	for (int i = 0; i < 16; i++){
		h_BRB_T0_tof_time[i] = new TH1D(Form("h_BRB_T0_tof_time_%i", i), 
				Form("Run %i TOF between T0 and TOF-%i (BRB);time TOF - T0[ns]; counts",run_number, i), 
				300, 0, 100);
	}
	vector<TH2D*> h_BRB_T0_tof_timediff(8);
	for (int i = 0; i < 8; i++){
		h_BRB_T0_tof_timediff[i] = new TH2D(Form("h_BRB_T0_tof_timediff_%i", i), 
				Form("Run %i Histogram of time hits in TOF-%i - T0 and TOF-%i - T0 (BRB);time TOF-%i - T0;time TOF-%i - T0;counts", run_number, i, i+8, i, i+8),
				70, 34, 44, 70, 34, 44);
	}
	vector<TH2D*> h_BRB_tof_avg_timediff(8);
	for (int i = 0; i < 8; i++){
		h_BRB_tof_avg_timediff[i] = new TH2D(Form("h_BRB_tof_avg_timediff_%i", i), 
				Form("Run %i Histogram of time hits in TOF-%i - TOF_avg and TOF-%i - TOF_avg (BRB);time TOF-%i - TOF_avg;time TOF-%i - TOF_avg;counts", run_number, i, i+8, i, i+8),
				100, -1, 1, 100, -1, 1);
	}
	vector<TH1D*> h_BRB_avg_tof_tof_time(16);
	for (int i = 0; i < 16; i++){
		h_BRB_avg_tof_tof_time[i] = new TH1D(Form("h_BRB_avg_tof_tof_time_%i", i), 
				Form("Run %i TOF between avg_tof and TOF-%i (BRB);time TOF - avg_tof[ns]; counts",run_number, i), 
				100, -5, 5);
	}
	vector<TH2D*> h_BRB_tof_avg_chargediff(8);
	for (int i = 0; i < 8; i++){
		h_BRB_tof_avg_chargediff[i] = new TH2D(Form("h_BRB_tof_avg_chargediff_%i", i), 
				Form("Run %i Histogram of charge hits in TOF-%i - TOF_avg and TOF-%i - TOF_avg (BRB);charge TOF-%i - TOF_avg;charge TOF-%i - TOF_avg;counts", run_number, i, i+8, i, i+8),
				100, -1000, 1000, 100, -1000, 1000);
	}
	TH2D* h_hits = new TH2D("h_hits", "Histogram of hits in TOF; x[mm]; y[mm]", 50, -scint_dimensions[4][0]/2, scint_dimensions[4][0]/2,
				8, start_y[0] - scint_dimensions[0][1]/2, start_y[7] + scint_dimensions[0][1]/2);

	vector<TH1D*> h_BRB_q_diff(8);
	for(int i = 0; i < 8; i++){
		h_BRB_q_diff[i] = new TH1D(Form("tof_q_diff_%i", i), Form("Run %i diferrence of charge hits in TOF-%i - TOF-%i; Charge TOF-%i - TOF %i ; counts", 
					run_number, i, i+8, i, i+8), 400, -4000, 4000);
	}
	vector<TH2D*> h_BRB_tof_chargediff(8);
	for (int i = 0; i < 8; i++){
		h_BRB_tof_chargediff[i] = new TH2D(Form("h_BRB_tof_chargediff_%i", i), 
				Form("Run %i Histogram of charge hits in TOF-%i and TOF-%i (BRB);charge TOF-%i;charge TOF-%i;counts", run_number, i, i+8, i, i+8),
				100, 0, 6000, 100, 0, 6000);
	}
	vector<TH2D*> h_BRB_total_avg_time_diff(8);
	for(int i = 0; i < 8; i++){
		h_BRB_total_avg_time_diff[i] = new TH2D(Form("tof_total_avg_time_diff_%i", i), Form("Run %i diferrence of time from average time hit in TOF-%i and TOF-%i; Time TOF-%i - avg; Time Tof-%i - avg; counts", 
					run_number, i, i+8, i, i+8), 100, -10, 10, 100, -10, 10);
	}
	vector<TH1D*> h_charge(16);
	for (int i = 0; i < 16; i++){
		h_charge[i] = new TH1D(Form("h_charge_%i", i), 
				Form("Run %i TOF charge on TOF-%i;charge; counts",run_number, i), 
				100, 0, 6000);
	}
	TH1D* h_sigma;
	h_sigma = new TH1D("h_sigma", Form("Run %i Time sigmas of scintillators; scint ## ; sigma[ns]",run_number), 
				16, 0, 8); 

	TH1D* h_speed;
	h_speed = new TH1D("h_speed", Form("Run %i Speed from sigma of shortest and longest scintillators; ; speed[ns]",run_number), 
				1, 0, 1); 

	TH1D* h_mean;
	h_mean = new TH1D("h_mean", Form("Run %i Mean value of t_{L} - t_{R} in a scintillator;scintillator;#mu_{#Deltat}[ns]",run_number), 
				8, 0, 8); 
	TH2D* h_hits_centered = new TH2D("h_hits_centered", "Histogram of hits in TOF; x[mm]; y[mm]", 50, -scint_dimensions[4][0]/2, scint_dimensions[4][0]/2,
				8, start_y[0] - scint_dimensions[0][1]/2, start_y[7] + scint_dimensions[0][1]/2);
	TH1D* h_npairs_hit = new TH1D("h_npairs_hit", Form("Run %i Number of hit pairs in one event;hits;counts", run_number),8, 0, 8);

	TH2D *h_ACTg1_t0tof_tof = new TH2D("h_ACTg1_t0tof_tof", Form("Run %i ACT charge group 1 and time of flight t0 - tof histogram;t0_tof_TOF [ns];ACT_g1 charge", run_number),
			200, 30, 50, 2000, 0, 8000);
	TH2D *h_ACTg2_t0tof_tof = new TH2D("h_ACTg2_t0tof_tof", Form("Run %i ACT charge group 2 and time of flight t0 - tof histogram;t0_tof_TOF [ns];ACT_g2 charge", run_number),
			200, 30, 50, 3000, 0, 14000);
	TH2D *h_ACTg2_t0tof_tof_g1_cut = new TH2D("h_ACTg2_t0tof_tof_g1_cut", Form("Run %i ACT charge group 2 and time of flight t0 - tof histogram after ACT G1 cut;t0_tof_TOF [ns];ACT_g2 charge", run_number),
			200, 30, 50, 3000, 0, 14000);
	vector<TH1D*> h_difftime_pi_cut(8);
	for (int i = 0; i < 8; i++){
		h_difftime_pi_cut[i] = new TH1D(Form("h_difftime_pi_cut_%i", i), 
				Form("Run %i Difference of time hits in TOF-%i and TOF-%i (BRB) of only pions;time TOF-%i - TOF-%i[ns];counts", run_number, i, i+8, i, i+8),
				500, -10, 10);
	}
	vector<TH1D*> h_difftime_mu_cut(8);
	for (int i = 0; i < 8; i++){
		h_difftime_mu_cut[i] = new TH1D(Form("h_difftime_mu_cut_%i", i), 
				Form("Run %i Difference of time hits in TOF-%i and TOF-%i (BRB) of only muons;time TOF-%i - TOF-%i[ns];counts", run_number, i, i+8, i, i+8),
				500, -10, 10);
	}
	vector<TH1D*> h_difftime_e_cut(8);
	for (int i = 0; i < 8; i++){
		h_difftime_e_cut[i] = new TH1D(Form("h_difftime_e_cut_%i", i), 
				Form("Run %i Difference of time hits in TOF-%i and TOF-%i (BRB) of only electrons;time TOF-%i - TOF-%i[ns];counts", run_number, i, i+8, i, i+8),
				500, -10, 10);
	}
	TH1D* h_T0_tof_time = new TH1D("h_T0_tof_time", "Histogram of time of flight from T0 to TOF;TOF - T0 time[ns]; counts", 200, 30, 50);
	TH1D* h_T1_tof_time = new TH1D("h_T1_tof_time", "Histogram of time of flight from T1 to TOF;TOF - T1 time[ns]; counts", 200, 15, 35);

	TH1D* h_pi_sigmas = new TH1D("h_pi_sigmas", "Histogram of sigma values in pi cuts;scintillator;sigma", 8, 0, 8);
	TH1D* h_mu_sigmas = new TH1D("h_mu_sigmas", "Histogram of sigma values in mu cuts;scintillator;sigma", 8, 0, 8);
	TH1D* h_e_sigmas = new TH1D("h_e_sigmas", "Histogram of sigma values in e cuts;scintillator;sigma", 8, 0, 8);


	TH1D* h_pi_charges = new TH1D("h_pi_charges", "Histogram of charge values in pi cuts;charge;counts", 3000, 0, 16000);
	TH1D* h_mu_charges = new TH1D("h_mu_charges", "Histogram of charge values in mu cuts;charge;counts", 3000, 0, 16000);
	TH1D* h_e_charges = new TH1D("h_e_charges", "Histogram of charge values in e cuts;charge;counts", 3000, 0, 16000);

	TH2D* h_timediff_charge_total = new TH2D("h_timediff_charge_total",
				"Histogram of charge in scintillators on time difference in opposite SiPMs;event charge sum;#sigma_{#Delta t}[ns]",
				2500, 0, 16000,
				100, -2, 2);
	
	vector<TH2D*> h_timediff_charge(8);
	for (int i = 0; i < 8; i++){
	h_timediff_charge[i] = new TH2D(Form("h_timediff_charge_%i", i), 
				Form("Histogram of charge in scintillator %i on time difference in opposite SiPMs;event charge sum;#sigma_{#Delta t}[ns]", i),
				2500, 0, 16000,
				100, -2, 2);
	}
	TH2D *h_time_charge_total = new TH2D("h_time_charge_total", 
				"Histogram of charge in all scintillators on time resolution in opposite SiPMs;event charge sum;#sigma_{#Delta t}[ns]",
				2500, 0, 16000,
				100, 0.1, 0.8);
	TH2D *h_timediff_charge_total_ecut = new TH2D("h_timediff_charge_total_ecut", 
				"Histogram of charge of electrons in all scintillators on time resolution in opposite SiPMs;event charge sum;#sigma_{#Delta t}[ns]",
				2500, 0, 16000,
				100, -2, 2);
	TH2D *h_timediff_charge_total_mucut = new TH2D("h_timediff_charge_total_mucut", 
				"Histogram of charge of muons in all scintillators on time resolution in opposite SiPMs;event charge sum;#sigma_{#Delta t}[ns]",
				2500, 0, 16000,
				100, -2, 2);
	TH2D *h_timediff_charge_total_picut = new TH2D("h_timediff_charge_total_picut", 
				"Histogram of charge of pions in all scintillators on time resolution in opposite SiPMs;event charge sum;#sigma_{#Delta t}[ns]",
				2500, 0, 16000,
				100, -2, 2);

	vector<TH2D*>h_time_charge(8);
	for (int i = 0; i < 8; i++){
		h_time_charge[i] = new TH2D(Form("h_time_charge_%i", i), 
				Form("Histogram of charge in scintillator %i on time resolution in opposite SiPMs;event charge sum;#sigma_{#Delta t}[ns]", i),
				2500, 0, 16000,
				50, 0, 1);
	}
	TH1D *h_T0_tof_vme = new TH1D("h_T0_tof_vme", "TOF between T0 and TOF (VME);tof-T0;count", 500, 0, 100);
	vector<TH2D*> h_timediff_charge_ecut(8);
	for (int i = 0; i < 8; i++){
	h_timediff_charge_ecut[i] = new TH2D(Form("h_timediff_charge_ecut_%i", i), 
				Form("Histogram of charge of electrons in scintillator %i on time difference in opposite SiPMs;event charge sum;#sigma_{#Delta t}[ns]", i),
				2500, 0, 16000,
				100, -2, 2);
	}
	vector<TH2D*> h_timediff_charge_mucut(8);
	for (int i = 0; i < 8; i++){
	h_timediff_charge_mucut[i] = new TH2D(Form("h_timediff_charge_mucut_%i", i), 
				Form("Histogram of charge of muons in scintillator %i on time difference in opposite SiPMs;event charge sum;#sigma_{#Delta t}[ns]", i),
				2500, 0, 16000,
				100, -2, 2);
	}
	vector<TH2D*> h_timediff_charge_picut(8);
	for (int i = 0; i < 8; i++){
	h_timediff_charge_picut[i] = new TH2D(Form("h_timediff_charge_picut_%i", i), 
				Form("Histogram of charge of pions in scintillator %i on time difference in opposite SiPMs;event charge sum;#sigma_{#Delta t}[ns]", i),
				2500, 0, 16000,
				100, -2, 2);
	}
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
	vector<vector<vector<double>>> beamline_ACT_charges(nentries, vector<vector<double>>(12));
	vector<vector<vector<double>>> beamline_T0_time(nentries, vector<vector<double>>(4));
	vector<vector<vector<double>>> beamline_T1_time(nentries, vector<vector<double>>(4));
	vector<vector<double>> TDCT0_hits(nentries);
	vector<vector<double>> TDCT1_hits(nentries);


	for(long ievent = 0; ievent < nentries; ievent++){
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

			// ######################### CUT ON A PAIR OF SiPMs HIT #################################

			bool pair_hit = false;
			bool SiPM_hit = false;
			int SiPM_id;	
			auto pt = corresponding_pmt.find(SiPM_id);
			auto pt2 = inverted_pmt.find(SiPM_id);
			for (int i = 0; i < hit_mpmt_card_ids->size(); i++){
				if (hit_mpmt_card_ids->at(i) == 132 && tof_pmt_ids.count(SiPM_id)){
					SiPM_id = hit_pmt_channel_ids->at(i);

					SiPM_hit = true;
					break;
				}
			}
			if(SiPM_hit){
				if(tof_pmt_ids.count(pt->second) || tof_pmt_ids.count(pt2->second)) pair_hit = true;
			} 


			// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
			for (const auto& [ID, was_hit] : Trigger_hits){
				if (!was_hit) pass_cut = false;			 //Do not analyse event if one of the trigger detectors was not hit
			}
			for (const auto& [ID, was_hit] : HC_hits){
				if (was_hit) pass_cut = false;			 //Do not analyse event if a Hole couter was hit
			}
			if(!pair_hit) pass_cut = false;

			//if(!pass_cut) cout << "Event did not pass cut" << endl;
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
				h_beamline_trigger_time_0->Fill(bm_trig_time_0);
				Ntrigs_0++;
			}
			else if (pmt_ID == 46){
				bm_trig_time_1 = beamline_pmt_tdc_times->at(i);
				Ntrigs_1++;
				h_beamline_trigger_time_1->Fill(bm_trig_time_1);
			}
		}
		if (Ntrigs_1 > 1 || Ntrigs_0 > 1) continue; 

		TDCT0_hits[ievent].push_back(bm_trig_time_0);
		TDCT1_hits[ievent].push_back(bm_trig_time_1);

		for (int i = 0; i < beamline_pmt_qdc_charge->size(); i++){
			pmt_ID = beamline_pmt_qdc_ids->at(i);
			if (pmt_ID > 11 && pmt_ID < 24){
				beamline_ACT_charges[ievent][pmt_ID - 12].push_back(beamline_pmt_qdc_charge->at(i));
			}
		}
		for (int i = 0; i < beamline_pmt_tdc_ids->size(); i++){
			pmt_ID = beamline_pmt_tdc_ids->at(i);
			double pmt_time = beamline_pmt_tdc_times->at(i);

			if (T0_beamline_pmt_ids.count(pmt_ID)){
				h_beamline_T0_times[pmt_ID]->Fill(pmt_time);
				if (pmt_time > 150 || pmt_time < 40) continue;
				beamline_T0_time[ievent][pmt_ID].push_back(pmt_time);
			}
			else if (T1_beamline_pmt_ids.count(pmt_ID)){
				h_beamline_T1_times[pmt_ID-4]->Fill(pmt_time);
				if (pmt_time > 100 || pmt_time < 50) continue;
				beamline_T1_time[ievent][pmt_ID-4].push_back(pmt_time);
			}

			// ############################ TOF analysis #############################
			else if (tof_beamline_pmt_ids.count(pmt_ID)){
				beamline_pmt_time[ievent][pmt_ID - 48].push_back(pmt_time);

			}

			// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
		}
		// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

		// ##########################  BRB data analysis ###############################

		for (int i = 0; i < hit_mpmt_card_ids->size(); i++){
			card_ID = hit_mpmt_card_ids->at(i);
			pmt_ID = hit_pmt_channel_ids->at(i);
			if(card_ID == 130 && pmt_ID == 19) trigger_time[0] = hit_pmt_times->at(i);
			else if(card_ID == 131 && pmt_ID == 0) trigger_time[1] = hit_pmt_times->at(i);
			else if(card_ID == 132 && pmt_ID == 19) trigger_time[2] = hit_pmt_times->at(i);
		}
		for(int i = 0; i < hit_mpmt_card_ids->size(); i++){
			pmt_ID = hit_pmt_channel_ids->at(i);
			card_ID = hit_mpmt_card_ids->at(i);
			if(card_ID == 132 && tof_pmt_ids.count(pmt_ID)){
				auto vect_id = std::distance(tof_pmt_ids.begin(), tof_pmt_ids.find(pmt_ID)); 
				BRB_pmt_charge[ievent][vect_id].push_back(hit_pmt_charges->at(i));
				BRB_pmt_time[ievent][vect_id].push_back(hit_pmt_times->at(i) - trigger_time[2]);
			}//do things for mpmt card 132
		} //end loop over mpmt card numbers


	}// end of loop over events

	// ######################### ANALYSIS OF OPPOSITE SiPMs ###################################
	cout << endl;
	cout << "Starting analysis of SiPM pairs" << endl;
	int cout_events = 0;

	for (int event = 0; event < beamline_pmt_time.size(); event++){
		for (int iscint = 0; iscint < 8; iscint++){
			for (int i = 0; i < beamline_pmt_time[event][iscint].size(); i++){
				double sipm_time_a = beamline_pmt_time[event][iscint][i];
				if(sipm_time_a < 70.0 || sipm_time_a > 120.0) continue;
				for (int j = 0; j < beamline_pmt_time[event][iscint+8].size(); j++){
					double sipm_time_b = beamline_pmt_time[event][iscint+8][j];
					if(sipm_time_b < 70.0 || sipm_time_b > 120.0) continue;
					h_beamline_diff_time[iscint]->Fill(sipm_time_a - sipm_time_b);
				}
			}
		}
	}
	for (int ievent = 0; ievent < beamline_pmt_time.size(); ievent++){ 
		for (int iPM = 0; iPM < beamline_pmt_time[ievent].size(); iPM++){
			for (int i = 0; i < beamline_pmt_time[ievent][iPM].size(); i++){
				h_beamline_tof_times[iPM]->Fill(beamline_pmt_time[ievent][iPM][i]);
			}
		}
	}

	for (int ievent = 0; ievent < beamline_T0_time.size(); ievent++){
		if (ievent < cout_events){
			cout << endl; 
			cout << "########################### event " << ievent << " ###################" << endl;
		}

		double avg_time_T0 = 0;
		double avg_time_T1 = 0;
		int npmt_hit_T0 = 0;
		int npmt_hit_T1 = 0;
		for (int i = 0; i < beamline_T0_time[ievent].size(); i++){
			for (int j = 0; j < beamline_T0_time[ievent][i].size(); j++){
				avg_time_T0 += beamline_T0_time[ievent][i][j];
				npmt_hit_T0++;
				if (ievent < cout_events) cout << " Time in T0-" << i << " is " << beamline_T0_time[ievent][i][j] << endl;
			}
			for (int j = 0; j < beamline_T1_time[ievent][i].size(); j++){

				avg_time_T1 += beamline_T1_time[ievent][i][j];
				npmt_hit_T1++;
				if (ievent < cout_events) cout << " Time in T1-" << i << " is " << beamline_T1_time[ievent][i][j] << endl;
			}

		}
		avg_time_T0 /= npmt_hit_T0;
		avg_time_T1 /= npmt_hit_T1;
		h_beamline_T0_T1_time->Fill(avg_time_T1 - avg_time_T0);
		if(ievent < cout_events){
			cout << "Number of counts is " << npmt_hit_T0 << endl;
			cout << "Average time is " << avg_time_T0 << endl;
		}
		h_beamline_T0_avgtime->Fill(avg_time_T0);
		h_beamline_T0_counts->Fill(npmt_hit_T0);
		for (int i = 0; i < beamline_pmt_time[ievent].size() ; i++){
			for(int j = 0; j < beamline_pmt_time[ievent][i].size() ; j++){
				h_beamline_T0_tof_time[i]->Fill(beamline_pmt_time[ievent][i][j] - TDCT1_hits[ievent][0] - (avg_time_T0 - TDCT0_hits[ievent][0]));

				h_beamline_T1_tof_time[i]->Fill(beamline_pmt_time[ievent][i][j] - TDCT1_hits[ievent][0] - (avg_time_T1 - TDCT0_hits[ievent][0]));
			}
		}
		for (int i = 0; i < beamline_pmt_time[ievent].size()/2; i++){
			for (int j = 0; j < beamline_pmt_time[ievent][i].size(); j++){
				if (beamline_pmt_time[ievent][i].size() > 1) continue;
				double sipm_time_a = beamline_pmt_time[ievent][i][j];
				for (int k = 0; k < beamline_pmt_time[ievent][i+8].size(); k++){
					if (beamline_pmt_time[ievent][i+8].size() > 1) continue;
					double sipm_time_b = beamline_pmt_time[ievent][i+8][k];
					double avg_time = (sipm_time_a + sipm_time_b) / 2;
					h_T0_tof_vme->Fill(avg_time - TDCT1_hits[ievent][0] - (avg_time_T0 - TDCT0_hits[ievent][0]));
				}

			}
		}
		for (int i = 0; i < beamline_pmt_time[ievent].size() /2 ; i++){
			for(int j = 0; j < beamline_pmt_time[ievent][i].size() ; j++){
				if (ievent < cout_events) cout << "SiPM " << i << " time: " << beamline_pmt_time[ievent][i][j] << endl;
				for(int k = 0; k < beamline_pmt_time[ievent][i + 8].size(); k++){
					if (ievent < cout_events) cout << "SiPM " << i+8 << " time: " << beamline_pmt_time[ievent][i+8][k] << endl;
					double tof0_trig_time = beamline_pmt_time[ievent][i][j] - TDCT1_hits[ievent][0];
					double tof1_trig_time = beamline_pmt_time[ievent][i+8][k] - TDCT1_hits[ievent][0];
					double T0_trig_time = avg_time_T0 - TDCT0_hits[ievent][0];
					double T1_trig_time = avg_time_T1 - TDCT0_hits[ievent][0];
					h_beamline_T0_tof_timediff[i]->Fill(tof0_trig_time - T0_trig_time, tof1_trig_time - T0_trig_time);
					h_beamline_T1_tof_timediff[i]->Fill(tof0_trig_time - T1_trig_time, tof1_trig_time - T1_trig_time);
					if (ievent < cout_events){
						cout << "Time of flight between T0 and TOF-"<< i << " is " << beamline_pmt_time[ievent][i][j] - avg_time_T0 << endl;
						cout << "Time of flight between T0 and TOF-"<< i+8 << " is " << beamline_pmt_time[ievent][i+8][k] - avg_time_T0 << endl;
					}
				}
			}
		}
	}

	cout << "Finished VME data analysis, starting BRB data analysis" << endl;


	// ######################## BRB opposite SIPMs analysis ######################################

	vector<vector<double>> scint_hits(8);

	for (int ievent = 0; ievent < beamline_ACT_charges.size(); ievent++){
		double ACT_g1_charge = 0;
		double ACT_g2_charge = 0;
		if (TDCT0_hits[ievent].size() < 1) continue;
		double avg_time_T0 = 0;
		int npmt_hit_T0 = 0;
		double avg_time_T1 = 0;
		int npmt_hit_T1 = 0;

		for (int i = 0; i < beamline_T0_time[ievent].size(); i++){
			for (int j = 0; j < beamline_T0_time[ievent][i].size(); j++){
				avg_time_T0 += beamline_T0_time[ievent][i][j];
				npmt_hit_T0++;
			}
		}
		for (int i = 0; i < beamline_T1_time[ievent].size(); i++){
			for (int j = 0; j < beamline_T1_time[ievent][i].size(); j++){
				avg_time_T1 += beamline_T1_time[ievent][i][j];
				npmt_hit_T1++;
			}
		}

		avg_time_T0 /= npmt_hit_T0;
		avg_time_T1 /= npmt_hit_T1;

		double T0_time = avg_time_T0 - TDCT0_hits[ievent][0];
		double T1_time = avg_time_T1 - TDCT0_hits[ievent][0];

		for (int i = 0; i < 6; i++){
			for (int icharge = 0; icharge < beamline_ACT_charges[ievent][i].size(); icharge++){
				ACT_g1_charge+=beamline_ACT_charges[ievent][i][icharge];
			}
		}
		for (int i = 6; i < 12; i++){
			for (int icharge = 0; icharge < beamline_ACT_charges[ievent][i].size(); icharge++){
				ACT_g2_charge+=beamline_ACT_charges[ievent][i][icharge];
			}
		}

		for (int iscint = 0; iscint < 8; iscint++){
			double avg_SiPM_time;
			for (int i = 0; i < BRB_pmt_time[ievent][iscint].size(); i++){
				double sipm_time_a = BRB_pmt_time[ievent][iscint][i];
				double sipm_charge_a = BRB_pmt_charge[ievent][iscint][i];
				bool pair_hit = true;
				if (sipm_time_a < -165 || sipm_time_a > -145) continue;
				for (int j = 0; j < BRB_pmt_time[ievent][iscint+8].size(); j++){
					double sipm_time_b = BRB_pmt_time[ievent][iscint+8][j];
					double sipm_charge_b = BRB_pmt_charge[ievent][iscint+8][j];
					if (sipm_time_b < -165 || sipm_time_b > -145) continue;
					double avg_time = (sipm_time_a + sipm_time_b)/2;
					double t0tof_time = avg_time - T0_time;
					double t1tof_time = avg_time - T1_time;
					if (pair_hit) {
						h_ACTg1_t0tof_tof->Fill(t0tof_time, ACT_g1_charge);
						h_ACTg2_t0tof_tof->Fill(t0tof_time, ACT_g2_charge);
						h_T0_tof_time->Fill(t0tof_time);
						h_T1_tof_time->Fill(t1tof_time);
						if (ACT_g1_charge < ACT_g1_cut) h_ACTg2_t0tof_tof_g1_cut->Fill(t0tof_time, ACT_g2_charge);
						if (ACT_g1_charge < ACT_g1_cut && ACT_g2_charge < ACT_g2_cut){
							h_difftime_pi_cut[iscint]->Fill(sipm_time_a - sipm_time_b);
						}
						if (ACT_g1_charge < ACT_g1_cut && ACT_g2_charge > ACT_g2_cut){
							h_difftime_mu_cut[iscint]->Fill(sipm_time_a - sipm_time_b);
						}
						if (ACT_g1_charge > ACT_g1_cut){
							h_difftime_e_cut[iscint]->Fill(sipm_time_a - sipm_time_b);
						}
					}
					pair_hit=false;
				}
			}
		}
	}


	for (int event = 0; event < BRB_pmt_time.size(); event++){
		for (int iscint = 0; iscint < 8; iscint++){
			for (int i = 0; i < BRB_pmt_time[event][iscint].size(); i++){
				double sipm_time_a = BRB_pmt_time[event][iscint][i];
				if (sipm_time_a < -165 || sipm_time_a > -145) continue;
				scint_hits[iscint].push_back(sipm_time_a);
				for (int j = 0; j < BRB_pmt_time[event][iscint+8].size(); j++){
					double sipm_time_b = BRB_pmt_time[event][iscint+8][j];
					if (sipm_time_b < -165 || sipm_time_b > -145) continue;
					scint_hits[iscint].push_back(sipm_time_b);
				}
			}
		}
	}

	for (int event = 0; event < BRB_pmt_time.size(); event++){
		int n_pairs_hit = 0;
		for (int iscint = 0; iscint < 8; iscint++){
			for (int i = 0; i < BRB_pmt_time[event][iscint].size(); i++){
				bool pair_hit = true;
				double sipm_time_a = BRB_pmt_time[event][iscint][i];
				if (sipm_time_a < -165 || sipm_time_a > -145) continue;
				for (int j = 0; j < BRB_pmt_time[event][iscint+8].size(); j++){
					double sipm_time_b = BRB_pmt_time[event][iscint+8][j];
					if (sipm_time_b < -165 || sipm_time_b > -145) continue;
					h_BRB_diff_time[iscint]->Fill(sipm_time_a - sipm_time_b);
					if (pair_hit) n_pairs_hit++;
					pair_hit = false;

				}
			}
		}
		h_npairs_hit->Fill(n_pairs_hit);
	}
	cout << "Number of no pair hits: "<< h_npairs_hit->GetBinContent(1) << " Number of 1 pairs hit: " << h_npairs_hit->GetBinContent(2) << endl;
	cout << "Ratio of no hits to 1 pair hit: " << h_npairs_hit->GetBinContent(2) / h_npairs_hit->GetBinContent(1) << endl;
	double total_avg_time[8];
	for (int iscint = 0; iscint < 8; iscint++){
		for (int i = 0; i < scint_hits[iscint].size(); i++){
			total_avg_time[iscint] += scint_hits[iscint][i];

		}
		total_avg_time[iscint] /= scint_hits[iscint].size();
	}
	for (int ievent = 0; ievent < BRB_pmt_time.size(); ievent++){
		if (verbose){
			cout << endl;
			cout << " ##################  EVENT " << ievent << "  ######################" << endl;
		}
		if (TDCT0_hits[ievent].size() < 1) continue;
		double avg_time_T0 = 0;
		int npmt_hit_T0 = 0;
		for (int i = 0; i < beamline_T0_time[ievent].size(); i++){
			for (int j = 0; j < beamline_T0_time[ievent][i].size(); j++){
				avg_time_T0 += beamline_T0_time[ievent][i][j];
				npmt_hit_T0++;
			}
		}
		avg_time_T0 /= npmt_hit_T0;
		double T0_time = avg_time_T0 - TDCT0_hits[ievent][0];
		double avg_sipm_time;
		for (int iscint = 0; iscint < 8; iscint++){
			for (int i = 0; i < BRB_pmt_time[ievent][iscint].size(); i++){
				if (BRB_pmt_time[ievent][iscint].size() < 1) continue;
				double sipm_time_a = BRB_pmt_time[ievent][iscint][i];
				for (int j = 0; j < BRB_pmt_time[ievent][iscint+8].size(); j++){
					double sipm_time_b = BRB_pmt_time[ievent][iscint+8][j];
					avg_sipm_time = (sipm_time_a + sipm_time_b)/2;
					//cout << "h_BRB_T0_tof_timediff content: " << sipm_time_a - T0_time << " " << sipm_time_b - T0_time << endl;
					h_BRB_T0_tof_timediff[iscint]->Fill(sipm_time_a - T0_time, sipm_time_b - T0_time);
					h_BRB_tof_avg_timediff[iscint]->Fill(sipm_time_a - avg_sipm_time, sipm_time_b - avg_sipm_time);
					h_BRB_total_avg_time_diff[iscint]->Fill(sipm_time_a - total_avg_time[iscint], sipm_time_b - total_avg_time[iscint]);
					double position = (sipm_time_a - sipm_time_b) * 185 / 2.0 ;
					h_hits->Fill(position, start_y[iscint]);
				}
			}
		}
		for (int iPM = 0; iPM < BRB_pmt_time[ievent].size(); iPM++){
			for (int i = 0; i < BRB_pmt_time[ievent][iPM].size(); i++){
				//cout << "h_BRB_T0_tof_times content: " << BRB_pmt_time[ievent][iPM][i] - T0_time << endl;
				h_BRB_tof_times[iPM]->Fill(BRB_pmt_time[ievent][iPM][i]);
				h_BRB_T0_tof_time[iPM]->Fill(BRB_pmt_time[ievent][iPM][i] - T0_time);
				h_BRB_avg_tof_tof_time[iPM]->Fill(BRB_pmt_time[ievent][iPM][i] - avg_sipm_time);
			}
		}
		if(verbose) cout << "Size of charge vector" << BRB_pmt_charge[ievent].size() << endl;
		for (int iscint = 0; iscint < BRB_pmt_charge[ievent].size() / 2; iscint++){
			if (verbose) cout << "accessing " << iscint << " scintillator" << endl;
			if (verbose) cout << "sipm_a size: " << BRB_pmt_charge[ievent][iscint].size() << endl;
			for (int sipm_a = 0; sipm_a < BRB_pmt_charge[ievent][iscint].size(); sipm_a++){
				float sipm_a_charge = BRB_pmt_charge[ievent][iscint][sipm_a];
				if (verbose) cout << "sipm_b size: " << BRB_pmt_charge[ievent][iscint + 8].size() << endl;
				for (int sipm_b = 0; sipm_b < BRB_pmt_charge[ievent][iscint + 8].size(); sipm_b++){
					float sipm_b_charge = BRB_pmt_charge[ievent][iscint+8][sipm_b];
					float avg_charge = (sipm_b_charge + sipm_a_charge) / 2;
					/*if (iscint + 8 == 10){
					  if(sipm_a_charge < 800 || sipm_b_charge < 400) continue;
					  }
					  else if (sipm_a_charge < 800 || sipm_b_charge < 600) continue;
					  */
					h_BRB_tof_avg_chargediff[iscint]->Fill(sipm_a_charge - avg_charge, sipm_b_charge - avg_charge);
					if (verbose) cout << "SiPM_A charge: " << sipm_a_charge << " SiPM_B charge: " << sipm_b_charge << endl;
					h_BRB_q_diff[iscint]->Fill(sipm_a_charge - sipm_b_charge);
					h_BRB_tof_chargediff[iscint]->Fill(sipm_a_charge, sipm_b_charge);
				}
			}
		}
		for (int iPM = 0; iPM < 16; iPM++){
			for (int i = 0; i < BRB_pmt_charge[ievent][iPM].size(); i++){
				h_charge[iPM]->Fill(BRB_pmt_charge[ievent][iPM][i]);
			}
		}
	}


	cout << "Analysis done, now drawing histograms" << endl;

	// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	gStyle->SetPalette(1);
	gStyle->SetOptFit(1);

	TCanvas *c_vme_tof = new TCanvas("c_T0_tof", "c_T0_tof", 1800, 900);
	h_T0_tof_vme->Draw("hist");
	c_vme_tof->Print("vme_T0_tof.png");

	TCanvas *c_T1_tof = new TCanvas("c_T1_tof", "c_T1_tof", 1800, 900);
	c_T1_tof->SetLogy();
	h_T1_tof_time->Draw("hist");
	c_T1_tof->Print("h_T1_tof_time.png");

	vector<double> fit_BRB_sigma(8);
	vector<double> fit_BRB_sigma_error(8);
	vector<double> fit_BRB_mean(8);

	vector<double> fit_pi_sigma(8);
	vector<double> fit_pi_sigma_error(8);
	vector<double> fit_pi_mean(8);

	vector<double> fit_mu_sigma(8);
	vector<double> fit_mu_sigma_error(8);
	vector<double> fit_mu_mean(8);

	vector<double> fit_e_sigma(8);
	vector<double> fit_e_sigma_error(8);
	vector<double> fit_e_mean(8);

	TCanvas* c_npairs_hit = new TCanvas("c_npairs_hit", "c_npairs_hit", 1800, 900);
	h_npairs_hit->Draw("hist");

	c_npairs_hit->Print("h_npairs_hit.pdf");


	TCanvas *c_T0tof_ACTg1 = new TCanvas("c_T0tof_ACTg1", "c_T0tof_ACTg1", 1800, 900);
	h_ACTg1_t0tof_tof->Draw("colz");
	c_T0tof_ACTg1->Print("h_ACTg1_T0tof_tof.png");

	TCanvas *c_T0tof_ACTg2 = new TCanvas("c_T0tof_ACTg2", "c_T0tof_ACTg2", 1800, 900);
	h_ACTg2_t0tof_tof->Draw("colz");
	c_T0tof_ACTg2->Print("h_ACTg2_T0tof_tof.png");

	TCanvas *c_T0tof_ACTg2_g1cut = new TCanvas("c_T0tof_ACTg2_g1cut", "c_T0tof_ACTg2_g1cut", 1800, 900);
	h_ACTg2_t0tof_tof_g1_cut->Draw("colz");
	c_T0tof_ACTg2_g1cut->Print("h_ACTg2_g1cut_T0tof_tof.png");

	c_T0tof_ACTg2_g1cut->Clear();
	h_ACTg2_t0tof_tof_g1_cut->GetYaxis()->SetRangeUser(0, 2000);
	h_ACTg2_t0tof_tof_g1_cut->Draw("colz");
	c_T0tof_ACTg2_g1cut->Print("h_ACTg2_g1cut_T0tof_tof_zoom.png");

	TCanvas* c_bm_diff_times = new TCanvas("c_bm_diff_times", "c_bm_diff_times", 1800, 900);
	c_bm_diff_times->Divide(4, 2);
	for (int i = 0; i < 8; i++){
		c_bm_diff_times->cd(i+1);
		double x_max = h_beamline_diff_time[i]->GetBinCenter(h_beamline_diff_time[i]->GetMaximumBin());
		h_beamline_diff_time[i]->Fit("gaus", "", "RS", x_max - 2, x_max + 2);
		TF1* fit_fun = h_beamline_diff_time [i]->GetFunction("gaus");
		h_beamline_diff_time[i]->Draw("hist");
		fit_fun->Draw("same");
	}
	c_bm_diff_times->Print("h_bm_diff_times.pdf");

	TCanvas* c_brb_diff_times = new TCanvas("c_brb_diff_times", "c_brb_diff_times", 1800, 900);
	c_brb_diff_times->Divide(4, 2);
	for (int i = 0; i < 8; i++){
		c_brb_diff_times->cd(i+1);
		double x_max = h_BRB_diff_time[i]->GetBinCenter(h_BRB_diff_time[i]->GetMaximumBin());
		h_BRB_diff_time[i]->Fit("gaus", "", "RS", x_max - 2, x_max + 2);
		TF1* fit_fun = h_BRB_diff_time [i]->GetFunction("gaus");
		h_BRB_diff_time[i]->Draw("hist");
		fit_fun->Draw("same");
		fit_BRB_sigma[i] = fit_fun->GetParameter(2);
		fit_BRB_sigma_error[i] = fit_fun->GetParError(2);
		fit_BRB_mean[i] = fit_fun->GetParameter(1);
	}
	c_brb_diff_times->Print("h_brb_diff_times.pdf");


	TCanvas* c_speed_fromavg = new TCanvas("c_speeds_fromavg", "c_speeds_fromavg", 1800, 900);
	double avg1 = (fit_BRB_sigma[0] + fit_BRB_sigma[7]) / 2;
	double avg2 = (fit_BRB_sigma[3] + fit_BRB_sigma[4]) / 2;
	double avg1_unc = sqrt(fit_BRB_sigma_error[0] * fit_BRB_sigma_error[0] + fit_BRB_sigma_error[7] * fit_BRB_sigma_error[7]) / 2;
	double avg2_unc = sqrt(fit_BRB_sigma_error[3] * fit_BRB_sigma_error[3] + fit_BRB_sigma_error[4] * fit_BRB_sigma_error[4]) / 2;
	double avg_speed = calculate_speed(avg1, avg2, scint_dimensions[0][0], scint_dimensions[3][0]);
	double v_eff_sigma = calculate_speed_uncertainty(avg1, avg2, scint_dimensions[0][0], scint_dimensions[3][0], avg1_unc, avg2_unc, avg_speed);

	h_speed->SetBinContent(1, avg_speed);
	h_speed->SetBinError(1, v_eff_sigma);

	h_speed->Draw("histe");

	TCanvas* c_brb_sigmas = new TCanvas("c_brb_sigmas", "c_brb_sigmas", 1800, 900);
	h_sigma->SetStats(0);
	for(int i = 0; i < 8; i+=2){
		h_sigma->SetBinContent(2*i + 1, fit_BRB_sigma[i/2]);
		h_sigma->SetBinContent(2*i + 3, fit_BRB_sigma[7 - i/2]);
		h_sigma->SetBinError(2*i + 1, fit_BRB_sigma_error[i/2]);
		h_sigma->SetBinError(2*i + 3, fit_BRB_sigma_error[7 - i/2]);
		TString label = Form("scintillator %i and %i", i/2, 7 - i/2);
		h_sigma->GetXaxis()->SetBinLabel(2*i + 1, label);
	}
	h_sigma->Draw("histe");
	c_brb_sigmas->Print("h_sigmas.pdf");

	cout << "Fit brb sigma 7: " << fit_BRB_sigma[7] << " fit brb sigma 4: " << fit_BRB_sigma[4] << " scint_dimensions 7: " << scint_dimensions[7][0] << " scint_dimensions 4: " << scint_dimensions[4][0] << endl;

	vector<TH1D*> h_speed2(8);
	int num = 1;
	for (int i = 0; i  < 8; i++){

		h_speed2[i] = new TH1D("", "", 10, 0, 10);
		h_speed2[i]->SetName(Form("Effective speed calculated from comparing scintillator %i with others", i));
		h_speed2[i]->SetTitle(Form("Effective speed calculated from comparing scintillator %i with others", i));
		h_speed2[i]->GetXaxis()->SetTitle("number");
		h_speed2[i]->GetYaxis()->SetTitle("speed [mm/ns]");
		num = 1;

		for (int j = 0; j < 8; j++){
			if(i==j || i== 7-j) continue;
			double value = calculate_speed(fit_BRB_sigma[i], fit_BRB_sigma[j], scint_dimensions[i][0], scint_dimensions[j][0]);
			if(std::isnan(value) || std::isinf(value) || value == 0) continue;
			h_speed2[i]->SetBinContent(num, value);
			num++;
			//cout << endl;
			//cout << "Check velocity value: " << value << endl;



		}
		cout << " NUM: " << num << endl;
		//	cout << "Check graph point number: N = " <<h_speed2[i]->GetN() << endl; 
	}
	cout << "ChecK: graphs done " << endl;
	TCanvas* c_speeds = new TCanvas("c_speeds", "c_speeds", 1800, 900);
	c_speeds->Divide(4, 2);
	for(int i = 0; i < 8; i++){
		c_speeds->cd(i+1);
		h_speed2[i]->Draw("APL");
	}
	c_speeds->Print("c_speeds.pdf");

	vector<TBox*> scints;
	for(int i = 0; i < scint_dimensions.size(); i++){
		TBox* help_box = new TBox(-scint_dimensions[i][0]/2, start_y[i] - scint_dimensions[i][1]/2, scint_dimensions[i][0]/2, start_y[i] + scint_dimensions[i][1]/2) ;
		scints.push_back(help_box);
	}//end create scintillator tiles
	 //

	TCanvas* c_BRB_T0_tof_difftimes = new TCanvas("c_BRB_T0_tof_difftimes", "c_BRB_T0_tof_difftimes", 1800, 900);
	c_BRB_T0_tof_difftimes->Divide(4, 2);
	for (int i = 0; i < 8; i++){
		c_BRB_T0_tof_difftimes->cd(i+1);
		h_BRB_T0_tof_timediff[i]->Draw("colz");
	}
	c_BRB_T0_tof_difftimes->Print("h_BRB_T0_tof_difftimes.pdf");


	for (int i = 0; i < 8; i++){
		h_mean->SetBinContent(i+1, fit_BRB_mean[i]);
		h_mean->SetBinError(i+1, fit_BRB_sigma[i]);
	}
	TCanvas* c_means = new TCanvas("c_means", "c_means", 1800, 900);
	h_mean->Draw("e1");
	TLine* line = new TLine(0, 0, 8, 0);
	line->SetLineColor(kRed);
	line->Draw("same");
	c_means->Print("c_means.pdf");

	TCanvas* c_hit_rec = new TCanvas("", "", 1800, 900);
	h_hits->Draw("colz");
	draw_boxes(scints);
	c_hit_rec->Print("c_position_reconstruction.pdf");


	for(int ievent = 0; ievent < BRB_pmt_time.size(); ievent++){
		for (int iscint = 0; iscint < 8; iscint++){
			double scint_bias = fit_BRB_mean[iscint];
			for (int i = 0; i < BRB_pmt_time[ievent][iscint].size(); i++){
				if (BRB_pmt_time[ievent][iscint].size() < 1) continue;
				double sipm_time_a = BRB_pmt_time[ievent][iscint][i];
				for (int j = 0; j < BRB_pmt_time[ievent][iscint+8].size(); j++){
					double sipm_time_b = BRB_pmt_time[ievent][iscint+8][j];
					//cout << "h_BRB_T0_tof_timediff content: " << sipm_time_a - T0_time << " " << sipm_time_b - T0_time << endl;
					double position = (sipm_time_a - sipm_time_b - scint_bias) * 185 / 2.0 ;
					h_hits_centered->Fill(position, start_y[iscint]);
				}
			}
		}
	}
	TCanvas* c_hit_rec_centered = new TCanvas("c_hit_rec_centered", "c_hit_rec_centered", 1800, 900);
	h_hits_centered->Draw("colz");
	draw_boxes(scints);
	c_hit_rec_centered->Print("c_position_reconstruction_centered.pdf");

	TCanvas *c_T0_tof_time = new TCanvas("c_T0_tof_time", "c_T0_tof_time", 1800, 900);
	c_T0_tof_time->SetLogy();
	h_T0_tof_time->Draw("hist");

	c_T0_tof_time->Print("h_T0_tof_time.png");

	TCanvas *c_difftimes_pi_cut = new TCanvas("c_difftimes_pi_cut", "c_difftimes_pi_cut", 1800, 900);
	c_difftimes_pi_cut->Divide(4, 2);
	for (int i = 0; i < 8; i++){
		c_difftimes_pi_cut->cd(i+1);
		double x_max = h_difftime_pi_cut[i]->GetBinCenter(h_difftime_pi_cut[i]->GetMaximumBin());
		h_difftime_pi_cut[i]->Fit("gaus", "", "RS", x_max - 2, x_max + 2);
		TF1* fit_fun = h_difftime_pi_cut[i]->GetFunction("gaus");
		h_difftime_pi_cut[i]->Draw("hist");
		fit_fun->Draw("same");

		fit_pi_sigma[i] = fit_fun->GetParameter(2);
		fit_pi_sigma_error[i] = fit_fun->GetParError(2);
		fit_pi_mean[i] = fit_fun->GetParameter(1);
		//if (h_difftime_pi_cut[0]->GetEntries() < 100 || h_difftime_pi_cut[7]->GetEntries() < 100) continue;
		/* fit_pi_sigma[i] = h_difftime_pi_cut[i]->GetStdDev();
		   fit_pi_sigma_error[i] = h_difftime_pi_cut[i]->GetStdDevError();
		   fit_pi_mean[i] = h_difftime_pi_cut[i]->GetMean(); */
	}
	c_difftimes_pi_cut->Print("h_difftime_pi_cut.png");



	TCanvas *c_difftimes_mu_cut = new TCanvas("c_difftimes_mu_cut", "c_difftimes_mu_cut", 1800, 900);
	c_difftimes_mu_cut->Divide(4, 2);
	for (int i = 0; i < 8; i++){
		c_difftimes_mu_cut->cd(i+1);
		double x_max = h_difftime_mu_cut[i]->GetBinCenter(h_difftime_mu_cut[i]->GetMaximumBin());
		h_difftime_mu_cut[i]->Fit("gaus", "", "RS", x_max - 2, x_max + 2);
		TF1* fit_fun = h_difftime_mu_cut[i]->GetFunction("gaus");
		h_difftime_mu_cut[i]->Draw("hist");
		fit_fun->Draw("same");

		fit_mu_sigma[i] = fit_fun->GetParameter(2);
		fit_mu_sigma_error[i] = fit_fun->GetParError(2);
		fit_mu_mean[i] = fit_fun->GetParameter(1);
		//if (h_difftime_mu_cut[0]->GetEntries() < 100 || h_difftime_mu_cut[7]->GetEntries() < 100) continue;
		/* fit_mu_sigma[i] = h_difftime_mu_cut[i]->GetStdDev();
		   fit_mu_sigma_error[i] = h_difftime_mu_cut[i]->GetStdDevError();
		   fit_mu_mean[i] = h_difftime_mu_cut[i]->GetMean(); */
	}
	c_difftimes_mu_cut->Print("h_difftime_mu_cut.png");

	TCanvas *c_difftimes_e_cut = new TCanvas("c_difftimes_e_cut", "c_difftimes_e_cut", 1800, 900);
	c_difftimes_e_cut->Divide(4, 2);
	for (int i = 0; i < 8; i++){
		c_difftimes_e_cut->cd(i+1);
		double x_max = h_difftime_e_cut[i]->GetBinCenter(h_difftime_e_cut[i]->GetMaximumBin());
		h_difftime_e_cut[i]->Fit("gaus", "", "RS", x_max - 2, x_max + 2);
		TF1* fit_fun = h_difftime_e_cut[i]->GetFunction("gaus");
		h_difftime_e_cut[i]->Draw("hist");
		fit_fun->Draw("same");

		fit_e_sigma[i] = fit_fun->GetParameter(2);
		fit_e_sigma_error[i] = fit_fun->GetParError(2);
		fit_e_mean[i] = fit_fun->GetParameter(1);
		if (h_difftime_e_cut[0]->GetEntries() < 100 || h_difftime_e_cut[7]->GetEntries() < 100) continue;

		/* fit_e_sigma[i] = h_difftime_e_cut[i]->GetStdDev();
		   fit_e_sigma_error[i] = h_difftime_e_cut[i]->GetStdDevError();
		   fit_e_mean[i] = h_difftime_e_cut[i]->GetMean(); */
	}
	c_difftimes_e_cut->Print("h_difftime_e_cut.png");
	gStyle->SetOptStat(0);
	TCanvas *c_pi_sigmas = new TCanvas("c_pi_sigmas", "c_pi_sigmas", 1800, 900);
	for (int i = 0; i < 8; i++){
		h_pi_sigmas->SetBinContent(i+1, fit_pi_sigma[i]);
		h_pi_sigmas->SetBinError(i+1, fit_pi_sigma_error[i]);
		h_mu_sigmas->SetBinContent(i+1, fit_mu_sigma[i]);
		h_mu_sigmas->SetBinError(i+1, fit_mu_sigma_error[i]);
		h_e_sigmas->SetBinContent(i+1, fit_e_sigma[i]);
		h_e_sigmas->SetBinError(i+1, fit_e_sigma_error[i]);
	}
	h_pi_sigmas->GetYaxis()->SetRangeUser(0, 0.65);
	h_pi_sigmas->Draw("histe");
	h_mu_sigmas->SetLineColor(kRed);
	h_e_sigmas->SetLineColor(kBlack);
	h_mu_sigmas->Draw("histesame");
	h_e_sigmas->Draw("histesame");
	TLegend *leg_comp = new TLegend(0.7, 0.8, 0.9, 0.9);
	leg_comp->AddEntry(h_e_sigmas, "Sigmas in electron cuts", "le");
	leg_comp->AddEntry(h_mu_sigmas, "Sigmas in muon cuts", "le");
	leg_comp->AddEntry(h_pi_sigmas, "Sigmas in pion cuts", "le");
	leg_comp->Draw("same");

	c_pi_sigmas->Print("h_sigmas_compare.png");

	for (int ievent = 0; ievent < beamline_ACT_charges.size(); ievent++){
		double ACT_g1_charge = 0;
		double ACT_g2_charge = 0;
		if (TDCT0_hits[ievent].size() < 1) continue;
		double avg_time_T0 = 0;
		int npmt_hit_T0 = 0;

		for (int i = 0; i < beamline_T0_time[ievent].size(); i++){
			for (int j = 0; j < beamline_T0_time[ievent][i].size(); j++){
				avg_time_T0 += beamline_T0_time[ievent][i][j];
				npmt_hit_T0++;
			}
		}

		avg_time_T0 /= npmt_hit_T0;

		double T0_time = avg_time_T0 - TDCT0_hits[ievent][0];

		for (int i = 0; i < 6; i++){
			for (int icharge = 0; icharge < beamline_ACT_charges[ievent][i].size(); icharge++){
				ACT_g1_charge+=beamline_ACT_charges[ievent][i][icharge];
			}
		}
		for (int i = 6; i < 12; i++){
			for (int icharge = 0; icharge < beamline_ACT_charges[ievent][i].size(); icharge++){
				ACT_g2_charge+=beamline_ACT_charges[ievent][i][icharge];
			}
		}

		for (int iscint = 0; iscint < 8; iscint++){
			double avg_SiPM_time;
			for (int i = 0; i < BRB_pmt_time[ievent][iscint].size(); i++){
				double sipm_time_a = BRB_pmt_time[ievent][iscint][i];
				double sipm_charge_a = BRB_pmt_charge[ievent][iscint][i];
				bool pair_hit = true;
				if (sipm_time_a < -165 || sipm_time_a > -145) continue;
				for (int j = 0; j < BRB_pmt_time[ievent][iscint+8].size(); j++){
					double sipm_time_b = BRB_pmt_time[ievent][iscint+8][j];
					double sipm_charge_b = BRB_pmt_charge[ievent][iscint+8][j];
					if (sipm_time_b < -165 || sipm_time_b > -145) continue;
				
					if (pair_hit) {
						h_timediff_charge_total->Fill(sipm_charge_a + sipm_charge_b, sipm_time_a - sipm_time_b);
						h_timediff_charge[iscint]->Fill(sipm_charge_a + sipm_charge_b, sipm_time_a - sipm_time_b);
						double tof_time = (sipm_time_a + sipm_time_b)/2 - T0_time;
						if (tof_time > 45) continue; // cut out protons
						if (ACT_g1_charge < ACT_g1_cut && ACT_g2_charge < ACT_g2_cut && ACT_g2_charge < 3000){ //pions
							h_pi_charges->Fill(sipm_charge_a + sipm_charge_b);
							h_time_charge[iscint]->Fill(sipm_charge_a + sipm_charge_b, fit_pi_sigma[iscint]);
							h_time_charge_total->Fill(sipm_charge_a + sipm_charge_b, fit_pi_sigma[iscint]);
							h_timediff_charge_picut[iscint]->Fill(sipm_charge_a + sipm_charge_b, sipm_time_a - sipm_time_b);
							h_timediff_charge_total_picut->Fill(sipm_charge_a + sipm_charge_b, sipm_time_a - sipm_time_b);

						}
						if (ACT_g1_charge < ACT_g1_cut && ACT_g2_charge > ACT_g2_cut && ACT_g2_charge < 3000){ //muons
							h_mu_charges->Fill(sipm_charge_a + sipm_charge_b);
							h_time_charge[iscint]->Fill(sipm_charge_a + sipm_charge_b, fit_mu_sigma[iscint]);
							h_time_charge_total->Fill(sipm_charge_a + sipm_charge_b, fit_mu_sigma[iscint]);
							h_timediff_charge_mucut[iscint]->Fill(sipm_charge_a + sipm_charge_b, sipm_time_a - sipm_time_b);
							h_timediff_charge_total_mucut->Fill(sipm_charge_a + sipm_charge_b, sipm_time_a - sipm_time_b);
						}
						if (ACT_g1_charge > ACT_g1_cut){ //electrons
							h_e_charges->Fill(sipm_charge_a + sipm_charge_b);
							h_time_charge[iscint]->Fill(sipm_charge_a + sipm_charge_b, fit_e_sigma[iscint]);
							h_time_charge_total->Fill(sipm_charge_a + sipm_charge_b, fit_e_sigma[iscint]);
							h_timediff_charge_ecut[iscint]->Fill(sipm_charge_a + sipm_charge_b, sipm_time_a - sipm_time_b);
							h_timediff_charge_total_ecut->Fill(sipm_charge_a + sipm_charge_b, sipm_time_a - sipm_time_b);

						}
					}
					pair_hit=false;
				}
			}
		}
	}
	TCanvas *c_time_charge = new TCanvas("c_time_charge", "c_time_charge", 1800, 900);
	c_time_charge->Divide(4,2);
	for (int i = 0; i < 8; i++){
		c_time_charge->cd(i+1);
		h_time_charge[i]->Draw("colz");
	}
	c_time_charge->Print("h_time_charge.pdf");
	c_time_charge->Print("h_time_charge.png");

	TCanvas *c_time_charge_total = new TCanvas("c_time_charge_total", "c_time_charge_total", 1800, 900);
	h_time_charge_total->Draw("colz");
	c_time_charge_total->Print("h_time_charge_total.png");

	TCanvas *c_timediff_charge = new TCanvas("c_timediff_charge", "c_timediff_charge", 1800, 900);
	c_timediff_charge->Divide(4,2);
	for (int i = 0; i < 8; i++){
		c_timediff_charge->cd(i+1);
		h_timediff_charge[i]->Draw("colz");
	}
	c_timediff_charge->Print("h_timediff_charge.pdf");
	c_timediff_charge->Print("h_timediff_charge.png");

	TCanvas *c_timediff_charge_total = new TCanvas("c_timediff_charge_total", "c_timediff_charge_total", 1800, 900);
	h_timediff_charge_total->Draw("colz");
	c_timediff_charge_total->Print("h_timediff_charge_total.png");
	TCanvas *c_timediff_charge_total_ecut = new TCanvas("c_timediff_charge_total_ecut", "c_timediff_charge_total_ecut", 1800, 900);
	h_timediff_charge_total_ecut->Draw("colz");
	c_timediff_charge_total_ecut->Print("h_timediff_charge_total_ecut.png");
	TCanvas *c_timediff_charge_total_mucut = new TCanvas("c_timediff_charge_total_mucut", "c_timediff_charge_total_mucut", 1800, 900);
	h_timediff_charge_total_mucut->Draw("colz");
	c_timediff_charge_total_mucut->Print("h_timediff_charge_total_mucut.png");
	TCanvas *c_timediff_charge_total_picut = new TCanvas("c_timediff_charge_total_picut", "c_timediff_charge_total_picut", 1800, 900);
	h_timediff_charge_total_picut->Draw("colz");
	c_timediff_charge_total_picut->Print("h_timediff_charge_total_picut.png");

	TCanvas *c_timediff_charge_ecut = new TCanvas("c_timediff_charge_ecut", "c_timediff_charge_ecut", 1800, 900);
	c_timediff_charge_ecut->Divide(4,2);
	for (int i = 0; i < 8; i++){
		c_timediff_charge_ecut->cd(i+1);
		h_timediff_charge_ecut[i]->Draw("colz");
	}
	c_timediff_charge_ecut->Print("h_timediff_charge.pdf");
	c_timediff_charge_ecut->Print("h_timediff_charge.png");

	TCanvas *c_timediff_charge_mucut = new TCanvas("c_timediff_charge_mucut", "c_timediff_charge_mucut", 1800, 900);
	c_timediff_charge_mucut->Divide(4,2);
	for (int i = 0; i < 8; i++){
		c_timediff_charge_mucut->cd(i+1);
		h_timediff_charge_mucut[i]->Draw("colz");
	}
	c_timediff_charge_mucut->Print("h_timediff_charge_mucut.pdf");
	c_timediff_charge_mucut->Print("h_timediff_charge_mucut.png");
	
	TCanvas *c_timediff_charge_picut = new TCanvas("c_timediff_charge_picut", "c_timediff_charge_picut", 1800, 900);
	c_timediff_charge_picut->Divide(4,2);
	for (int i = 0; i < 8; i++){
		c_timediff_charge_picut->cd(i+1);
		h_timediff_charge_picut[i]->Draw("colz");
	}
	c_timediff_charge_picut->Print("h_timediff_charge_picut.pdf");
	c_timediff_charge_picut->Print("h_timediff_charge_picut.png");

	TCanvas *c_charge_picut = new TCanvas("c_charge_picut", "c_charge_picut", 1800, 900);
	h_pi_charges->Draw("hist");
	c_charge_picut->Print("h_charge_picut.pdf");
	c_charge_picut->Print("h_charge_picut.png");

	TCanvas *c_charge_mucut = new TCanvas("c_charge_mucut", "c_charge_mucut", 1800, 900);
	h_mu_charges->Draw("hist");
	c_charge_mucut->Print("h_charge_mucut.pdf");
	c_charge_mucut->Print("h_charge_mucut.png");

	TCanvas *c_charge_ecut = new TCanvas("c_charge_ecut", "c_charge_ecut", 1800, 900);
	h_e_charges->Draw("hist");
	c_charge_ecut->Print("h_charge_ecut.pdf");
	c_charge_ecut->Print("h_charge_ecut.png");
	// ################# end drawing histograms ###############################

	gSystem->cd("../../..");
	std::ofstream out(out_file, std::ios::app);
	for (int i = 0; i < fit_BRB_sigma.size(); i++){
		if (energy == 0) out << run_number <<"\t" << fit_BRB_sigma[i] << "\t" << fit_BRB_sigma_error[i] << "\t" << scint_dimensions[i][0] <<"\n";
		else out << run_number << "\t" << energy << "\t" << fit_BRB_sigma[i] << "\t" << fit_BRB_sigma_error[i] << "\t" << scint_dimensions[i][0] <<"\n";
	}
	out.close();

	// ####################         PARTICLE TOF DETECTOR SIGMA COMPARTISON START         ##################################




	// XXXXXXXXXXXXXXXXXXXXXXX FINISH COMPARING SIGMAS OF DIFFERENT PARITCLES XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX




	// ###################### OLD GRAPHS, HISTOGRAMS ##############################

	/*
	   TCanvas* c_bm_trig_0 = new TCanvas("", "", 1800, 900);
	   h_beamline_trigger_time_0->Draw("hist");

	   TCanvas* c_bm_trig_1 = new TCanvas("", "", 1800, 900);
	   h_beamline_trigger_time_1->Draw("hist");

	   TCanvas* c_bm_tof_times = new TCanvas("", "", 1800, 900);
	   c_bm_tof_times->Divide(4, 4);
	   for (int i = 0; i < 16; i++){
	   c_bm_tof_times->cd(i+1);
	   h_beamline_tof_times[i]->Draw("hist");
	   }

	   TCanvas* c_brb_tof_times = new TCanvas("", "", 1800, 900);
	   c_brb_tof_times->Divide(4, 4);
	   for (int i = 0; i < 16; i++){
	   c_brb_tof_times->cd(i+1);
	   h_BRB_tof_times[i]->Draw("hist");
	   }

	   TCanvas* c_T0_avg_times = new TCanvas("", "", 900, 900);
	   h_beamline_T0_avgtime->Draw("hist");

	   TCanvas* c_T0_counts = new TCanvas("", "", 900, 900);
	   h_beamline_T0_counts->Draw("hist");
	   c_T0_counts->Print("h_T0_counts.pdf");

	   TCanvas* c_T0_times = new TCanvas("", "", 1800, 900);
	   c_T0_times->Divide(2,2);
	   for (int i = 0; i < 4; i++){
	   c_T0_times->cd(i+1);
	   h_beamline_T0_times[i]->Draw("hist");
	   }
	   c_T0_times->Print("h_T0_times.pdf");

	   TCanvas* c_T1_times = new TCanvas("", "", 1800, 900);
	   c_T1_times->Divide(2,2);
	   for (int i = 0; i < 4; i++){
	   c_T1_times->cd(i+1);
	   h_beamline_T1_times[i]->Draw("hist");
	   }
	   c_T1_times->Print("h_T1_times.pdf");

	   TCanvas* c_T0_tof_difftimes = new TCanvas("", "", 1800, 900);
	   c_T0_tof_difftimes->Divide(4,2);
	   for (int i = 0; i < 8; i++){
	   c_T0_tof_difftimes->cd(i+1);
	   h_beamline_T0_tof_timediff[i]->Draw("colz");
	   }
	   c_T0_tof_difftimes->Print("h_T0_tof_difftimes.pdf");

	   TCanvas* c_T1_tof_difftimes = new TCanvas("", "", 1800, 900);
	   c_T1_tof_difftimes->Divide(4,2);
	   for (int i = 0; i < 8; i++){
	   c_T1_tof_difftimes->cd(i+1);
	   h_beamline_T1_tof_timediff[i]->Draw("colz");
	   }
	   c_T1_tof_difftimes->Print("h_T1_tof_difftimes.pdf");

	   TCanvas* c_charges = new TCanvas("c_charges", "c_charges", 1800, 900);
	   c_charges->Divide(4, 4);
	   for (int i = 0; i < 16; i++){
	   c_charges->cd(i+1);
	   h_charge[i]->Draw("hist");
	   }

	   c_charges->Print("h_charges.pdf");

	   TCanvas* c_BRB_tof_avg_diffcharges = new TCanvas("c_BRB_tof_diffcharges", "c_BRB_tof_diffcharges", 1800, 900);
	   c_BRB_tof_avg_diffcharges->Divide(4, 2);
	for (int i = 0; i < 8; i++){
		c_BRB_tof_avg_diffcharges->cd(i+1);
		h_BRB_tof_avg_chargediff[i]->Draw("colz");
	}
	c_BRB_tof_avg_diffcharges->Print("h_BRB_tof_diffcharges.pdf");

	TCanvas* c_BRB_q_diff= new TCanvas("", "", 1800, 900);
	c_BRB_q_diff->Divide(4, 2);
	for (int i = 0; i < 8; i++){
		c_BRB_q_diff->cd(i+1);
		h_BRB_q_diff[i]->Draw("colz");
	}
	c_BRB_q_diff->Print("h_BRB_tof_diffcharges.pdf");

	TCanvas* c_BRB_chargediff= new TCanvas("", "", 1800, 900);
	c_BRB_chargediff->Divide(4, 2);
	for (int i = 0; i < 8; i++){
		c_BRB_chargediff->cd(i+1);
		h_BRB_tof_chargediff[i]->Draw("colz");
	}
	c_BRB_chargediff->Print("h_BRB_tof_chargediffs.pdf");


	TCanvas* c_BRB_total_avg_time_diff= new TCanvas("", "", 1800, 900);
	c_BRB_total_avg_time_diff->Divide(4, 2);
	for (int i = 0; i < 8; i++){
		c_BRB_total_avg_time_diff->cd(i+1);
		h_BRB_total_avg_time_diff[i]->Draw("colz");
	}
	c_BRB_total_avg_time_diff->Print("h_BRB_tof_total_avg_time_diffs.pdf");



	double e_time = 6.5/3e8 * std::sqrt(1.0 + 0.511/260.0) * pow(10, 9);
	cout << "e_time is " << e_time << endl;
	double u_time= 6.5/3e8 * std::sqrt(1.0 + 105.7/260.0) * pow(10, 9);
	cout << "u_time is " << u_time << endl;
	double pi_time = 6.5/3e8 * std::sqrt(1.0 + 139.57/260.0) * pow(10, 9);
	cout << "pi_time is " << pi_time << endl;
	TCanvas* c_T0_tof_times = new TCanvas("", "", 1800, 900);
	c_T0_tof_times->Divide(4,4);
	for (int i = 0; i < 16; i++){
		c_T0_tof_times->cd(i+1);
		h_beamline_T0_tof_time[i]->Draw("hist");
		double y_max =  h_beamline_T0_tof_time[i]->GetMaximum();
		TLine* line_e = new TLine(e_time, 0, e_time, y_max);
		line_e->SetLineColor(kBlue);
		TLine* line_u = new TLine(u_time, 0 , u_time, y_max);
		line_u->SetLineColor(kRed);
		TLine* line_pi = new TLine(pi_time, 0, pi_time, y_max);
		line_pi->SetLineColor(kGreen);
		line_e->Draw("same");
		line_u->Draw("same");
		line_pi->Draw("same");
	}
	c_T0_tof_times->Print("h_T0_tof_times.pdf");

	cout << "pi_time is " << pi_time << endl;
	TCanvas* c_T1_tof_times = new TCanvas("", "", 1800, 900);
	c_T1_tof_times->Divide(4,4);
	for (int i = 0; i < 16; i++){
		c_T1_tof_times->cd(i+1);
		h_beamline_T1_tof_time[i]->Draw("hist");
	}
	c_T1_tof_times->Print("h_T1_tof_times.pdf");

	TCanvas* c_T0_T1_times = new TCanvas("", "", 900, 900);
	h_beamline_T0_T1_time->Draw("hist");
	c_T0_T1_times->Print("h_T0_T1_time.pdf");

	TCanvas* c_BRB_T0_tof_times = new TCanvas("c_BRB_T0_tof_times", "c_BRB_T0_tof_times", 1800, 900);
	c_BRB_T0_tof_times->Divide(4, 4);
	for (int i = 0; i < 16; i++){
		c_BRB_T0_tof_times->cd(i+1);
		h_BRB_T0_tof_time[i]->Draw("hist");
	}

	c_BRB_T0_tof_times->Print("h_BRB_T0_tof_times.pdf");


	TCanvas* c_BRB_tof_avg_difftimes = new TCanvas("c_BRB_tof_difftimes", "c_BRB_tof_difftimes", 1800, 900);
	c_BRB_tof_avg_difftimes->Divide(4, 2);
	for (int i = 0; i < 8; i++){
		c_BRB_tof_avg_difftimes->cd(i+1);
		h_BRB_tof_avg_timediff[i]->Draw("colz");
	}
	c_BRB_tof_avg_difftimes->Print("h_BRB_tof_difftimes.pdf");


	TCanvas* c_BRB_avg_tof_tof_times = new TCanvas("c_BRB_avg_tof_tof_times", "c_BRB_avg_tof_tof_times", 1800, 900);
	c_BRB_avg_tof_tof_times->Divide(4, 4);
	for (int i = 0; i < 16; i++){
		c_BRB_avg_tof_tof_times->cd(i+1);
		h_BRB_avg_tof_tof_time[i]->Draw("hist");
	}

	c_BRB_avg_tof_tof_times->Print("h_BRB_avg_tof_tof_times.pdf");
	*/


} //end of
