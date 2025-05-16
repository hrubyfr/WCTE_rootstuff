#include "math.h"

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

#include "TTree.h"
#include <TVirtualPad.h>

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

//################################################ MAIN FUNCTION ###########################################


void read_matched_data(TString filename = "/media/frantisek/T7\ Shield/WCTE_data/WCTE_offline_R1361S0.root"){

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

	TString plots_folder = Form("analysis_plots/matching/file_%i", run_number);
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
				50, 0, 100, 50, 0, 100);
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


		/*
		   for (int i = 0; i < hit_mpmt_card_ids->size(); i++){
		   if(hit_mpmt_card_ids->at(i) == 132 && (tof_pmt_ids.find(hit_pmt_channel_ids->at(i)) != tof_pmt_ids.end())){
		   double diff_time = hit_pmt_times->at(i) - trigger_time[2];
		   tof_minus_trigger->Fill(diff_time);
		   double T0_tof_time = diff_time - T0_avg_time;
		   double T1_tof_time = diff_time - T1_avg_time;
		   T0_TOF_times->Fill(T0_tof_time);
		   T1_TOF_times->Fill(T1_tof_time);
		   }
		   }

		   for(int i = 0; i < hit_mpmt_card_ids->size(); i++){
		   int pmt_id = hit_pmt_channel_ids->at(i);
		   int mpmt_id = hit_mpmt_card_ids->at(i);
		   if (beam_cards.count(mpmt_id)){
		   if(mpmt_id == 130){
		   card1_hist->Fill(pmt_id);
		   } //do things for mpmt card 130
		   else if(mpmt_id == 131){
		   card2_hist->Fill(pmt_id);
		   }//do things for mpmt card 131
		   else if(hit_mpmt_card_ids->at(i) == 132){
		   card3_hist->Fill(pmt_id);


		   pmt_charge[hit_pmt_channel_ids->at(i)].push_back(hit_pmt_charges->at(i));
		   pmt_time[hit_pmt_channel_ids->at(i)].push_back(hit_pmt_times->at(i) - trigger_time[2]);
		   if (pmt_charge[hit_pmt_channel_ids->at(i)].size() != pmt_time[hit_pmt_channel_ids->at(i)].size()) cout << "Charge and times sizes are not equal!" << endl;
		   }//do things for mpmt card 132
		   } //end if over cards 130 through 132
		   } //end loop over mpmt card numbers
		     // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

		     // ################################## WAVEFORM ANALYSIS ############################# 
		     double waveform_trigger;

		     int N_trigs = 0;
		     for (int i = 0; i < pmt_waveform_times->size(); i++){
		     int pmt = pmt_waveform_pmt_channel_ids->at(i);
		     if (pmt_waveform_mpmt_card_ids->at(i) == 132 && pmt_waveform_pmt_channel_ids->at(i) == 19){
		     waveform_trigger = pmt_waveform_times->at(i);
		     h_132_trigger_wftime->Fill(pmt_waveform_times->at(i));
		     N_trigs++;
		     }
		     }
		     if(N_trigs < 1) continue;

		     for (int i = 0; i < pmt_waveforms->size(); i++){
		     int mpmt_id = pmt_waveform_mpmt_card_ids->at(i);
		     int pmt_id = pmt_waveform_pmt_channel_ids->at(i);
		     double waveform_time = pmt_waveform_times->at(i);
		     if (mpmt_id == 132 && tof_pmt_ids.count(pmt_id)){
		     int hist_id = std::distance(tof_pmt_ids.begin(), tof_pmt_ids.find(pmt_id));
		     double diff_time = waveform_time - waveform_trigger;

		     tof_waveform_times[hist_id]->Fill(pmt_waveform_times->at(i) - waveform_trigger);
		     if (diff_time > -155 && diff_time < -150){
		     auto& waveform = pmt_waveforms->at(i);
		     for (int itime = 0; itime < waveform.size(); itime++){
		     tof_heatmaps[hist_id]->Fill(itime, waveform.at(itime));
		     }
		     }
		     }
		     } // end loop over waveforms
		       // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
		       */
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

	for (int event = 0; event < BRB_pmt_time.size(); event++){
		for (int iscint = 0; iscint < 8; iscint++){
			for (int i = 0; i < BRB_pmt_time[event][iscint].size(); i++){
				double sipm_time_a = BRB_pmt_time[event][iscint][i];
				if (sipm_time_a < -165 || sipm_time_a > -145) continue;
				for (int j = 0; j < BRB_pmt_time[event][iscint+8].size(); j++){
					double sipm_time_b = BRB_pmt_time[event][iscint+8][j];
					if (sipm_time_b < -165 || sipm_time_b > -145) continue;
					h_BRB_diff_time[iscint]->Fill(sipm_time_a - sipm_time_b);
				}
			}
		}
	}
	for (int ievent = 0; ievent < BRB_pmt_time.size(); ievent++){
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
	}
	cout << "Analysis done, now drawing histograms" << endl;

	// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	gStyle->SetPalette(1);
	gStyle->SetOptFit(1);
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


	TCanvas* c_bm_diff_times = new TCanvas("", "", 1800, 900);
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
	c_bm_diff_times->Close();

	TCanvas* c_brb_diff_times = new TCanvas("", "", 1800, 900);
	c_brb_diff_times->Divide(4, 2);
	for (int i = 0; i < 8; i++){
		c_brb_diff_times->cd(i+1);
		double x_max = h_BRB_diff_time[i]->GetBinCenter(h_beamline_diff_time[i]->GetMaximumBin());
		h_BRB_diff_time[i]->Fit("gaus", "", "RS", x_max - 2, x_max + 2);
		TF1* fit_fun = h_BRB_diff_time [i]->GetFunction("gaus");
		h_BRB_diff_time[i]->Draw("hist");
		fit_fun->Draw("same");
	}
	c_brb_diff_times->Print("h_brb_diff_times.pdf");
	c_brb_diff_times->Close();

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

	TCanvas* c_BRB_T0_tof_difftimes = new TCanvas("c_BRB_T0_tof_difftimes", "c_BRB_T0_tof_difftimes", 1800, 900);
	c_BRB_T0_tof_difftimes->Divide(4, 2);
	for (int i = 0; i < 8; i++){
		c_BRB_T0_tof_difftimes->cd(i+1);
		h_BRB_T0_tof_timediff[i]->Draw("colz");
	}
	c_BRB_T0_tof_difftimes->Print("h_BRB_T0_tof_difftimes.pdf");

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

} //end of code
