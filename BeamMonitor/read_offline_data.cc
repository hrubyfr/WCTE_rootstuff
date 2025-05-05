#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "math.h"


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
#include <algorithm>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <ios>
#include <ostream>
#include <vector>
#include <array>
#include <set>


using std::array;
using std::string;
using std::vector;
using std::cout;
using std::endl;
using ROOT::Math::XYVector;
using std::set;

XYVector construct_SiPM(double x, double y){
	return XYVector(x, y);
}

double distance(double x, double y){
	return sqrt(pow(x, 2) + pow(y, 2));
}



vector<vector<XYVector>> sort_points(vector<XYVector> points, vector<array<double, 2>> scint_dimensions, vector<double> start_y){
	vector<vector<XYVector>> sorted_points(9);
	for (int i = 0; i < points.size(); i++){
		if (abs(points[i].X()) < abs(scint_dimensions[0][0]/2) && points[i].Y() < start_y[0] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[0] - scint_dimensions[0][1]/2) sorted_points[0].push_back(points[i]);
		else if (abs(points[i].X()) < abs(scint_dimensions[1][0]/2) && points[i].Y() < start_y[1] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[1] - scint_dimensions[0][1]/2) sorted_points[1].push_back(points[i]);
		else if (abs(points[i].X()) < abs(scint_dimensions[2][0]/2) && points[i].Y() < start_y[2] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[2] - scint_dimensions[0][1]/2) sorted_points[2].push_back(points[i]);
		else if (abs(points[i].X()) < abs(scint_dimensions[3][0]/2) && points[i].Y() < start_y[3] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[3] - scint_dimensions[0][1]/2) sorted_points[3].push_back(points[i]);
		else if (abs(points[i].X()) < abs(scint_dimensions[4][0]/2) && points[i].Y() < start_y[4] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[4] - scint_dimensions[0][1]/2) sorted_points[4].push_back(points[i]);
		else if (abs(points[i].X()) < abs(scint_dimensions[5][0]/2) && points[i].Y() < start_y[5] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[5] - scint_dimensions[0][1]/2) sorted_points[5].push_back(points[i]);
		else if (abs(points[i].X()) < abs(scint_dimensions[6][0]/2) && points[i].Y() < start_y[6] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[6] - scint_dimensions[0][1]/2) sorted_points[6].push_back(points[i]);
		else if (abs(points[i].X()) < abs(scint_dimensions[7][0]/2) && points[i].Y() < start_y[7] + scint_dimensions[0][1]/2 && points[i].Y() > start_y[7] - scint_dimensions[0][1]/2) sorted_points[7].push_back(points[i]);
		else sorted_points[8].push_back(points[i]);
	}//end for loop over points
	return sorted_points;
}//end sort point function



double Det_response(XYVector SiPM_position, XYVector position){
	XYVector Det_position(SiPM_position.X(), SiPM_position.Y());
	auto rel_position = Det_position.X() - position.X();
	// std::cout << "Relative position of point to detector is: " << rel_position << endl;
	// std::cout << "Distance of point from detector is " << distance(rel_position.X(), rel_position.Y()) << " cm" << endl;
	double r = rel_position;
	double lambda = 1200.0;
	double atten = exp(-abs(rel_position)/lambda);
	auto response = atten / pow(r, 2);
	return response;
}

void draw_boxes(vector<TBox*> boxes){
	for(int i = 0; i < boxes.size(); i++){
		boxes[i]->SetFillStyle(0);
		boxes[i]->SetLineColor(kBlack);
		boxes[i]->Draw("same");
	}
}


double GetMinimum_chi2(const int niter, const double tolerance, const float charge1, const float charge2, const XYVector SiPM1, const XYVector SiPM2){
	int iter = 0;
	double x = (SiPM1.X() * charge1 + SiPM2.X() * charge2)/(charge1 + charge2); //Start point guess from weighted charge		
	double step_size = 10.0;
	double chi2 = 0;
	double chi2_history = 0;
	auto evaluate = [&x, SiPM1, SiPM2, charge1, charge2, &chi2](){
		double dx1 = x - SiPM1.X();
		double dx2 = x - SiPM2.X();

		double theory_sig1 = 1.0 / (dx1 * dx1);
		double theory_sig2 = 1.0 / (dx2 * dx2);

		double residual1, residual2;					//Get chi2 at start point
		residual1 = charge1 - theory_sig1;
		residual2 = charge2 - theory_sig2;

		chi2 = residual1 * residual1 + residual2 * residual2;
		return chi2;
	};
	evaluate();
	int i_step = 0;
	int direction;
	if (charge1 > charge2) direction = -1;
	else direction = 1;
	while (chi2 > tolerance){
		x += step_size * direction;
		chi2_history = chi2;
		chi2 = evaluate();


		if (chi2 > chi2_history) {
			direction *=-1;
			step_size /= 10.0;
		}
		if (abs(x) > abs(SiPM1.X())) {
			std::cerr << "ERROR: X OUT OF BOUNDS" << endl;
			return -1;
		}
		if(iter == niter) break;
		iter++;
	}
	return x;
}

//################################################ MAIN FUNCTION ###########################################


void read_offline_data(TString filename = "data/offline_data/WCTE_offline_R1308S0.root"){

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

	TString plots_folder = Form("analysis_plots/file_%i", run_number);
	if(gSystem->AccessPathName(plots_folder)) gSystem->Exec("mkdir " + plots_folder);  // make a folder to save plots to
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
	set<int> T0_beamline_pmt_ids = {0, 1, 2, 3};
	set<int> T1_beamline_pmt_ids = {4, 5, 6, 7};

	// ################ Initialise histograms #########################

	TH1D* card1_hist  = new TH1D("card1_hist", "Histogram of hit pmts on card 130; pmt_id; counts", 20, 0, 20); 
	TH1D* card2_hist  = new TH1D("card2_hist", "Histogram of hit pmts on card 131; pmt_id; counts", 20, 0, 20); 
	TH1D* card3_hist  = new TH1D("card3_hist", "Histogram of hit pmts on card 132; pmt_id; counts", 20, 0, 20); 

	vector<TH1D*> tof_diff(8);
	vector<TH2D*> tof_v_tof(8);
	vector<TH1D*> tof_q_diff(8);
	vector<TH2D*> tof_qvq(8);
	vector<TH2D*> tof_heatmaps(16);
	vector<TH1D*> tof_waveform_times(16);
	vector<TH1D*> tof_charges(16);
	vector<TH1D*> tof_times(16);

	TH1D* tof_minus_trigger = new TH1D("tof_minus_trigger", "Tof times - Trigger time; time; count", 500, -400, 0);
	TH1D* BM_T0_times = new TH1D("BM_T0_times", "T0 detection time; time; count", 500, -500, 500);
	TH1D* BM_T1_times = new TH1D("BM_T1_times", "T1 detection time; time; count", 500, -500, 500);
	TH1D* T0_times = new TH1D("T0_times", "T0 detection time; time; count", 500, -500, 500);
	TH1D* T1_times = new TH1D("T1_times", "T1 detection time; time; count", 500, -500, 500);

	for (int i = 0; i < tof_diff.size(); i++){
		tof_diff[i] = new TH1D(Form("tof_diff_%i", i), Form("Run %i diferrence of time hits in TOF-%i - TOF-%i; time TOF-%i - TOF %i [ns];counts", 
					run_number, i, i+8, i, i+8), 1000, -100, 100);
		tof_v_tof[i] = new TH2D(Form("tof_v_tof_%i", i), Form("Run %i time hits in opposite SiPMs - TOF-%i and TOF-%i; time TOF-%i; time TOF %i; counts", 
					run_number, i, i+8, i, i+8), 100, -500, 500, 100, -500, 500);
		tof_q_diff[i] = new TH1D(Form("tof_q_diff_%i", i), Form("Run %i diferrence of charge hits in TOF-%i - TOF-%i; Charge TOF-%i - TOF %i ; counts", 
					run_number, i, i+8, i, i+8), 1000, -8000, 8000);
		tof_qvq[i] = new TH2D(Form("tof_qvq_%i", i), Form("Run %i charge hits in opposite SiPMs - TOF-%i and TOF-%i; charge TOF-%i; charge TOF %i; counts", 
					run_number, i, i+8, i, i+8), 100, 0, 8000, 100, 0, 8000);
	}


	TH1D* T0_TOF_times = new TH1D("T0_TOF_times", Form("Run %i T0 to TOF flight times; T0-TOF time; counts", run_number),500, -10000, 10000);
	TH1D* T1_TOF_times = new TH1D("T1_TOF_times", Form("Run %i T1 to TOF flight times; T1-TOF time; counts", run_number),500, -500, 500);
	
	TH1D* h_132_trigger_wftime = new TH1D("h_132_trigger_wftime", Form("Run %i board 132 waveform trigger time; trigger time; counts", run_number), 500, -500, 4000);


	for (int i = 0; i < 16; i++){
		tof_times[i] = new TH1D(Form("tof_times_%i", i), Form("Run %i TOF-%i time", run_number, i), 300, -500, 100);
		tof_charges[i] = new TH1D(Form("tof_charges_%i", i), Form("Run %i TOF-%i charge", run_number, i), 500, 0, 8000);
		tof_heatmaps[i] = new TH2D(Form("tof_heatmaps_%i", i), Form("Run %i TOF-%i heatmap", run_number, i), 128, 0, 128, 100, 0, 1000);
		tof_waveform_times[i] = new TH1D(Form("tof_waveform_times_%i", i), Form("Run %i TOF-%i Waveform times", run_number, i), 300, -500, 0);

	}

	int verb = 1000;


	std::map<int, bool> Trigger_hits;

	for (int i = 0; i < 8; i++){
		Trigger_hits[i] = false;		//Fill map with detctor IDs we want to do cuts on (0, 1, 2, 3 - T0; 4, 5, 6, 7 - T1; 43, 44 - T4)
	}
	Trigger_hits[42] = false;
	Trigger_hits[43] = false;

	std::map<int, bool> HC_hits;
	HC_hits[9] = false;				//Fill map with HC IDs - we want to exclude all event when HCs were hit (9 - HC-0; 10 - HC-1)
	HC_hits[10] = false;

	std::set<int> T0_IDs = {12, 13, 14, 15};
	std::set<int> T1_IDs = {13, 14, 15, 16};
	vector<vector<float>> pmt_charge(20);
	vector<vector<float>> pmt_time(20);
	std::map<int, int> corresponding_pmt = {{0, 8}, {1, 10},{2, 11},{3, 12}, {4, 13},{5, 14},{6, 15},{7, 16},}; 		//Map of opposite SiPMs
	std::map<int, int> inverted_pmt = {{8, 0}, {10, 1},{11, 2},{12, 3}, {13, 4},{14, 5},{15, 6},{16, 7},}; 		//Map of opposite SiPMs
	for(long ievent = 0; ievent < nentries; ievent++){
		tree->GetEntry(ievent);

		// ################## console print out progress ########################

		if (ievent % verb == 0){
			printf("\rProcessed %ld / %ld events (%.1f%%)", ievent, nentries, 100.0 * ievent / nentries);		//output progress on console
			fflush(stdout);
		}

		// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



		/*
		   cout << "Beamline PMT qcd IDs ("; 
		   for (int i = 0; i < beamline_pmt_qdc_ids->size(); i++){
		   cout << beamline_pmt_qdc_ids->at(i) << ", ";
		   }
		   cout << ")" << endl;

		   cout << "Beamline PMT qcd charges: ("; 
		   for (int i = 0; i < beamline_pmt_qdc_charge->size(); i++){
		   cout << beamline_pmt_qdc_charge->at(i) << ", ";
		   }
		   cout << ")" << endl;
		   */
		//cout << "Trigger_times size: " <<  trigger_times->size() << ", hit_mpmt_card_ids size: " << hit_mpmt_card_ids->size() << endl;


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




		if (pass_cut){
			// Set trigger times for the beam monitors cards

			double trigger_time[3];
			double T0_diff_time;
			int card_ID;
			int pmt_ID;
			int T0_hits = 0;
			int T1_hits = 0;
			double T0_avg_time;
			double T1_avg_time;
			for (int i = 0; i < hit_mpmt_card_ids->size(); i++){
				card_ID = hit_mpmt_card_ids->at(i);
				pmt_ID = hit_pmt_channel_ids->at(i);
				if(card_ID == 130 && pmt_ID == 19) trigger_time[0] = hit_pmt_times->at(i);
				else if(card_ID == 131 && pmt_ID == 0) trigger_time[1] = hit_pmt_times->at(i);
				else if(card_ID == 132 && pmt_ID == 19) trigger_time[2] = hit_pmt_times->at(i);
			}

			for (int i = 0; i < hit_mpmt_card_ids->size(); i++){
				if(card_ID == 130 && T0_IDs.count(pmt_ID)){
					T0_times->Fill(hit_pmt_times->at(i) - trigger_time[0]);
					T0_avg_time += hit_pmt_times->at(i) - trigger_time[0];
					T0_hits++;
				}

				else if(card_ID == 131 && T1_IDs.count(pmt_ID)){
					T1_times->Fill(hit_pmt_times->at(i) - trigger_time[1]);
					T1_avg_time += hit_pmt_times->at(i) - trigger_time[1];
					T1_hits++;
				}

			}
			T0_avg_time /= T0_hits;
			T1_avg_time /= T1_hits;
			// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
			double bm_trig_time_0;
			double bm_trig_time_1;
			for (int i = 0; i < beamline_pmt_tdc_ids->size(); i++){
				pmt_ID = beamline_pmt_tdc_ids->at(i);
				if (pmt_ID == 31) bm_trig_time_0 = beamline_pmt_tdc_times->at(i); 
				else if (pmt_ID == 46) bm_trig_time_1 = beamline_pmt_tdc_times->at(i);
			}


			for (int i = 0; i < beamline_pmt_tdc_ids->size(); i++){
				pmt_ID = beamline_pmt_tdc_ids->at(i);

				if (T0_beamline_pmt_ids.count(pmt_ID)){
					BM_T0_times->Fill(beamline_pmt_tdc_times->at(i));
				}
				else if (T1_beamline_pmt_ids.count(pmt_ID)){
					BM_T1_times->Fill(beamline_pmt_tdc_times->at(i));
				}
			}
			// ##########################  BRB data analysis ###############################

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
		}	//end if condition (check for cuts)
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

	}// end of loop over events
	TCanvas* c_T0 = new TCanvas("", "", 1800, 900);
	c_T0->Divide(2);
	c_T0->cd(1);
	BM_T0_times->Draw("hist");
	c_T0->cd(2);
	T0_times->Draw("hist");
	TCanvas* c_T1 = new TCanvas("", "", 1800, 900);
	c_T1->Divide(2);
	c_T1->cd(1);
	BM_T1_times->Draw("hist");
	c_T1->cd(2);
	T1_times->Draw("hist");

	TCanvas* can = new TCanvas("", "", 900, 900);
	can->SetLeftMargin(0.15);
	tof_minus_trigger->Draw("hist");
	can->Print("c_tof_trigger.pdf");

	for (int i = 0; i < 16; i++){
		auto pmt = tof_pmt_ids.begin();
		std::advance(pmt, i);
		for (const auto& time : pmt_time[*pmt]){
			tof_times[i]->Fill(time);
		}
		for (const auto& charge : pmt_charge[*pmt]){
			tof_charges[i]->Fill(charge);
		}

	}

	gStyle->SetPalette(1);

	TCanvas* c_wf_trig = new TCanvas("", "", 900, 900);
	h_132_trigger_wftime->Draw("hist");

	TCanvas* c_wf_times = new TCanvas("", "", 1800, 900);
	c_wf_times->Divide(4, 4);
	for (int i = 0; i < 16; i++){
		c_wf_times->cd(i + 1);
		tof_waveform_times[i]->Draw("hist");
	}
	TCanvas* c_T0_tof = new TCanvas("", "", 900, 900);
	T0_TOF_times->Draw("hist");
	TCanvas* c_T1_tof = new TCanvas("", "", 900, 900);
	T1_TOF_times->Draw("hist");


	TCanvas* c_charge = new TCanvas("", "", 1800, 900);
	c_charge->Divide(4, 4);
	for (int i = 0; i < 16; i++){
		c_charge->cd(i+1);
		tof_charges[i]->Draw("hist");
	}
	c_charge->Print("c_charge.png");
	c_charge->Print("c_charge.pdf");

	TCanvas* c_time = new TCanvas("", "", 1800, 900);
	c_time->Divide(4, 4);
	for (int i = 0; i < 16; i++){
		c_time->cd(i+1);
		tof_times[i]->Draw("hist");
	}
	c_time->Print("c_time.png");
	c_time->Print("c_time.pdf");

	TCanvas* c_heatmaps = new TCanvas("", "", 1800, 900);
	c_heatmaps->Divide(4, 4);
	for (int i = 0; i < 16; i++){
		c_heatmaps->cd(i+1);
		tof_heatmaps[i]->Draw("hist");
	}
	c_heatmaps->Print("c_heatmaps.png");
	c_heatmaps->Print("c_heatmaps.pdf");


	for (const auto& pair : corresponding_pmt){
		for(int i = 0; i < pmt_time[pair.first].size(); i++){
			for (int j = 0; j < pmt_time[pair.second].size(); j++){
				tof_diff[pair.first]->Fill(pmt_time[pair.first][i] - pmt_time[pair.second][j]);
				tof_v_tof[pair.first]->Fill(pmt_time[pair.first][i], pmt_time[pair.second][j]);
				if(abs(pmt_time[pair.first][i] - pmt_time[pair.second][j]) < 100){
					tof_q_diff[pair.first]->Fill(pmt_charge[pair.first][i] - pmt_charge[pair.second][j]);
					tof_qvq[pair.first]->Fill(pmt_charge[pair.first][i], pmt_charge[pair.second][j]);


				}

			} 
		}
	}



	vector<TCanvas*> c_tvt(4);
	for(int ican = 0; ican < c_tvt.size(); ican++){
		c_tvt[ican] = new TCanvas(Form("can_%i", ican), "", 1800, 900);
		c_tvt[ican]->Divide(4, 2);
	}
	gStyle->SetOptFit(1);
	for(int i = 0; i < tof_diff.size(); i++){

		c_tvt[0]->cd(i+1);
		gPad->SetLeftMargin(0.15);
		double x_max = tof_diff[i]->GetBinCenter(tof_diff[i]->GetMaximumBin());
		tof_diff[i]->Fit("gaus", "RS", "", x_max - 10, x_max + 10);
		TF1* fit_fun = tof_diff[i]->GetFunction("gaus");
		tof_diff[i]->Draw("hist");
		fit_fun->Draw("same");
		c_tvt[1]->cd(i+1);
		gPad->SetLogz();
		gPad->SetRightMargin(0.15);
		gPad->SetLeftMargin(0.15);
		tof_v_tof[i]->SetStats(kFALSE);
		tof_v_tof[i]->Draw("colz");
		c_tvt[2]->cd(i+1);
		gPad->SetLeftMargin(0.15);
		x_max = tof_q_diff[i]->GetBinCenter(tof_q_diff[i]->GetMaximumBin());
		tof_q_diff[i]->Fit("gaus", "RS", "", x_max - 500, x_max + 500);
		fit_fun = tof_q_diff[i]->GetFunction("gaus");
		tof_q_diff[i]->Draw("hist");
		fit_fun->Draw("same");
		c_tvt[3]->cd(i+1);
		gPad->SetLogz();
		gPad->SetRightMargin(0.15);
		gPad->SetLeftMargin(0.15);
		tof_qvq[i]->SetStats(kFALSE);
		tof_qvq[i]->Draw("colz");

	}
	if(!apply_cuts){
		c_tvt[0]->Print("c_tof_diff_nocuts.png");
		c_tvt[1]->Print("c_tof_tvt_nocuts.png");
		c_tvt[2]->Print("c_tof_q_diff_nocuts.png");
		c_tvt[3]->Print("c_tof_qvq_nocuts.png");
	}
	else{
		c_tvt[0]->Print("c_tof_diff.png");
		c_tvt[1]->Print("c_tof_tvt.png");
		c_tvt[2]->Print("c_tof_q_diff.png");
		c_tvt[3]->Print("c_tof_qvq.png");
	}

	TCanvas* can1 = new TCanvas("", "", 1800, 900);
	card1_hist->Draw("hist");

	TCanvas* can2 = new TCanvas("", "", 1800, 900);
	card2_hist->Draw("hist");

	TCanvas* can3 = new TCanvas("", "", 1800, 900);
	card3_hist->Draw("hist");


	TCanvas* can_zoom = new TCanvas("", "", 900, 900);
	tof_q_diff[4]->GetXaxis()->SetRangeUser(-400, 400);
	tof_q_diff[4]->Draw("hist");
	can_zoom->Print("tof_zoom.pdf");
	// ############   LOOK AT POSITION RECONSTRUCTION    #####################################
	/*
	   int NDet = 16;
	   vector<array<double, 2>> scint_dimensions = {{51, 16.25}, 
	   {94, 16.25},
	   {112, 16.25},
	   {123, 16.25},
	   {123, 16.25},
	   {112, 11.25},
	   {94, 16.25},
	   {51, 16.25}
	   };
	   vector<double> start_y = {-58.875, -42.125, -25.375, -8.625, 8.625, 25.375, 42.125, 58.875};
	   vector<array<double, 2>> det_positions;
	   for (int i = 0; i < 2; i++){
	   for (int j = 0; j < NDet/2; j++){
	   array<double, 2> position;
	   if(i == 0) position = {scint_dimensions[j][0]/2, start_y[j]};
	   if(i == 1) position = {-scint_dimensions[j][0]/2, start_y[j]};

	   det_positions.push_back(position);
	   }
	   }//end for loop

	   double radius = 75;
	   std::vector<XYVector> SiPMs;
	   for (int i = 0; i  < NDet; i++){
	   SiPMs.push_back(construct_SiPM(det_positions[i][0], det_positions[i][1]));
	   }

	   vector<vector<array<double, 2>>> rec_coordinates(8);  //data arrays definition
	   vector<double> chi2_values;
	   for(int iscint = 0; iscint < 8; iscint++){
	   for (int ievent = 0; ievent < pmt_charge[iscint].size(); ievent++){

	   double tolerance = 10e-12;
	   int max_iter = 1000;
	   double x_rec = GetMinimum_chi2(max_iter, tolerance, pmt_charge.at(iscint).at(ievent), pmt_charge.at(corresponding_pmt.find(iscint)->second).at(ievent), SiPMs[iscint], SiPMs[iscint+8]);

	   array<double, 2> help_arr = {x_rec, start_y[iscint]};
	   rec_coordinates[iscint].push_back(help_arr);

	   cout << "End point: (" << x_rec << ", " << start_y[iscint] << ")" << endl; 
	   } //end loop over event in scintillator
	   }	//end loop over scintillators

	//##############################################
	vector<TBox*> scints;
	for(int i = 0; i < scint_dimensions.size(); i++){
	TBox* help_box = new TBox(-scint_dimensions[i][0]/2, start_y[i] - scint_dimensions[i][1]/2, scint_dimensions[i][0]/2, start_y[i] + scint_dimensions[i][1]/2) ;
	scints.push_back(help_box);
	}//end create scintillator tiles
	 //
	 //#################################################

	 TH2D* h_hits = new TH2D("h_hist_final", "Histogram of hits of reconstructed points+ x[mm]; y[mm]", 41, -scint_dimensions[4][0]/2, scint_dimensions[4][0]/2,
	 8, start_y[0] - scint_dimensions[0][1]/2, start_y[7] + scint_dimensions[0][1]/2);
	 for (int i = 0; i < SiPMs.size() / 2; i++){
	 for(int j = 0; j < rec_coordinates[i].size(); j++){
	 h_hits->Fill(rec_coordinates[i][j][0], rec_coordinates[i][j][1]);
	 }
	 }

	 TCanvas* c_hist = new TCanvas("", "", 1800, 900);
	 draw_boxes(scints);
	 c_hist->cd(2);
	 h_hits->Draw("colz");
	 draw_boxes(scints);
	 */

	/*

	   TRandom3 rand;






	   int n_enum = 50;
	   auto points = sort_points(get_grid(radius, n_enum), scint_dimensions, start_y);		//Generating grid



	   long npoints = 0;
	   int n_crash = 0;
	   double blur = 0.001;			//defining parameters used in minimization and signal blurring

	   std::vector<std::vector<double>> signal_memory;

	   TCanvas* can = new TCanvas("", "", 1800, 900);

	   can->Divide(2);
	   can->cd(1);
	   TGraph* g_start = new TGraph(npoints);
	   TGraph* g_rec = new TGraph(npoints);
	   int ip = 0;
	   for (int i = 0; i < SiPMs.size() / 2; i++){
	   for(int j = 0; j < points[i].size(); j++){
	   g_rec->SetPoint(ip, rec_coordinates[i][j][0], rec_coordinates[i][j][1]);
	   g_start->SetPoint(ip, points[i][j].X(), points[i][j].Y());
	   ip++;
	   }
	   }

	   g_start->Draw("AP*");
	   draw_boxes(scints);
	   can->cd(2);
	   g_rec->Draw("AP*");
	   draw_boxes(scints);

	//############################################################################
	*/
} //end of code
