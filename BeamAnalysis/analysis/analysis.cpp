#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>   
#include <string>
#include "nlohmann/json.hpp" 
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>

#include "event.h"
#include "constants.h"

int main(int argc, char *argv[]) {

     if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <run number>" << std::endl;
        return 1;
    }

    // Convert the argument to an integer (run number)
    int runNumber = std::stoi(argv[1]);

    std::string run_str = std::to_string(runNumber);
    // EventBuffer buffer = EventBuffer(run_str);

    //name the file
    std::string fileName = "../Data/events_run0" + std::to_string(runNumber) + "_beamline_newFormat.root";
    //std::string fileName = "/Users/bsmithers/software/data/tuples/events_run0" + std::to_string(runNumber) + "_beamline_newFormat.root";


    // Open the ROOT file
    TFile *file = TFile::Open(fileName.c_str());

    // Check if the file opened successfully
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open the ROOT file!" << std::endl;
        return 1;
    }

    // List contents of the file
    file->ls();

    // retrieve the beam info
    TTree *tree = (TTree*)file->Get("beam_monitor_info");
    if (tree) {
        // Do something with the tree, like print its structure
        tree->Print();
    }

    // Set up to read the branch with a std::vector<double> type
    std::vector<double> *beamline_tdc_time = nullptr;
    tree->SetBranchAddress("beamline_tdc_time", &beamline_tdc_time);

    std::vector<double> *beamline_tdc_time_id = nullptr;
    tree->SetBranchAddress("beamline_tdc_time_id", &beamline_tdc_time_id);


    std::vector<double> *beamline_qdc_charge = nullptr;
    tree->SetBranchAddress("beamline_qdc_charge", &beamline_qdc_charge);

    std::vector<double> *beamline_qdc_charge_id = nullptr;
    tree->SetBranchAddress("beamline_qdc_charge_id", &beamline_qdc_charge_id);

    int beamline_tdc_nhits = 0;
    tree->SetBranchAddress("beamline_tdc_nhits", &beamline_tdc_nhits);

    //config file with the channel names and info
    nlohmann::json configData;
    std::ifstream inFile("TDC_and_QDC_channels.json");
    if (inFile.is_open()) {
        inFile >> configData;
        inFile.close();

    }// config file is open




    double offset = 0.0;
    int charge_index_log = 0;

    bool eventMismatch = false;
    int last_good_event = 0;

	double run = 0.0;
    // Loop over the entries in the tree
    Long64_t nEntries = tree->GetEntries();

    int TDC_shift = 0; //shifting dynamically by the number of event needed to recover from event mismatches 
    int QDC_shift = 0; //shifting dynamically by the number of event needed to recover from event mismatches 

    /*
        ==================== Read in Data, build event vector ====================
    */

    std::vector<Event> event_vector = {};

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        std::vector<std::vector<double>> TDC_values_perPMT= std::vector<std::vector<double>>(nPMTs);
        std::vector<double> QDC_values_perPMT = std::vector<double>(nPMTs);
        std::vector<int> nTDC_values_perPMT(nPMTs);


        for (int TDC_entry = 0; TDC_entry < beamline_tdc_time_id->size(); TDC_entry++) {            

            int PMT_number = (*beamline_tdc_time_id)[TDC_entry];
            std::string detector = configData[Form("PMT_%i", PMT_number)]["detector"];


            //check if the previous entry has already been filled to get warning out
            TDC_values_perPMT[PMT_number].push_back((*beamline_tdc_time)[TDC_entry]*25e-3);
            //insrease the counter 

            nTDC_values_perPMT[PMT_number] += 1;

        }//read TDC

        // std::cout << QDC_shift << std::endl;
        //do not apply the shift in event number to the QDC entries, only the TDC ones are fixed
        // tree->GetEntry(i + QDC_shift);

        for (int QDC_entry = 0; QDC_entry < beamline_qdc_charge->size(); QDC_entry++) {
            int PMT_number = (*beamline_qdc_charge_id)[QDC_entry];
            
            //PMT numbers: 60-62 # TDC refs, 63+: HD elements
            if(PMT_number < 60){
                QDC_values_perPMT[PMT_number] = (*beamline_qdc_charge)[QDC_entry];
            }
        };//read in QDC

        event_vector.push_back(Event(
            TDC_values_perPMT, 
            QDC_values_perPMT,
            nTDC_values_perPMT
        ));

    }
    std::cout << "filling event_vector complete" << std::endl;
    while(true){
    	int argument;
    	std::cin >> argument;
    	std::cout << std::endl;
    	if (std::cin.fail()){break;}
    	else{
    		std::cout<< "charge of the second event on " << argument + 1 <<". PMT is " << event_vector[1].QDC(argument)<<std::endl;
    	}
    }

}