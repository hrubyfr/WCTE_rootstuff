//this is a simple code to open events_run0XXX_beamline.root
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>   
#include <string>
#include <nlohmann/json.hpp> 
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>

#include "constants.h"
#include "buffer.h"
#include "event.h"
#include "utils.h"

bool LEAD_GLASS_RUN = true; //change only this if LEAD 


int main(int argc, char *argv[]) {

     if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <run number>" << std::endl;
        return 1;
    }

    // Convert the argument to an integer (run number)
    int runNumber = std::stoi(argv[1]);

    std::string run_str = std::to_string(runNumber);
    EventBuffer buffer = EventBuffer(run_str);

    //name the file
    std::string fileName = "/home/wcte/Desktop/tuples/events_run0" + std::to_string(runNumber) + "_beamline_newFormat.root";
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
            
            // TODO verify one QDC per QDC

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

    /*
            ==================== Run over data once - check time windows ======================
    */    

    
    std::vector<double> tof_bins = {};
    for(double bin=-0.5; bin<500.5; bin++){
        tof_bins.push_back(bin);
    }
    std::vector<int> tof0_counts(500);
    std::vector<int> tof1_counts(500);

    int insert_index = 0;
    
    // bin the events. I don't like root... that's why it's like this 
    for (int i=0; i<event_vector.size(); i++){
        // only use the hits on the first two, those seem to be more well-behaved
        for (int pmtid=28; pmtid<28+nTOF_PMTs; pmtid++){
            // first for the TOF 0 
            for(int hit=0; hit<event_vector[i].TDC(pmtid).size(); hit++){
                insert_index = binary_search(event_vector[i].TDC(pmtid)[hit], tof_bins);
                tof0_counts[insert_index]++;
            }
            for(int hit=0; hit<event_vector[i].TDC(pmtid+nTOF_PMTs).size(); hit++){
                insert_index = binary_search(event_vector[i].TDC(pmtid+nTOF_PMTs)[hit], tof_bins);
                tof1_counts[insert_index]++;
            }   
        }
    }

    // scan over, get the bin with the highest value 
    int _tof0_peak = -1;
    int _tof0_peak_index = -1;
    int _tof1_peak = -1;
    int _tof1_peak_index = -1;
    for (int b_index=0; b_index<tof0_counts.size(); b_index++){
        if (tof0_counts[b_index]>_tof0_peak){
            _tof0_peak = tof0_counts[b_index];
            _tof0_peak_index = b_index;
        }
        if (tof1_counts[b_index]>_tof1_peak){
            _tof1_peak = tof1_counts[b_index];
            _tof1_peak_index = b_index;
        }
    }
    
    double tof0_peak_time = 0.5 + _tof0_peak_index;
    double tof1_peak_tiem = 0.5 + _tof1_peak_index;

    // accept hits in a 50ns window centered here
    std::vector<double> tof0_range = { tof0_peak_time -25, tof0_peak_time +25 };
    std::vector<double> tof1_range = { tof1_peak_tiem -25 , tof1_peak_tiem + 25};

    std::cout << "TOF0 peak " <<tof0_peak_time << std::endl;
    std::cout << "TOF1 peak " <<tof1_peak_tiem << std::endl;

    std::cout <<"TOF0 range "<< tof0_range[0] <<", "<<tof0_range[1] <<std::endl;
    std::cout <<"TOF1 range "<< tof1_range[0] <<", "<<tof1_range[1] <<std::endl;
    /*
            ==================== Fill Event Buffers on Second Pass======= ======================
    */    

    for (int i = 0; i < event_vector.size(); ++i) {

        

        int nT0_tdc_hits = 0;
        // int nT1_tdc_hits = 0;

        int nMisMatch = 0;

        for (int T0_PMTs = 0; T0_PMTs < 8; T0_PMTs ++){
            //append 1 if there is a hit in that PMT
            nT0_tdc_hits += (event_vector[i].nTDC(T0_PMTs)!=0);
            if (event_vector[i].nTDC(T0_PMTs)==0 && event_vector[i].QDC(T0_PMTs) > configData[Form("PMT_%i", T0_PMTs)]["QDC_pedestal"]){
                // std::cout << "Event " << i << ": There is mismatch for PMT " << T0_PMTs <<  std::endl;
                nMisMatch += 1;
            }
        }

        if (nMisMatch >= 2){
            QDC_shift -= 1;
        }

        double actg1 =0.0;
        double actg2 =0.0;

        for(int index=0; index<6; index++){
            actg1 += event_vector[i].QDC(10+index);
            actg2 += event_vector[i].QDC(16+index);
        }

        int nTDC_hits_hc0_hc1 = 0; 
        for(int PMT=22; PMT<24; PMT++){
            
            if (event_vector[i].TDC(PMT).size()>0){
                if (event_vector[i].TDC(PMT)[0] > -1) nTDC_hits_hc0_hc1 +=1;
            }
        }

        int nTDC_MUT_hits = 0; 
        for(int PMT=26; PMT<28; PMT++){
            if (event_vector[i].TDC(PMT).size()>0){
                if (event_vector[i].TDC(PMT)[0] > -1) nTDC_MUT_hits +=1;
            }
        }

        bool rejectEvent = false;

        //check the event mismatch where there are large QDC values but no TDC values, in one of the PMTs (take the OR between T0-0L and T0-0R)

        int PMT_to_check_matching = 7;
        double pedestal_max_QDC_value =  configData[Form("PMT_%i", PMT_to_check_matching)]["QDC_pedestal"] ; // this is the value that we should not be above of in case there are no TDC events
        double pedestal_max_QDC_value_tight =  configData[Form("PMT_%i", PMT_to_check_matching)]["QDC_pedestal"] ; // this is the value that we should not be above of in case there are no TDC events

        if (false) {// (event_vector[i].nTDC(PMT_to_check_matching) == 0 && event_vector[i].QDC(PMT_to_check_matching) > pedestal_max_QDC_value && (!eventMismatch)) ){
            if (!eventMismatch){
                std::cout << "The event for which there is a mismatch (T0-0L or T0-0R: nTDC==0, QDC>" << pedestal_max_QDC_value << ") is " << i << std::endl;
                std::cout << "The last event for which there was no mismatch was " << last_good_event << std::endl;
                std::cout << "QDC_shift " << QDC_shift << std::endl;
            }

            //buffer.Clear();
            eventMismatch = true;  
            

            
        } //trying to fix the shift

        if ((event_vector[i].nTDC(PMT_to_check_matching) == 0 && event_vector[i].QDC(PMT_to_check_matching) < pedestal_max_QDC_value)){
            //keepign track of the last well matched event where there is low QDC value and no TDC value
            last_good_event = i;
            buffer.Flush();
            eventMismatch = false;  
        }         

        //hole counter veto
        if (hc_veto_on){
            if (nTDC_hits_hc0_hc1 > 0){
                rejectEvent = true;
            } 
        }
        //muon tagging
        // if (nTDC_MUT_hits == 0) continue;
        //electron veto
        if (ele_veto_on){
            if (actg1 > 900){
                rejectEvent = true;
            }
        }

        // //first look at implemeting cut on TS charge  
        // if (QDC_values_perPMT[6]< 300) continue;
        // if (QDC_values_perPMT[7]> 1500) continue;

        //tof difference
        double tof0_count = 0.0;
        double tof1_count = 0.0;
        int n_tof0 = 0;
        int n_tof1 = 0;
        bool is_tof0 = false;
        bool is_tof1 = false;
        double tof_diff = 0.0;

        //t0, t1 difference 
        bool is_t0 = false;
        bool is_t1 = false;
        double t_diff = 0.0;
        double t0_count = 0.0;
        double t1_count = 0.0;
        int n_t0 = 0;
        int n_t1 = 0;
        double shifted_time = 0.0;
		
        bool keep = false; 

        

        //read in the TDCT0i reference values, if there is no reference, will not save the time for that event (soemthing must have gone wrong with the trigger?)
        double TDCT00_time = 0.;
        double TDCT01_time = 0.;
        double TDCT02_time = 0.;
        //to substract the reference time, for now seems to be some hardware disagreement, not implementing it  
				
        for (int PMT=0; PMT < 63; PMT++){
             if (PMT == TDCT00_PMT_id){
                 if (event_vector[i].TDC(PMT).size() > 0){
                    //why are we substreacting by 290 ??
                     TDCT00_time = event_vector[i].TDC(PMT)[0] - 290;
                     // std::cout << "The one TDC entry is empty for reference PMT " << PMT << std::endl;
                 }
             }
             if (PMT == TDCT01_PMT_id ){
                 if (event_vector[i].TDC(PMT).size() > 0){
                     TDCT01_time = event_vector[i].TDC(PMT)[0] -290;
                     // std::cout << "The one TDC entry is empty for reference PMT " << PMT << std::endl;
                 }
             }
             if (PMT == TDCT02_PMT_id){
                if (event_vector[i].TDC(PMT).size() > 0){
                     TDCT02_time = event_vector[i].TDC(PMT)[0] -290;
                     // std::cout << "The one TDC entry is empty for reference PMT " << PMT << std::endl;
                 }
             }
         }

	    // std::cout <<" tdct0 "<< TDCT00_time <<", 1 is "<<TDCT01_time <<", 2 is "<<TDCT02_time << std::endl;
            
        
        if (!eventMismatch && !rejectEvent){

            // std::cout << "Event of interest is " << i << " nTDC values of PMT of interest is " << nTDC_values_perPMT[ACT_PMT_of_interest] << " QDC value is " << QDC_values_perPMT[ACT_PMT_of_interest]<< std::endl;

            for (int PMT=0; PMT < 63; PMT++){                                                
                is_t0 = PMT>=0 && PMT<4; 
                is_t1 = PMT>=4 && PMT<8;

                is_tof0 = PMT>27 && PMT<44; 
                is_tof1 = PMT>=44 && PMT<60;

                buffer.Buffer(event_vector[i].QDC(PMT), BufferT::qdc_per_pmt, PMT);
                buffer.Buffer(event_vector[i].nTDC(PMT), BufferT::ntdc_per_pmt, PMT);
                charge_index_log++;

                /*
                        38 
                */
                // Iadjust to TDCT02 for the TOF PMTs 
                if (PMT==38){
                    offset = TDCT01_time;
                }else if(PMT==59){
                    offset = TDCT00_time;
                }else if (PMT>=28 && PMT <=59){
                    offset = TDCT02_time;
                }

                int this_tof_counter = 0; 
                double this_tof_accumulate = 0.0;
                
                for (int hit = 0 ; hit < event_vector[i].TDC(PMT).size(); hit ++){
                    buffer.Buffer(event_vector[i].TDC(PMT)[hit], BufferT::ht_all_bufferkind);
                    shifted_time =  event_vector[i].TDC(PMT)[hit] - offset; 
                    // std::cout << hit << std::endl;                    
                    if (is_tof0 || is_tof1){
                        if(is_tof0 && event_vector[i].TDC(PMT)[hit] > tof0_range[0] && event_vector[i].TDC(PMT)[hit]<tof0_range[1]){
                            //tof0_count +=shifted_time;
                            this_tof_counter++;
                            this_tof_accumulate += shifted_time;
                            buffer.Buffer(shifted_time, BufferT::tdc_per_pmt, PMT);
                            buffer.Buffer(charge_index_log, BufferT::event_no);
                        }else if(is_tof1 && event_vector[i].TDC(PMT)[hit] > tof1_range[0] && event_vector[i].TDC(PMT)[hit] < tof1_range[1]){
                            //tof1_count += shifted_time;
                            this_tof_counter ++;
                            this_tof_accumulate += shifted_time;
                            buffer.Buffer(shifted_time, BufferT::tdc_per_pmt, PMT); 
                            buffer.Buffer(charge_index_log, BufferT::event_no);
                        }else{
                            buffer.Buffer(NAN, BufferT::tdc_per_pmt, PMT); 
                            buffer.Buffer(charge_index_log, BufferT::event_no);
                        }
                    }else if(is_t0 || is_t1){
                        // no offset 
                        buffer.Buffer(event_vector[i].TDC(PMT)[hit], BufferT::tdc_per_pmt, PMT);
                        buffer.Buffer(charge_index_log, BufferT::event_no); 
                        
                        if (event_vector[i].TDC(PMT)[hit] > min_TS0_time_window && event_vector[i].TDC(PMT)[hit] < max_TS0_time_window){
                            if (is_t0){
                                t0_count += event_vector[i].TDC(PMT)[hit];
                                n_t0++;
                            }else if(is_t1){
                                t1_count += event_vector[i].TDC(PMT)[hit];
                                n_t1++;
                            }
                        }
                }else{
                    buffer.Buffer(event_vector[i].TDC(PMT)[hit], BufferT::tdc_per_pmt, PMT); 
                    buffer.Buffer(charge_index_log, BufferT::event_no);
                }        
            }
            if (is_tof0 && this_tof_counter>0){
                n_tof0++;
                tof0_count += this_tof_accumulate/this_tof_counter;
            }else if(is_tof1 && this_tof_counter>0){
                n_tof1++;
                tof1_count += this_tof_accumulate/this_tof_counter;
            }



        }
        //here, require that there is no more than 16 TDC entries within the window per TOF module, to remove 2 particle events
        if (LEAD_GLASS_RUN){
            keep = n_tof0 > 0 && n_tof1 > 0; // && n_tof0 < 16 && n_tof1 < 16;
        }{
            
            keep = n_tof0>10 && n_tof1>10; 
        }
        if (keep){
                buffer.Buffer(tof0_count /n_tof0, BufferT::tof0);
                buffer.Buffer(tof1_count /n_tof1, BufferT::tof1);

                tof_diff =  (tof1_count/n_tof1) - (tof0_count/n_tof0);
                t_diff = (t1_count/n_t1) - (t0_count/n_t0);

                buffer.Buffer(tof_diff, BufferT::tof_diff);
                buffer.Buffer(t_diff, BufferT::ts_diff);
                
                buffer.Buffer(actg1, BufferT::act_g1); // 
                buffer.Buffer(actg2, BufferT::act_g2);
                buffer.Buffer(n_tof0, BufferT::n_tof0);
                buffer.Buffer(n_tof1, BufferT::n_tof1);
            }else{
                // keeps these out of 2D and 1D arrays, but lets us keep everything as the right length
                buffer.Buffer(NAN, BufferT::tof0);
                buffer.Buffer(NAN, BufferT::tof1);
                buffer.Buffer(NAN, BufferT::tof_diff);
                buffer.Buffer(NAN, BufferT::ts_diff);

                // bin the number of hits 
                buffer.Buffer(actg1, BufferT::act_g1); 
                buffer.Buffer(actg2, BufferT::act_g2);
                buffer.Buffer(n_tof0, BufferT::n_tof0);
                buffer.Buffer(n_tof1, BufferT::n_tof1);
            }
        }
        
    }
    buffer.Flush();
    buffer.MakePlots();
    return 0;
}
