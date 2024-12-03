
#ifndef EVENTBUFF_H
#define EVENTBUFF_H

#include <TTree.h>
#include "constants.h"
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TVectorD.h>
#include <vector>
#include <nlohmann/json.hpp> 

enum BufferT{
    act_g1,
    act_g2,
    ntdc_per_pmt, //2d
    qdc_per_pmt, //2d
    tdc_per_pmt, //2d
    ht_all_bufferkind,
    tof0,
    tof1,
    tof_diff,
    ts_diff,
    n_tof0,
    n_tof1,  
    event_no
};

class EventBuffer{
    /*
        Rather than directly add events into the histograms, we buffer them temporarily.
        When we notice a miscount, we clear the buffer.
        When we notice another clear proton event, we flush the buffer.
        When we finish, we flush the buffer and fill histograms.
    */
    private:
        bool initialized = false; 

        int buf_ct = 0;

        std::string run_str;
        nlohmann::json configData;

        TH2D * h_q_nt[nPMTs];
        TH2D * h_q_t[nPMTs][nHits];
        TH1D * h_t[nPMTs];
        TH1D * h_q[nPMTs];
        TH1D * h_t_all;

        TH1D* ACT_g1;
        TH1D* ACT_g2;

        TH1D* tof0;
        TH1D* tof1;
        TH1D* TOF_tof_diff;
        TH1D* TOF_ts_diff;
        TH1D* tof0_nhit;
        TH1D* tof1_nhit;

        TH2D* ACTg1_TOF;
        TH2D* ACTg2_TOF;

        TH2D* ACT0L_ACT0R;
        TH2D* ACT1L_ACT1R;
        TH2D* ACTg2L_R;
        
        TH2D* TOFdiff_Tdiff;

        TH2D* ACT1_Lead;
        TH2D* ACT2_Lead;

        TH2D* Lead_ACT1;
        TH2D* Lead_ACT2;

        TH2D* Lead_Time;

        // now the vector buffers
        std::vector<std::vector<double>> ntdc_per_pmt_buffer = std::vector<std::vector<double>>(nPMTs);
        std::vector<std::vector<double>> tdc_per_pmt_buffer = std::vector<std::vector<double>>(nPMTs);
        std::vector<std::vector<double>> qdc_per_pmt_buffer = std::vector<std::vector<double>>(nPMTs);

        std::vector<double> htall_buf;

        std::vector<double> act_g1_buf;
        std::vector<double> act_g2_buf;
        std::vector<double> tof0_buf;
        std::vector<double> tof1_buf;
        std::vector<double> tof_diff_buf;
        std::vector<double> ts_diff_buf;
        std::vector<double> ntof0_buf;
        std::vector<double> ntof1_buf;
        std::vector<int> event_number_buffer;

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

        int tdct_ref=0;

        double offset = 0;
        double TDCT00_time = 0;
        double TDCT01_time = 0;
        double TDCT02_time = 0;

        bool keep;

    public:
        EventBuffer(std::string run_no);

        void UpdateOffsets(double TDCT00_time, double TDCT01_time, double TDCT02_time);

        // called per-entry 
        void Buffer(double what, BufferT kind);
        void Buffer(double what, BufferT kind, int pmt_id);
        //void Buffer(std::vector<double> what, BufferT kind);


        /*
            Clears the buffer deques of all contents.
        */
        void Clear();

        /*
            Flushes the buffer deques and fills the histograms. 
        */
        void Flush();

        void MakePlots();

};


#endif