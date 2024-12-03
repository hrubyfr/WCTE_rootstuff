#include "buffer.h"
#include <TCanvas.h>
#include <TStyle.h>

#include <iostream>
#include <fstream> 
#include <exception>

EventBuffer::EventBuffer(std::string runno){

    std::ifstream inFile("TDC_and_QDC_channels.json");
    if (inFile.is_open()) {
        inFile >> configData;
        inFile.close();
    }// config file is open

    run_str= runno;

    h_t_all = new TH1D("t_allPMTs", "", 300, 150, 300);

    ACT_g1 = new TH1D("act group 1", "", 2000, 0, 4000);
    ACT_g2 = new TH1D("act group 2", "", 2000, 0, 4000);

    tof0 = new TH1D("tof0 Group", "", 250, 0, 500);
    tof1 = new TH1D("tof1 Group", "", 250, 0, 500);
    TOF_tof_diff = new TH1D("TOF_tof_diff", "", 120, 40, 80);
    TOF_ts_diff = new TH1D("TOF_ts_diff", "", 100, -100, 100);
    tof0_nhit = new TH1D("tof0 Group nHits", "", 17, -0.5, 16.5);
    tof1_nhit = new TH1D("tof1 Group nHits", "", 17, -0.5, 16.5); 

    ACTg1_TOF = new TH2D("ACTg1_TOF", "", 60, 30, 90, 100, 0, 4000);
    ACTg2_TOF = new TH2D("ACTg2_TOF", "", 60, 30, 90, 100, 0, 4000);

    ACT0L_ACT0R = new TH2D("ACT0L_ACT0R", "", 100, 0, 4000, 100, 0, 4000);
    ACT1L_ACT1R = new TH2D("ACT1L_ACT1R", "", 100, 0, 4000, 100, 0, 4000);
    ACTg2L_R = new TH2D("ACTg2L_R", "", 200, 0, 4000, 200, 0, 4000);


    TOFdiff_Tdiff = new TH2D("TOFdiffTdiff", "", 100, 30, 90, 100, -40, 40);

    ACT1_Lead = new TH2D("ACT02 Lead", "", 100, 0, 8000, 100, 220, 2000);
    ACT2_Lead = new TH2D("ACT35 Lead", "", 100, 0, 8000, 100, 220, 2000);

    Lead_ACT1 = new TH2D("Lead ACT02 ", "", 100, 0, 2000, 100, 0, 8000);
    Lead_ACT2 = new TH2D("Lead ACT35", "", 100, 0, 2000, 100, 0, 8000);

    Lead_Time = new TH2D("tof_diff vs Lead Charge", "", 100, 40, 65, 100, 0, 4000);

        //initialise them
    for (int PMT=0; PMT < nPMTs; PMT++){
        h_q_nt[PMT] = new TH2D(("q_nt_PMT" + std::to_string(PMT)).c_str(), "", 5, 0, 5, 100, 0, 4200);
        h_t[PMT] = new TH1D(("t_PMT" + std::to_string(PMT)).c_str(), "", 500, 0,700);
        if (PMT>9 && PMT<22){
            h_q[PMT] = new TH1D(("q_PMT" + std::to_string(PMT)).c_str(), "", 200, 0, 2000);
                }else if(PMT==25){
            h_q[PMT] = new TH1D(("q_PMT" + std::to_string(PMT)).c_str(), "", 200, 0, 2000);
                }else{
            h_q[PMT] = new TH1D(("q_PMT" + std::to_string(PMT)).c_str(), "", 100, 0, 4200);
        }
        for (int hit=0; hit < nHits; hit ++){
            h_q_t[PMT][hit] = new TH2D(("q_nt_PMT" + std::to_string(PMT)+"_hit"+ std::to_string(hit)).c_str(), "", 100, 0, 500, 100, 0, 4200);
        }
        
    }
}

void EventBuffer::UpdateOffsets(double _TDCT00_time, double _TDCT01_time, double _TDCT02_time){
    TDCT00_time = _TDCT00_time;
    TDCT01_time = _TDCT01_time;
    TDCT02_time = _TDCT02_time;
}

void EventBuffer::Buffer(double what, BufferT kind){
    if(kind==BufferT::ntdc_per_pmt || kind==BufferT::qdc_per_pmt || kind==BufferT::tdc_per_pmt){
        throw std::runtime_error("Specify pmt_id!");
    }
    buf_ct ++;

    switch(kind){
        case BufferT::act_g1:
            act_g1_buf.push_back(what);
            break;
        case BufferT::act_g2:
            act_g2_buf.push_back(what);
            break;
        case BufferT::tof0:
            tof0_buf.push_back(what);
            break;
        case BufferT::tof1:
            tof1_buf.push_back(what);
            break;
        case BufferT::tof_diff:
            tof_diff_buf.push_back(what);
            break;
        case BufferT::ts_diff:
            ts_diff_buf.push_back(what);
            break;
        case BufferT::n_tof0:
            ntof0_buf.push_back(what);
            break;
        case BufferT::n_tof1:
            ntof1_buf.push_back(what);
            break;
        case BufferT::event_no:
            event_number_buffer.push_back(what);
            break;
        case BufferT::ht_all_bufferkind:
            htall_buf.push_back(what);
            break;
        default:
            std::runtime_error("Unreachable state");
            break;
    }
}

void EventBuffer::Buffer(double what, BufferT kind, int pmt_id){
    if(kind!=BufferT::ntdc_per_pmt && kind!=BufferT::qdc_per_pmt && kind!=BufferT::tdc_per_pmt){
        throw std::runtime_error("Do NOT specify pmt_id!");
    }
    if(pmt_id<0 || pmt_id>nPMTs-1){
        throw std::runtime_error("Invalid pmt id");
    }
    buf_ct++;
    switch (kind){
        case BufferT::ntdc_per_pmt:
            ntdc_per_pmt_buffer[pmt_id].push_back(what);
            break;
        case BufferT::qdc_per_pmt:
            qdc_per_pmt_buffer[pmt_id].push_back(what);
            break;
        case BufferT::tdc_per_pmt:
            tdc_per_pmt_buffer[pmt_id].push_back(what);
            break;
        default : 
            std::runtime_error("Unreachable state");
            break;
    }

}

void EventBuffer::Clear(){
    std::cout << "Cleaning event buffer" <<std::endl;
    ntdc_per_pmt_buffer.clear();
    tdc_per_pmt_buffer.clear();
    qdc_per_pmt_buffer.clear();

    act_g1_buf.clear();
    act_g2_buf.clear();
    
    tof0_buf.clear();
    tof1_buf.clear();
    tof_diff_buf.clear();
    ts_diff_buf.clear();

    ntof0_buf.clear();
    ntof1_buf.clear();

    ntdc_per_pmt_buffer = std::vector<std::vector<double>>(nPMTs);
    tdc_per_pmt_buffer =  std::vector<std::vector<double>>(nPMTs);
    qdc_per_pmt_buffer = std::vector<std::vector<double>>(nPMTs);

    act_g1_buf = {};
    act_g2_buf = {};
    
    tof0_buf  ={};
    tof1_buf  ={};
    tof_diff_buf ={};
    ts_diff_buf ={};

    ntof0_buf ={};
    ntof1_buf ={};

}

void EventBuffer::Flush(){    
    // These have ONE entry in the buffer per entry in root file 
    for(int i=0; i<act_g1_buf.size();i++){

        // n hits! 
        tof0_nhit->Fill(ntof0_buf[i]);
        tof1_nhit->Fill(ntof1_buf[i]);

        ACT0L_ACT0R->Fill(qdc_per_pmt_buffer[10][i], qdc_per_pmt_buffer[11][i]);
        ACT1L_ACT1R->Fill(qdc_per_pmt_buffer[12][i], qdc_per_pmt_buffer[13][i]);
        ACTg2L_R->Fill(qdc_per_pmt_buffer[16][i] + qdc_per_pmt_buffer[18][i] + qdc_per_pmt_buffer[20][i], qdc_per_pmt_buffer[17][i] + qdc_per_pmt_buffer[19][i] + qdc_per_pmt_buffer[21][i]);
        TOFdiff_Tdiff->Fill(tof_diff_buf[i], ts_diff_buf[i]);

        ACT_g1->Fill(act_g1_buf[i]);
        ACT_g2->Fill(act_g2_buf[i]);
        ACTg1_TOF->Fill(tof_diff_buf[i], act_g1_buf[i]);
        ACTg2_TOF->Fill(tof_diff_buf[i], act_g2_buf[i]);
        if (!std::isnan(tof_diff_buf[i])){

            TOF_tof_diff->Fill(tof_diff_buf[i]);
            TOF_ts_diff->Fill(ts_diff_buf[i]);
            tof0->Fill(tof0_buf[i]);
            tof1->Fill(tof1_buf[i]);

        }
        ACT1_Lead->Fill( act_g1_buf[i], qdc_per_pmt_buffer[25][i]);
        ACT2_Lead->Fill( act_g2_buf[i], qdc_per_pmt_buffer[25][i]);
        Lead_ACT1->Fill( qdc_per_pmt_buffer[25][i], act_g1_buf[i]);
        Lead_ACT2->Fill( qdc_per_pmt_buffer[25][i], act_g2_buf[i]);
        Lead_Time->Fill(tof_diff_buf[i], qdc_per_pmt_buffer[25][i]);
    }

    // At least one entry per PMT 
    for (int pmt=0; pmt<tdc_per_pmt_buffer.size(); pmt++){

        // multiple hits per PMT         
        for (int i=0; i<tdc_per_pmt_buffer[pmt].size(); i++){
            h_t[pmt]->Fill(tdc_per_pmt_buffer[pmt][i]);
            h_q_t[pmt][0]->Fill(tdc_per_pmt_buffer[pmt][i], qdc_per_pmt_buffer[pmt][event_number_buffer[i]]);
        }

        // just one entry per PMT 
        for (int i=0; i<qdc_per_pmt_buffer[pmt].size(); i++){
            h_q[pmt]->Fill(qdc_per_pmt_buffer[pmt][i]);
            h_q_nt[pmt]->Fill(ntdc_per_pmt_buffer[pmt][i], qdc_per_pmt_buffer[pmt][i]);

        }
    }
    // for every event, but not every pmt, many hits per event
    for (int i=0; i<htall_buf.size();i++){
        h_t_all->Fill(htall_buf[i]);
    }

    Clear();
}

void EventBuffer::MakePlots(){

    TCanvas * c;
    TCanvas * d;
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);


    std::cout << "making plots"<<std::endl;
    std::cout << "buffer called "<<buf_ct<<" times"<<std::endl;

    c = new TCanvas();

    ACT0L_ACT0R->SetTitle(("Run " + run_str + " ACT0R vs ACT0L (HC veto = " + hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
    ACT0L_ACT0R->GetXaxis()->SetTitle("QDC ACT0L");
    ACT0L_ACT0R->GetYaxis()->SetTitle("QDC ACT0R");
    ACT0L_ACT0R->Draw();
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_ACT0LR_comparison.png").c_str());

    c = new TCanvas();
    ACT1L_ACT1R->SetTitle(("Run " + run_str + " ACT1R vs ACT1L (HC veto = " + hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
    ACT1L_ACT1R->GetXaxis()->SetTitle("QDC ACT1L");
    ACT1L_ACT1R->GetYaxis()->SetTitle("QDC ACT1R");
    ACT1L_ACT1R->Draw();
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_ACT1LR_comparison.png").c_str());

    c = new TCanvas();
    ACTg2L_R->SetTitle(("Run " + run_str + " ACT3-5 Charge Sum Left vs Right (HC veto = "+ hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
    ACTg2L_R->GetXaxis()->SetTitle("Sum QDC left");
    ACTg2L_R->GetYaxis()->SetTitle("Sum QDC right");
    ACTg2L_R->Draw("COLZ");
    ACTg2L_R->SetMaximum(150);
    ACTg2L_R->SetMinimum(0);
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_ACTgroup2LR_comparison.png").c_str());

    c = new TCanvas();
    h_t_all->SetTitle(("Run " + run_str + " all PMTs Time (HC veto = " + hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
    h_t_all->GetXaxis()->SetTitle("Time (ns)");
    h_t_all->GetYaxis()->SetTitle("Entries");
    h_t_all->Draw();
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_allPMTs_time.png").c_str());

    c = new TCanvas();
    ACT_g1->SetTitle(("Run " + run_str + " ACT first three boxes (HC veto = "+ hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
    ACT_g1->GetXaxis()->SetTitle("QDC");
    ACT_g1->GetYaxis()->SetTitle("Entries");
    ACT_g1->Draw();
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_ACT_group1.png").c_str());

    c = new TCanvas();
    ACT_g2->SetTitle(("Run " + run_str + " ACT last three boxes (HC veto = "+ hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
    ACT_g2->GetXaxis()->SetTitle("QDC");
    ACT_g2->GetYaxis()->SetTitle("Entries");
    ACT_g2->GetYaxis()->SetRange(0, 800);
    ACT_g2->Draw();
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_ACT_group2.png").c_str());

    c = new TCanvas();
    tof0->SetTitle(("Run " + run_str + " tof0 Mean On-Time Time").c_str());
    tof0->GetXaxis()->SetTitle("Time [ns]");
    tof0->GetYaxis()->SetTitle("Entries");
    tof0->Draw();
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_tof0_mean_ontime_time.png").c_str());

    c = new TCanvas();
    tof1->SetTitle(("tof1 Mean On-Time Time (HC veto = "+ hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
    tof1->GetXaxis()->SetTitle("Time [ns]");
    tof1->GetYaxis()->SetTitle("Entries");
    tof1->Draw();
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_tof1_mean_ontime_time.png").c_str());

    TF1 *gauss = new TF1 ("gauss","[0]*TMath::Exp(-0.5*pow((x-[1])/[2],2))+[3]",40,80);
    gauss->SetParameter(0,100);
    gauss->SetParameter(1,59);
    gauss->SetParameter(2,1.0);
    gauss->SetParameter(3,20);
    TOF_tof_diff->Fit("gauss");

    c = new TCanvas();
    TOF_tof_diff->SetTitle(("Run " + run_str + " On-Time Difference TOF (HC veto = "+ hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
    TOF_tof_diff->GetXaxis()->SetTitle("Time [ns]");
    TOF_tof_diff->GetYaxis()->SetTitle("Entries");
    TOF_tof_diff->Draw();
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_TOF0_TOF1_mean_ontime_time_diff.png").c_str());

    c = new TCanvas();
    TOF_ts_diff->SetTitle(("Run " + run_str + " On-Time Difference TS (wHCveto)").c_str());
    TOF_ts_diff->GetXaxis()->SetTitle("Time [ns]");
    TOF_ts_diff->GetYaxis()->SetTitle("Entries");
    TOF_ts_diff->Draw();
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_T0_T1_mean_ontime_time_diff.png").c_str());

    c = new TCanvas();
    TOFdiff_Tdiff->SetTitle(("Run " + run_str + " On-Time Difference TOF vs TS time of flight (HC veto = "+ hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
    TOFdiff_Tdiff->GetXaxis()->SetTitle("TOF diff Time [ns]");
    TOFdiff_Tdiff->GetYaxis()->SetTitle("TS diff Time [ns]");
    TOFdiff_Tdiff->Draw();
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_TvsTOF_mean_ontime_time_diff.png").c_str());


    

    c = new TCanvas();
    ACTg1_TOF->SetTitle(("Run " + run_str + " ACT0-2 Charge Sum vs tof_diff (HC veto = "+ hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
    ACTg1_TOF->GetXaxis()->SetTitle("Time [ns]");
    ACTg1_TOF->GetYaxis()->SetTitle("Sum QDC");
    ACTg1_TOF->GetZaxis()->SetTitle("Counts");
    ACTg1_TOF->SetStats(0);
    ACTg1_TOF->Draw("COLZ");
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_ACTgroup1vsTOF.png").c_str());

    c = new TCanvas();
    ACTg2_TOF->SetTitle(("Run " + run_str + " ACT3-5 Charge Sum vs tof_diff (HC veto = "+ hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
    ACTg2_TOF->GetXaxis()->SetTitle("Time [ns]");
    ACTg2_TOF->GetYaxis()->SetTitle("Sum QDC");
    ACTg2_TOF->GetZaxis()->SetTitle("Counts");
    ACTg2_TOF->SetStats(0);
    ACTg2_TOF->Draw("COLZ");
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_ACTgroup2vTOF.png").c_str());

    c = new TCanvas();
    ACT1_Lead->SetTitle(("Run " + run_str +  " ACT0-2 vs Lead Charge (HC veto = "+ hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +" )").c_str());
    ACT1_Lead->GetXaxis()->SetTitle("ACT0-2 QDC");
    ACT1_Lead->GetYaxis()->SetTitle("PbGlass QDC");
    ACT1_Lead->GetZaxis()->SetTitle("Counts");
    ACT1_Lead->SetStats(1);
    ACT1_Lead->Draw("COLZ");
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_ACT_Lead_1.png").c_str());

    c = new TCanvas();
    Lead_ACT1->SetTitle(("Run " + run_str +  " ACT0-2 vs Lead Charge (HC veto = "+ hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +" )").c_str());
    Lead_ACT1->GetYaxis()->SetTitle("ACT0-2 QDC");
    Lead_ACT1->GetXaxis()->SetTitle("PbGlass QDC");
    Lead_ACT1->GetZaxis()->SetTitle("Counts");
    Lead_ACT1->SetStats(1);
    Lead_ACT1->Draw("COLZ");
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_Lead_ACTg1.png").c_str());

    c = new TCanvas();
    Lead_ACT2->SetTitle(("Run " + run_str +  " ACT3-5 vs Lead Charge (HC veto = "+ hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +" )").c_str());
    Lead_ACT2->GetYaxis()->SetTitle("ACT3-5 QDC");
    Lead_ACT2->GetXaxis()->SetTitle("PbGlass QDC");
    Lead_ACT2->GetZaxis()->SetTitle("Counts");
    Lead_ACT2->SetStats(1);
    Lead_ACT2->Draw("COLZ");
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_Lead_ACTg2.png").c_str());

    c = new TCanvas(); 
    ACT2_Lead->SetTitle(("Run " + run_str +  " ACT3-5 vs Lead Charge (HC veto = "+ hc_veto_on_char + "ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
    ACT2_Lead->GetXaxis()->SetTitle("ACT3-5 QDC");
    ACT2_Lead->GetYaxis()->SetTitle("PbGlass QDC");
    ACT2_Lead->GetZaxis()->SetTitle("Counts");
    ACT2_Lead->SetStats(1);
    ACT2_Lead->Draw("COLZ");
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_ACT_Lead_2.png").c_str());

    c = new TCanvas();
    //tof0_nhit
    tof0_nhit->SetTitle(("Run " + run_str +  "TOF0 N Hits").c_str());
    tof0_nhit->GetXaxis()->SetTitle("N Hits");
    tof0_nhit->GetYaxis()->SetTitle("N Entries");
    tof0_nhit->SetStats(0);
    tof0_nhit->Draw();
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_tof0_nhits.png").c_str());

    c = new TCanvas();
    //tof0_nhit
    tof1_nhit->SetTitle(("Run " + run_str +  "TOF1 N Hits").c_str());
    tof1_nhit->GetXaxis()->SetTitle("N Hits");
    tof1_nhit->GetYaxis()->SetTitle("N Entries");
    tof1_nhit->SetStats(0);
    tof1_nhit->Draw();
    c->Print(("plots/run"+run_str+ "/Run" + run_str + "_tof1_nhits.png").c_str());


    for (int PMT=0; PMT < nPMTs; PMT++){
                    std::cout <<"Making plots for "<< PMT<<std::endl;
        //2d nTDChits vs QDC
        c = new TCanvas();
        std::string PMT_name = configData[Form("PMT_%i", PMT)]["detector"];
        h_q_nt[PMT]->SetTitle(("Run " + run_str + " PMT " + std::to_string(PMT) + " " + PMT_name + " QDC vs n TDC hits").c_str());
        h_q_nt[PMT]->GetXaxis()->SetTitle("Number of TDC hits");
        h_q_nt[PMT]->GetYaxis()->SetTitle("QDC");
        h_q_nt[PMT]->Draw();
        // c->Print(("plots/run"+run_str+ "/Run" + run_str + "_PMT" + std::to_string(PMT) + "_QDCvsnTDC.png").c_str());

        //1D TDC
        c = new TCanvas();
        // std::string PMT_name = configData[Form("PMT_%i", PMT)]["detector"];
        h_t[PMT]->SetTitle(("Run " + run_str + " PMT " + std::to_string(PMT) + " " + PMT_name + " Time \n (HC veto = " + hc_veto_on_char + " ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
        h_t[PMT]->GetXaxis()->SetTitle("Time (ns)");
        h_t[PMT]->GetYaxis()->SetTitle("Entries");
        h_t[PMT]->Draw();
        c->Print(("plots/run"+run_str+ "/Run" + run_str + "PMT"+ std::to_string(PMT)+ "_time.png").c_str());
        
        //1D QDC
        c = new TCanvas();
        // std::string PMT_name = configData[Form("PMT_%i", PMT)]["detector"];
        h_q[PMT]->SetTitle(("Run " + run_str + " PMT " + std::to_string(PMT) + " " + PMT_name + " Charge \n (HC veto = "+ hc_veto_on_char + " ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
        h_q[PMT]->GetXaxis()->SetTitle("QDC ");
        h_q[PMT]->GetYaxis()->SetTitle("Entries");
        h_q[PMT]->Draw();
        c->Print(("plots/run"+run_str+ "/Run" + run_str + "PMT"+ std::to_string(PMT)+ "_charge.png").c_str());



        for (int hit = 0 ; hit < 1; hit ++){
            d = new TCanvas();
            h_q_t[PMT][hit]->SetTitle(("Run " + run_str + " PMT " + std::to_string(PMT) + " " + PMT_name + " QDC vs TDC \n (HC veto = "+ hc_veto_on_char + " ele veto = " + ele_veto_on_char + " " + act_hit_veto_on_char +")").c_str());
            h_q_t[PMT][hit]->GetXaxis()->SetTitle("TDC (ns)");
            h_q_t[PMT][hit]->GetYaxis()->SetTitle("QDC");
            h_q_t[PMT][hit]->Draw();
            d->Print(("plots/run"+run_str+ "/Run" + run_str + "_PMT" + std::to_string(PMT) +"_QDCvsTDC_allHits.png").c_str());
        }
    }
}