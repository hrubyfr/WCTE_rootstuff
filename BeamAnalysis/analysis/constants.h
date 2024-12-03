#ifndef H_CONST
#define H_CONST

#include <string>

const int nPMTs = 64;
const int nHits = 6;

const double nTOF_cut = 10;

const int ACT_PMT_of_interest = 5;

const int hc_veto_on = 0;
const int ele_veto_on = 0;
const bool requireTDC_HitInPMTx = true;
const std::string hc_veto_on_char = std::to_string(hc_veto_on);
const std::string ele_veto_on_char = std::to_string(ele_veto_on);

const std::string act_hit_veto_on_char = "";

const int nTOF_PMTs = 16;


const double min_TS0_time_window = 50.;
const double max_TS0_time_window = 200.;

const int TDCT00_PMT_id = 60;
const int TDCT01_PMT_id = 61;
const int TDCT02_PMT_id = 62;




#endif