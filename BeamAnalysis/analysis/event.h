#ifndef EVENT_H
#define EVENT_H

#include "constants.h"

class Event{
    public:
        Event(){}
        Event(std::vector<std::vector<double>> _TDC_values_perPMT, std::vector<double> _QDC_values_perPMT,std::vector<int> _nTDC_values_perPMT) :
                TDC_values_perPMT(_TDC_values_perPMT), 
                QDC_values_perPMT(_QDC_values_perPMT),
                nTDC_values_perPMT(_nTDC_values_perPMT){}

        std::vector<double> TDC(int pmt_index){return TDC_values_perPMT[pmt_index];}
        double QDC(int pmt_index){return QDC_values_perPMT[pmt_index];}
        int nTDC(int pmt_index){return nTDC_values_perPMT[pmt_index];}

    private:
        std::vector<std::vector<double>> TDC_values_perPMT;
        std::vector<double> QDC_values_perPMT;
        std::vector<int> nTDC_values_perPMT;
};

#endif