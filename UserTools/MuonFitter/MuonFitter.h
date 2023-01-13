#ifndef MuonFitter_H
#define MuonFitter_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TTree.h"

#include "Hit.h"

/**
 * \class MuonFitter
 *
 * $ Author: J.He              $
 * $ Date: 2023/01/02          $
 * Contact: juhe@ucdavis.edu
*/

class MuonFitter: public Tool {


  public:

    MuonFitter(); ///< Simple constructor
    bool Initialise(std::string configfile,DataModel &data);
    bool Execute();
    bool Finalise();


  private:

    //constants
    double LUX_AREA = 470.;          //10-inch R7081 Hamamatsu
    double ETEL_AREA = 613.;         //11-inch D784UKFLB ETEL
    double HAMAMATSU_AREA = 330.;    //8-inch R5912-100 Hamamatsu
    double WATCHBOY_AREA = 470.;     //10-inch R7081 Hamamatsu
    double WATCHMAN_AREA = 470.;     //10-inch R7081-100 Hamamatsu

    double PHOTON_DENSITY = 328.;    // [#photons/cm]
    double CHERENKOV_ANGLE = 42.*TMath::Pi()/180.;   //or is it 41? [rads]
    double ALPHA_FACTOR = PHOTON_DENSITY / (4.*TMath::Pi()*TMath::Sin(CHERENKOV_ANGLE)*TMath::Tan(CHERENKOV_ANGLE));


    //logging variables
    int verbosity;
    int v_error = 0;
    int v_warning = 1;
    int v_message = 2;
    int v_debug = 3;
    std::string logmessage;

    //config variables
    bool isData;

    TFile *root_outp = nullptr;

    //histograms
    TH1D *h_alpha = nullptr;
    TH1D *h_exp_PE = nullptr;

    //geometry variables
    Geometry *geom = nullptr;
    double tank_height;
    double tank_radius;
    int n_tank_pmts;
    int n_mrd_pmts;
    int n_veto_pmts;

    //maps
    std::map<int, double> ChannelKeyToSPEMap;
    std::map<int, double> m_pmt_eff;
    std::map<int, double> m_pmt_alpha;
    std::map<double, std::vector<Hit>> *m_all_clusters = nullptr;
    std::map<double, std::vector<MCHit>> *m_all_clusters_MC = nullptr;
    std::map<double, std::vector<unsigned long>> *m_all_clusters_detkeys =nullptr;

};


#endif
