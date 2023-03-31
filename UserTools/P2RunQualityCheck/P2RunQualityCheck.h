#ifndef P2RunQualityCheck_H
#define P2RunQualityCheck_H

#include <string>
#include <iostream>
#include <deque>

#include "Tool.h"
#include "ADCPulse.h"
#include "Hit.h"
#include "Particle.h"
#include "BeamStatus.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"


/**
 * \class P2RunQualityCheck
 *
 *
 * $ Author: Julie He
 * $ Date: 16 Sept 2021
 * $ Contact: juhe@ucdavis.edu
**/

class P2RunQualityCheck: public Tool {


 public:

  P2RunQualityCheck(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  void InitTree();
  void InitHist();
  void WriteHist();
  void InitGraph();
  void WriteGraph();

 private:

  // config variables
  int verbosity;
  std::string outpfile_prefix;
  int maxEntries;
  uint32_t globalRunNumber;

  uint32_t trigword;
  int trigext;
  int prevpart;
  int globalEntry;
  int beamEntry;
  int n_beam_evts;
  int n_beam_pot;
  int n_ok_evts;

  std::ofstream beampot_file;

  Geometry *fGeo = nullptr;

  std::map<int, double> map_chankey2spe;
  std::map<unsigned long, vector<Hit>> *tdcdata = nullptr;
	std::vector<BoostStore> *theMrdTracks; //the actual tracks
  std::map<double, std::vector<Hit>> *m_all_clusters = nullptr;

  // ROOT 
  TFile *p2rqc_root_outp = nullptr;

  TLine *l = nullptr;

  // ROOT histograms
  //1D
  TH1D *h_beampot_ok_all = nullptr;
  TH1D *h_beampot_nok_all = nullptr;
  TH1D *h_pot_beam = nullptr;
  TH1D *h_pot_led = nullptr;
  TH1D *h_pot_cosmics = nullptr;
  TH1D *h_pot_ok_beam = nullptr;
  TH1D *h_pot_ok_led = nullptr;
  TH1D *h_pot_ok_cosmics = nullptr;
  TH1D *h_pot_nok_beam = nullptr;
  TH1D *h_pot_nok_led = nullptr;
  TH1D *h_pot_nok_cosmics = nullptr;

  //2D
  TH2D *h2_pot_all = nullptr;
  TH2D *h2_pot_all_zoomx = nullptr;
  TH2D *h2_pot_all_zoomxy = nullptr;
  TH2D *h2_pot_ok_beam = nullptr;

  // ROOT graphs
  TCanvas *c_short_part = nullptr;
  TGraph *g_short_part = nullptr;
  TCanvas *c_short_all = nullptr;
  TGraph *g_short_all = nullptr;

  // ROOT trees
  TTree *t_tank = nullptr;
  uint32_t fRunNumber;
  uint32_t fSubrunNumber;
  uint32_t fRunType;
  uint64_t fStartTime;
  uint32_t fEventNumber;
  uint64_t fEventTimeTank;
  int fPartNumber;
  int fVetoHit;
  std::vector<double> fHitPE;

  double fBeamPOT;
  int fBeamOK;

  TTree *t_cluster = nullptr;
  double fClusterCharge;
  double fClusterTime;
  double fClusterPE;
  int fClusterHits;
  double fClusterChargeBalance;
  std::vector<double> fClusterHitPE;

  // verbosity variables
  int v_error = 0;
  int v_warning = 1;
  int v_message = 2;
  int v_debug = 3;

};


#endif
