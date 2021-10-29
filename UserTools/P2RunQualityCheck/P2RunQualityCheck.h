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

  Geometry *fGeo = nullptr;

  std::map<int, double> map_chankey2spe;
  std::map<unsigned long, vector<Hit>> *tdcdata = nullptr;
	std::vector<BoostStore> *theMrdTracks; //the actual tracks
  std::map<double, std::vector<Hit>> *m_all_clusters = nullptr;

  uint32_t trigword;
  int trigext;
  int prevpart;
  int globalEntry_i;
  int beam_entry_i;
  int n_beam_evts;
  int n_beam_pot;
  int n_ok_evts;
  int e;
  int gi;
  int pi;
	int numsubevs;

  size_t max_deq_size;
  std::deque<double> pot_deq;
  std::vector<std::pair<int, double> > gradient_vec;

  // ROOT 
  TFile *p2rqc_root_outp = nullptr;

  // ROOT histograms
  int th2_xlim;
  TH2D *h_beam_ok_all = nullptr;
  TH2D *h_beam_pot_all = nullptr;
  TH2D *h_beam_pot_all_zoomx = nullptr;
  TH2D *h_beam_pot_all_zoomxy = nullptr;
  TH2D *h_beam_ok_beam = nullptr;
  TH2D *h_beam_pot_beam = nullptr;
  TH2D *h_tankcharge_part = nullptr;
  TH2D *h_tankcharge_part_zoom = nullptr;
  TH2D *h_tankcharge_all = nullptr;
  TH2D *h_tankcharge_all_zoom = nullptr;
  TH1D *h_clusterTime_all = nullptr;
  TH1D *h_clusterTime_short = nullptr;
  TH1D *h_clusterCharge_short = nullptr;
  TH1D *h_clusterPE_short = nullptr;
  TH1D *h_clusterTime_long = nullptr;
  TH1D *h_clusterCharge_long = nullptr;
  TH1D *h_clusterPE_long = nullptr;

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
