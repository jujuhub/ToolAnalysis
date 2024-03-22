#ifndef MuonFitter_H
#define MuonFitter_H

#include <algorithm>
#include <string>
#include <iostream>
//#include <unordered_map>
#include <boost/algorithm/string.hpp>

#include "Tool.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TMarker.h"
#include "TPolyMarker.h"
#include "TTree.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TPolyLine3D.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoSphere.h"
#include "TVirtualGeoTrack.h"
#include "TLatex.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"

#include "Hit.h"

/**
 * \class MuonFitter
 *
 * $ Author: J.He              $
 * $ Last Updated: 2023/09/18  $
 * Contact: juhe@ucdavis.edu
*/

class MuonFitter: public Tool {


  public:

    MuonFitter(); ///< Simple constructor
    bool Initialise(std::string configfile,DataModel &data);
    bool Execute();
    bool Finalise();
    Position Line3D(double x1, double y1, double z1, double x2, double y2, double z2, double C, std::string coord);
    void reset_3d();
    bool FileExists(std::string fname);
    void LoadTankTrackFits();
    double CalcTankdEdx(double input_E);
    double CalcMRDdEdx(double input_E);


  private:

    //constants
    double LUX_AREA = 470.;          //10-inch R7081 Hamamatsu
    double ETEL_AREA = 613.;         //11-inch D784UKFLB ETEL
    double HAMAMATSU_AREA = 330.;    //8-inch R5912-100 Hamamatsu
    double WATCHBOY_AREA = 470.;     //10-inch R7081 Hamamatsu
    double WATCHMAN_AREA = 470.;     //10-inch R7081-100 Hamamatsu

    double PHOTON_DENSITY = 328.;    // [#photons/cm]
    double CHER_ANGLE_DEG = 41.;    //[deg] 41 or 42?
    double CHER_ANGLE_RAD = CHER_ANGLE_DEG*TMath::Pi()/180.;   //[rads]
    double ALPHA_FACTOR = PHOTON_DENSITY / (4.*TMath::Pi()*TMath::Sin(CHER_ANGLE_RAD)*TMath::Tan(CHER_ANGLE_RAD));
    double I_O = 9.200e-5;   //=11.5*8*0.000001
    double I_H = 1.920e-5;   //=19.2*1*0.000001
    double I_GD = 5.888e-4;   //=19.2*1*0.000001
    double I_S = 1.632e-4;    //=10.6*16*0.000001
    double I_FE = 2.86e-4;    //=11*26*0.000001


    //logging
    int verbosity;
    int v_error = 0;
    int v_warning = 1;
    int v_message = 2;
    int v_debug = 3;
    std::string logmessage;

    //config variables
    bool isData;
    bool plot3d = false;
    bool draw3d_fmv = false;
    bool draw3d_mrd = false;
    bool save_hists = false;
    double PMTMRDOffset;
    double deltaL = 0.;
    double insideAngle = -5.;
    double outsideAngle = 5.;
    double PMTQCut = 3.;
    double EtaThreshold = 500.;
    bool display_truth = false;
    bool fit_mode = false;
    std::string tankTrackFitFile;

    //text files
    std::ofstream pos_file;
    std::ofstream cpp_file;
    std::ofstream pehits_file;
    std::ofstream truetrack_file;
    std::ofstream nhits_trlen_file;
    std::ofstream lg_ediff_file;

    //root,plots
    TFile *root_outp = nullptr;
    TCanvas *canvas_3d;
    TCanvas *canvas_vtxq;
    TCanvas *c_vtx_charge;
    TCanvas *c_vtx_detkey;
    TCanvas *c_charge_per_pmt;
    TCanvas *c_effarea_detkey;
    TCanvas *c_fpmt_detkey;
    TCanvas *c_effarea_ai;
    TCanvas *c_fpmt_ai;
    TCanvas *c_eta_ai;

    TCanvas *canvas_ev_display; //test saving of ev displays
    std::map<uint32_t, TCanvas *> m_evdisplays;

    TCanvas *c_h_tzero;
    TH1D *h_tzero = nullptr;

    TGraph *gr_vtx_charge_in = nullptr;
    TGraph *gr_vtx_charge_out = nullptr;
    TGraph *gr_qdensity_in = nullptr;
    TGraph *gr_qdensity_out = nullptr;
    TGraph *gr_vtx_detkey_in = nullptr;
    TGraph *gr_vtx_detkey_out = nullptr;
    TGraph *gr_effarea_detkey = nullptr;
    TGraph *gr_effarea_ai = nullptr;
    TGraph *gr_fpmt_detkey = nullptr;
    TGraph *gr_fpmt_ai = nullptr;
    TGraph *gr_eta_ai = nullptr;
    TGraph *gr_qincone_ai = nullptr;
    TGraph *gr_qoutcone_ai = nullptr;
    TGraph *gr_running_avg = nullptr;
    TH1D *h_alpha = nullptr;
    TH1D *h_expected_PE = nullptr;
    TH1D *h_phot_inc_angle = nullptr;
    TH1D *h_hit_angles = nullptr;
    TH1D *h_fitted_tank_track_len = nullptr;
    TH1D *h_closest_approach = nullptr;
    TH1D *h_num_mrd_layers = nullptr;
    TH1D *h_lastvtx_z = nullptr;
    TH1D *h_clusterhit_x = nullptr;
    TH1D *h_clusterhit_y = nullptr;
    TH1D *h_clusterhit_z = nullptr;
    TH1D *h_clusterhit_detkey = nullptr;
    TH1D *h_truevtx_angle = nullptr;
    TH1D *h_tankexit_to_pmt = nullptr;
    TH1D *h_tanktrack_ai = nullptr;
    TH1D *h_eff_area_pmt = nullptr;
    TH1D *h_fpmt = nullptr;
    TH2D *h_eta_ai = nullptr;
    TH1D *h_clusterhit_timespread = nullptr;
    TH1D *h_clusterhit_time = nullptr;
    TH1D *h_qincone_truevtx = nullptr;
    TH1D *h_qoutcone_truevtx = nullptr;
    TH2D *h_total_pe_hits = nullptr;
    TH2D *h_charge_detkey = nullptr;
    TH1D *h_truevtx_trueexit_track = nullptr;
    TH1D *h_pmt_charge = nullptr;
    TH1D *h_lr_avg_eta = nullptr;
    TH1D *h_avg_eta = nullptr;
    TH1D *h_truefitdiff_x = nullptr;
    TH1D *h_truefitdiff_y = nullptr;
    TH1D *h_truefitdiff_z = nullptr;
    TH1D *h_tdiff = nullptr;
    TH1D *h_uber_t0widths = nullptr;
    TH1D *h_true_tanktrack_len = nullptr;
    TH1D *h_fitted_tank_track = nullptr;
    TH1D *h_truefit_len_diff = nullptr;
    TH1D *h_vtxfit_x = nullptr;
    TH1D *h_vtxfit_y = nullptr;
    TH1D *h_vtxfit_z = nullptr;
    TH1D *h_truevtx_x = nullptr;
    TH1D *h_truevtx_y = nullptr;
    TH1D *h_truevtx_z = nullptr;
    TH2D *h_topview_fit = nullptr;
    TH2D *h_topview_truth = nullptr;
    TH2D *h_sideview_fit = nullptr;
    TH2D *h_sideview_truth = nullptr;
    TH1D *h_deltaR = nullptr;
    TH1D *h_transverse = nullptr;
    TH1D *h_parallel = nullptr;
    TH1D *h_deltaR_4pi = nullptr;
    TH2D *h_true_reco_E = nullptr;
    TH1D *h_true_reco_Ediff = nullptr;
    TH2D *h_tank_mrd_E = nullptr;
    TH2D *h_Ediff_frac_tank = nullptr;
    TH2D *h_Ediff_frac_mrd = nullptr;
    TH1D *h_mrd_eloss_diff = nullptr;
    TH2D *h_mrd_eloss_SBvANNIE = nullptr;
    TH1D *h_true_reco_Ediff_sb = nullptr;
    TH1D *h_true_reco_Ediff_sb_outerE = nullptr;
    TH1D *h_tank_track_diff_small = nullptr;
    TH1D *h_tank_track_diff_large = nullptr;
    TH1D *h_mrd_track_diff_small = nullptr;
    TH1D *h_mrd_track_diff_large = nullptr;
    TH1D *h_deltaR_small = nullptr;
    TH1D *h_deltaR_large = nullptr;
    TH2D *h_mrd_nlyrs_reco = nullptr;
    TH1D *h_mrd_track_diff = nullptr;
    TH1D *h_mrd_angle_diff = nullptr;
    TH1D *h_mrd_angle = nullptr;
    TH1D *h_clusterPE = nullptr;
    TH1D *h_clusterPE_fit = nullptr;
    TH1D *h_clusterPE_fit_haspion = nullptr;
    TH1D *h_clusterPE_lrg_ediff = nullptr;
    TH1D *h_clusterPE_lrg_ediff_haspion = nullptr;
    TH1D *h_pca_angle = nullptr;
    TH1D *h_pca_reco_angle = nullptr;
    TH1D *h_pca_true_angle = nullptr;

    //event variables
    int partnumber;
    uint32_t evnum;
    int runnumber;
    int subrunnumber;
    int mcevnum;
    uint16_t mctrignum;
    std::string mcFile;
    double trueVtxTime;
    double trueVtxX, trueVtxY, trueVtxZ;
    double trueDirX, trueDirY, trueDirZ;
    double trueStopVtxX, trueStopVtxY, trueStopVtxZ;
    double trueAngle;
    double trueTrackLengthInWater;
    double trueTrackLengthInMRD;
    double trueMuonEnergy;
    int nrings;
    std::vector<unsigned int> particles_ring;
    std::map<unsigned long, std::vector<Hit>>* m_hits = nullptr;

    //geometry
    Geometry *geom = nullptr;
    double tank_height;
    double tank_radius;
    double tank_center_x;
    double tank_center_y;
    double tank_center_z;
    double detector_version;
    std::string detector_config;
    int n_tank_pmts;
    int n_mrd_pmts;
    int n_veto_pmts;

    //3d geometry
    TGeoManager *ageom = nullptr;
    TGeoMaterial *vacuum = nullptr;
    TGeoMaterial *Fe = nullptr;
    TGeoMedium *Air = nullptr;
    TGeoMedium *Iron = nullptr;
    TGeoVolume *EXPH = nullptr;  //experimental hall
    char blockName[100];  //char array to form strings w/ sprintf
    int N = 0, maxN = 0;  //Nth node
    TGeoVolume *bBlock = nullptr;   //building block of PMTs
    std::map<unsigned long, int> detkey_to_node;
    std::map<unsigned long, double> x_pmt, y_pmt, z_pmt;
    std::map<unsigned long, double> x_pmt_dir, y_pmt_dir, z_pmt_dir;
    Int_t track_index = 0;
    Int_t ntracks = 0;
    TVirtualGeoTrack *track = nullptr;

    //maps,vectors
    std::map<int, double> ChannelKeyToSPEMap;
    std::map<unsigned long, vector<Hit>> *tdcdata = nullptr;
    std::map<int, double> m_pmt_area;
    std::map<int, double> m_pmt_eff;
    std::map<int, double> m_pmt_alpha;
    std::map<double, std::vector<Hit>> *m_all_clusters = nullptr;
    std::map<double, std::vector<MCHit>> *m_all_clusters_MC = nullptr;
    std::map<double, std::vector<unsigned long>> *m_all_clusters_detkeys = nullptr;
    std::map<std::string, double> m_tank_track_fits;
    std::vector<MCParticle> *mcParticles = nullptr;
    std::vector<int> LayersHit;
    std::vector<int> PMTsHit;

    //mrd
    int numTracksInEv = 0;
    std::vector<BoostStore>* mrdTracks;   //the reconstructed tracks
    Position mrdStartVertex;
    Position mrdStopVertex;
    Position tankExitPoint;
    double trackAngle = -69.;
    double trackAngleError = -69.;
    double penetrationDepth = -69.;
    Position mrdEntryPoint;
    int numLayersHit = 0;
    double energyLoss = 0.;
    double energyLossError = 0.;
    double mrdStartTime = -69.;

    //misc
    double avg_eta;
    double left_avg_eta;
    double right_avg_eta;
    int num_left_eta;
    int num_right_eta;
    double bestFitAi = 0;
};


#endif
