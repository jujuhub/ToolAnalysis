#include "MuonFitter.h"

MuonFitter::MuonFitter():Tool(){}


bool MuonFitter::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // Default values
  verbosity = 3;
  isData = true;
  PMTMRDOffset = 745.;

  // --------------------------------------------------
  // --- Retrieve config variables --------------------
  // --------------------------------------------------
  m_variables.Get("verbosity", verbosity);
  m_variables.Get("IsData", isData);
  m_variables.Get("LuxArea", LUX_AREA);
  m_variables.Get("EtelArea", ETEL_AREA);
  m_variables.Get("HamamatsuArea", HAMAMATSU_AREA);
  m_variables.Get("WatchboyArea", WATCHBOY_AREA);
  m_variables.Get("WatchmanArea", WATCHMAN_AREA);
  m_variables.Get("PMTMRDOffset", PMTMRDOffset);
  m_variables.Get("Plot3D", plot3d);
  m_variables.Get("Draw3DFMV", draw3d_fmv);
  m_variables.Get("Draw3DMRD", draw3d_mrd);
  m_variables.Get("SaveHistograms", save_hists);
  m_variables.Get("DeltaL", deltaL);
  m_variables.Get("InsideAngle", insideAngle);
  m_variables.Get("OutsideAngle", outsideAngle);
  m_variables.Get("PMTChargeThreshold", PMTQCut);
  m_variables.Get("EtaThreshold", EtaThreshold);
  m_variables.Get("DisplayTruth", display_truth);
  m_variables.Get("FitMode", fit_mode);
  m_variables.Get("TankTrackFitFile", tankTrackFitFile);


  // Output ROOT file
  std::string outfile;
  m_variables.Get("OutputFile", outfile);
  Log("MuonFitter Tool: Saving output into " + outfile, v_message, verbosity);
  root_outp = new TFile(outfile.c_str(), "RECREATE");


  // --------------------------------------------------
  // --- Initialize canvases, graphs, histograms ------
  // --------------------------------------------------
  // Canvases
  c_vtx_charge = new TCanvas("c_vtx_charge", "Total Charge Seen at Each Vertex", 800, 600);
  c_vtx_detkey = new TCanvas("c_vtx_detkey", "Num of PMTs In/Out Cone at Each Vertex", 800, 600);
  c_charge_per_pmt = new TCanvas("c_charge_per_pmt", "Charge Per PMT at Each Vertex", 800, 600);
  c_effarea_detkey = new TCanvas("c_effarea_detkey", "Effective PMT Area vs detkey", 800, 600);
  c_fpmt_detkey = new TCanvas("c_fpmt_detkey", "Effective PMT Area as Fraction of Frustum Area", 800, 600);
  c_effarea_ai = new TCanvas("c_effarea_ai", "Effective PMT Area vs ai (tank track)", 800, 600);
  c_fpmt_ai = new TCanvas("c_fpmt_ai", "Effective PMT Area as Fraction of Frustum Area", 800, 600);
  c_eta_ai = new TCanvas("c_eta_ai", "Charge per cm (eta)", 800, 600);
  c_h_tzero = new TCanvas("c_h_tzero", "t0", 800, 600);

  // Histograms
  h_alpha = new TH1D("h_alpha", "Alpha Value Distribution", 10, 0., 10.);
  h_alpha->GetXaxis()->SetTitle("alpha_i (x1e3)");   //TODO:latex

  h_expected_PE = new TH1D("h_expected_PE", "Expected Num of PE", 100, 0., 100.);
  h_expected_PE->GetXaxis()->SetTitle("expected # of PE");

  h_phot_inc_angle = new TH1D("h_phot_inc_angle", "Photon Incident Angle", 180, 0., 180.);
  h_phot_inc_angle->GetXaxis()->SetTitle("incident angle [deg]");

  h_hit_angles = new TH1D("h_hit_angles", "Hit Angles wrt Vertex and MRD Track Dir", 180, 0., 180.);
  h_hit_angles->GetXaxis()->SetTitle("hit angle [deg]");

  h_fitted_tank_track_len = new TH1D("h_fitted_tank_track_len", "Length of Fitted Tank Tracks", 500, 0., 500.);
  h_fitted_tank_track_len->GetXaxis()->SetTitle("track length [cm]");

  h_closest_approach = new TH1D("h_closest_approach", "Distance of Closest Approach (Btwn Vtx & Tank Ctr)", 500, 0., 500.);
  h_closest_approach->GetXaxis()->SetTitle("distance to tank center [cm]");

  h_num_mrd_layers = new TH1D("h_num_mrd_layers", "Num MRD Layers Hit", 50, 0, 50.);
  h_num_mrd_layers->GetXaxis()->SetTitle("# layers");
  
  //h_truevtx_z = new TH1D("h_truevtx_z", "Z coordinate of True Vertex", 300, -150., 150.);
  //h_truevtx_z->GetXaxis()->SetTitle("Z [cm]");

  h_lastvtx_z = new TH1D("h_lastvtx_z", "Z coordinate of Last Candidate Vertex", 300, -150., 150.);
  h_lastvtx_z->GetXaxis()->SetTitle("Z [cm]");

  h_clusterhit_x = new TH1D("h_clusterhit_x", "X coordinate of Single Cluster Hits", 300, -150., 150.);
  h_clusterhit_x->GetXaxis()->SetTitle("X [cm]");

  h_clusterhit_y = new TH1D("h_clusterhit_y", "Y coordinate of Single Cluster Hits", 400, -200., 200.);
  h_clusterhit_y->GetXaxis()->SetTitle("Y [cm]");

  h_clusterhit_z = new TH1D("h_clusterhit_z", "Z coordinate of Single Cluster Hits", 300, -150., 150.);
  h_clusterhit_z->GetXaxis()->SetTitle("Z [cm]");

  h_clusterhit_detkey = new TH1D("h_clusterhit_detkey", "Detkeys of Single Cluster Hits", 135, 330, 465);
  h_clusterhit_detkey->GetXaxis()->SetTitle("detkey");

  h_truevtx_angle = new TH1D("h_truevtx_angle", "Angle of True Muon Vertex", 360, -180., 180.);
  h_truevtx_angle->GetXaxis()->SetTitle("angle [deg]");

  h_tankexit_to_pmt = new TH1D("h_tankexit_to_pmt", "Distance from Tank Exit Point to PMT Position", 350, 0., 350.);
  h_tankexit_to_pmt->GetXaxis()->SetTitle("Ri [cm]");

  h_tanktrack_ai = new TH1D("h_tanktrack_ai", "Tank Track (ai)", 200, 0., 200.);
  h_tanktrack_ai->GetXaxis()->SetTitle("ai [cm]");

  h_eff_area_pmt = new TH1D("h_eff_area_pmt", "Effective Area of PMT Seen by Photons", 650, 0., 650.);
  h_eff_area_pmt->GetXaxis()->SetTitle("effective area [cm^2]");

  h_fpmt = new TH1D("h_fpmt", "area fraction f", 100, 0., 1.);
  h_fpmt->GetXaxis()->SetTitle("f");

  h_eta_ai = new TH2D("h_eta_ai", "#eta vs. Tank Track Length", 50, 0., 500., 500, 0., 5000.);
  h_eta_ai->GetXaxis()->SetTitle("ai [cm]");
  h_eta_ai->GetYaxis()->SetTitle("#eta [nPE/f]");

  h_clusterhit_timespread = new TH1D("h_clusterhit_timespread", "Time Spread of Clusters (Latest - Earliest)", 200, 0, 200);
  h_clusterhit_timespread->GetXaxis()->SetTitle("latest t - earliest t [ns]");

  h_clusterhit_time = new TH1D("h_clusterhit_time", "Hit Time of Hits in Cluster", 70000, 0, 70000);
  h_clusterhit_time->GetXaxis()->SetTitle("hit time [ns]");

  h_qincone_truevtx = new TH1D("h_qincone_truevtx", "Hits That Fall Inside Cone of True Vertex (<42 deg)", 500, 0, 500);
  h_qincone_truevtx->GetXaxis()->SetTitle("hit PE [#PE]");

  h_qoutcone_truevtx = new TH1D("h_qoutcone_truevtx", "Hits That Fall Outside Cone of True Vertex (>50 deg)", 500, 0, 500);
  h_qoutcone_truevtx->GetXaxis()->SetTitle("hit PE [#PE]");

  h_total_pe_hits = new TH2D("h_total_pe_hits", "Total PE vs Total Hits", 500, 0., 500., 8000, 0., 8000.);
  h_total_pe_hits->GetXaxis()->SetTitle("# hits");
  h_total_pe_hits->GetYaxis()->SetTitle("# PE");

  h_charge_detkey = new TH2D("h_charge_detkey", "Total PE vs PMT detkey", 134, 331, 465, 500, 0., 10000.);
  h_charge_detkey->GetXaxis()->SetTitle("detkey");
  h_charge_detkey->GetYaxis()->SetTitle("nPE");

  h_truevtx_trueexit_track = new TH1D("h_truevtx_trueexit_track", "Tank Track Length (true vertex to true tank exit)", 350, 0, 350);
  h_truevtx_trueexit_track->GetXaxis()->SetTitle("track length [cm]");

  h_pmt_charge = new TH1D("h_pmt_charge", "Total Charge Seen By PMT In A Cluster", 500, 0, 1000);
  h_pmt_charge->GetXaxis()->SetTitle("charge [pe]");

  h_lr_avg_eta = new TH1D("h_lr_avg_eta", "Avg #eta to Left/Right of True Tank Track Length", 100, 0, 5000);
  h_lr_avg_eta->GetXaxis()->SetTitle("eta [PE/cm]");

  h_avg_eta = new TH1D("h_avg_eta", "Overall Average Eta", 100, 0, 5000);
  h_avg_eta->GetXaxis()->SetTitle("eta [PE/cm]");

  h_truefitdiff_x = new TH1D("h_truefitdiff_x", "Fitted Vertex - True Vertex (X)", 300, -150., 150.);
  h_truefitdiff_x->GetXaxis()->SetTitle("diff (fit-true) [cm]");

  h_truefitdiff_y = new TH1D("h_truefitdiff_y", "Fitted Vertex - True Vertex (Y)", 300, -150., 150.);
  h_truefitdiff_y->GetXaxis()->SetTitle("diff (fit-true) [cm]");

  h_truefitdiff_z = new TH1D("h_truefitdiff_z", "Fitted Vertex - True Vertex (Z)", 300, -150., 150.);
  h_truefitdiff_z->GetXaxis()->SetTitle("diff (fit-true) [cm]");

  h_tdiff = new TH1D("h_tdiff", "Time Residuals", 1000, -1000., 1000.);
  h_tdiff->GetXaxis()->SetTitle("time residual [ns]");

  h_uber_t0widths = new TH1D("h_uber_t0widths", "Distribution of Timing Distribution Widths", 31, -1., 30.);
  h_uber_t0widths->GetXaxis()->SetTitle("std devs [ns]");

  h_true_tanktrack_len = new TH1D("h_true_tanktrack_len", "Distribution of True Track Lengths", 400, 0, 2000);
  h_true_tanktrack_len->GetXaxis()->SetTitle("track length [cm]");

  h_fitted_track_len = new TH1D("h_fitted_track_len", "Distribution of Fitted Track Length", 200, 0, 400);
  h_fitted_track_len->GetXaxis()->SetTitle("track length [cm]");

  h_truefit_len_diff = new TH1D("h_truefit_len_diff", "Difference in Fitted and True Track Lengths", 100, 0, 200);
  h_truefit_len_diff->GetXaxis()->SetTitle("abs(fit-truth) [cm]");

  h_vtxfit_x = new TH1D("h_vtxfit_x", "Fitted Vertex X", 60, -150., 150.);
  h_vtxfit_x->GetXaxis()->SetTitle("x [cm]");

  h_vtxfit_y = new TH1D("h_vtxfit_y", "Fitted Vertex Y", 60, -150., 150.);
  h_vtxfit_y->GetXaxis()->SetTitle("y [cm]");

  h_vtxfit_z = new TH1D("h_vtxfit_z", "Fitted Vertex Z", 60, -150., 150.);
  h_vtxfit_z->GetXaxis()->SetTitle("z [cm]");

  h_truevtx_x = new TH1D("h_truevtx_x", "True Vertex X", 60, -150., 150.);
  h_truevtx_x->GetXaxis()->SetTitle("x [cm]");

  h_truevtx_y = new TH1D("h_truevtx_y", "True Vertex Y", 60, -150., 150.);
  h_truevtx_y->GetXaxis()->SetTitle("y [cm]");

  h_truevtx_z = new TH1D("h_truevtx_z", "True Vertex Z", 60, -150., 150.);
  h_truevtx_z->GetXaxis()->SetTitle("z [cm]");

  h_topview_fit = new TH2D("h_topview_fit", "Distributions of FITTED Vertices in Tank (Top View)", 60, -150., 150., 60, -150., 150.);
  h_topview_fit->GetXaxis()->SetTitle("Z [cm]");
  h_topview_fit->GetYaxis()->SetTitle("X [cm]");

  h_topview_truth = new TH2D("h_topview_truth", "Distributions of TRUTH Vertices in Tank (Top View)", 60, -150., 150., 60, -150., 150.);
  h_topview_truth->GetXaxis()->SetTitle("Z [cm]");
  h_topview_truth->GetYaxis()->SetTitle("X [cm]");

  h_sideview_fit = new TH2D("h_sideview_fit", "Distributions of FITTED Vertices in Tank (Side View)", 60, -150., 150., 60, -150., 150.);
  h_sideview_fit->GetXaxis()->SetTitle("Z [cm]");
  h_sideview_fit->GetYaxis()->SetTitle("Y [cm]");

  h_sideview_truth = new TH2D("h_sideview_truth", "Distributions of TRUTH Vertices in Tank (Side View)", 60, -150., 150., 60, -150., 150.);
  h_sideview_truth->GetXaxis()->SetTitle("Z [cm]");
  h_sideview_truth->GetYaxis()->SetTitle("Y [cm]");

  h_deltaR = new TH1D("h_deltaR", "#Deltar", 75, 0., 150.);
  h_deltaR->GetXaxis()->SetTitle("#Deltar [cm]");

  h_transverse = new TH1D("h_transverse", "transverse distance", 75, 0., 150.);
  h_transverse->GetXaxis()->SetTitle("[cm]");

  h_parallel = new TH1D("h_parallel", "parallel distance", 75, 0., 150.);
  h_parallel->GetXaxis()->SetTitle("[cm]");

  h_deltaR_4pi = new TH1D("h_deltaR_4pi", "#Deltar (density corrected)", 75, 0., 150.);
  h_deltaR_4pi->GetXaxis()->SetTitle("#Deltar");
  h_deltaR_4pi->GetYaxis()->SetTitle("#frac{n}{4#pir^{2}dr} [cm^{-3}]");

  h_true_reco_E = new TH2D("h_true_reco_E", "Reco vs True #mu Energy (in tank?)", 250, 0., 5000., 250, 0., 5000.);
  h_true_reco_E->GetXaxis()->SetTitle("TrueE [MeV]");
  h_true_reco_E->GetYaxis()->SetTitle("RecoE (tank+mrd) [MeV]");

  h_true_reco_Ediff = new TH1D("h_true_reco_Ediff", "Diff in Reco and True #mu Energy", 100, -500, 500);
  h_true_reco_Ediff->GetXaxis()->SetTitle("recoE - trueE [MeV]");

  h_tank_mrd_E = new TH2D("h_tank_mrd_E", "Reco #mu Energy in Tank vs in MRD", 250, 0., 5000., 250, 0., 5000.);
  h_tank_mrd_E->GetXaxis()->SetTitle("Reco MrdE [MeV]");
  h_tank_mrd_E->GetYaxis()->SetTitle("Reco TankE [MeV]");

  h_Ediff_frac_tank = new TH2D("h_Ediff_frac_tank", "#mu Energy Difference (reco - true vs Tank Track Fraction) ", 50, 0., 1., 100, -500, 500);
  h_Ediff_frac_tank->GetXaxis()->SetTitle("tank track fraction");
  h_Ediff_frac_tank->GetYaxis()->SetTitle("E diff (reco-true) [MeV]");

  h_Ediff_frac_mrd = new TH2D("h_Ediff_frac_mrd", "#mu Energy Difference (reco - true vs MRD Track Fraction) ", 50, 0., 1., 100, -500, 500);
  h_Ediff_frac_mrd->GetXaxis()->SetTitle("mrd track fraction");
  h_Ediff_frac_mrd->GetYaxis()->SetTitle("E diff (reco-true) [MeV]");

  h_mrd_eloss_diff = new TH1D("h_mrd_eloss_diff", "diff in MRD energy loss (SB-ANNIE)", 200, -500, 1500);
  h_mrd_eloss_diff->GetXaxis()->SetTitle("sb-annie [MeV]");

  h_mrd_eloss_SBvANNIE = new TH2D("h_mrd_eloss_SBvANNIE", "SB vs ANNIE MRD Energy Loss", 250, 0, 5000, 250, 0, 5000);
  h_mrd_eloss_SBvANNIE->GetXaxis()->SetTitle("annie eloss [MeV]");
  h_mrd_eloss_SBvANNIE->GetYaxis()->SetTitle("sb eloss [MeV]");

  h_true_reco_Ediff_sb = new TH1D("h_true_reco_Ediff_sb", "#mu Energy Difference (reco-true, SB)", 200, -500, 1500);
  h_true_reco_Ediff_sb->GetXaxis()->SetTitle("recoE - trueE [MeV]");

  h_true_reco_Ediff_sb_outerE = new TH1D("h_true_reco_Ediff_sb_outerE", "#mu Energy Difference (reco-true, SB)", 200, -500, 1500); 
  h_true_reco_Ediff_sb_outerE->GetXaxis()->SetTitle("recoE - trueE [MeV]");

  // Graphs
  // total charge at ea vtx
  gr_vtx_charge_in = new TGraph();
  //gr_vtx_charge_in->SetName("gr_vtx_charge_in");
  //gr_vtx_charge_in->SetTitle("Charge as function of vertex Z (inside cone)");
  gr_vtx_charge_in->SetLineColor(8);
  gr_vtx_charge_in->SetLineWidth(2);
  gr_vtx_charge_in->SetMarkerStyle(20);
  gr_vtx_charge_in->SetMarkerColor(8);
  gr_vtx_charge_in->SetFillColor(0);
  gr_vtx_charge_in->GetXaxis()->SetTitle("Z [cm]");
  gr_vtx_charge_in->GetYaxis()->SetTitle("charge [#PE]");

  gr_vtx_charge_out = new TGraph();
  gr_vtx_charge_out->SetLineColor(9);
  gr_vtx_charge_out->SetLineWidth(2);
  gr_vtx_charge_out->SetMarkerStyle(20);
  gr_vtx_charge_out->SetMarkerColor(9);
  gr_vtx_charge_out->SetFillColor(0);
  gr_vtx_charge_out->GetXaxis()->SetTitle("Z [cm]");
  gr_vtx_charge_out->GetYaxis()->SetTitle("charge [#PE]");

  // charge per pmt at ea vtx
  gr_qdensity_in = new TGraph();
  gr_qdensity_in->SetLineColor(8);
  gr_qdensity_in->SetLineWidth(2);
  gr_qdensity_in->SetMarkerStyle(20);
  gr_qdensity_in->SetMarkerColor(8);
  gr_qdensity_in->SetFillColor(0);
  gr_qdensity_in->GetXaxis()->SetTitle("Z [cm]");
  gr_qdensity_in->GetYaxis()->SetTitle("charge per PMT [#PE]");

  gr_qdensity_out = new TGraph();
  gr_qdensity_out->SetLineColor(9);
  gr_qdensity_out->SetLineWidth(2);
  gr_qdensity_out->SetMarkerStyle(20);
  gr_qdensity_out->SetMarkerColor(9);
  gr_qdensity_out->SetFillColor(0);
  gr_qdensity_out->GetXaxis()->SetTitle("Z [cm]");
  gr_qdensity_out->GetYaxis()->SetTitle("charge per PMT [#PE]");

  // total num of pmts at ea vtx
  gr_vtx_detkey_in = new TGraph();
  gr_vtx_detkey_in->SetLineColor(8);
  gr_vtx_detkey_in->SetLineWidth(2);
  gr_vtx_detkey_in->SetMarkerStyle(20);
  gr_vtx_detkey_in->SetMarkerColor(8);
  gr_vtx_detkey_in->SetFillColor(0);
  gr_vtx_detkey_in->GetXaxis()->SetTitle("Z [cm]");
  gr_vtx_detkey_in->GetYaxis()->SetTitle("num PMTs");

  gr_vtx_detkey_out = new TGraph();
  gr_vtx_detkey_out->SetLineColor(9);
  gr_vtx_detkey_out->SetLineWidth(2);
  gr_vtx_detkey_out->SetMarkerStyle(20);
  gr_vtx_detkey_out->SetMarkerColor(9);
  gr_vtx_detkey_out->SetFillColor(0);
  gr_vtx_detkey_out->GetXaxis()->SetTitle("Z [cm]");
  gr_vtx_detkey_out->GetYaxis()->SetTitle("num PMTs");

  gr_effarea_detkey = new TGraph();
  gr_effarea_detkey->SetLineWidth(0);
  gr_effarea_detkey->SetMarkerStyle(24);
  gr_effarea_detkey->SetMarkerColor(38);
  gr_effarea_detkey->SetFillColor(0);
  gr_effarea_detkey->GetXaxis()->SetTitle("detkey");
  gr_effarea_detkey->GetYaxis()->SetTitle("effective PMT area [cm^{2}]");

  gr_fpmt_detkey = new TGraph();
  gr_fpmt_detkey->SetLineWidth(0);
  gr_fpmt_detkey->SetMarkerStyle(24);
  gr_fpmt_detkey->SetMarkerColor(8);
  gr_fpmt_detkey->SetFillColor(0);
  gr_fpmt_detkey->GetXaxis()->SetTitle("detkey");
  gr_fpmt_detkey->GetYaxis()->SetTitle("f");

  gr_effarea_ai = new TGraph();
  gr_effarea_ai->SetLineWidth(0);
  gr_effarea_ai->SetMarkerStyle(25);
  gr_effarea_ai->SetMarkerColor(40);
  gr_effarea_ai->SetFillColor(0);
  gr_effarea_ai->GetXaxis()->SetTitle("a_{i} [cm]");
  gr_effarea_ai->GetYaxis()->SetTitle("effective PMT area [cm^{2}]");

  gr_fpmt_ai = new TGraph();
  gr_fpmt_ai->SetLineWidth(0);
  gr_fpmt_ai->SetMarkerStyle(27);
  gr_fpmt_ai->SetMarkerColor(30);
  gr_fpmt_ai->SetFillColor(0);
  gr_fpmt_ai->GetXaxis()->SetTitle("a_{i} [cm]");
  gr_fpmt_ai->GetYaxis()->SetTitle("f");

  gr_eta_ai = new TGraph();
  gr_eta_ai->SetLineWidth(0);
  gr_eta_ai->SetMarkerStyle(25);
  gr_eta_ai->SetMarkerColor(46);
  gr_eta_ai->SetFillColor(0);
  gr_eta_ai->GetXaxis()->SetTitle("a_{i} [cm]");
  gr_eta_ai->GetYaxis()->SetTitle("#eta (nPE/f)");

  gr_running_avg = new TGraph();
  gr_running_avg->SetLineWidth(0);
  gr_running_avg->SetMarkerStyle(8);
  gr_running_avg->SetMarkerColor(38);
  gr_running_avg->SetFillColor(0);
  gr_running_avg->GetXaxis()->SetTitle("a_{i} [cm]");
  gr_running_avg->GetYaxis()->SetTitle("#eta (nPE/f)");

  gr_qincone_ai = new TGraph();
  gr_qincone_ai->SetLineWidth(0);
  gr_qincone_ai->SetMarkerStyle(26);
  gr_qincone_ai->SetMarkerColor(8);
  gr_qincone_ai->SetFillColor(0);

  gr_qoutcone_ai = new TGraph();
  gr_qoutcone_ai->SetLineWidth(0);
  gr_qoutcone_ai->SetMarkerStyle(32);
  gr_qoutcone_ai->SetMarkerColor(9);
  gr_qoutcone_ai->SetFillColor(0);


  // --------------------------------------------------
  // --- Retrieve Store variables ---------------------
  // --------------------------------------------------
  m_data->CStore.Get("ChannelNumToTankPMTSPEChargeMap", ChannelKeyToSPEMap);  //same for data and mc


  // --------------------------------------------------
  // --- Get Geometry ---------------------------------
  // --------------------------------------------------
  auto get_geo = m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry", geom);
  if (!get_geo)
  {
    Log("MuonFitter Tool: Error retrieving Geometry from ANNIEEvent!", v_error, verbosity);
    return false;
  }
  tank_radius = geom->GetTankRadius();
  tank_height = geom->GetTankHalfheight();
  tank_height /= 2;
  detector_version = geom->GetVersion();
  Log("MuonFitter Tool: Using detector version " + std::to_string(detector_version), v_message, verbosity);
  double barrel_compression = 0.82;
  if (detector_config == "ANNIEp2v6" && !isData) { tank_height *= barrel_compression; }
  else if (isData) { tank_height = 1.2833; }
  if (tank_radius < 1. || isData) { tank_radius = 1.37504; }

  n_tank_pmts = geom->GetNumDetectorsInSet("Tank");
  n_mrd_pmts = geom->GetNumDetectorsInSet("MRD");
  n_veto_pmts = geom->GetNumDetectorsInSet("Veto");
  std::cout<<"  [debug] num Tank/MRD/Veto pmts: "<<n_tank_pmts<<"/"<<n_mrd_pmts<<"/"<<n_veto_pmts<<std::endl;

  Position detector_center = geom->GetTankCentre();
  tank_center_x = detector_center.X();  //[m]
  tank_center_y = detector_center.Y();
  tank_center_z = detector_center.Z();
  std::cout<<"  [debug] tank_center: "<<tank_center_x<<","<<tank_center_y<< ","<<tank_center_z<<std::endl;
  std::cout<<"  [debug] tank_radius: "<<tank_radius<<" m"<<std::endl;
  std::cout<<"  [debug] tank_height: "<<tank_height<<" m"<<std::endl;


  // --------------------------------------------------
  // --- ANNIE in 3D ----------------------------------
  // --------------------------------------------------
  Log("MuonFitter Tool: Creating 3D Geometry", v_debug, verbosity);
  ageom = new TGeoManager("ageom", "ANNIE in 3D");
  TGeoNode *node;

  //material
  vacuum = new TGeoMaterial("vacuum",0,0,0);
  Fe = new TGeoMaterial("Fe",55.845,26,7.87);

  //create media
  Air = new TGeoMedium("Vacuum",0,vacuum);
  Iron = new TGeoMedium("Iron",1,Fe);

  //create volume
  EXPH = ageom->MakeBox("EXPH",Air,300,300,300);
  ageom->SetTopVolume(EXPH);
  ageom->SetTopVisible(0);
  // If you want to see the boundary, input the number 1 instead of 0:
  // geom->SetTopVisible(1);

  //draw CANVAS/TOP CENTER
  bBlock = ageom->MakeSphere("EXPH_vol_center", Iron, 0,3,0,180,0,360);
  bBlock->SetLineColor(1);
  EXPH->AddNodeOverlap(bBlock,N++,new TGeoTranslation(0,0,0));

  //draw TANK
  TGeoVolume *annietank = ageom->MakeTubs("annietank", Iron, 0,tank_radius*100.,tank_height*100.,0,360);  //convert to cm
  annietank->SetLineColor(38);
  //EXPH->AddNodeOverlap(annietank,N++,new TGeoCombiTrans(tank_center_x*100.,tank_center_y*100.,tank_center_z*100.,new TGeoRotation("annietank",0,90,0)));
  EXPH->AddNodeOverlap(annietank,N++,new TGeoCombiTrans(0,0,0,new TGeoRotation("annietank",0,90,0)));
  node = EXPH->GetNode(N-1);


  // --------------------------------------------------
  // --- Read in TANK PMTs ----------------------------
  // --------------------------------------------------
  std::map<std::string, std::map<unsigned long, Detector *>> *Detectors = geom->GetDetectors();

  Log("MuonFitter Tool: Adding tank PMTs to 3D geometry", v_debug, verbosity);
  for (std::map<unsigned long, Detector *>::iterator it = Detectors->at("Tank").begin(); it != Detectors->at("Tank").end(); ++it)
  {
    Detector *apmt = it->second;
    unsigned long detkey = it->first;
    std::string det_type = apmt->GetDetectorType();
    Position position_PMT = apmt->GetDetectorPosition();
    // pmt xyz corrected to (0,0,0)
    x_pmt.insert(std::pair<int,double>(detkey, 100.*(position_PMT.X()-tank_center_x)));
    y_pmt.insert(std::pair<int,double>(detkey, 100.*(position_PMT.Y()-tank_center_y)));
    z_pmt.insert(std::pair<int,double>(detkey, 100.*(position_PMT.Z()-tank_center_z)));
    // pmt xyz w/o tank center correction
    //x_pmt.insert(std::pair<unsigned long, double>(detkey, position_PMT.X()*100.));
    //y_pmt.insert(std::pair<unsigned long, double>(detkey, position_PMT.Y()*100.));
    //z_pmt.insert(std::pair<unsigned long, double>(detkey, position_PMT.Z()*100.));
    // pmt orientation
    Direction direction_PMT = apmt->GetDetectorDirection();
    x_pmt_dir.insert(std::pair<unsigned long, double>(detkey, direction_PMT.X()*100.));
    y_pmt_dir.insert(std::pair<unsigned long, double>(detkey, direction_PMT.Y()*100.));
    z_pmt_dir.insert(std::pair<unsigned long, double>(detkey, direction_PMT.Z()*100.));

    // ANNIE in 3D: Drawing TANK PMTs
    sprintf(blockName,"tank_pmt%lu",detkey);
    bBlock = ageom->MakeSphere(blockName, Iron, 0,1.5,0,180,0,360); //TODO:include PMT type and location condition
    bBlock->SetLineColor(41);
    EXPH->AddNodeOverlap(bBlock,1,new TGeoTranslation(x_pmt[detkey], y_pmt[detkey], z_pmt[detkey]));   //these pos coords work w/ tank center (0,0,0)
    //EXPH->AddNodeOverlap(bBlock,detkey,new TGeoTranslation(position_PMT.X()*100., position_PMT.Y()*100., position_PMT.Z()*100.));   //these pos coords work w/ tank_center_x/y/z*100.
    detkey_to_node.insert(std::pair<unsigned long, int>(detkey, N++));

    // Get PMT areas & efficiency values
    double det_area = 0.;
    if (det_type == "LUX" || det_type == "ANNIEp2v7-glassFaceWCPMT_R7081") { det_area = LUX_AREA; }
    else if (det_type == "ETEL" || det_type == "ANNIEp2v7-glassFaceWCPMT_D784KFLB") { det_area = ETEL_AREA; }
    else if (det_type == "Hamamatsu" || det_type == "ANNIEp2v7-glassFaceWCPMT_R5912HQE") { det_area = HAMAMATSU_AREA; }
    else if (det_type == "Watchboy" || det_type == "ANNIEp2v7-glassFaceWCPMT_R7081") { det_area = WATCHBOY_AREA; }
    else if (det_type == "Watchman" || det_type == "ANNIEp2v7-glassFaceWCPMT_R7081HQE") { det_area = WATCHMAN_AREA; }
    else { Log("MuonFitter Tool: Unrecognized detector type! Setting det_area to 0.", v_error, verbosity); }
    m_pmt_area.insert(std::pair<int, double>(detkey, det_area));

    double eff = 0.25;  //apmt->GetEfficiency(); 
                        //TODO:get effective eff for Cherenkov spectrum
    m_pmt_eff.insert(std::pair<int, double>(detkey, eff));  //unnecessary?

    double alpha = ALPHA_FACTOR * eff * det_area;
    m_pmt_alpha.insert(std::pair<int, double>(detkey, alpha));
    //std::cout<<" [debug] detkey, alpha: "<<detkey<<",  "<<alpha<<std::endl;
    h_alpha->Fill(alpha/1000.);
  }


  // --------------------------------------------------
  // --- Read in MRD PMTs -----------------------------
  // --------------------------------------------------
  Log("MuonFitter Tool: Adding MRD PMTs to 3D geoemtry", v_debug, verbosity);
  for (std::map<unsigned long, Detector *>::iterator it = Detectors->at("MRD").begin(); it != Detectors->at("MRD").end(); ++it)
  {
    Detector *amrdpmt = it->second;
    unsigned long detkey = it->first;
    unsigned long chankey = amrdpmt->GetChannels()->begin()->first;
    Paddle *mrdpaddle = (Paddle *) geom->GetDetectorPaddle(detkey);

    double xmin = mrdpaddle->GetXmin();
    double xmax = mrdpaddle->GetXmax();
    double ymin = mrdpaddle->GetYmin();
    double ymax = mrdpaddle->GetYmax();
    double zmin = mrdpaddle->GetZmin();
    double zmax = mrdpaddle->GetZmax();

    //ANNIE in 3D: drawing MRD
    if (draw3d_mrd)
    {
      Position position_MRD = mrdpaddle->GetOrigin();
      sprintf(blockName,"mrd_pmt%lu",detkey);
      bBlock = ageom->MakeBox(blockName, Iron, (xmax-xmin)*100., (ymax-ymin)*100., (zmax-zmin)*100.); //TODO:include PMT orientation
      bBlock->SetLineColor(16);
      EXPH->AddNodeOverlap(bBlock,detkey,new TGeoTranslation((position_MRD.X()-tank_center_x)*100., (position_MRD.Y()-tank_center_y)*100., (position_MRD.Z()-tank_center_z)*100.));
      detkey_to_node.insert(std::pair<unsigned long, int>(detkey, N++));

      if (verbosity > v_debug)
      {
        std::cout << " [debug] position_MRD (no tank_center correction): " << position_MRD.X() << "," << position_MRD.Y() << "," << position_MRD.Z() << std::endl;
        std::cout << " [debug] blockName: " << blockName << std::endl;
        std::cout << " [debug] detkey_to_node[detkey]: " << detkey_to_node[detkey] << std::endl;
      }
    }
  }

  // --------------------------------------------------
  // --- Read in FMV PMTs -----------------------------
  // --------------------------------------------------
  Log("MuonFitter Tool: Adding FMV PMTs to 3D geoemtry", v_debug, verbosity);
  for (std::map<unsigned long, Detector *>::iterator it = Detectors->at("Veto").begin(); it != Detectors->at("Veto").end(); ++it)
  {
    Detector *avetopmt = it->second;
    unsigned long detkey = it->first;
    unsigned long chankey = avetopmt->GetChannels()->begin()->first;
    Paddle *vetopaddle = (Paddle *) geom->GetDetectorPaddle(detkey);

    double xmin = vetopaddle->GetXmin();
    double xmax = vetopaddle->GetXmax();
    double ymin = vetopaddle->GetYmin();
    double ymax = vetopaddle->GetYmax();
    double zmin = vetopaddle->GetZmin();
    double zmax = vetopaddle->GetZmax();

    //ANNIE in 3D: drawing FMV
    if (draw3d_fmv)
    {
      Position position_FMV = vetopaddle->GetOrigin();
      sprintf(blockName,"fmv_pmt%lu",detkey);
      bBlock = ageom->MakeBox(blockName, Iron, (xmax-xmin)*100., (ymax-ymin)*100., (zmax-zmin)*100.); //TODO:include PMT orientation
      bBlock->SetLineColor(20);
      EXPH->AddNodeOverlap(bBlock,detkey,new TGeoTranslation((position_FMV.X()-tank_center_x)*100., (position_FMV.Y()-tank_center_y)*100., (position_FMV.Z()-tank_center_z)*100.));
      detkey_to_node.insert(std::pair<unsigned long, int>(detkey, N++));

      if (verbosity > v_debug)
      {
        std::cout << " [debug] position_FMV (x,y,z): " << position_FMV.X() << "," << position_FMV.Y() << "," << position_FMV.Z() << std::endl;
        std::cout << " [debug] blockName: " << blockName << std::endl;
        std::cout << " [debug] detkey_to_node[detkey]: " << detkey_to_node[detkey] << std::endl;
      }
    }
  }

  //ANNIE in 3D: set max nodes
  maxN = N;
  Log("MuonFitter Tool: Number of nodes in 3D geometry: " + std::to_string(maxN), v_debug, verbosity);


  // Save 3D plots
  ageom->CloseGeometry();   //close 3d geometry
  EXPH->SetVisibility(0);
  if (plot3d)
  {
    canvas_3d = new TCanvas("canvas_3d", "3D Event Display", 800, 600);
    canvas_3d->cd(1);
    EXPH->Draw();
    canvas_3d->Modified();
    canvas_3d->Update();
  }

  // Save start & stop vertices to file; Save (ai,eta) values
  //std::string pos_fname = "posFile.txt";
  //pos_file.open(pos_fname.c_str());
  //pos_file << "##evnum,startX,startY,startZ,stopX,stopY,stopZ" << std::endl;
  std::string pos_fname = "ev_ai_eta.txt";
  pos_file.open(pos_fname.c_str(), std::ostream::app);
  if (!pos_file)
  {
    std::cout << "File " << pos_fname << " does not exist! Creating now.." << std::endl;
    pos_file.open(pos_fname.c_str());
    pos_file << "##ev_id,ai,eta" << std::endl;
  }

  // Save charge data to file
  //std::string cpp_fname = "charge_per_pmt.txt";
  //cpp_file.open(cpp_fname.c_str());
  //cpp_file << "##evnum,nVtx,Qin,Qout,avgSumIn,avgSumOut,avgQin,avgQout" << std::endl;

  std::string pehits_fname = "tot_pe_hits.txt";
  pehits_file.open(pehits_fname.c_str(), std::ostream::app);
  if (!pehits_file)
  {
    std::cout << "File " << pehits_fname << " does not exist! Creating now..." << std::endl;
    pehits_file.open(pehits_fname.c_str());
    pehits_file << "##partfile,main_cluster_hits,main_cluster_charge" << std::endl;
  }

  // Save avg eta to left and right of true tank track length
  if (!isData && display_truth)
  {
    truetrack_file.open("true_track_len.txt", std::ostream::app);
    if (!truetrack_file)
    {
      std::cout << " [debug] true_track_len.txt does not exist! Creating now..." << std::endl;
      truetrack_file.open("true_track_len.txt");
      truetrack_file << "#event_id,true_track_len,left_avg,right_avg" << std::endl;
    }
  }

  // Save fitted_track_len, nhits, nhits_incone to file
  nhits_trlen_file.open("nhits_trlen.txt", std::ostream::app);
  if (!nhits_trlen_file)
  {
    std::cout << "File nhits_trlen.txt does not exist! Creating now..." << std::endl;
    nhits_trlen_file.open("nhits_trlen.txt");
    nhits_trlen_file << "##event_id,track_fit,nhits,nhits_incone,totalpe,totalpe_incone" << std::endl;
  }

  // Load vertex fits
  if (fit_mode)
  {
    if (this->FileExists(tankTrackFitFile))
    {
      // Load event part, number and track fit
      std::cout << "MuonFitter Tool: Loading tank track fits" << std::endl;
      this->LoadTankTrackFits();
    }
    else
    {
      Log("MuonFitter Tool: File for tank track fit does not exist! Continuing tool not in fit mode.", v_error, verbosity);
      fit_mode = false;
    }
  }

  if (verbosity > 2) std::cout << "MuonFitter Tool: Initialization complete" << std::endl;

  return true;
}


bool MuonFitter::Execute(){
  Log("MuonFitter Tool: Executing", v_debug, verbosity);

  int eventExists = m_data->Stores.count("ANNIEEvent");
  if (!eventExists)
  {
    Log("MuonFitter Tool: No ANNIEEvent store!", v_error, verbosity);
    return false;
  }

  // --------------------------------------------------
  // --- Reset variables ------------------------------
  // --------------------------------------------------
  //this->ResetVariables();   //TODO:write ResetVariables func
  bool drawEvent = false;
  if (plot3d) { reset_3d(); }

  gr_vtx_charge_in->Set(0);
  gr_vtx_charge_out->Set(0);
  gr_qdensity_in->Set(0);
  gr_qdensity_out->Set(0);
  gr_effarea_detkey->Set(0);
  gr_fpmt_detkey->Set(0);
  gr_effarea_ai->Set(0);
  gr_fpmt_ai->Set(0);
  gr_eta_ai->Set(0);
  gr_running_avg->Set(0);

  if (h_tzero) h_tzero->Delete();
  h_tzero = new TH1D("h_tzero", "t0?", 60, -15, 15);
  h_tzero->GetXaxis()->SetTitle("[ns]");

  avg_eta = 0;
  left_avg_eta = 0.;
  right_avg_eta = 0.;
  num_left_eta = 0;
  num_right_eta = 0;

  m_data->CStore.Set("FittedTrackLengthInWater", -888.);
  Position dummy_vtx(-888, -888, -888);
  m_data->CStore.Set("FittedMuonVertex", dummy_vtx);

  int get_ok;

  // --------------------------------------------------
  // --- Get event info (data) ------------------------
  // --------------------------------------------------
  if (isData)
  {
    get_ok = m_data->Stores["ANNIEEvent"]->Get("EventNumber", evnum);
    std::cout << "MuonFitter Tool: Working on event " << evnum << std::endl;
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving EventNumber from ANNIEEvent!", v_error, verbosity); return false; }
    get_ok = m_data->Stores["ANNIEEvent"]->Get("RunNumber", runnumber);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving RunNumber from ANNIEEvent!", v_error, verbosity); return false; }
    get_ok = m_data->Stores["ANNIEEvent"]->Get("PartNumber", partnumber);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving PartNumber from ANNIEEvent!", v_error, verbosity); partnumber = -1; }
  }

  // --------------------------------------------------
  // --- Get event info (MC) ------------------------
  // --------------------------------------------------
  TVector3 trueTrackDir;
  if (!isData)
  {
    get_ok = m_data->Stores["ANNIEEvent"]->Get("MCParticles", mcParticles);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving MCParticles from ANNIEEvent!", v_error, verbosity); return false; }  //needed to retrieve true vertex and direction
    get_ok = m_data->Stores["ANNIEEvent"]->Get("MCEventNum", mcevnum);
    std::cout << "MuonFitter Tool: Working on event " << mcevnum << std::endl;
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving MCEventNum from ANNIEEvent!", v_error, verbosity); return false; }
    get_ok = m_data->Stores["ANNIEEvent"]->Get("MCTriggernum", mctrignum);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving MCTriggernum from ANNIEEvent!", v_error, verbosity); return false; }
    get_ok = m_data->Stores["ANNIEEvent"]->Get("MCFile", mcFile);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving MCFile from ANNIEEvent!", v_error, verbosity); mcFile = "-1"; }
    std::string delim = ".";
    std::string tmp_str = mcFile.erase(0, mcFile.find(delim) + delim.length());
    partnumber = stoi(tmp_str.substr(0, tmp_str.find(delim)));
    //std::cout << " [debug] mcFile, partnumber: " << mcFile << ", " << partnumber << std::endl;

    // Get RecoEvent variables
    RecoVertex *truevtx = 0;
    get_ok = m_data->Stores["RecoEvent"]->Get("TrueVertex", truevtx);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving TrueVertex from RecoEvent Store!", v_error, verbosity); return false; }
    trueVtxX = truevtx->GetPosition().X();
    trueVtxY = truevtx->GetPosition().Y();
    trueVtxZ = truevtx->GetPosition().Z();
    trueVtxTime = truevtx->GetTime();
    trueDirX = truevtx->GetDirection().X();
    trueDirY = truevtx->GetDirection().Y();
    trueDirZ = truevtx->GetDirection().Z();
    trueTrackDir = TVector3(trueDirX,trueDirY,trueDirZ).Unit();
    //TODO:double check true track direction doesn't change when .Unit()
    //is it only the muon direction? what about neutrino, pion directions/vertices?

    double trueAngleRad = TMath::ACos(trueDirZ);
    trueAngle = trueAngleRad/(TMath::Pi()/180.);
    h_truevtx_angle->Fill(trueAngle);

    RecoVertex *truestopvtx = 0;
    get_ok = m_data->Stores["RecoEvent"]->Get("TrueStopVertex", truestopvtx);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving TrueStopVertex from RecoEvent Store!", v_error, verbosity); return false; }
    trueStopVtxX = truestopvtx->GetPosition().X();
    trueStopVtxY = truestopvtx->GetPosition().Y();
    trueStopVtxZ = truestopvtx->GetPosition().Z();

    get_ok = m_data->Stores["RecoEvent"]->Get("TrueTrackLengthInWater", trueTrackLengthInWater);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving TrueTrackLengthInWater from RecoEvent Store!", v_error, verbosity); trueTrackLengthInWater = -1; }
    trueTrackLengthInWater = trueTrackLengthInWater*100.;   //convert to cm
    h_true_tanktrack_len->Fill(trueTrackLengthInWater);

    get_ok = m_data->Stores["RecoEvent"]->Get("TrueTrackLengthInMRD", trueTrackLengthInMRD);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving TrueTrackLengthInMRD from RecoEvent Store!", v_error, verbosity); trueTrackLengthInMRD = -1; }

    get_ok = m_data->Stores["RecoEvent"]->Get("TrueMuonEnergy", trueMuonEnergy);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving TrueMuonEnergy from RecoEvent Store!", v_error, verbosity); trueMuonEnergy = -1; }

    get_ok = m_data->Stores["RecoEvent"]->Get("NRings", nrings);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving NRings, true from RecoEvent!", v_error, verbosity); }
    get_ok = m_data->Stores["RecoEvent"]->Get("IndexParticlesRing", particles_ring);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving IndexParticlesRing, true from RecoEvent!", v_error, verbosity); }
  }

  // --------------------------------------------------
  // --- Get event info (both) ------------------------
  // --------------------------------------------------
  uint32_t trigword;
  get_ok = m_data->Stores.at("ANNIEEvent")->Get("TriggerWord",trigword);
  std::cout << "trigword: " << trigword << std::endl;

  // use cosmic triggers to get charge limit on PMTs
  if (trigword == 36)
  {
    //load hits
    if (isData)
    {
      get_ok = m_data->Stores["ANNIEEvent"]->Get("Hits", m_hits);
      if (!get_ok) std::cout << "MuonFitter Tool: Unable to retrieve hits for this cosmic trigger. "<< std::endl;
      std::cout << " m_hits size: " << m_hits->size() << std::endl;

      for (std::pair<unsigned long, std::vector<Hit>> &&apair : *m_hits)
      {
        //for each detkey and its vector of hits
        unsigned long chankey = apair.first;
        std::vector<Hit> &v_hits = apair.second;

        std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(chankey);
        if (it != ChannelKeyToSPEMap.end())
        {
          Detector* this_detector = geom->ChannelToDetector(chankey);
          unsigned long detkey = this_detector->GetDetectorID();  //chankey same as detkey
          for (Hit &ahit : v_hits)
          {
            double hit_charge = ahit.GetCharge();
            double hit_PE = hit_charge / ChannelKeyToSPEMap.at(chankey);
            h_charge_detkey->Fill(detkey, hit_PE);

            //look into large PE hits
            if (hit_PE > 2000.) 
            {
              std::cout << " [debug] found high PE hit! p" << partnumber;
              if (isData) std::cout << evnum;
              else std::cout << mcevnum;
              std::cout << " detkey " << detkey << " with " << hit_PE << "p.e." << std::endl;
            }
          }
        }
      } //end loop thru m_hits
    } //end if isData
    else  //isMC
    {
      std::cout << "not cosmic!" << std::endl;
    }
  }

  // skip entry if not beam (5)
  //if (trigword != 5) return true;


  // --------------------------------------------------
  // --- Check for particles other than muon (MC) -----
  // --------------------------------------------------
  bool hasPion = false;
  int n_rings = 0;
  if (!isData)
  {
    //if (mcParticles->size() > 1) return false;
    for (unsigned int mcp_i = 0; mcp_i < mcParticles->size(); mcp_i++)
    {
      MCParticle aparticle = mcParticles->at(mcp_i);
      //std::cout << " [debug] Ev " << mcevnum << ", particle pdg code: " << aparticle.GetPdgCode() << std::endl;
      if (std::find(particles_ring.begin(), particles_ring.end(), mcp_i) != particles_ring.end())
      {
        ++n_rings;
        //if (aparticle.GetPdgCode() == 22 | aparticle.GetPdgCode() == -22) hasPion = true;
      }
    }
    std::cout << " [debug] n_rings (my counter): " << n_rings << std::endl;
    if (n_rings > 1) hasPion = true;
  }
  //if (hasPion) return true;   //skip event if there is pion
  if (hasPion) std::cout << " [debug] has pion!" << std::endl;    //just indicate that there is pion


  // --------------------------------------------------
  // --- Check for FMV hits ---------------------------
  // --------------------------------------------------
  bool hasVeto = false;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("TDCData", tdcdata);   //named the same whether it's data or MC
  if (!get_ok) { Log("MuonFitter Tool: No TDCData object in ANNIEEvent! Abort!", v_error, verbosity); return true; }

  if (tdcdata->size() > 0)
  {
    Log("MuonFitter Tool: Looping over FMV/MRD hits... looking for Veto activity", v_debug, verbosity);
    for (auto&& anmrdpmt : (*tdcdata))
    {
      unsigned long chankey = anmrdpmt.first;
      Detector *thedetector = geom->ChannelToDetector(chankey);
      unsigned long detkey = thedetector->GetDetectorID();
      if (thedetector->GetDetectorElement() == "Veto") { hasVeto = true;}
    }
  }
  if (hasVeto)
  {
    Log("MuonFitter Tool: Found FMV/Veto hit!", v_debug, verbosity);
    m_data->CStore.Set("DrawEventDisplay", drawEvent);
    return true;
  }


  // --------------------------------------------------
  // --- Check for MRD tracks -------------------------
  // --------------------------------------------------
  //get_ok = m_data->Stores["ANNIEEvent"]->Get("MRDTriggerType", mrdTriggerType);  //XXX:care about this?
  get_ok = m_data->Stores["MRDTracks"]->Get("NumMrdTracks", numTracksInEv);
  get_ok = m_data->Stores["MRDTracks"]->Get("MRDTracks", mrdTracks);    //XXX:might need MC version

  if (!get_ok) { Log("MuonFitter Tool: Couldn't retrieve MRD tracks info. Did you run TimeClustering/FindMRDTracks first?", v_debug, verbosity); return false; }

  //skip if more than 1 track in MRD
  if (numTracksInEv != 1)
  {
    Log("MuonFitter Tool: More than 1 reconstructed track found!", v_debug, verbosity);
    m_data->CStore.Set("DrawEventDisplay", drawEvent);
    return true;
  }


  // --------------------------------------------------
  // --- Get MRD track params -------------------------
  // --------------------------------------------------
  double reco_mrd_track = 0.;
  bool isMrdStopped;
  bool isMrdPenetrating;
  bool isMrdSideExit;
  for(int track_i = 0; track_i < numTracksInEv; track_i++)
  {
    BoostStore* thisTrackAsBoostStore = &(mrdTracks->at(track_i));
    
    //get track properties needed for through-going muon selection
    thisTrackAsBoostStore->Get("StartVertex", mrdStartVertex);
    thisTrackAsBoostStore->Get("StopVertex", mrdStopVertex);
    thisTrackAsBoostStore->Get("TrackAngle", trackAngle);
    thisTrackAsBoostStore->Get("TrackAngleError", trackAngleError);
    thisTrackAsBoostStore->Get("PenetrationDepth", penetrationDepth);
    thisTrackAsBoostStore->Get("MrdEntryPoint", mrdEntryPoint);
    //thisTrackAsBoostStore->Get("NumLayersHit", numLayersHit);
    thisTrackAsBoostStore->Get("LayersHit", numLayersHit);
    thisTrackAsBoostStore->Get("EnergyLoss", energyLoss);
    thisTrackAsBoostStore->Get("EnergyLossError", energyLossError);
    thisTrackAsBoostStore->Get("StartTime", mrdStartTime);
    thisTrackAsBoostStore->Get("TankExitPoint", tankExitPoint);
    thisTrackAsBoostStore->Get("IsMrdStopped", isMrdStopped);
    thisTrackAsBoostStore->Get("IsMrdPenetrating", isMrdPenetrating);
    thisTrackAsBoostStore->Get("IsMrdSideExit", isMrdSideExit);

    reco_mrd_track = sqrt(pow((mrdStopVertex.X()-mrdStartVertex.X()),2)+pow(mrdStopVertex.Y()-mrdStartVertex.Y(),2)+pow(mrdStopVertex.Z()-mrdStartVertex.Z(),2));
    reco_mrd_track = reco_mrd_track*100.;       //convert to cm
    penetrationDepth = penetrationDepth*100.;   //convert to cm

    std::cout << "  [debug] mrdStartVertex: " << mrdStartVertex.X() << "," << mrdStartVertex.Y() << "," << mrdStartVertex.Z() << std::endl;
    std::cout << "  [debug] mrdStopVertex: " << mrdStopVertex.X() << "," << mrdStopVertex.Y() << "," << mrdStopVertex.Z() << std::endl;
    std::cout << "  [debug] tankExitPoint: " << tankExitPoint.X() << "," << tankExitPoint.Y() << "," << tankExitPoint.Z() << std::endl;
    std::cout << "  [debug] mrdEntryPoint: " << mrdEntryPoint.X() << "," << mrdEntryPoint.Y() << "," << mrdEntryPoint.Z() << std::endl;
    std::cout << "  [debug] trackAngle: " << trackAngle << std::endl;
    std::cout << "  [debug] numLayersHit: " << numLayersHit << std::endl;
    std::cout << "  [debug] penetrationDepth: " << penetrationDepth << std::endl;
    std::cout << "  [debug] reco_mrd_track: " << reco_mrd_track << std::endl;
    std::cout << "  [debug] energyLoss: " << energyLoss << std::endl;
    std::cout << "  [debug] mrdStartTime: " << mrdStartTime << std::endl;
    std::cout << "  [debug] isMrdStopped: " << isMrdStopped << std::endl;
    std::cout << "  [debug] isMrdPenetrating: " << isMrdPenetrating << std::endl;
    std::cout << "  [debug] isMrdSideExit: " << isMrdSideExit << std::endl;

    //save MRD start/stop coordinates to file
    //pos_file << evnum << "," << mrdStartVertex.X() << ","
    //                         << mrdStartVertex.Y() << ","
    //                         << mrdStartVertex.Z() << ","
    //                         << mrdStopVertex.X() << ","
    //                         << mrdStopVertex.Y() << ","
    //                         << mrdStopVertex.Z() << std::endl;
    h_num_mrd_layers->Fill(numLayersHit);
  }
  //XXX:check whether MRD track falls into the fiducial volume

  double tankExitPointX = 100.*(tankExitPoint.X()-tank_center_x);
  double tankExitPointY = 100.*(tankExitPoint.Y()-tank_center_y);
  double tankExitPointZ = 100.*(tankExitPoint.Z()-tank_center_z);
  TVector3 tankExit(tankExitPointX, tankExitPointY, tankExitPointZ);

  //#TODO:or use mrdStartVertex?
  double mrdEntryPointX = 100.*(mrdEntryPoint.X()-tank_center_x);
  double mrdEntryPointY = 100.*(mrdEntryPoint.Y()-tank_center_y);
  double mrdEntryPointZ = 100.*(mrdEntryPoint.Z()-tank_center_z);

  // Create vectors for MRD start/stop
  TVector3 mrdStart((mrdStartVertex.X()-tank_center_x)*100., (mrdStartVertex.Y()-tank_center_y)*100., (mrdStartVertex.Z()-tank_center_z)*100.);   //[cm]
  TVector3 mrdStop((mrdStopVertex.X()-tank_center_x)*100., (mrdStopVertex.Y()-tank_center_y)*100., (mrdStopVertex.Z()-tank_center_z)*100.);
  TVector3 mrdTrackDir = (mrdStop - mrdStart).Unit();

  //add mrdStartVertex to 3D geo
  if (isData) sprintf(blockName,"mrdStart_%u",evnum);
  else sprintf(blockName,"mrdStart_%u",mcevnum);
  bBlock = ageom->MakeSphere(blockName, Iron, 0,1.5,0,180,0,360);
  bBlock->SetLineColor(8);
  EXPH->AddNodeOverlap(bBlock,69,new TGeoTranslation(mrdStart.X(), mrdStart.Y(), mrdStart.Z()));
  N++;

  //add mrdStopVertex to 3D geo
  if (isData) sprintf(blockName,"mrdStop_%u",evnum);
  else sprintf(blockName,"mrdStop_%u",mcevnum);
  bBlock = ageom->MakeSphere(blockName, Iron, 0,1.5,0,180,0,360);
  bBlock->SetLineColor(2);
  EXPH->AddNodeOverlap(bBlock,69,new TGeoTranslation(mrdStop.X(), mrdStop.Y(), mrdStop.Z()));
  N++;

  // Select for tracks that end in MRD
  if (!isMrdStopped)
  {
    std::cout << " [debug] MUON DID NOT STOP IN MRD. Skipping..." << std::endl;
    return true;
  }

  if (isMrdSideExit)
  {
    std::cout << " [debug] MUON EXITTED SIDE OF MRD. Skipping..." << std::endl;
    //return true;
  }



  // --------------------------------------------------
  // --- Calculate several vertex candidates ----------
  // --------------------------------------------------
  std::vector<TVector3> vtxCandidates;
  int c = 1;
  bool outsideTank = false;
  bool inFV = false;
  while (!outsideTank)
  {
    TVector3 v = mrdStart - c*10.*mrdTrackDir;  //[cm]
    if (c <= 5) { vtxCandidates.push_back(v); }
    else if (c > 5)
    {
      //check that vtx is contained in tank cylinder
      if ((pow(v.X(),2) + pow(v.Z(),2) >= pow(tank_radius*100.,2)) || (v.Y() > tank_height*100.) || (v.Y() < -tank_height*100.))
      {
        outsideTank = true;
      }
      else { vtxCandidates.push_back(v); }

      //check if any vtx intersects with Fiducial Volume around center
      if (pow(v.X(),2) + pow(v.Y(),2) + pow(v.Z(),2) <= pow(tank_radius*100.-50.,2))
      {
        //std::cout << " [debug] (distance from tank center)^2: " << pow(v.X(),2) + pow(v.Y(),2) + pow(v.Z(),2) << std::endl;
        inFV = true;
      }
    }

    // add candidate vertex to 3D geo
/*    sprintf(blockName,"tankVtx_%d", c);
    bBlock = ageom->MakeSphere(blockName, Iron, 0,1.5,0,180,0,360);
    bBlock->SetLineColor(6);
    EXPH->AddNodeOverlap(bBlock,69,new TGeoTranslation(v.X(), v.Y(), v.Z()));
    N++; */
    c++;
  } //done calculating vtx candidates
  // make sure track originates in tank
  if (!inFV) return true;



  std::cout << "MuonFitter Tool: Num of vtxCandidates: " << vtxCandidates.size() << std::endl;
  double l_zmax = vtxCandidates.at(0).Z();
  double l_zmin = vtxCandidates.at(vtxCandidates.size()-1).Z();

  //check distance btwn vertices
  for (int i = 1; i < vtxCandidates.size(); ++i)
  {
    TVector3 v0(vtxCandidates.at(i-1));
    TVector3 v1(vtxCandidates.at(i));
    double d = TMath::Sqrt(pow(v1.X()-v0.X(),2) + pow(v1.Y()-v0.Y(),2) + pow(v1.Z()-v0.Z(),2));
    //std::cout << "  [debug] distance btwn candidate vtx: " << d << std::endl;
  }
  if (!isData)
  {
    h_lastvtx_z->Fill(vtxCandidates.at(vtxCandidates.size()-1).Z());
    //h_truevtx_z->Fill(trueVtxZ);
  }


  // --------------------------------------------------
  // --- Check short tracks ---------------------------
  // --------------------------------------------------
  if (!isData)
  {
    // skip if true vertex is downstream of tank
    if (trueVtxZ > 0.)
    {
      std::cout << " SKIPPING event. True vertex is in downstream half of tank. " << mcevnum << std::endl;
      //return true;
    }

    TVector3 trueTankTrack = TVector3(tankExitPointX, tankExitPointY, tankExitPointZ) - TVector3(trueVtxX, trueVtxY, trueVtxZ);
    // skip if tank track is short
    if (trueTrackLengthInWater < 100.)
    {
      std::cout << " Track less than 100 cm! SKIPPING! " << mcevnum << std::endl;
      return true;
    }

    double trueTankTrackAngle = trueTankTrack.Angle(TVector3(0,0,1));
    std::cout << " [debug] trueTankTrackAngle: " << trueTankTrackAngle*180./TMath::Pi() << std::endl;
  }


  // For each event, load clusters and their hits
  //   (code from PhaseIITreeMaker)
  //   (this contains x,y,z info of hit (PMT) location
  // XXX: maybe load event hits instead?
  // XXX: make cuts in num of clusters?

  // --------------------------------------------------
  // --- Load Clusters --------------------------------
  // --------------------------------------------------
  bool has_clusters = false;
  if (isData) { has_clusters = m_data->CStore.Get("ClusterMap", m_all_clusters); }
  else { has_clusters = m_data->CStore.Get("ClusterMapMC", m_all_clusters_MC); }
  if (!has_clusters)
  {
    std::cout << "MuonFitter Tool: No clusters found in CStore! Did you run ClusterFinder tool before this?" << std::endl;
    return false;
  }

  has_clusters = m_data->CStore.Get("ClusterMapDetkey", m_all_clusters_detkeys);  //same for data and mc
  if (!has_clusters)
  {
    std::cout << "MuonFitter Tool: No cluster detkeys found in CStore! Did you run ClusterFinder tool before this?" << std::endl;
    return false;
  }

  Log("MuonFitter Tool: Accessing pairs in all_clusters map", v_debug, verbosity);
  int cluster_num = 0;
  int cluster_size = 0;
  if (isData) { cluster_size = (int) m_all_clusters->size(); }
  else { cluster_size = (int) m_all_clusters_MC->size(); }
  std::cout << " [debug] cluster_size (num clusters in event): " << cluster_size << std::endl;
  //if (cluster_size <= 1) { return false; }    // XXX: test single cluster hit positions (!=1)

  std::map<double, std::vector<Hit>>::iterator it_cluster_pair;
  std::map<double, std::vector<MCHit>>::iterator it_cluster_pair_MC;
  if (isData) { it_cluster_pair = (*m_all_clusters).begin(); }
  else { it_cluster_pair_MC = (*m_all_clusters_MC).begin(); }


  double avg_charge_per_pmt_in = 0., avg_charge_per_pmt_out = 0.;
  int max_vertex = -1;
  double max_vtx_charge = 0.;

  // Find the main cluster (max charge and in [0..2000ns] time window)
  double main_cluster_time = 0;
  double main_cluster_charge = 0;
  double main_cluster_hits = 0;

  bool found_muon = false;
  double earliest_hittime = 0;
  std::vector<double> v_cluster_times;
  std::map<unsigned long, double> m_cluster_charge;
  if (isData)
  {
    for (std::pair<double, std::vector<Hit>>&& apair : *m_all_clusters)
    {
      m_cluster_charge.clear();
      for (unsigned long detkey = 332 ; detkey < 464; detkey++) { m_cluster_charge.emplace(detkey, 0.); }

      std::vector<Hit>& cluster_hits = apair.second;
      double temp_time = 0;
      double temp_charge = 0;
      int temp_hits = 0;
      std::cout << " [debug] cluster_hits.size (num hits in cluster): " << cluster_hits.size() << std::endl;
      for (int ihit = 0; ihit < cluster_hits.size(); ihit++)
      {
        temp_hits++;
        temp_time += cluster_hits.at(ihit).GetTime();
        temp_charge += cluster_hits.at(ihit).GetCharge();

        //XXX: check single cluster hit positions
        int chankey = cluster_hits.at(ihit).GetTubeId();
        std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(chankey);
        if (it != ChannelKeyToSPEMap.end())
        {
          Detector* this_detector = geom->ChannelToDetector(chankey);
          unsigned long detkey = this_detector->GetDetectorID();  //chankey same as detkey
          Position det_pos = this_detector->GetDetectorPosition();
          h_clusterhit_x->Fill(100.*(det_pos.X()-tank_center_x));
          h_clusterhit_y->Fill(100.*(det_pos.Y()-tank_center_y));
          h_clusterhit_z->Fill(100.*(det_pos.Z()-tank_center_z));
          h_clusterhit_detkey->Fill(detkey);
          //m_cluster_charge[detkey] += cluster_hits.at(ihit).GetCharge();
        }
      }
      if (temp_hits > 0) temp_time /= temp_hits;  //mean time
      std::cout << " [debug] temp_time: " << temp_time << std::endl;
      if (temp_time > 2000. || temp_time < 0.) continue;  // not in time window
      if (temp_charge > main_cluster_charge)
      {
        found_muon = true;
        main_cluster_charge = temp_charge;
        main_cluster_time = apair.first;    //same as mean time above
        main_cluster_hits = cluster_hits.size();
      }
    }
  }
  else  //is MC
  {
    for (std::pair<double, std::vector<MCHit>>&& apair : *m_all_clusters_MC)
    {
      m_cluster_charge.clear();
      for (unsigned long detkey = 332 ; detkey < 464; detkey++) { m_cluster_charge.emplace(detkey, 0.); }
      std::vector<MCHit>& cluster_hits_MC = apair.second;
      double temp_time = 0;
      double temp_charge = 0;
      int temp_hits = 0;
      std::vector<double> v_cluster_hits;

      std::cout << " [debug] cluster_hits_MC.size (num hits in cluster): " << cluster_hits_MC.size() << std::endl;

      for (int ihit = 0; ihit < cluster_hits_MC.size(); ihit++)
      {
        temp_hits++;
        //std::cout << " [debug] hit time: " << cluster_hits_MC.at(ihit).GetTime() << std::endl;
        v_cluster_hits.push_back(cluster_hits_MC.at(ihit).GetTime());
        h_clusterhit_time->Fill(cluster_hits_MC.at(ihit).GetTime());
        temp_time += cluster_hits_MC.at(ihit).GetTime();
        temp_charge += cluster_hits_MC.at(ihit).GetCharge();

        //XXX: check single cluster hit positions
        int chankey = cluster_hits_MC.at(ihit).GetTubeId();
        std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(chankey);
        if (it != ChannelKeyToSPEMap.end())
        {
          Detector* this_detector = geom->ChannelToDetector(chankey);
          unsigned long detkey = this_detector->GetDetectorID();  //chankey same as detkey
          Position det_pos = this_detector->GetDetectorPosition();
          h_clusterhit_x->Fill(100.*(det_pos.X()-tank_center_x));
          h_clusterhit_y->Fill(100.*(det_pos.Y()-tank_center_y));
          h_clusterhit_z->Fill(100.*(det_pos.Z()-tank_center_z));
          h_clusterhit_detkey->Fill(detkey);
          m_cluster_charge[detkey] += cluster_hits_MC.at(ihit).GetCharge();
        }
      }
      sort(v_cluster_hits.begin(), v_cluster_hits.end());
      std::cout << " [debug] all cluster hit times (sorted): ";
      for (int ct = 0; ct < v_cluster_hits.size(); ++ct)
      {
        std::cout << v_cluster_hits.at(ct) << ",";
      }
      std::cout << std::endl;
      h_clusterhit_timespread->Fill(v_cluster_hits.at(v_cluster_hits.size()-1)-v_cluster_hits.at(0));
      if (temp_hits > 0) temp_time /= temp_hits;  //mean time
      v_cluster_times.push_back(temp_time);
      if (temp_time > 2000. || temp_time < 0.) continue;  // not in time window
      if (temp_charge > main_cluster_charge)
      {
        found_muon = true;
        main_cluster_charge = temp_charge;
        main_cluster_time = apair.first;    //same as mean time above
        //main_cluster_hits = cluster_hits_MC.size();
        earliest_hittime = v_cluster_hits.at(0);
        std::cout << " [debug] earliest_hittime: " << earliest_hittime << std::endl;

        main_cluster_hits = 0;
        for (unsigned long detkey = 332; detkey < 464; detkey++)
        {
          if (m_cluster_charge[detkey] > 0.)
          {
            main_cluster_hits++;
          }
        }

      }
    }
  }
  std::cout << " [debug] all cluster times: ";
  for (int ct = 0; ct < v_cluster_times.size(); ++ct)
  {
    std::cout << v_cluster_times.at(ct) << ",";
  }
  std::cout << std::endl;
  std::cout << " [debug] main_cluster_time: " << main_cluster_time << std::endl;
  std::cout << " [debug] main_cluster_charge: " << main_cluster_charge << std::endl;
  std::cout << " [debug] main_cluster_hits: " << main_cluster_hits << std::endl;

  // --------------------------------------------------
  // --- Check coincidence with MRD -------------------
  // --------------------------------------------------
  double tankmrd_tdiff = mrdStartTime - main_cluster_time;
  std::cout << " [debug] tankmrd_tdiff: " << tankmrd_tdiff << std::endl;
  if ((tankmrd_tdiff > PMTMRDOffset-50) && (tankmrd_tdiff < PMTMRDOffset+50))
  {
    std::cout << "MuonFitter tool: Main cluster is coincident with MRD cluster!" << std::endl;
  }
  else
  {
    found_muon = false;
    std::cout << "MuonFitter tool: Main cluster is NOT coincident with MRD cluster! "
              << "Tank cluster and MRD cluster times are too different." << std::endl;
  }

  // save info about main/max cluster
  h_total_pe_hits->Fill(main_cluster_hits, main_cluster_charge);


  // --------------------------------------------------
  // --- FOUND MUON CANDIDATE -------------------------
  // --------------------------------------------------
  double max_eta = 0.;
  double fitted_track_len = -999.;
  TVector3 fitted_vtx(-999,-999,-999);
  if (found_muon)
  {
    std::cout << "MuonFitter Tool: Found muon candidate!" << std::endl;
    std::cout << "MuonFitter Tool: p" << partnumber << "_";
    if (isData) std::cout << evnum;
    else std::cout << mcevnum;
    std::cout << std::endl;

    // save to txt files
    pehits_file << "p" << partnumber << "_";
    if (isData) pehits_file << evnum;
    else pehits_file << mcevnum;
    pehits_file << "," << main_cluster_hits << "," << main_cluster_charge << std::endl;

    // Load cluster hits & detkeys
    std::vector<Hit> cluster_hits;
    std::vector<MCHit> cluster_hits_MC;
    if (isData) { cluster_hits = m_all_clusters->at(main_cluster_time); }
    else { cluster_hits_MC = m_all_clusters_MC->at(main_cluster_time); }
    std::vector<unsigned long> cluster_detkeys = m_all_clusters_detkeys->at(main_cluster_time);
    std::vector<double> x_hits, y_hits, z_hits;

    // --------------------------------------------------
    // --- Fitted tank track ----------------------------
    // --------------------------------------------------
    //double fitted_track_len = -999.;
    //TVector3 fitted_vtx(-999,-999,-999);
    if (fit_mode)
    {
      std::stringstream ev_id;
      ev_id << "p" << partnumber << "_";
      if (isData) ev_id << evnum;
      else ev_id << mcevnum;

      // Find the fitted track length for this event
      std::map<std::string, double>::iterator it = m_tank_track_fits.find(ev_id.str());
      if (it != m_tank_track_fits.end())
      {
        fitted_track_len = m_tank_track_fits.at(ev_id.str());
        std::cout << " [debug] found track for " << ev_id.str() << ": " << fitted_track_len << endl;
        h_fitted_track_len->Fill(fitted_track_len);

        //calculate initial MIP energy loss in tank from tank track length
        double tank_track = fitted_track_len;
        if (!isData && display_truth)
        { //check energy calc with true track length in water
          tank_track = trueTrackLengthInWater;
        }
        double dEdx = 2.000;    //1.992 MeV/cm: #TODO:get more accurate conversion
        double intank_edep = tank_track*dEdx;
        double outer_track_len = TMath::Sqrt(pow(tankExitPointX-mrdEntryPointX,2) + pow(tankExitPointY-mrdEntryPointY,2) + pow(tankExitPointZ-mrdEntryPointZ,2));
        double outtank_edep = outer_track_len*dEdx;
        std::cout << " [debug] outer_track_len: " << outer_track_len << std::endl;
        std::cout << " [debug] outtank_edep: " << outtank_edep << std::endl;

        // Fine-tune reco muon energy by finding best dE/dx
        double reco_mu_e = intank_edep + energyLoss + outtank_edep;
        double T_before = reco_mu_e;
        bool eres_ok = false;
        int n_tries = 0;
        std::cout << " [debug] initial dEdx: " << dEdx << std::endl;
/*        while (!eres_ok)
        {
          double new_dEdx = Calc_dEdx(T_before);
          std::cout << " [debug] new_dEdx: " << new_dEdx << std::endl;
          double T_after = tank_track*new_dEdx + energyLoss + outer_track_len*new_dEdx;
          double E_res = TMath::Abs(T_after-T_before)/T_before;
          std::cout << " [debug] T_before: " << T_before << std::endl;
          std::cout << " [debug] T_after: " << T_after << std::endl;
          std::cout << " [debug] Eres: " << E_res << std::endl;
          if (E_res <= 0.02)
          {
            reco_mu_e = T_after;  //set as new reco muon energy
            eres_ok = true;
          }
          else
          {
            T_before = T_after;
            n_tries++;
            if (n_tries > 10)
            {
              std::cout << " [debug] Did 10 tries!" << std::endl;
              eres_ok = true;
            }
          }
        } //end while eres_ok loop
*/

        // Fine-tune reco muon energy by finding dEdx for changing energy after traveling some distance dx
        double d_tanktrack = fitted_track_len;
        double dx = 10;   //cm
        double input_Emu = reco_mu_e;
        double sum_tank_Emu = 0.;
        std::cout << " [debug] Calculating muon energy iteratively..." << std::endl;
        std::cout << " [debug] Initial d_tanktrack: " << d_tanktrack << std::endl;
        while (d_tanktrack > 0)
        {
          if (d_tanktrack < dx)
          {
            double deltaE = Calc_dEdx(input_Emu)*d_tanktrack;
            sum_tank_Emu += deltaE;
            d_tanktrack -= d_tanktrack;
            std::cout << " [debug] d_tanktrack: " << d_tanktrack << std::endl;
          }
          else
          {
            double deltaE = Calc_dEdx(input_Emu)*dx;
            sum_tank_Emu += deltaE;
            std::cout << " [debug] dEdx: " << deltaE/dx << std::endl; 
            std::cout << " [debug] sum_tank_Emu: " << sum_tank_Emu << std::endl; 
            
            //calculate new input energy & remaining tank track length
            input_Emu -= deltaE;
            d_tanktrack -= dx;
          }
        }
        std::cout << " [debug] Initial reconstructed Emu: " << reco_mu_e << std::endl;
        reco_mu_e = sum_tank_Emu + energyLoss + outtank_edep;
        std::cout << " [debug] New reconstructed Emu: " << reco_mu_e << std::endl;


        // Use SciBooNE method to determine Edep in MRD
        double mrd_track = reco_mrd_track;    //using reco mrd track
        if (!isData && display_truth)
        {
          mrd_track = trueTrackLengthInMRD;   //using true mrd track
        }
        double mrd_edep_sb = Calc_MrdEnergyLoss(mrd_track);
        std::cout << " [debug] mrd_edep_sb: " << mrd_edep_sb << std::endl;
        h_mrd_eloss_diff->Fill(mrd_edep_sb - energyLoss);       //diff btwn SciBooNE method Edep and ANNIE method Edep
        h_mrd_eloss_SBvANNIE->Fill(energyLoss, mrd_edep_sb);

        // Get the vertex coordinates from fitted track length
        fitted_vtx = tankExit - fitted_track_len*mrdTrackDir;  //TODO:draw best fit and true vtx in 3D geo
        std::cout << " [debug] fitted_vtx (x,y,z): " << fitted_vtx.X() << "," << fitted_vtx.Y() << "," << fitted_vtx.Z() << std::endl;
        h_vtxfit_x->Fill(fitted_vtx.X());
        h_vtxfit_y->Fill(fitted_vtx.Y());
        h_vtxfit_z->Fill(fitted_vtx.Z());
        h_topview_fit->Fill(fitted_vtx.Z(),fitted_vtx.X());
        h_sideview_fit->Fill(fitted_vtx.Z(),fitted_vtx.Y());

        if (!isData)  //truth info
        {
          std::cout << " [debug] true vertex (x,y,z): " << trueVtxX << "," << trueVtxY << "," << trueVtxZ << std::endl;
          h_truefitdiff_x->Fill(fitted_vtx.X()-trueVtxX);
          h_truefitdiff_y->Fill(fitted_vtx.Y()-trueVtxY);
          h_truefitdiff_z->Fill(fitted_vtx.Z()-trueVtxZ);
          h_truefit_len_diff->Fill(TMath::Abs(fitted_track_len-trueTrackLengthInWater));
          h_truevtx_x->Fill(trueVtxX);
          h_truevtx_y->Fill(trueVtxY);
          h_truevtx_z->Fill(trueVtxZ);
          h_topview_truth->Fill(trueVtxZ,trueVtxX);
          h_sideview_truth->Fill(trueVtxZ,trueVtxY);

          // calculcate spatial resolution
          double deltaR = TMath::Sqrt(pow(fitted_vtx.X()-trueVtxX,2) + pow(fitted_vtx.Y()-trueVtxY,2) + pow(fitted_vtx.Z()-trueVtxZ,2));
          std::cout << " [debug] deltaR: " << deltaR << std::endl;
          h_deltaR->Fill(deltaR);
          TVector3 true_vtx = TVector3(trueVtxX,trueVtxY,trueVtxZ);
          double transverse = deltaR*TMath::Sin((true_vtx-fitted_vtx).Angle(mrdTrackDir));
          std::cout << " [debug] transverse: " << transverse << std::endl;
          h_transverse->Fill(TMath::Abs(transverse));
          double parallel = deltaR*TMath::Cos((true_vtx-fitted_vtx).Angle(mrdTrackDir));
          std::cout << " [debug] parallel: " << parallel << std::endl;
          h_parallel->Fill(TMath::Abs(parallel));

          // compare with truth info
          double true_mu_e = trueMuonEnergy - 105.66;
          double ediff = reco_mu_e - true_mu_e;
          h_true_reco_E->Fill(true_mu_e, reco_mu_e);
          h_true_reco_Ediff->Fill(ediff);
          // compare SciBooNE method of MRD energy deposition
          h_true_reco_Ediff_sb->Fill((intank_edep + mrd_edep_sb) - true_mu_e);   //using SB method of determining MRD energy dep
          h_true_reco_Ediff_sb_outerE->Fill((intank_edep + mrd_edep_sb + outtank_edep) - true_mu_e);   //same as above, but with outer track

          // look at fractions of the tank and mrd track lengths
          double f_tanktrack = (tank_track + outer_track_len)/ (tank_track + outer_track_len + reco_mrd_track);
          double f_mrdtrack = (reco_mrd_track) / (tank_track + outer_track_len + reco_mrd_track);
          h_Ediff_frac_tank->Fill(f_tanktrack, ediff);
          h_Ediff_frac_mrd->Fill(f_mrdtrack, ediff);

        }
        h_tank_mrd_E->Fill(energyLoss, intank_edep);
      }
      else std::cout << " no fitted track found for: " << ev_id.str() << std::endl;

/*      std::cout << " [debug] Setting FittedTrackLengthInWater to CStore" << std::endl;
      m_data->CStore.Set("FittedTrackLengthInWater", fitted_track_len);
      std::cout << " [debug] Setting FittedMuonVertex to CStore" << std::endl;
      Position fitted_muon_vtx(fitted_vtx.X(), fitted_vtx.Y(), fitted_vtx.Z());
      m_data->CStore.Set("FittedMuonVertex", fitted_muon_vtx);

      double test_tracklen = -999.;
      m_data->CStore.Get("FittedTrackLengthInWater", test_tracklen);
      std::cout << " [debug] test retrieval of track length from store: " << test_tracklen << std::endl;
      Position test_vtx(-999, -999, -999);
      m_data->CStore.Get("FittedMuonVertex", test_vtx);
      std::cout << " [debug] test retrieval of fitted vertex from store: " << test_vtx.X() << "," << test_vtx.Y() << "," << test_vtx.Z() << std::endl;
*/
    }

    // --------------------------------------------------
    // METHOD: CALCULATE ai FOR EACH PMT
    // --------------------------------------------------
    std::map<unsigned long, double> charge;
    std::map<unsigned long, std::vector<double>> hittime;
    std::map<int, std::vector<double>> m_PE_ai;
    std::map<int, std::vector<double>> m_fpmt_ai;

    charge.clear();
    hittime.clear();
    for (unsigned long detkey = 332; detkey < 464; ++detkey)
    {
      charge.emplace(detkey, 0.);
      hittime.emplace(detkey, std::vector<double>{-99.});
    }

    if (isData)
    {
      for (int i = 0; i < (int)cluster_hits.size(); ++i)
      {
        //for each hit
        int chankey = cluster_hits.at(i).GetTubeId();
        std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(chankey);
        if (it != ChannelKeyToSPEMap.end())
        {
          Detector* this_detector = geom->ChannelToDetector(chankey);
          unsigned long detkey = this_detector->GetDetectorID();  //chankey same as detkey
          Position det_pos = this_detector->GetDetectorPosition();
          double hit_charge = cluster_hits.at(i).GetCharge();
          double hit_PE = hit_charge / ChannelKeyToSPEMap.at(chankey);
          double hit_time = cluster_hits.at(i).GetTime();
          hittime[detkey].push_back(hit_time);
          //h_charge_detkey->Fill(detkey, hit_PE);

          // keep track of total charge seen by each PMT
          charge[detkey] += hit_PE;
        } //end if ChannelKeyToSPEMap
      } //end cluster_hits loop
    }
    else  //is MC
    {
      for (int i = 0; i < (int)cluster_hits_MC.size(); ++i)
      {
        //for each hit
        int chankey = cluster_hits_MC.at(i).GetTubeId();
        std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(chankey);
        if (it != ChannelKeyToSPEMap.end())
        {
          Detector* this_detector = geom->ChannelToDetector(chankey);
          unsigned long detkey = this_detector->GetDetectorID();  //chankey same as detkey
          Position det_pos = this_detector->GetDetectorPosition();
          double hit_PE = cluster_hits_MC.at(i).GetCharge();  //in PE
          double hit_charge = hit_PE * ChannelKeyToSPEMap.at(chankey);
          double hit_time = cluster_hits_MC.at(i).GetTime();
          hittime[detkey].push_back(hit_time);
          h_charge_detkey->Fill(detkey, hit_PE);

          // keep track of total charge seen by each PMT
          charge[detkey] += hit_PE;
        } //end if ChannelKetToSPEMap
      } //end cluster_hits_MC loop
    }

    // Go through charge, hittime maps
    int gi = 0;    //counter for TGraph
    int nhits = 0, nhits_incone = 0;
    double totalpe = 0, totalpe_incone = 0;
    std::vector<double> v_tzero;
    for (unsigned long detkey = 332; detkey < 464; ++detkey)
    {
      double pmt_PE = charge[detkey];
      h_pmt_charge->Fill(pmt_PE);
      
      if (pmt_PE == 0) continue;    //NOTE: there will be PMTs w/ 0 p.e. from charge map initialization

      std::vector<double> v_hit_times = hittime[detkey];
      double pmt_t = 0;
      int nhits_pmt = v_hit_times.size();
      std::cout << " [debug] all hit times for detkey " << detkey << ": ";
      for (int ht = 0; ht < nhits_pmt; ++ht) { std::cout << v_hit_times.at(ht) << ","; }
      std::cout << std::endl;

      //get earliest hit time bc Cherenkov radiation is usually earliest time
      sort(v_hit_times.begin(), v_hit_times.end());
      pmt_t = v_hit_times.at(0);
      if (pmt_t == -99.) pmt_t = v_hit_times.at(1);
      std::cout << " [debug] final pmt_t: " << pmt_t << std::endl;
      

      //PE CUT
      if (pmt_PE < PMTQCut)
      {
        std::cout << " [debug] SKIPPING charge[" << detkey << "] < " << PMTQCut << "pe: " << charge[detkey] << std::endl;
        pmt_PE = 0.;    //set PMT charge to 0
        continue;       //skip PMT entirely
      }

      nhits += 1;   //keep track of total hits that meet PE cut in this event
      totalpe += pmt_PE;

      double hitX = x_pmt[detkey];
      double hitY = y_pmt[detkey];
      double hitZ = z_pmt[detkey];
      double dirPMTX = x_pmt_dir[detkey];
      double dirPMTY = y_pmt_dir[detkey];
      double dirPMTZ = z_pmt_dir[detkey];
      
      TVector3 pmt_dir = TVector3(dirPMTX, dirPMTY, dirPMTZ).Unit();

      // get vector from tankExitPoint to PMT (Ri)
      TVector3 vec_Ri = TVector3(hitX,hitY,hitZ) - TVector3(tankExitPointX,tankExitPointY,tankExitPointZ);
      double Ri = vec_Ri.Mag();
      h_tankexit_to_pmt->Fill(Ri);
      std::cout << " [debug] vec_Ri: " << vec_Ri.X() << "," << vec_Ri.Y() << "," << vec_Ri.Z() << std::endl;
      std::cout << " [debug] vec_Ri direction: " << vec_Ri.Unit().X() << "," << vec_Ri.Unit().Y() << "," << vec_Ri.Unit().Z() << std::endl;

      // get angle btwn Ri and muon direction (ai)
      double ang_alpha = vec_Ri.Angle(-mrdTrackDir);  //rad
      std::cout << " [debug] tankExitPoint: " << tankExitPointX << "," << tankExitPointY << "," << tankExitPointZ << std::endl;

      double ai = Ri * TMath::Sin(ang_alpha) / TMath::Tan(CHER_ANGLE_RAD) + Ri * TMath::Cos(ang_alpha);
      h_tanktrack_ai->Fill(ai);

      // get the vector from the vertex of ai to the PMT (bi)
      //TVector3 vec_ai = TVector3(tankExitPointX, tankExitPointY, tankExitPointZ) - ai*mrdTrackDir;  //in tank along mrd track (opposite mrdTrackDir)
      TVector3 vec_ai = tankExit - ai*mrdTrackDir;    //vector to where ai starts in tank, along mrd track
      TVector3 vec_bi = TVector3(hitX,hitY,hitZ) - vec_ai;
      double bi = TMath::Sqrt(pow(hitX-vec_ai.X(),2) + pow(hitY-vec_ai.Y(),2) + pow(hitZ-vec_ai.Z(),2));
      std::cout << " [debug] vec_ai: " << vec_ai.X() << "," << vec_ai.Y() << "," << vec_ai.Z() << std::endl;
      std::cout << " [debug] vec_bi: " << vec_bi.X() << "," << vec_bi.Y() << "," << vec_bi.Z() << std::endl;

      // do something with time
      double t0 = -999.;
      if (fit_mode)
      {
        // check if hit falls inside cone of FITTED vertex
        TVector3 vtx2tankExit = tankExit - fitted_vtx;   //similar to ai vector
        TVector3 vtx2pmt = TVector3(hitX,hitY,hitZ) - fitted_vtx;

        double ang = vtx2pmt.Angle(mrdTrackDir);

        if (!isData && display_truth)
        {
          std::cout << " [debug] vtx2tankExit (before truth): " << vtx2tankExit.X() << "," << vtx2tankExit.Y() << "," << vtx2tankExit.Z() << std::endl;
          vtx2tankExit = TVector3(trueStopVtxX,trueStopVtxY,trueStopVtxZ) - TVector3(trueVtxX,trueVtxY,trueVtxZ);
          vtx2pmt = TVector3(hitX,hitY,hitZ) - TVector3(trueVtxX,trueVtxY,trueVtxZ);
          std::cout << " [debug] vtx2tankExit (after truth): " << vtx2tankExit.X() << "," << vtx2tankExit.Y() << "," << vtx2tankExit.Z() << std::endl;
          //ang = vtx2pmt.Angle(vtx2tankExit);
          ang = vtx2pmt.Angle(trueTrackDir);
          std::cout << " [debug] ang (w/ trueTrackDir, w/ vtx2tankExit): " << ang << "," << vtx2pmt.Angle(vtx2tankExit) << std::endl;
        }

        if (ang <= CHER_ANGLE_RAD)
        {
          //t0 = pmt_t - (TVector3(hitX,hitY,hitZ).Dot(mrdTrackDir) + ai) - TMath::Sqrt(1.33*1.33-1)*TVector3(hitX,hitY,hitZ).Cross(mrdTrackDir).Mag();
          t0 = pmt_t - (vtx2pmt.Dot(mrdTrackDir) + TMath::Sqrt(1.33*1.33-1)*vtx2pmt.Cross(mrdTrackDir).Mag()) / 30.;

          if (!isData && display_truth)
          {
            t0 = pmt_t - (vtx2pmt.Dot(trueTrackDir) + TMath::Sqrt(1.33*1.33-1)*vtx2pmt.Cross(trueTrackDir).Mag()) / 30.;
          }
          std::cout << " [debug] pmt_t made it inside cone! " << pmt_t << "," << t0 << std::endl;
          v_tzero.push_back(t0);
          nhits_incone += 1;
          totalpe_incone += pmt_PE;
        }
        else { std::cout << " [debug] pmt_t did not make it inside cone! " << pmt_t << std::endl; }
      }


      // find angle btwn true vertex and hit
      if (!isData)
      {
        TVector3 vec_truevtx_hit = TVector3(hitX,hitY,hitZ) - TVector3(trueVtxX,trueVtxY,trueVtxZ);
        double anglePmtTrueVtx = vec_truevtx_hit.Angle(mrdTrackDir)*180./TMath::Pi();
      }

      // get angle btwn vector bi and pmt dir
      double psi = vec_bi.Angle(-pmt_dir);    //used only to get angle; don't use for magnitude
      if (psi > TMath::Pi()/2.) psi = TMath::Pi()-psi;
      //std::cout << " [debug] psi(+), psi(-): " << psi*180./TMath::Pi() << ", " << vec_bi.Angle(pmt_dir)*180./TMath::Pi() << std::endl;

      // calculate the area of PMT & frustum
      double eff_area_pmt = 0.5*m_pmt_area[detkey]*(1.+TMath::Cos(psi));
      h_eff_area_pmt->Fill(eff_area_pmt);
      double area_frustum = 2.*TMath::Pi()*(deltaL)*Ri;
      double f_pmt = eff_area_pmt / area_frustum;
      h_fpmt->Fill(f_pmt);

      gr_effarea_detkey->SetPoint(gi, detkey, eff_area_pmt);
      gr_fpmt_detkey->SetPoint(gi, detkey, f_pmt);
      gr_effarea_ai->SetPoint(gi, ai, eff_area_pmt);
      gr_fpmt_ai->SetPoint(gi, ai, f_pmt);

      h_eta_ai->Fill(ai, pmt_PE / f_pmt);

      for (int a = 55; a < 500; a+=deltaL)
      {
        if (ai >= a-deltaL/2 && ai < a+deltaL/2)
        {
          m_PE_ai[a].push_back(pmt_PE);
          m_fpmt_ai[a].push_back(f_pmt);
        }
      }


      //XXX: checking PMT direction by drawing sphere in front of PMT
/*      TVector3 Ri_dir = TVector3(hitX, hitY, hitZ) + 10.*vec_Ri.Unit();
      if (detkey == 400)
      {
        // tankExit
        sprintf(blockName,"tankExit_%u",N);
        bBlock = ageom->MakeSphere(blockName, Iron, 0,2,0,180,0,360);
        bBlock->SetLineColor(3);
        EXPH->AddNodeOverlap(bBlock,69,new TGeoTranslation(tankExitPointX, tankExitPointY, tankExitPointZ));
        N++;

        // direction of Ri
        sprintf(blockName,"RiDir_%u",N);
        bBlock = ageom->MakeSphere(blockName, Iron, 0,1.5,0,180,0,360);
        bBlock->SetLineColor(5);
        EXPH->AddNodeOverlap(bBlock,69,new TGeoTranslation(Ri_dir.X(), Ri_dir.Y(), Ri_dir.Z()));
        N++;

        // bi vertex
        sprintf(blockName, "bi_%u",N);
        bBlock = ageom->MakeSphere(blockName, Iron, 0,2,0,180,0,360);
        bBlock->SetLineColor(9);
        EXPH->AddNodeOverlap(bBlock,69,new TGeoTranslation(vec_bi.X()+vec_ai.X(), vec_bi.Y()+vec_ai.Y(), vec_bi.Z()+vec_ai.Z()));
        N++;
      } //end if detkey==400 */

      ++gi;

    } //end charge[detkey] loop

    if (fit_mode)
    {
      nhits_trlen_file << "p" << partnumber << "_";
      if (isData) nhits_trlen_file << evnum;
      else nhits_trlen_file << mcevnum;
      nhits_trlen_file << "," << fitted_track_len << "," << nhits << "," << nhits_incone << "," << totalpe << "," << totalpe_incone << std::endl;
    }


    // Find mean t0 and subtract to center hist around 0
    double mean_t0 = 0.;
    for (int t = 0; t < v_tzero.size(); ++t) { mean_t0 += v_tzero.at(t); }
    if (v_tzero.size() != 0) mean_t0 /= v_tzero.size();
    std::cout << " [debug] mean_t0: " << mean_t0 << std::endl;

    for (int t = 0; t < v_tzero.size(); ++t) { h_tzero->Fill(v_tzero.at(t)-mean_t0); }

    h_uber_t0widths->Fill(h_tzero->GetStdDev());


    // find which pmt charges are inside/outside cone of true track direction
    if (!isData)
    {
      for (unsigned long detkey = 332; detkey < 464; ++detkey)
      {
        if (charge[detkey] != 0.)
        {
          // find vector from true vertex to PMT
          TVector3 vec_pmt = TVector3(x_pmt[detkey],y_pmt[detkey],z_pmt[detkey]) - TVector3(trueVtxX,trueVtxY,trueVtxZ);
          std::cout << " [debug] xyz_pmt[detkey] " << detkey << ": " << x_pmt[detkey] << "," << y_pmt[detkey] << "," << z_pmt[detkey] << std::endl;
          // find angle between pmt vector (vec_pmt) and true track direction
          double pmt_angle = vec_pmt.Angle(trueTrackDir)*180./TMath::Pi();
          std::cout << " [debug] vec_pmt: " << vec_pmt.X() << "," << vec_pmt.Y() << "," << vec_pmt.Z() << std::endl;
          std::cout << " [debug] pmt_angle: " << pmt_angle << std::endl;
          std::cout << " [debug] trueTrackDir: " << trueTrackDir.X() << "," << trueTrackDir.Y() << "," << trueTrackDir.Z() << std::endl;
          std::cout << " [debug] trueVtx: " << trueVtxX << "," << trueVtxY << "," << trueVtxZ << std::endl;
          std::cout << " [debug] trueStopVtx: " << trueStopVtxX << "," << trueStopVtxY << "," << trueStopVtxZ << std::endl;

          if (pmt_angle < CHER_ANGLE_DEG+insideAngle) h_qincone_truevtx->Fill(charge[detkey]);
          if (pmt_angle > CHER_ANGLE_DEG+outsideAngle) h_qoutcone_truevtx->Fill(charge[detkey]);
        }
      }
    }

    // Make graph of eta vs ai (tank track length)
    int j = 0;
    double running_avg = 0;
    bool found_vtx = false;
    for (auto const &pair: m_fpmt_ai)
    {
      double total_PE_ai = 0;
      double total_fpmt_ai = 0;
      for (int e = 0; e < pair.second.size(); ++e)
      {
        total_fpmt_ai += pair.second.at(e);
        total_PE_ai += m_PE_ai[pair.first].at(e);
      }
      double total_eta = total_PE_ai / total_fpmt_ai;
      if (total_eta > max_eta) max_eta = total_eta;     //for plotting
      
      gr_eta_ai->SetPoint(j, pair.first, total_eta);
      // save (ai,eta) data to txt file
      pos_file << "p" << partnumber << "_";
      if (isData) pos_file << evnum;
      else pos_file << mcevnum;
      pos_file << "," << pair.first << "," << total_eta << std::endl;
      avg_eta += total_eta;
      
      if (j == 0) running_avg = avg_eta;
      else
      {
        std::cout << " [debug] avg_eta (before dividing j, so it is total eta): " << avg_eta << std::endl;
        running_avg = avg_eta / (j+1);

        if (running_avg < EtaThreshold && !found_vtx)
        {
          // first instance the avg dips below threshold
          std::cout << " [debug] BEST VERTEX (ai): " << pair.first << std::endl;
          found_vtx = true;
          bestFitAi = pair.first;
        }
      }

      std::cout << " [debug] ai, running_avg: " << pair.first << ", " << running_avg << std::endl;
      gr_running_avg->SetPoint(j, pair.first, running_avg);

      if (!isData)
      {
        if (pair.first < trueTrackLengthInWater)   //left of true track length
        {
          std::cout << " [debug] total_eta (left): " << total_eta << std::endl;
          num_left_eta++;
          left_avg_eta += total_eta;
        }
        if (pair.first > trueTrackLengthInWater)   //right of true track length
        {
          std::cout << " [debug] total_eta (right): " << total_eta << std::endl;
          num_right_eta++;
          right_avg_eta += total_eta;
        }
      }
      j+=1;
    }
    avg_eta /= j;
    h_avg_eta->Fill(avg_eta);
    std::cout << " [debug] avg eta: " << avg_eta << ", j: " << j << std::endl;

    if (!isData)
    {
      std::cout << " [debug] num pmts (left, right): " << num_left_eta << ", " << num_right_eta << std::endl;
      if (num_left_eta != 0) left_avg_eta /= num_left_eta;
      if (num_right_eta != 0) right_avg_eta /= num_right_eta;
      std::cout << " [debug] left avg: " << left_avg_eta << ", right avg: " << right_avg_eta << std::endl;
      h_lr_avg_eta->Fill(left_avg_eta);
      h_lr_avg_eta->Fill(right_avg_eta);
    }

    // running avg of entire graph
    // this might just be smoothing...
    for (int e = 1; e < j-2; ++e)
    {
      double three_pt_avg = (gr_eta_ai->GetY()[e-1] + gr_eta_ai->GetY()[e] + gr_eta_ai->GetY()[e+1]) / 3.;
      //gr_running_avg->SetPoint(e-1, gr_eta_ai->GetX(e), three_pt_avg);
    }

    // save to txt files
    if (!isData && display_truth)
    {
      // Save truth information
      truetrack_file << "p" << partnumber << "_" << mcevnum << ",";
      truetrack_file << trueTrackLengthInWater << "," << left_avg_eta << "," << right_avg_eta << std::endl;
    }

  
    // --------------------------------------------------
    // METHOD: STEPPING BACK FROM TANKEXITPT FOR VTX CANDIDATES
    // --------------------------------------------------
    std::map<int, double> m_vtx_charge_incone;   //map of total charge seen at ea vtx candidate
    std::map<int, double> m_vtx_charge_per_pmt_incone;
    std::map<int, double> m_vtx_detkey_incone;
    std::map<int, std::vector<double>> m_vtx_vec_charge_incone;   //map of charge vec seen at ea vtx candidate
    std::map<int, std::vector<unsigned long>> m_vtx_vec_detkey_incone;
    std::map<int, std::vector<double>> m_vtx_vec_charge_outcone;
    std::map<int, std::vector<unsigned long>> m_vtx_vec_detkey_outcone;
    for (int c = 0; c < vtxCandidates.size(); ++c)
    {
      //for each vtx candidate
      TVector3 vtx_to_mrd = TVector3(mrdStart.X(), mrdStart.Y(), mrdStart.Z()) - vtxCandidates.at(c);

      // Find all pmts that fall in cone, hit or not
      std::vector<unsigned long> v_all_detkeys_incone;
      for (unsigned long detkey = 332; detkey < 465; ++detkey)
      {
        TVector3 vtx_to_pmt = TVector3(x_pmt[detkey], y_pmt[detkey], z_pmt[detkey]) - vtxCandidates.at(c);
        double pmtAngle = vtx_to_pmt.Angle(vtx_to_mrd)*180./TMath::Pi();
        if (pmtAngle < CHER_ANGLE_DEG)
        {
          v_all_detkeys_incone.push_back(detkey);
        }
      }
      //std::cout << " [debug] total detkeys in cone: " << v_all_detkeys_incone.size() << std::endl;
      gr_vtx_detkey_in->SetPoint(c, vtxCandidates.at(c).Z(), v_all_detkeys_incone.size());
      gr_vtx_detkey_out->SetPoint(c, vtxCandidates.at(c).Z(), 132-v_all_detkeys_incone.size());

      std::vector<unsigned long> v_detkeys_incone;
      std::vector<unsigned long> v_detkeys_outcone;
      std::vector<double> v_charges_incone;
      std::vector<double> v_charges_outcone;
      std::cout << " [debug] WORKING ON VTX c = " << c << std::endl;
      if (isData)
      {
        for (int i = 0; i < (int)cluster_hits.size(); ++i)
        {
          //for each hit
          int chankey = cluster_hits.at(i).GetTubeId();
          std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(chankey);
          if (it != ChannelKeyToSPEMap.end())
          {
            Detector* this_detector = geom->ChannelToDetector(chankey);
            unsigned long detkey = this_detector->GetDetectorID();  //chankey same as detkey
            Position det_pos = this_detector->GetDetectorPosition();
            double hit_charge = cluster_hits.at(i).GetCharge();
            double hit_PE = hit_charge / ChannelKeyToSPEMap.at(chankey);
            double hit_time = cluster_hits.at(i).GetTime();
            double hitX = 100.*(det_pos.X()-tank_center_x);   //[cm]
            double hitY = 100.*(det_pos.Y()-tank_center_y);
            double hitZ = 100.*(det_pos.Z()-tank_center_z);

            Direction dirPMT = this_detector->GetDetectorDirection();

            //angle of vector from vtx to PMT wrt to track
            TVector3 vtx_to_hit = TVector3(hitX, hitY, hitZ) - vtxCandidates.at(c);
            double pmtAngleRad = vtx_to_hit.Angle(vtx_to_mrd);
            double pmtAngleDeg = pmtAngleRad*180./TMath::Pi();
            h_hit_angles->Fill(pmtAngleDeg);

            //TVector3 pmt_dir = TVector3(hitX, hitY, hitZ) - TVector3(dirPMT.X()*100., dirPMT.Y()*100., dirPMT.Z()*100.);
            TVector3 pmt_dir = TVector3(dirPMT.X()*100., dirPMT.Y()*100., dirPMT.Z()*100.).Unit();

            //XXX: checking PMT direction by drawing sphere in front PMT
            TVector3 p2 = TVector3(hitX, hitY, hitZ) - 10.*pmt_dir;
            if (detkey == 400)
            {
              sprintf(blockName,"pmtDir_%u",N);
              bBlock = ageom->MakeSphere(blockName, Iron, 0,5,0,180,0,360);
              bBlock->SetLineColor(6);
              EXPH->AddNodeOverlap(bBlock,69,new TGeoTranslation(p2.X(), p2.Y(), p2.Z()));
              N++;
            }

            if (pmtAngleDeg > CHER_ANGLE_DEG)
            {
              v_detkeys_outcone.push_back(detkey);
              //v_charges_outcone.push_back(hit_charge);
              v_charges_outcone.push_back(hit_PE);
            }
            else
            {
              v_detkeys_incone.push_back(detkey);
              //v_charges_incone.push_back(hit_charge);
              v_charges_incone.push_back(hit_PE);
              std::cout << " [debug] Ev " << evnum << ": detkey " << detkey << " charge made it in cone (nC, PE): " << hit_charge << ", " << hit_PE << std::endl;

              // Find the incident angle of photon wrt pmt
              double photInAngleRad = pmt_dir.Angle(vtx_to_hit);
              double photInAngleDeg = photInAngleRad*180./TMath::Pi();
              h_phot_inc_angle->Fill(photInAngleDeg);
              //std::cout << "  [debug] incident photon angle (rad, deg): " << photInAngleRad << ", " << photInAngleDeg << std::endl;
            }
          } //exit ChannelKeyToSPEMap if statement
        } //done looping through cluster hits for this vtx
      }
      else  //is MC
      {
        for (int i = 0; i < (int)cluster_hits_MC.size(); ++i)
        {
          //for each hit
          int chankey = cluster_hits_MC.at(i).GetTubeId();
          std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(chankey);
          if (it != ChannelKeyToSPEMap.end())
          {
            Detector* this_detector = geom->ChannelToDetector(chankey);
            unsigned long detkey = this_detector->GetDetectorID();  //chankey same as detkey
            Position det_pos = this_detector->GetDetectorPosition();
            double hit_PE = cluster_hits_MC.at(i).GetCharge();  //in PE
            double hit_charge = hit_PE * ChannelKeyToSPEMap.at(chankey);
            double hit_time = cluster_hits_MC.at(i).GetTime();
            double hitX = 100.*(det_pos.X()-tank_center_x);   //[cm]
            double hitY = 100.*(det_pos.Y()-tank_center_y);
            double hitZ = 100.*(det_pos.Z()-tank_center_z);

            // individual PMT charge cut
/*            if (hit_PE < 10) continue;
            std::cout << " [debug] more than 10 PE!" << std::endl;
            if (hit_time > earliest_hittime+40) continue;
            std::cout << " [debug] hit is within 40 ns of earliest hit time!" << std::endl;
*/

            Direction dirPMT = this_detector->GetDetectorDirection();

            //angle of vector from vtx to PMT wrt to track
            TVector3 vtx_to_hit = TVector3(hitX, hitY, hitZ) - vtxCandidates.at(c);
            double pmtAngleRad = vtx_to_hit.Angle(vtx_to_mrd);
            double pmtAngleDeg = pmtAngleRad*180./TMath::Pi();
            h_hit_angles->Fill(pmtAngleDeg);

            TVector3 pmt_dir = TVector3(dirPMT.X()*100., dirPMT.Y()*100., dirPMT.Z()*100.).Unit();

            //XXX: checking PMT direction by drawing sphere in front PMT
/*            TVector3 p2 = TVector3(hitX, hitY, hitZ) + 10.*pmt_dir;
            if (detkey == 400)
            {
              sprintf(blockName,"pmtDir_%u",N);
              bBlock = ageom->MakeSphere(blockName, Iron, 0,1.5,0,180,0,360);
              bBlock->SetLineColor(6);
              EXPH->AddNodeOverlap(bBlock,69,new TGeoTranslation(p2.X(), p2.Y(), p2.Z()));
              N++;
            }*/

            if (pmtAngleDeg > CHER_ANGLE_DEG)
            {
              v_detkeys_outcone.push_back(detkey);
              v_charges_outcone.push_back(hit_PE);
            }
            else
            {
              v_detkeys_incone.push_back(detkey);
              v_charges_incone.push_back(hit_PE);
              std::cout << " [debug] Ev " << evnum << ": detkey " << detkey << " charge made it in cone (nC, PE): " << hit_charge << ", " << hit_PE << std::endl;

              // Find the incident angle of photon wrt pmt
              double photInAngleRad = pmt_dir.Angle(vtx_to_hit);
              double photInAngleDeg = photInAngleRad*180./TMath::Pi();
              h_phot_inc_angle->Fill(photInAngleDeg);
            }


            // --------------------------------------------------
            // --- PMT and cone/frustum area correction ---------
            // --------------------------------------------------
            if (c==0)   //only want to do this with cluster once
            {
            } // end if c==0

          } //exit ChannelKeyToSPEMap if statement 
        } //done looping through cluster hits for this vtx c
      } //end isData/isMC

      // INSIDE cone at each vtx
      double tot_charge_incone = 0;
      for (int q = 0; q < v_charges_incone.size(); ++q) { tot_charge_incone += v_charges_incone.at(q); }
      gr_vtx_charge_in->SetPoint(c, vtxCandidates.at(c).Z(), tot_charge_incone);
      m_vtx_charge_incone[c] = tot_charge_incone;   //same as graph

      m_vtx_vec_charge_incone[c] = v_charges_incone;
      m_vtx_vec_detkey_incone[c] = v_detkeys_incone;
      m_vtx_vec_charge_outcone[c] = v_charges_outcone;
      m_vtx_vec_detkey_outcone[c] = v_detkeys_outcone;

      //calculcate charge per pmt
      //cpp_file << evnum << "," << vtxCandidates.size() <<",";
      int tot_pmts_incone = v_all_detkeys_incone.size();
      if (tot_pmts_incone == 0)
      {
        gr_qdensity_in->SetPoint(c, vtxCandidates.at(c).Z(), 0.);
        avg_charge_per_pmt_in += 0.;
        m_vtx_charge_per_pmt_incone[c] = 0.;
        //std::cout << " ev" << evnum << " charge per pmt (in): 0." << std::endl;
        //cpp_file << "0.,";
      }
      else
      {
        gr_qdensity_in->SetPoint(c, vtxCandidates.at(c).Z(), tot_charge_incone / tot_pmts_incone);
        avg_charge_per_pmt_in += tot_charge_incone / tot_pmts_incone;
        m_vtx_charge_per_pmt_incone[c] = tot_charge_incone / tot_pmts_incone;
        //std::cout << " ev" << evnum << " charge per pmt (in): " << tot_charge_incone/tot_pmts_incone << std::endl;
        //cpp_file << tot_charge_incone/tot_pmts_incone << ",";
      }

      //debug
      //std::cout << "  [debug] detkeys that were hit (incone): ";
      //for (int d = 0; d < v_detkeys_incone.size(); ++d) { std::cout << v_detkeys_incone.at(d) << ", "; }
      //std::cout << std::endl;

      //consolidate hit pmts (find duplicates)
      sort(v_detkeys_incone.begin(), v_detkeys_incone.end());
      std::map<unsigned long, int> m_hit_pmts_incone;
      for (int d = 0; d < v_detkeys_incone.size(); ++d)
      {
        m_hit_pmts_incone[v_detkeys_incone.at(d)]++;
      }
      //std::cout << "  [debug] checking duplicate detkeys...total: " << m_hit_pmts_incone.size() << std::endl;
      //for (auto const &pair: m_it_pmts_incone) { std::cout << "{" << pair.first << ":" << pair.second <<"}" << std::endl; }


      // OUTSIDE cone at each vtx
      double tot_charge_outcone = 0;
      for (int q = 0; q < v_charges_outcone.size(); ++q) { tot_charge_outcone += v_charges_outcone.at(q); }
      gr_vtx_charge_out->SetPoint(c, vtxCandidates.at(c).Z(), tot_charge_outcone);

/*        std::cout << "  [debug] detkeys (outcone): ";
      for (int d = 0; d < v_detkeys_outcone.size(); ++d) { std::cout << v_detkeys_outcone.at(d) << ", "; }
      std::cout << std::endl;*/

      int tot_pmts_outcone = 132-v_all_detkeys_incone.size();
      if (tot_pmts_outcone == 0)
      {
        gr_qdensity_out->SetPoint(c, vtxCandidates.at(c).Z(), 0.);
        avg_charge_per_pmt_out += 0.;
        //std::cout << " ev" << evnum << " charge per pmt (out): 0." << std::endl;
        //cpp_file << "0.,";
      }
      else
      {
        gr_qdensity_out->SetPoint(c, vtxCandidates.at(c).Z(), tot_charge_outcone / tot_pmts_outcone);
        avg_charge_per_pmt_out += tot_charge_outcone / tot_pmts_outcone;
        //std::cout << " ev" << evnum << " charge per pmt (out): " << tot_charge_outcone / tot_pmts_outcone << std::endl;
        //cpp_file << tot_charge_outcone/tot_pmts_outcone << ",";
      }

      //cpp_file << avg_charge_per_pmt_in << "," << avg_charge_per_pmt_out << ",," << std::endl;

      //ANNIE in 3D: scale size of PMT with amount of charge seen in ENTIRE TANK
      // XXX:scale size of PMT w amt of charge seen in entire tank or in/out of cone?
      // TODO:this should also be done with just the selected vtx, not the final
      //IN CONE
/*        std::map<unsigned long, double> m_pmt_charge_incone;
      for (int p = 0; p < v_charges_incone.size(); ++p)
      {
        unsigned long detkey = v_detkeys_incone.at(p);
        m_pmt_charge_incone[detkey] += v_charges_incone.at(p);
        //std::cout << "  [debug] m_pmt_charge_incone[" << detkey << "]: " << m_pmt_charge_incone[detkey] << std::endl;
      }
      for (auto const &pair: m_pmt_charge_incone)
      {
        TGeoNode *node = EXPH->GetNode(detkey_to_node[pair.first]);
        //if (pair.first == 400) node ->GetVolume()->SetLineColor(6);
        node->GetVolume()->SetLineColor(4);
        TGeoSphere *sphere = (TGeoSphere *)(node->GetVolume()->GetShape());
        sphere->SetSphDimensions(0,20.*(pair.second/(tot_charge_incone+tot_charge_outcone))+5.,0,180,0,360);
      }

      //OUT CONE
        node->GetVolume()->SetLineColor(30);
        TGeoSphere *sphere = (TGeoSphere *)(node->GetVolume()->GetShape());
        sphere->SetSphDimensions(0,20.*(pair.second/(tot_charge_outcone+tot_charge_incone))+5.,0,180,0,360);
      } */
    } //done going thru vtxCandidates

      //compare avg charge per pmt
      avg_charge_per_pmt_in = avg_charge_per_pmt_in / vtxCandidates.size();
      avg_charge_per_pmt_out = avg_charge_per_pmt_out / vtxCandidates.size();

      // --------------------------------------------------
      // --- Find best vertex -----------------------------
      // --------------------------------------------------
      max_vtx_charge = m_vtx_charge_per_pmt_incone[0];
      for (auto const &pair: m_vtx_charge_per_pmt_incone)
      {
        if (pair.second > max_vtx_charge)
        {
          max_vtx_charge = pair.second;
          max_vertex = pair.first;
        }
      }
      auto X = std::max_element(m_vtx_charge_per_pmt_incone.begin(), m_vtx_charge_per_pmt_incone.end(), [] (const pair<int, double>& p1, const pair<int, double>&p2) { return p1.second < p2.second; });
      max_vertex = X->first;
      max_vtx_charge = X->second;
      std::cout << " Max charge using max_element method: " << max_vertex << ", " << max_vtx_charge << std::endl;

      //TODO: find bestVtx
      TVector3 bestVtx(vtxCandidates.at(max_vertex).X(), vtxCandidates.at(max_vertex).Y(), vtxCandidates.at(max_vertex).Z());
      double tankTrackLength = TMath::Sqrt(pow(mrdStart.X()-bestVtx.X(),2) + pow(mrdStart.Y()-bestVtx.Y(),2) + pow(mrdStart.Z()-bestVtx.Z(),2));
      //std::cout << " [debug] tankTrackLength: " << tankTrackLength << std::endl;
      h_fitted_tank_track_len->Fill(tankTrackLength);
      double distClosestApproach = TMath::Sqrt(pow(bestVtx.X()-tank_center_x,2) + pow(bestVtx.Y()-tank_center_y,2) + pow(bestVtx.Z()-tank_center_z,2));
      h_closest_approach->Fill(distClosestApproach);

      double best_vtx_tot_charge_in = 0., best_vtx_tot_charge_out = 0.;
      std::cout << " [debug] best vtx charge vector (in): ";
      for (int i = 0; i < m_vtx_vec_charge_incone[max_vertex].size(); ++i)
      {
        best_vtx_tot_charge_in += m_vtx_vec_charge_incone[max_vertex].at(i);
        std::cout << m_vtx_vec_charge_incone[max_vertex].at(i) << ",";
      }
      std::cout << " total: " << best_vtx_tot_charge_in << std::endl;

      std::cout << " [debug] best vtx charge vector (out): ";
      for (int i = 0; i < m_vtx_vec_charge_outcone[max_vertex].size(); ++i)
      {
        best_vtx_tot_charge_out += m_vtx_vec_charge_outcone[max_vertex].at(i);
        std::cout << m_vtx_vec_charge_outcone[max_vertex].at(i) << ",";
      }
      std::cout << " total: " << best_vtx_tot_charge_out << std::endl;

      //ANNIE in 3D
/*      sprintf(blockName,"maxTankVtx_%d", c);
      bBlock = ageom->MakeSphere(blockName, Iron, 0,1.5,0,180,0,360);
      bBlock->SetLineColor(6);
      EXPH->AddNodeOverlap(bBlock,69,new TGeoTranslation(vtxCandidates.at(max_vertex).X(), vtxCandidates.at(max_vertex).Y(), vtxCandidates.at(max_vertex).Z()));
      N++; */

      std::map<unsigned long, double> m_pmt_charge_incone;
      for (int p = 0; p < m_vtx_vec_detkey_incone[max_vertex].size(); ++p)
      {
        unsigned long detkey = m_vtx_vec_detkey_incone[max_vertex].at(p);
        m_pmt_charge_incone[detkey] += m_vtx_vec_charge_incone[max_vertex].at(p);
      }
      for (auto const &pair: m_pmt_charge_incone)
      {
        TGeoNode *node = EXPH->GetNode(detkey_to_node[pair.first]);
        node->GetVolume()->SetLineColor(4);
        TGeoSphere *sphere = (TGeoSphere *)(node->GetVolume()->GetShape());
        sphere->SetSphDimensions(0,20.*(pair.second/(best_vtx_tot_charge_in+best_vtx_tot_charge_out))+5.,0,180,0,360);
      }

      std::map<unsigned long, double> m_pmt_charge_outcone;
      for (int p = 0; p < m_vtx_vec_detkey_outcone[max_vertex].size(); ++p)
      {
        unsigned long detkey = m_vtx_vec_detkey_outcone[max_vertex].at(p);
        m_pmt_charge_outcone[detkey] += m_vtx_vec_charge_outcone[max_vertex].at(p);
      }
      for (auto const &pair: m_pmt_charge_outcone)
      {
        TGeoNode *node = EXPH->GetNode(detkey_to_node[pair.first]);
        if (pair.first == 400) node->GetVolume()->SetLineColor(0);
        else node->GetVolume()->SetLineColor(30);
        TGeoSphere *sphere = (TGeoSphere *)(node->GetVolume()->GetShape());
        sphere->SetSphDimensions(0,20.*(pair.second/(best_vtx_tot_charge_out+best_vtx_tot_charge_in))+5.,0,180,0,360);
      }

/*
          double det_eff = m_pmt_eff[detkey];
          double alpha_i = -99.;
          alpha_i = m_pmt_alpha[detkey];

          // Calculate expected num of PE for each hit
          double exp_PE = alpha_i * (1.+TMath::Cos(0.)) / 2.;  //TODO:figure out angle that goes into cos
          h_expected_PE->Fill(exp_PE);

          //do something w charge and PE?
        }
      }
*/
      // check map of vtx and charge
      // XXX: there might be multiple clusters...
/*      std::cout << "  [debug] checking vertex and charge map" << std::endl;
      int pt = 0;
      for (std::map<int, double>::iterator it = m_charge_per_vtx.begin(); it != m_charge_per_vtx.end(); ++it)
      {
        std::cout << " [debug] vertex " << it->first << ", charge= " << it->second << std::endl;
        gr_vtx_charge_in->SetPoint(pt, vtxCandidates.at(it->first).Z(), it->second);
        ++pt;
      }*/

  } //end if found_muon


  std::cout << " [debug] Setting FittedTrackLengthInWater to CStore " << fitted_track_len << std::endl;
  m_data->CStore.Set("FittedTrackLengthInWater", fitted_track_len);
  std::cout << " [debug] Setting FittedMuonVertex to CStore" << std::endl;
  Position fitted_muon_vtx(fitted_vtx.X(), fitted_vtx.Y(), fitted_vtx.Z());
  m_data->CStore.Set("FittedMuonVertex", fitted_muon_vtx);

  double test_tracklen = -999.;
  m_data->CStore.Get("FittedTrackLengthInWater", test_tracklen);
  std::cout << " [debug] test retrieval of track length from store: " << test_tracklen << std::endl;
  Position test_vtx(-999, -999, -999);
  m_data->CStore.Get("FittedMuonVertex", test_vtx);
  std::cout << " [debug] test retrieval of fitted vertex from store: " << test_vtx.X() << "," << test_vtx.Y() << "," << test_vtx.Z() << std::endl;

  // --------------------------------------------------
  // --- Draw event -----------------------------------
  // --------------------------------------------------
  root_outp->cd();
  if (gr_vtx_charge_in->GetMaxSize() != 0) drawEvent = true;
  // Tell EventDisplay to draw event
  m_data->CStore.Set("DrawEventDisplay", drawEvent);

  // ANNIE in 3D
  if (plot3d && drawEvent)
  {
    if (verbosity > v_debug) { std::cout << "MuonFitter tool: Saving 3D plots..." << std::endl; }

    std::stringstream ss_evdisplay3d_title, ss_evdisplay3d_name;
    ss_evdisplay3d_title << "evdisplay3d_ev";
    ss_evdisplay3d_name << "canvas_evdisplay3d_p" << partnumber << "_ev";
    if (!isData)
    {
      ss_evdisplay3d_title << mcevnum;
      ss_evdisplay3d_name << mcevnum;
    }
    else
    {
      ss_evdisplay3d_title << evnum;
      ss_evdisplay3d_name << evnum;
    }
    canvas_3d->SetTitle(ss_evdisplay3d_title.str().c_str());
    canvas_3d->SetName(ss_evdisplay3d_name.str().c_str());
    ageom->DrawTracks();
    canvas_3d->Modified();
    canvas_3d->Update();
    canvas_3d->Write();
  }

  //saving graph
  //if (gr_vtx_charge_in->GetMaxSize() != 0)
  if (found_muon)
  {
    // Total charge
    c_vtx_charge->cd();
    std::stringstream ss_cvtxq_grtitle, ss_cvtxq_grname;
    ss_cvtxq_grtitle << "Total Charge Seen at Each Vertex Ev_";
    ss_cvtxq_grname << "c_vtx_charge_p" << partnumber << "_ev";
    if (isData)
    {
      ss_cvtxq_grtitle << evnum;
      ss_cvtxq_grname << evnum;
    }
    else
    {
      ss_cvtxq_grtitle << mcevnum;
      ss_cvtxq_grname << mcevnum;
    }
    gr_vtx_charge_in->Draw("alp");
    gr_vtx_charge_out->Draw("lp");
    c_vtx_charge->Modified();
    c_vtx_charge->Update();
    c_vtx_charge->SetTitle(ss_cvtxq_grtitle.str().c_str());
    c_vtx_charge->SetName(ss_cvtxq_grname.str().c_str());
    if (!isData)
    {
      TLine *lTrueVtxZ_total = new TLine(trueVtxZ, 0, trueVtxZ, gr_vtx_charge_in->GetPointY(vtxCandidates.size()-1));
      lTrueVtxZ_total->SetLineColor(4);
      lTrueVtxZ_total->SetLineWidth(2);
      lTrueVtxZ_total->Draw();
    }
    TLegend *legend1 = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend1->AddEntry(gr_vtx_charge_in, "inside cone", "lp");
    legend1->AddEntry(gr_vtx_charge_out, "outside cone", "lp");
    legend1->Draw();
    //c_vtx_charge->Write();

    // Charge per PMT
    c_charge_per_pmt->cd();
    std::stringstream ss_cv_qdensity_grtitle, ss_cv_qdensity_grname;
    ss_cv_qdensity_grtitle << "Charge Per PMT at Each Vertex Ev_";
    ss_cv_qdensity_grname << "c_charge_per_pmt_p" << partnumber << "_ev";
    if (isData)
    {
      ss_cv_qdensity_grtitle << evnum;
      ss_cv_qdensity_grname << evnum;
    }
    else
    {
      ss_cv_qdensity_grtitle << mcevnum;
      ss_cv_qdensity_grname << mcevnum;
    }
    gr_qdensity_in->Draw("alp");
    gr_qdensity_out->Draw("lp");
    //draw marker over largest charge per pmt value
    TMarker *mkr_max = new TMarker(vtxCandidates.at(max_vertex).Z(), max_vtx_charge, 43);
    mkr_max->SetMarkerColor(2);
    mkr_max->SetMarkerSize(1.5);
    mkr_max->Draw();
    c_charge_per_pmt->Modified();
    c_charge_per_pmt->Update();
    c_charge_per_pmt->SetTitle(ss_cv_qdensity_grtitle.str().c_str());
    c_charge_per_pmt->SetName(ss_cv_qdensity_grname.str().c_str());
    TLine *l1 = new TLine(l_zmin, avg_charge_per_pmt_in, l_zmax, avg_charge_per_pmt_in);
    l1->SetLineColor(8);
    l1->SetLineWidth(2);
    //l1->Draw();
    if (!isData)
    {
      TLine *lTrueVtxZ = new TLine(trueVtxZ, 0, trueVtxZ, max_vtx_charge);
      lTrueVtxZ->SetLineColor(4);
      lTrueVtxZ->SetLineWidth(2);
      lTrueVtxZ->Draw();
    }
    TLegend *legend2 = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend2->AddEntry(gr_qdensity_in, "inside cone", "lp");
    legend2->AddEntry(gr_qdensity_out, "outside cone", "lp");
    legend2->Draw();
    c_charge_per_pmt->Write();
    //cpp_file << ",,,,,," << avg_charge_per_pmt_in << "," << avg_charge_per_pmt_out << std::endl;

    // total detkeys
    c_vtx_detkey->cd();
    std::stringstream ss_vtx_detkey_title, ss_vtx_detkey_name;
    ss_vtx_detkey_title << "Num PMTs Seen at Each Vertex Ev_";
    ss_vtx_detkey_name << "c_vtx_detkey_p" << partnumber << "_ev";
    if (isData)
    {
      ss_vtx_detkey_title << evnum;
      ss_vtx_detkey_name << evnum;
    }
    else
    {
      ss_vtx_detkey_title << mcevnum;
      ss_vtx_detkey_name << mcevnum;
    }
    gr_vtx_detkey_out->Draw("alp");
    gr_vtx_detkey_in->Draw("lp");
    c_vtx_detkey->Modified();
    c_vtx_detkey->Update();
    c_vtx_detkey->SetTitle(ss_vtx_detkey_title.str().c_str());
    c_vtx_detkey->SetName(ss_vtx_detkey_name.str().c_str());
    TLegend *legend3 = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend3->AddEntry(gr_vtx_detkey_in, "inside cone", "lp");
    legend3->AddEntry(gr_vtx_detkey_out, "outside cone", "lp");
    legend3->Draw();
    c_vtx_detkey->Write();

    // area/frustum correction
    c_effarea_detkey->cd();
    std::stringstream ss_effarea_title, ss_effarea_name;
    ss_effarea_title << "Effective PMT Area Ev_";
    ss_effarea_name << "c_effarea_detkey_p" << partnumber << "_ev";
    if (isData)
    {
      ss_effarea_title << evnum;
      ss_effarea_name << evnum;
    }
    else
    {
      ss_effarea_title << mcevnum;
      ss_effarea_name << mcevnum;
    }
    gr_effarea_detkey->Draw("alp");
    c_effarea_detkey->Modified();
    c_effarea_detkey->Update();
    c_effarea_detkey->SetTitle(ss_effarea_title.str().c_str());
    c_effarea_detkey->SetName(ss_effarea_name.str().c_str());
    //c_effarea_detkey->Write();

    c_fpmt_detkey->cd();
    std::stringstream ss_fpmt_title, ss_fpmt_name;
    ss_fpmt_title << "Effective PMT Area as Fraction of Frustrum Area Ev_";
    ss_fpmt_name << "c_fpmt_detkey_p" << partnumber << "_ev";
    if (isData)
    {
      ss_fpmt_title << evnum;
      ss_fpmt_name << evnum;
    }
    else
    {
      ss_fpmt_title << mcevnum;
      ss_fpmt_name << mcevnum;
    }
    gr_fpmt_detkey->Draw("alp");
    c_fpmt_detkey->Modified();
    c_fpmt_detkey->Update();
    c_fpmt_detkey->SetTitle(ss_fpmt_title.str().c_str());
    c_fpmt_detkey->SetName(ss_fpmt_name.str().c_str());
    //c_fpmt_detkey->Write();

    c_effarea_ai->cd();
    std::stringstream ss_effarea_ai_name;
    ss_effarea_ai_name << "c_effarea_ai_p" << partnumber << "_ev";
    if (isData)
    {
      ss_effarea_ai_name << evnum;
    }
    else
    {
      ss_effarea_ai_name << mcevnum;
    }
    gr_effarea_ai->Draw("alp");
    c_effarea_ai->Modified();
    c_effarea_ai->Update();
    c_effarea_ai->SetTitle(ss_effarea_title.str().c_str());
    c_effarea_ai->SetName(ss_effarea_ai_name.str().c_str());
    c_effarea_ai->Write();

    c_fpmt_ai->cd();
    std::stringstream ss_fpmt_ai_name;
    ss_fpmt_ai_name << "c_fpmt_ai_p" << partnumber << "_ev";
    if (isData)
    {
      ss_fpmt_ai_name << evnum;
    }
    else
    {
      ss_fpmt_ai_name << mcevnum;
    }
    gr_fpmt_ai->Draw("alp");
    c_fpmt_ai->Modified();
    c_fpmt_ai->Update();
    c_fpmt_ai->SetTitle(ss_fpmt_title.str().c_str());
    c_fpmt_ai->SetName(ss_fpmt_ai_name.str().c_str());
    c_fpmt_ai->Write();

    c_eta_ai->cd();
    std::stringstream ss_eta_title, ss_eta_ai_name, ss_truetrack;
    ss_eta_title << "Eta: #PE Divided by f (tot PMT charge > " << PMTQCut << ") Ev_";
    ss_eta_ai_name << "c_eta_ai_p" << partnumber << "_ev";
    if (isData)
    {
      ss_eta_title << evnum;
      ss_eta_ai_name << evnum;
    }
    else
    {
      ss_eta_title << mcevnum;
      ss_eta_ai_name << mcevnum;
      ss_truetrack << "true track length in water: " << trueTrackLengthInWater << " cm";
      if (display_truth) gr_eta_ai->SetTitle(ss_truetrack.str().c_str());
    }
    gr_eta_ai->Draw("alp");
    gr_eta_ai->GetXaxis()->SetLimits(45., 505.);
    gr_eta_ai->GetHistogram()->SetMinimum(0.);
    gr_running_avg->Draw("lp");
    //gr_qincone_ai->Draw("lp");
    //gr_qoutcone_ai->Draw("lp");
    c_eta_ai->Modified();
    c_eta_ai->Update();
    c_eta_ai->SetTitle(ss_eta_title.str().c_str());
    c_eta_ai->SetName(ss_eta_ai_name.str().c_str());

    TLegend *legend4 = new TLegend(0.55, 0.65, 0.89, 0.89);
    legend4->AddEntry(gr_eta_ai, "#eta", "lp");
    //legend4->AddEntry(gr_qincone_ai, "inside cone", "lp");
    //legend4->AddEntry(gr_qoutcone_ai, "outside cone", "lp");
    legend4->AddEntry(gr_running_avg, "running avg #eta", "lp");

    TLine *lBestFitAi = new TLine(bestFitAi, 0, bestFitAi, max_eta+1);
    lBestFitAi->SetLineColor(2);
    lBestFitAi->SetLineWidth(2);
    //lBestFitAi->Draw();
    //legend4->AddEntry(lBestFitAi, "best fit", "l");
    if (!isData)
    { //indicate distance btwn trueVtx and tankExitPoint
      TLine *lTrueTankTrack = new TLine(trueTrackLengthInWater, 0, trueTrackLengthInWater, max_eta+1);
      lTrueTankTrack->SetLineColor(4);
      lTrueTankTrack->SetLineWidth(2);
      if (display_truth)
      {
        lTrueTankTrack->Draw();
        legend4->AddEntry(lTrueTankTrack, "true tank track length", "l");
      }

      TLine *lLeftAvgEta = new TLine(50, left_avg_eta, trueTrackLengthInWater, left_avg_eta);
      lLeftAvgEta->SetLineColor(2);
      lLeftAvgEta->SetLineWidth(2);
      if (display_truth)
      {
        lLeftAvgEta->Draw();
        legend4->AddEntry(lLeftAvgEta, "avg #eta to left of true track length", "l");
      }

      TLine *lRightAvgEta = new TLine(trueTrackLengthInWater, right_avg_eta, 455, right_avg_eta);
      lRightAvgEta->SetLineColor(8);
      lRightAvgEta->SetLineWidth(2);
      if (display_truth)
      {
        lRightAvgEta->Draw();
        legend4->AddEntry(lRightAvgEta, "avg #eta to right of true track length", "l");
      }
    }
    TLine *lAvgEta = new TLine(55, avg_eta, 500-deltaL/2, avg_eta);
    lAvgEta->SetLineColor(46);
    lAvgEta->SetLineWidth(2);
    lAvgEta->Draw();
    legend4->AddEntry(lAvgEta, "avg #eta", "l");
    legend4->Draw();
    c_eta_ai->Write();
    if (!isData && display_truth) ss_eta_ai_name << "_truth";
    ss_eta_ai_name << ".png";
    c_eta_ai->SaveAs(ss_eta_ai_name.str().c_str());

    
    c_h_tzero->cd();
    std::stringstream ss_t0_title, ss_t0_name;
    ss_t0_title << "t0" << " p" << partnumber << "_ev";
    ss_t0_name << "c_h_tzero_p" << partnumber << "_ev";
    if (isData)
    {
      ss_t0_title << evnum;
      ss_t0_name << evnum;
    }
    else
    {
      ss_t0_title << mcevnum;
      ss_t0_name << mcevnum;
    }
    h_tzero->Draw();
    c_h_tzero->Modified();
    c_h_tzero->Update();
    c_h_tzero->SetTitle(ss_t0_title.str().c_str());
    if (!isData && display_truth) ss_t0_name << "_truth";
    c_h_tzero->SetName(ss_t0_name.str().c_str());
    c_h_tzero->Write();
    ss_t0_name << ".png";
    if (fit_mode) c_h_tzero->SaveAs(ss_t0_name.str().c_str());
  }

  return true;
}


bool MuonFitter::Finalise(){
  // Save output
  root_outp->cd();
  h_alpha->Write();
  h_expected_PE->Write();
  h_phot_inc_angle->Write();
  //h_hit_angles->Write();
  //h_fitted_tank_track_len->Write();
  //h_closest_approach->Write();
  //h_num_mrd_layers->Write();
  //h_clusterhit_x->Write();
  //h_clusterhit_y->Write();
  //h_clusterhit_z->Write();
  //h_clusterhit_detkey->Write();
  //h_clusterhit_timespread->Write();
  //h_clusterhit_time->Write();
  if (!isData)
  {
    //h_lastvtx_z->Write();
    //h_truevtx_z->Write();
    //h_truevtx_angle->Write();
    h_true_tanktrack_len->Write();
    h_truevtx_x->Write();
    h_truevtx_y->Write();
    h_truevtx_z->Write();
    h_topview_truth->Write();
    h_sideview_truth->Write();
  }
  //h_tankexit_to_pmt->Write();
  //h_tanktrack_ai->Write();
  h_eff_area_pmt->Write();
  h_fpmt->Write();
  //h_eta_ai->Write();
  //h_qincone_truevtx->Write();
  //h_qoutcone_truevtx->Write();
  //h_total_pe_hits->Write();
  h_pmt_charge->Write();
  h_charge_detkey->Write();
  //h_tdiff->Write();
  if (fit_mode)
  {
    h_uber_t0widths->Write();
    h_fitted_track_len->Write();
    h_vtxfit_x->Write();
    h_vtxfit_y->Write();
    h_vtxfit_z->Write();
    h_topview_fit->Write();
    h_sideview_fit->Write();
    h_tank_mrd_E->Write();
    if (!isData)
    {
      h_truefit_len_diff->Write();
      h_truefitdiff_x->Write();
      h_truefitdiff_y->Write();
      h_truefitdiff_z->Write();
      h_deltaR->Write();
      h_transverse->Write();
      h_parallel->Write();
      for (int i = 1; i <= 75; ++i)
      {
        h_deltaR_4pi->SetBinContent(i, h_deltaR->GetBinContent(i)/(4.*TMath::Pi()*pow(2*i-1,2)*2.));
      }
      h_deltaR_4pi->Write();
      h_true_reco_E->Write();
      h_true_reco_Ediff->Write();
      h_Ediff_frac_tank->Write();
      h_Ediff_frac_mrd->Write();
      h_mrd_eloss_diff->Write();
      h_mrd_eloss_SBvANNIE->Write();
      h_true_reco_Ediff_sb->Write();
      h_true_reco_Ediff_sb_outerE->Write();
    }
  }

  // Close 3D canvas
  if (plot3d)
  {
    canvas_3d->Clear();
    canvas_3d->Close();
    delete canvas_3d;
  }

  delete c_vtx_charge;
  delete c_vtx_detkey;
  delete c_charge_per_pmt;
  delete c_effarea_detkey;
  delete c_fpmt_detkey;
  delete c_effarea_ai;
  delete c_fpmt_ai;
  delete c_eta_ai;
  delete c_h_tzero;

  root_outp->Close();
  pos_file.close();
  //cpp_file.close();
  pehits_file.close();
  if (!isData && display_truth) truetrack_file.close();

  Log("MuonFitter Tool: Exiting", v_message, verbosity);
  return true;
}

Position MuonFitter::Line3D(double x1, double y1, double z1, double x2, double y2, double z2, double C, std::string coord)
{
  // returns the XYZ coordinates of a new point in a 3D line
  // given two initial points and one X,Y,Z value
  double a = x2-x1;
  double b = y2-y1;
  double c = z2-z1;
  double x = -69., y = -69., z = -69., L = -69.;

  if (coord == "x" || coord == "X")
  {
    L = (C-x1)/a;
    x = C;
    y = b*L+y1;
    z = c*L+z1;
  }
  else if (coord == "y" || coord == "Y")
  {
    L = (C-y1)/b;
    x = a*L+x1;
    y = C;
    z = c*L+z1;
  }
  else if (coord == "z" || coord == "Z")
  {
    L = (C-z1)/c;
    x = a*L+x1;
    y = b*L+y1;
    z = C;
  }
  else { Log("MuonFitter Tool: Unable to recognize coordinate!", v_debug, verbosity); }

  return(Position(x,y,z));
}

void MuonFitter::reset_3d()
{
  //reset TANK PMTs
  TGeoNode *node;
  TGeoSphere *sphere;
  for (unsigned long detkey=332; detkey<464; detkey++)
  {
    node = EXPH->GetNode(detkey_to_node[detkey]);
    sphere = (TGeoSphere *)(node->GetVolume()->GetShape());
    sphere->SetSphDimensions(0,1.5,0,180,0,180);
    node->GetVolume()->SetLineColor(41);
  }

  //remove extra node(s)
  if (N > maxN)
  {
    while (N > maxN)
    {
      node = EXPH->GetNode(maxN);
      EXPH->RemoveNode(node);
      --N;
    }
  }

  //remove tracks
  if(ageom->GetNtracks() > 0)
  {
    ageom->ClearTracks();
    if (verbosity > v_debug) { std::cout << "Clearing tracks... current number: " << ageom->GetNtracks() << std::endl; }
  }

  canvas_3d->Modified();
  canvas_3d->Update();

}

bool MuonFitter::FileExists(std::string fname)
{
  ifstream myfile(fname.c_str());
  return myfile.good();
}

void MuonFitter::LoadTankTrackFits()
{
  ifstream myfile(tankTrackFitFile.c_str());
  std::cout << "opened file" << std::endl;
  std::string line;
  if (myfile.is_open())
  {
    // Loop over lines, get event id, fit (should only be one line here)
    while (getline(myfile, line))
    {
      //if (verbosity > 3) std::cout << line << std::endl; //has our stuff
      if (line.find("#") != std::string::npos) continue;
      std::vector<std::string> fileLine;
      boost::split(fileLine, line, boost::is_any_of(","), boost::token_compress_on);
      std::string event_id = "";
      double tank_track_fit = -9999.;
      event_id = fileLine.at(0);
      tank_track_fit = std::stod(fileLine.at(1));
      //if (verbosity > 4) std::cout << "event_id, tank_track_fit: " << event_id << " , " << tank_track_fit << std::endl;
      m_tank_track_fits.insert(std::pair<std::string,double>(event_id, tank_track_fit));
    }
  }
  return;
}

double MuonFitter::Calc_dEdx(double muon_T0)
{
  double gamma = muon_T0/105.66;
  double beta = TMath::Sqrt(gamma*gamma-1.)/gamma;
  double Tmax = (2.*0.511*beta*beta*gamma*gamma)/(1.+2.*0.511*gamma/105.66 + pow(0.511/105.66,2));
  double ox_ln = TMath::Log(2.*0.511*beta*beta*gamma*gamma*Tmax/(I_O*I_O));
  double ox_dEdx = (0.999)*(16./18.)*1*0.307*0.5*pow(1./beta,2)*(0.5*ox_ln-beta*beta-0.5);
  double h_ln = TMath::Log(2.*0.511*beta*beta*gamma*gamma*Tmax/(I_H*I_H));
  double h_dEdx = (0.999)*(2./18.)*1*0.307*0.5*pow(1./beta,2)*(0.5*h_ln-beta*beta-0.5);

  // including Gd2(SO4)3
  double gd_ln = TMath::Log(2.*0.511*beta*beta*gamma*gamma*Tmax/(I_GD*I_GD));
  double gd_dEdx = (0.001)*(314./602.)*1*0.307*0.5*pow(1./beta,2)*(0.5*gd_ln-beta*beta-0.5);
  double s_ln = TMath::Log(2.*0.511*beta*beta*gamma*gamma*Tmax/(I_S*I_S));
  double s_dEdx = (0.001)*(96./602.)*1*0.307*0.5*pow(1./beta,2)*(0.5*s_ln-beta*beta-0.5);
  double gdox_dEdx = (0.001)*(192./602.)*1*0.307*0.5*pow(1./beta,2)*(0.5*ox_ln-beta*beta-0.5);

  return (ox_dEdx + h_dEdx + gd_dEdx + s_dEdx + gdox_dEdx);
}

double MuonFitter::Calc_MrdEnergyLoss(double rng)
{
  // SciBooNE method of determining energy loss in MRD
  //assumes rng is in cm
  //returns energy in MeV
  return ((0.0982739 + 0.013045*rng)*1.0e3);
}
