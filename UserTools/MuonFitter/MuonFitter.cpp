#include "MuonFitter.h"

/* *****************************************************************
 * tool name: MuonFitter
 * author: Julie He
 *
 * desc: This tool takes
 *
 * versions:
 * ...
 * 240331v1JH: fixed line of code defining reco_mu_e (prev defined 
 *              using mrdEnergyLoss)
 * 240401v2JH: removed unnecessary code, old algorithms
 * 240405v1JH: add option to use simple energy reco (just add eloss)
 * 240407v1JH: reco mode no longer includes finding muon, added 
 *              cluster_time to m_tank_track_fits
 * 240425v1JH: incorporated scattering effects into MRD track length
 * *****************************************************************
 */

MuonFitter::MuonFitter():Tool(){}


bool MuonFitter::Initialise(std::string configfile, DataModel &data){

  ///////////////////////// Useful header /////////////////////////
  //-- Load config file
  if(configfile!="") m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  //-- Set default values
  verbosity = 3;
  isData = true;
  PMTMRDOffset = 745.;

  // -------------------------------------------------------------
  // --- Retrieve config variables -------------------------------
  // -------------------------------------------------------------
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
  m_variables.Get("StepSizeAi", step_size_ai);
  m_variables.Get("InsideAngle", insideAngle);
  m_variables.Get("OutsideAngle", outsideAngle);
  m_variables.Get("PMTChargeThreshold", PMTQCut);
  m_variables.Get("EtaThreshold", EtaThreshold);
  m_variables.Get("DisplayTruth", display_truth);
  m_variables.Get("RecoMode", reco_mode);
  m_variables.Get("AiEtaFile", aiEtaFile);
  m_variables.Get("TankTrackFitFile", tankTrackFitFile);
  m_variables.Get("UseNumLayers", use_nlyrs);
  m_variables.Get("UsePCA", use_pca);
  m_variables.Get("UseConnDots", use_conn_dots);
  m_variables.Get("UseELoss", use_eloss);
  m_variables.Get("UseSimpleEReco", use_simple_ereco);
  m_variables.Get("RecoEnergyShift", ERECO_SHIFT);

  if (use_nlyrs) Log("MuonFitter Tool: Using num layers to determine MRD track length", v_message, verbosity);
  if (use_pca) Log("MuonFitter Tool: Using PCA to determine MRD track angle", v_message, verbosity);
  if (use_conn_dots) Log("MuonFitter Tool: Using connect the dots method to determine MRD track angle", v_message, verbosity);
  if (use_eloss) Log("MuonFitter Tool: Using current ANNIE tools to determine MRD energy loss", v_message, verbosity);
  if (use_simple_ereco) Log("MuonFitter Tool: Just add tank and MRD energy depositions (don't update dEdx values)", v_message, verbosity);


  // Output ROOT file
  std::string outfile;
  m_variables.Get("OutputFile", outfile);
  Log("MuonFitter Tool: Saving output into " + outfile, v_message, verbosity);
  root_outp = new TFile(outfile.c_str(), "RECREATE");


  // -------------------------------------------------------------
  // --- Initialize canvases, graphs, histograms -----------------
  // -------------------------------------------------------------
  //-- Canvases
  c_effarea_detkey = new TCanvas("c_effarea_detkey", "Effective PMT Area vs detkey", 800, 600);
  c_fpmt_detkey = new TCanvas("c_fpmt_detkey", "Effective PMT Area as Fraction of Frustum Area", 800, 600);
  c_effarea_ai = new TCanvas("c_effarea_ai", "Effective PMT Area vs a_{i} (tank track segment)", 800, 600);
  c_fpmt_ai = new TCanvas("c_fpmt_ai", "Effective PMT Area as Fraction of Frustum Area", 800, 600);
  c_eta_ai = new TCanvas("c_eta_ai", "Charge per cm (#eta)", 800, 600);
  c_h_tzero = new TCanvas("c_h_tzero", "t0", 800, 600);

  //-- Histograms
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

  h_tanktrack_ai = new TH1D("h_tanktrack_ai", "Segment of Tank Track (a_{i})", 250, 0., 500.);
  h_tanktrack_ai->GetXaxis()->SetTitle("a_{i} [cm]");

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

  h_fitted_tank_track = new TH1D("h_fitted_tank_track", "Distribution of Fitted Track Length", 200, 0, 400);
  h_fitted_tank_track->GetXaxis()->SetTitle("track length [cm]");

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

  h_true_reco_E = new TH2D("h_true_reco_E", "Reco vs True #mu Energy", 250, 0., 5000., 250, 0., 5000.);
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

  h_mrd_eloss_diff = new TH1D("h_mrd_eloss_diff", "diff in MRD energy loss", 200, -500, 1500);
  h_mrd_eloss_diff->GetXaxis()->SetTitle("[MeV]");

  h_tank_track_diff_small = new TH1D("h_tank_track_diff_small", "Difference in TANK track lengths (Reco-True) When E_{diff} < 200 MeV", 150, -150, 150);
  h_tank_track_diff_small->GetXaxis()->SetTitle("reco-truth [cm]");

  h_tank_track_diff_large = new TH1D("h_tank_track_diff_large", "Difference in TANK track lengths (Reco-True) When E_{diff} > 200 MeV", 150, -150, 150);
  h_tank_track_diff_large->GetXaxis()->SetTitle("reco-truth [cm]");

  h_mrd_track_diff_small = new TH1D("h_mrd_track_diff_small", "Difference in MRD track lengths (Reco-True) When E_{diff} < 200 MeV", 150, -150, 150);
  h_mrd_track_diff_small->GetXaxis()->SetTitle("reco-truth [cm]");

  h_mrd_track_diff_large = new TH1D("h_mrd_track_diff_large", "Difference in MRD track lengths (Reco-True) When E_{diff} > 200 MeV", 150, -150, 150);
  h_mrd_track_diff_large->GetXaxis()->SetTitle("reco-truth [cm]");

  h_deltaR_small = new TH1D("h_deltaR_small", "#DeltaR When E_{diff} < 200 MeV", 76, 0, 152);
  h_deltaR_small->GetXaxis()->SetTitle("#DeltaR [cm]");

  h_deltaR_large = new TH1D("h_deltaR_large", "#DeltaR When E_{diff} > 200 MeV", 76, 0, 152);
  h_deltaR_large->GetXaxis()->SetTitle("#DeltaR [cm]");

  h_mrd_nlyrs_reco = new TH2D("h_mrd_nlyrs_reco", "MRD Track Length: Num Layers vs Current Reco", 40, 0., 200., 20, 0., 100.);
  h_mrd_nlyrs_reco->GetXaxis()->SetTitle("current mrd reco [cm]");
  h_mrd_nlyrs_reco->GetYaxis()->SetTitle("num layers [cm]");

  h_mrd_track_diff = new TH1D("h_mrd_track_diff", "Difference in MRD Track Lengths (True - Reco)", 100, -250, 250);
  h_mrd_track_diff->GetXaxis()->SetTitle("reco-true [cm]");

  h_mrd_track_diff_nlyrs = new TH1D("h_mrd_track_diff_nlyrs", "Difference in MRD Track Lengths (True - Reco nlyrs)", 100, -250, 250);
  h_mrd_track_diff_nlyrs->GetXaxis()->SetTitle("reco nlyrs-true [cm]");

  h_mrd_angle_diff = new TH1D("h_mrd_angle_diff", "Difference in Muon Direction Angle (Reco-True)", 300, -30, 30);
  h_mrd_angle_diff->GetXaxis()->SetTitle("reco-true [deg]");

  h_mrd_angle = new TH1D("h_mrd_angle", "Reconstructed Muon Direction Angle Close to True Angle", 200, -50, 50);
  h_mrd_angle->GetXaxis()->SetTitle("reco angle [deg]");

  h_clusterPE = new TH1D("h_clusterPE", "Cluster PE of Events with Found Muon", 100, 0, 10000);
  h_clusterPE->GetXaxis()->SetTitle("cluster PE [p.e.]");

  h_clusterPE_fit = new TH1D("h_clusterPE_fit", "Cluster PE of Events with Tank Track Fits", 100, 0, 10000);
  h_clusterPE_fit->GetXaxis()->SetTitle("cluster PE [p.e.]");

  h_clusterPE_fit_haspion = new TH1D("h_clusterPE_fit_haspion", "Cluster PE of Events with Tank Track Fits", 100, 0, 10000);
  h_clusterPE_fit_haspion->GetXaxis()->SetTitle("cluster PE [p.e.]");

  h_clusterPE_lrg_ediff = new TH1D("h_clusterPE_lrg_ediff", "Cluster PE of Events with Large E Diff", 100, 0, 10000);
  h_clusterPE_lrg_ediff->GetXaxis()->SetTitle("cluster PE [p.e.]");

  h_clusterPE_lrg_ediff_haspion = new TH1D("h_clusterPE_lrg_ediff_haspion", "Cluster PE of Events with Large E Diff and Has Pion", 100, 0, 10000);
  h_clusterPE_lrg_ediff_haspion->GetXaxis()->SetTitle("cluster PE [p.e.]");

  h_pca_angle = new TH1D("h_pca_angle", "Reconstructed Muon Direction Angle Using PCA", 200, -50, 50);
  h_pca_angle->GetXaxis()->SetTitle("pca angle [deg]");

  h_pca_reco_angle = new TH1D("h_pca_reco_angle", "PCA vs Reco Angle", 200, -50, 50);
  h_pca_reco_angle->GetXaxis()->SetTitle("pca-reco angle [deg]");

  h_pca_true_angle = new TH1D("h_pca_true_angle", "PCA vs True Angle", 200, -50, 50);
  h_pca_true_angle->GetXaxis()->SetTitle("pca-true angle [deg]");

  h_total_track_diff = new TH1D("h_total_track_diff", "Difference in Total Track Lengths (True - Reco)", 100, -250, 250);
  h_total_track_diff->GetXaxis()->SetTitle("reco-true [cm]");

  h_total_track_diff_nlyrs = new TH1D("h_total_track_diff_nlyrs", "Difference in MRD Track Lengths (True - Reco nlyrs)", 100, -250, 250);
  h_total_track_diff_nlyrs->GetXaxis()->SetTitle("reco nlyrs-true [cm]");

  h_remainder_track_last20MeV = new TH1D("h_remainder_track_last20MeV", "Remaining Track Length in MRD When Input E_{#mu} < 20 MeV", 25, 0, 50);
  h_remainder_track_last20MeV->GetXaxis()->SetTitle("remaining MRD track [cm]");

  h_true_reco_Ediff_last20MeV = new TH1D("h_true_reco_Ediff_last20MeV", "Energy Difference Btwn True and Reco Energy When Input E_{#mu} < 20 MeV", 100, -500, 500);
  h_true_reco_Ediff_last20MeV->GetXaxis()->SetTitle("E diff (reco-true) [MeV]");

  h_remainder_track_Ediff_last20MeV = new TH2D("h_remainder_track_Ediff_last20MeV", "Remaining MRD Track Length vs True/Reco Energy Diff When Input E_{#mu} < 20 MeV", 100, -500, 500, 25, 0, 50);
  h_remainder_track_Ediff_last20MeV->GetXaxis()->SetTitle("E diff (reco-true) [MeV]");
  h_remainder_track_Ediff_last20MeV->GetYaxis()->SetTitle("remaining MRD track [cm]");

  //-- Graphs
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
  gr_eta_ai->GetYaxis()->SetTitle("#eta (n_{PE}/f)");

  gr_running_avg = new TGraph();
  gr_running_avg->SetLineWidth(0);
  gr_running_avg->SetMarkerStyle(8);
  gr_running_avg->SetMarkerColor(38);
  gr_running_avg->SetFillColor(0);
  gr_running_avg->GetXaxis()->SetTitle("a_{i} [cm]");
  gr_running_avg->GetYaxis()->SetTitle("#eta (nPE/f)");


  // -------------------------------------------------------------
  // --- Retrieve Store variables --------------------------------
  // -------------------------------------------------------------
  m_data->CStore.Get("ChannelNumToTankPMTSPEChargeMap", ChannelKeyToSPEMap);  //same for data and mc


  // -------------------------------------------------------------
  // --- Get Geometry --------------------------------------------
  // -------------------------------------------------------------
  //-- Code based on UserTools/EventDisplay/EventDisplay.cpp
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
  //-- Use compressed barrel radius for ANNIEp2v6 detector configuration (only MC)
  if (detector_config == "ANNIEp2v6" && !isData) { tank_height *= barrel_compression; }
  else if (isData) { tank_height = 1.2833; }
  //-- Set tank radius to standard value of old anniev2 configuration (v4/v6 seems to have very different radius?)
  if (tank_radius < 1. || isData) { tank_radius = 1.37504; }

  Position detector_center = geom->GetTankCentre();
  tank_center_x = 100.*detector_center.X();  //convert to cm
  tank_center_y = 100.*detector_center.Y();
  tank_center_z = 100.*detector_center.Z();
  //-- QA: Check tank center, radius, and height being used
  std::cout << " MuonFitter Tool: tank_center xyz [cm]: " << tank_center_x << "," << tank_center_y << "," << tank_center_z << std::endl;
  std::cout << " MuonFitter Tool: tank_radius [m]: " << tank_radius << std::endl;
  std::cout << " MuonFitter Tool: tank_height [m]: " << tank_height << std::endl;

  //-- QA: Check number of PMTs in each subdetector
  n_tank_pmts = geom->GetNumDetectorsInSet("Tank");
  n_mrd_pmts = geom->GetNumDetectorsInSet("MRD");
  n_veto_pmts = geom->GetNumDetectorsInSet("Veto");
  std::cout << " MuonFitter Tool: Number of Tank / MRD / Veto PMTs in this geometry: " << n_tank_pmts << " / " << n_mrd_pmts << " / " << n_veto_pmts << std::endl;


  // -------------------------------------------------------------
  // --- ANNIE in 3D ---------------------------------------------
  // -------------------------------------------------------------
  // --- Set up 3D geometry for viewing 
  Log("MuonFitter Tool: Creating 3D Geometry", v_debug, verbosity);
  ageom = new TGeoManager("ageom", "ANNIE in 3D");
  TGeoNode *node;

  //-- material
  vacuum = new TGeoMaterial("vacuum",0,0,0);
  Fe = new TGeoMaterial("Fe",55.845,26,7.87);

  //-- create media
  Air = new TGeoMedium("Vacuum",0,vacuum);
  Iron = new TGeoMedium("Iron",1,Fe);

  //-- create volume
  EXPH = ageom->MakeBox("EXPH",Air,300,300,300);
  ageom->SetTopVolume(EXPH);
  ageom->SetTopVisible(0);
  //-- If you want to see the boundary, input the number 1 instead of 0:
  // geom->SetTopVisible(1);

  //-- draw CANVAS/TOP CENTER
  bBlock = ageom->MakeSphere("EXPH_vol_center", Iron, 0,3,0,180,0,360);
  bBlock->SetLineColor(1);
  EXPH->AddNodeOverlap(bBlock,N++,new TGeoTranslation(0,0,0));

  //-- draw TANK
  TGeoVolume *annietank = ageom->MakeTubs("annietank", Iron, 0, tank_radius*100., tank_height*100., 0, 360);  //convert to cm
  annietank->SetLineColor(38);
  EXPH->AddNodeOverlap(annietank,N++,new TGeoCombiTrans(0,0,0,new TGeoRotation("annietank",0,90,0)));
  node = EXPH->GetNode(N-1);


  // -------------------------------------------------------------
  // --- Read in TANK PMTs ---------------------------------------
  // -------------------------------------------------------------
  //-- Code based on UserTools/EventDisplay/EventDisplay.cpp
  std::map<std::string, std::map<unsigned long, Detector *>> *Detectors = geom->GetDetectors();

  Log("MuonFitter Tool: Adding tank PMTs to 3D geometry", v_debug, verbosity);
  for (std::map<unsigned long, Detector *>::iterator it = Detectors->at("Tank").begin(); it != Detectors->at("Tank").end(); ++it)
  {
    Detector *apmt = it->second;
    unsigned long detkey = it->first;
    std::string det_type = apmt->GetDetectorType();
    Position position_PMT = apmt->GetDetectorPosition();  //in meters
    //-- PMT xyz corrected by tank center, in cm
    x_pmt.insert(std::pair<int,double>(detkey, 100.*position_PMT.X()-tank_center_x));
    y_pmt.insert(std::pair<int,double>(detkey, 100.*position_PMT.Y()-tank_center_y));
    z_pmt.insert(std::pair<int,double>(detkey, 100.*position_PMT.Z()-tank_center_z));
    //-- PMT orientation
    Direction direction_PMT = apmt->GetDetectorDirection();
    x_pmt_dir.insert(std::pair<unsigned long, double>(detkey, direction_PMT.X()));
    y_pmt_dir.insert(std::pair<unsigned long, double>(detkey, direction_PMT.Y()));
    z_pmt_dir.insert(std::pair<unsigned long, double>(detkey, direction_PMT.Z()));

    //-- Get PMT areas for ring area overlap later
    double det_area = 0.;
    if (det_type == "LUX" || det_type == "ANNIEp2v7-glassFaceWCPMT_R7081") { det_area = LUX_AREA; }
    else if (det_type == "ETEL" || det_type == "ANNIEp2v7-glassFaceWCPMT_D784KFLB") { det_area = ETEL_AREA; }
    else if (det_type == "Hamamatsu" || det_type == "ANNIEp2v7-glassFaceWCPMT_R5912HQE") { det_area = HAMAMATSU_AREA; }
    else if (det_type == "Watchboy" || det_type == "ANNIEp2v7-glassFaceWCPMT_R7081") { det_area = WATCHBOY_AREA; }
    else if (det_type == "Watchman" || det_type == "ANNIEp2v7-glassFaceWCPMT_R7081HQE") { det_area = WATCHMAN_AREA; }
    else { Log("MuonFitter Tool: Unrecognized detector type! Setting det_area to 0.", v_error, verbosity); }
    m_pmt_area.insert(std::pair<int, double>(detkey, det_area));

    //-- ANNIE in 3D: Drawing TANK PMTs
    sprintf(blockName,"tank_pmt%lu",detkey);
    bBlock = ageom->MakeSphere(blockName, Iron, 0,1.5,0,180,0,360);
    bBlock->SetLineColor(41);
    EXPH->AddNodeOverlap(bBlock,1,new TGeoTranslation(x_pmt[detkey], y_pmt[detkey], z_pmt[detkey]));
    detkey_to_node.insert(std::pair<unsigned long, int>(detkey, N++));
  }


  // -------------------------------------------------------------
  // --- Read in MRD PMTs ----------------------------------------
  // -------------------------------------------------------------
  //-- Code based on UserTools/EventDisplay/EventDisplay.cpp
  Log("MuonFitter Tool: Adding MRD PMTs to 3D geoemtry", v_debug, verbosity);
  for (std::map<unsigned long, Detector *>::iterator it = Detectors->at("MRD").begin(); it != Detectors->at("MRD").end(); ++it)
  {
    Detector *amrdpmt = it->second;
    unsigned long detkey = it->first;
    unsigned long chankey = amrdpmt->GetChannels()->begin()->first;
    Paddle *mrdpaddle = (Paddle *)geom->GetDetectorPaddle(detkey);

    //-- Retrieve xyz bounds of paddle, orientation, etc
    //-- NOTE: so far this has not been used
    double xmin = mrdpaddle->GetXmin();
    double xmax = mrdpaddle->GetXmax();
    double ymin = mrdpaddle->GetYmin();
    double ymax = mrdpaddle->GetYmax();
    double zmin = mrdpaddle->GetZmin();
    double zmax = mrdpaddle->GetZmax();
    int orientation = mrdpaddle->GetOrientation();    //0:horizontal,1:vertical
    int half = mrdpaddle->GetHalf();    //0 or 1
    int side = mrdpaddle->GetSide();

    std::vector<double> xdim{xmin,xmax};
    std::vector<double> ydim{ymin,ymax};
    std::vector<double> zdim{zmin,zmax};

    //-- Store the detkey and the bounds of the MRD paddle
    //-- NOTE: mrd_{xyz} is type std::vector<unsigned long, std::vector<double>>
    mrd_x.emplace(detkey,xdim);
    mrd_y.emplace(detkey,ydim);
    mrd_z.emplace(detkey,zdim);

    //-- Store detkey and paddle center xyz (tank center corrected)
    mrd_center_x.emplace(detkey,100.*(mrdpaddle->GetOrigin()).X()-tank_center_x);
    mrd_center_y.emplace(detkey,100.*(mrdpaddle->GetOrigin()).Y()-tank_center_y);
    mrd_center_z.emplace(detkey,100.*(mrdpaddle->GetOrigin()).Z()-tank_center_z);


    //-- ANNIE in 3D: drawing MRD
    if (draw3d_mrd)
    {
      Position position_MRD = mrdpaddle->GetOrigin();
      sprintf(blockName,"mrd_pmt%lu",detkey);
      bBlock = ageom->MakeBox(blockName, Iron, (xmax-xmin)*100., (ymax-ymin)*100., (zmax-zmin)*100.);
      bBlock->SetLineColor(16);
      EXPH->AddNodeOverlap(bBlock,detkey,new TGeoTranslation(100.*position_MRD.X()-tank_center_x, 100.*position_MRD.Y()-tank_center_y, 100.*position_MRD.Z()-tank_center_z));
      detkey_to_node.insert(std::pair<unsigned long, int>(detkey, N++));

      //-- QA: Check detkey_to_node map
      if (verbosity > v_debug)
      {
        std::cout << " [debug] blockName: " << blockName << std::endl;
        std::cout << " [debug] detkey_to_node[detkey]: " << detkey_to_node[detkey] << std::endl;
      }
    }
  }

  // ------------------------------------------------------------
  // --- Read in FMV PMTs ---------------------------------------
  // ------------------------------------------------------------
  //-- Code based on UserTools/EventDisplay/EventDisplay.cpp
  Log("MuonFitter Tool: Adding FMV PMTs to 3D geoemtry", v_debug, verbosity);
  for (std::map<unsigned long, Detector *>::iterator it = Detectors->at("Veto").begin(); it != Detectors->at("Veto").end(); ++it)
  {
    Detector *avetopmt = it->second;
    unsigned long detkey = it->first;
    unsigned long chankey = avetopmt->GetChannels()->begin()->first;
    Paddle *vetopaddle = (Paddle *)geom->GetDetectorPaddle(detkey);

    double xmin = vetopaddle->GetXmin();
    double xmax = vetopaddle->GetXmax();
    double ymin = vetopaddle->GetYmin();
    double ymax = vetopaddle->GetYmax();
    double zmin = vetopaddle->GetZmin();
    double zmax = vetopaddle->GetZmax();

    //-- ANNIE in 3D: drawing FMV
    if (draw3d_fmv)
    {
      Position position_FMV = vetopaddle->GetOrigin();
      sprintf(blockName,"fmv_pmt%lu",detkey);
      bBlock = ageom->MakeBox(blockName, Iron, (xmax-xmin)*100., (ymax-ymin)*100., (zmax-zmin)*100.);
      bBlock->SetLineColor(20);
      EXPH->AddNodeOverlap(bBlock,detkey,new TGeoTranslation(100.*position_FMV.X()-tank_center_x, 100.*position_FMV.Y()-tank_center_y, 100.*position_FMV.Z()-tank_center_z));
      detkey_to_node.insert(std::pair<unsigned long, int>(detkey, N++));

      //-- QA: Check detkey_to_node map
      if (verbosity > v_debug)
      {
        std::cout << " [debug] blockName: " << blockName << std::endl;
        std::cout << " [debug] detkey_to_node[detkey]: " << detkey_to_node[detkey] << std::endl;
      }
    }
  }

  //-- ANNIE in 3D: Set max number of nodes
  maxN = N;
  Log("MuonFitter Tool: Number of nodes in 3D geometry: " + std::to_string(maxN), v_debug, verbosity);


  //-- Save 3D plots
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

  // ------------------------------------------------------------
  // --- Create txt files for separate analyses -----------------
  // ------------------------------------------------------------
  //-- Save start & stop vertices to file; Save (ai,eta) values
  //std::string pos_fname = "posFile.txt";
  //pos_file.open(pos_fname.c_str());
  //pos_file << "##evnum,startX,startY,startZ,stopX,stopY,stopZ" << std::endl;
  //std::string pos_fname = "ev_ai_eta.txt";
  std::string pos_fname = aiEtaFile.c_str();
  pos_file.open(pos_fname.c_str(), std::ostream::app);
  if (!pos_file)
  {
    std::cout << " [pos_file] File " << pos_fname << " does not exist! Creating now.." << std::endl;
    pos_file.open(pos_fname.c_str());
    pos_file << "##ev_id,cluster_time,ai,eta" << std::endl;
  }

  //-- Save charge data to file
  //std::string cpp_fname = "charge_per_pmt.txt";
  //cpp_file.open(cpp_fname.c_str());
  //cpp_file << "##evnum,nVtx,Qin,Qout,avgSumIn,avgSumOut,avgQin,avgQout" << std::endl;

  /*std::string pehits_fname = "tot_pe_hits.txt";
  pehits_file.open(pehits_fname.c_str(), std::ostream::app);
  if (!pehits_file)
  {
    std::cout << "File " << pehits_fname << " does not exist! Creating now..." << std::endl;
    pehits_file.open(pehits_fname.c_str());
    pehits_file << "##partfile,main_cluster_hits,main_cluster_charge" << std::endl;
  }*/

  //-- Save avg eta to left and right of true tank track length
  if (!isData)
  {
    truetrack_file.open("true_track_len.txt", std::ostream::app);
    if (!truetrack_file)
    {
      std::cout << " [debug] true_track_len.txt does not exist! Creating now..." << std::endl;
      truetrack_file.open("true_track_len.txt");
      truetrack_file << "#event_id,true_track_len,left_avg,right_avg" << std::endl;
    }
  }

  //-- Save fitted_tank_track, nhits, nhits_incone to file
  //-- TODO: Get rid of this
  /*nhits_trlen_file.open("nhits_trlen.txt", std::ostream::app);
  if (!nhits_trlen_file)
  {
    std::cout << "File nhits_trlen.txt does not exist! Creating now..." << std::endl;
    nhits_trlen_file.open("nhits_trlen.txt");
    nhits_trlen_file << "##event_id,track_fit,nhits,nhits_incone,totalpe,totalpe_incone" << std::endl;
  }*/

  //-- Save info about events with Ediff > 200 MeV
  if (!isData)
  {
    lg_ediff_file.open("lg_ediff.txt", std::ostream::app);
    if (!lg_ediff_file)
    {
      std::cout << "File lg_ediff.txt does not exist! Creating now..." << std::endl;
      lg_ediff_file.open("lg_ediff.txt");
      lg_ediff_file << "##event_id,ediff,pions,tanktrackF,tanktrackT,mrdtrackF,mrdtrackT,muonEF,muonET" << std::endl;
    }
  }

  // ------------------------------------------------------------
  // --- Load vertex fits (RECO MODE ONLY) ----------------------
  // ------------------------------------------------------------
  if (reco_mode)
  {
    if (this->FileExists(tankTrackFitFile))
    {
      //-- Load event part, number and track fit
      std::cout << "MuonFitter Tool: Loading tank track fits from " << tankTrackFitFile << std::endl;
      this->LoadTankTrackFits();
    }
    else
    {
      Log("MuonFitter Tool: File for tank track fit does not exist! Continuing tool NOT in RECO MODE.", v_error, verbosity);
      reco_mode = false;
    }
  }

  Log("MuonFitter Tool: Initialization complete", v_debug, verbosity);

  return true;
}


bool MuonFitter::Execute(){
  Log("MuonFitter Tool: Executing", v_debug, verbosity);

  //-- Retrieve ANNIEEvent
  int eventExists = m_data->Stores.count("ANNIEEvent");
  if (!eventExists)
  {
    Log("MuonFitter Tool: No ANNIEEvent store!", v_error, verbosity);
    return false;
  }

  // ------------------------------------------------------------
  // --- Reset variables ----------------------------------------
  // ------------------------------------------------------------
  //this->ResetVariables();   //TODO:write ResetVariables func
  bool drawEvent = false;   //-- Make graphs for this event (e.g. ev displays, eta vs ai)
  if (plot3d) { reset_3d(); }

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

  //-- Initialize tank track length and muon vertex to be saved in
  //-- CStore for downstream tools
  m_data->CStore.Set("FittedTrackLengthInWater", -888.);
  Position dummy_vtx(-888, -888, -888);
  m_data->CStore.Set("FittedMuonVertex", dummy_vtx);
  m_data->CStore.Set("RecoMuonKE", -888.);


  // -------------------------------------------------------------
  // --- Get event info (DATA) -----------------------------------
  // -------------------------------------------------------------
  int get_ok = false;
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

  // -------------------------------------------------------------
  // --- Get event info (MC) -------------------------------------
  // -------------------------------------------------------------
  TVector3 trueTrackDir;
  RecoVertex *truevtx = 0;
  if (!isData)
  {
    get_ok = m_data->Stores["ANNIEEvent"]->Get("MCParticles", mcParticles);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving MCParticles from ANNIEEvent!", v_error, verbosity); return false; }  //<< needed to retrieve true vertex and direction
    get_ok = m_data->Stores["ANNIEEvent"]->Get("MCEventNum", mcevnum);
    std::cout << "MuonFitter Tool: Working on event " << mcevnum << std::endl;
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving MCEventNum from ANNIEEvent!", v_error, verbosity); return false; }
    get_ok = m_data->Stores["ANNIEEvent"]->Get("MCTriggernum", mctrignum);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving MCTriggernum from ANNIEEvent!", v_error, verbosity); return false; }
    get_ok = m_data->Stores["ANNIEEvent"]->Get("MCFile", mcFile);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving MCFile from ANNIEEvent!", v_error, verbosity); mcFile = "-1"; }
 
    //-- Extract file part number
    //-- NOTE: event numbers weren't incrementing properly...
    std::string delim = ".";
    std::string tmp_str = mcFile.erase(0, mcFile.find(delim) + delim.length());
    partnumber = stoi(tmp_str.substr(0, tmp_str.find(delim)));

    //-- Get RecoEvent variables
    get_ok = m_data->Stores["RecoEvent"]->Get("TrueVertex", truevtx);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving TrueVertex from RecoEvent Store!", v_error, verbosity); return false; }
    trueVtxX = truevtx->GetPosition().X();  //already in cm
    trueVtxY = truevtx->GetPosition().Y();
    trueVtxZ = truevtx->GetPosition().Z();
    trueVtxTime = truevtx->GetTime();
    trueDirX = truevtx->GetDirection().X();
    trueDirY = truevtx->GetDirection().Y();
    trueDirZ = truevtx->GetDirection().Z();
    trueTrackDir = TVector3(trueDirX,trueDirY,trueDirZ).Unit();

    trueAngleRad = TMath::ACos(trueDirZ);  //calculated like this in other tools
    trueAngleDeg = trueAngleRad/(TMath::Pi()/180.);
    h_truevtx_angle->Fill(trueAngleDeg);

    RecoVertex *truestopvtx = 0;
    get_ok = m_data->Stores["RecoEvent"]->Get("TrueStopVertex", truestopvtx);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving TrueStopVertex from RecoEvent Store!", v_error, verbosity); return false; }
    trueStopVtxX = truestopvtx->GetPosition().X();  //already in cm
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

  //-- Get ev_id for matching
  std::stringstream ev_id;
  ev_id << "p" << partnumber << "_";
  if (isData) ev_id << evnum;
  else ev_id << mcevnum;
  std::cout << "MuonFitter Tool: Working on event " << ev_id.str() << std::endl;

  if (reco_mode)
  {
    //-- Skip events that weren't fitted to make processing faster
    std::map<std::string, std::vector<double>>::iterator it = m_tank_track_fits.find(ev_id.str());
    if (it == m_tank_track_fits.end()) return true;
  }

  // ------------------------------------------------------------
  // --- Get event info (both) ----------------------------------
  // ------------------------------------------------------------
  uint32_t trigword;
  get_ok = m_data->Stores.at("ANNIEEvent")->Get("TriggerWord",trigword);
  std::cout << " trigword: " << trigword << std::endl;


  // ------------------------------------------------------------
  // --- CUT: Select beam only events ---------------------------
  // ------------------------------------------------------------
  if (trigword != 5) return true;


  // ------------------------------------------------------------
  // --- Check for particles other than muon (MC ONLY) ----------
  // ------------------------------------------------------------
  bool hasPion = false;
  int n_rings = 0;
  if (!isData)
  {
    for (unsigned int mcp_i = 0; mcp_i < mcParticles->size(); mcp_i++)
    {
      MCParticle aparticle = mcParticles->at(mcp_i);
      //-- QA: Check MC evnum and particle PDG code
      //std::cout << " [debug] Ev " << mcevnum << ", particle pdg code: " << aparticle.GetPdgCode() << std::endl;
      if (std::find(particles_ring.begin(), particles_ring.end(), mcp_i) != particles_ring.end())
      {
        //-- Use number of rings to determine if there are particles other than muon
        ++n_rings;
      }
    }
    std::cout << " [debug] n_rings (my counter): " << n_rings << std::endl;
    if (n_rings > 1) hasPion = true;  //might include other particles
  }
  if (hasPion) std::cout << " [debug] has pion / other particle!" << std::endl;


  // ------------------------------------------------------------
  // --- CUT: Check for FMV hits --------------------------------
  // ------------------------------------------------------------
  bool hasVeto = false;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("TDCData", tdcdata);   //'TDCData' for both data and MC
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
  //-- Skip event if there is a veto hit
  if (hasVeto)
  {
    Log("MuonFitter Tool: Found FMV/Veto hit!", v_debug, verbosity);
    return true;
  }


  // ------------------------------------------------------------
  // --- CUT: Check for MRD tracks ------------------------------
  // ------------------------------------------------------------
  //get_ok = m_data->Stores["ANNIEEvent"]->Get("MRDTriggerType", mrdTriggerType);  //XXX:care about this? not defined yet
  get_ok = m_data->Stores["MRDTracks"]->Get("NumMrdTracks", numTracksInEv);
  get_ok = m_data->Stores["MRDTracks"]->Get("MRDTracks", mrdTracks);    //XXX:might need MC version

  if (!get_ok) { Log("MuonFitter Tool: Couldn't retrieve MRD tracks info. Did you run TimeClustering/FindMRDTracks first?", v_debug, verbosity); return false; }

  //-- Skip event if num tracks not equal to 1 in MRD
  if (numTracksInEv != 1)
  {
    Log("MuonFitter Tool: More than 1 reconstructed track found!", v_debug, verbosity);
    return true;
  }


  // ------------------------------------------------------------
  // --- Get MRD track params -----------------------------------
  // ------------------------------------------------------------
  double reco_mrd_track = 0.;
  double nlyrs_mrd_track = 0.;
  double conn_dots_mrd_track = 0.;
  bool isMrdStopped;
  bool isMrdPenetrating;
  bool isMrdSideExit;

  for(int track_i = 0; track_i < numTracksInEv; track_i++)
  {
    BoostStore* thisTrackAsBoostStore = &(mrdTracks->at(track_i));
    
    thisTrackAsBoostStore->Get("StartVertex", mrdStartVertex);  //m
    thisTrackAsBoostStore->Get("StopVertex", mrdStopVertex);    //m
    thisTrackAsBoostStore->Get("MrdEntryPoint", mrdEntryPoint); //m
    thisTrackAsBoostStore->Get("TankExitPoint", tankExitPoint); //m
    thisTrackAsBoostStore->Get("StartTime", mrdStartTime);
    thisTrackAsBoostStore->Get("TrackAngle", trackAngleRad);
    thisTrackAsBoostStore->Get("TrackAngleError", trackAngleError);
    thisTrackAsBoostStore->Get("PenetrationDepth", penetrationDepth);
    thisTrackAsBoostStore->Get("NumLayersHit", numLayersHit);
    thisTrackAsBoostStore->Get("LayersHit", LayersHit);
    thisTrackAsBoostStore->Get("PMTsHit", MrdPMTsHit);
    thisTrackAsBoostStore->Get("EnergyLoss", mrdEnergyLoss);
    thisTrackAsBoostStore->Get("EnergyLossError", mrdEnergyLossError);
    thisTrackAsBoostStore->Get("IsMrdStopped", isMrdStopped);
    thisTrackAsBoostStore->Get("IsMrdPenetrating", isMrdPenetrating);
    thisTrackAsBoostStore->Get("IsMrdSideExit", isMrdSideExit);

    //-- Calculate MRD track length using MRD track start and stop
    //-- NOTE: this was done in other tools
    reco_mrd_track = sqrt(pow((mrdStopVertex.X()-mrdStartVertex.X()),2)+pow(mrdStopVertex.Y()-mrdStartVertex.Y(),2)+pow(mrdStopVertex.Z()-mrdStartVertex.Z(),2));

    //-- Do some conversions
    reco_mrd_track = reco_mrd_track*100.;       //convert to cm
    penetrationDepth = penetrationDepth*100.;
    trackAngleDeg = trackAngleRad*180./TMath::Pi();

    //-- QA: Check MRD track variables
    std::cout << "  [debug] mrdStartVertex: " << mrdStartVertex.X() << "," << mrdStartVertex.Y() << "," << mrdStartVertex.Z() << std::endl;
    std::cout << "  [debug] mrdStopVertex: " << mrdStopVertex.X() << "," << mrdStopVertex.Y() << "," << mrdStopVertex.Z() << std::endl;
    std::cout << "  [debug] tankExitPoint: " << tankExitPoint.X() << "," << tankExitPoint.Y() << "," << tankExitPoint.Z() << std::endl;
    std::cout << "  [debug] mrdEntryPoint: " << mrdEntryPoint.X() << "," << mrdEntryPoint.Y() << "," << mrdEntryPoint.Z() << std::endl;
    std::cout << "  [debug] mrdStartTime: " << mrdStartTime << std::endl;
    std::cout << "  [debug] trackAngleRad: " << trackAngleRad << std::endl;
    std::cout << "  [debug] penetrationDepth: " << penetrationDepth << std::endl;
    std::cout << "  [debug] reco_mrd_track (ANNIE tools): " << reco_mrd_track << std::endl;
    std::cout << "  [debug] mrdEnergyLoss: " << mrdEnergyLoss << std::endl;
    std::cout << "  [debug] isMrdStopped: " << isMrdStopped << std::endl;
    std::cout << "  [debug] isMrdPenetrating: " << isMrdPenetrating << std::endl;
    std::cout << "  [debug] isMrdSideExit: " << isMrdSideExit << std::endl;
    std::cout << "  [debug] numLayersHit: " << numLayersHit << std::endl;
    std::cout << "  [debug] LayersHit.size(): " << LayersHit.size() << std::endl;
    std::cout << "  [debug] MrdPMTsHit.size(): " << MrdPMTsHit.size() << std::endl;

  } //-- Done retrieving MRD track params

  //TODO:Check whether MRD track falls into the fiducial volume

  //-- QA: Check the values of LayersHit vector
  for (int i = 0; i < LayersHit.size(); ++i)
  {
    std::cout << " [debug] LayersHit at " << i << ": " << LayersHit.at(i) << std::endl;
  } //-- returns numbers from 0-11 (layer number ordered)

  //-- Histogram number of layers traversed
  h_num_mrd_layers->Fill(LayersHit.size());


  //-- QA: Check the values of MrdPMTsHit vector
  for (int i = 0; i < MrdPMTsHit.size(); ++i)
  {
    std::cout << " [debug] MrdPMTsHit at " << i << ": " << MrdPMTsHit.at(i) << std::endl;
  } //-- returns unsorted list of MRD detkeys


  //-- Use num layers to determine amount of track length in iron
  //-- NOTE: this is an "effective" track length
  //-- (5cm)*(num layers)/cos(theta) << 5cm is thickness of iron slab
  nlyrs_mrd_track = 5.*LayersHit.size()/abs(TMath::Cos(trackAngleRad));

  //-- Compare MRD tracks reconstructed w/ ANNIE methods and num layers
  //h_mrd_nlyrs_reco->Fill(reco_mrd_track, nlyrs_mrd_track);


  // ------------------------------------------------------------
  // --- Use PCA to determine track angle in MRD ----------------
  // ------------------------------------------------------------
  //-- Code based on https://github.com/JaydipSingh/EnergyStudies/blob/main/trackAna_module.cc#L468
  //-- First define vectors to store xyz of each hit
  std::vector<double> fRec_SpacePoint_X;
  std::vector<double> fRec_SpacePoint_Y;
  std::vector<double> fRec_SpacePoint_Z;

  for (int i = 0; i < MrdPMTsHit.size(); ++i)
  {
    unsigned long detkey = (unsigned long)MrdPMTsHit.at(i);
    std::cout << " [debug] PCA detkey: " << detkey << std::endl;
    if (detkey < 26) continue;  //skip FMV hits
    Paddle *mrdpaddle = (Paddle *)geom->GetDetectorPaddle(detkey);
    Position position_MRD = mrdpaddle->GetOrigin();
    
    fRec_SpacePoint_X.push_back(mrd_center_x[detkey]);
    fRec_SpacePoint_Y.push_back(mrd_center_y[detkey]);
    fRec_SpacePoint_Z.push_back(mrd_center_z[detkey]);
  }
    
  //-- Put all the xyz info into one array
  double *xyz[3];
  TVector3 pos(0,0,0);
  TVector3 dir;
  xyz[0] = fRec_SpacePoint_X.data();
  xyz[1] = fRec_SpacePoint_Y.data();
  xyz[2] = fRec_SpacePoint_Z.data();
  size_t npts = fRec_SpacePoint_X.size();

  if (npts < 2)
  {
    //throw std::exception("tpcvechitfinder2_module.cc: too few TPCClusters to fit a line in linefit");
    std::cout << " MuonFitter Tool: PCA Method: Too few TPCClusters to fit a line in linefit" << std::endl;
  }

  TMatrixDSym covmat(3);  //covariance matrix (use symmetric version)

  // position is just the average of the coordinates
  double psum[3] = {0,0,0};
  for (size_t ipoint=0; ipoint<npts; ++ipoint)
  {
    for(size_t j=0; j<3; j++)
    {
      psum[j] += xyz[j][ipoint];
    }
  }
  for (size_t j=0; j<3; ++j)
  {
    psum[j] /= npts;
  }
  pos.SetXYZ(psum[0],psum[1],psum[2]);

  for(size_t i=0; i<3; ++i)
  {
    for (size_t j=0; j<= i; ++j)
    {
      double csum=0;
      for (size_t ipoint=0; ipoint<npts; ++ipoint)
      {
        csum += (xyz[i][ipoint] - psum[i]) * (xyz[j][ipoint] - psum[j]);
      }
      csum /= (npts-1);
      covmat[i][j] = csum;
      covmat[j][i] = csum;
    }
  }
  TVectorD eigenvalues(3);
  TMatrixD eigenvectors = covmat.EigenVectors(eigenvalues);

  double dirv[3] = {0,0,0};
  for (size_t i=0; i<3; ++i)
  {
    dirv[i] = eigenvectors[i][0];
  }
  dir.SetXYZ(dirv[0],dirv[1],dirv[2]);

  double pcaAngleRad = dir.Angle(TVector3(0,0,1));
  double pcaAngleDeg = pcaAngleRad*180./TMath::Pi();

  std::cout << " [debug] PCA dir XYZ: " << dir.X() << "," << dir.Y() << "," << dir.Z() << std::endl;
  std::cout << " [debug] PCA angle (rad): " << pcaAngleRad << std::endl;
  std::cout << " [debug] PCA angle (deg): " << pcaAngleDeg << std::endl;
  h_pca_angle->Fill(pcaAngleDeg);
  h_pca_true_angle->Fill(pcaAngleDeg-trueAngleDeg);
  h_pca_reco_angle->Fill(pcaAngleDeg-trackAngleDeg);

  //-- Update nlyrs_mrd_track so that it uses the PCA-reconstructed track angle
  double pca_mrd_track = 5.*LayersHit.size()/abs(TMath::Cos(pcaAngleRad));
  std::cout << " [debug] PCA nlyrs_mrd_track: " << pca_mrd_track << std::endl;
  if (use_pca) { nlyrs_mrd_track = pca_mrd_track; }

  //cout<<"Ave Space Point Position P  =(\t "<<pos[0]<<",\t"<<pos[1]<<",\t"<<pos[2]<<"\t)"<<endl;
  // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ PCA END ################################################
  // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$         Track Length Estimation from Space Point  $$$$$


  //-- Determine MRD track length by "connecting the dots"
  //-- First sort MRD PMTs IDs; MrdPMTsHit is type std::vector<int>
  std::sort(MrdPMTsHit.begin(), MrdPMTsHit.end());

  //-- QA: Check that MRD pmts have been sorted
  std::cout << " [debug] sorted MRD PMTs: ";
  for (int i = 0; i < MrdPMTsHit.size(); ++i) { std::cout << MrdPMTsHit.at(i) << ","; }
  std::cout << std::endl;

  //-- Now connect the dots
  for (int i = 1; i < MrdPMTsHit.size(); ++i)
  {
    unsigned long detkey1 = (unsigned long)MrdPMTsHit.at(i-1);
    unsigned long detkey2 = (unsigned long)MrdPMTsHit.at(i);

    //-- XXX:For some reason, there are FMV PMTs in these clusters
    if (detkey1 < 26) continue;

    double dist_btwn_lyrs = TMath::Sqrt(pow((mrd_center_x[detkey2]-mrd_center_x[detkey1]),2) + pow((mrd_center_y[detkey2]-mrd_center_y[detkey1]),2) + pow((mrd_center_z[detkey2]-mrd_center_z[detkey1]),2));

    conn_dots_mrd_track += dist_btwn_lyrs;

    //-- QA: Check calculations
    std::cout << " [debug] detkey 1,2: " << detkey1 << "," << detkey2 << std::endl;
    std::cout << " [debug] dist_btwn_lyrs: " << dist_btwn_lyrs << std::endl;
    std::cout << " [debug] conn_dots_mrd_track: " << conn_dots_mrd_track << std::endl;
  } //-- Done connecting the dots

  std::cout << " [debug] FINAL conn_dots_mrd_track: " << conn_dots_mrd_track << std::endl;
  std::cout << " [debug] FINAL nlyrs_mrd_track: " << nlyrs_mrd_track << std::endl;
  //-- Repurpose hist to compare nlyrs and conn dots
  h_mrd_nlyrs_reco->Fill(nlyrs_mrd_track, conn_dots_mrd_track);



  // ------------------------------------------------------------
  // --- CUT: Select min MRD track length / Track angle ---------
  // ------------------------------------------------------------
  //-- Select for tracks that end in MRD
  if (!isMrdStopped)
  {
    std::cout << " MuonFitter Tool: MUON DID NOT STOP IN MRD. Skipping..." << std::endl;
    return true;
  }

  if (isMrdSideExit)  //<< Not working
  {
    std::cout << " MuonFitter Tool: MUON EXITTED SIDE OF MRD. Skipping..." << std::endl;
    //return true;
  }

  //TODO: do something with angles
/*  if (reco_mrd_track < 50.)
  //if (reco_mrd_track < 20.)   //for nlyrs method
  {
    std::cout << " [debug] MRD track too short!" << std::endl;
    return true;
  }*/

  // MRD angle cut
/*  if (!isData)
  {
    double angle_diff = trackAngleRad*180./TMath::Pi()-trueAngleDeg;
    if (abs(angle_diff) > 5.) return true;  //select events where reco angle is close to true angle
  }*/
  //if (abs(trackAngleRad*180./TMath::Pi()) > 20.) return true;  //select events that are straight going

  // ------------------------------------------------------------
  // --- Now look at tank exit and MRD points -------------------
  // ------------------------------------------------------------
  double tankExitPointX = 100.*tankExitPoint.X()-tank_center_x;
  double tankExitPointY = 100.*tankExitPoint.Y()-tank_center_y;
  double tankExitPointZ = 100.*tankExitPoint.Z()-tank_center_z;
  TVector3 tankExit(tankExitPointX, tankExitPointY, tankExitPointZ);

  double mrdEntryPointX = 100.*mrdEntryPoint.X()-tank_center_x;
  double mrdEntryPointY = 100.*mrdEntryPoint.Y()-tank_center_y;
  double mrdEntryPointZ = 100.*mrdEntryPoint.Z()-tank_center_z;

  //-- Create vectors for MRD start/stop
  TVector3 mrdStart(100.*mrdStartVertex.X()-tank_center_x, 100.*mrdStartVertex.Y()-tank_center_y, 100.*mrdStartVertex.Z()-tank_center_z);   //cm
  TVector3 mrdStop(100.*mrdStopVertex.X()-tank_center_x, 100.*mrdStopVertex.Y()-tank_center_y, 100.*mrdStopVertex.Z()-tank_center_z);
  //-- Get track direction from MRD start/stop
  TVector3 mrdTrackDir = (mrdStop - mrdStart).Unit();
  //-- QA:Compare angle from mrdTrackDir to trackAngleRad
  std::cout << " [debug] mrdTrackDir, trackAngleRad: " << mrdTrackDir.Angle(TVector3(0,0,1)) << "," << trackAngleRad << std::endl;

  //-- ANNIE in 3D: Add mrdStart and mrdStop to 3D geo
  if (isData) sprintf(blockName,"mrdStart_%u",evnum);
  else sprintf(blockName,"mrdStart_%u",mcevnum);
  bBlock = ageom->MakeSphere(blockName, Iron, 0,1.5,0,180,0,360);
  bBlock->SetLineColor(8);
  EXPH->AddNodeOverlap(bBlock,69,new TGeoTranslation(mrdStart.X(), mrdStart.Y(), mrdStart.Z()));
  N++;

  if (isData) sprintf(blockName,"mrdStop_%u",evnum);
  else sprintf(blockName,"mrdStop_%u",mcevnum);
  bBlock = ageom->MakeSphere(blockName, Iron, 0,1.5,0,180,0,360);
  bBlock->SetLineColor(2);
  EXPH->AddNodeOverlap(bBlock,69,new TGeoTranslation(mrdStop.X(), mrdStop.Y(), mrdStop.Z()));
  N++;



  // ------------------------------------------------------------
  // --- Calculate several vertex candidates --------------------
  // ------------------------------------------------------------
  //-- This is kind of like an angle cut
  std::vector<TVector3> vtxCandidates;
  int c = 1;
  bool outsideTank = false;
  bool inSFV = false;   //spherical FV around tank center
  while (!outsideTank)
  {
    TVector3 v = mrdStart - c*10.*mrdTrackDir;  //cm
    if (c <= 5) { vtxCandidates.push_back(v); }   //get into tank first
    else if (c > 5)
    {
      //-- Check that vtx is contained in tank cylinder
      if ((pow(v.X(),2) + pow(v.Z(),2) >= pow(tank_radius*100.,2)) || (v.Y() > tank_height*100.) || (v.Y() < -tank_height*100.))
      {
        outsideTank = true;
      }
      else { vtxCandidates.push_back(v); }

      //-- Check if any vtx intersects 1-m sphere FV around center
      if (pow(v.X(),2) + pow(v.Y(),2) + pow(v.Z(),2) <= pow(tank_radius*100.-50.,2))
      {
        inSFV = true;
      }
    }
    c++;
  } //-- Done calculating vtx candidates

  //-- CUT: Make sure track originates in tank
  if (!inSFV) return true;

  std::cout << "MuonFitter Tool: Num of vtxCandidates: " << vtxCandidates.size() << std::endl;

  //-- Check distance btwn vertices to make sure it's 10cm
  for (int i = 1; i < vtxCandidates.size(); ++i)
  {
    TVector3 v0(vtxCandidates.at(i-1));
    TVector3 v1(vtxCandidates.at(i));
    double d = TMath::Sqrt(pow(v1.X()-v0.X(),2) + pow(v1.Y()-v0.Y(),2) + pow(v1.Z()-v0.Z(),2));
    //std::cout << "  [debug] distance btwn candidate vtx: " << d << std::endl;
  } //-- yes, it's 10cm


  // ------------------------------------------------------------
  // --- Check for short tracks ---------------------------------
  // ------------------------------------------------------------
  if (!isData)
  {
    //-- Flag if true vertex is downstream of tank
    if (trueVtxZ > 0.)
    {
      std::cout << " MuonFitter Tool: True vertex is in downstream half of tank. p" << partnumber << "_" << mcevnum << std::endl;
    }

    TVector3 trueTankTrack = tankExit - TVector3(trueVtxX, trueVtxY, trueVtxZ);   //<< before finding trueTrackLengthInWater
    //-- QA:Compare track legnths
    std::cout << " [debug] trueTankTrack, trueTrackLengthInWater: " << trueTankTrack.Mag() << "," << trueTrackLengthInWater << std::endl;

    //-- Flag if tank track is short
    if (trueTrackLengthInWater < 100.)
    {
      std::cout << " MuonFitter Tool: True tank track length is less than 100 cm! p" << partnumber << "_" << mcevnum << std::endl;
      //return true;
    }
  }


  // ------------------------------------------------------------
  // --- Load Tank Clusters -------------------------------------
  // ------------------------------------------------------------
  bool has_clusters = false;
  if (isData) { has_clusters = m_data->CStore.Get("ClusterMap", m_all_clusters); }
  else { has_clusters = m_data->CStore.Get("ClusterMapMC", m_all_clusters_MC); }
  if (!has_clusters)
  {
    std::cout << " MuonFitter Tool: No clusters found in CStore! Did you run ClusterFinder tool before this?" << std::endl;
    return false;
  }

  has_clusters = m_data->CStore.Get("ClusterMapDetkey", m_all_clusters_detkeys);  //same for data and mc
  if (!has_clusters)
  {
    std::cout << " MuonFitter Tool: No cluster detkeys found in CStore! Did you run ClusterFinder tool before this?" << std::endl;
    return false;
  }

  Log("MuonFitter Tool: Accessing pairs in all_clusters map", v_debug, verbosity);
  int cluster_num = 0;
  int cluster_size = 0;
  if (isData) { cluster_size = (int)m_all_clusters->size(); }
  else { cluster_size = (int)m_all_clusters_MC->size(); }
  std::cout << " [debug] cluster_size (num clusters in event): " << cluster_size << std::endl;


  // ------------------------------------------------------------
  // --- Find main cluster when more than one cluster -----------
  // ------------------------------------------------------------
  //-- Max charge and in [0..2000ns] time window
  //-- Code based on UserTools/EventDisplay/EventDisplay.cpp
  bool found_muon = false;
  double main_cluster_time = 0;
  double main_cluster_charge = 0;
  double main_cluster_hits = 0;

  //bool found_muon = false;
  double earliest_hittime = 0;
  std::vector<double> v_cluster_times;
  if (isData)
  {
    for (std::pair<double, std::vector<Hit>>&& apair : *m_all_clusters)
    {
      std::vector<Hit>& cluster_hits = apair.second;
      double temp_time = 0;
      double temp_charge = 0;
      int temp_hits = 0;

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
          h_clusterhit_x->Fill(x_pmt[detkey]);
          h_clusterhit_y->Fill(y_pmt[detkey]);
          h_clusterhit_z->Fill(z_pmt[detkey]);
          h_clusterhit_detkey->Fill(detkey);
        }
      }
      if (temp_hits > 0) temp_time /= temp_hits;  //mean time
      std::cout << " [debug] temp_time [ns]: " << temp_time << std::endl;
      if (temp_time > 2000. || temp_time < 0.) continue;  // not in time window
      if (temp_charge > main_cluster_charge)
      {
        found_muon = true;
        main_cluster_charge = temp_charge;  //sum of all charges
        main_cluster_time = apair.first;    //same as mean time above
        main_cluster_hits = cluster_hits.size();
      }
    }
  }
  else  //is MC
  {
    for (std::pair<double, std::vector<MCHit>>&& apair : *m_all_clusters_MC)
    {
      std::vector<MCHit>& cluster_hits_MC = apair.second;
      double temp_time = 0;
      double temp_charge = 0;
      int temp_hits = 0;
      std::vector<double> v_cluster_hits;

      for (int ihit = 0; ihit < cluster_hits_MC.size(); ihit++)
      {
        temp_hits++;
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
          h_clusterhit_x->Fill(x_pmt[detkey]);
          h_clusterhit_y->Fill(y_pmt[detkey]);
          h_clusterhit_z->Fill(z_pmt[detkey]);
          h_clusterhit_detkey->Fill(detkey);

          //-- QA: Check if coordinates are the same. If so, just use existing map
          Position det_pos = this_detector->GetDetectorPosition();
        }
      }

      //-- Sort cluster hit times to determine diff btwn earliest and latest times
      sort(v_cluster_hits.begin(), v_cluster_hits.end());
      std::cout << " [debug] all MC cluster hit times (sorted): ";
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
        main_cluster_hits = cluster_hits_MC.size();
        earliest_hittime = v_cluster_hits.at(0);
        std::cout << " [debug] earliest_hittime: " << earliest_hittime << std::endl;
      }
    }
  }
  //-- QA: Check cluster times, main cluster time, charge, num hits
  std::cout << " [debug] all cluster times: ";
  for (int ct = 0; ct < v_cluster_times.size(); ++ct)
  {
    std::cout << v_cluster_times.at(ct) << ",";
  }
  std::cout << std::endl;
  std::cout << " [debug] main_cluster_time [ns]: " << main_cluster_time << std::endl;
  std::cout << " [debug] main_cluster_charge [nC/pe]: " << main_cluster_charge << std::endl;
  std::cout << " [debug] main_cluster_hits [#]: " << main_cluster_hits << std::endl;


  // ------------------------------------------------------------
  // --- Check coincidence tank and MRD activity ----------------
  // ------------------------------------------------------------
  //-- Select events that are within +/- 50ns of PMTMRDOffset
  //-- in data: PMTMRDOffset ~745ns; in MC: PMTMRDOffset = 0
  //-- Code based on UserTools/EventSelector/EventSelector.cpp
  double tankmrd_tdiff = mrdStartTime - main_cluster_time;
  std::cout << " MuonFitter Tool: Time difference between tank and MRD activity [ns]: " << tankmrd_tdiff << std::endl;
  if ((tankmrd_tdiff > PMTMRDOffset-50) && (tankmrd_tdiff < PMTMRDOffset+50))
  {
    std::cout << "MuonFitter Tool: Main cluster is coincident with MRD cluster!" << std::endl;
  }
  else
  {
    found_muon = false;
    std::cout << "MuonFitter Tool: Main cluster is NOT coincident with MRD cluster! "
              << "Tank cluster and MRD cluster times are too different." << std::endl;
  }

  //-- Save num hits and charge of main/max cluster
  //-- NOTE: in MC, charge is given in PE
  h_total_pe_hits->Fill(main_cluster_hits, main_cluster_charge);


  // ------------------------------------------------------------
  // --- FOUND MUON CANDIDATE -----------------------------------
  // ------------------------------------------------------------
  double max_eta = 0.;  //for plotting

  if (!reco_mode)
  {
  if (found_muon)
  {
    std::cout << "MuonFitter Tool: Found muon candidate! Event: p" << partnumber << "_";
    if (isData) std::cout << evnum;
    else std::cout << mcevnum;
    std::cout << std::endl;
    drawEvent = true;

    //-- Save nhits, cluster_charge to txt files
    //-- TODO: get rid of this
    /*pehits_file << "p" << partnumber << "_";
    if (isData) pehits_file << evnum;
    else pehits_file << mcevnum;
    pehits_file << "," << main_cluster_hits << "," << main_cluster_charge << std::endl;*/

    //-- Save main cluster charge
    h_clusterPE->Fill(main_cluster_charge);


    // ------------------------------------------------------------
    // --- Load MAIN cluster hits and detkeys ---------------------
    // ------------------------------------------------------------
    std::vector<Hit> cluster_hits;
    std::vector<MCHit> cluster_hits_MC;
    if (isData) { cluster_hits = m_all_clusters->at(main_cluster_time); }
    else { cluster_hits_MC = m_all_clusters_MC->at(main_cluster_time); }
    std::vector<unsigned long> cluster_detkeys = m_all_clusters_detkeys->at(main_cluster_time);
    std::vector<double> x_hits, y_hits, z_hits;

    // ------------------------------------------------------------
    // --- Ring Imaging: Calculate ai for each PMT ----------------
    // ------------------------------------------------------------
    std::map<unsigned long, double> charge;   //-- used for making charge cuts at PMT level
    std::map<unsigned long, std::vector<double>> hittime;   //-- used for t0 calc
    std::map<int, std::vector<double>> m_PE_ai;
    std::map<int, std::vector<double>> m_fpmt_ai;

    charge.clear();
    hittime.clear();
    for (unsigned long detkey = 332; detkey < 464; ++detkey)
    {
      charge.emplace(detkey, 0.);
      hittime.emplace(detkey, std::vector<double>{-99.});
    }

    //-- First collect all the charges seen by the PMT in the cluster
    if (isData)
    {
      for (int i = 0; i < (int)cluster_hits.size(); ++i)
      {
        //for each cluster hit
        int chankey = cluster_hits.at(i).GetTubeId();
        std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(chankey);
        if (it != ChannelKeyToSPEMap.end())
        {
          Detector* this_detector = geom->ChannelToDetector(chankey);
          unsigned long detkey = this_detector->GetDetectorID();  //chankey same as detkey
          double hit_charge = cluster_hits.at(i).GetCharge();
          double hit_PE = hit_charge / ChannelKeyToSPEMap.at(chankey);
          double hit_time = cluster_hits.at(i).GetTime();
          hittime[detkey].push_back(hit_time);

          // keep track of total charge seen by each PMT
          charge[detkey] += hit_PE;

        } //end if ChannelKeyToSPEMap
      } //end cluster_hits loop
    }
    else  //is MC
    {
      for (int i = 0; i < (int)cluster_hits_MC.size(); ++i)
      {
        //for each cluster hit
        int chankey = cluster_hits_MC.at(i).GetTubeId();
        std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(chankey);
        if (it != ChannelKeyToSPEMap.end())
        {
          Detector* this_detector = geom->ChannelToDetector(chankey);
          unsigned long detkey = this_detector->GetDetectorID();  //chankey same as detkey
          double hit_PE = cluster_hits_MC.at(i).GetCharge();  //charge in MC is in PE
          double hit_charge = hit_PE * ChannelKeyToSPEMap.at(chankey);
          double hit_time = cluster_hits_MC.at(i).GetTime();
          hittime[detkey].push_back(hit_time);

          // keep track of total charge seen by each PMT
          charge[detkey] += hit_PE;

        } //end if ChannelKetToSPEMap
      } //end cluster_hits_MC loop
    }

    //-- Now go through charge, hittime maps
    int gi = 0;    //counter for TGraph
    int nhits = 0, nhits_incone = 0;
    double totalpe = 0, totalpe_incone = 0;
    std::vector<double> v_tzero;

    for (unsigned long detkey = 332; detkey < 464; ++detkey)
    {
      //-- Get the charge for this PMT
      double pmt_PE = charge[detkey];
      h_pmt_charge->Fill(pmt_PE);   //save charge to check what PMTs see
      
      //-- NOTE: there are PMTs w/ 0pe from charge map initialization
      if (pmt_PE == 0) continue;

      std::vector<double> v_hittimes = hittime[detkey];
      double pmt_t = 0;
      int nhits_pmt = v_hittimes.size();
      std::cout << " [debug] all hit times for detkey " << detkey << ": ";
      for (int ht = 0; ht < nhits_pmt; ++ht) { std::cout << v_hittimes.at(ht) << ","; }
      std::cout << std::endl;

      //-- Get earliest hit time bc Cherenkov radiation is usu earliest
      sort(v_hittimes.begin(), v_hittimes.end());
      pmt_t = v_hittimes.at(0);
      if (pmt_t == -99.) pmt_t = v_hittimes.at(1);
      std::cout << " MuonFitter Tool: Earliest time for this PMT: " << pmt_t << std::endl;
      

      //-- CUT: Skip PMTs in this cluster that don't see enough light
      if (pmt_PE < PMTQCut)
      {
        std::cout << " MuonFitter Tool: SKIPPING charge[" << detkey << "] < " << PMTQCut << "pe: " << charge[detkey] << std::endl;
        pmt_PE = 0.;    //set PMT charge to 0
        continue;       //skip PMT entirely
      }

      //-- Keep track of total hits, PE that meet cut in this cluster
      nhits += 1;
      totalpe += pmt_PE;

      //-- Get position and direction of PMT
      double hitX = x_pmt[detkey];
      double hitY = y_pmt[detkey];
      double hitZ = z_pmt[detkey];
      double dirPMTX = x_pmt_dir[detkey];
      double dirPMTY = y_pmt_dir[detkey];
      double dirPMTZ = z_pmt_dir[detkey];
      
      TVector3 pmt_dir = TVector3(dirPMTX, dirPMTY, dirPMTZ).Unit();

      //-- Get vector from tankExitPoint to PMT (Ri)
      TVector3 vec_Ri = TVector3(hitX,hitY,hitZ) - tankExit;
      double Ri = vec_Ri.Mag();
      h_tankexit_to_pmt->Fill(Ri);
      //-- QA: Check that vector was calculated correctly
      std::cout << " [debug] hitXYZ: " << hitX << "," << hitY << "," << hitZ << std::endl;
      std::cout << " [debug] tankExitXYZ: " << tankExit.X() << "," << tankExit.Y() << "," << tankExit.Z() << std::endl;
      std::cout << " [debug] vec_Ri: " << vec_Ri.X() << "," << vec_Ri.Y() << "," << vec_Ri.Z() << std::endl;
      std::cout << " [debug] vec_Ri direction: " << vec_Ri.Unit().X() << "," << vec_Ri.Unit().Y() << "," << vec_Ri.Unit().Z() << std::endl;

      //-- Get angle btwn Ri and muon direction (ai)
      double ang_alpha = vec_Ri.Angle(-mrdTrackDir);  //rad
      double ai = Ri * TMath::Sin(ang_alpha) / TMath::Tan(CHER_ANGLE_RAD) + Ri * TMath::Cos(ang_alpha);
      h_tanktrack_ai->Fill(ai);

      //-- Get the vector from vertex of ai to the PMT (bi)
      TVector3 vec_ai = tankExit - ai*mrdTrackDir;
      TVector3 vec_bi = TVector3(hitX,hitY,hitZ) - vec_ai;
      double bi = TMath::Sqrt(pow(hitX-vec_ai.X(),2) + pow(hitY-vec_ai.Y(),2) + pow(hitZ-vec_ai.Z(),2));
      //-- QA: Check that vector was calculated correctly
      std::cout << " [debug] vec_ai: " << vec_ai.X() << "," << vec_ai.Y() << "," << vec_ai.Z() << std::endl;
      std::cout << " [debug] ai: " << ai << std::endl;
      std::cout << " [debug] vec_bi: " << vec_bi.X() << "," << vec_bi.Y() << "," << vec_bi.Z() << std::endl;
      std::cout << " [debug] bi: " << bi << std::endl;


      //-- Find angle btwn true vertex and hit (MC ONLY)
      //-- not sure what is being done w this rn
      if (!isData)
      {
        TVector3 vec_truevtx_hit = TVector3(hitX,hitY,hitZ) - TVector3(trueVtxX,trueVtxY,trueVtxZ);
        double anglePmtTrueVtx = vec_truevtx_hit.Angle(mrdTrackDir)*180./TMath::Pi();
      }

      //-- Get angle btwn vector bi and PMT direction
      double psi = vec_bi.Angle(-pmt_dir);    //used only to get angle; don't use for magnitude
      if (psi > TMath::Pi()/2.) psi = TMath::Pi()-psi;
      //-- QA: Check that angle was calculated correctly
      //std::cout << " [debug] psi(+), psi(-): " << psi*180./TMath::Pi() << ", " << vec_bi.Angle(pmt_dir)*180./TMath::Pi() << std::endl;

      //-- Calculate the area of PMT & frustum (ring)
      double eff_area_pmt = 0.5 * m_pmt_area[detkey] * (1. + TMath::Cos(psi));  // effective area seen by photons from emission point
      h_eff_area_pmt->Fill(eff_area_pmt);
      double area_frustum = 2. * TMath::Pi() * step_size_ai * Ri;
      double f_pmt = eff_area_pmt / area_frustum;   //fraction of photocathode of total frustrum (ring) area
      h_fpmt->Fill(f_pmt);

      gr_effarea_detkey->SetPoint(gi, detkey, eff_area_pmt);
      gr_fpmt_detkey->SetPoint(gi, detkey, f_pmt);
      gr_effarea_ai->SetPoint(gi, ai, eff_area_pmt);
      gr_fpmt_ai->SetPoint(gi, ai, f_pmt);

      h_eta_ai->Fill(ai, pmt_PE / f_pmt);

      //-- For each distance (ai) from tank exit point bin, collect the
      //-- charges associated w/ this distance
      for (int a = 55; a < 500; a+=step_size_ai)
      {
        if (ai >= a-step_size_ai/2 && ai < a+step_size_ai/2)
        {
          m_PE_ai[a].push_back(pmt_PE);
          m_fpmt_ai[a].push_back(f_pmt);
        }
      }

      ++gi;

    } //-- Done going through charge[detkey] map

    // ------------------------------------------------------------
    // -- Make eta vs ai (tank track segment) graphs for fitting --
    // ------------------------------------------------------------
    //-- eta represents photon density. In principle, it should be 
    //-- a high, constant value inside Cherenkov disc, and should 
    //-- drop off once outside the disc
    int j = 0;
    double running_avg = 0;

    //-- Go thru ai and PMT area fraction maps and add up fractions
    //-- and charges seen at each track segment (ai)
    for (auto const &pair: m_fpmt_ai)
    {
      double total_PE_ai = 0;
      double total_fpmt_ai = 0;
      for (int e = 0; e < pair.second.size(); ++e)
      {
        total_fpmt_ai += pair.second.at(e);
        total_PE_ai += m_PE_ai[pair.first].at(e);
      }
      //-- Calculate photon density at this track segment (ai)
      double total_eta = total_PE_ai / total_fpmt_ai;

      //-- Used to indicate where highest photon density occurs
      if (total_eta > max_eta) max_eta = total_eta;     //for plotting
      
      //-- Plot (ai,eta) pair
      gr_eta_ai->SetPoint(j, pair.first, total_eta);

      //-- Save (ev_id,cluster_time,ai,eta) data to txt file for ML scripts
      pos_file << "p" << partnumber << "_";
      if (isData) pos_file << evnum;
      else pos_file << mcevnum;
      pos_file << "," << main_cluster_time << "," << pair.first << "," << total_eta << std::endl;

      //-- Get the overall avg photon density (eta)
      avg_eta += total_eta;
      
      if (j == 0) running_avg = avg_eta;
      else
      {
        std::cout << " [debug] avg_eta (before dividing j, so right now it is total eta): " << avg_eta << std::endl;
        running_avg = avg_eta / (j+1);
      }

      std::cout << " [debug] ai, running_avg: " << pair.first << ", " << running_avg << std::endl;
      gr_running_avg->SetPoint(j, pair.first, running_avg);

      //-- Keep track of avg eta to the left and right of the trueTrackLegnthInWater:
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

    //-- Keep track of avg eta to the left and right of the trueTrackLegnthInWater:
    if (!isData)
    {
      std::cout << " [debug] num pmts (left, right): " << num_left_eta << ", " << num_right_eta << std::endl;
      if (num_left_eta != 0) left_avg_eta /= num_left_eta;
      if (num_right_eta != 0) right_avg_eta /= num_right_eta;
      std::cout << " [debug] left avg: " << left_avg_eta << ", right avg: " << right_avg_eta << std::endl;
      h_lr_avg_eta->Fill(left_avg_eta);
      h_lr_avg_eta->Fill(right_avg_eta);
    }

    //-- Save truth info to txt files
    if (!isData)
    {
      truetrack_file << "p" << partnumber << "_" << mcevnum << ",";
      truetrack_file << trueTrackLengthInWater << "," << trueTrackLengthInMRD << "," << trueMuonEnergy << std::endl;
      //truetrack_file << trueTrackLengthInWater << "," << left_avg_eta << "," << right_avg_eta << std::endl;
    }
  } //-- End if found_muon
  }//-- End !reco_mode


  // ------------------------------------------------------------
  // --- RECO MODE: Fit tank track, muon energy, vertex ---------
  // ------------------------------------------------------------
  double fitted_tank_track = -999.;
  TVector3 fitted_vtx(-999,-999,-999);
  double reco_muon_ke = -999.;
  if (reco_mode)
  {
    double t0 = -999.;
    bool save_t0 = false;

    //-- Get ev_id for matching << this has been moved further up in script
    /*std::stringstream ev_id;
    ev_id << "p" << partnumber << "_";
    if (isData) ev_id << evnum;
    else ev_id << mcevnum;*/

    if (abs(trackAngleRad*180./TMath::Pi()) > 5.)
    {
      std::cout << " [debug] angle more than 5 degrees! Event: p" << partnumber << "_";
      if (isData) std::cout << evnum;
      else std::cout << mcevnum;
      std::cout << std::endl;
    }

    //-- Get the fitted tank track length for this event from file
    std::map<std::string, std::vector<double>>::iterator it = m_tank_track_fits.find(ev_id.str());
    if (it != m_tank_track_fits.end())
    {
      std::vector<double> v_fit_ctime = m_tank_track_fits.at(ev_id.str());
      double fit_cluster_time = (double)v_fit_ctime.at(0);
      fitted_tank_track = (double)v_fit_ctime.at(1);
      std::cout << " MuonFitter Tool: Found track, cluster time for " << ev_id.str() << ": " << fitted_tank_track << ", " << fit_cluster_time << endl;
      h_fitted_tank_track->Fill(fitted_tank_track);

      //-- Skip the bad fits
      //-- NOTE: This may no longer be needed if using RNN to fit; might need different cuts
      if (fitted_tank_track < 0) return false;

      if (!isData)
      {
        //-- Check diff btwn reco and true info
        //h_mrd_angle_diff->Fill(trackAngleDeg - trueAngleDeg);
        h_mrd_track_diff->Fill(reco_mrd_track - trueTrackLengthInMRD);
        //h_mrd_track_diff->Fill(conn_dots_mrd_track - trueTrackLengthInMRD);  //connect dots
        //h_mrd_track_diff_nlyrs->Fill(nlyrs_mrd_track - trueTrackLengthInMRD);
        h_mrd_track_diff_nlyrs->Fill(conn_dots_mrd_track - nlyrs_mrd_track);
        h_total_track_diff->Fill((reco_mrd_track+fitted_tank_track) - (trueTrackLengthInMRD+trueTrackLengthInWater));
        //h_total_track_diff->Fill((conn_dots_mrd_track+fitted_tank_track) - (trueTrackLengthInMRD+trueTrackLengthInWater));
        h_total_track_diff_nlyrs->Fill((5.*LayersHit.size()/abs(TMath::Cos(trackAngleRad))+fitted_tank_track) - (trueTrackLengthInMRD+trueTrackLengthInWater));
      }

      //-- Load cluster
      //-- NOTE: Must use main_cluster_time from above since precision is lost when saving to file
      std::vector<Hit> cluster_hits;
      std::vector<MCHit> cluster_hits_MC;
      if (isData)
      {
        cluster_hits = m_all_clusters->at(main_cluster_time);
        //main_cluster_nhits = cluster_hits.size();
        for (int ihit = 0; ihit < cluster_hits.size(); ihit++)
        {
          main_cluster_charge += cluster_hits.at(ihit).GetCharge();
        }
      }
      else
      {
        cluster_hits_MC = m_all_clusters_MC->at(main_cluster_time);
        //main_cluster_nhits = cluster_hits_MC.size();
        for (int ihit = 0; ihit < cluster_hits.size(); ihit++)
        {
          main_cluster_charge += cluster_hits_MC.at(ihit).GetCharge();
        }
      }
      std::vector<unsigned long> cluster_detkeys = m_all_clusters_detkeys->at(main_cluster_time);

      //-- Save main cluster charge for fitted events
      h_clusterPE_fit->Fill(main_cluster_charge);
      if (hasPion) h_clusterPE_fit_haspion->Fill(main_cluster_charge);  //look at charge for events with pions

      // ------------------------------------------------------------
      // --- Fit muon vertex from fitted tank track length ----------
      // ------------------------------------------------------------
      //-- Start from tank exit point and move backwards into tank
      fitted_vtx = tankExit - fitted_tank_track * mrdTrackDir;
      std::cout << " MuonFitter Tool: Fitted vtx xyz: " << fitted_vtx.X() << "," << fitted_vtx.Y() << "," << fitted_vtx.Z() << std::endl;
      h_vtxfit_x->Fill(fitted_vtx.X());
      h_vtxfit_y->Fill(fitted_vtx.Y());
      h_vtxfit_z->Fill(fitted_vtx.Z());
      h_topview_fit->Fill(fitted_vtx.Z(),fitted_vtx.X());
      h_sideview_fit->Fill(fitted_vtx.Z(),fitted_vtx.Y());

      //-- Check if hit falls inside cone of FITTED vertex
      TVector3 vtx2tankExit = tankExit - fitted_vtx;   //similar to ai vector
/*      TVector3 vtx2pmt = TVector3(hitX,hitY,hitZ) - fitted_vtx;

      double ang = vtx2pmt.Angle(mrdTrackDir);

      if (!isData && display_truth)
      {
        std::cout << " [debug] vtx2tankExit (before truth): " << vtx2tankExit.X() << "," << vtx2tankExit.Y() << "," << vtx2tankExit.Z() << std::endl;
        vtx2tankExit = TVector3(trueStopVtxX,trueStopVtxY,trueStopVtxZ) - TVector3(trueVtxX,trueVtxY,trueVtxZ);
        vtx2pmt = TVector3(hitX,hitY,hitZ) - TVector3(trueVtxX,trueVtxY,trueVtxZ);
        std::cout << " [debug] vtx2tankExit (after truth): " << vtx2tankExit.X() << "," << vtx2tankExit.Y() << "," << vtx2tankExit.Z() << std::endl;
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
*/

      // ------------------------------------------------------------
      // --- Reconstruct Muon Energy --------------------------------
      // ------------------------------------------------------------
      //-- Calculate initial MIP energy loss in TANK from tank track length
      double tank_track = fitted_tank_track;
      if (!isData && display_truth)
      {
        // Check energy calc with trueTrackLengthInWater
        tank_track = trueTrackLengthInWater;
      }
      double tank_dEdx = 2.000;    //1.992 MeV/cm
      double tank_edep = tank_track * tank_dEdx;
      double tank_mrd_dist = TMath::Sqrt(pow(tankExitPointX-mrdEntryPointX,2) + pow(tankExitPointY-mrdEntryPointY,2) + pow(tankExitPointZ-mrdEntryPointZ,2));
      double outtank_edep = tank_mrd_dist * 1.05;  //this should be dE/dx in air..; currently not used
      std::cout << " [debug] tank_mrd_dist [cm]: " << tank_mrd_dist << std::endl;
      std::cout << " [debug] outtank_edep: " << outtank_edep << std::endl;

      //-- Calculate initial MIP energy loss in MRD from MRD track length
      //-- TODO: Use another dEdx value for iron?
      double mrd_track = reco_mrd_track;
      if (use_nlyrs) mrd_track = nlyrs_mrd_track;
      if (use_conn_dots) mrd_track = conn_dots_mrd_track;

      double mrd_dEdx = 11.3;     //MeV/cm
      double mrd_edep = mrd_track * mrd_dEdx;
      if (use_eloss) mrd_edep = mrdEnergyLoss;  //use official ANNIE tool reco MRD energy loss as the starting point
      std::cout << " [debug] mrd_edep: " << mrd_edep << std::endl;  //check which MRD energy is used

      if (use_nlyrs)
      {
        std::cout << " [debug] Adjust MRD track length based on initial mrd_edep" << std::endl;
        if (mrd_edep >= 400.)
        {
          mrd_track = mrd_track / 0.95;
        }
        else if (mrd_edep < 400.)
        {
          mrd_track = mrd_track / 0.80;
        }

        // Update with new track length
        mrd_edep = mrd_track * mrd_dEdx;
        std::cout << " [debug] nlyrs upd mrd_track: " << mrd_track << std::endl;
        std::cout << " [debug] nlyrs upd mrd_edep: " << mrd_edep << std::endl;
      }

      //-- Get initial reco muon energy by adding tank and MRD energy depositions
      double reco_mu_e = tank_edep + mrd_edep;


      //-- Fine-tune reconstructed muon energy by finding best dE/dx
      //-- and updating the muon energy until it doesn't change by 2%
      double T_before = reco_mu_e;
      bool eres_ok = false;
      int n_tries = 0;
      std::cout << " [debug] initial tank_dEdx: " << tank_dEdx << std::endl;

/*      while (!eres_ok)
      {
        double new_dEdx = CalcTankdEdx(T_before);
        std::cout << " [debug] new_dEdx: " << new_dEdx << std::endl;
        double T_after = tank_track*new_dEdx + mrdEnergyLoss + tank_mrd_dist*new_dEdx;
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

        // ------------------------------------------------------------
        // --- Fine-tune reco muon energy by finding dEdx for changing 
        // --- muon energy after traveling some distance dx -----------
        // ------------------------------------------------------------
        std::cout << " MuonFitter Tool: Calculating muon energy iteratively..." << std::endl;
        double dx = 10;   //cm; step size to update dE/dx
        double input_Emu = reco_mu_e;   //input energy for dE/dx calculation

        //-- Use true energy minus muon's rest mass as initial input energy
        //if (!isData && display_truth) { input_Emu = trueMuonEnergy-105.66; }

        //-- Start with energy deposition in tank
        double sum_tank_Emu = 0.;   //keep track of new reco muon energy in tank
        double d_tanktrack = fitted_tank_track;   //this will be used to keep track (haha get it) of remaining tank track length
        if (!isData && display_truth)
        {
          //-- Check energy calc with true track length in water
          d_tanktrack = trueTrackLengthInWater;
        }
        std::cout << " [debug] Initial d_tanktrack: " << d_tanktrack << std::endl;

        while (d_tanktrack > 0)
        {
          if (d_tanktrack < dx)
          {
            std::cout << " [debug] remaining d_tanktrack: " << d_tanktrack << std::endl;

            double deltaE = CalcTankdEdx(input_Emu) * d_tanktrack;
            sum_tank_Emu += deltaE;
            //-- Calculate new input energy for MRD energy calc
            input_Emu -= deltaE;
            d_tanktrack -= d_tanktrack;
            std::cout << " [debug] input_Emu: " << input_Emu << std::endl; 
            std::cout << " [debug] sum_tank_Emu: " << sum_tank_Emu << std::endl; 
            std::cout << " [debug] d_tanktrack: " << d_tanktrack << std::endl;  //should be 0
            std::cout << " [debug] final dEdx: " << deltaE/dx << std::endl; 
          }
          else
          {
            double deltaE = CalcTankdEdx(input_Emu) * dx;
            sum_tank_Emu += deltaE;
            //-- Calculate new input energy and remaining tank track length
            input_Emu -= deltaE;
            d_tanktrack -= dx;
            std::cout << " [debug] input_Emu: " << input_Emu << std::endl; 
            std::cout << " [debug] sum_tank_Emu: " << sum_tank_Emu << std::endl; 
            std::cout << " [debug] d_tanktrack: " << d_tanktrack << std::endl;
            std::cout << " [debug] dEdx: " << deltaE/dx << std::endl; 
          }
        }
        std::cout << " [debug] input_Emu (for MRD): " << input_Emu << "," << reco_mu_e-sum_tank_Emu << std::endl; //should be equal
        std::cout << " [debug] Initial reconstructed Emu: " << reco_mu_e << std::endl;
        //-- Update reco muon energy with the new tank energy (same initial MRD energy deposition)
        reco_mu_e = sum_tank_Emu + mrd_edep;
        std::cout << " [debug] New reconstructed Emu (updated tank energy): " << reco_mu_e << std::endl;

        //-- Now update energy deposition in MRD
        std::cout << " [debug] Now calculating muon energy in MRD..." << std::endl;
        double d_mrdtrack = mrd_track;
        double sum_mrd_Emu = 0.;
        //-- Use trueTrackLengthInMRD to check energy calculation
        //if (!isData && display_truth) { d_mrdtrack = trueTrackLengthInMRD; }

        std::cout << " [debug] Initial d_mrdtrack: " << d_mrdtrack << std::endl;
        while (d_mrdtrack > 0)
        {
          if (input_Emu < 0) break;

          //-- TODO: Need condition on remaining track length
          if (input_Emu < 20.)
          {
            std::cout << " [debug] input_Emu is LESS THAN 20 MeV: " << input_Emu << std::endl;
            std::cout << " [debug] Remaining d_mrdtrack: " << d_mrdtrack << std::endl;

            /*if (use_nlyrs && (d_mrdtrack > dx))
            {
              // Up the remaining energy since we still have a lot of track left
              input_Emu = d_mrdtrack * 10.1;
            }*/

            //-- If input energy < 20 MeV, just add remaining energy
            sum_mrd_Emu += input_Emu;

            //-- If input energy < 20, add (remaining track)*(previous dE/dx)
            //sum_mrd_Emu += d_mrdtrack * mrd_dEdx;

            std::cout << " [debug] sum_mrd_Emu: " << sum_mrd_Emu << std::endl; 
            input_Emu -= input_Emu;

            //-- hist: How much track length do we have left?
            h_remainder_track_last20MeV->Fill(d_mrdtrack);
            //-- hist: What is the diff btwn reco_mu_e now and true_mu_e?
            h_true_reco_Ediff_last20MeV->Fill((sum_tank_Emu+sum_mrd_Emu)-(trueMuonEnergy-105.66));
            h_remainder_track_Ediff_last20MeV->Fill((sum_tank_Emu+sum_mrd_Emu)-(trueMuonEnergy-105.66), d_mrdtrack);
            break;
          }

          if (d_mrdtrack < dx)
          {
            std::cout << " [debug] remaining d_mrdtrack: " << d_mrdtrack << std::endl;
            double deltaE = CalcMRDdEdx(input_Emu) * d_mrdtrack;
            sum_mrd_Emu += deltaE;
            input_Emu -= deltaE;
            d_mrdtrack -= d_mrdtrack;
            mrd_dEdx = deltaE/dx; //save for remaining track length calc

            std::cout << " [debug] input_Emu: " << input_Emu << std::endl;
            std::cout << " [debug] sum_mrd_Emu: " << sum_mrd_Emu << std::endl; 
            std::cout << " [debug] d_mrdtrack: " << d_mrdtrack << std::endl;
            std::cout << " [debug] final dEdx: " << deltaE/dx << std::endl; 
          }
          else
          {
            double deltaE = CalcMRDdEdx(input_Emu) * dx;
            sum_mrd_Emu += deltaE;
            input_Emu -= deltaE;
            d_mrdtrack -= dx;
            mrd_dEdx = deltaE/dx;

            std::cout << " [debug] input_Emu: " << input_Emu << std::endl;
            std::cout << " [debug] sum_mrd_Emu: " << sum_mrd_Emu << std::endl; 
            std::cout << " [debug] d_mrdtrack: " << d_mrdtrack << std::endl;
            std::cout << " [debug] dEdx: " << deltaE/dx << std::endl; 
          }
        } //-- Done updating MRD energy
        std::cout << " [debug] Updated reconstructed MRD energy (mrd_edep) vs mrdEnergyLoss: " << mrd_edep << "," << mrdEnergyLoss << std::endl;
        //-- Update reconstructed muon energy with new MRD energy deposition
        reco_mu_e = sum_tank_Emu + sum_mrd_Emu;
        std::cout << " [debug] New reconstructed Emu (after updating MRD): " << reco_mu_e << std::endl;

        if (use_simple_ereco) 
        {
          reco_mu_e = tank_edep + mrd_edep;     //-- mrd_edep diff if use_eloss is 1
          //reco_mu_e = sum_tank_Emu + mrd_edep;   //-- TODO:later do this where tank gets dynamically updated but MRD doesn't
          std::cout << " [debug] reconstructed Emu (Simple E Reco): " << reco_mu_e << std::endl;
        }

        reco_muon_ke = reco_mu_e;   //-- store for CStore

        h_mrd_angle->Fill(trackAngleRad*180./TMath::Pi());    //-- Save the angle of reconstructed event
        h_mrd_eloss_diff->Fill(sum_mrd_Emu - mrdEnergyLoss);  //-- Diff btwn const dE/dx and ANNIE method

        //-- Calculated deltaR and other variables (MC ONLY -- requires truth info)
        if (!isData)
        {
          //-- Look at xyz diff btwn truth and fit/reco
          h_truefitdiff_x->Fill(fitted_vtx.X()-trueVtxX);
          h_truefitdiff_y->Fill(fitted_vtx.Y()-trueVtxY);
          h_truefitdiff_z->Fill(fitted_vtx.Z()-trueVtxZ);
          h_truefit_len_diff->Fill(TMath::Abs(fitted_tank_track-trueTrackLengthInWater));
          //-- Look at trueVtxXYZ
          h_truevtx_x->Fill(trueVtxX);
          h_truevtx_y->Fill(trueVtxY);
          h_truevtx_z->Fill(trueVtxZ);
          h_topview_truth->Fill(trueVtxZ,trueVtxX);
          h_sideview_truth->Fill(trueVtxZ,trueVtxY);

          //-- Calculate spatial resolution (deltaR)
          double deltaR = TMath::Sqrt(pow(fitted_vtx.X()-trueVtxX,2) + pow(fitted_vtx.Y()-trueVtxY,2) + pow(fitted_vtx.Z()-trueVtxZ,2));
          h_deltaR->Fill(deltaR);

          TVector3 true_vtx = TVector3(trueVtxX,trueVtxY,trueVtxZ);
          double transverse = deltaR * TMath::Sin((true_vtx-fitted_vtx).Angle(mrdTrackDir));
          h_transverse->Fill(TMath::Abs(transverse));

          double parallel = deltaR * TMath::Cos((true_vtx-fitted_vtx).Angle(mrdTrackDir));
          h_parallel->Fill(TMath::Abs(parallel));

          //-- Compare with true muon energy
          double true_mu_e = trueMuonEnergy - 105.66;   //minus rest mass
          reco_mu_e = reco_mu_e + ERECO_SHIFT;
          double ediff = reco_mu_e - true_mu_e;
          h_true_reco_E->Fill(true_mu_e, reco_mu_e);
          h_true_reco_Ediff->Fill(ediff);

          //-- Use timing to select good vertex fits
          /*if (h_tzero->GetStdDev() < 2.) 
          {
            std::cout << " [debug] h_t0 width: " << h_tzero->GetStdDev() << std::endl;
            h_true_reco_E->Fill(true_mu_e, reco_mu_e);
            h_true_reco_Ediff->Fill(ediff);
            //if (main_cluster_charge < 3500.) h_true_reco_Ediff->Fill(ediff);  //make a charge cut
          }*/

          if (abs(ediff) > 100)
          {
            std::cout << " !!HEY! THERE IS A HUGE ENERGY DIFFERENCE!! Event: " << ev_id.str() <<", Ediff: " << ediff << std::endl;
            std::cout << " >> TANK track lengths (reco; true): " << tank_track << "; " << trueTrackLengthInWater << std::endl;
            std::cout << " >> MRD track lengths (reco; true): " << mrd_track << "; " << trueTrackLengthInMRD << std::endl;
            std::cout << " >> Muon energy (reco; true): " << reco_mu_e << "; " << true_mu_e << std::endl;
            std::cout << " >> Muon vertex (reco; true): " << fitted_vtx.X() << "," << fitted_vtx.Y() << "," << fitted_vtx.Z() << "; " << trueVtxX << "," << trueVtxY << "," << trueVtxZ << std::endl;

            h_tank_track_diff_large->Fill(tank_track-trueTrackLengthInWater);
            h_mrd_track_diff_large->Fill(mrd_track-trueTrackLengthInMRD);
            h_deltaR_large->Fill(deltaR);

            //-- Check diff btwn reco and true track angles to check for bias in angle
            double angle_diff = trackAngleDeg-trueAngleDeg;
            h_mrd_angle_diff->Fill(angle_diff);

            //-- Save info to file
            //save_t0 = true;
            lg_ediff_file << ev_id.str() << ","
                          << ediff << ","
                          << hasPion << ","
                          << tank_track << ","
                          << trueTrackLengthInWater << ","
                          << mrd_track << ","
                          << trueTrackLengthInMRD << ","
                          << reco_mu_e << ","
                          << true_mu_e << std::endl;

            //-- Save main_cluster_charge of events with ediff > 100 and/or has pion
            h_clusterPE_lrg_ediff->Fill(main_cluster_charge);
            if (hasPion) h_clusterPE_lrg_ediff_haspion->Fill(main_cluster_charge);
          }
          else
          {
            //-- Save info for events with ediff < 100
            h_tank_track_diff_small->Fill(tank_track-trueTrackLengthInWater);
            h_mrd_track_diff_small->Fill(mrd_track-trueTrackLengthInMRD);
            h_deltaR_small->Fill(deltaR);
          }


          //-- Look at fractions of the tank and MRD track lengths out of entire track length
          double f_tanktrack = (tank_track + tank_mrd_dist) / (tank_track + tank_mrd_dist + reco_mrd_track);
          double f_mrdtrack = (reco_mrd_track) / (tank_track + tank_mrd_dist + reco_mrd_track);
          h_Ediff_frac_tank->Fill(f_tanktrack, ediff);
          h_Ediff_frac_mrd->Fill(f_mrdtrack, ediff);

        }
        //-- histogram: Compare energy deposition in MRD and tank
        h_tank_mrd_E->Fill(mrdEnergyLoss, tank_edep);
      }
      else std::cout << " MuonFitter Tool: womp womp. No fitted track found for: " << ev_id.str() << std::endl;


    // Find mean t0 and subtract to center hist around 0
    //double mean_t0 = 0.;
/*    for (int t = 0; t < v_tzero.size(); ++t) { mean_t0 += v_tzero.at(t); }
    if (v_tzero.size() != 0) mean_t0 /= v_tzero.size();
    std::cout << " [debug] mean_t0: " << mean_t0 << std::endl;

    for (int t = 0; t < v_tzero.size(); ++t) { h_tzero->Fill(v_tzero.at(t)-mean_t0); }

    h_uber_t0widths->Fill(h_tzero->GetStdDev());
    if (!isData && reco_mode && save_t0)
    {
      lg_ediff_file << "," << h_tzero->GetStdDev() << std::endl;
    }
*/
    //save_t0 = false;  //reset


    // find which pmt charges are inside/outside cone of true track direction
/*    if (!isData)
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
    }*/

  } //-- End reco_mode


  // ------------------------------------------------------------
  // --- Store variables to CStore for downstream tools ---------
  // ------------------------------------------------------------
  std::cout << " MuonFitter Tool: Setting FittedTrackLengthInWater to CStore: " << fitted_tank_track << std::endl;
  m_data->CStore.Set("FittedTrackLengthInWater", fitted_tank_track);  //could be -888 or -999 if no fit
  std::cout << " MuonFitter Tool: Setting FittedMuonVertex to CStore" << std::endl;
  Position fitted_muon_vtx(fitted_vtx.X(), fitted_vtx.Y(), fitted_vtx.Z());   //could be -888 or -999 if no fit
  m_data->CStore.Set("FittedMuonVertex", fitted_muon_vtx);
  std::cout << " MuonFitter Tool: Setting RecoMuonKE to CStore: " << reco_muon_ke << std::endl;
  m_data->CStore.Set("RecoMuonKE", reco_muon_ke);   //could be -888 or -999 if no fit


  // ------------------------------------------------------------
  // --- Draw event ---------------------------------------------
  // ------------------------------------------------------------
  root_outp->cd();

  //-- ANNIE in 3D
  //-- Code based on UserTools/EventDisplay/EventDisplay.cpp
  if (plot3d && drawEvent)
  {
    if (verbosity > v_debug) { std::cout << "MuonFitter Tool: Saving 3D plots..." << std::endl; }

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

  //-- Save graphs for fitting
  if (!reco_mode && found_muon)
  {
    //------------------------------------------------------------
    //-- Area/frustum correction ---------------------------------
    //------------------------------------------------------------

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


    //------------------------------------------------------------
    //-- Fraction of PMT overlap with ring/frustrum area ---------
    //------------------------------------------------------------

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


    //------------------------------------------------------------
    //-- Photon density (eta) vs tank track segment (ai) ---------
    //------------------------------------------------------------

    c_eta_ai->cd();
    std::stringstream ss_eta_title, ss_eta_ai_name, ss_truetrack;
    ss_eta_title << "#eta: #PE Divided By f (tot PMT charge > " << PMTQCut << ") Ev_";
    ss_eta_ai_name << "c_eta_ai";
    if (isData) ss_eta_ai_name << "_r" << runnumber;
    else ss_eta_ai_name << "_wcsim";
    ss_eta_ai_name << "_p" << partnumber << "_ev";
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
    c_eta_ai->Modified();
    c_eta_ai->Update();
    c_eta_ai->SetTitle(ss_eta_title.str().c_str());
    c_eta_ai->SetName(ss_eta_ai_name.str().c_str());

    TLegend *legend4 = new TLegend(0.55, 0.65, 0.89, 0.89);
    legend4->AddEntry(gr_eta_ai, "#eta", "lp");
    legend4->AddEntry(gr_running_avg, "running avg #eta", "lp");

    if (!isData)
    { //indicate distance btwn trueVtx and tankExitPoint
      /*TLine *lTrueTankTrack = new TLine(trueTrackLengthInWater, 0, trueTrackLengthInWater, max_eta+1);
      lTrueTankTrack->SetLineColor(4);
      lTrueTankTrack->SetLineWidth(2);
      if (display_truth)
      {
        lTrueTankTrack->Draw();
        legend4->AddEntry(lTrueTankTrack, "true tank track length", "l");
      }*/

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
    TLine *lAvgEta = new TLine(55, avg_eta, 500-step_size_ai/2, avg_eta);
    lAvgEta->SetLineColor(46);
    lAvgEta->SetLineWidth(2);
    lAvgEta->Draw();
    legend4->AddEntry(lAvgEta, "avg #eta", "l");
    legend4->Draw();
    c_eta_ai->Write();
    if (!isData && display_truth) ss_eta_ai_name << "_truth";
    ss_eta_ai_name << ".png";
    c_eta_ai->SaveAs(ss_eta_ai_name.str().c_str());
  } //-- End !reco_mode

    
  if (reco_mode)
  {
    c_h_tzero->cd();
    std::stringstream ss_t0_title, ss_t0_name;
    ss_t0_title << "t0" << " p" << partnumber << "_ev";
    ss_t0_name << "c_h_tzero";
    if (isData) ss_t0_name << "_r" << runnumber;
    else ss_t0_name << "_wcsim";
    ss_t0_name << "_p" << partnumber << "_ev";
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
    c_h_tzero->SaveAs(ss_t0_name.str().c_str());
  }

  return true;
}


bool MuonFitter::Finalise(){
  // Save output
  root_outp->cd();
  //h_expected_PE->Write();
  //h_phot_inc_angle->Write();
  //h_hit_angles->Write();
  //h_fitted_tank_track_len->Write();
  //h_closest_approach->Write();
  h_num_mrd_layers->Write();
  //h_clusterhit_x->Write();
  //h_clusterhit_y->Write();
  //h_clusterhit_z->Write();
  //h_clusterhit_detkey->Write();
  //h_clusterhit_timespread->Write();
  //h_clusterhit_time->Write();
  if (!isData)
  {
    //h_truevtx_z->Write();
    //h_truevtx_angle->Write();
    h_true_tanktrack_len->Write();
    h_truevtx_x->Write();
    h_truevtx_y->Write();
    h_truevtx_z->Write();
    h_topview_truth->Write();
    h_sideview_truth->Write();
    h_mrd_angle_diff->Write();
    h_mrd_angle->Write();
    h_clusterPE->Write();
    h_clusterPE_fit->Write();
    h_clusterPE_fit_haspion->Write();
    h_clusterPE_lrg_ediff->Write();
    h_clusterPE_lrg_ediff_haspion->Write();
    h_pca_angle->Write();
    h_pca_true_angle->Write();
    h_pca_reco_angle->Write();
    h_remainder_track_last20MeV->Write();
    h_true_reco_Ediff_last20MeV->Write();
    h_remainder_track_Ediff_last20MeV->Write();
  }
  h_tankexit_to_pmt->Write();
  h_tanktrack_ai->Write();
  h_eff_area_pmt->Write();
  h_fpmt->Write();
  //h_eta_ai->Write();
  //h_qincone_truevtx->Write();
  //h_qoutcone_truevtx->Write();
  //h_total_pe_hits->Write();
  h_pmt_charge->Write();
  //h_tdiff->Write();
  if (reco_mode)
  {
    h_uber_t0widths->Write();
    h_fitted_tank_track->Write();
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
      h_tank_track_diff_small->Write();
      h_tank_track_diff_large->Write();
      h_mrd_track_diff_small->Write();
      h_mrd_track_diff_large->Write();
      h_deltaR_small->Write();
      h_deltaR_large->Write();
      h_mrd_nlyrs_reco->Write();
      h_mrd_track_diff_nlyrs->Write();
      h_mrd_track_diff->Write();
      h_total_track_diff_nlyrs->Write();
      h_total_track_diff->Write();
    }
  }

  // Close 3D canvas
  if (plot3d)
  {
    canvas_3d->Clear();
    canvas_3d->Close();
    delete canvas_3d;
  }

  delete c_effarea_detkey;
  delete c_fpmt_detkey;
  delete c_effarea_ai;
  delete c_fpmt_ai;
  delete c_eta_ai;
  delete c_h_tzero;

  root_outp->Close();
  pos_file.close();
  //cpp_file.close();
  //pehits_file.close();
  if (!isData) truetrack_file.close();
  if (!isData) lg_ediff_file.close();

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
      double cluster_time = -9999.;
      event_id = fileLine.at(0);
      cluster_time = std::stod(fileLine.at(1));
      tank_track_fit = std::stod(fileLine.at(2));
      std::vector<double> v_fit_ctime{cluster_time, tank_track_fit};
      if (verbosity > 5) std::cout << "Loading tank track fit: " << event_id << ", " << cluster_time << ", " << tank_track_fit << std::endl;
      m_tank_track_fits.insert(std::pair<std::string,std::vector<double>>(event_id, v_fit_ctime));
    }
  }
  return;
}

double MuonFitter::CalcTankdEdx(double input_E)
{
  double gamma = (input_E+105.66)/105.66;
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

double MuonFitter::CalcMRDdEdx(double input_E)
{
  // Using iterative dE/dx method
  double gamma = (input_E+105.66)/105.66;
  double beta = TMath::Sqrt(gamma*gamma-1.)/gamma;
  double Tmax = (2.*0.511*beta*beta*gamma*gamma)/(1.+2.*0.511*gamma/105.66 + pow(0.511/105.66,2));
  double fe_ln = TMath::Log(2.*0.511*beta*beta*gamma*gamma*Tmax/(I_FE*I_FE));
  double fe_dEdx = 7.121*0.307*0.464*pow(1./beta,2)*(0.5*fe_ln-beta*beta-0.5);

  return fe_dEdx;
}
