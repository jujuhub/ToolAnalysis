//TODO:consider renaming tool to TankEnergyEstimator/EnergyEstimator?
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

  // Files
  std::string outfile;
  m_variables.Get("OutputFile", outfile);
  Log("MuonFitter Tool: Saving output into " + outfile, v_message, verbosity);
  root_outp = new TFile(outfile.c_str(), "RECREATE");

  // --------------------------------------------------
  // --- Initialize canvases, graphs, histograms ------
  // --------------------------------------------------
  c_vtx_charge = new TCanvas("c_vtx_charge", "Total Charge Seen at Each Vertex", 800, 600);
  c_vtx_detkey = new TCanvas("c_vtx_detkey", "Num of PMTs In/Out Cone at Each Vertex", 800, 600);
  c_charge_per_pmt = new TCanvas("c_charge_per_pmt", "Charge Per PMT at Each Vertex", 800, 600);
  c_effarea_detkey = new TCanvas("c_effarea_detkey", "Effective PMT Area vs detkey", 800, 600);
  c_fpmt_detkey = new TCanvas("c_fpmt_detkey", "Effective PMT Area as Fraction of Frustum Area", 800, 600);
  c_effarea_ai = new TCanvas("c_effarea_ai", "Effective PMT Area vs ai (tank track)", 800, 600);
  c_fpmt_ai = new TCanvas("c_fpmt_ai", "Effective PMT Area as Fraction of Frustum Area", 800, 600);
  c_eta_ai = new TCanvas("c_eta_ai", "Charge per cm (eta)", 800, 600);

  h_alpha = new TH1D("h_alpha", "Alpha Value Distribution", 10, 0., 10.);
  h_alpha->GetXaxis()->SetTitle("alpha_i (x1e3)");   //TODO:latex

  h_expected_PE = new TH1D("h_expected_PE", "Expected Num of PE", 100, 0., 100.);
  h_expected_PE->GetXaxis()->SetTitle("expected # of PE");

  h_phot_inc_angle = new TH1D("h_phot_inc_angle", "Photon Incident Angle", 180, 0., 180.);
  h_phot_inc_angle->GetXaxis()->SetTitle("incident angle [deg]");

  h_hit_angles = new TH1D("h_hit_angles", "Hit Angles wrt Vertex and MRD Track Dir", 180, 0., 180.);
  h_hit_angles->GetXaxis()->SetTitle("hit angle [deg]");

  h_fitted_tank_track_len = new TH1D("h_fitted_tank_track_len", "Length of Tank Tracks", 500, 0., 500.);
  h_fitted_tank_track_len->GetXaxis()->SetTitle("track length [cm]");

  h_closest_approach = new TH1D("h_closest_approach", "Distance of Closest Approach (Btwn Vtx & Tank Ctr)", 500, 0., 500.);
  h_closest_approach->GetXaxis()->SetTitle("distance to tank center [cm]");

  h_num_mrd_layers = new TH1D("h_num_mrd_layers", "Num MRD Layers Hit", 50, 0, 50.);
  h_num_mrd_layers->GetXaxis()->SetTitle("# layers");
  
  h_truevtx_z = new TH1D("h_truevtx_z", "Z coordinate of True Vertex", 300, -150., 150.);
  h_truevtx_z->GetXaxis()->SetTitle("Z [cm]");

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

  h_eta_ai = new TH2D("h_eta_ai", "#PE/f vs. Tank Track Length", 50, 0., 500., 500, 0., 5000.);
  h_eta_ai->GetXaxis()->SetTitle("ai [cm]");
  h_eta_ai->GetYaxis()->SetTitle("#PE/f");

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

  h_truevtx_recoexit_track = new TH1D("h_truevtx_recoexit_track", "Tank Track Length (true vertex to reco tank exit)", 350, 0, 350);
  h_truevtx_recoexit_track->GetXaxis()->SetTitle("track length [cm]");

  h_truevtx_trueexit_track = new TH1D("h_truevtx_trueexit_track", "Tank Track Length (true vertex to true tank exit)", 350, 0, 350);
  h_truevtx_trueexit_track->GetXaxis()->SetTitle("track length [cm]");

  h_pmt_charge = new TH1D("h_pmt_charge", "Total Charge Seen By PMT", 500, 0, 500);
  h_pmt_charge->GetXaxis()->SetTitle("charge [pe]");

  h_lr_avg_eta = new TH1D("h_lr_avg_eta", "Avg Eta to Left/Right of True Tank Track Length", 100, 0, 5000);
  h_lr_avg_eta->GetXaxis()->SetTitle("eta [PE/cm]");

  h_avg_eta = new TH1D("h_avg_eta", "Overall Average Eta", 100, 0, 5000);
  h_avg_eta->GetXaxis()->SetTitle("eta [PE/cm]");


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
  gr_effarea_detkey->GetYaxis()->SetTitle("effective PMT area [cm2]");

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
  gr_effarea_ai->GetXaxis()->SetTitle("ai [cm]");
  gr_effarea_ai->GetYaxis()->SetTitle("effective PMT area [cm2]");

  gr_fpmt_ai = new TGraph();
  gr_fpmt_ai->SetLineWidth(0);
  gr_fpmt_ai->SetMarkerStyle(27);
  gr_fpmt_ai->SetMarkerColor(30);
  gr_fpmt_ai->SetFillColor(0);
  gr_fpmt_ai->GetXaxis()->SetTitle("ai [cm]");
  gr_fpmt_ai->GetYaxis()->SetTitle("f");

  gr_eta_ai = new TGraph();
  gr_eta_ai->SetLineWidth(0);
  gr_eta_ai->SetMarkerStyle(25);
  gr_eta_ai->SetMarkerColor(46);
  gr_eta_ai->SetFillColor(0);
  gr_eta_ai->GetXaxis()->SetTitle("ai [cm]");
  gr_eta_ai->GetYaxis()->SetTitle("#PE/f (eta)");

  gr_running_avg = new TGraph();
  gr_running_avg->SetLineWidth(0);
  gr_running_avg->SetMarkerStyle(8);
  gr_running_avg->SetMarkerColor(38);
  gr_running_avg->SetFillColor(0);
  gr_running_avg->GetXaxis()->SetTitle("ai [cm]");
  gr_running_avg->GetYaxis()->SetTitle("#PE/f (eta)");

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
  //m_data->CStore.Get("EventDisplay", canvas_ev_display);
  //m_data->CStore.Get("EventDisplays", m_evdisplays);

  // --------------------------------------------------
  // --- Get Geometry ---------------------------------
  // --------------------------------------------------
  // --- Retrieve Store variables ---------------------
  // --------------------------------------------------
  m_data->CStore.Get("ChannelNumToTankPMTSPEChargeMap", ChannelKeyToSPEMap);  //same for data and mc
  //m_data->CStore.Get("EventDisplay", canvas_ev_display);
  //m_data->CStore.Get("EventDisplays", m_evdisplays);

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

  // Save start & stop vertices to file
  std::string pos_fname = "posFile.txt";
  //pos_file.open(pos_fname.c_str());
  //pos_file << "##evnum,startX,startY,startZ,stopX,stopY,stopZ" << std::endl;

  std::string cpp_fname = "charge_per_pmt.txt";
  //cpp_file.open(cpp_fname.c_str());
  //cpp_file << "##evnum,nVtx,Qin,Qout,avgSumIn,avgSumOut,avgQin,avgQout" << std::endl;

  pehits_file.open("tot_pe_hits.txt", std::ostream::app);
  if (!pehits_file)
  {
    std::cout << " [debug] tot_pe_hits.txt does not exist! Creating now..." << std::endl;
    pehits_file.open("tot_pe_hits.txt");
    pehits_file << "##partfile,max_cluster_hits,max_cluster_charge" << std::endl;
  }

  lravg_file.open("left_right_avg.txt", std::ostream::app);
  if (!lravg_file)
  {
    std::cout << " [debug] left_right_avg.txt does not exist! Creating now..." << std::endl;
    lravg_file.open("left_right_avg.txt");
    lravg_file << "##partfile,left_avg,right_avg" << std::endl;
  }

/*  lravg_file.open("true_track_length.txt", std::ostream::app);
  if (!lravg_file)
  {
    std::cout << " [debug] left_right_avg.txt does not exist! Creating now..." << std::endl;
    lravg_file.open("true_track_length.txt");
    lravg_file << "##partfile,truetracklength" << std::endl;
  }*/

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

  avg_eta = 0;
  left_avg_eta = 0.;
  right_avg_eta = 0.;
  num_left_eta = 0;
  num_right_eta = 0;


  // --------------------------------------------------
  // --- Get event info (data) ------------------------
  // --------------------------------------------------
  int get_ok;
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

    double trueAngleRad = TMath::ACos(trueDirZ);
    trueAngle = trueAngleRad/(TMath::Pi()/180.);
    h_truevtx_angle->Fill(trueAngle);

    RecoVertex *truestopvtx = 0;
    get_ok = m_data->Stores["RecoEvent"]->Get("TrueStopVertex", truestopvtx);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving TrueStopVertex from RecoEvent Store!", v_error, verbosity); return false; }
    trueStopVtxX = truestopvtx->GetPosition().X();
    trueStopVtxY = truestopvtx->GetPosition().Y();
    trueStopVtxZ = truestopvtx->GetPosition().Z();

    get_ok = m_data->Stores["RecoEvent"]->Get("NRings", nrings);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving NRings, true from RecoEvent!", v_error, verbosity); }
    get_ok = m_data->Stores["RecoEvent"]->Get("IndexParticlesRing", particles_ring);
    if (not get_ok) { Log("MuonFitter Tool: Error retrieving IndexParticlesRing, true from RecoEvent!", v_error, verbosity); }
  }


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
  if (hasPion) return true;


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
  double tracklength = 0.;
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
    tracklength = sqrt(pow((mrdStopVertex.X()-mrdStartVertex.X()),2)+pow(mrdStopVertex.Y()-mrdStartVertex.Y(),2)+pow(mrdStopVertex.Z()-mrdStartVertex.Z(),2));

    std::cout << "  [debug] mrdStartVertex: " << mrdStartVertex.X() << "," << mrdStartVertex.Y() << "," << mrdStartVertex.Z() << std::endl;
    std::cout << "  [debug] mrdStopVertex: " << mrdStopVertex.X() << "," << mrdStopVertex.Y() << "," << mrdStopVertex.Z() << std::endl;
    std::cout << "  [debug] tankExitPoint: " << tankExitPoint.X() << "," << tankExitPoint.Y() << "," << tankExitPoint.Z() << std::endl;
    std::cout << "  [debug] mrdEntryPoint: " << mrdEntryPoint.X() << "," << mrdEntryPoint.Y() << "," << mrdEntryPoint.Z() << std::endl;
    std::cout << "  [debug] trackAngle: " << trackAngle << std::endl;
    std::cout << "  [debug] numLayersHit: " << numLayersHit << std::endl;
    std::cout << "  [debug] penetrationDepth: " << penetrationDepth << std::endl;
    std::cout << "  [debug] tracklength: " << tracklength << std::endl;
    std::cout << "  [debug] energyLoss: " << energyLoss << std::endl;
    std::cout << "  [debug] mrdStartTime: " << mrdStartTime << std::endl;

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

      //check if any vtx intersects with FV
      if (pow(v.X(),2) + pow(v.Y(),2) + pow(v.Z(),2) <= pow(tank_radius*100.-50.,2))
      {
        std::cout << " [debug] (distance from tank center)^2: " << pow(v.X(),2) + pow(v.Y(),2) + pow(v.Z(),2) << std::endl;
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
  if (!inFV) { return true; }

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
    h_truevtx_z->Fill(trueVtxZ);
  }


  // --------------------------------------------------
  // --- Check short tracks ---------------------------
  // --------------------------------------------------
  // get distance btwn trueVtx and reco tankExitPoint
  double trueRecoTankTrackLength = 0.;
  if (!isData)
  {
    trueRecoTankTrackLength = TMath::Sqrt(pow(tankExitPointX-trueVtxX,2) + pow(tankExitPointY-trueVtxY,2) + pow(tankExitPointZ-trueVtxZ,2));
    std::cout << " [debug] trueRecoTankTrackLength: " << trueRecoTankTrackLength << std::endl;

    // skip if true vertex is downstream of tank
    if (trueVtxZ > 0.)
    {
      std::cout << " SKIPPING event. True vertex is in downstream half of tank. " << mcevnum << std::endl;
      //return true;
    }

    TVector3 trueTankTrack = TVector3(tankExitPointX, tankExitPointY, tankExitPointZ) - TVector3(trueVtxX, trueVtxY, trueVtxZ);
    // skip if tank track is short
    if (trueRecoTankTrackLength < 100.)
    {
      std::cout << " Track less than 100 cm! SKIPPING! " << mcevnum << std::endl;
      return true;
    }
    h_truevtx_recoexit_track->Fill(trueRecoTankTrackLength);

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
  double max_cluster = 0;
  double max_cluster_charge = 0;
  double max_cluster_hits = 0;

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
      if (temp_charge > max_cluster_charge)
      {
        found_muon = true;
        max_cluster_charge = temp_charge;
        max_cluster = apair.first;    //actually the cluster time (same as mean time above)
        max_cluster_hits = cluster_hits.size();
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
        std::cout << " [debug] hit time: " << cluster_hits_MC.at(ihit).GetTime() << std::endl;
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
      std::cout << " [debug] all cluster hit times: ";
      for (int ct = 0; ct < v_cluster_hits.size(); ++ct)
      {
        std::cout << v_cluster_hits.at(ct) << ",";
      }
      std::cout << std::endl;
      h_clusterhit_timespread->Fill(v_cluster_hits.at(v_cluster_hits.size()-1)-v_cluster_hits.at(0));
      if (temp_hits > 0) temp_time /= temp_hits;  //mean time
      v_cluster_times.push_back(temp_time);
      if (temp_time > 2000. || temp_time < 0.) continue;  // not in time window
      if (temp_charge > max_cluster_charge)
      {
        found_muon = true;
        max_cluster_charge = temp_charge;
        max_cluster = apair.first;    //actually the cluster time (same as mean time above)
        //max_cluster_hits = cluster_hits_MC.size();
        earliest_hittime = v_cluster_hits.at(0);
        std::cout << " [debug] earliest_hittime: " << earliest_hittime << std::endl;

        max_cluster_hits = 0;
        for (unsigned long detkey = 332; detkey < 464; detkey++)
        {
          if (m_cluster_charge[detkey] > 0.)
          {
            max_cluster_hits++;
            h_pmt_charge->Fill(m_cluster_charge[detkey]);
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
  std::cout << " [debug] max_cluster (chosen cluster time): " << max_cluster << std::endl;
  std::cout << " [debug] max_cluster_charge: " << max_cluster_charge << std::endl;
  std::cout << " [debug] max_cluster_hits: " << max_cluster_hits << std::endl;

  // --------------------------------------------------
  // --- Check coincidence with MRD -------------------
  // --------------------------------------------------
  double tankmrd_tdiff = mrdStartTime - max_cluster;
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
  h_total_pe_hits->Fill(max_cluster_hits, max_cluster_charge);

  // save to txt files
  pehits_file << "p" << partnumber << "_";
  if (isData) pehits_file << evnum;
  else pehits_file << mcevnum;
  pehits_file << "," << max_cluster_hits << "," << max_cluster_charge << std::endl;


  double max_eta = 0.;

  if (found_muon)
  {
    std::cout << "MuonFitter Tool: Found muon candidate!" << std::endl;

    // Load cluster hits & detkeys
    std::vector<Hit> cluster_hits;
    std::vector<MCHit> cluster_hits_MC;
    if (isData) { cluster_hits = m_all_clusters->at(max_cluster); }
    else { cluster_hits_MC = m_all_clusters_MC->at(max_cluster); }
    std::vector<unsigned long> cluster_detkeys = m_all_clusters_detkeys->at(max_cluster);
    std::vector<double> x_hits, y_hits, z_hits;
    double cluster_charge = 0.;   //TODO:make cluster or individual PMT charge cut?

    // METHOD: CALCULATE ai FOR EACH PMT
    std::map<unsigned long, double> charge;
    std::map<int, std::vector<double>> m_PE_ai;
    std::map<int, std::vector<double>> m_fpmt_ai;

    if (isData)
    {
      charge.clear();
      for (unsigned long detkey = 332; detkey < 464; ++detkey) { charge.emplace(detkey, 0.); }

      std::cout << " Need to fill in code for data case (ai method)" << std::endl;
    }
    else  //is MC
    {
      charge.clear();
      for (unsigned long detkey = 332; detkey < 464; ++detkey) { charge.emplace(detkey, 0.); }

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
          x_hits.push_back(hitX);
          y_hits.push_back(hitY);
          z_hits.push_back(hitZ);
          //Direction dirPMT = this_detector->GetDetectorDirection();
          //TVector3 pmt_dir = TVector3(dirPMT.X()*100., dirPMT.Y()*100., dirPMT.Z()*100.).Unit();

          // keep track of total charge seen by each PMT
          charge[detkey] += hit_PE;
        } //end if ChannelKetToSPEMap; TEST PE CUT
      } //end cluster_hits_MC loop; TEST PE CUT


      int i = 0;    //index for TGraph
      for (unsigned long detkey = 332; detkey < 464; ++detkey)
      {
          //PE CUT
          //NOTE: there will be PMTs w/ 0 charge from initialization of charge map
          if (charge[detkey] == 0) continue;

          double hit_PE = charge[detkey];
          if (charge[detkey] < PMTQCut)
          {
            std::cout << " [debug] SKIPPING charge[" << detkey << "] < " << PMTQCut << "pe: " << charge[detkey] << std::endl;
            //std::cout << " [debug] JK setting charge to 0" << std::endl;
            hit_PE = 0.;
            continue;
          }

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
          std::cout << " [debug] ang_alpha, detkey: " << ang_alpha*180./TMath::Pi() << ", " << detkey << std::endl;
          std::cout << " [debug] tankExitPoint: " << tankExitPointX << "," << tankExitPointY << "," << tankExitPointZ << std::endl;

          double ai = Ri * TMath::Sin(ang_alpha) / TMath::Tan(CHER_ANGLE_RAD) + Ri * TMath::Cos(ang_alpha);
          h_tanktrack_ai->Fill(ai);

          // get the vector from the vertex of ai to the PMT (bi)
          std::cout << " [debug] hitXYZ for detkey " << detkey << ": " << hitX << "," << hitY << "," << hitZ << std::endl;
          TVector3 vec_ai = TVector3(tankExitPointX, tankExitPointY, tankExitPointZ) - ai*mrdTrackDir;  //this is along the track
          std::cout << " [debug] vec_ai: " << vec_ai.X() << "," << vec_ai.Y() << "," << vec_ai.Z() << std::endl;
          TVector3 vec_bi = TVector3(hitX,hitY,hitZ) - vec_ai;
          std::cout << " [debug] vec_bi: " << vec_bi.X() << "," << vec_bi.Y() << "," << vec_bi.Z() << std::endl;

          //std::cout << " [debug] angle check: " << vec_bi.Angle(mrdTrackDir)*180./TMath::Pi() << std::endl;

          // find angle btwn true vertex and hit
          TVector3 vec_truevtx_hit = TVector3(hitX,hitY,hitZ) - TVector3(trueVtxX,trueVtxY,trueVtxZ);
          double anglePmtTrueVtx = vec_truevtx_hit.Angle(mrdTrackDir)*180./TMath::Pi();
          //if (anglePmtTrueVtx < CHER_ANGLE_DEG+insideAngle) h_qincone_truevtx->Fill(it_PE);
          //if (anglePmtTrueVtx > CHER_ANGLE_DEG+outsideAngle) h_qoutcone_truevtx->Fill(hit_PE);

          // loop through all hits for this ai and determine whether in/out cone
/*          double qInCone = 0., qOutCone = 0.;
          for (int j = 0; j < (int)cluster_hits_MC.size(); ++j)
          {
            //for each hit
            int chankey = cluster_hits_MC.at(j).GetTubeId();
            std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(chankey);
            if (it != ChannelKeyToSPEMap.end())
            {
              Detector* this_detector = geom->ChannelToDetector(chankey);
              unsigned long detkey = this_detector->GetDetectorID();  //chankey same as detkey
              Position det_pos = this_detector->GetDetectorPosition();
              double hit_PE = cluster_hits_MC.at(j).GetCharge();  //in PE
              //double hit_charge = hit_PE * ChannelKeyToSPEMap.at(chankey);
              double hit_time = cluster_hits_MC.at(j).GetTime();
              double hitX = 100.*(det_pos.X()-tank_center_x);   //[cm]
              double hitY = 100.*(det_pos.Y()-tank_center_y);   //[cm]
              double hitZ = 100.*(det_pos.Z()-tank_center_z);   //[cm]

              TVector3 vec_ai_pmt = TVector3(hitX,hitY,hitZ) - vec_ai;
              double anglePmtAi = vec_ai_pmt.Angle(mrdTrackDir)*180./TMath::Pi();
              //std::cout << " [debug] anglePmtAi: " << anglePmtAi << std::endl;
              if (anglePmtAi < CHER_ANGLE_DEG-5.) qInCone += hit_PE;
              if (anglePmtAi > CHER_ANGLE_DEG+5.) qOutCone += hit_PE;
            }
          }
          gr_qincone_ai->SetPoint(i, ai, qInCone);
          gr_qoutcone_ai->SetPoint(i, ai, qOutCone); */

          // get angle btwn vector bi and pmt dir
          double psi = vec_bi.Angle(-pmt_dir);
          if (psi > TMath::Pi()/2.) psi = TMath::Pi()-psi;
          //std::cout << " [debug] psi(+), psi(-): " << psi*180./TMath::Pi() << ", " << vec_bi.Angle(pmt_dir)*180./TMath::Pi() << std::endl;

          // calculate the area of PMT & frustum
          //double area_pmt = 0.5*TMath::Pi()*pow(4.*2.54, 2)*(1. + TMath::Cos(psi));
          double eff_area_pmt = 0.5*m_pmt_area[detkey]*(1.+TMath::Cos(psi));
          h_eff_area_pmt->Fill(eff_area_pmt);
          double area_frustum = 2.*TMath::Pi()*(deltaL)*Ri;
          double f_pmt = eff_area_pmt / area_frustum;
          h_fpmt->Fill(f_pmt);

          gr_effarea_detkey->SetPoint(i, detkey, eff_area_pmt);
          gr_fpmt_detkey->SetPoint(i, detkey, f_pmt);
          gr_effarea_ai->SetPoint(i, ai, eff_area_pmt);
          gr_fpmt_ai->SetPoint(i, ai, f_pmt);
          //gr_eta_ai->SetPoint(i, ai, hit_PE / f_pmt);

          h_eta_ai->Fill(ai, hit_PE / f_pmt);

          for (int i = 55; i < 500; i+=deltaL)
          {
            if (ai >= i-deltaL/2 && ai < i+deltaL/2)
            {
              m_PE_ai[i].push_back(hit_PE);
              m_fpmt_ai[i].push_back(f_pmt);
            }
          }


          //XXX: checking PMT direction by drawing sphere in front of PMT
/*          TVector3 Ri_dir = TVector3(hitX, hitY, hitZ) + 10.*vec_Ri.Unit();
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

          ++i;
        } //end charge[detkey] loop; TEST PE CUT

    //    } //end if ChannelKeytoSPEMap; COMMENTED OUT FOR TEST PE CUT
    //  } //end for loop thru cluster hits; COMMENTED OUT FOR TEST PE CUT
    } //end is MC

    // find which pmt charges are inside/outside cone of true track direction
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
      if (total_eta > max_eta) max_eta = total_eta;
      
      gr_eta_ai->SetPoint(j, pair.first, total_eta);
      avg_eta += total_eta;
      
      if (j == 0) running_avg = avg_eta;
      else
      {
        running_avg = avg_eta / j;

        if (running_avg < EtaThreshold && !found_vtx)
        {
          // first instance the avg dips below threshold
          std::cout << " [debug] BEST VERTEX (ai): " << pair.first << std::endl;
          found_vtx = true;
          bestFitAi = pair.first;
        }
      }

      std::cout << " [debug] ai, running_avg: " << pair.first << ", " << running_avg << std::endl;
      std::cout << " [debug] EtaThreshold: " << EtaThreshold << std::endl;
      gr_running_avg->SetPoint(j, pair.first, running_avg);

      if (pair.first < trueRecoTankTrackLength)   //left of true track length
      {
        std::cout << " [debug] total_eta (left): " << total_eta << std::endl;
        num_left_eta++;
        left_avg_eta += total_eta;
      }
      if (pair.first > trueRecoTankTrackLength)   //right of true track length
      {
        std::cout << " [debug] total_eta (right): " << total_eta << std::endl;
        num_right_eta++;
        right_avg_eta += total_eta;
      }
      j+=1;
    }
    avg_eta /= j;
    std::cout << " [debug] num pmts (left, right): " << num_left_eta << ", " << num_right_eta << std::endl;
    if (num_left_eta != 0) left_avg_eta /= num_left_eta;
    if (num_right_eta != 0) right_avg_eta /= num_right_eta;
    std::cout << " [debug] left avg: " << left_avg_eta << ", right avg: " << right_avg_eta << std::endl;
    std::cout << " [debug] avg eta: " << avg_eta << ", j: " << j << std::endl;
    h_avg_eta->Fill(avg_eta);
    h_lr_avg_eta->Fill(left_avg_eta);
    h_lr_avg_eta->Fill(right_avg_eta);

    // running avg of entire graph
    // this might just be smoothing...
    for (int e = 1; e < j-2; ++e)
    {
      double three_pt_avg = (gr_eta_ai->GetY()[e-1] + gr_eta_ai->GetY()[e] + gr_eta_ai->GetY()[e+1]) / 3.;
      //gr_running_avg->SetPoint(e-1, gr_eta_ai->GetX(e), three_pt_avg);
    }

    // save to txt files
    lravg_file << "p" << partnumber << "_";
    if (isData) lravg_file << evnum;
    else lravg_file << mcevnum;
    lravg_file << "," << left_avg_eta << "," << right_avg_eta << std::endl;
    //lravg_file << "," << trueRecoTankTrackLength << std::endl;


  
    // METHOD: STEPPING BACK FROM TANKEXITPT FOR VTX CANDIDATES
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
      //TODO:FIGURE OUT WHY AVG OUT IS HIGHER THAN ITS VALUES
      //std::cout << " AVG CHARGE PER PMT IN/OUT: " << avg_charge_per_pmt_in << ", " << avg_charge_per_pmt_out << std::endl;
      //cpp_file << ",,,,,," << avg_charge_per_pmt_in << "," << avg_charge_per_pmt_out << "loop" << std::endl;

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
  if (gr_vtx_charge_in->GetMaxSize() != 0)
  {
    if (isData) std::cout << " [debug] Saving graphs to canvas... Ev " << evnum << std::endl;
    else std::cout << " [debug] Saving graphs to canvas... Ev " << mcevnum << std::endl;

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
      //ss_effarea_title << evnum;
      ss_effarea_ai_name << evnum;
    }
    else
    {
      //ss_effarea_title << mcevnum;
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
      //ss_fpmt_title << evnum;
      ss_fpmt_ai_name << evnum;
    }
    else
    {
      //ss_fpmt_title << mcevnum;
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
      ss_truetrack << "true track: " << trueRecoTankTrackLength << " cm";
      //gr_eta_ai->SetTitle(ss_truetrack.str().c_str());
    }
    gr_eta_ai->Draw("alp");
    gr_running_avg->Draw("lp");
    //gr_qincone_ai->Draw("lp");
    //gr_qoutcone_ai->Draw("lp");
    c_eta_ai->Modified();
    c_eta_ai->Update();
    c_eta_ai->SetTitle(ss_eta_title.str().c_str());
    c_eta_ai->SetName(ss_eta_ai_name.str().c_str());

    TLegend *legend4 = new TLegend(0.55, 0.65, 0.89, 0.89);
    legend4->AddEntry(gr_eta_ai, "eta", "lp");
    //legend4->AddEntry(gr_qincone_ai, "inside cone", "lp");
    //legend4->AddEntry(gr_qoutcone_ai, "outside cone", "lp");
    legend4->AddEntry(gr_running_avg, "running avg eta", "lp");

    TLine *lBestFitAi = new TLine(bestFitAi, 0, bestFitAi, max_eta+1);
    lBestFitAi->SetLineColor(2);
    lBestFitAi->SetLineWidth(2);
    //lBestFitAi->Draw();
    //legend4->AddEntry(lBestFitAi, "best fit", "l");
    if (!isData)
    { //indicate distance btwn trueVtx and tankExitPoint
      TLine *lTrueTankTrack = new TLine(trueRecoTankTrackLength, 0, trueRecoTankTrackLength, max_eta+1);
      lTrueTankTrack->SetLineColor(4);
      lTrueTankTrack->SetLineWidth(2);
      //lTrueTankTrack->Draw();
      //legend4->AddEntry(lTrueTankTrack, "true tank track length", "l");

      TLine *lLeftAvgEta = new TLine(55, left_avg_eta, trueRecoTankTrackLength, left_avg_eta);
      lLeftAvgEta->SetLineColor(46);
      lLeftAvgEta->SetLineWidth(2);
      //lLeftAvgEta->Draw();
      //legend4->AddEntry(lLeftAvgEta, "avg eta left of true track length", "l");

      TLine *lRightAvgEta = new TLine(trueRecoTankTrackLength, right_avg_eta, 450, right_avg_eta);
      lRightAvgEta->SetLineColor(8);
      lRightAvgEta->SetLineWidth(2);
      //lRightAvgEta->Draw();
      //legend4->AddEntry(lRightAvgEta, "avg eta right of true track length", "l");

      TLine *lAvgEta = new TLine(55, avg_eta, 450, avg_eta);
      lAvgEta->SetLineColor(5);
      lAvgEta->SetLineWidth(2);
      lAvgEta->Draw();
      legend4->AddEntry(lAvgEta, "avg eta", "l");
    }
    legend4->Draw();
    c_eta_ai->Write();
    ss_eta_ai_name << ".png";
    c_eta_ai->SaveAs(ss_eta_ai_name.str().c_str());
  }

  return true;
}


bool MuonFitter::Finalise(){
  // Save output
  root_outp->cd();
  h_alpha->Write();
  h_expected_PE->Write();
  h_phot_inc_angle->Write();
  h_hit_angles->Write();
  h_fitted_tank_track_len->Write();
  h_closest_approach->Write();
  //h_num_mrd_layers->Write();
  //h_clusterhit_x->Write();
  //h_clusterhit_y->Write();
  //h_clusterhit_z->Write();
  //h_clusterhit_detkey->Write();
  //h_clusterhit_timespread->Write();
  //h_clusterhit_time->Write();
  if (!isData)
  {
    h_lastvtx_z->Write();
    h_truevtx_z->Write();
    h_truevtx_angle->Write();
  }
  h_tankexit_to_pmt->Write();
  h_tanktrack_ai->Write();
  h_eff_area_pmt->Write();
  h_fpmt->Write();
  //h_eta_ai->Write();
  h_qincone_truevtx->Write();
  h_qoutcone_truevtx->Write();
  //h_total_pe_hits->Write();
  h_truevtx_recoexit_track->Write();
  h_pmt_charge->Write();

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

  root_outp->Close();
  //pos_file.close();
  //cpp_file.close();
  pehits_file.close();
  lravg_file.close();

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
