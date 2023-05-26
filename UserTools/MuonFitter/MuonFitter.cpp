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
  isData = false;

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
  m_variables.Get("Plot3D", plot3d);
  m_variables.Get("Draw3DFMV", draw3d_fmv);
  m_variables.Get("Draw3DMRD", draw3d_mrd);
  m_variables.Get("SaveHistograms", save_hists);

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

  if (save_hists) { canvas_vtxq = new TCanvas("canvas_vtxq", "Total charge vs. track length", 800, 600); }

  h_alpha = new TH1D("h_alpha", "Alpha Value Distribution", 10, 0., 10.);
  h_alpha->GetXaxis()->SetTitle("alpha_i (x1e3)");   //TODO:latex

  h_expected_PE = new TH1D("h_expected_PE", "Expected Num of PE", 100, 0., 100.);
  h_expected_PE->GetXaxis()->SetTitle("expected # of PE");

  h_phot_inc_angle = new TH1D("h_phot_inc_angle", "Photon Incident Angle", 180, 0., 180.);
  h_phot_inc_angle->GetXaxis()->SetTitle("incident angle [deg]");

  h_hit_angles = new TH1D("h_hit_angles", "Hit Angles wrt Vertex and MRD Track Dir", 180, 0., 180.);
  h_hit_angles->GetXaxis()->SetTitle("hit angle [deg]");

  h_tank_track_len = new TH1D("h_tank_track_len", "Length of Tank Tracks", 500, 0., 500.);
  h_tank_track_len->GetXaxis()->SetTitle("track length [cm]");

  h_closest_approach = new TH1D("h_closest_approach", "Distance of Closest Approach (Btwn Vtx & Tank Ctr)", 500, 0., 500.);
  h_closest_approach->GetXaxis()->SetTitle("distance to tank center [cm]");

  h_num_mrd_layers= new TH1D("h_num_mrd_layers", "Num MRD Layers Hit", 50, 0, 50.);
  h_num_mrd_layers->GetXaxis()->SetTitle("# layers");


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

  // --------------------------------------------------
  // --- Retrieve Store variables ---------------------
  // --------------------------------------------------
  m_data->CStore.Get("ChannelNumToTankPMTSPEChargeMap", ChannelKeyToSPEMap);
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
  bBlock->SetLineColor(6);
  EXPH->AddNodeOverlap(bBlock,N++,new TGeoTranslation(0,0,0));

  //draw TANK
  TGeoVolume *annietank = ageom->MakeTubs("annietank", Iron, 0,tank_radius*100.,tank_height*100.,0,360);  //convert to cm
  annietank->SetLineColor(38);
  EXPH->AddNodeOverlap(annietank,N++,new TGeoCombiTrans(tank_center_x*100.,tank_center_y*100.,tank_center_z*100.,new TGeoRotation("annietank",0,90,0)));
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
    //x_pmt.insert(std::pair<int,double>(detkey,position_PMT.X()-tank_center_x));
    //y_pmt.insert(std::pair<int,double>(detkey,position_PMT.Y()-tank_center_y));
    //z_pmt.insert(std::pair<int,double>(detkey,position_PMT.Z()-tank_center_z));
    // pmt xyz
    x_pmt.insert(std::pair<unsigned long, double>(detkey, position_PMT.X()*100.));
    y_pmt.insert(std::pair<unsigned long, double>(detkey, position_PMT.Y()*100.));
    z_pmt.insert(std::pair<unsigned long, double>(detkey, position_PMT.Z()*100.));
    // pmt orientation
    Direction direction_PMT = apmt->GetDetectorDirection();
    x_pmt_dir.insert(std::pair<unsigned long, double>(detkey, direction_PMT.X()*100.));
    y_pmt_dir.insert(std::pair<unsigned long, double>(detkey, direction_PMT.Y()*100.));
    z_pmt_dir.insert(std::pair<unsigned long, double>(detkey, direction_PMT.Z()*100.));

    // ANNIE in 3D: Drawing TANK PMTs
    sprintf(blockName,"tank_pmt%lu",detkey);
    bBlock = ageom->MakeSphere(blockName, Iron, 0,1.5,0,180,0,360); //TODO:include PMT type and location condition
    bBlock->SetLineColor(41);
    //EXPH->AddNodeOverlap(bBlock,1,new TGeoTranslation(x_pmt[detkey]*100., y_pmt[detkey]*100., z_pmt[detkey]*100.));   //these pos coords work w/ tank center (0,0,0)
    EXPH->AddNodeOverlap(bBlock,detkey,new TGeoTranslation(position_PMT.X()*100., position_PMT.Y()*100., position_PMT.Z()*100.));   //these pos coords work w/ tank_center_x/y/z*100.
    detkey_to_node.insert(std::pair<unsigned long, int>(detkey, N++));

    // Get PMT areas & efficiency values
    double det_area = 0.;
    if (det_type == "LUX") { det_area = LUX_AREA; }
    else if (det_type == "ETEL") { det_area = ETEL_AREA; }
    else if (det_type == "Hamamatsu") { det_area = HAMAMATSU_AREA; }
    else if (det_type == "Watchboy") { det_area = WATCHBOY_AREA; }
    else if (det_type == "Watchman") { det_area = WATCHMAN_AREA; }
    else { Log("MuonFitter Tool: Unrecognized detectory type! Setting det_area to 0.", v_error, verbosity); }
    //std::cout<<" [debug] det_type: "<<det_type<<", det_area: "<<det_area<<std::endl;

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
      EXPH->AddNodeOverlap(bBlock,detkey,new TGeoTranslation(position_MRD.X()*100., position_MRD.Y()*100., position_MRD.Z()*100.));
      detkey_to_node.insert(std::pair<unsigned long, int>(detkey, N++));

      if (verbosity > v_debug)
      {
        std::cout << " [debug] position_MRD (x,y,z): " << position_MRD.X() << "," << position_MRD.Y() << "," << position_MRD.Z() << std::endl;
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
      EXPH->AddNodeOverlap(bBlock,detkey,new TGeoTranslation(position_FMV.X()*100., position_FMV.Y()*100., position_FMV.Z()*100.));
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
  pos_file.open(pos_fname.c_str());
  pos_file << "##evnum,startX,startY,startZ,stopX,stopY,stopZ" << std::endl;

  std::string cpp_fname = "charge_per_pmt.txt";
  cpp_file.open(cpp_fname.c_str());
  cpp_file << "##evnum,nVtx,Qin,Qout,avgSumIn,avgSumOut,avgQin,avgQout" << std::endl;

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
  // --- Get event info -------------------------------
  // --------------------------------------------------
  int get_ok;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("EventNumber", evnum);
  if (not get_ok) { Log("MuonFitter Tool: Error retrieving EventNumber, true from ANNIEEvent!", v_error, verbosity); return false; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("RunNumber", runnumber);
  if (not get_ok) { Log("EventDisplay tool: Error retrieving RunNumber, true from ANNIEEvent!", v_error, verbosity); return false; }


  // --------------------------------------------------
  // --- Reset variables ------------------------------
  // --------------------------------------------------
  //this->ResetVariables();   //TODO:write ResetVariables func
  bool drawEvent = false;
  if (plot3d) { reset_3d(); }
  if (save_hists) { canvas_vtxq->Clear(); }

  if (h_vtx_charge) { h_vtx_charge->Delete(); }
  h_vtx_charge = new TH2D("h_vtx_charge", "charge as function of vertex Z", 300, 100., 400., 100, 0., 5.);
  h_vtx_charge->GetXaxis()->SetTitle("vertex Z [cm]");
  h_vtx_charge->GetYaxis()->SetTitle("charge [#PE]");

  if (h_charge_per_pmt) { h_charge_per_pmt->Delete(); }
  h_charge_per_pmt = new TH2D("h_charge_per_pmt", "charge per pmt as function of vertex Z", 300, 100., 400., 100, 0., 0.1);
  h_charge_per_pmt->GetXaxis()->SetTitle("vertex Z [cm]");
  h_charge_per_pmt->GetYaxis()->SetTitle("charge [#PE]");

  gr_vtx_charge_in->Set(0);
  gr_vtx_charge_out->Set(0);
  gr_qdensity_in->Set(0);
  gr_qdensity_out->Set(0);


  // --------------------------------------------------
  // --- Check for FMV hits ---------------------------
  // --------------------------------------------------
  bool hasVeto = false;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("TDCData", tdcdata);
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
  get_ok = m_data->Stores["MRDTracks"]->Get("MRDTracks", mrdTracks);

  if (!get_ok)
  {
    Log("MuonFitter Tool: Couldn't retrieve MRD tracks info. Did you run TimeClustering/FindMRDTracks first?", v_debug, verbosity);
    return false;
  }

  //skip if more than 1 track in MRD
  if (numTracksInEv != 1)
  {
    Log("MuonFitter Tool: More than 1 reconstructed track found!", v_debug, verbosity);
    m_data->CStore.Set("DrawEventDisplay", drawEvent);
    return true;
  }

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
    tracklength = sqrt(pow((mrdStopVertex.X()-mrdStartVertex.X()),2)+pow(mrdStopVertex.Y()-mrdStartVertex.Y(),2)+pow(mrdStopVertex.Z()-mrdStartVertex.Z(),2));

    std::cout << "  [debug] mrdStartVertex: " << mrdStartVertex.X() << "," << mrdStartVertex.Y() << "," << mrdStartVertex.Z() << std::endl;
    std::cout << "  [debug] mrdStopVertex: " << mrdStopVertex.X() << "," << mrdStopVertex.Y() << "," << mrdStopVertex.Z() << std::endl;
    std::cout << "  [debug] trackAngle: " << trackAngle << std::endl;
    std::cout << "  [debug] numLayersHit: " << numLayersHit << std::endl;
    std::cout << "  [debug] penetrationDepth: " << penetrationDepth << std::endl;
    std::cout << "  [debug] tracklength: " << tracklength << std::endl;
    std::cout << "  [debug] energyLoss: " << energyLoss << std::endl;

    //save MRD start/stop coordinates to file
    pos_file << evnum << "," << mrdStartVertex.X() << ","
                             << mrdStartVertex.Y() << ","
                             << mrdStartVertex.Z() << ","
                             << mrdStopVertex.X() << ","
                             << mrdStopVertex.Y() << ","
                             << mrdStopVertex.Z() << std::endl;
    h_num_mrd_layers->Fill(numLayersHit);
  }
  //XXX:check whether MRD track falls into the fiducial volume


  // Create vectors for MRD start/stop
  TVector3 mrdStart(mrdStartVertex.X()*100., mrdStartVertex.Y()*100., mrdStartVertex.Z()*100.);   //[cm]
  TVector3 mrdStop(mrdStopVertex.X()*100., mrdStopVertex.Y()*100., mrdStopVertex.Z()*100.);
  TVector3 mrdTrackDir = (mrdStop - mrdStart).Unit();

  //add mrdStartVertex to 3D geo
  sprintf(blockName,"mrdStart_%u",evnum);
  bBlock = ageom->MakeSphere(blockName, Iron, 0,1.5,0,180,0,360);
  bBlock->SetLineColor(8);
  EXPH->AddNodeOverlap(bBlock,69,new TGeoTranslation(mrdStart.X(), mrdStart.Y(), mrdStart.Z()));
  N++;

  //add mrdStopVertex to 3D geo
  sprintf(blockName,"mrdStop_%u",evnum);
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
    TVector3 v = mrdStart - c*10.*mrdTrackDir;
    if (c <= 5) { vtxCandidates.push_back(v); }
    else if (c > 5)
    {
      //check that vtx is contained in tank cylinder
      if ((pow(v.X()-tank_center_x*100.,2) + pow(v.Z()-tank_center_z*100.,2) >= pow(tank_radius*100.,2)) || (v.Y() > tank_center_y*100.+tank_height*100.) || (v.Y() < tank_center_y*100.-tank_height*100.))
      {
        outsideTank = true;
      }
      else { vtxCandidates.push_back(v); }

      //check if any vtx is in FV
      if (pow(v.X()-tank_center_x*100.,2) + pow(v.Y()-tank_center_y*100.,2) + pow(v.Z()-tank_center_z*100.,2) <= pow(tank_radius*100.-50.,2))
      {
        std::cout << " [debug] (distance from tank center)^2: " << pow(v.X()-tank_center_x*100.,2) + pow(v.Y()-tank_center_y*100.,2) + pow(v.Z()-tank_center_z*100.,2) << std::endl;
        inFV = true;
      }
    }

    // add vertex to 3D geo
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

  has_clusters = m_data->CStore.Get("ClusterMapDetkey", m_all_clusters_detkeys);
  if (!has_clusters)
  {
    std::cout << "MuonFitter Tool: No cluster detkeys found in CStore! Did you run ClusterFinder tool before this?" << std::endl;
    return false;
  }

  Log("MuonFitter Tool: Accessing pairs in all_clusters map", v_debug, verbosity);
  int cluster_num = 0;
  int cluster_size = 0;
  if (isData) { cluster_size = (int) m_all_clusters->size(); std::cout << " [debug] cluster_size: " << cluster_size << std::endl; }
  else { cluster_size = (int) m_all_clusters_MC->size(); }

  std::map<double, std::vector<Hit>>::iterator it_cluster_pair;
  std::map<double, std::vector<MCHit>>::iterator it_cluster_pair_MC;
  if (isData) { it_cluster_pair = (*m_all_clusters).begin(); }
  else { it_cluster_pair_MC = (*m_all_clusters_MC).begin(); }

  bool loop_map = true;
  if (cluster_size == 0) { loop_map = false; }  //TODO:require cluster_size = 1?

  double avg_charge_per_pmt_in = 0., avg_charge_per_pmt_out = 0.;
  int max_vertex = -1;
  double max_charge = 0.;
  while (loop_map)
  {
    if (isData)
    {
      double cluster_time = it_cluster_pair->first;
      std::vector<Hit> cluster_hits = it_cluster_pair->second;
      //std::cout<<"  [debug] cluster_time: "<<cluster_time<<std::endl;

      // Make charge and time cuts
      if (cluster_time > 2000.)   //if not prompt, move on to next event
      {
        it_cluster_pair++;
        if (it_cluster_pair == (*m_all_clusters).end()) { loop_map = false; }
        continue;
      }

      //do something with hits
      std::vector<double> x_hits, y_hits, z_hits;
      double cluster_charge = 0.;   //TODO:make cluster or individual PMT charge cut?
  
      std::map<int, double> m_vtx_charge_incone;   //map of total charge seen at ea vtx candidate
      std::map<int, double> m_vtx_charge_per_pmt_incone;
      std::map<int, double> m_vtx_detkey_incone;
      std::map<int, std::vector<double>> m_vtx_vec_charge_incone;   //map of charge vec seen at ea vtx candidate
      std::map<int, std::vector<unsigned long>> m_vtx_vec_detkey_incone;
      std::map<int, std::vector<double>> m_vtx_vec_charge_outcone;
      std::map<int, std::vector<unsigned long>> m_vtx_vec_detkey_outcone;
      //double avg_charge_per_pmt_in = 0., avg_charge_per_pmt_out = 0.;
      for (int c = 0; c < vtxCandidates.size(); ++c)
      {
        //for each vtx candidate
        TVector3 vtx_to_mrd = TVector3(mrdStart.X(), mrdStart.Y(), mrdStart.Z()) - vtxCandidates.at(c);

        // Find all pmts that fall in cone, hit or not
        //std::cout << "  [debug] all detkeys in cone (hit or not) for vtx c" << c << ": ";
        std::vector<unsigned long> v_all_detkeys_incone;
        for (unsigned long detkey = 332; detkey < 465; ++detkey)
        {
          TVector3 vtx_to_pmt = TVector3(x_pmt[detkey], y_pmt[detkey], z_pmt[detkey]) - vtxCandidates.at(c);
          double pmtAngle = vtx_to_pmt.Angle(vtx_to_mrd)*180./TMath::Pi();
          if (pmtAngle < CHER_ANGLE_DEG)
          {
            v_all_detkeys_incone.push_back(detkey);
            //std::cout << detkey << ", ";
          }
        }
        //std::cout << std::endl;
        //std::cout << " [debug] total detkeys in cone: " << v_all_detkeys_incone.size() << std::endl;
        gr_vtx_detkey_in->SetPoint(c, vtxCandidates.at(c).Z(), v_all_detkeys_incone.size());
        gr_vtx_detkey_out->SetPoint(c, vtxCandidates.at(c).Z(), 132-v_all_detkeys_incone.size());

        std::vector<unsigned long> v_detkeys_incone;
        std::vector<unsigned long> v_detkeys_outcone;
        std::vector<double> v_charges_incone;
        std::vector<double> v_charges_outcone;
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
            double hitX = det_pos.X(); //- tank_center_x;   //[m]
            double hitY = det_pos.Y(); //- tank_center_y;
            double hitZ = det_pos.Z(); //- tank_center_z;
            x_hits.push_back(hitX);
            y_hits.push_back(hitY);
            z_hits.push_back(hitZ);

            Direction dirPMT = this_detector->GetDetectorDirection();

            //angle of vector from vtx to PMT wrt to track
            TVector3 vtx_to_hit = TVector3(hitX*100., hitY*100., hitZ*100.) - vtxCandidates.at(c);
            double pmtAngleRad = vtx_to_hit.Angle(vtx_to_mrd);
            double pmtAngleDeg = pmtAngleRad*180./TMath::Pi();
            h_hit_angles->Fill(pmtAngleDeg);

            TVector3 pmt_dir = TVector3(hitX*100., hitY*100., hitZ*100.) - TVector3(dirPMT.X()*100., dirPMT.Y()*100., dirPMT.Z()*100.);

            //debug: testing pmt dir
/*            TVector3 p1 = pmt_dir.Unit();
            TVector3 p2 = TVector3(hitX*100., hitY*100., hitZ*100.) - 10.*p1;
            if (detkey == 400)
            {
              sprintf(blockName,"pmtDir%u",N);
              bBlock = ageom->MakeSphere(blockName, Iron, 0,1.5,0,180,0,360);
              bBlock->SetLineColor(2);
              EXPH->AddNodeOverlap(bBlock,69,new TGeoTranslation(p2.X(), p2.Y(), p2.Z()));
              N++;
            }
*/
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

              // Find the incident angle of photon wrt pmt
              double photInAngleRad = pmt_dir.Angle(vtx_to_hit);
              double photInAngleDeg = photInAngleRad*180./TMath::Pi();
              h_phot_inc_angle->Fill(photInAngleDeg);
              //std::cout << "  [debug] incident photon angle (rad, deg): " << photInAngleRad << ", " << photInAngleDeg << std::endl;
            }
          }
        } //done going through cluster hits for this vtx

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
        cpp_file << evnum << "," << vtxCandidates.size() <<",";
        int tot_pmts_incone = v_all_detkeys_incone.size();
        if (tot_pmts_incone == 0)
        {
          gr_qdensity_in->SetPoint(c, vtxCandidates.at(c).Z(), 0.);
          avg_charge_per_pmt_in += 0.;
          m_vtx_charge_per_pmt_incone[c] = 0.;
          //std::cout << " ev" << evnum << " charge per pmt (in): 0." << std::endl;
          cpp_file << "0.,";
        }
        else
        {
          gr_qdensity_in->SetPoint(c, vtxCandidates.at(c).Z(), tot_charge_incone / tot_pmts_incone);
          avg_charge_per_pmt_in += tot_charge_incone / tot_pmts_incone;
          m_vtx_charge_per_pmt_incone[c] = tot_charge_incone / tot_pmts_incone;
          //std::cout << " ev" << evnum << " charge per pmt (in): " << tot_charge_incone/tot_pmts_incone << std::endl;
          cpp_file << tot_charge_incone/tot_pmts_incone << ",";
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
        //for (auto const &pair: m_hit_pmts_incone) { std::cout << "{" << pair.first << ":" << pair.second <<"}" << std::endl; }


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
          cpp_file << "0.,";
        }
        else
        {
          gr_qdensity_out->SetPoint(c, vtxCandidates.at(c).Z(), tot_charge_outcone / tot_pmts_outcone);
          avg_charge_per_pmt_out += tot_charge_outcone / tot_pmts_outcone;
          //std::cout << " ev" << evnum << " charge per pmt (out): " << tot_charge_outcone / tot_pmts_outcone << std::endl;
          cpp_file << tot_charge_outcone/tot_pmts_outcone << ",";
        }

        //std::cout << " [debug] ev" << evnum << " avg_charge_per_pmt_in/out, num vtx: " << avg_charge_per_pmt_in << ", " << avg_charge_per_pmt_out << ", " << vtxCandidates.size() << std::endl;
        cpp_file << avg_charge_per_pmt_in << "," << avg_charge_per_pmt_out << ",," << std::endl;

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
        std::map<unsigned long, double> m_pmt_charge_outcone;
        for (int p = 0; p < v_detkeys_outcone.size(); ++p)
        {
          unsigned long detkey = v_detkeys_outcone.at(p);
          m_pmt_charge_outcone[detkey] += v_charges_outcone.at(p);
          //std::cout << "  [debug] m_pmt_charge_outcone[" << detkey << "]: " << m_pmt_charge_outcone[detkey] << std::endl;
        }
        for (auto const &pair: m_pmt_charge_outcone)
        {
          TGeoNode *node = EXPH->GetNode(detkey_to_node[pair.first]);
          node->GetVolume()->SetLineColor(30);
          TGeoSphere *sphere = (TGeoSphere *)(node->GetVolume()->GetShape());
          sphere->SetSphDimensions(0,20.*(pair.second/(tot_charge_outcone+tot_charge_incone))+5.,0,180,0,360);
        } */
      } //done going thru vtxCandidates

      //compare avg charge per pmt
      avg_charge_per_pmt_in = avg_charge_per_pmt_in / vtxCandidates.size();
      avg_charge_per_pmt_out = avg_charge_per_pmt_out / vtxCandidates.size();
      //TODO:FIGURE OUT WHY AVG OUT IS HIGHER THAN ITS VALUES
      std::cout << " AVG CHARGE PER PMT IN/OUT: " << avg_charge_per_pmt_in << ", " << avg_charge_per_pmt_out << std::endl;
      cpp_file << ",,,,,," << avg_charge_per_pmt_in << "," << avg_charge_per_pmt_out << "loop" << std::endl;

      // --------------------------------------------------
      // --- Find best vertex -----------------------------
      // --------------------------------------------------
      max_charge = m_vtx_charge_per_pmt_incone[0];
      for (auto const &pair: m_vtx_charge_per_pmt_incone)
      {
        //std::cout << " [debug] current charge, vtx: " << pair.second <<", " << pair.first << std::endl;
        if (pair.second > max_charge)
        {
          max_charge = pair.second;
          max_vertex = pair.first;
        }
        //std::cout << " [debug] max charge, vtx: " << max_charge <<", " << max_vertex << std::endl;
      }
      auto X = std::max_element(m_vtx_charge_per_pmt_incone.begin(), m_vtx_charge_per_pmt_incone.end(), [] (const pair<int, double>& p1, const pair<int, double>&p2) { return p1.second < p2.second; });
      max_vertex = X->first;
      max_charge = X->second;
      std::cout << " Max charge using max_element method: " << max_vertex << ", " << max_charge << std::endl;

      TVector3 bestVtx(vtxCandidates.at(max_vertex).X(), vtxCandidates.at(max_vertex).Y(), vtxCandidates.at(max_vertex).Z());
      double tankTrackLength = TMath::Sqrt(pow(mrdStart.X()-bestVtx.X(),2) + pow(mrdStart.Y()-bestVtx.Y(),2) + pow(mrdStart.Z()-bestVtx.Z(),2));
      //std::cout << " [debug] tankTrackLength: " << tankTrackLength << std::endl;
      h_tank_track_len->Fill(tankTrackLength);
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
      sprintf(blockName,"tankVtx_%d", c);
      bBlock = ageom->MakeSphere(blockName, Iron, 0,1.5,0,180,0,360);
      bBlock->SetLineColor(6);
      EXPH->AddNodeOverlap(bBlock,69,new TGeoTranslation(vtxCandidates.at(max_vertex).X(), vtxCandidates.at(max_vertex).Y(), vtxCandidates.at(max_vertex).Z()));
      N++;

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
        node->GetVolume()->SetLineColor(30);
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
        h_vtx_charge->Fill(vtxCandidates.at(it->first).Z(), it->second);
        gr_vtx_charge_in->SetPoint(pt, vtxCandidates.at(it->first).Z(), it->second);
        ++pt;
      }*/
    }
    else  //is MC
    {
      double cluster_time = it_cluster_pair_MC->first;
      std::vector<MCHit> cluster_hits = it_cluster_pair_MC->second;
    }

    // Move on to next cluster or exit
    if (isData)
    {
      it_cluster_pair++;
      if (it_cluster_pair == (*m_all_clusters).end()) { loop_map = false; }
    }
    else
    {
      it_cluster_pair_MC++;
      if (it_cluster_pair_MC == (*m_all_clusters_MC).end()) { loop_map = false; }
    }
  } //end while loop_map

  //make some cuts? need to select muon candidates
  // -beam conditions,timing,charge,nhits,etc
  //calculate alpha for each hit

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
    if (!isData)
    {
      int dummyev = 69;
      ss_evdisplay3d_title << "evdisplay3d_ev" << dummyev;
      ss_evdisplay3d_name << "canvas_evdisplay3d_ev" << dummyev;
    }
    else
    {
      ss_evdisplay3d_title << "evdisplay3d_ev" << evnum;
      ss_evdisplay3d_name << "canvas_evdisplay3d_ev" << evnum;
    }
    canvas_3d->SetTitle(ss_evdisplay3d_title.str().c_str());
    canvas_3d->SetName(ss_evdisplay3d_name.str().c_str());
    ageom->DrawTracks();
    canvas_3d->Modified();
    canvas_3d->Update();
    canvas_3d->Write();
  }

  // Write histograms
  /*
  std::stringstream ss_hvtxq_title, ss_hvtxq_name;
  ss_hvtxq_title << "Histogram: Vtx v Q - Ev_" << evnum;
  ss_hvtxq_name << "h_vtx_q_ev" << evnum;
  h_vtx_charge->SetTitle(ss_hvtxq_title.str().c_str());
  h_vtx_charge->SetName(ss_hvtxq_name.str().c_str());
  std::stringstream ss_hqpmt_title, ss_hqpmt_name;
  ss_hqpmt_title << "Histogram: Charge Per PMT - Ev_" << evnum;
  ss_hqpmt_name << "h_charge_per_pmt_ev" << evnum;
  h_charge_per_pmt->SetTitle(ss_hqpmt_title.str().c_str());
  h_charge_per_pmt->SetName(ss_hqpmt_name.str().c_str());

  if (save_hists)
  {
    std::cout << "  [debug] Saving histograms to canvas..." << std::endl;
    std::stringstream ss_cvtxq_title, ss_cvtxq_name;
    ss_cvtxq_title << "Histogram: Vtx v Q - Ev_" << evnum;
    ss_cvtxq_name << "canvas_vtxq_ev" << evnum;
    canvas_vtxq->cd(1);
    h_vtx_charge->Draw("P");
    canvas_vtxq->Modified();
    canvas_vtxq->Update();
    canvas_vtxq->SetTitle(ss_cvtxq_title.str().c_str());
    canvas_vtxq->SetName(ss_cvtxq_name.str().c_str());
  }

  std::cout << "  [debug] hist entries: " << h_vtx_charge->GetEntries() << std::endl;
  if (h_vtx_charge->GetEntries() != 0)
  {
    if (save_hists) { canvas_vtxq->Write(); }
    if (!save_hists) { h_vtx_charge->Draw("P"); }
    h_vtx_charge->Write();
  }
  if (h_charge_per_pmt->GetEntries() != 0)
  {
    h_charge_per_pmt->Draw("P");
    h_charge_per_pmt->Write();
  }
*/
  //saving graph
  if (gr_vtx_charge_in->GetMaxSize() != 0)
  {
    std::cout << " [debug] Saving graphs to canvas... Ev " << evnum << std::endl;
    // total charge
    c_vtx_charge->cd();
    std::stringstream ss_cvtxq_grtitle, ss_cvtxq_grname;
    ss_cvtxq_grtitle << "Total Charge Seen at Each Vertex Ev_" << evnum;
    ss_cvtxq_grname << "c_vtx_charge_ev" << evnum;
    gr_vtx_charge_in->Draw("alp");
    gr_vtx_charge_out->Draw("lp");
    c_vtx_charge->Modified();
    c_vtx_charge->Update();
    c_vtx_charge->SetTitle(ss_cvtxq_grtitle.str().c_str());
    c_vtx_charge->SetName(ss_cvtxq_grname.str().c_str());
    TLegend *legend1 = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend1->AddEntry(gr_vtx_charge_in, "inside cone", "lp");
    legend1->AddEntry(gr_vtx_charge_out, "outside cone", "lp");
    legend1->Draw();
    c_vtx_charge->Write();

    // charge per PMT
    c_charge_per_pmt->cd();
    std::stringstream ss_cv_qdensity_grtitle, ss_cv_qdensity_grname;
    ss_cv_qdensity_grtitle << "Charge Per PMT at Each Vertex Ev_" << evnum;
    ss_cv_qdensity_grname << "c_charge_per_pmt_ev" << evnum;
    gr_qdensity_in->Draw("alp");
    gr_qdensity_out->Draw("lp");
    //draw marker over largest charge per pmt value
    TMarker *mkr_max = new TMarker(vtxCandidates.at(max_vertex).Z(), max_charge, 43);
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
    l1->Draw();
    TLine *l2 = new TLine(l_zmin, avg_charge_per_pmt_out, l_zmax, avg_charge_per_pmt_out);
    l2->SetLineColor(9);
    l2->SetLineWidth(2);
    l2->Draw();
    TLegend *legend2 = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend2->AddEntry(gr_qdensity_in, "inside cone", "lp");
    legend2->AddEntry(gr_qdensity_out, "outside cone", "lp");
    legend2->Draw();
    c_charge_per_pmt->Write();
    cpp_file << ",,,,,," << avg_charge_per_pmt_in << "," << avg_charge_per_pmt_out << std::endl;

    // total detkeys
    c_vtx_detkey->cd();
    std::stringstream ss_vtx_detkey_title, ss_vtx_detkey_name;
    ss_vtx_detkey_title << "Num PMTs Seen at Each Vertex Ev_" << evnum;
    ss_vtx_detkey_name << "c_vtx_detkey_ev" << evnum;
    gr_vtx_detkey_in->Draw("alp");
    gr_vtx_detkey_out->Draw("lp");
    c_vtx_detkey->Modified();
    c_vtx_detkey->Update();
    c_vtx_detkey->SetTitle(ss_vtx_detkey_title.str().c_str());
    c_vtx_detkey->SetName(ss_vtx_detkey_name.str().c_str());
    c_vtx_detkey->Write();
  }

  return true;
}


bool MuonFitter::Finalise(){
  // Save output
  root_outp->cd();
  h_alpha->Write();
  h_expected_PE->Write();
  //h_vtx_charge->Write();
  h_phot_inc_angle->Write();
  h_hit_angles->Write();
  h_tank_track_len->Write();
  h_closest_approach->Write();
  h_num_mrd_layers->Write();

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

  root_outp->Close();
  pos_file.close();
  cpp_file.close();

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
