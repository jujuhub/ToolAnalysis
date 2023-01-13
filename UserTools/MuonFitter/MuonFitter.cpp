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

  // Retrieve config variables
  m_variables.Get("verbosity", verbosity);
  m_variables.Get("IsData", isData);

  // Files
  std::string outfile;
  m_variables.Get("OutputFile", outfile);
  std::cout << "MuonFitter Tool: Saving output into " << outfile << std::endl;
  root_outp = new TFile(outfile.c_str(), "RECREATE");

  // Initialize histograms
  h_alpha = new TH1D("h_alpha", "Alpha Value Distribution", 100, 0., 50.);
  h_alpha->GetXaxis()->SetTitle("alpha_i");   //TODO:latex
  h_exp_PE = new TH1D("h_exp_PE", "Expected Num of PE", 100, 0., 100.);
  h_exp_PE->GetXaxis()->SetTitle("expected # of PE");

  // Retrieve Store variables
  m_data->CStore.Get("ChannelNumToTankPMTSPEChargeMap", ChannelKeyToSPEMap);
  auto get_geometry = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry", geom);
  if (!get_geometry)
  {
    Log("MuonFitter Tool: Error retrieving Geometry from ANNIEEvent!", v_error, verbosity);
    return false;
  }
  //tank_radius = geom->GetTankRadius();
  //tank_height = geom->GetTankHalfheight();
  //tank_height /= 2;

  n_tank_pmts = geom->GetNumDetectorsInSet("Tank");
  n_mrd_pmts = geom->GetNumDetectorsInSet("MRD");
  n_veto_pmts = geom->GetNumDetectorsInSet("Veto");
  std::cout<<" [debug] num tank,mrd,veto pmts: "<<n_tank_pmts<<","<<n_mrd_pmts<<","<<n_veto_pmts<<std::endl;

  std::map<std::string, std::map<unsigned long, Detector *>> *Detectors = geom->GetDetectors();

  for (std::map<unsigned long, Detector *>::iterator it = Detectors->at("Tank").begin(); it != Detectors->at("Tank").end(); ++it)
  {
    //loop might be unneeded unless we get diff efficiencies for ea PMT
    Detector *apmt = it->second;
    unsigned long detkey = it->first;
    std::string det_type = apmt->GetDetectorType();

    double det_area = 0.;
    if (det_type == "LUX") { det_area = LUX_AREA; }
    else if (det_type == "ETEL") { det_area = ETEL_AREA; }
    else if (det_type == "Hamamatsu") { det_area = HAMAMATSU_AREA; }
    else if (det_type == "Watchboy") { det_area = WATCHBOY_AREA; }
    else if (det_type == "Watchman") { det_area = WATCHMAN_AREA; }
    else { Log("MuonFitter Tool: Unrecognized detectory type! Setting det_area to 0.", v_error, verbosity); }
    std::cout<<" [debug] det_type: "<<det_type<<", det_area: "<<det_area<<std::endl;

    double eff = 0.25;  //apmt->GetEfficiency(); //TODO:get effective efficiency for Cherenkov spectrum
    m_pmt_eff.insert(std::pair<int, double>(detkey, eff));

    double alpha = ALPHA_FACTOR * eff * det_area;
    m_pmt_alpha.insert(std::pair<int, double>(detkey, alpha));
    std::cout<<" [debug] detkey,alpha: "<<detkey<<","<<alpha<<std::endl;
    h_alpha->Fill(alpha);
  }

  if (verbosity>2) std::cout << "MuonFitter Tool: Initialization complete" << std::endl;

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

  // Reset variables
  //this->ResetVariables();

  //for each event, load clusters and their individual hits
  // (this contains x,y,z info of the hit (PMT) location
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
  if (isData) { cluster_size = (int) m_all_clusters->size(); }
  else { cluster_size = (int) m_all_clusters_MC->size(); }

  std::map<double, std::vector<Hit>>::iterator it_cluster_pair;
  std::map<double, std::vector<MCHit>>::iterator it_cluster_pair_MC;
  if (isData) { it_cluster_pair = (*m_all_clusters).begin(); }
  else { it_cluster_pair_MC = (*m_all_clusters_MC).begin(); }

  bool loop_map = true;
  if (cluster_size == 0) { loop_map = false; }

  while (loop_map)
  {
    //double alpha_i = -99.;

    if (isData)
    {
      double cluster_time = it_cluster_pair->first;
      std::vector<Hit> cluster_hits = it_cluster_pair->second;
      std::cout<<" [debug] cluster_time: "<<cluster_time<<std::endl;

      //do something with hits and cluster time?
      for (int i = 0; i < (int)cluster_hits.size(); i++)
      {
        int chankey = cluster_hits.at(i).GetTubeId();
        std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(chankey);
        if (it != ChannelKeyToSPEMap.end())
        {
          Detector* this_detector = geom->ChannelToDetector(chankey);
          unsigned long detkey = this_detector->GetDetectorID();
          Position det_pos = this_detector->GetDetectorPosition();
          std::string det_type = this_detector->GetDetectorType();
          double hit_charge = cluster_hits.at(i).GetCharge();
          double hit_PE = hit_charge / ChannelKeyToSPEMap.at(chankey);
          double hitX = det_pos.X(); //- tank_center_x;  //TODO:need 2 define
          double hitY = det_pos.Y(); //- tank_center_y;
          double hitZ = det_pos.Z(); //- tank_center_z;

/**          double det_area = 0.;
          if (det_type == "LUX") { det_area = LUX_AREA; }
          else if (det_type == "ETEL") { det_area = ETEL_AREA; }
          else if (det_type == "Hamamatsu") { det_area = HAMAMATSU_AREA; }
          else if (det_type == "Watchboy") { det_area = WATCHBOY_AREA; }
          else if (det_type == "Watchman") { det_area = WATCHMAN_AREA; }
          else { Log("MuonFitter Tool: Unrecognized detectory type!", v_error, verbosity); }
          std::cout<<" [debug] det_type: "<<det_type<<", det_area: "<<det_area<<std::endl; **/

          double det_eff = m_pmt_eff[detkey];
          //alpha_i = ALPHA_FACTOR * det_eff * det_area;
          double alpha_i = -99.;
          alpha_i = m_pmt_alpha[detkey];
          //h_alpha->Fill(alpha_i);

          //expected num of PE
          double exp_PE = alpha_i * (1.+TMath::Cos(0.)) / 2.;  //TOOD:figure out angle that goes into cos
          h_exp_PE->Fill(exp_PE);

          //do something w charge and PE?
        }
      }
    }
    else
    {
      double cluster_time = it_cluster_pair_MC->first;
      std::vector<MCHit> cluster_hits = it_cluster_pair_MC->second;
    }

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

  //test
  //h_alpha->Fill(1.);

  return true;
}


bool MuonFitter::Finalise(){
  // Save output
  root_outp->cd();
  std::cout<<" [debug] saving histogram..."<<std::endl;
  h_alpha->Write();
  h_exp_PE->Write();

  root_outp->Close();

  Log("MuonFitter Tool: Exiting", v_message, verbosity);
  return true;
}
