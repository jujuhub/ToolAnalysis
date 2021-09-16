#include "PhaseIINeutronBG.h"

PhaseIINeutronBG::PhaseIINeutronBG():Tool(){}


bool PhaseIINeutronBG::Initialise(std::string configfile, DataModel &data){

  ///////////////////////// Useful header /////////////////////////
  if (configfile != "") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; // assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // default variable values
  verbosity = 1;
  outpfile_prefix = "PhaseIINBG_default_output";
  min_clusterPE = 5;
  max_clusterPE = 100;

  // load user-defined variables from config file 
  m_variables.Get("verbosity", verbosity);
  m_variables.Get("OutputFilePrefix", outpfile_prefix);
  m_variables.Get("ClusterPEMin", min_clusterPE);
  m_variables.Get("ClusterPEMax", max_clusterPE);
  if (verbosity > 2)
  {
    std::cout << " PhaseIINeutronBG tool: config variables loaded. Moving on.." << std::endl;
  }

  // must load geo in Initialise or else it gets overwritten
  bool get_ok = m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry", fGeo);
  if (!get_ok)
  {
    std::cout << " PhaseIINeutronBG tool: UH OH! Could not load AnnieGeometry!" << std::endl;
    return false;
  }

  // initialize ROOT stuff
  std::string root_outpfile_ext = ".root";
  std::string root_outpfile_name = outpfile_prefix + root_outpfile_ext;

  if (verbosity > 2) { std::cout << " PhaseIINeutronBG tool: root file will be saved as " << root_outpfile_name << std::endl; }

  p2nbg_root_outp = new TFile(root_outpfile_name.c_str(), "RECREATE");

  this->InitTree();
  this->InitHist();

  m_data->CStore.Get("ChannelNumToTankPMTSPEChargeMap", map_chankey2spe);

  std::cout << " PhaseIINeutronBG tool: Initalization complete." << std::endl;
  if (verbosity > 3)
  {
    std::cout << "   Moving on to world destru--I mean, neutron background analysis........" << std::endl;
  }

  return true;
}


bool PhaseIINeutronBG::Execute(){
  if (verbosity > 5) { std::cout << "MUCH VERBOSE! SO WORDS! MANY TALKS!" << std::endl; }

  //  load event info  //
  int annieeventexists = m_data->Stores.count("ANNIEEvent");
  if (!annieeventexists)
  {
    std::cerr << "PhaseIINeutronBG tool: No ANNIEEvent store!" << std::endl;
    return false;
  }

  //  get ptr to ANNIEEvent Store  //
  //auto* annieEvt = m_data->Stores["ANNIEEvent"];
  //annieEvt->Get("TriggerWord", trigword); //this syntax also werks

  bool get_ok;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("RunNumber", fRunNumber);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No RunNumber object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("SubrunNumber", fSubrunNumber);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No SubrunNumber object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("RunType", fRunType);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No RunType object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("RunStartTime", fStartTime);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No RunStartTime object in ANNIEEvent! Abort!", v_error, verbosity); return true; }

  get_ok = m_data->Stores["ANNIEEvent"]->Get("EventNumber", fEventNumber);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No EventNumber object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("EventTimeTank", fEventTimeTank);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No EventTimeTank object in ANNIEEvent! Abort!", v_error, verbosity); return true; }

  bool isBeam = false;
  bool isExtTrig = false;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("TriggerWord", trigword);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No TriggerWord object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("TriggerExtended", trigext);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No TriggerExtended object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  if (trigword == 5) { isBeam = true; }
  if (trigext == 2) { isExtTrig = true; }


  //  check for Veto hits  //
  fVetoHit = 0;
  bool hasVeto = false;

  get_ok = m_data->Stores["ANNIEEvent"]->Get("TDCData", tdcdata);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No TDCData object in ANNIEEvent! Abort!", v_error, verbosity); return true; }

  if (tdcdata->size() > 0)
  {
    Log("PhaseIINeutronBG tool: Looping over FMV/MRD hits... looking for Veto activity", v_debug, verbosity);
    for (auto&& anmrdpmt : (*tdcdata))
    {
      unsigned long chankey = anmrdpmt.first;
      Detector *thedetector = fGeo->ChannelToDetector(chankey);
      unsigned long detkey = thedetector->GetDetectorID();
      if (thedetector->GetDetectorElement() == "Veto") { fVetoHit = 1; hasVeto = true;}
    }
  }


  //  check for MRD tracks  //
  bool hasMRDTracks = false;
  int n_tracks_evt = -99;

  get_ok = m_data->Stores["MRDTracks"]->Get("NumMrdTracks", n_tracks_evt);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No NumMrdTracks object in ANNIEEvent! Abort!", v_error, verbosity); return true; }

  if (n_tracks_evt != 0)
  {
    hasMRDTracks = true;
    std::cout << "  [[ DEBUG ]] NumMrdTracks = " << n_tracks_evt << std::endl;
  }


  //  looking at tank hits  //
  std::map<unsigned long, std::vector<Hit>> *tank_hits = nullptr;

  get_ok = m_data->Stores["ANNIEEvent"]->Get("Hits", tank_hits);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No Hits object in ANNIEEvent! Abort!", v_error, verbosity); return true; }

  bool t_hitPE_ok = true;
  fHitPE.clear();
  for (std::pair<unsigned long, std::vector<Hit>>&& apair : *tank_hits)
  {
    unsigned long chankey = apair.first;
    std::map<int, double>::iterator it = map_chankey2spe.find(chankey);
    if (it != map_chankey2spe.end())
    {
      std::vector<Hit>& thisPMTHits = apair.second;
      for (Hit &ahit : thisPMTHits)
      {
        double hit_charge = ahit.GetCharge();
        double t_hitPE = hit_charge / map_chankey2spe.at(chankey);
        fHitPE.push_back(t_hitPE);

        if (t_hitPE > 5.) { t_hitPE_ok = false; std::cout << "  [[ DEBUG ]] t_hitPE = " << t_hitPE << std::endl; }   //if ANY PMT PE > 5
      }
    }
  }


  //  working w/ tank clusters  //
  get_ok = m_data->CStore.Get("ClusterMap", m_all_clusters);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No ClusterMap object in ANNIEEvent! Abort!", v_error, verbosity); return true; }

  int n_cluster = 0; //NOT the number of clusters in evt; issa label
  int n_prompt_cluster = 0;
  bool hasPromptCluster = false;

  for (std::pair<double, std::vector<Hit>>&& cluster_pair : *m_all_clusters)
  {
    //std::cout << "  [[ DEBUG ]] picklerick" << std::endl;

    // reset variables
    bool isPrompt = false;
    bool isNonCC = false;
    bool pc_hitPE_ok = true;
    double cluster_charge = 0.;
    double cluster_time = cluster_pair.first;
    double cluster_PE = 0;
    std::vector<Hit> cluster_hits = cluster_pair.second;
    fClusterTime = -9999.;
    fClusterHitPE.clear();

    if (cluster_time < 2000.)
    {
      isPrompt = true;  // this is for individual cluster; resets later
      hasPromptCluster = true;  // this is for the entire evt
      n_prompt_cluster += 1;
      //std::cout << "  [[ DEBUG ]] n_prompt_cluster = " << n_prompt_cluster << std::endl;
    }

    // calculate cluster charge
    for (int i = 0; i < cluster_hits.size(); i++)
    {
      bool lonely_PMT = false;  // this will not happen bc min hits to form clusters = 5
      if (cluster_hits.size() < 2) { lonely_PMT = true; std::cout << " [[ DEBUG ]] FOUND A LONELY BOI!" << std::endl; }

      int hit_ID = cluster_hits.at(i).GetTubeId();
      std::map<int, double>::iterator it = map_chankey2spe.find(hit_ID);

      if (it != map_chankey2spe.end())
      {
        double hit_charge = cluster_hits.at(i).GetCharge();
        double c_hitPE = hit_charge / map_chankey2spe.at(hit_ID);
        fClusterHitPE.push_back(c_hitPE);
        cluster_charge += hit_charge;
        cluster_PE += c_hitPE;

        if (isPrompt && c_hitPE > 5.) { pc_hitPE_ok = false; }   //if PROMPT cluster & HIT PE > 5 in cluster
      }
      else
      {
        if (verbosity > 2) { std::cout << "PhaseIINeutronBG tool: FOUND A HIT FOR CHANKEY " << hit_ID << " BUT NO CONVERSION TO PE AVAILABLE. SKIPPING PE..." << std::endl; }
      }
    }
    if (verbosity > 3 && !pc_hitPE_ok) { std::cout << "PhaseIINeutronBG tool: Found a PROMPT cluster with cluster hit PE > 5" << std::endl; }

    // get charge balance
    bool good_class = this->LoadTankClusterClassifiers(cluster_time);
    if (!good_class) { Log("PhaseIINeutronBG tool: NO cluster classifiers..", v_debug, verbosity); }


    fClusterTime = cluster_pair.first;
    fClusterNumber = n_cluster;
    fClusterCharge = cluster_charge;
    fClusterPE = cluster_PE;
    fClusterHits = cluster_hits.size();

    // ALL cluster events
    h_clusterCharge->Fill(cluster_charge);
    h_clusterTime->Fill(cluster_time);
    h_clusterPE->Fill(cluster_PE);


    // labels
    if (pc_hitPE_ok) { isNonCC = true; }


    if (isBeam)
    {
      h_clusterTime_beam->Fill(cluster_time);
      h_clusterPE_beam->Fill(cluster_PE);
    }

    if (isBeam && isNonCC)
    {
      h_clusterTime_nonCCbeam->Fill(cluster_time);
      h_clusterPE_nonCCbeam->Fill(cluster_PE);
    }

    if (cluster_time < 2000.) //ALL prompt
    {
      h_clusterTime_prompt->Fill(cluster_time);
      h_clusterPE_prompt->Fill(cluster_PE);

      if (isBeam && isNonCC)
      {
        h_clusterTime_nonCCbeam_prompt->Fill(cluster_time);
        h_clusterPE_nonCCbeam_prompt->Fill(cluster_PE);

        if (!hasVeto)
        {
          h_clusterTime_nonCCbeam_prompt_noVeto->Fill(cluster_time);
          h_clusterPE_nonCCbeam_prompt_noVeto->Fill(cluster_PE);
        }
      }
    } else if (cluster_time > 2000.) //ALL delayed
    {
      h_clusterTime_delayed->Fill(cluster_time);
      h_clusterPE_delayed->Fill(cluster_PE);

      if (isBeam && isNonCC)
      {
        h_clusterTime_nonCCbeam_delayed->Fill(cluster_time);
        h_clusterPE_nonCCbeam_delayed->Fill(cluster_PE);

        if (!hasVeto)
        {
          h_clusterTime_nonCCbeam_delayed_noVeto->Fill(cluster_time);
          h_clusterPE_nonCCbeam_delayed_noVeto->Fill(cluster_PE);
        }
      }
    }

    t_TankCluster->Fill();
    n_cluster += 1;

  }



  // get run type info - ignore events that aren't beam? (EventSelector)
  // CIT bg calculation - is this per event?
  // identification of skyshine, dirt neutrons


  return true;
}


bool PhaseIINeutronBG::Finalise(){

  p2nbg_root_outp->cd();
  t_TankCluster->Write("", TObject::kOverwrite);
  this->WriteHist();

  p2nbg_root_outp->Close();
  delete p2nbg_root_outp;

  std::cout << "PhaseIINeutronBG tool exitting" << std::endl;
  if (verbosity > 5)
  {
    std::cout << " ...I'm going to go now...." << std::endl;
  }

  return true;
}

void PhaseIINeutronBG::InitTree()
{
  /***********************************
   ********** Define TTree ***********
   ***********************************/
  
  p2nbg_root_outp->cd();

  t_TankCluster = new TTree("phaseIITankClusterTree", "ANNIE Phase II Tank Cluster Tree");

  t_TankCluster->Branch("runNumber", &fRunNumber, "runNumber/I");
  t_TankCluster->Branch("subrunNumber", &fSubrunNumber, "subrunNumber/I");
  t_TankCluster->Branch("runType", &fRunType, "runType/I");
  t_TankCluster->Branch("startTime", &fStartTime, "startTime/l");

  t_TankCluster->Branch("eventNumber", &fEventNumber, "eventNumber/I");
  t_TankCluster->Branch("eventTimeTank", &fEventTimeTank, "eventTimeTank/l");
  t_TankCluster->Branch("clusterNumber", &fClusterNumber, "clusterNumber/I");
  t_TankCluster->Branch("clusterCharge", &fClusterCharge, "clusterCharge/D");
  t_TankCluster->Branch("clusterTime", &fClusterTime, "clusterTime/D");
  t_TankCluster->Branch("clusterPE", &fClusterPE, "clusterPE/D");
  t_TankCluster->Branch("clusterMaxPE", &fClusterMaxPE, "clusterMaxPE/D");
  t_TankCluster->Branch("clusterHits", &fClusterHits, "clusterHits/I");
  t_TankCluster->Branch("clusterChargeBalance", &fClusterChargeBalance, "clusterChargeBalance/D");

  t_TankCluster->Branch("triggerWord", &trigword, "triggerWord/I");
  t_TankCluster->Branch("triggerWordExtended", &trigext, "triggerWordExtended/I");

  t_TankCluster->Branch("vetoHit", &fVetoHit, "vetoHit/I");
  t_TankCluster->Branch("hitPE", &fHitPE);
  t_TankCluster->Branch("cluster_hitPE", &fClusterHitPE);

  gROOT->cd();
}

void PhaseIINeutronBG::InitHist()
{
  /***********************************
   ******** Define Histograms ********
   ***********************************/

  p2nbg_root_outp->cd();

  h_clusterCharge = new TH1F("h_clusterCharge", "clusterCharge", 1000, 0, 1);

  h_clusterTime = new TH1F("h_clusterTime", "clusterTime", 80000, 0, 80000);
  h_clusterTime_beam = new TH1F("h_clusterTime_beam", "clusterTime for BEAM evts", 80000, 0, 80000);
  h_clusterTime_prompt = new TH1F("h_clusterTime_prompt", "clusterTime for ALL PROMPT evts", 2100, 0, 2100);
  h_clusterTime_delayed = new TH1F("h_clusterTime_delayed", "clusterTime for ALL DELAYED evts", 80000, 0, 80000);
  h_clusterTime_nonCCbeam = new TH1F("h_clusterTime_nonCCbeam", "clusterTime for NON-CC BEAM evts", 80000, 0, 80000);
  h_clusterTime_nonCCbeam_prompt = new TH1F("h_clusterTime_nonCCbeam_prompt", "clusterTime for NON-CC BEAM evts - prompt", 2100, 0, 2100);
  h_clusterTime_nonCCbeam_prompt_noVeto = new TH1F("h_clusterTime_nonCCbeam_prompt_noVeto", "clusterTime for NON-CC BEAM evts w NO VETO - prompt", 2100, 0, 2100);
  h_clusterTime_nonCCbeam_delayed = new TH1F("h_clusterTime_nonCCbeam_delayed", "clusterTime for NON-CC BEAM evts - delayed", 80000, 0, 80000);
  h_clusterTime_nonCCbeam_delayed_noVeto = new TH1F("h_clusterTime_nonCCbeam_delayed_noVeto", "clusterTime for NON-CC BEAM evts w NO VETO - delayed", 80000, 0, 80000);

  h_clusterPE = new TH1F("h_clusterPE", "clusterPE", 300, 0, 300);
  h_clusterPE_beam = new TH1F("h_clusterPE_beam", "clusterPE for BEAM evts", 300, 0, 300);
  h_clusterPE_prompt = new TH1F("h_clusterPE_prompt", "clusterPE for ALL PROMPT evts", 300, 0, 300);
  h_clusterPE_delayed = new TH1F("h_clusterPE_delayed", "clusterPE for ALL DELAYED evts", 300, 0, 300);
  h_clusterPE_nonCCbeam = new TH1F("h_clusterPE_nonCCbeam", "clusterPE for NON-CC BEAM evts", 300, 0, 300);
  h_clusterPE_nonCCbeam_prompt = new TH1F("h_clusterPE_nonCCbeam_prompt", "clusterPE for NON-CC BEAM evts - prompt", 300, 0, 300);
  h_clusterPE_nonCCbeam_prompt_noVeto = new TH1F("h_clusterPE_nonCCbeam_prompt_noVeto", "clusterPE for NON-CC BEAM evts w NO VETO - prompt", 300, 0, 300);
  h_clusterPE_nonCCbeam_delayed = new TH1F("h_clusterPE_nonCCbeam_delayed", "clusterPE for NON-CC BEAM evts - delayed", 300, 0, 300);
  h_clusterPE_nonCCbeam_delayed_noVeto = new TH1F("h_clusterPE_nonCCbeam_delayed_noVeto", "clusterPE for NON-CC BEAM evts w NO VETO - delayed", 300, 0, 300);

  gROOT->cd();
}

void PhaseIINeutronBG::WriteHist()
{
  /***********************************
   ******** Write Histograms *********
   ***********************************/

  p2nbg_root_outp->cd();

  TDirectory *dir_allhist = p2nbg_root_outp->mkdir("Histograms");
  dir_allhist->cd();

  h_clusterCharge->Write();

  h_clusterTime->Write();
  h_clusterTime_beam->Write();
  h_clusterTime_prompt->Write();
  h_clusterTime_delayed->Write();
  h_clusterTime_nonCCbeam->Write();
  h_clusterTime_nonCCbeam_prompt->Write();
  h_clusterTime_nonCCbeam_prompt_noVeto->Write();
  h_clusterTime_nonCCbeam_delayed->Write();
  h_clusterTime_nonCCbeam_delayed_noVeto->Write();

  h_clusterPE->Write();
  h_clusterPE_beam->Write();
  h_clusterPE_prompt->Write();
  h_clusterPE_delayed->Write();
  h_clusterPE_nonCCbeam->Write();
  h_clusterPE_nonCCbeam_prompt->Write();
  h_clusterPE_nonCCbeam_prompt_noVeto->Write();
  h_clusterPE_nonCCbeam_delayed->Write();
  h_clusterPE_nonCCbeam_delayed_noVeto->Write();

  gROOT->cd();
}

bool PhaseIINeutronBG::LoadTankClusterClassifiers(double cluster_time)
{
  bool got_ccp = m_data->Stores["ANNIEEvent"]->Get("ClusterChargePoints", cluster_CP);
  bool got_ccb = m_data->Stores["ANNIEEvent"]->Get("ClusterChargeBalances", cluster_CB);
  bool got_cmpe = m_data->Stores["ANNIEEvent"]->Get("ClusterMaxPEs", cluster_maxPEs);
  bool good_class = got_ccp && got_ccb && got_cmpe;
  if (!good_class) { Log("PhaseIINeutronBG tool: One of the charge cluster classifiers is not available", v_debug, verbosity); }
  else
  {
    Log("PhaseIINeutronBG tool: Setting fCluster variables to classifier parameters", v_debug, verbosity);
    fClusterMaxPE = cluster_maxPEs.at(cluster_time);
    fClusterChargeBalance = cluster_CB.at(cluster_time);
  }

  return good_class;
}
