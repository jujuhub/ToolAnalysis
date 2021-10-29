#include "P2RunQualityCheck.h"

P2RunQualityCheck::P2RunQualityCheck():Tool(){}


bool P2RunQualityCheck::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // default variable values
  verbosity = 3;
  outpfile_prefix = "PhaseIIRQC_DefaultOutput";
  maxEntries = 10;
  th2_xlim = 10000;
  globalRunNumber = 1;
  max_deq_size = 100;
  globalEntry_i = 0;
  beam_entry_i = 0;
  n_beam_evts = 0;
  n_beam_pot = 0;
  prevpart = 0;
  n_ok_evts = 0;
  e = 0;
  gi = 0;
  pi = 0;

  // Load user-defined variables from config file 
  m_variables.Get("verbosity", verbosity);
  m_variables.Get("OutputFilePrefix", outpfile_prefix);
  m_variables.Get("EntriesPerFile", maxEntries);
  m_variables.Get("TH2XUpperLim", th2_xlim);
  m_variables.Get("RunNumber", globalRunNumber);
  m_variables.Get("MaxDequeSize", max_deq_size);

  if (verbosity > 2) { std::cout << " P2RunQualityCheck tool: config variables loaded. Moving on.." << std::endl; }

  // Load geo in Initialise (or else it gets overwritten)
  bool get_ok = m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry", fGeo);
  if (!get_ok)
  {
    std::cout << " P2RunQualityCheck tool: RUH ROH! Could not load AnnieGeometry!" << std::endl;
    return false;
  }

  // Initialize ROOT stuff
  std::string root_outpfile_ext = ".root";
  std::string root_outpfile_name = outpfile_prefix + root_outpfile_ext;

  if (verbosity > 2) { std::cout << " P2RunQualityCheck tool: root file will be saved as " << root_outpfile_name << std::endl; }

  p2rqc_root_outp = new TFile(root_outpfile_name.c_str(), "RECREATE");

  this->InitTree();
  this->InitHist();
  this->InitGraph();

  // BoostStores stuff
  m_data->CStore.Get("ChannelNumToTankPMTSPEChargeMap", map_chankey2spe);

  std::cout << " P2RunQualityCheck tool: Initalization complete." << std::endl;
  if (verbosity > 3) { std::cout << "   Moving on to world destru--I mean, run quality check........" << std::endl; }

  return true;
}


bool P2RunQualityCheck::Execute(){
  if (verbosity > 5) { std::cout << "MUCH VERBOSE! SO WORDS! MANY TALKS!" << std::endl; }

  //***** Load annie event *****//
  int annieeventexists = m_data->Stores.count("ANNIEEvent");
  if (!annieeventexists)
  {
    std::cerr << "P2RunQualityCheck tool: No ANNIEEvent store!" << std::endl;
    return false;
  }
  globalEntry_i++;


  //***** Retrieve variables from store *****//
  bool get_ok;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("RunNumber", fRunNumber);
  if (!get_ok) { Log("P2RunQualityCheck tool: No RunNumber object in ANNIEEvent! Abort!", v_error, verbosity); fRunNumber = 0;/* return true;*/ }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("SubrunNumber", fSubrunNumber);
  if (!get_ok) { Log("P2RunQualityCheck tool: No SubrunNumber object in ANNIEEvent! Abort!", v_error, verbosity); fSubrunNumber = 0;/* return true;*/ }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("RunType", fRunType);
  if (!get_ok) { Log("P2RunQualityCheck tool: No RunType object in ANNIEEvent! Abort!", v_error, verbosity); fRunType = 0;/* return true;*/ }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("RunStartTime", fStartTime);
  if(!get_ok) { Log("P2RunQualityCheck tool: No RunStartTime object in ANNIEEvent! Abort!", v_error, verbosity); fStartTime = 0;/* return true;*/ }

  get_ok = m_data->Stores["ANNIEEvent"]->Get("EventNumber", fEventNumber);
  if (!get_ok) { Log("P2RunQualityCheck tool: No EventNumber object in ANNIEEvent! Abort!", v_error, verbosity); fEventNumber = 0;/* return true;*/ }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("EventTimeTank", fEventTimeTank);
  if (!get_ok) { Log("P2RunQualityCheck tool: No EventTimeTank object in ANNIEEvent! Abort!", v_error, verbosity); fEventTimeTank = 0;/* return true;*/ }

  get_ok = m_data->Stores["ANNIEEvent"]->Get("PartNumber", fPartNumber);
  if (!get_ok) { Log("P2RunQualityCheck tool: No PartNumber object in ANNIEEvent! Abort!", v_error, verbosity); fPartNumber = -999;/* return true;*/ }

  //this is for looking at 100 (or N) evts per file part
  if ((fPartNumber > prevpart) && (fPartNumber != -999))
  {
    //std::cout << "  [[ DEBUG ]] did this werk??" << std::endl;
    //std::cout << "  [[ DEBUG ]] e = " << e << std::endl;
    e = 0;  //reset entry/evt counter
    prevpart = fPartNumber; //set to current partnumber
  }

  //trigger type
  bool isBeam = false;
  bool isExtTrig = false;
  bool isCCTrig = false;
  bool isNonCCTrig = false;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("TriggerWord", trigword);
  if (!get_ok) { Log("P2RunQualityCheck tool: No TriggerWord object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("TriggerExtended", trigext);
  if (!get_ok) { Log("P2RunQualityCheck tool: No TriggerExtended object in ANNIEEvent! Abort!", v_error, verbosity); trigword = -99;/* return true;*/ }
  if (trigword == 5) { isBeam = true; beam_entry_i++; std::cout << "  [DEBUG] globalEntry_i: " << globalEntry_i << ", beam_entry_i: " << beam_entry_i << std::endl; }
  if (trigext == 0) { isExtTrig = false; }
  if (trigext == 1) { isExtTrig = true; isCCTrig = true; }
  if (trigext == 2) { isExtTrig = true; isNonCCTrig = true; }

  //***** Beam info *****//
  BeamStatus beamstat;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("BeamStatus", beamstat);
  if (get_ok)
  {
    //std::cout << "  [DEBUG] hello beam!" << std::endl;
    if (beamstat.ok())
    {
      h_beam_pot_all->Fill(globalEntry_i, beamstat.pot());
      h_beam_pot_all_zoomx->Fill(globalEntry_i, beamstat.pot());
      h_beam_pot_all_zoomxy->Fill(globalEntry_i, beamstat.pot());
      pot_deq.push_back(beamstat.pot());
      h_beam_ok_all->Fill(globalEntry_i, 1);
      if (isBeam)
      {
        h_beam_ok_beam->Fill(beam_entry_i, 1);
        h_beam_pot_beam->Fill(globalEntry_i, beamstat.pot());
      }
    }
    else
    {
      std::cout << "  [DEBUG] not ok beam pot: " << beamstat.pot() << std::endl;
      h_beam_pot_all->Fill(globalEntry_i, 0.);
      h_beam_pot_all_zoomx->Fill(globalEntry_i, 0.); 
      h_beam_pot_all_zoomxy->Fill(globalEntry_i, 0.);
      pot_deq.push_back(0.);
      h_beam_ok_all->Fill(globalEntry_i, 0);
      if (isBeam)
      {
        h_beam_ok_beam->Fill(beam_entry_i, 0);
        //h_beam_pot_beam->Fill(globalEntry_i, 0.);
      }
    }
    //h_beam_pot_all->Fill(globalEntry_i, beamstat.pot());
    if (isBeam)
    {
      //h_beam_pot_beam->Fill(beam_entry_i, beamstat.pot());
      if (beamstat.pot() < 2.5e12) { n_beam_pot++; }
    }

    //pot_deq.push_back(beamstat.pot());
  }
  else if (!get_ok) { Log("P2RunQualityCheck tool: BeamStatus object in ANNIEEvent! Abort!", v_error, verbosity); return false; }

  if (pot_deq.size() > max_deq_size)
  {
    pot_deq.pop_front();
    std::cout << "  [DEBUG] pot_deq.size(): " << pot_deq.size() << std::endl;
  }

  double deq_mean = 0.;
  double front_pot;
  double back_pot;
  if (pot_deq.size() == max_deq_size)
  {
    //find mean of deque
    double deq_sum = 0.;
    for (size_t i = 0; i < pot_deq.size(); i++)
    {
      deq_sum += pot_deq[i];
    }
    deq_mean = deq_sum / pot_deq.size();
    std::cout << "  [DEBUG] deq_mean: " << deq_mean << std::endl;

    //find gradient
    std::cout << "  [DEBUG] globalEntry: " << globalEntry_i << std::endl;
    front_pot = pot_deq.front();
    std::cout << "  [DEBUG] front_pot: " << front_pot << std::endl;
    back_pot = pot_deq.back();
    std::cout << "  [DEBUG] back_pot: " << back_pot << std::endl;

    gradient_vec.push_back(std::make_pair(globalEntry_i, (back_pot - front_pot) / max_deq_size));
  }


  //***** Check for Veto hits *****//
  fVetoHit = 0;
  bool hasVeto = false;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("TDCData", tdcdata);
  if (!get_ok) { Log("P2RunQualityCheck tool: No TDCData object in ANNIEEvent! Abort!", v_error, verbosity); return true; }

  if (tdcdata->size() > 0)
  {
    Log("P2RunQualityCheck tool: Looping over FMV/MRD hits... looking for Veto activity", v_debug, verbosity);
    for (auto&& anmrdpmt : (*tdcdata))
    {
      unsigned long chankey = anmrdpmt.first;
      Detector *thedetector = fGeo->ChannelToDetector(chankey);
      unsigned long detkey = thedetector->GetDetectorID();
      if (thedetector->GetDetectorElement() == "Veto") { fVetoHit = 1; hasVeto = true;}
    }
  }


  //***** Check for MRD tracks *****//
  bool hasMRDTracks = false;
  int n_tracks_evt = -99;
  get_ok = m_data->Stores["MRDTracks"]->Get("NumMrdTracks", n_tracks_evt);
  if (!get_ok) { Log("P2RunQualityCheck tool: No NumMrdTracks object in ANNIEEvent! Abort!", v_error, verbosity); return true; }

  if (n_tracks_evt != 0)  //should this be > 0? can it ever be negative?
  {
    hasMRDTracks = true;
    if (verbosity > 3) { std::cout << " P2RunQualityCheck tool: NumMrdTracks = " << n_tracks_evt << std::endl; }
  }

  // Get the actual tracks //TODO:UNFINISHED
/*	get_ok = m_data->Stores["MRDTracks"]->Get("NumMrdSubEvents", numsubevs);
	if (!get_ok)
  {
		Log("P2RunQualityCheck Tool: No NumMrdSubEvents in ANNIEEvent!", v_error, verbosity);
		return false;
	}
	get_ok = m_data->Stores["MRDTracks"]->Get("MRDTracks", theMrdTracks);
	if (!get_ok)
  {
		Log("P2RunQualityCheck Tool: No MRDTracks object in MRDTracks BoostStore!", v_error, verbosity);
		Log("MRDTracks store contents:", v_error, verbosity);
		m_data->Stores["MRDTracks"]->Print(false);
		return false;
	}

	for(int tracki=0; tracki<numtracksinev; tracki++){
		BoostStore* thisTrackAsBoostStore = &(theMrdTracks->at(tracki));
		if(verbosity>3) cout<<"track "<<tracki<<" at "<<thisTrackAsBoostStore<<endl;
		// Get the track details from the BoostStore
		thisTrackAsBoostStore->Get("MrdTrackID",MrdTrackID);                    // int
		thisTrackAsBoostStore->Get("MrdSubEventID",MrdSubEventID);              // int
*/


  //***** Looking at individual (not cluster) tank hits *****//
  // for physics, doesn't make sense to look at individual hits
  std::map<unsigned long, std::vector<Hit>> *tank_hits = nullptr;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("Hits", tank_hits);
  if (!get_ok) { Log("P2RunQualityCheck tool: No Hits object in ANNIEEvent! Abort!", v_error, verbosity); return true; }

  fHitPE.clear();
  double total_pmt_charge = 0.;   //for entire entry/event
  for (std::pair<unsigned long, std::vector<Hit>>&& apair : *tank_hits)
  {
    unsigned long chankey = apair.first;
    std::map<int, double>::iterator it = map_chankey2spe.find(chankey);
    if (it != map_chankey2spe.end())
    {
      std::vector<Hit>& thisPMTsHits = apair.second;
      for (Hit &ahit : thisPMTsHits)
      {
        double hit_charge = ahit.GetCharge();
        double t_hitPE = hit_charge / map_chankey2spe.at(chankey);
        fHitPE.push_back(t_hitPE);
        total_pmt_charge += t_hitPE; //or hit_charge
      }
    }
    else { std::cout << " P2RunQualityCheck tool: Couldn't find chankey in chankey2spe map!" << std::endl; }
  }

  if (isBeam && !isExtTrig && !hasVeto)
  {
    n_ok_evts += 1;
    h_tankcharge_all->Fill(globalEntry_i, total_pmt_charge);
    h_tankcharge_all_zoom->Fill(globalEntry_i, total_pmt_charge);
    g_short_all->SetPoint(pi, globalEntry_i, total_pmt_charge);
    pi += 1;

    if (e < maxEntries)
    {
      //std::cout << "  [[ DEBUG ]] what is e? >> " << e << std::endl;
      h_tankcharge_part->Fill(gi, total_pmt_charge);
      h_tankcharge_part_zoom->Fill(gi, total_pmt_charge);
      g_short_part->SetPoint(gi, gi, total_pmt_charge);
      e += 1;
      gi += 1;
      std::cout << "  [[ DEBUG ]] gi = " << gi << " | total_pmt_charge = " << total_pmt_charge << std::endl;
    }
  }


  //***** Looking at tank clusters *****//
  get_ok = m_data->CStore.Get("ClusterMap", m_all_clusters);
  if (!get_ok) { Log("P2RunQualityCheck tool: No ClusterMap object in ANNIEEvent! Did you run the ClusterFinder tool first? Abort!", v_error, verbosity); return false; }

  double total_cluster_charge = 0.;
  if (m_all_clusters->size() != 0)
  {
  for (std::pair<double, std::vector<Hit>>&& cluster_pair : *m_all_clusters)
  {
    double cluster_time = cluster_pair.first;
    std::vector<Hit> cluster_hits = cluster_pair.second;
    double cluster_charge = 0.;
    double cluster_PE = 0.;

    //calculate cluster charge
    for (int i = 0; i < (int)cluster_hits.size(); i++)
    {
      int hitID = cluster_hits.at(i).GetTubeId();
      std::map<int, double>::iterator it =  map_chankey2spe.find(hitID);

      if (it != map_chankey2spe.end())
      {
        double cluster_hitcharge = cluster_hits.at(i).GetCharge();
        double cluster_hitPE = cluster_hitcharge / map_chankey2spe.at(hitID);
        cluster_charge += cluster_hitcharge;
        cluster_PE += cluster_hitPE;
      }
      else { if (verbosity > 2) std::cout << " P2RunQualityCheck tool: FOUND A HIT FOR CHANKEY " << hitID << " BUT NO CONVERSION TO PE AVAILABLE. SKIPPING PE..." << std::endl; }
    }//end loop over cluster hits

    fClusterTime = cluster_pair.first;
    fClusterCharge = cluster_charge;
    fClusterPE = cluster_PE;
    fClusterHits = cluster_hits.size(); //number of hits in cluster

    if (isBeam) { h_clusterTime_all->Fill(cluster_time); }

    if (isBeam && !hasVeto)
    {
      if (isExtTrig) //could be cc or non-cc
      {
        h_clusterTime_long->Fill(cluster_time);
        h_clusterCharge_long->Fill(cluster_charge);
        h_clusterPE_long->Fill(cluster_PE);
      }
      if (!isExtTrig)
      {
        h_clusterTime_short->Fill(cluster_time);
        h_clusterCharge_short->Fill(cluster_charge);
        h_clusterPE_short->Fill(cluster_PE);
      }
    }
    total_cluster_charge += cluster_charge;
    t_cluster->Fill();
  }//end cluster map loop
  }

  t_tank->Fill();
  //t_cluster->Fill();

  return true;
}


bool P2RunQualityCheck::Finalise(){

  //find mean of gradient values
  double avg_grad = 0.;
  for (size_t i = 0; i < gradient_vec.size(); ++i)
  {
    std::pair<int, double> tmpPair;
    tmpPair = gradient_vec.at(i);
    avg_grad += tmpPair.second;
  }
  avg_grad /= (int)gradient_vec.size();

  std::cout << "  [DEBUG] avg_grad: " << avg_grad << " and gradient vector size: " << gradient_vec.size() << std::endl;

  p2rqc_root_outp->cd();
  t_tank->Write("", TObject::kOverwrite);
  t_cluster->Write("", TObject::kOverwrite);
  this->WriteHist();
  this->WriteGraph();

  p2rqc_root_outp->Close();
  delete p2rqc_root_outp;

  std::cout << " P2RunQualityCheck tool: \n  Global number of entries = " << globalEntry_i << "\n   Number of ok evts = " << n_ok_evts << std::endl;
  std::cout << " P2RunQualityCheck tool: \n  Total number of beam entries = " << beam_entry_i << "\n   Number of beam entries with POT < 2.5e12 = " << n_beam_pot << std::endl;
  std::cout << "P2RunQualityCheck tool exitting" << std::endl;

  return true;
}


void P2RunQualityCheck::InitTree()
{
  /***********************************
   ********** Define TTree ***********
   ***********************************/
  
  p2rqc_root_outp->cd();

  t_tank = new TTree("TankTree", "ANNIE Phase II Tank Tree");

  t_tank->Branch("runNumber", &fRunNumber, "runNumber/I");
  t_tank->Branch("subrunNumber", &fSubrunNumber, "subrunNumber/I");
  t_tank->Branch("runType", &fRunType, "runType/I");
  t_tank->Branch("startTime", &fStartTime, "startTime/l");

  t_tank->Branch("eventNumber", &fEventNumber, "eventNumber/I");
  t_tank->Branch("eventTimeTank", &fEventTimeTank, "eventTimeTank/l");
  t_tank->Branch("partNumber", &fPartNumber, "partNumber/I");

  t_tank->Branch("triggerWord", &trigword, "triggerWord/I");
  t_tank->Branch("triggerWordExtended", &trigext, "triggerWordExtended/I");

  t_tank->Branch("vetoHit", &fVetoHit, "vetoHit/I");
  t_tank->Branch("hitPE", &fHitPE);

  //clusters
  t_cluster = new TTree("TankClusterTree", "Tank Cluster Tree");

  t_cluster->Branch("clusterCharge", &fClusterCharge, "clusterCharge/D");
  t_cluster->Branch("clusterTime", &fClusterTime, "clusterTime/D");
  t_cluster->Branch("clusterPE", &fClusterPE, "clusterPE/D");
  t_cluster->Branch("clusterHits", &fClusterHits, "clusterHits/I");
  t_cluster->Branch("clusterChargeBalance", &fClusterChargeBalance, "clusterChargeBalance/D");
  t_cluster->Branch("cluster_hitPE", &fClusterHitPE);

  gROOT->cd();
}

void P2RunQualityCheck::InitHist()
{
  /***********************************
   ******** Define Histograms ********
   ***********************************/

  p2rqc_root_outp->cd();

  h_beam_ok_all = new TH2D("h_beam_ok_all", "beam OK (all evts)", 50000, 0, 50000, 5, 0, 5);
  h_beam_pot_all = new TH2D("h_beam_pot_all", "beam POT (all evts)", 80000, 0, 80000, 1000, 0, 5.0e12);
  h_beam_pot_all_zoomx = new TH2D("h_beam_pot_all_zoomx", "beam POT (all evts, zoomed in on x)", 1000, 0, 1000, 1000, 0, 5.0e12);
  h_beam_pot_all_zoomxy = new TH2D("h_beam_pot_all_zoomxy", "beam POT (all evts, zoomed in on x,y)", 1000, 0, 1000, 1000, 0, 1.0e12);

  h_beam_ok_beam = new TH2D("h_beam_ok_beam", "beam OK (beam-only evts)", 50000, 0, 50000, 5, 0, 5);
  h_beam_pot_beam = new TH2D("h_beam_pot_beam", "beam POT (beam-only evts)", 80000, 0, 80000, 1000, 0, 5.0e12);

  h_tankcharge_part = new TH2D("h_tankcharge_part", "total charge in tank (100 evts)", 1000, 0, 1000, 1000, 0, 300);
  h_tankcharge_part_zoom = new TH2D("h_tankcharge_part_zoom", "total charge in tank (100 evts)", 1000, 0, 1000, 1000, 0, 50);
  h_tankcharge_all = new TH2D("h_tankcharge_all", "total charge in tank (all events)", th2_xlim, 0, th2_xlim, 1000, 0, 300);
  h_tankcharge_all_zoom = new TH2D("h_tankcharge_all_zoom", "total charge in tank (all events)", th2_xlim, 0, th2_xlim, 1000, 0, 50);

  h_clusterTime_all = new TH1D("h_clusterTime_all", "clusterTime of beam events", 1000, 0, 70000);

  h_clusterTime_short = new TH1D("h_clusterTime_short", "clusterTime of beam events (no FMV)", 1000, 0, 75000);
  h_clusterCharge_short = new TH1D("h_clusterCharge_short", "clusterCharge of beam events (no FMV)", 1000, 0, 1);
  h_clusterPE_short = new TH1D("h_clusterPE_short", "clusterPE of beam events (no FMV)", 300, 0, 300);

  h_clusterTime_long = new TH1D("h_clusterTime_long", "clusterTime of extended beam events (no FMV)", 1000, 0, 75000);
  h_clusterCharge_long = new TH1D("h_clusterCharge_long", "clusterCharge of extended beam events (no FMV)", 1000, 0, 1);
  h_clusterPE_long = new TH1D("h_clusterPE_long", "clusterPE of extended beam events (no FMV)", 300, 0, 300);

  gROOT->cd();
}

void P2RunQualityCheck::WriteHist()
{
  /***********************************
   ******** Write Histograms *********
   ***********************************/

  p2rqc_root_outp->cd();

  TDirectory *dir_allhist = p2rqc_root_outp->mkdir("histograms");
  dir_allhist->cd();

  // 2D histograms
  h_beam_ok_all->Write();
  h_beam_pot_all->GetXaxis()->SetTitle("entry no.");
  h_beam_pot_all->GetYaxis()->SetTitle("POT");
  h_beam_pot_all->Write();
  h_beam_pot_all_zoomx->GetXaxis()->SetTitle("entry no.");
  h_beam_pot_all_zoomx->GetYaxis()->SetTitle("POT");
  h_beam_pot_all_zoomx->Write();
  h_beam_pot_all_zoomxy->GetXaxis()->SetTitle("entry no.");
  h_beam_pot_all_zoomxy->GetYaxis()->SetTitle("POT");
  h_beam_pot_all_zoomxy->Write();
  h_beam_ok_beam->Write();
  h_beam_pot_beam->GetXaxis()->SetTitle("entry no.");
  h_beam_pot_beam->GetYaxis()->SetTitle("POT");
  h_beam_pot_beam->Write();

  h_tankcharge_part->GetXaxis()->SetTitle("evt no. [100 per file]");
  h_tankcharge_part->GetYaxis()->SetTitle("total charge (from Hits)");
  h_tankcharge_part->Write();
  h_tankcharge_part_zoom->GetXaxis()->SetTitle("evt no. [100 per file]");
  h_tankcharge_part_zoom->GetYaxis()->SetTitle("total charge (from Hits)");
  h_tankcharge_part_zoom->Write();

  h_tankcharge_all->GetXaxis()->SetTitle("evt no.");
  h_tankcharge_all->GetYaxis()->SetTitle("total charge (from Hits)");
  h_tankcharge_all->Write();
  h_tankcharge_all_zoom->GetXaxis()->SetTitle("evt no.");
  h_tankcharge_all_zoom->GetXaxis()->SetTitle("total charge (from Hits)");
  h_tankcharge_all_zoom->Write();

  // 1D histograms
  h_clusterTime_all->GetXaxis()->SetTitle("clusterTime [ns]");
  h_clusterTime_all->Write();
  h_clusterTime_short->GetXaxis()->SetTitle("clusterTime [ns]");
  h_clusterTime_short->Write();
  h_clusterCharge_short->GetXaxis()->SetTitle("clusterCharge []");
  h_clusterCharge_short->Write();
  h_clusterPE_short->GetXaxis()->SetTitle("clusterPE [#]");
  h_clusterPE_short->Write();

  h_clusterTime_long->GetXaxis()->SetTitle("clusterTime [ns]");
  h_clusterTime_long->Write();
  h_clusterCharge_long->GetXaxis()->SetTitle("clusterCharge []");
  h_clusterCharge_long->Write();
  h_clusterPE_long->GetXaxis()->SetTitle("clusterPE [#]");
  h_clusterPE_long->Write();

  gROOT->cd();
}

void P2RunQualityCheck::InitGraph()
{
  /***********************************
   *********** Init Graphs ***********
   ***********************************/

  c_short_part = new TCanvas("c_short_part", "100 Short Entries", 900, 600);
  g_short_part = new TGraph();
  g_short_part->SetTitle("total tank charge vs. event time");

  c_short_all = new TCanvas("c_short_all", "All Short Entries", 900, 600);
  g_short_all = new TGraph();
  g_short_all->SetTitle("total tank charge vs. event time");

  gROOT->cd();
}

void P2RunQualityCheck::WriteGraph()
{
  /***********************************
   ********** Write Graphs ***********
   ***********************************/

  p2rqc_root_outp->cd();
  TDirectory *dir_allgr = p2rqc_root_outp->mkdir("graphs");
  dir_allgr->cd();

  g_short_part->GetXaxis()->SetTitle("event time");
  g_short_part->GetYaxis()->SetTitle("total tank charge");
  g_short_part->Write("g_short_part");

  c_short_part->cd();
  g_short_part->Draw("AP");
  c_short_part->Write();
  //c_short_part->SaveAs("testcanvas.png");

  g_short_all->GetXaxis()->SetTitle("event time");
  g_short_all->GetYaxis()->SetTitle("total tank charge");
  g_short_all->Write("g_short_all");

  c_short_all->cd();
  g_short_all->Draw("AP");
  c_short_all->Write();

  gROOT->cd();
}

