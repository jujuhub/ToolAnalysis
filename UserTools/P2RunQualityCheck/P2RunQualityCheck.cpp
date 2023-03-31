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
  globalRunNumber = 1337;
  globalEntry = 0;
  beamEntry = 0;
  n_beam_evts = 0;
  n_beam_pot = 0;
  prevpart = 0;
  n_ok_evts = 0;

  // Load user-defined variables from config file 
  m_variables.Get("verbosity", verbosity);
  m_variables.Get("EntriesPerFile", maxEntries);
  m_variables.Get("OutputFilePrefix", outpfile_prefix);
  m_variables.Get("RunNumber", globalRunNumber);

  if (verbosity > 4) { std::cout << " P2RunQualityCheck tool: config variables loaded. Moving on.." << std::endl; }

  std::string pot_fname = "R" + std::to_string(globalRunNumber) + "_beampot.txt";
  std::cout << "P2RunQualityCheck tool: Saving POT file as... " << pot_fname << std::endl;
  beampot_file.open(pot_fname.c_str());
  beampot_file << "#EventTimeTank BeamPOT" << std::endl;

  // BoostStore stuff
  bool get_ok;
  //load geo in Initialise (or else it gets overwritten)
  get_ok = m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry", fGeo);
  if (!get_ok) { std::cout << " P2RunQualityCheck tool: RUH ROH! Could not load AnnieGeometry!" << std::endl; return false; }

  get_ok = m_data->CStore.Get("ChannelNumToTankPMTSPEChargeMap", map_chankey2spe);
  if (!get_ok) { Log("P2RunQualityCheck tool: Could not find ChannelNumToTankPMTSPEChargeMap!", v_error, verbosity); }

  // Initialize ROOT stuff
  std::string root_outpfile_ext = ".root";
  std::string root_outpfile_name = outpfile_prefix + root_outpfile_ext;
  if (verbosity > 2) { std::cout << " P2RunQualityCheck tool: root file will be saved as " << root_outpfile_name << std::endl; }

  p2rqc_root_outp = new TFile(root_outpfile_name.c_str(), "RECREATE");

  this->InitTree();
  this->InitHist();
  this->InitGraph();

  std::cout << " P2RunQualityCheck tool: Initalization complete." << std::endl;
  if (verbosity > 4) { std::cout << "   Moving on to world destru--I mean, run quality check........" << std::endl; }

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
  globalEntry++;


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
/*  if ((fPartNumber > prevpart) && (fPartNumber != -999))
  {
    //std::cout << "  [[ DEBUG ]] did this werk??" << std::endl;
    //std::cout << "  [[ DEBUG ]] e = " << e << std::endl;
    e = 0;  //reset entry/evt counter
    prevpart = fPartNumber; //set to current partnumber
  } */

  //trigger type
  bool isBeam = false;
  bool isExtTrig = false;
  bool isCCTrig = false;
  bool isNonCCTrig = false;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("TriggerWord", trigword);
  if (!get_ok) { Log("P2RunQualityCheck tool: No TriggerWord object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("TriggerExtended", trigext);
  if (!get_ok) { Log("P2RunQualityCheck tool: No TriggerExtended object in ANNIEEvent! Abort!", v_error, verbosity); trigword = -99;/* return true;*/ }
  if (trigword == 5) { isBeam = true; beamEntry++; }
  if (trigext == 0) { isExtTrig = false; }
  if (trigext == 1) { isExtTrig = true; isCCTrig = true; }
  if (trigext == 2) { isExtTrig = true; isNonCCTrig = true; }

  //***** Beam info *****//
  BeamStatus beamstat;

  //reset some variables
  fBeamPOT = -8888;
  fBeamOK = -1;

  get_ok = m_data->Stores["ANNIEEvent"]->Get("BeamStatus", beamstat);
  if (get_ok)
  {
    //std::cout << "  [DEBUG] hello beam!" << std::endl;

    //all entries regardless of beam ok
    fBeamPOT = beamstat.pot();
    if (beamstat.ok()) { fBeamOK = 1; }
    else if (!beamstat.ok()) { fBeamOK = 0; }

    if (trigword == 5) h_pot_beam->Fill(beamstat.pot());
    if (trigword == 31) h_pot_led->Fill(beamstat.pot());
    if (trigword == 36) h_pot_cosmics->Fill(beamstat.pot());

    h2_pot_all->Fill(globalEntry, beamstat.pot());
    h2_pot_all_zoomx->Fill(globalEntry, beamstat.pot());
    h2_pot_all_zoomxy->Fill(globalEntry, beamstat.pot());

    if (beamstat.ok())
    {
      h_beampot_ok_all->Fill(beamstat.pot());
      if (trigword == 5) h_pot_ok_beam->Fill(beamstat.pot());
      if (trigword == 31) h_pot_ok_led->Fill(beamstat.pot());
      if (trigword == 36) h_pot_ok_cosmics->Fill(beamstat.pot());

      if (isBeam)
      {
        h2_pot_ok_beam->Fill(globalEntry, beamstat.pot());

        beampot_file << fEventTimeTank << " " << beamstat.pot() << std::endl;
      }
    }
    else  //!beamstat.ok()
    {
      std::cout << "  [DEBUG] not ok beam pot: " << beamstat.pot() << std::endl;

      h_beampot_nok_all->Fill(beamstat.pot());
      if (trigword == 5) h_pot_nok_beam->Fill(beamstat.pot());
      if (trigword == 31) h_pot_nok_led->Fill(beamstat.pot());
      if (trigword == 36) h_pot_nok_cosmics->Fill(beamstat.pot());
    }

    if (isBeam)
    {
      if (beamstat.pot() < 2.5e12) { n_beam_pot++; }
    }
  }
  else if (!get_ok) { Log("P2RunQualityCheck tool: BeamStatus object in ANNIEEvent! Abort!", v_error, verbosity); return false; }


  //***** Check for Veto hits *****//
/*  fVetoHit = 0;
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
*/


  //***** Check for MRD tracks *****//
/*  bool hasMRDTracks = false;
  int n_tracks_evt = -99;
  get_ok = m_data->Stores["MRDTracks"]->Get("NumMrdTracks", n_tracks_evt);
  if (!get_ok) { Log("P2RunQualityCheck tool: No NumMrdTracks object in ANNIEEvent! Abort!", v_error, verbosity); return true; }

  if (n_tracks_evt != 0)  //should this be > 0? can it ever be negative?
  {
    hasMRDTracks = true;
    if (verbosity > 3) { std::cout << " P2RunQualityCheck tool: NumMrdTracks = " << n_tracks_evt << std::endl; }
  }
*/


  t_tank->Fill();

  return true;
}


bool P2RunQualityCheck::Finalise(){

  beampot_file.close();

  p2rqc_root_outp->cd();
  t_tank->Write("", TObject::kOverwrite);
  this->WriteHist();
  this->WriteGraph();

  p2rqc_root_outp->Close();
  delete p2rqc_root_outp;

  std::cout << " P2RunQualityCheck tool: \n  Global number of entries = " << globalEntry << "\n   Number of ok evts = " << n_ok_evts << std::endl;
  std::cout << " P2RunQualityCheck tool: \n  Total number of beam entries = " << beamEntry << "\n   Number of beam entries with POT < 2.5e12 = " << n_beam_pot << std::endl;
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

  t_tank->Branch("beamok", &fBeamOK, "beamok/I");
  t_tank->Branch("beampot", &fBeamPOT, "beampot/D");


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

  h_beampot_ok_all = new TH1D("h_beampot_ok_all", "POT, beam_ok=true (all evts)", 1000, -10.e10, 10.e10);
  h_beampot_nok_all = new TH1D("h_beampot_nok_all", "POT, beam_ok=false (all evts)", 1000, -10.e10, 10.e10);
  h_pot_beam = new TH1D("h_pot_beam", "POT (trigword=5, ok&nok)", 5000, -10.e10, 5.5e12);
  h_pot_led = new TH1D("h_pot_led", "POT (trigword=31, ok&nok)", 5000, -10.e10, 5.5e12);
  h_pot_cosmics = new TH1D("h_pot_cosmics", "POT (trigword=36, ok&nok)", 5000, -10.e10, 5.5e12);
  h_pot_ok_beam = new TH1D("h_pot_ok_beam", "POT (trigword=5, beam_ok=true)", 5000, -10.e10, 5.5e12);
  h_pot_ok_led = new TH1D("h_pot_ok_led", "POT (trigword=31, beam_ok=true)", 5000, -10.e10, 5.5e12);
  h_pot_ok_cosmics = new TH1D("h_pot_ok_cosmics", "POT (trigword=36, beam_ok=true)", 5000, -10.e10, 5.5e12);
  h_pot_nok_beam = new TH1D("h_pot_nok_beam", "POT (trigword=5, beam_ok=false)", 5000, -10.e10, 5.5e12);
  h_pot_nok_led = new TH1D("h_pot_nok_led", "POT (trigword=31, beam_ok=false)", 5000, -10.e10, 5.5e12);
  h_pot_nok_cosmics = new TH1D("h_pot_nok_cosmics", "POT (trigword=36, beam_ok=false)", 5000, -10.e10, 5.5e12);

  h2_pot_all = new TH2D("h2_pot_all", "beam POT (all evts)", 80000, 0, 80000, 1000, 0, 5.3e12);
  h2_pot_all_zoomx = new TH2D("h2_pot_all_zoomx", "beam POT (all evts, zoomed in on x)", 1000, 0, 1000, 1000, 0, 5.3e12);
  h2_pot_all_zoomxy = new TH2D("h2_pot_all_zoomxy", "beam POT (all evts, zoomed in on x,y)", 1000, 0, 1000, 1000, 0, 1.0e12);
  h2_pot_ok_beam = new TH2D("h2_pot_ok_beam", "POT (trigword=5, beam_ok=true)", 50000, 0, 50000, 1000, 0, 5.3e12);

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

  h_beampot_ok_all->GetXaxis()->SetTitle("POT");
  h_beampot_ok_all->Write();
  h_beampot_nok_all->GetXaxis()->SetTitle("POT");
  h_beampot_nok_all->Write();
  h_pot_beam->GetXaxis()->SetTitle("POT");
  h_pot_beam->Write();
  h_pot_led->GetXaxis()->SetTitle("POT");
  h_pot_led->Write();
  h_pot_cosmics->GetXaxis()->SetTitle("POT");
  h_pot_cosmics->Write();
  h_pot_ok_beam->GetXaxis()->SetTitle("POT");
  h_pot_ok_beam->Write();
  h_pot_ok_led->GetXaxis()->SetTitle("POT");
  h_pot_ok_led->Write();
  h_pot_ok_cosmics->GetXaxis()->SetTitle("POT");
  h_pot_ok_cosmics->Write();
  h_pot_nok_beam->GetXaxis()->SetTitle("POT");
  h_pot_nok_beam->Write();
  h_pot_nok_led->GetXaxis()->SetTitle("POT");
  h_pot_nok_led->Write();
  h_pot_nok_cosmics->GetXaxis()->SetTitle("POT");
  h_pot_nok_cosmics->Write();

  // 2D histograms
  h2_pot_all->GetXaxis()->SetTitle("entry no.");
  h2_pot_all->GetYaxis()->SetTitle("POT");
  h2_pot_all->Write();
  h2_pot_all_zoomx->GetXaxis()->SetTitle("entry no.");
  h2_pot_all_zoomx->GetYaxis()->SetTitle("POT");
  h2_pot_all_zoomx->Write();
  h2_pot_all_zoomxy->GetXaxis()->SetTitle("entry no.");
  h2_pot_all_zoomxy->GetYaxis()->SetTitle("POT");
  h2_pot_all_zoomxy->Write();
  h2_pot_ok_beam->GetXaxis()->SetTitle("entry no.");
  h2_pot_ok_beam->GetYaxis()->SetTitle("POT");
  h2_pot_ok_beam->Write();

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

