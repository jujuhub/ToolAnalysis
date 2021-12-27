#include "RunValidation.h"

RunValidation::RunValidation():Tool(){}


bool RunValidation::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  createImages = false;

  m_variables.Get("verbosity",verbosity);
  m_variables.Get("OutputPath",outfile_path);
  m_variables.Get("InvertMRDTimes",invert_mrd_times); 
  m_variables.Get("RunNumber",user_runnumber);
  m_variables.Get("SubRunNumber",user_subrunnumber);
  m_variables.Get("RunType",user_runtype);
  m_variables.Get("FilenameSuffix",filename_suffix);
  //InvertMRDTimes should only be necessary if the times have not been saved correctly from the common stop mark in common stop mode
  m_variables.Get("SinglePEGains",singlePEgains);
  m_variables.Get("CreateImages",createImages);
  m_variables.Get("SavePath",savePath);

  if (createImages != 0 && createImages != 1){
    createImages = 0;
  }

  Log("RunValidation tool: Initialise",v_message,verbosity);

  //Get ANNIE Geometry
  m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);

  //Read in single p.e. gains //TODO: replace with global gain file in LoadGeometry
  ifstream file_singlepe(singlePEgains.c_str());
  unsigned long temp_chankey;
  double temp_gain;
  while (!file_singlepe.eof()){
    file_singlepe >> temp_chankey >> temp_gain;
    if (file_singlepe.eof()) break;
    pmt_gains.emplace(temp_chankey,temp_gain);
  }
  file_singlepe.close();

  //Initialize counting variables
  n_entries=0;
  n_pmt_clusters=0;
  n_pmt_clusters_threshold=0;
  n_mrd_clusters=0;
  n_pmt_mrd_clusters=0;
  n_pmt_mrd_nofacc=0;
  n_pmt_mrd_time=0;
  n_pmt_mrd_time_nofacc=0;
  n_pmt_mrd_time_facc=0;
  n_facc=0;
  n_facc_pmt=0;
  n_facc_mrd=0;
  n_facc_pmt_mrd=0;

  first_entry = true;

  return true;
}


bool RunValidation::Execute(){

  Log("RunValidation tool: Execute",v_message,verbosity);
  
  //-------------------------------------------------------------------------
  //Get necessary objects related to PMT clusters (from ClusterFinder tool)
  //-------------------------------------------------------------------------

  int get_ok;
  bool has_clusters = m_data->CStore.Get("ClusterMap",m_all_clusters);
  if (not has_clusters) { Log("RunValidation Tool: Error retrieving ClusterMap from CStore, did you run ClusterFinder beforehand?",v_error,verbosity);/* return false; */}
  has_clusters = m_data->CStore.Get("ClusterMapDetkey",m_all_clusters_detkey);
  if (not has_clusters) { Log("RunValidation Tool: Error retrieving ClusterMapDetkey from CStore, did you run ClusterFinder beforehand?",v_error,verbosity);/* return false;*/ }

  // Get ClusterClassifiers-related objects
  get_ok = m_data->Stores["ANNIEEvent"]->Get("ClusterChargeBalances", ClusterChargeBalances);
  if (not get_ok) { Log("RunValidation Tool: Error retrieving ClusterChargeBalances from ANNIEEvent, did you run ClusterClassifiers beforehand?", v_error, verbosity); /*return false;*/ }


  //-------------------------------------------------------------------------
  //Get necessary objects related to MRD clusters (from TimeClustering tool)
  //-------------------------------------------------------------------------

  get_ok = m_data->CStore.Get("MrdTimeClusters",MrdTimeClusters);
  if (not get_ok) { Log("RunValidation Tool: Error retrieving MrdTimeClusters map from CStore, did you run TimeClustering beforehand?",v_error,verbosity); return false; }
  if (MrdTimeClusters.size()!=0){
    get_ok = m_data->CStore.Get("MrdDigitTimes",MrdDigitTimes);
    if (not get_ok) { Log("RunValidation Tool: Error retrieving MrdDigitTimes map from CStore, did you run TimeClustering beforehand?",v_error,verbosity); return false; }
    get_ok = m_data->CStore.Get("MrdDigitChankeys",mrddigitchankeysthisevent);
    if (not get_ok) { Log("RunValidation Tool: Error retrieving MrdDigitChankeys, did you run TimeClustering beforehand",v_error,verbosity); return false;}
  }

  //-------------------------------------------------------------------------
  //--------------------- Get ANNIEEvent objects ----------------------------
  //-------------------------------------------------------------------------
   
  get_ok = m_data->Stores["ANNIEEvent"]->Get("RunNumber",RunNumber);
  if (!get_ok){
    Log("RunValidation tool: Did not find RunNumber in ANNIEEvent! Abort",v_error,verbosity);
    return false;
  }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("SubrunNumber",SubRunNumber);
  if (!get_ok){
    Log("RunValidation tool: Did not find SubrunNumber in ANNIEEvent! Abort",v_error,verbosity);
    return false;
  }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("RunType",RunType);
  if (!get_ok){
    Log("RunValidation tool: Did not find RunType in ANNIEEvent! Abort",v_error,verbosity);
    return false;
  }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("RunStartTime",RunStartTime);
  if (!get_ok){
    Log("RunValidation tool: Did not find RunStartTime in ANNIEEvent! Abort",v_error,verbosity);
    return false;
  }
  
  //Check if stored run information is sensible and replace with user input, if not
  if (RunNumber == -1) {
    RunNumber = user_runnumber;
    RunStartTime = 0;
  }
  if (SubRunNumber == -1) SubRunNumber = user_subrunnumber;
  if (RunType == -1) RunType = user_runtype;

  get_ok = m_data->Stores["ANNIEEvent"]->Get("EventNumber",EventNumber);
  if (!get_ok) {
    Log("RunValidation tool: Did not find EventNumber in ANNIEEvent! Abort",v_error,verbosity);
    return false;
  }
  EventTimeTank=0;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("EventTimeTank",EventTimeTank);
  if (!get_ok){
    Log("RunValidation tool: Did not find EventTimeTank in ANNIEEvent! Abort",v_error,verbosity);
  //  return false;
  }
  EventTimeTank/=1.e6;
  bool has_tdcdata = m_data->Stores["ANNIEEvent"]->Get("TDCData",TDCData);
  Log("RunValidation tool: RunNumber: "+std::to_string(RunNumber)+", SubrunNumber: "+std::to_string(SubRunNumber)+", RunType: "+std::to_string(RunType)+", RunStartTime: "+std::to_string(RunStartTime)+", EventNumber: "+std::to_string(EventNumber)+", EventTimeTank: "+std::to_string(EventTimeTank),v_message,verbosity);

  if (start_time == 0) start_time = EventTimeTank;

  std::string MrdTriggerType;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("MRDTriggerType",MrdTriggerType);
  if (!get_ok) {
    Log("RunValidation tool: Did not find MRDTriggerType in ANNIEEvent! Abort",v_error,verbosity);
    /*return false;*/	//Some events don't have MRD data -> ok
  }
  if (MrdTriggerType == "Cosmic") Log("RunValidation tool: Cosmic MRD trigger at event "+std::to_string(EventNumber),v_message,verbosity);

  int triggerword;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("TriggerWord",triggerword);
  if (!get_ok) {
    Log("RunValidation tool: Did not find TriggerWord in ANNIEEvent! Assume TriggerWord = 5",v_warning,verbosity);
    triggerword = 5;
  }

  int fExtended;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("TriggerExtended",fExtended);
  if (!get_ok){
    Log("RunValidation tool: Did not find TriggerExtended in ANNIEEvent! Assume TriggerExtended = 0",v_warning,verbosity);
    fExtended = 0;
  }

  BeamStatus beamstat;
  m_data->Stores["ANNIEEvent"]->Get("BeamStatus",beamstat);
  //std::cout <<"BeamStatus: "<<beamstat.Print()<<std::endl;
  double fPot = beamstat.pot();

  //-------------------------------------------------------------------------
  //--------------------- Initialise histograms------------------------------
  //-------------------------------------------------------------------------
  
  if (first_entry) this->DefineHistograms();
 
  Triggerwords->Fill(triggerword);
  
  // Define some booleans
  if (first_entry) current_time = EventTimeTank;
  else if (EventTimeTank > 0) current_time = EventTimeTank;

  bool mrd_hit = false;
  bool facc_hit = false;
  bool tank_hit = false;  
  bool coincident_tank_mrd = false;
  int n_delayed_clusters = 0;
  int n_delayed_clusters_CB = 0;
  bool coincident_tank_veto = false;
  bool coincident_mrd_veto = false;
  std::vector<double> delayed_clusters_times;
  std::vector<double> delayed_clusters_times_CB;  

  //-------------------------------------------------------------------------
  //--------------------- Fill size of raw waveforms-------------------------
  //-------------------------------------------------------------------------

  std::map<unsigned long, std::vector<Waveform<unsigned short>>> raw_waveform_map;
  bool has_raw = m_data->Stores["ANNIEEvent"]->Get("RawADCData",raw_waveform_map);
  /*if (!has_raw) {
    Log("RunValidation tool: Did not find RawADCData in ANNIEEvent! Abort",v_message,verbosity);
  }*/

  if (has_raw){
    for (auto& temp_pair : raw_waveform_map) {
      const auto& achannel_key = temp_pair.first;
      auto& araw_waveforms = temp_pair.second;
      for (unsigned int i=0; i< araw_waveforms.size(); i++){
        auto samples = araw_waveforms.at(i).GetSamples();
        int size_sample = samples->size();
        ADCWaveform_Samples->Fill(size_sample);
      }
    }
  } else {
    std::map<unsigned long, std::vector<int>> raw_acqsize_map;
    has_raw = m_data->Stores["ANNIEEvent"]->Get("RawAcqSize",raw_acqsize_map);
    int size_of_window = 2000;
    if (!has_raw) {
      Log("RunValidation tool: Did not find RawAcqSize in ANNIEEvent! Abort",v_warning,verbosity);
    } else {
      for (auto& temp_pair : raw_acqsize_map) {
        const auto& achannel_key = temp_pair.first;
        auto& araw_acqsize = temp_pair.second;
        for (unsigned int i=0; i< araw_acqsize.size(); i++){
          int size_sample = araw_acqsize.at(i);
          ADCWaveform_Samples->Fill(size_sample);
        }
      }
    }
  }

  //-------------------------------------------------------------------------
  //----------------Fill cluster time histograms-----------------------------
  //-------------------------------------------------------------------------

  double PMT_prompt_time=0.;
  int n_pmt_hits=0;
  double max_charge=0.;
  double max_charge_pe=0.;
  int clustersize = 0;
  if (has_clusters){
    clustersize = m_all_clusters->size();
    if (clustersize != 0){
      for(std::pair<double,std::vector<Hit>>&& apair : *m_all_clusters){
        double cluster_time = apair.first;
        std::vector<Hit>&Hits = apair.second;
        std::vector<unsigned long> detkeys = m_all_clusters_detkey->at(cluster_time);
        double charge_balance = ClusterChargeBalances.at(cluster_time);
        if (cluster_time > 12000){
          n_delayed_clusters++; //count delayed clusters in window 12us ... 67us
          delayed_clusters_times.push_back(cluster_time);
          PMT_DelayedTime->Fill(cluster_time); 
          if (charge_balance < 0.4) {
            n_delayed_clusters_CB++;
            delayed_clusters_times_CB.push_back(cluster_time);
            PMT_DelayedTime_CB->Fill(cluster_time);
          }
        }
        double global_time=0.;
        double global_charge=0.;
        double global_chargeperpmt=0.;
        int nhits=0;
        for (unsigned int i_hit = 0; i_hit < Hits.size(); i_hit++){
          unsigned long detkey = detkeys.at(i_hit);
          double time = Hits.at(i_hit).GetTime();
          double charge = Hits.at(i_hit).GetCharge()/pmt_gains[detkey];
          if (charge > 2) {
            if (triggerword == 5) PMT_t_clusters_2pe->Fill(time);
            else if (triggerword == 36) PMT_t_clusters_2pe_cosmics->Fill(time);
            else if (triggerword == 31) PMT_t_clusters_2pe_led->Fill(time);
   	    PMT_Channelkeys->Fill(detkey);	//Only fill channelkey histogram for hits above 2p.e.         
            if (time < 10080 || time > 10400) PMT_t_clusters_2pe_full->Fill(time);
            if (charge > 5) {
              if (triggerword == 5) PMT_t_clusters_5pe->Fill(time);
              else if (triggerword == 36) PMT_t_clusters_5pe_cosmics->Fill(time);
              else if (triggerword == 31) PMT_t_clusters_5pe_led->Fill(time);
              if (time < 10080 || time > 10400) PMT_t_clusters_5pe_full->Fill(time);
              if (charge > 10) {
                if (time < 10080 || time > 10400) PMT_t_clusters_10pe_full->Fill(time);
                if (triggerword == 5) PMT_t_clusters_10pe->Fill(time);
                else if (triggerword == 36) PMT_t_clusters_10pe_cosmics->Fill(time);
                else if (triggerword == 31) PMT_t_clusters_10pe_led->Fill(time);
                if (charge > 30) {
                  if (triggerword == 5) PMT_t_clusters_30pe->Fill(time);
                  else if (triggerword == 36) PMT_t_clusters_30pe_cosmics->Fill(time);
                  else if (triggerword == 31) PMT_t_clusters_30pe_led->Fill(time);
                  if (time < 10080 || time > 10400) PMT_t_clusters_30pe_full->Fill(time);
                }
              }
            }
	  }
	  global_time+=time;
          global_charge+=charge;
          nhits++;
          if (time < 10080 || time > 10400) PMT_t_clusters_full->Fill(time);
          if (time < 2000.) {
            if (triggerword == 5) PMT_t_clusters->Fill(time);
            else if (triggerword == 36) PMT_t_clusters_cosmics->Fill(time);
            else if (triggerword == 31) PMT_t_clusters_led->Fill(time);
          }
        }
	
        if (nhits>0) {
          global_time/=nhits;
          global_chargeperpmt=global_charge/double(nhits);
        }
        
        if (global_time < 2000.) {
          PMT_prompt_charge->Fill(global_charge);
          PMT_prompt_charge_zoom->Fill(global_charge);
          if (triggerword == 5) PMT_prompt_charge_Beam->Fill(global_charge);
          else if (triggerword == 31) PMT_prompt_charge_LED->Fill(global_charge);
          else if (triggerword == 36) PMT_prompt_charge_Cosmic->Fill(global_charge);
          if (nhits>=10) PMT_prompt_charge_10hits->Fill(global_charge);
          PMT_chargeperpmt->Fill(global_chargeperpmt);
          if (global_charge>100) PMT_chargeperpmt_100pe->Fill(global_chargeperpmt);
          if (global_charge > max_charge_pe) {
            max_charge_pe = global_charge;
            PMT_prompt_time = global_time;
          }
          PMT_prompt_charge_CB->Fill(global_charge,charge_balance);
          n_pmt_hits++;
          tank_hit = true;
        }
        else {
          if (global_time < 10080 || global_time > 10400) {
            PMT_delayed_charge->Fill(global_charge);
            PMT_delayed_charge_zoom->Fill(global_charge);
            if (nhits>=10) PMT_delayed_charge_10hits->Fill(global_charge);
            PMT_delayed_charge_CB->Fill(global_charge,charge_balance);
          }
	}
      }
    }
  }

  std::vector<double> mrd_cluster_times;
  double global_mrd_time=0;
  int n_mrd_hits=0;
  int max_mrd_hits = 0;
  bool cosmic_region = false;
  for (unsigned int i_cluster = 0; i_cluster < MrdTimeClusters.size(); i_cluster++){
    double temp_mrd_time = 0;
    int temp_mrd_hits = 0;
    std::vector<int> single_mrdcluster = MrdTimeClusters.at(i_cluster);
    int numdigits = single_mrdcluster.size();
    if (numdigits > 50) continue; 	//unphysical, noise event
    for (int thisdigit = 0; thisdigit < numdigits; thisdigit++){
      int digit_value = single_mrdcluster.at(thisdigit);
      unsigned long chankey = mrddigitchankeysthisevent.at(digit_value);
      Detector *thedetector = geom->ChannelToDetector(chankey);
      unsigned long detkey = thedetector->GetDetectorID();
      if (thedetector->GetDetectorElement()!="MRD") continue;
      else {
        MRD_Channelkeys->Fill(chankey);
        mrd_hit = true;
        int time = MrdDigitTimes.at(digit_value);
 	if (invert_mrd_times) time = 4000.-time;
        if (time > 1700 && time < 1750) cosmic_region = true;
        MRD_t_clusters->Fill(time);
        if (triggerword == 5) MRD_t_clusters_beam->Fill(time);
        if (triggerword == 36) MRD_t_clusters_cosmic->Fill(time);
        //MRD_t_clusters->Fill(time);
        temp_mrd_time+=time;
        temp_mrd_hits++;
      }
    }
    if (temp_mrd_hits > 0) temp_mrd_time /= temp_mrd_hits;
    if (temp_mrd_hits > max_mrd_hits){
      max_mrd_hits = temp_mrd_hits;
      global_mrd_time = temp_mrd_time;
      n_mrd_hits = temp_mrd_hits;
    }
    mrd_cluster_times.push_back(temp_mrd_time);
  }

  //-------------------------------------------------------------------------
  //---------------Loop over FACC entries------------------------------------
  //-------------------------------------------------------------------------
  
  if(has_tdcdata){
    if (TDCData->size()==0){
      Log("RunValidation tool: TDC data is empty in this event.",v_message,verbosity);
    } else {
      for (auto&& anmrdpmt : (*TDCData)){
        unsigned long chankey = anmrdpmt.first;
        FMV_Channelkeys->Fill(chankey);
        Detector* thedetector = geom->ChannelToDetector(chankey);
        unsigned long detkey = thedetector->GetDetectorID();
        if (thedetector->GetDetectorElement()=="Veto") {
          //facc_hit = true;
          std::vector<Hit> fmv_hits = anmrdpmt.second;
          for (int i_hit=0; i_hit < (int) fmv_hits.size(); i_hit++){
            Hit fmv_hit = fmv_hits.at(i_hit);
            double time_diff = fmv_hit.GetTime()-PMT_prompt_time;
            if (time_diff > (740) && time_diff < (840)){
              coincident_tank_veto = true;
            }
            FMV_PMT_Deltat->Fill(time_diff);
            if (max_charge_pe>100) FMV_PMT_Deltat_100pe->Fill(time_diff);
            FMV_t_clusters->Fill(fmv_hit.GetTime());
            if (PMT_prompt_time >0) FMV_PMT_t->Fill(fmv_hit.GetTime(),PMT_prompt_time);
            if (max_charge_pe > 100) FMV_PMT_t_100pe->Fill(fmv_hit.GetTime(),PMT_prompt_time);
            if (triggerword == 5) {
              FMV_t_clusters_beam->Fill(fmv_hit.GetTime());
              FMV_PMT_Deltat_Beam->Fill(time_diff);
              if (PMT_prompt_time > 0) FMV_PMT_t_Beam->Fill(fmv_hit.GetTime(),PMT_prompt_time);
              if (max_charge_pe > 100){
                FMV_PMT_Deltat_Beam_100pe->Fill(time_diff);
                FMV_PMT_t_100pe_Beam->Fill(fmv_hit.GetTime(),PMT_prompt_time);
              }
            }
            else if (triggerword == 36) {
              FMV_t_clusters_cosmic->Fill(fmv_hit.GetTime());
              FMV_PMT_Deltat_Cosmic->Fill(time_diff);
              if (PMT_prompt_time > 0) FMV_PMT_t_Cosmic->Fill(fmv_hit.GetTime(),PMT_prompt_time);
              if (max_charge_pe>100){
                FMV_PMT_Deltat_Cosmic_100pe->Fill(time_diff);
                FMV_PMT_t_100pe_Cosmic->Fill(fmv_hit.GetTime(),PMT_prompt_time);
              }
            }
            for (int i_mrd=0; i_mrd < (int) mrd_cluster_times.size(); i_mrd++){
              double diff = fmv_hit.GetTime()-mrd_cluster_times.at(i_mrd);
              MRD_FMV_Deltat->Fill(diff);
              if (triggerword == 5) MRD_FMV_Deltat_Beam->Fill(diff);
              else if (triggerword == 36) MRD_FMV_Deltat_Cosmic->Fill(diff);
              if (fabs(diff)<100) coincident_mrd_veto = true;
            }
          }
        }
      }
    }
  } else {
    Log("RunValidation tool: No TDC data available in this event.",v_message,verbosity);
  }
  
  if (coincident_tank_veto || coincident_mrd_veto) facc_hit = true;

  //-------------------------------------------------------------------------
  //---------------Increment event type counters-----------------------------
  //-------------------------------------------------------------------------
  
  bool pmt_mrd_time_coinc=false;
  if (n_mrd_hits>0) global_mrd_time/=n_mrd_hits;
  if (n_mrd_hits>0 && n_pmt_hits>0){
    for (int i_cluster = 0; i_cluster < (int) mrd_cluster_times.size(); i_cluster++){
    double i_mrd_time = mrd_cluster_times.at(i_cluster);
    MRD_PMT_t->Fill(i_mrd_time,PMT_prompt_time);
    MRD_PMT_Deltat->Fill(i_mrd_time-PMT_prompt_time);
    if (triggerword == 5){
      MRD_PMT_t_Beam->Fill(i_mrd_time,PMT_prompt_time);
      MRD_PMT_Deltat_Beam->Fill(i_mrd_time-PMT_prompt_time);
    } else if (triggerword == 36){
      MRD_PMT_t_Cosmic->Fill(i_mrd_time,PMT_prompt_time);
      MRD_PMT_Deltat_Cosmic->Fill(i_mrd_time-PMT_prompt_time);
    }
    if ((i_mrd_time-PMT_prompt_time)>700 && (i_mrd_time-PMT_prompt_time)<800) pmt_mrd_time_coinc = true;
    if (max_charge_pe > 100){
      MRD_PMT_t_100pe->Fill(i_mrd_time,PMT_prompt_time);
      MRD_PMT_Deltat_100pe->Fill(i_mrd_time-PMT_prompt_time);
      if (triggerword == 5){
        MRD_PMT_t_100pe_Beam->Fill(i_mrd_time,PMT_prompt_time);
        MRD_PMT_Deltat_100pe_Beam->Fill(i_mrd_time-PMT_prompt_time);
      } else if (triggerword == 36){
        MRD_PMT_t_100pe_Cosmic->Fill(i_mrd_time,PMT_prompt_time);
        MRD_PMT_Deltat_100pe_Cosmic->Fill(i_mrd_time-PMT_prompt_time);
      }
    }
    }
    coincident_tank_mrd = true;
  }

  if (clustersize > 0){ //Only count delayed clusters if there was a cluster at all
    PMT_DelayedMult->Fill(n_delayed_clusters);
  }

  //Increment event counters only for beam triggers

  if (triggerword == 5) {
    timestamps_beam.push_back(EventTimeTank);
    if (fExtended > 0) timestamps_beam_extended.push_back(EventTimeTank);
    if (fPot > 0) pot_values.push_back(fPot);
  } else if (triggerword == 31) {
    timestamps_led.push_back(EventTimeTank);
  } else if (triggerword == 36) {
    timestamps_cosmic.push_back(EventTimeTank);
  }

  if (triggerword == 5){
  if (tank_hit) n_pmt_clusters++;
  if (tank_hit && max_charge_pe > 100.) n_pmt_clusters_threshold++;	//100 p.e. threshold arbitrary at this point
  if (mrd_hit) n_mrd_clusters++;
  if (pmt_mrd_time_coinc) n_pmt_mrd_clusters++;
  if (pmt_mrd_time_coinc && !facc_hit) n_pmt_mrd_nofacc++;
  if (pmt_mrd_time_coinc) {
    n_pmt_mrd_time++;
    PMT_DelayedMult_Coinc->Fill(n_delayed_clusters);
    PMT_prompt_charge_MRDCoinc->Fill(max_charge_pe);
    for (int i_del = 0; i_del < (int) delayed_clusters_times.size(); i_del++){
      PMT_DelayedTime_Coinc->Fill(delayed_clusters_times.at(i_del));
    }
  }
  if (pmt_mrd_time_coinc && !facc_hit){
    PMT_DelayedMult_Coinc_NoFMV->Fill(n_delayed_clusters);
    PMT_DelayedMult_Coinc_NoFMV_CB->Fill(n_delayed_clusters_CB);
    PMT_prompt_charge_MRDCoinc_NoFMV->Fill(max_charge_pe);
    n_pmt_mrd_time_nofacc++;
    for (int i_del = 0; i_del < (int) delayed_clusters_times.size(); i_del++){
      PMT_DelayedTime_Coinc_NoFMV->Fill(delayed_clusters_times.at(i_del));
    }
    for (int i_del = 0; i_del < (int) delayed_clusters_times_CB.size(); i_del++){
      PMT_DelayedTime_Coinc_NoFMV_CB->Fill(delayed_clusters_times_CB.at(i_del));
    }
  }
  if (pmt_mrd_time_coinc && facc_hit) n_pmt_mrd_time_facc++;
  if (facc_hit) n_facc++;
  if (tank_hit && facc_hit) {
    n_facc_pmt++;
    PMT_prompt_charge_FMV->Fill(max_charge_pe);
  }
  if (mrd_hit && facc_hit) n_facc_mrd++;
  if (pmt_mrd_time_coinc && facc_hit) n_facc_pmt_mrd++;
  n_entries++;

  bool facc_mrd = (facc_hit && mrd_hit);
  if (coincident_tank_mrd){
    Log("RunValidation tool: Coincident tank + mrd clusters found. facc_hit && mrd_hit: "+std::to_string(facc_mrd)+", n_pmt_mrd_clusters: "+std::to_string(n_pmt_mrd_clusters)+", n_facc_mrd: "+std::to_string(n_facc_mrd),v_message,verbosity);
  }
  if (facc_hit && mrd_hit){
    Log("RunValidation tool: FACC hit found. coincident_tank_mrd: "+std::to_string(coincident_tank_mrd)+", n_pmt_mrd_clusters: "+std::to_string(n_pmt_mrd_clusters)+", n_facc_mrd: "+std::to_string(n_facc_mrd),v_message,verbosity);
  }
  }

  return true;
}


bool RunValidation::Finalise(){

  Log("RunValidation tool: Finalise",v_message,verbosity);
 
  //-------------------------------------------------------------------------
  //---------------------------Calculate rates-------------------------------
  //-------------------------------------------------------------------------

  double rate_entries=0;
  double rate_tank=0;
  double rate_tank_threshold=0;
  double rate_mrd=0;
  double rate_coincident=0;
  double rate_coincident_nofacc=0;
  double rate_coincident_time=0;
  double rate_coincident_time_nofacc=0;
  double rate_facc=0;
  double rate_facc_pmt=0;
  double rate_facc_mrd=0;
  double rate_facc_pmt_mrd=0;
  double fraction_entries=0;
  double fraction_tank=0;
  double fraction_tank_threshold=0;
  double fraction_mrd=0;
  double fraction_coincident=0;
  double fraction_coincident_nofacc=0;
  double fraction_facc=0;
  double fraction_facc_pmt=0;
  double fraction_facc_mrd=0;
  double fraction_facc_pmt_mrd=0;
  double fraction_coincident_time = 0;
  double fraction_coincident_time_nofacc = 0;
  double fraction_coincident_time_facc = 0;


  double time_diff = (current_time-start_time)/1000.;
  if (time_diff > 0.){
    rate_entries = n_entries/time_diff;
    rate_tank = n_pmt_clusters/time_diff;
    rate_tank_threshold = n_pmt_clusters_threshold/time_diff;
    rate_mrd = n_mrd_clusters/time_diff;
    rate_coincident = n_pmt_mrd_clusters/time_diff;
    rate_coincident_nofacc = n_pmt_mrd_nofacc/time_diff;
    rate_coincident_time = n_pmt_mrd_time/time_diff;
    rate_coincident_time_nofacc = n_pmt_mrd_time_nofacc/time_diff;
    rate_facc = n_facc/time_diff;
    rate_facc_pmt = n_facc_pmt/time_diff;
    rate_facc_mrd = n_facc_mrd/time_diff;
    rate_facc_pmt_mrd = n_facc_pmt_mrd/time_diff;
  }
  if (n_entries > 0){
    fraction_entries = 100;
    fraction_tank = double(n_pmt_clusters)/n_entries*100;
    fraction_tank_threshold = double(n_pmt_clusters_threshold)/n_entries*100;
    fraction_mrd = double(n_mrd_clusters)/n_entries*100;
    fraction_coincident = double(n_pmt_mrd_clusters)/n_entries*100;
    fraction_coincident_nofacc = double(n_pmt_mrd_nofacc)/n_entries*100;
    fraction_facc = double(n_facc)/n_entries*100;
    fraction_facc_pmt = double(n_facc_pmt)/n_entries*100;
    fraction_facc_mrd = double(n_facc_mrd)/n_entries*100;
    fraction_facc_pmt_mrd = double(n_facc_pmt_mrd)/n_entries*100;
    fraction_coincident_time = double(n_pmt_mrd_time)/n_entries*100;
    fraction_coincident_time_nofacc = double(n_pmt_mrd_time_nofacc)/n_entries*100;
  }

  ANNIE_rates->SetBinContent(1,rate_entries);
  ANNIE_rates->SetBinContent(2,rate_tank);
  ANNIE_rates->SetBinContent(3,rate_tank_threshold);
  ANNIE_rates->SetBinContent(4,rate_mrd);
  ANNIE_rates->SetBinContent(5,rate_coincident_time);
  ANNIE_rates->SetBinContent(6,rate_coincident_time_nofacc);
  ANNIE_rates->SetBinContent(7,rate_facc);
  ANNIE_rates->SetBinContent(8,rate_facc_pmt);
  ANNIE_rates->SetBinContent(9,rate_facc_mrd);
  ANNIE_rates->SetBinContent(10,rate_facc_pmt_mrd);
  
  ANNIE_counts->SetBinContent(1,n_entries);
  ANNIE_counts->SetBinContent(2,n_pmt_clusters);
  ANNIE_counts->SetBinContent(3,n_pmt_clusters_threshold);
  ANNIE_counts->SetBinContent(4,n_mrd_clusters);
  ANNIE_counts->SetBinContent(5,n_pmt_mrd_time);
  ANNIE_counts->SetBinContent(6,n_pmt_mrd_time_nofacc);
  ANNIE_counts->SetBinContent(7,n_facc);
  ANNIE_counts->SetBinContent(8,n_facc_pmt);
  ANNIE_counts->SetBinContent(9,n_facc_mrd);
  ANNIE_counts->SetBinContent(10,n_facc_pmt_mrd);
  
  ANNIE_fractions->SetBinContent(1,fraction_entries);
  ANNIE_fractions->SetBinContent(2,fraction_tank);
  ANNIE_fractions->SetBinContent(3,fraction_tank_threshold);
  ANNIE_fractions->SetBinContent(4,fraction_mrd);
  ANNIE_fractions->SetBinContent(5,fraction_coincident_time);
  ANNIE_fractions->SetBinContent(6,fraction_coincident_time_nofacc);
  ANNIE_fractions->SetBinContent(7,fraction_facc);
  ANNIE_fractions->SetBinContent(8,fraction_facc_pmt);
  ANNIE_fractions->SetBinContent(9,fraction_facc_mrd);
  ANNIE_fractions->SetBinContent(10,fraction_facc_pmt_mrd);

  const char *category_label[10] = {"All Events","PMT Clusters","PMT Clusters > 100p.e.","MRD Clusters","PMT+MRD Clusters","PMT+MRD, No FMV","FMV","FMV+PMT","FMV+MRD","FMV+PMT+MRD"};
  for (int i_label=0;i_label<10;i_label++){
    ANNIE_rates->GetXaxis()->SetBinLabel(i_label+1,category_label[i_label]);
    ANNIE_counts->GetXaxis()->SetBinLabel(i_label+1,category_label[i_label]);
    ANNIE_fractions->GetXaxis()->SetBinLabel(i_label+1,category_label[i_label]);
  }

  std::stringstream title_trigger_beam, title_trigger_cosmic, title_trigger_led, title_trigger_beam_extended;
  title_trigger_beam << "Trigger Rates Beam - Run "<<GlobalRunNumber;
  title_trigger_beam_extended << "Trigger Rates Beam Extended - Run "<<GlobalRunNumber;
  title_trigger_cosmic << "Trigger Rates Cosmic - Run "<<GlobalRunNumber;
  title_trigger_led << "Trigger Rates LED - Run "<<GlobalRunNumber;

  TriggerRate_Beam = new TH1D("TriggerRate_Beam",title_trigger_beam.str().c_str(),100,start_time/1000.,current_time/(1000.));
  TriggerRate_Beam_Extended = new TH1D("TriggerRate_Beam_Extended",title_trigger_beam_extended.str().c_str(),100,start_time/1000.,current_time/(1000.));
  TriggerRate_Cosmic = new TH1D("TriggerRate_Cosmic",title_trigger_cosmic.str().c_str(),100,start_time/1000.,current_time/(1000.));
  TriggerRate_LED = new TH1D("TriggerRate_LED",title_trigger_led.str().c_str(),100,start_time/1000.,current_time/(1000.));

  std::cout <<"start_time: "<<start_time<<", current_time: "<<current_time<<std::endl;

  for (int i_trigger = 0; i_trigger < (int) timestamps_beam.size(); i_trigger++){
    TriggerRate_Beam->Fill(timestamps_beam.at(i_trigger)/1000.);
    if (i_trigger == 0) std::cout <<"Beam timestamp: "<<timestamps_beam.at(i_trigger)/1000.<<std::endl;
    else if (i_trigger == (timestamps_beam.size()-1)) std::cout <<"Last Beam timestamp: "<<timestamps_beam.at(i_trigger)/1000.<<std::endl;
  }
  for (int i_trigger = 0; i_trigger < (int) timestamps_beam_extended.size(); i_trigger++){
    TriggerRate_Beam_Extended->Fill(timestamps_beam_extended.at(i_trigger)/1000.);
  }
  for (int i_trigger = 0; i_trigger < (int) timestamps_cosmic.size(); i_trigger++){
    TriggerRate_Cosmic->Fill(timestamps_cosmic.at(i_trigger)/1000.);
  }
  for (int i_trigger = 0; i_trigger < (int) timestamps_led.size(); i_trigger++){
    TriggerRate_LED->Fill(timestamps_led.at(i_trigger)/1000.);
  }

  double scale_factor = 100./time_diff;
  TriggerRate_Beam->Scale(scale_factor);
  TriggerRate_Beam_Extended->Scale(scale_factor);
  TriggerRate_Cosmic->Scale(scale_factor);
  TriggerRate_LED->Scale(scale_factor);

  double rate_trigger_beam = 0;
  double rate_trigger_beam_extended = 0;
  double rate_trigger_cosmic = 0;
  double rate_trigger_led = 0;
  for (int i_bin = 0; i_bin < (int) TriggerRate_Beam->GetXaxis()->GetNbins(); i_bin++){
    rate_trigger_beam += TriggerRate_Beam->GetBinContent(i_bin+1);
    rate_trigger_beam_extended += TriggerRate_Beam_Extended->GetBinContent(i_bin+1);
    rate_trigger_cosmic += TriggerRate_Cosmic->GetBinContent(i_bin+1);
    rate_trigger_led += TriggerRate_LED->GetBinContent(i_bin+1);
  }
  rate_trigger_beam /= 100.;
  rate_trigger_beam_extended /= 100.;
  rate_trigger_cosmic /= 100.;
  rate_trigger_led /= 100.;

  double pot_total=0;
  for (int i_pot=0; i_pot < pot_values.size(); i_pot++){
    pot_total += pot_values.at(i_pot);
  }
  pot_total /= time_diff;

  TriggerRate_Beam->GetXaxis()->SetTimeDisplay(1);
  TriggerRate_Beam->GetXaxis()->SetLabelSize(0.03);
  TriggerRate_Beam->GetXaxis()->SetLabelOffset(0.03);
  TriggerRate_Beam->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  TriggerRate_Beam->GetYaxis()->SetTitleSize(0.045);
  TriggerRate_Beam->GetYaxis()->SetTitle("R_{beam} [Hz]");
  TriggerRate_Beam_Extended->GetXaxis()->SetTimeDisplay(1);
  TriggerRate_Beam_Extended->GetXaxis()->SetLabelSize(0.03);
  TriggerRate_Beam_Extended->GetXaxis()->SetLabelOffset(0.03);
  TriggerRate_Beam_Extended->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  TriggerRate_Beam->GetYaxis()->SetTitleSize(0.045);
  TriggerRate_Beam->GetYaxis()->SetTitle("R_{beam,extended} [Hz]");
  TriggerRate_Cosmic->GetXaxis()->SetTimeDisplay(1);
  TriggerRate_Cosmic->GetXaxis()->SetLabelSize(0.03);
  TriggerRate_Cosmic->GetXaxis()->SetLabelOffset(0.03);
  TriggerRate_Cosmic->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  TriggerRate_Cosmic->GetYaxis()->SetTitleSize(0.045);
  TriggerRate_Cosmic->GetYaxis()->SetTitle("R_{cosmic} [Hz]");
  TriggerRate_LED->GetXaxis()->SetTimeDisplay(1);
  TriggerRate_LED->GetXaxis()->SetLabelSize(0.03);
  TriggerRate_LED->GetXaxis()->SetLabelOffset(0.03);
  TriggerRate_LED->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
  TriggerRate_LED->GetYaxis()->SetTitleSize(0.045);
  TriggerRate_LED->GetYaxis()->SetTitle("R_{LED} [Hz]");

  tree_stats = new TTree("tree_stats","Run Stats");
  tree_stats->Branch("start_time",&start_time);
  tree_stats->Branch("end_time",&current_time);
  tree_stats->Branch("rate_trigger_beam",&rate_trigger_beam);
  tree_stats->Branch("rate_trigger_beam_extended",&rate_trigger_beam_extended);
  tree_stats->Branch("rate_trigger_cosmic",&rate_trigger_cosmic);
  tree_stats->Branch("rate_trigger_led",&rate_trigger_led);
  tree_stats->Branch("rate_entries",&rate_entries);
  tree_stats->Branch("rate_tank",&rate_tank);
  tree_stats->Branch("rate_tank_threshold",&rate_tank_threshold);
  tree_stats->Branch("rate_mrd",&rate_mrd);
  tree_stats->Branch("rate_coincident_time",&rate_coincident_time);
  tree_stats->Branch("rate_coincident_time_nofacc",&rate_coincident_time_nofacc);
  tree_stats->Branch("rate_facc",&rate_facc);
  tree_stats->Branch("rate_facc_pmt",&rate_facc_pmt);
  tree_stats->Branch("rate_facc_mrd",&rate_facc_mrd);
  tree_stats->Branch("rate_facc_pmt_mrd",&rate_facc_pmt_mrd);
  tree_stats->Branch("pot_total",&pot_total);
  tree_stats->Fill();

  Log("RunValidation tool: Finished analysing run. Summary:",v_warning,verbosity);
  Log("RunValidation tool: Duration: "+std::to_string(int(time_diff/3600))+"h:"+std::to_string(int(time_diff/60)%60)+"min:"+std::to_string(int(time_diff)%60)+"sec",v_warning,verbosity);
  Log("RunValidation tool: Num Entries: "+std::to_string(n_entries)+", PMT Clusters: "+std::to_string(n_pmt_clusters)+", PMT Clusters > 100p.e.: "+std::to_string(n_pmt_clusters_threshold)+", MRD Clusters: "+std::to_string(n_mrd_clusters)+", PMT+MRD clusters: "+std::to_string(n_pmt_mrd_clusters)+", FACC Events: "+std::to_string(n_facc)+", FACC+PMT: "+std::to_string(n_facc_pmt)+", FACC+MRD: "+std::to_string(n_facc_mrd)+", FACC+PMT+MRD: "+std::to_string(n_facc_pmt_mrd),v_warning,verbosity); 
  Log("RunValidation tool: Total Rate: "+std::to_string(rate_entries)+", Rate PMT Clusters: "+std::to_string(rate_tank)+", Rate PMT Clusters > 100p.e.: "+std::to_string(rate_tank_threshold)+", Rate MRD Clusters: "+std::to_string(rate_mrd)+", Rate PMT+MRD clusters: "+std::to_string(rate_coincident)+", FACC Events: "+std::to_string(rate_facc)+", Rate FACC+PMT: "+std::to_string(rate_facc_pmt)+", Rate FACC+MRD: "+std::to_string(rate_facc_mrd)+", Rate FACC+PMT+MRD: "+std::to_string(rate_facc_pmt_mrd),v_warning,verbosity); 
  Log("RunValidation tool: Fraction PMT Clusters: "+std::to_string(fraction_tank)+", Fraction PMT Clusters > 100p.e.: "+std::to_string(fraction_tank_threshold)+", Fraction MRD Clusters: "+std::to_string(fraction_mrd)+", Fraction PMT+MRD clusters: "+std::to_string(fraction_coincident)+", Fraction FACC: "+std::to_string(fraction_facc)+", Fraction FACC+PMT: "+std::to_string(fraction_facc_pmt)+", Fraction FACC+MRD: "+std::to_string(fraction_facc_mrd)+", Fraction FACC+PMT+MRD: "+std::to_string(fraction_facc_pmt_mrd),v_warning,verbosity); 
  Log("RunValidation tool: Fraction Coincident PMT/MRD clusters: "+std::to_string(fraction_coincident_time)+", Fraction Coincident PMT/MRD clusters with FMV hit: "+std::to_string(fraction_facc_pmt_mrd),v_warning,verbosity);
  Log("RunValidation tool: End of Summary",v_warning,verbosity);


  //-------------------------------------------------------------------------
  //--------------------Write histograms to file-----------------------------
  //-------------------------------------------------------------------------
  
  outfile->cd();
  tree_stats->Write();
  MRD_t_clusters->Write();
  MRD_t_clusters_beam->Write();
  MRD_t_clusters_cosmic->Write();
  FMV_t_clusters->Write();
  FMV_t_clusters_beam->Write();
  FMV_t_clusters_cosmic->Write();
  PMT_t_clusters->Write();
  PMT_t_clusters_2pe->Write();
  PMT_t_clusters_5pe->Write();
  PMT_t_clusters_10pe->Write();
  PMT_t_clusters_30pe->Write();
  PMT_t_clusters_cosmics->Write();
  PMT_t_clusters_2pe_cosmics->Write();
  PMT_t_clusters_5pe_cosmics->Write();
  PMT_t_clusters_10pe_cosmics->Write();
  PMT_t_clusters_30pe_cosmics->Write();
  PMT_t_clusters_led->Write();
  PMT_t_clusters_2pe_led->Write();
  PMT_t_clusters_5pe_led->Write();
  PMT_t_clusters_10pe_led->Write();
  PMT_t_clusters_30pe_led->Write();
  PMT_t_clusters_full->Write();
  PMT_t_clusters_2pe_full->Write();
  PMT_t_clusters_5pe_full->Write();
  PMT_t_clusters_10pe_full->Write();
  PMT_t_clusters_30pe_full->Write();
  MRD_PMT_t->Write();
  MRD_PMT_t_100pe->Write();
  MRD_PMT_t_Beam->Write();
  MRD_PMT_t_100pe_Beam->Write();
  MRD_PMT_t_Cosmic->Write();
  MRD_PMT_t_100pe_Cosmic->Write();
  FMV_PMT_t->Write();
  FMV_PMT_t_100pe->Write();
  FMV_PMT_t_Beam->Write();
  FMV_PMT_t_100pe_Beam->Write();
  FMV_PMT_t_Cosmic->Write();
  FMV_PMT_t_100pe_Cosmic->Write();
  MRD_PMT_Deltat->Write();
  MRD_PMT_Deltat_100pe->Write();
  MRD_PMT_Deltat_Beam->Write();
  MRD_PMT_Deltat_100pe_Beam->Write();
  MRD_PMT_Deltat_Cosmic->Write();
  MRD_PMT_Deltat_100pe_Cosmic->Write();
  FMV_PMT_Deltat->Write();
  FMV_PMT_Deltat_100pe->Write();
  FMV_PMT_Deltat_Beam->Write();
  FMV_PMT_Deltat_Beam_100pe->Write();
  FMV_PMT_Deltat_Cosmic->Write();
  FMV_PMT_Deltat_Cosmic_100pe->Write();
  MRD_FMV_Deltat->Write();
  MRD_FMV_Deltat_Beam->Write();
  MRD_FMV_Deltat_Cosmic->Write();
  PMT_DelayedMult->Write();
  PMT_DelayedMult_Coinc->Write();
  PMT_DelayedMult_Coinc_NoFMV->Write();
  PMT_DelayedMult_Coinc_NoFMV_CB->Write();
  PMT_DelayedTime->Write();
  PMT_DelayedTime_CB->Write();
  PMT_DelayedTime_Coinc->Write();
  PMT_DelayedTime_Coinc_NoFMV->Write();
  PMT_DelayedTime_Coinc_NoFMV_CB->Write();
  PMT_prompt_charge->Write();
  PMT_prompt_charge_zoom->Write();
  PMT_prompt_charge_10hits->Write();
  PMT_prompt_charge_MRDCoinc->Write();
  PMT_prompt_charge_MRDCoinc_NoFMV->Write();
  PMT_prompt_charge_FMV->Write();
  PMT_prompt_charge_Beam->Write();
  PMT_prompt_charge_Cosmic->Write();
  PMT_prompt_charge_LED->Write();
  PMT_prompt_charge_CB->Write();
  PMT_delayed_charge->Write();
  PMT_delayed_charge_zoom->Write();
  PMT_delayed_charge_10hits->Write();
  PMT_delayed_charge_CB->Write();
  PMT_chargeperpmt->Write();
  PMT_chargeperpmt_100pe->Write();
  ADCWaveform_Samples->Write();
  Triggerwords->Write();
  PMT_Channelkeys->Write();
  MRD_Channelkeys->Write();
  FMV_Channelkeys->Write();
  ANNIE_counts->Write();
  ANNIE_rates->Write();
  ANNIE_fractions->Write();
  TriggerRate_Beam->Write();
  TriggerRate_Beam_Extended->Write();
  TriggerRate_Cosmic->Write();
  TriggerRate_LED->Write();

  canvas_beamspill->Divide(2,1);
  canvas_beamspill->cd(1);
  PMT_t_clusters_2pe->Rebin(5);
  PMT_t_clusters_2pe->Draw();
  canvas_beamspill->cd(2);
  MRD_t_clusters_beam->Rebin(5);
  MRD_t_clusters_beam->Draw();
  canvas_beamspill->Write();

  canvas_mrd_pmt->Divide(2,2);
  canvas_mrd_pmt->cd(1);
  MRD_PMT_t_100pe->Draw("colz");
  canvas_mrd_pmt->cd(2);
  MRD_PMT_Deltat_100pe->Draw();
  canvas_mrd_pmt->cd(3);
  FMV_PMT_Deltat_100pe->Draw();  
  canvas_mrd_pmt->cd(4);
  MRD_FMV_Deltat->Draw();
  canvas_mrd_pmt->Write();

  canvas_triggerwords->Divide(2,1);
  canvas_triggerwords->cd(1);
  Triggerwords->Draw();
  TPad *p = (TPad*) canvas_triggerwords->cd(2);
  p->SetLogy();
  ADCWaveform_Samples->Draw();  
  canvas_triggerwords->Write();

  canvas_prompt_charge->Divide(2,2);
  TPad *pq1 = (TPad*) canvas_prompt_charge->cd(1);
  pq1->SetLogy();
  PMT_prompt_charge->Draw();
  TPad *pq2 = (TPad*) canvas_prompt_charge->cd(2);
  pq2->SetLogy();
  PMT_prompt_charge_MRDCoinc_NoFMV->Draw();
  TPad *pq3 = (TPad*) canvas_prompt_charge->cd(3);
  pq3->SetLogy();
  PMT_delayed_charge_zoom->Draw();
  canvas_prompt_charge->cd(4);
  PMT_delayed_charge_CB->Draw("colz");
  canvas_prompt_charge->Write();

  canvas_neutron_mult->Divide(2,2);
  TPad *pn1 = (TPad*) canvas_neutron_mult->cd(1);
  pn1->SetLogy();
  PMT_DelayedMult->Draw();
  TPad *pn2 = (TPad*) canvas_neutron_mult->cd(2);
  pn2->SetLogy();
  PMT_DelayedMult_Coinc_NoFMV_CB->Draw();
  TPad *pn3 = (TPad*) canvas_neutron_mult->cd(3);
  pn3->SetLogy();
  PMT_DelayedTime->Draw();
  canvas_neutron_mult->cd(4);
  PMT_DelayedTime_Coinc_NoFMV_CB->Draw("colz");
  canvas_neutron_mult->Write();

  canvas_channelkey->Divide(3,1);
  canvas_channelkey->SetLogy();
  canvas_channelkey->cd(1);
  PMT_Channelkeys->Draw();
  canvas_channelkey->cd(2);
  MRD_Channelkeys->Draw();
  canvas_channelkey->cd(3);
  FMV_Channelkeys->Draw();
  canvas_channelkey->Write();

  canvas_rates->Divide(2,1);
  canvas_rates->SetLogy();
  TPad *pr1 = (TPad*) canvas_rates->cd(1);
  pr1->SetLogy();
  pr1->SetGridy();
  ANNIE_rates->Draw();
  TPad *pr2 = (TPad*) canvas_rates->cd(2); 
  pr2->SetLogy();
  pr2->SetGridy();
  ANNIE_fractions->Draw();
  canvas_rates->Write();

  if (createImages){

    std::stringstream save_mrd_t, save_mrd_t_beam, save_mrd_t_cosmic, save_fmv_t, save_fmv_t_beam, save_fmv_t_cosmic, save_pmt_t, save_pmt_t_2pe, save_pmt_t_5pe, save_pmt_t_10pe, save_pmt_t_30pe, save_pmt_tfull, save_pmt_tfull_2pe, save_pmt_tfull_5pe, save_pmt_tfull_10pe, save_pmt_tfull_30pe;
    std::stringstream save_mrd_pmt_t, save_mrd_pmt_t_100pe, save_mrd_pmt_deltat, save_mrd_pmt_deltat_100pe, save_mrd_pmt_t_cosmic, save_mrd_pmt_t_100pe_cosmic, save_mrd_pmt_t_beam, save_mrd_pmt_t_100pe_beam, save_mrd_pmt_deltat_cosmic, save_mrd_pmt_deltat_100pe_cosmic, save_mrd_pmt_deltat_beam, save_mrd_pmt_deltat_100pe_beam, save_fmv_pmt_deltat, save_fmv_pmt_deltat_100pe, save_fmv_pmt_deltat_beam, save_fmv_pmt_deltat_cosmic, save_fmv_pmt_deltat_beam_100pe, save_fmv_pmt_deltat_cosmic_100pe, save_mrd_fmv_deltat, save_mrd_fmv_deltat_beam, save_mrd_fmv_deltat_cosmic;
    std::stringstream save_pmt_prompt, save_pmt_prompt_zoom, save_pmt_prompt_10, save_pmt_prompt_mrdcoinc, save_pmt_prompt_mrdcoinc_nofmv, save_pmt_prompt_fmv, save_pmt_prompt_beam, save_pmt_prompt_cosmic, save_pmt_prompt_led, save_pmt_prompt_cb, save_pmt_del, save_pmt_del_zoom, save_pmt_del_10, save_pmt_del_cb, save_pmt_avgq, save_pmt_avgq_100pe;
    std::stringstream save_pmt_chkey, save_mrd_chkey, save_fmv_chkey, save_annie_counts, save_annie_rates, save_annie_fractions;
    std::stringstream save_del_mult, save_del_mult_coinc, save_del_mult_coinc_nofmv, save_del_mult_coinc_nofmv_cb;
    std::stringstream save_adcwaveform, save_triggerwords;
    std::stringstream save_fmv_pmt_t, save_fmv_pmt_t_100pe, save_fmv_pmt_t_cosmic, save_fmv_pmt_t_100pe_cosmic, save_fmv_pmt_t_beam, save_fmv_pmt_t_100pe_beam;
    std::stringstream save_trigger_beam, save_trigger_beam_extended, save_trigger_cosmic, save_trigger_led;
    std::stringstream save_canvas_beamspill, save_canvas_mrd_pmt, save_canvas_triggerwords, save_canvas_prompt_charge, save_canvas_neutron_mult, save_canvas_channelkey, save_canvas_rates;

    save_mrd_t << savePath << "R" << GlobalRunNumber << "_MRDTimes.png";
    save_mrd_t_beam << savePath << "R" << GlobalRunNumber << "_MRDTimesBeam.png";
    save_mrd_t_cosmic << savePath << "R" << GlobalRunNumber << "_MRDTimesCosmic.png";
    save_fmv_t << savePath << "R" << GlobalRunNumber << "_FMVTimes.png";
    save_fmv_t_beam << savePath << "R" << GlobalRunNumber << "_FMVTimesBeam.png";
    save_fmv_t_cosmic << savePath << "R" << GlobalRunNumber << "_FMVTimesCosmic.png";
    save_pmt_t << savePath << "R" << GlobalRunNumber << "_PMTTimes.png";
    save_pmt_t_2pe << savePath << "R" << GlobalRunNumber << "_PMTTimes_2pe.png";
    save_pmt_t_5pe << savePath << "R" << GlobalRunNumber << "_PMTTimes_5pe.png";
    save_pmt_t_10pe << savePath << "R" << GlobalRunNumber << "_PMTTimes_10pe.png";
    save_pmt_t_30pe << savePath << "R" << GlobalRunNumber << "_PMTTimes_30pe.png";
    save_pmt_tfull << savePath << "R" << GlobalRunNumber << "_PMTTimes_FullTimeWindow.png";
    save_pmt_tfull_2pe << savePath << "R" << GlobalRunNumber << "_PMTTimes_FullWindow_2pe.png";
    save_pmt_tfull_5pe << savePath << "R" << GlobalRunNumber << "_PMTTimes_FullWindow_5pe.png";
    save_pmt_tfull_10pe << savePath << "R" << GlobalRunNumber << "_PMTTimes_FullWindow_10pe.png";
    save_pmt_tfull_30pe << savePath << "R" << GlobalRunNumber << "_PMTTimes_FullWindow_30pe.png";
    save_mrd_pmt_t << savePath << "R" << GlobalRunNumber << "_MRDPMT_Times.png";
    save_mrd_pmt_t_100pe << savePath << "R" << GlobalRunNumber << "_MRDPMT_Times_100pe.png";
    save_mrd_pmt_t_cosmic << savePath << "R" << GlobalRunNumber << "_MRDPMT_Times_Cosmic.png";
    save_mrd_pmt_t_100pe_cosmic << savePath << "R" << GlobalRunNumber << "_MRDPMT_Times_100pe_Cosmic.png";
    save_mrd_pmt_t_beam << savePath << "R" << GlobalRunNumber << "_MRDPMT_Times_Beam.png";
    save_mrd_pmt_t_100pe_beam << savePath << "R" << GlobalRunNumber << "_MRDPMT_Times_100pe_Beam.png";
    save_fmv_pmt_t << savePath << "R" << GlobalRunNumber << "_FMVPMT_Times.png";
    save_fmv_pmt_t_100pe << savePath << "R" << GlobalRunNumber << "_FMVPMT_Times_100pe.png";
    save_fmv_pmt_t_cosmic << savePath << "R" << GlobalRunNumber << "_FMVPMT_Times_Cosmic.png";
    save_fmv_pmt_t_100pe_cosmic << savePath << "R" << GlobalRunNumber << "_FMVPMT_Times_100pe_Cosmic.png";
    save_fmv_pmt_t_beam << savePath << "R" << GlobalRunNumber << "_FMVPMT_Times_Beam.png";
    save_fmv_pmt_t_100pe_beam << savePath << "R" << GlobalRunNumber << "_FMVPMT_Times_100pe_Beam.png";
    save_mrd_pmt_deltat << savePath << "R" << GlobalRunNumber << "_MRDPMT_DeltaTimes.png";
    save_mrd_pmt_deltat_100pe << savePath << "R" << GlobalRunNumber << "_MRDPMT_DeltaTimes_100pe.png";
    save_mrd_pmt_deltat_cosmic << savePath << "R" << GlobalRunNumber << "_MRDPMT_DeltaTimes_Cosmic.png";
    save_mrd_pmt_deltat_100pe_cosmic << savePath << "R" << GlobalRunNumber << "_MRDPMT_DeltaTimes_100pe_Cosmic.png";
    save_mrd_pmt_deltat_beam << savePath << "R" << GlobalRunNumber << "_MRDPMT_DeltaTimes_Beam.png";
    save_mrd_pmt_deltat_100pe_beam << savePath << "R" << GlobalRunNumber << "_MRDPMT_DeltaTimes_100pe_Beam.png";
    save_fmv_pmt_deltat << savePath << "R" << GlobalRunNumber << "_FMVPMT_DeltaTimes.png";
    save_fmv_pmt_deltat_100pe << savePath << "R" << GlobalRunNumber << "_FMVPMT_DeltaTimes_100pe.png";
    save_fmv_pmt_deltat_beam << savePath << "R" << GlobalRunNumber << "_FMVPMT_DeltaTimes_Beam.png";
    save_fmv_pmt_deltat_beam_100pe << savePath << "R" << GlobalRunNumber << "_FMVPMT_DeltaTimes_Beam_100pe.png";
    save_fmv_pmt_deltat_cosmic << savePath << "R" << GlobalRunNumber << "_FMVPMT_DeltaTimes_Cosmic.png";
    save_fmv_pmt_deltat_cosmic_100pe << savePath << "R" << GlobalRunNumber << "_FMVPMT_DeltaTimes_Cosmic_100pe.png";
    save_mrd_fmv_deltat << savePath << "R" << GlobalRunNumber << "_MRDFMV_DeltaTimes.png";
    save_mrd_fmv_deltat_beam << savePath << "R" << GlobalRunNumber << "_MRDFMV_DeltaTimes_Beam.png";
    save_mrd_fmv_deltat_cosmic << savePath << "R" << GlobalRunNumber << "_MRDFMV_DeltaTimes_Cosmic.png";
    save_pmt_prompt << savePath << "R" << GlobalRunNumber << "_PMTCharge_Prompt.png";
    save_pmt_prompt_zoom << savePath << "R" << GlobalRunNumber << "_PMTCharge_Prompt_zoom.png";
    save_pmt_prompt_10 << savePath << "R" << GlobalRunNumber << "_PMTCharge_Prompt_10hits.png";
    save_pmt_prompt_mrdcoinc << savePath << "R" << GlobalRunNumber << "_PMTCharge_Prompt_MRDCoinc.png";
    save_pmt_prompt_mrdcoinc_nofmv << savePath << "R" << GlobalRunNumber << "_PMTCharge_Prompt_MRDCoinc_NoFMV.png";
    save_pmt_prompt_fmv << savePath << "R" << GlobalRunNumber << "_PMTCharge_Prompt_FMV.png";
    save_pmt_prompt_beam << savePath << "R" << GlobalRunNumber << "_PMTCharge_Prompt_Beam.png";
    save_pmt_prompt_cosmic << savePath << "R" << GlobalRunNumber << "_PMTCharge_Prompt_Cosmic.png";
    save_pmt_prompt_led << savePath << "R" << GlobalRunNumber << "_PMTCharge_Prompt_LED.png";
    save_pmt_prompt_cb << savePath << "R" << GlobalRunNumber << "_PMTCharge_Prompt_CB.png";
    save_pmt_del << savePath << "R" << GlobalRunNumber << "_PMTCharge_Delayed.png";
    save_pmt_del_10 << savePath << "R" << GlobalRunNumber << "_PMTCharge_Delayed_10hits.png";
    save_pmt_del_cb << savePath << "R" << GlobalRunNumber << "_PMTCharge_Delayed_CB.png";
    save_pmt_del_zoom << savePath << "R" << GlobalRunNumber << "_PMTCharge_Delayed_zoom.png";
    save_pmt_avgq << savePath << "R" << GlobalRunNumber << "_PMTCharge_AvgQ.png";
    save_pmt_avgq_100pe << savePath << "R" << GlobalRunNumber << "_PMTCharge_AvgQ_100pe.png";
    save_pmt_chkey << savePath << "R" << GlobalRunNumber << "_PMTChannelkeys.png";
    save_mrd_chkey << savePath << "R" << GlobalRunNumber << "_MRDChannelkeys.png";
    save_fmv_chkey << savePath << "R" << GlobalRunNumber << "_FMVChannelkeys.png";
    save_annie_counts << savePath << "R" << GlobalRunNumber << "_ANNIECounts.png";
    save_annie_rates << savePath << "R" << GlobalRunNumber << "_ANNIERates.png";
    save_annie_fractions << savePath << "R" << GlobalRunNumber << "_ANNIEFractions.png";
    save_del_mult << savePath << "R" << GlobalRunNumber << "_DelayedMult.png";
    save_del_mult_coinc << savePath << "R" << GlobalRunNumber << "_DelayedMult_Coinc.png";
    save_del_mult_coinc_nofmv << savePath << "R" << GlobalRunNumber << "_DelayedMult_Coinc_NoFMV.png";
    save_del_mult_coinc_nofmv_cb << savePath << "R" << GlobalRunNumber << "_DelayedMult_Coinc_NoFMV_CB.png";
    save_adcwaveform << savePath << "R" << GlobalRunNumber << "_ADCWaveform.png";
    save_triggerwords << savePath << "R" << GlobalRunNumber << "_Triggerwords.png";
    save_trigger_beam << savePath << "R" << GlobalRunNumber << "_TriggerRatesBeam.png";
    save_trigger_beam_extended << savePath << "R" << GlobalRunNumber << "_TriggerRatesBeamExtended.png";
    save_trigger_cosmic << savePath << "R" << GlobalRunNumber << "_TriggerRatesCosmic.png";
    save_trigger_led << savePath << "R" << GlobalRunNumber << "_TriggerRatesLED.png";
    save_canvas_beamspill << savePath << "R" << GlobalRunNumber << "_CanvasBeamSpill.png";
    save_canvas_mrd_pmt << savePath << "R" << GlobalRunNumber << "_CanvasMRDPMT.png";
    save_canvas_triggerwords << savePath << "R" << GlobalRunNumber << "_CanvasTriggerWords.png";
    save_canvas_prompt_charge << savePath << "R" << GlobalRunNumber << "_CanvasPromptCharge.png";
    save_canvas_neutron_mult << savePath << "R" << GlobalRunNumber << "_CanvasNeutronMult.png";
    save_canvas_channelkey << savePath << "R" << GlobalRunNumber << "_CanvasChannelkey.png";
    save_canvas_rates << savePath << "R" << GlobalRunNumber << "_CanvasRates.png";


    TCanvas *canvas = new TCanvas("canvas","Canvas",900,600);
    canvas->cd();
    MRD_t_clusters->Draw();
    canvas->SaveAs(save_mrd_t.str().c_str());
    canvas->Clear();
    MRD_t_clusters_beam->Draw();
    canvas->SaveAs(save_mrd_t_beam.str().c_str());
    canvas->Clear();
    MRD_t_clusters_cosmic->Draw();
    canvas->SaveAs(save_mrd_t_cosmic.str().c_str());
    canvas->Clear();
    FMV_t_clusters->Draw();
    canvas->SaveAs(save_fmv_t.str().c_str());
    canvas->Clear();
    FMV_t_clusters_beam->Draw();
    canvas->SaveAs(save_fmv_t_beam.str().c_str());
    canvas->Clear();
    FMV_t_clusters_cosmic->Draw();
    canvas->SaveAs(save_fmv_t_cosmic.str().c_str());
    canvas->Clear();
    PMT_t_clusters->Draw();
    canvas->SaveAs(save_pmt_t.str().c_str());
    canvas->Clear();
    PMT_t_clusters_2pe->Draw();
    canvas->SaveAs(save_pmt_t_2pe.str().c_str());
    canvas->Clear();
    PMT_t_clusters_5pe->Draw();
    canvas->SaveAs(save_pmt_t_5pe.str().c_str());
    canvas->Clear();
    PMT_t_clusters_10pe->Draw();
    canvas->SaveAs(save_pmt_t_10pe.str().c_str());
    canvas->Clear();
    PMT_t_clusters_30pe->Draw();
    canvas->SaveAs(save_pmt_t_30pe.str().c_str());
    canvas->Clear();
    PMT_t_clusters_full->Draw();
    canvas->SaveAs(save_pmt_tfull.str().c_str());
    canvas->Clear();
    PMT_t_clusters_2pe_full->Draw();
    canvas->SaveAs(save_pmt_tfull_2pe.str().c_str());
    canvas->Clear();
    PMT_t_clusters_5pe_full->Draw();
    canvas->SaveAs(save_pmt_tfull_5pe.str().c_str());
    canvas->Clear();
    PMT_t_clusters_10pe_full->Draw();
    canvas->SaveAs(save_pmt_tfull_10pe.str().c_str());
    canvas->Clear();
    PMT_t_clusters_30pe_full->Draw();
    canvas->SaveAs(save_pmt_tfull_30pe.str().c_str());
    canvas->Clear();
    MRD_PMT_t->Draw();
    canvas->SaveAs(save_mrd_pmt_t.str().c_str());
    canvas->Clear();
    MRD_PMT_t_100pe->Draw();
    canvas->SaveAs(save_mrd_pmt_t_100pe.str().c_str());
    canvas->Clear();
    MRD_PMT_t_Beam->Draw();
    canvas->SaveAs(save_mrd_pmt_t_beam.str().c_str());
    canvas->Clear();
    MRD_PMT_t_100pe_Beam->Draw();
    canvas->SaveAs(save_mrd_pmt_t_100pe_beam.str().c_str());
    canvas->Clear();
    MRD_PMT_t_Cosmic->Draw();
    canvas->SaveAs(save_mrd_pmt_t_cosmic.str().c_str());
    canvas->Clear();
    MRD_PMT_t_100pe_Cosmic->Draw();
    canvas->SaveAs(save_mrd_pmt_t_100pe_cosmic.str().c_str());
    canvas->Clear();
    FMV_PMT_t->Draw();
    canvas->SaveAs(save_fmv_pmt_t.str().c_str());
    canvas->Clear();
    FMV_PMT_t_100pe->Draw();
    canvas->SaveAs(save_fmv_pmt_t_100pe.str().c_str());
    canvas->Clear();
    FMV_PMT_t_Beam->Draw();
    canvas->SaveAs(save_fmv_pmt_t_beam.str().c_str());
    canvas->Clear();
    FMV_PMT_t_100pe_Beam->Draw();
    canvas->SaveAs(save_fmv_pmt_t_100pe_beam.str().c_str());
    canvas->Clear();
    FMV_PMT_t_Cosmic->Draw();
    canvas->SaveAs(save_fmv_pmt_t_cosmic.str().c_str());
    canvas->Clear();
    FMV_PMT_t_100pe_Cosmic->Draw();
    canvas->SaveAs(save_fmv_pmt_t_100pe_cosmic.str().c_str());
    canvas->Clear();
    MRD_PMT_Deltat->Draw();
    canvas->SaveAs(save_mrd_pmt_deltat.str().c_str());
    canvas->Clear();
    MRD_PMT_Deltat_100pe->Draw();
    canvas->SaveAs(save_mrd_pmt_deltat_100pe.str().c_str());
    canvas->Clear();
    MRD_PMT_Deltat_Cosmic->Draw();
    canvas->SaveAs(save_mrd_pmt_deltat_cosmic.str().c_str());
    canvas->Clear();
    MRD_PMT_Deltat_100pe_Cosmic->Draw();
    canvas->SaveAs(save_mrd_pmt_deltat_100pe_cosmic.str().c_str());
    canvas->Clear();
    MRD_PMT_Deltat_Beam->Draw();
    canvas->SaveAs(save_mrd_pmt_deltat_beam.str().c_str());
    canvas->Clear();
    MRD_PMT_Deltat_100pe_Beam->Draw();
    canvas->SaveAs(save_mrd_pmt_deltat_100pe_beam.str().c_str());
    canvas->Clear();
    FMV_PMT_Deltat->Draw();
    canvas->SaveAs(save_fmv_pmt_deltat.str().c_str());
    canvas->Clear();
    FMV_PMT_Deltat_100pe->Draw();
    canvas->SaveAs(save_fmv_pmt_deltat_100pe.str().c_str());
    canvas->Clear();
    FMV_PMT_Deltat_Beam->Draw();
    canvas->SaveAs(save_fmv_pmt_deltat_beam.str().c_str());
    canvas->Clear();
    FMV_PMT_Deltat_Beam_100pe->Draw();
    canvas->SaveAs(save_fmv_pmt_deltat_beam_100pe.str().c_str());
    canvas->Clear();
    FMV_PMT_Deltat_Cosmic->Draw();
    canvas->SaveAs(save_fmv_pmt_deltat_cosmic.str().c_str());
    canvas->Clear();
    FMV_PMT_Deltat_Cosmic_100pe->Draw();
    canvas->SaveAs(save_fmv_pmt_deltat_cosmic_100pe.str().c_str());
    canvas->Clear();
    MRD_FMV_Deltat->Draw();
    canvas->SaveAs(save_mrd_fmv_deltat.str().c_str());
    canvas->Clear();
    PMT_DelayedMult->Draw();
    canvas->SaveAs(save_del_mult.str().c_str());
    canvas->Clear();
    PMT_DelayedMult_Coinc->Draw();
    canvas->SaveAs(save_del_mult_coinc.str().c_str());
    canvas->Clear();
    PMT_DelayedMult_Coinc_NoFMV->Draw();
    canvas->SaveAs(save_del_mult_coinc_nofmv.str().c_str());
    canvas->Clear();
    PMT_DelayedMult_Coinc_NoFMV_CB->Draw();
    canvas->SaveAs(save_del_mult_coinc_nofmv_cb.str().c_str());
    canvas->Clear();
    PMT_prompt_charge->Draw();
    canvas->SaveAs(save_pmt_prompt.str().c_str());
    canvas->Clear();
    PMT_prompt_charge_zoom->Draw();
    canvas->SaveAs(save_pmt_prompt_zoom.str().c_str());
    canvas->Clear();
    PMT_prompt_charge_10hits->Draw();
    canvas->SaveAs(save_pmt_prompt_10.str().c_str());
    canvas->Clear();
    PMT_prompt_charge_MRDCoinc->Draw();
    canvas->SaveAs(save_pmt_prompt_mrdcoinc.str().c_str());
    canvas->Clear();
    PMT_prompt_charge_MRDCoinc_NoFMV->Draw();
    canvas->SaveAs(save_pmt_prompt_mrdcoinc_nofmv.str().c_str());
    canvas->Clear();
    PMT_prompt_charge_FMV->Draw();
    canvas->SaveAs(save_pmt_prompt_fmv.str().c_str());
    canvas->Clear();
    PMT_prompt_charge_Beam->Draw();
    canvas->SaveAs(save_pmt_prompt_beam.str().c_str());
    canvas->Clear();
    PMT_prompt_charge_Cosmic->Draw();
    canvas->SaveAs(save_pmt_prompt_cosmic.str().c_str());
    canvas->Clear();
    PMT_prompt_charge_LED->Draw();
    canvas->SaveAs(save_pmt_prompt_led.str().c_str());
    canvas->Clear();
    PMT_prompt_charge_CB->Draw();
    canvas->SaveAs(save_pmt_prompt_cb.str().c_str());
    canvas->Clear();
    PMT_delayed_charge->Draw();
    canvas->SaveAs(save_pmt_del.str().c_str());
    canvas->Clear();    
    PMT_delayed_charge_zoom->Draw();
    canvas->SaveAs(save_pmt_del_zoom.str().c_str());
    canvas->Clear();    
    PMT_delayed_charge_10hits->Draw();
    canvas->SaveAs(save_pmt_del_10.str().c_str());
    canvas->Clear();    
    PMT_delayed_charge_CB->Draw();
    canvas->SaveAs(save_pmt_del_cb.str().c_str());
    canvas->Clear();    
    PMT_chargeperpmt->Draw();
    canvas->SaveAs(save_pmt_avgq.str().c_str());
    canvas->Clear();    
    PMT_chargeperpmt_100pe->Draw();
    canvas->SaveAs(save_pmt_avgq_100pe.str().c_str());
    canvas->Clear();    
    ANNIE_counts->Draw();
    canvas->SaveAs(save_annie_counts.str().c_str());
    canvas->Clear();    
    ANNIE_rates->Draw();
    canvas->SaveAs(save_annie_rates.str().c_str());
    canvas->Clear();    
    ANNIE_fractions->Draw();
    canvas->SaveAs(save_annie_fractions.str().c_str());
    canvas->Clear();    
    ADCWaveform_Samples->Draw();
    canvas->SaveAs(save_adcwaveform.str().c_str());
    canvas->Clear();    
    PMT_Channelkeys->Draw();
    canvas->SaveAs(save_pmt_chkey.str().c_str());
    canvas->Clear();
    MRD_Channelkeys->Draw();
    canvas->SaveAs(save_mrd_chkey.str().c_str());
    canvas->Clear();
    FMV_Channelkeys->Draw();
    canvas->SaveAs(save_fmv_chkey.str().c_str());
    canvas->Clear();
    Triggerwords->Draw();
    canvas->SaveAs(save_triggerwords.str().c_str());
    canvas->Clear();
    TriggerRate_Beam->Draw();
    canvas->SaveAs(save_trigger_beam.str().c_str());
    canvas->Clear();
    TriggerRate_Beam_Extended->Draw();
    canvas->SaveAs(save_trigger_beam_extended.str().c_str());
    canvas->Clear();
    TriggerRate_Cosmic->Draw();
    canvas->SaveAs(save_trigger_cosmic.str().c_str());
    canvas->Clear();
    TriggerRate_LED->Draw();
    canvas->SaveAs(save_trigger_led.str().c_str());
    canvas->Clear();
    canvas_beamspill->SaveAs(save_canvas_beamspill.str().c_str());
    canvas_mrd_pmt->SaveAs(save_canvas_mrd_pmt.str().c_str());
    canvas_triggerwords->SaveAs(save_canvas_triggerwords.str().c_str());
    canvas_prompt_charge->SaveAs(save_canvas_prompt_charge.str().c_str());
    canvas_neutron_mult->SaveAs(save_canvas_neutron_mult.str().c_str());
    canvas_channelkey->SaveAs(save_canvas_channelkey.str().c_str());
    canvas_rates->SaveAs(save_canvas_rates.str().c_str());
  }

  outfile->Close();
  delete outfile;

  return true;
} 

void RunValidation::DefineHistograms(){
    
    GlobalRunNumber = RunNumber;
    GlobalSubRunNumber = SubRunNumber;
    start_time = EventTimeTank;

    std::stringstream file_name;
    file_name << outfile_path << "RunValidation_R"<<GlobalRunNumber<<"S"<<GlobalSubRunNumber<<"T"<<RunType<<".root";
    outfile = new TFile(file_name.str().c_str(),"RECREATE");
    outfile->cd();

    std::stringstream title_mrd, title_mrd_beam, title_mrd_cosmic, title_fmv, title_fmv_beam, title_fmv_cosmic, title_pmt, title_pmt_2pe, title_pmt_5pe, title_pmt_10pe, title_pmt_30pe, title_mrd_pmt, title_mrd_pmt_100pe, title_mrd_pmt_delta, title_mrd_pmt_delta_100pe, title_mrd_pmt_cosmic, title_mrd_pmt_100pe_cosmic, title_mrd_pmt_delta_cosmic, title_mrd_pmt_delta_100pe_cosmic, title_fmv_pmt_delta, title_fmv_pmt_delta_100pe, title_mrd_fmv_delta, title_pmt_prompt, title_pmt_prompt_10, title_pmt_prompt_MRDCoinc, title_pmt_prompt_MRDCoinc_NoFMV, title_pmt_prompt_FMV, title_pmt_prompt_beam, title_pmt_prompt_cosmic, title_pmt_prompt_led, title_pmt_prompt_CB, title_chargeperpmt, title_chargeperpmt_100pe, title_pmt_delayed, title_pmt_delayed_10, title_pmt_delayed_CB, title_counts, title_rates, title_fractions, title_delayed_mult, title_delayed_mult_coinc, title_delayed_mult_coinc_nofmv, title_delayed_mult_coinc_nofmv_cb, title_delayed_time, title_delayed_time_cb, title_delayed_time_coinc, title_delayed_time_coinc_nofmv, title_delayed_time_coinc_nofmv_cb, title_adcwaveform_samples, title_pmt_chankeys, title_mrd_chankeys, title_fmv_chankeys, title_triggerwords;

    std::stringstream title_pmt_led, title_pmt_2pe_led, title_pmt_5pe_led, title_pmt_10pe_led, title_pmt_30pe_led;
    std::stringstream title_pmt_cosmics, title_pmt_2pe_cosmics, title_pmt_5pe_cosmics, title_pmt_10pe_cosmics, title_pmt_30pe_cosmics;

    std::stringstream title_mrd_pmt_beam, title_mrd_pmt_100pe_beam, title_mrd_pmt_delta_beam, title_mrd_pmt_delta_100pe_beam;
    std::stringstream title_fmv_pmt_delta_beam, title_fmv_pmt_delta_beam_100pe, title_fmv_pmt_delta_cosmic, title_fmv_pmt_delta_cosmic_100pe, title_mrd_fmv_delta_beam, title_mrd_fmv_delta_cosmic;
    std::stringstream title_fmv_pmt, title_fmv_pmt_100pe, title_fmv_pmt_beam, title_fmv_pmt_100pe_beam, title_fmv_pmt_cosmic, title_fmv_pmt_100pe_cosmic;

    title_mrd << "MRD Cluster Times - Run "<<GlobalRunNumber;
    title_mrd_beam << "MRD Cluster Times (Beam) - Run "<<GlobalRunNumber;
    title_mrd_cosmic << "MRD Cluster Times (Cosmic) - Run "<<GlobalRunNumber;
    title_fmv << "FMV Cluster Times - Run "<<GlobalRunNumber;
    title_fmv_beam << "FMV Cluster Times (Beam) - Run "<<GlobalRunNumber;
    title_fmv_cosmic << "FMV Cluster Times (Cosmic) - Run "<<GlobalRunNumber;
    title_pmt << "PMT Cluster Times - Run "<<GlobalRunNumber;
    title_pmt_2pe << "PMT Cluster Times (>2 p.e.) - Run "<<GlobalRunNumber;
    title_pmt_5pe << "PMT Cluster Times (>5 p.e.) - Run "<<GlobalRunNumber;
    title_pmt_10pe << "PMT Cluster Times (>10 p.e.) - Run "<<GlobalRunNumber;
    title_pmt_30pe << "PMT Cluster Times (>30 p.e.) - Run "<<GlobalRunNumber;
    title_pmt_led << "PMT Cluster Times - LED triggerword - Run "<<GlobalRunNumber;
    title_pmt_2pe_led << "PMT Cluster Times - LED triggerword (>2 p.e.) - Run "<<GlobalRunNumber;
    title_pmt_5pe_led << "PMT Cluster Times - LED triggerword (>5 p.e.) - Run "<<GlobalRunNumber;
    title_pmt_10pe_led << "PMT Cluster Times - LED triggerword (>10 p.e.) - Run "<<GlobalRunNumber;
    title_pmt_30pe_led << "PMT Cluster Times - LED triggerword (>30 p.e.) - Run "<<GlobalRunNumber;
    title_pmt_cosmics << "PMT Cluster Times - Cosmic triggerword - Run "<<GlobalRunNumber;
    title_pmt_2pe_cosmics << "PMT Cluster Times - Cosmic triggerword (>2 p.e.) - Run "<<GlobalRunNumber;
    title_pmt_5pe_cosmics << "PMT Cluster Times - Cosmic triggerword (>5 p.e.) - Run "<<GlobalRunNumber;
    title_pmt_10pe_cosmics << "PMT Cluster Times - Cosmic triggerword (>10 p.e.) - Run "<<GlobalRunNumber;
    title_pmt_30pe_cosmics << "PMT Cluster Times - Cosmic triggerword (>30 p.e.) - Run "<<GlobalRunNumber;
    title_mrd_pmt << "MRD vs. PMT Cluster Times - Run "<<GlobalRunNumber;
    title_mrd_pmt_100pe << "MRD vs. PMT Cluster Times (>100 p.e.) - Run "<<GlobalRunNumber;
    title_mrd_pmt_beam << "MRD vs. PMT Cluster Times Beam - Run "<<GlobalRunNumber;
    title_mrd_pmt_100pe_beam << "MRD vs. PMT Cluster Times Beam (>100 p.e.) - Run "<<GlobalRunNumber;
    title_mrd_pmt_delta << "Difference MRD/PMT Cluster Times - Run "<<GlobalRunNumber;
    title_mrd_pmt_delta_100pe << "Difference MRD/PMT Cluster Times (>100 p.e.) - Run "<<GlobalRunNumber;
    title_mrd_pmt_cosmic << "MRD vs. PMT Cluster Times Cosmics - Run "<<GlobalRunNumber;
    title_mrd_pmt_100pe_cosmic << "MRD vs. PMT Cluster Times Cosmics (>100 p.e.) - Run "<<GlobalRunNumber;
    title_mrd_pmt_delta_cosmic << "Difference MRD/PMT Cluster Times Cosmics - Run "<<GlobalRunNumber;
    title_mrd_pmt_delta_100pe_cosmic << "Difference MRD/PMT Cluster Times Cosmics (>100 p.e.) - Run "<<GlobalRunNumber;
    title_fmv_pmt << "FMV vs. PMT Cluster Times - Run "<<GlobalRunNumber;
    title_fmv_pmt_100pe << "FMV vs. PMT Cluster Times (>100 p.e.) - Run "<<GlobalRunNumber;
    title_fmv_pmt_beam << "FMV vs. PMT Cluster Times - Beam - Run "<<GlobalRunNumber;
    title_fmv_pmt_100pe_beam << "FMV vs. PMT Cluster Times - Beam (>100 p.e.) - Run "<<GlobalRunNumber;
    title_fmv_pmt_cosmic << "FMV vs. PMT Cluster Times - Cosmic - Run "<<GlobalRunNumber;
    title_fmv_pmt_100pe_cosmic << "FMV vs. PMT Cluster Times - Cosmic (>100 p.e.) - Run "<<GlobalRunNumber;
    title_fmv_pmt_delta << "Difference FMV/PMT Cluster Times - Run "<<GlobalRunNumber;
    title_fmv_pmt_delta_100pe << "Difference FMV/PMT Cluster Times (>100 p.e.) - Run "<<GlobalRunNumber;
    title_fmv_pmt_delta_beam << "Difference FMV/PMT Cluster Times - Beam - Run "<<GlobalRunNumber;
    title_fmv_pmt_delta_beam_100pe << "Difference FMV/PMT Cluster Times - Beam (>100 p.e.) - Run "<<GlobalRunNumber;
    title_fmv_pmt_delta_cosmic << "Difference FMV/PMT Cluster Times - Cosmic - Run "<<GlobalRunNumber;
    title_fmv_pmt_delta_cosmic_100pe << "Difference FMV/PMT Cluster Times - Cosmic (>100 p.e.) - Run "<<GlobalRunNumber;
    title_mrd_fmv_delta << "Difference MRD/FMV Cluster Times - Run "<<GlobalRunNumber;
    title_mrd_fmv_delta_beam << "Difference MRD/FMV Cluster Times - Beam - Run "<<GlobalRunNumber;
    title_mrd_fmv_delta_cosmic << "Difference MRD/FMV Cluster Times - Cosmic - Run "<<GlobalRunNumber;
    title_pmt_prompt << "PMT Prompt Charge - Run "<<GlobalRunNumber;
    title_pmt_prompt_10 << "PMT Prompt Charge (>10 hits) - Run "<<GlobalRunNumber;
    title_pmt_prompt_MRDCoinc << "PMT Prompt Charge (MRD Coinc) - Run "<<GlobalRunNumber;
    title_pmt_prompt_MRDCoinc_NoFMV << "PMT Prompt Charge (MRD Coinc, No FMV) - Run "<<GlobalRunNumber;
    title_pmt_prompt_FMV << "PMT Prompt Charge (FMV hit) - Run "<<GlobalRunNumber;
    title_pmt_prompt_beam << "PMT Prompt Charge (Beam trigger) - Run "<<GlobalRunNumber;
    title_pmt_prompt_cosmic << "PMT Prompt Charge (Cosmic trigger) - Run "<<GlobalRunNumber;
    title_pmt_prompt_led << "PMT Prompt Charge (LED trigger) - Run "<<GlobalRunNumber;
    title_pmt_prompt_CB << "PMT Charge Balance vs prompt Q - Run "<<GlobalRunNumber;
    title_chargeperpmt << "PMT Charge/PMT - Run "<<GlobalRunNumber;
    title_chargeperpmt_100pe << "PMT Charge/PMT (>100 p.e.) - Run "<<GlobalRunNumber;
    title_pmt_delayed << "PMT Delayed Charge - Run "<<GlobalRunNumber;
    title_pmt_delayed_10 << "PMT Delayed Charge (>10 hits)- Run "<<GlobalRunNumber;
    title_pmt_delayed_CB << "PMT Charge Balance vs delayed Q - Run "<<GlobalRunNumber;
    title_pmt_chankeys << "PMT Channelkeys - Run "<<GlobalRunNumber;
    title_mrd_chankeys << "MRD Channelkeys - Run "<<GlobalRunNumber;
    title_fmv_chankeys << "FMV Channelkeys - Run "<<GlobalRunNumber;
    title_counts << "ANNIE Beam Trigger Counts - Run "<<GlobalRunNumber;
    title_rates << " ANNIE Beam Trigger Rates - Run "<<GlobalRunNumber;
    title_fractions << " ANNIE Beam Trigger Event Fractions - Run "<<GlobalRunNumber;
    title_delayed_mult << "Neutron Candidate Multiplicity Distribution - Run "<<GlobalRunNumber;
    title_delayed_mult_coinc << "Neutron Candidate Multiplicity Distribution (MRD/Tank coincidence) - Run "<<GlobalRunNumber;
    title_delayed_mult_coinc_nofmv << "Neutron Candidate Multiplicity Distribution (MRD/Tank coincidence, No FMV) - Run "<<GlobalRunNumber;
    title_delayed_mult_coinc_nofmv_cb << "Neutron Candidate Multiplicity Distribution (MRD/Tank coincidence, No FMV, CB < 0.4) - Run "<<GlobalRunNumber;
    title_delayed_time << "Neutron Candidate Time Distribution - Run "<<GlobalRunNumber;
    title_delayed_time_cb << "Neutron Candidate Time Distribution - Run "<<GlobalRunNumber;
    title_delayed_time_coinc << "Neutron Candidate Time Distribution (MRD/Tank coincidence) - Run "<<GlobalRunNumber;
    title_delayed_time_coinc_nofmv << "Neutron Candidate Time Distribution (MRD/Tank coincidence, No FMV) - Run "<<GlobalRunNumber;
    title_delayed_time_coinc_nofmv_cb << "Neutron Candidate Time Distribution (MRD/Tank coincidence, No FMV, CB < 0.4) - Run "<<GlobalRunNumber;
    title_adcwaveform_samples << "Number of samples in ADC waveforms - Run "<<GlobalRunNumber;
    title_triggerwords << "Triggerwords - Run "<<GlobalRunNumber;   

    MRD_t_clusters = new TH1D("MRD_t_clusters",title_mrd.str().c_str(),250,0,4000);
    MRD_t_clusters_beam = new TH1D("MRD_t_clusters_beam",title_mrd_beam.str().c_str(),250,0,4000);
    MRD_t_clusters_cosmic = new TH1D("MRD_t_clusters_cosmic",title_mrd_cosmic.str().c_str(),250,0,4000);
    FMV_t_clusters = new TH1D("FMV_t_clusters",title_fmv.str().c_str(),250,0,4000);
    FMV_t_clusters_beam = new TH1D("FMV_t_clusters_beam",title_fmv_beam.str().c_str(),250,0,4000);
    FMV_t_clusters_cosmic = new TH1D("FMV_t_clusters_cosmic",title_fmv_cosmic.str().c_str(),250,0,4000);
    PMT_t_clusters = new TH1D("PMT_t_clusters",title_pmt.str().c_str(),250,0,2000);
    PMT_t_clusters_2pe = new TH1D("PMT_t_clusters_2pe",title_pmt_2pe.str().c_str(),250,0,2000);
    PMT_t_clusters_5pe = new TH1D("PMT_t_clusters_5pe",title_pmt_5pe.str().c_str(),250,0,2000);
    PMT_t_clusters_10pe = new TH1D("PMT_t_clusters_10pe",title_pmt_10pe.str().c_str(),250,0,2000);
    PMT_t_clusters_30pe = new TH1D("PMT_t_clusters_30pe",title_pmt_30pe.str().c_str(),250,0,2000);
    PMT_t_clusters_cosmics = new TH1D("PMT_t_clusters_cosmics",title_pmt.str().c_str(),250,0,2000);
    PMT_t_clusters_2pe_cosmics = new TH1D("PMT_t_clusters_2pe_cosmics",title_pmt_2pe.str().c_str(),250,0,2000);
    PMT_t_clusters_5pe_cosmics = new TH1D("PMT_t_clusters_5pe_cosmics",title_pmt_5pe.str().c_str(),250,0,2000);
    PMT_t_clusters_10pe_cosmics = new TH1D("PMT_t_clusters_10pe_cosmics",title_pmt_10pe.str().c_str(),250,0,2000);
    PMT_t_clusters_30pe_cosmics = new TH1D("PMT_t_clusters_30pe_cosmics",title_pmt_30pe.str().c_str(),250,0,2000);
    PMT_t_clusters_led = new TH1D("PMT_t_clusters_led",title_pmt_led.str().c_str(),250,0,2000);
    PMT_t_clusters_2pe_led = new TH1D("PMT_t_clusters_2pe_led",title_pmt_2pe_led.str().c_str(),250,0,2000);
    PMT_t_clusters_5pe_led = new TH1D("PMT_t_clusters_5pe_led",title_pmt_5pe_led.str().c_str(),250,0,2000);
    PMT_t_clusters_10pe_led = new TH1D("PMT_t_clusters_10pe_led",title_pmt_10pe_led.str().c_str(),250,0,2000);
    PMT_t_clusters_30pe_led = new TH1D("PMT_t_clusters_30pe_led",title_pmt_30pe_led.str().c_str(),250,0,2000);
    PMT_t_clusters_full = new TH1D("PMT_t_clusters_full",title_pmt.str().c_str(),500,0,75000);
    PMT_t_clusters_2pe_full = new TH1D("PMT_t_clusters_2pe_full",title_pmt_2pe.str().c_str(),500,0,75000);
    PMT_t_clusters_5pe_full = new TH1D("PMT_t_clusters_5pe_full",title_pmt_5pe.str().c_str(),500,0,75000);
    PMT_t_clusters_10pe_full = new TH1D("PMT_t_clusters_10pe_full",title_pmt_10pe.str().c_str(),500,0,75000);
    PMT_t_clusters_30pe_full = new TH1D("PMT_t_clusters_30pe_full",title_pmt_30pe.str().c_str(),500,0,75000);
    MRD_PMT_t = new TH2D("MRD_PMT_t",title_mrd_pmt.str().c_str(),50,0,4000,50,0,2000);
    MRD_PMT_t_100pe = new TH2D("MRD_PMT_t_100pe",title_mrd_pmt_100pe.str().c_str(),50,0,4000,50,0,2000);
    MRD_PMT_t_Beam = new TH2D("MRD_PMT_t_Beam",title_mrd_pmt_beam.str().c_str(),50,0,4000,50,0,2000);
    MRD_PMT_t_100pe_Beam = new TH2D("MRD_PMT_t_100pe_Beam",title_mrd_pmt_100pe_beam.str().c_str(),50,0,4000,50,0,2000);
    MRD_PMT_t_Cosmic = new TH2D("MRD_PMT_t_Cosmic",title_mrd_pmt_cosmic.str().c_str(),50,0,4000,50,0,2000);
    MRD_PMT_t_100pe_Cosmic = new TH2D("MRD_PMT_t_100pe_Cosmic",title_mrd_pmt_100pe_cosmic.str().c_str(),50,0,4000,50,0,2000);
    FMV_PMT_t = new TH2D("FMV_PMT_t",title_fmv_pmt.str().c_str(),50,0,4000,50,0,2000);
    FMV_PMT_t_100pe = new TH2D("FMV_PMT_t_100pe",title_fmv_pmt_100pe.str().c_str(),50,0,4000,50,0,2000);
    FMV_PMT_t_Beam = new TH2D("FMV_PMT_t_Beam",title_fmv_pmt_beam.str().c_str(),50,0,4000,50,0,2000);
    FMV_PMT_t_100pe_Beam = new TH2D("FMV_PMT_t_100pe_Beam",title_fmv_pmt_100pe_beam.str().c_str(),50,0,4000,50,0,2000);
    FMV_PMT_t_Cosmic = new TH2D("FMV_PMT_t_Cosmic",title_fmv_pmt_cosmic.str().c_str(),50,0,4000,50,0,2000);
    FMV_PMT_t_100pe_Cosmic = new TH2D("FMV_PMT_t_100pe_Cosmic",title_fmv_pmt_100pe_cosmic.str().c_str(),50,0,4000,50,0,2000);
    MRD_PMT_Deltat = new TH1D("MRD_PMT_Deltat",title_mrd_pmt_delta.str().c_str(),500,-2000,4000);
    MRD_PMT_Deltat_100pe = new TH1D("MRD_PMT_Deltat_100pe",title_mrd_pmt_delta_100pe.str().c_str(),500,-2000,4000);
    MRD_PMT_Deltat_Beam = new TH1D("MRD_PMT_Deltat_Beam",title_mrd_pmt_delta_beam.str().c_str(),500,-2000,4000);
    MRD_PMT_Deltat_100pe_Beam = new TH1D("MRD_PMT_Deltat_100pe_Beam",title_mrd_pmt_delta_100pe_beam.str().c_str(),500,-2000,4000);
    MRD_PMT_Deltat_Cosmic = new TH1D("MRD_PMT_Deltat_Cosmic",title_mrd_pmt_delta_cosmic.str().c_str(),500,-2000,4000);
    MRD_PMT_Deltat_100pe_Cosmic = new TH1D("MRD_PMT_Deltat_100pe_Cosmic",title_mrd_pmt_delta_100pe_cosmic.str().c_str(),500,-2000,4000);
    FMV_PMT_Deltat = new TH1D("FMV_PMT_Deltat",title_fmv_pmt_delta.str().c_str(),500,-2000,4000);
    FMV_PMT_Deltat_100pe = new TH1D("FMV_PMT_Deltat_100pe",title_fmv_pmt_delta_100pe.str().c_str(),500,-2000,4000);
    FMV_PMT_Deltat_Beam = new TH1D("FMV_PMT_Deltat_Beam",title_fmv_pmt_delta_beam.str().c_str(),500,-2000,4000);
    FMV_PMT_Deltat_Beam_100pe = new TH1D("FMV_PMT_Deltat_Beam_100pe",title_fmv_pmt_delta_beam_100pe.str().c_str(),500,-2000,4000);
    FMV_PMT_Deltat_Cosmic = new TH1D("FMV_PMT_Deltat_Cosmic",title_fmv_pmt_delta_cosmic.str().c_str(),500,-2000,4000);
    FMV_PMT_Deltat_Cosmic_100pe = new TH1D("FMV_PMT_Deltat_Cosmic_100pe",title_fmv_pmt_delta_cosmic_100pe.str().c_str(),500,-2000,4000);
    MRD_FMV_Deltat = new TH1D("MRD_FMV_Deltat",title_mrd_fmv_delta.str().c_str(),500,-2000,4000);
    MRD_FMV_Deltat_Beam = new TH1D("MRD_FMV_Deltat_Beam",title_mrd_fmv_delta_beam.str().c_str(),500,-2000,4000);
    MRD_FMV_Deltat_Cosmic = new TH1D("MRD_FMV_Deltat_Cosmic",title_mrd_fmv_delta_cosmic.str().c_str(),500,-2000,4000);
    PMT_prompt_charge = new TH1D("PMT_prompt_charge",title_pmt_prompt.str().c_str(),200,0,5000);
    PMT_delayed_charge = new TH1D("PMT_delayed_charge",title_pmt_delayed.str().c_str(),200,0,5000);
    PMT_prompt_charge_zoom = new TH1D("PMT_prompt_charge_zoom",title_pmt_prompt.str().c_str(),200,0,500);
    PMT_delayed_charge_zoom = new TH1D("PMT_delayed_charge_zoom",title_pmt_delayed.str().c_str(),200,0,200);
    PMT_prompt_charge_10hits = new TH1D("PMT_prompt_charge_10hits",title_pmt_prompt_10.str().c_str(),200,0,500);
    PMT_delayed_charge_10hits = new TH1D("PMT_delayed_charge_10hits",title_pmt_delayed_10.str().c_str(),200,0,200);
    PMT_chargeperpmt = new TH1D("PMT_chargeperpmt",title_chargeperpmt.str().c_str(),50,0,30);
    PMT_chargeperpmt_100pe = new TH1D("PMT_chargeperpmt_100pe",title_chargeperpmt_100pe.str().c_str(),50,0,30);
    PMT_prompt_charge_MRDCoinc = new TH1D("PMT_prompt_charge_MRDCoinc",title_pmt_prompt_MRDCoinc.str().c_str(),200,0,5000);
    PMT_prompt_charge_MRDCoinc_NoFMV = new TH1D("PMT_prompt_charge_MRDCoinc_NoFMV",title_pmt_prompt_MRDCoinc_NoFMV.str().c_str(),200,0,5000);
    PMT_prompt_charge_FMV = new TH1D("PMT_prompt_charge_FMV",title_pmt_prompt_FMV.str().c_str(),200,0,5000);
    PMT_prompt_charge_Beam = new TH1D("PMT_prompt_charge_Beam",title_pmt_prompt_beam.str().c_str(),200,0,5000);
    PMT_prompt_charge_Cosmic = new TH1D("PMT_prompt_charge_Cosmic",title_pmt_prompt_cosmic.str().c_str(),200,0,5000);
    PMT_prompt_charge_LED = new TH1D("PMT_prompt_charge_LED",title_pmt_prompt_led.str().c_str(),200,0,5000);
    PMT_prompt_charge_CB = new TH2D("PMT_prompt_charge_CB",title_pmt_prompt_CB.str().c_str(),200,0,200,200,0,1);
    PMT_delayed_charge_CB = new TH2D("PMT_delayed_charge_CB",title_pmt_delayed_CB.str().c_str(),200,0,200,200,0,1);
    PMT_Channelkeys = new TH1D("PMT_Channelkeys",title_pmt_chankeys.str().c_str(),132,332,446);
    FMV_Channelkeys = new TH1D("FMV_Channelkeys",title_fmv_chankeys.str().c_str(),26,0,26);
    MRD_Channelkeys = new TH1D("MRD_Channelkeys",title_mrd_chankeys.str().c_str(),306,26,332);
    ANNIE_counts = new TH1D("ANNIE_counts",title_counts.str().c_str(),10,0,10);  
    ANNIE_rates = new TH1D("ANNIE_rates",title_rates.str().c_str(),10,0,10);  
    ANNIE_fractions = new TH1D("ANNIE_fractions",title_fractions.str().c_str(),10,0,10);  
    PMT_DelayedMult = new TH1D("PMT_DelayedMult",title_delayed_mult.str().c_str(),20,0,20);
    PMT_DelayedMult_Coinc = new TH1D("PMT_DelayedMult_Coinc",title_delayed_mult_coinc.str().c_str(),20,0,20);
    PMT_DelayedMult_Coinc_NoFMV = new TH1D("PMT_DelayedMult_Coinc_NoFMV",title_delayed_mult_coinc_nofmv.str().c_str(),20,0,20);
    PMT_DelayedMult_Coinc_NoFMV_CB = new TH1D("PMT_DelayedMult_Coinc_NoFMV_CB",title_delayed_mult_coinc_nofmv_cb.str().c_str(),20,0,20);
    PMT_DelayedTime = new TH1D("PMT_DelayedTime",title_delayed_time.str().c_str(),100,12000,67000);
    PMT_DelayedTime_CB = new TH1D("PMT_DelayedTime_CB",title_delayed_time_cb.str().c_str(),100,12000,67000);
    PMT_DelayedTime_Coinc = new TH1D("PMT_DelayedTime_Coinc",title_delayed_time_coinc.str().c_str(),100,12000,67000);
    PMT_DelayedTime_Coinc_NoFMV = new TH1D("PMT_DelayedTime_Coinc_NoFMV",title_delayed_time_coinc_nofmv.str().c_str(),100,12000,67000);
    PMT_DelayedTime_Coinc_NoFMV_CB = new TH1D("PMT_DelayedTime_Coinc_NoFMV_CB",title_delayed_time_coinc_nofmv_cb.str().c_str(),100,12000,67000);
    ADCWaveform_Samples = new TH1D("ADCWaveform_Samples",title_adcwaveform_samples.str().c_str(),5000,0,50000);  
    Triggerwords = new TH1D("Triggerwords",title_triggerwords.str().c_str(),60,0,60);

    MRD_t_clusters->GetXaxis()->SetTitle("t_{MRD} [ns]");
    MRD_t_clusters_beam->GetXaxis()->SetTitle("t_{MRD} [ns]");
    MRD_t_clusters_cosmic->GetXaxis()->SetTitle("t_{MRD} [ns]");
    FMV_t_clusters->GetXaxis()->SetTitle("t_{FMV} [ns]");
    FMV_t_clusters_beam->GetXaxis()->SetTitle("t_{FMV} [ns]");
    FMV_t_clusters_cosmic->GetXaxis()->SetTitle("t_{FMV} [ns]");
    PMT_t_clusters->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_2pe->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_5pe->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_10pe->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_30pe->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_cosmics->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_2pe_cosmics->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_5pe_cosmics->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_10pe_cosmics->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_30pe_cosmics->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_led->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_2pe_led->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_5pe_led->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_10pe_led->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_30pe_led->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_full->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_2pe_full->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_5pe_full->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_10pe_full->GetXaxis()->SetTitle("t_{PMT} [ns]");
    PMT_t_clusters_30pe_full->GetXaxis()->SetTitle("t_{PMT} [ns]");
    MRD_PMT_t->GetXaxis()->SetTitle("t_{MRD} [ns]");
    MRD_PMT_t->GetYaxis()->SetTitle("t_{PMT} [ns]");
    MRD_PMT_t_100pe->GetXaxis()->SetTitle("t_{MRD} [ns]");
    MRD_PMT_t_100pe->GetYaxis()->SetTitle("t_{PMT} [ns]");
    MRD_PMT_t_Beam->GetXaxis()->SetTitle("t_{MRD} [ns]");
    MRD_PMT_t_Beam->GetYaxis()->SetTitle("t_{PMT} [ns]");
    MRD_PMT_t_100pe_Beam->GetXaxis()->SetTitle("t_{MRD} [ns]");
    MRD_PMT_t_100pe_Beam->GetYaxis()->SetTitle("t_{PMT} [ns]");
    MRD_PMT_t_Cosmic->GetXaxis()->SetTitle("t_{MRD} [ns]");
    MRD_PMT_t_Cosmic->GetYaxis()->SetTitle("t_{PMT} [ns]");
    MRD_PMT_t_100pe_Cosmic->GetXaxis()->SetTitle("t_{MRD} [ns]");
    MRD_PMT_t_100pe_Cosmic->GetYaxis()->SetTitle("t_{PMT} [ns]");
    FMV_PMT_t->GetXaxis()->SetTitle("t_{FMV} [ns]");
    FMV_PMT_t->GetYaxis()->SetTitle("t_{PMT} [ns]");
    FMV_PMT_t_100pe->GetXaxis()->SetTitle("t_{FMV} [ns]");
    FMV_PMT_t_100pe->GetYaxis()->SetTitle("t_{PMT} [ns]");
    FMV_PMT_t_Beam->GetXaxis()->SetTitle("t_{FMV} [ns]");
    FMV_PMT_t_Beam->GetYaxis()->SetTitle("t_{PMT} [ns]");
    FMV_PMT_t_100pe_Beam->GetXaxis()->SetTitle("t_{FMV} [ns]");
    FMV_PMT_t_100pe_Beam->GetYaxis()->SetTitle("t_{PMT} [ns]");
    FMV_PMT_t_Cosmic->GetXaxis()->SetTitle("t_{FMV} [ns]");
    FMV_PMT_t_Cosmic->GetYaxis()->SetTitle("t_{PMT} [ns]");
    FMV_PMT_t_100pe_Cosmic->GetXaxis()->SetTitle("t_{FMV} [ns]");
    FMV_PMT_t_100pe_Cosmic->GetYaxis()->SetTitle("t_{PMT} [ns]");
    MRD_PMT_Deltat->GetXaxis()->SetTitle("t_{MRD}-t_{PMT} [ns]");
    MRD_PMT_Deltat_100pe->GetXaxis()->SetTitle("t_{MRD}-t_{PMT} [ns]");
    MRD_PMT_Deltat_Beam->GetXaxis()->SetTitle("t_{MRD}-t_{PMT} [ns]");
    MRD_PMT_Deltat_100pe_Beam->GetXaxis()->SetTitle("t_{MRD}-t_{PMT} [ns]");
    MRD_PMT_Deltat_Cosmic->GetXaxis()->SetTitle("t_{MRD}-t_{PMT} [ns]");
    MRD_PMT_Deltat_100pe_Cosmic->GetXaxis()->SetTitle("t_{MRD}-t_{PMT} [ns]");
    FMV_PMT_Deltat->GetXaxis()->SetTitle("t_{FMV}-t_{PMT} [ns]");
    FMV_PMT_Deltat_100pe->GetXaxis()->SetTitle("t_{FMV}-t_{PMT} [ns]");
    FMV_PMT_Deltat_Beam->GetXaxis()->SetTitle("t_{FMV}-t_{PMT} [ns]");
    FMV_PMT_Deltat_Beam_100pe->GetXaxis()->SetTitle("t_{FMV}-t_{PMT} [ns]");
    FMV_PMT_Deltat_Cosmic->GetXaxis()->SetTitle("t_{FMV}-t_{PMT} [ns]");
    FMV_PMT_Deltat_Cosmic_100pe->GetXaxis()->SetTitle("t_{FMV}-t_{PMT} [ns]");
    MRD_FMV_Deltat->GetXaxis()->SetTitle("t_{MRD}-t_{FMV} [ns]");
    MRD_FMV_Deltat_Beam->GetXaxis()->SetTitle("t_{MRD}-t_{FMV} [ns]");
    MRD_FMV_Deltat_Cosmic->GetXaxis()->SetTitle("t_{MRD}-t_{FMV} [ns]");
    PMT_prompt_charge->GetXaxis()->SetTitle("q_{prompt} [p.e.]");
    PMT_prompt_charge_10hits->GetXaxis()->SetTitle("q_{prompt} [p.e.]");
    PMT_prompt_charge_zoom->GetXaxis()->SetTitle("q_{prompt} [p.e.]");
    PMT_prompt_charge_MRDCoinc->GetXaxis()->SetTitle("q_{prompt} [p.e.]");
    PMT_prompt_charge_MRDCoinc_NoFMV->GetXaxis()->SetTitle("q_{prompt} [p.e.]");
    PMT_prompt_charge_FMV->GetXaxis()->SetTitle("q_{prompt} [p.e.]");
    PMT_prompt_charge_Beam->GetXaxis()->SetTitle("q_{prompt} [p.e.]");
    PMT_prompt_charge_Cosmic->GetXaxis()->SetTitle("q_{prompt} [p.e.]");
    PMT_prompt_charge_LED->GetXaxis()->SetTitle("q_{prompt} [p.e.]");
    PMT_delayed_charge->GetXaxis()->SetTitle("q_{delayed} [p.e.]");
    PMT_delayed_charge_10hits->GetXaxis()->SetTitle("q_{delayed} [p.e.]");
    PMT_delayed_charge_zoom->GetXaxis()->SetTitle("q_{delayed} [p.e.]");
    PMT_chargeperpmt->GetXaxis()->SetTitle("q/n_{pmt} [p.e.]");
    PMT_chargeperpmt_100pe->GetXaxis()->SetTitle("q/n_{pmt} [p.e.]");
    PMT_prompt_charge_CB->GetXaxis()->SetTitle("q_{prompt} [p.e.]");
    PMT_prompt_charge_CB->GetYaxis()->SetTitle("Charge Balance");
    PMT_prompt_charge_CB->SetStats(0);
    PMT_delayed_charge_CB->GetXaxis()->SetTitle("q_{delayed} [p.e.]");
    PMT_delayed_charge_CB->GetYaxis()->SetTitle("Charge Balance");
    PMT_delayed_charge_CB->SetStats(0);
    PMT_Channelkeys->GetXaxis()->SetTitle("Channelkey");
    MRD_Channelkeys->GetXaxis()->SetTitle("Channelkey");
    FMV_Channelkeys->GetXaxis()->SetTitle("Channelkey");
    ANNIE_rates->GetYaxis()->SetTitle("Rate [Hz]");
    ANNIE_counts->GetYaxis()->SetTitle("#");
    ANNIE_fractions->GetYaxis()->SetTitle("Fraction of events [%]");
    ANNIE_rates->SetStats(0);
    ANNIE_counts->SetStats(0);
    ANNIE_fractions->SetStats(0);
    PMT_DelayedMult->GetXaxis()->SetTitle("Neutron candidate multiplicity");
    PMT_DelayedMult_Coinc->GetXaxis()->SetTitle("Neutron candidate multiplicity");
    PMT_DelayedMult_Coinc_NoFMV->GetXaxis()->SetTitle("Neutron candidate multiplicity");
    PMT_DelayedMult_Coinc_NoFMV_CB->GetXaxis()->SetTitle("Neutron candidate multiplicity");
    PMT_DelayedMult->GetYaxis()->SetTitle("#");
    PMT_DelayedMult_Coinc->GetYaxis()->SetTitle("#");
    PMT_DelayedMult_Coinc_NoFMV->GetYaxis()->SetTitle("#");
    PMT_DelayedMult_Coinc_NoFMV_CB->GetYaxis()->SetTitle("#");
    PMT_DelayedTime->GetXaxis()->SetTitle("Cluster time [ns]");
    PMT_DelayedTime_CB->GetXaxis()->SetTitle("Cluster time [ns]");
    PMT_DelayedTime_Coinc->GetXaxis()->SetTitle("Cluster time [ns]");
    PMT_DelayedTime_Coinc_NoFMV->GetXaxis()->SetTitle("Cluster time [ns]");
    PMT_DelayedTime_Coinc_NoFMV_CB->GetXaxis()->SetTitle("Cluster time [ns]");
    PMT_DelayedTime->GetYaxis()->SetTitle("#");
    PMT_DelayedTime_Coinc->GetYaxis()->SetTitle("#");
    PMT_DelayedTime_Coinc_NoFMV->GetYaxis()->SetTitle("#");
    PMT_DelayedTime_Coinc_NoFMV_CB->GetYaxis()->SetTitle("#");
    ADCWaveform_Samples->GetXaxis()->SetTitle("Number of samples in waveform (ADC)");
    ADCWaveform_Samples->GetYaxis()->SetTitle("#");
    Triggerwords->SetStats(0);
    Triggerwords->GetXaxis()->SetTitle("Triggerword");
    Triggerwords->GetYaxis()->SetTitle("#");

    canvas_beamspill = new TCanvas("canvas_beamspill","Beamspill",900,600);
    canvas_mrd_pmt = new TCanvas("canvas_mrd_pmt","MRD-PMT alignment",900,600);
    canvas_triggerwords = new TCanvas("canvas_triggerwords","Triggerwords",900,600);
    canvas_prompt_charge = new TCanvas("canvas_prompt_charge","Prompt charge",900,600);
    canvas_neutron_mult = new TCanvas("canvas_neutron_mult","Neutron multiplicity",900,600);
    canvas_channelkey = new TCanvas("canvas_channelkey","Channelkeys",900,600);
    canvas_rates = new TCanvas("canvas_rates","Rates",900,600);

    gROOT->cd();

    first_entry = false;

}
