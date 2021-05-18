#include "PhaseIITreeMaker.h"

PhaseIITreeMaker::PhaseIITreeMaker():Tool(){}


bool PhaseIITreeMaker::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////
  
  m_variables.Get("verbose", verbosity);
  m_variables.Get("TankHitInfo_fill", TankHitInfo_fill);
  m_variables.Get("MRDHitInfo_fill", MRDHitInfo_fill);
  m_variables.Get("fillCleanEventsOnly", fillCleanEventsOnly);
  m_variables.Get("MCTruth_fill", MCTruth_fill);
  m_variables.Get("MRDReco_fill", MRDReco_fill);
  m_variables.Get("TankReco_fill", TankReco_fill);
  m_variables.Get("RecoDebug_fill", RecoDebug_fill);
  m_variables.Get("muonTruthRecoDiff_fill", muonTruthRecoDiff_fill);

  m_variables.Get("SiPMPulseInfo_fill",SiPMPulseInfo_fill);
  m_variables.Get("TankClusterProcessing",TankClusterProcessing);
  m_variables.Get("MRDClusterProcessing",MRDClusterProcessing);
  m_variables.Get("TriggerProcessing",TriggerProcessing);

  std::string output_filename;
  m_variables.Get("OutputFile", output_filename);
  fOutput_tfile = new TFile(output_filename.c_str(), "recreate");
  fPhaseIITankClusterTree = new TTree("phaseIITankClusterTree", "ANNIE Phase II Tank Cluster Tree");
  fPhaseIIMRDClusterTree = new TTree("phaseIIMRDClusterTree", "ANNIE Phase II MRD Cluster Tree");
  fPhaseIITrigTree = new TTree("phaseIITriggerTree", "ANNIE Phase II Ntuple Trigger Tree");

  m_data->CStore.Get("AuxChannelNumToTypeMap",AuxChannelNumToTypeMap);
  m_data->CStore.Get("ChannelNumToTankPMTSPEChargeMap",ChannelKeyToSPEMap);

  auto get_geometry= m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",geom);
  if(!get_geometry){
  	Log("PhaseIITreeMaker Tool: Error retrieving Geometry from ANNIEEvent!",v_error,verbosity); 
  	return false; 
  }
 
  if(TankClusterProcessing){
    fPhaseIITankClusterTree->Branch("runNumber",&fRunNumber,"runNumber/I");
    fPhaseIITankClusterTree->Branch("subrunNumber",&fSubrunNumber,"subrunNumber/I");
    fPhaseIITankClusterTree->Branch("runType",&fRunType,"runType/I");
    fPhaseIITankClusterTree->Branch("startTime",&fStartTime_Tree,"startTime/l");

    //Some lower level information to save
    fPhaseIITankClusterTree->Branch("eventNumber",&fEventNumber,"eventNumber/I");
    fPhaseIITankClusterTree->Branch("eventTimeTank",&fEventTimeTank_Tree,"eventTimeTank/l");
    fPhaseIITankClusterTree->Branch("clusterNumber",&fClusterNumber,"clusterNumber/I");
    fPhaseIITankClusterTree->Branch("clusterTime",&fClusterTime,"clusterTime/D");
    fPhaseIITankClusterTree->Branch("clusterCharge",&fClusterCharge,"clusterCharge/D");
    fPhaseIITankClusterTree->Branch("clusterPE",&fClusterPE,"clusterPE/D");
    fPhaseIITankClusterTree->Branch("clusterMaxPE",&fClusterMaxPE,"clusterMaxPE/D");
    fPhaseIITankClusterTree->Branch("clusterChargePointX",&fClusterChargePointX,"clusterChargePointX/D");
    fPhaseIITankClusterTree->Branch("clusterChargePointY",&fClusterChargePointY,"clusterChargePointY/D");
    fPhaseIITankClusterTree->Branch("clusterChargePointZ",&fClusterChargePointZ,"clusterChargePointZ/D");
    fPhaseIITankClusterTree->Branch("clusterChargeBalance",&fClusterChargeBalance,"clusterChargeBalance/D");
    fPhaseIITankClusterTree->Branch("clusterHits",&fClusterHits,"clusterHits/i");
    fPhaseIITankClusterTree->Branch("trigword",&fTriggerword);
    fPhaseIITankClusterTree->Branch("TankMRDCoinc",&fTankMRDCoinc);
    fPhaseIITankClusterTree->Branch("NoVeto",&fNoVeto);
    fPhaseIITankClusterTree->Branch("ADCSamples",&fADCWaveformSamples);
    fPhaseIITankClusterTree->Branch("ADCChankeys",&fADCWaveformChankeys);
    fPhaseIITankClusterTree->Branch("Extended",&fExtended);
    fPhaseIITankClusterTree->Branch("beam_pot",&fPot);
    fPhaseIITankClusterTree->Branch("beam_ok",&fBeamok);
    if(TankHitInfo_fill){
      fPhaseIITankClusterTree->Branch("filter",&fIsFiltered);
      fPhaseIITankClusterTree->Branch("hitX",&fHitX);
      fPhaseIITankClusterTree->Branch("hitY",&fHitY);
      fPhaseIITankClusterTree->Branch("hitZ",&fHitZ);
      fPhaseIITankClusterTree->Branch("hitT",&fHitT);
      fPhaseIITankClusterTree->Branch("hitQ",&fHitQ);
      fPhaseIITankClusterTree->Branch("hitPE",&fHitPE);
      fPhaseIITankClusterTree->Branch("hitType", &fHitType);
      fPhaseIITankClusterTree->Branch("hitDetID", &fHitDetID);
      fPhaseIITankClusterTree->Branch("hitChankey",&fHitChankey);
    }
    //SiPM Pulse Info; load into both trees for now...
    if(SiPMPulseInfo_fill){
      fPhaseIITankClusterTree->Branch("SiPMhitQ",&fSiPMHitQ);
      fPhaseIITankClusterTree->Branch("SiPMhitT",&fSiPMHitT);
      fPhaseIITankClusterTree->Branch("SiPMhitAmplitude",&fSiPMHitAmplitude);
      fPhaseIITankClusterTree->Branch("SiPMNum",&fSiPMNum);
      fPhaseIITankClusterTree->Branch("SiPM1NPulses",&fSiPM1NPulses,"SiPM1NPulses/I");
      fPhaseIITankClusterTree->Branch("SiPM2NPulses",&fSiPM2NPulses,"SiPM2NPulses/I");
    }
  } 

  if(MRDClusterProcessing){
    fPhaseIIMRDClusterTree->Branch("runNumber",&fRunNumber,"runNumber/I");
    fPhaseIIMRDClusterTree->Branch("subrunNumber",&fSubrunNumber,"subrunNumber/I");
    fPhaseIIMRDClusterTree->Branch("runType",&fRunType,"runType/I");
    fPhaseIIMRDClusterTree->Branch("startTime",&fStartTime_Tree,"startTime/l");


    //Some lower level information to save
    fPhaseIIMRDClusterTree->Branch("eventNumber",&fEventNumber,"eventNumber/I");
    fPhaseIIMRDClusterTree->Branch("eventTimeMRD",&fEventTimeMRD_Tree,"eventTimeMRD/l");
    fPhaseIIMRDClusterTree->Branch("eventTimeTank",&fEventTimeTank_Tree,"eventTimeTank/l");
    fPhaseIIMRDClusterTree->Branch("clusterNumber",&fMRDClusterNumber,"clusterNumber/I");
    fPhaseIIMRDClusterTree->Branch("clusterTime",&fMRDClusterTime,"clusterTime/D");
    fPhaseIIMRDClusterTree->Branch("clusterTimeSigma",&fMRDClusterTimeSigma,"clusterTimeSigma/D");
    fPhaseIIMRDClusterTree->Branch("clusterHits",&fMRDClusterHits,"clusterHits/i");
    fPhaseIIMRDClusterTree->Branch("trigword",&fTriggerword);
    fPhaseIIMRDClusterTree->Branch("TankMRDCoinc",&fTankMRDCoinc);
    fPhaseIIMRDClusterTree->Branch("NoVeto",&fNoVeto);
    fPhaseIIMRDClusterTree->Branch("Extended",&fExtended);
    fPhaseIIMRDClusterTree->Branch("beam_pot",&fPot);
    fPhaseIIMRDClusterTree->Branch("beam_ok",&fBeamok);

    if(MRDHitInfo_fill){
      fPhaseIIMRDClusterTree->Branch("MRDhitT",&fMRDHitT);
      fPhaseIIMRDClusterTree->Branch("MRDhitDetID", &fMRDHitDetID);
      fPhaseIIMRDClusterTree->Branch("MRDhitChankey", &fMRDHitChankey);
      fPhaseIIMRDClusterTree->Branch("FMVhitT",&fFMVHitT);
      fPhaseIIMRDClusterTree->Branch("FMVhitDetID",&fFMVHitDetID);
      fPhaseIIMRDClusterTree->Branch("FMVhitChankey",&fFMVHitChankey);
      fPhaseIIMRDClusterTree->Branch("vetoHit",&fVetoHit,"vetoHit/I");
    }
    if(MRDReco_fill){
      fPhaseIIMRDClusterTree->Branch("numClusterTracks",&fNumClusterTracks,"numClusterTracks/I");
    //Push back some properties
      fPhaseIIMRDClusterTree->Branch("MRDTrackAngle",&fMRDTrackAngle);
      fPhaseIIMRDClusterTree->Branch("MRDTrackAngleError",&fMRDTrackAngleError);
      fPhaseIIMRDClusterTree->Branch("MRDPenetrationDepth",&fMRDPenetrationDepth);
      fPhaseIIMRDClusterTree->Branch("MRDTrackLength",&fMRDTrackLength);
      fPhaseIIMRDClusterTree->Branch("MRDEntryPointRadius",&fMRDEntryPointRadius);
      fPhaseIIMRDClusterTree->Branch("MRDEnergyLoss",&fMRDEnergyLoss);
      fPhaseIIMRDClusterTree->Branch("MRDEnergyLossError",&fMRDEnergyLossError);
      fPhaseIIMRDClusterTree->Branch("MRDTrackStartX",&fMRDTrackStartX);
      fPhaseIIMRDClusterTree->Branch("MRDTrackStartY",&fMRDTrackStartY);
      fPhaseIIMRDClusterTree->Branch("MRDTrackStartZ",&fMRDTrackStartZ);
      fPhaseIIMRDClusterTree->Branch("MRDTrackStopX",&fMRDTrackStopX);
      fPhaseIIMRDClusterTree->Branch("MRDTrackStopY",&fMRDTrackStopY);
      fPhaseIIMRDClusterTree->Branch("MRDTrackStopZ",&fMRDTrackStopZ);
      fPhaseIIMRDClusterTree->Branch("MRDSide",&fMRDSide);
      fPhaseIIMRDClusterTree->Branch("MRDStop",&fMRDStop);
      fPhaseIIMRDClusterTree->Branch("MRDThrough",&fMRDThrough);
    }
  }

  if(TriggerProcessing){
    //Metadata for standard Events
    fPhaseIITrigTree->Branch("runNumber",&fRunNumber,"runNumber/I");
    fPhaseIITrigTree->Branch("subrunNumber",&fSubrunNumber,"subrunNumber/I");
    fPhaseIITrigTree->Branch("runType",&fRunType,"runType/I");
    fPhaseIITrigTree->Branch("startTime",&fStartTime_Tree,"startTime/l");

    //Some lower level information to save
    fPhaseIITrigTree->Branch("eventNumber",&fEventNumber,"eventNumber/I");
    fPhaseIITrigTree->Branch("eventTimeTank",&fEventTimeTank_Tree,"eventTimeTank/l");
    fPhaseIITrigTree->Branch("eventTimeMRD",&fEventTimeMRD_Tree,"eventTimeMRD/l");
    fPhaseIITrigTree->Branch("nhits",&fNHits,"nhits/I");

    //CTC & event selection related variables
    fPhaseIITrigTree->Branch("trigword",&fTriggerword);
    fPhaseIITrigTree->Branch("HasTank",&fHasTank);
    fPhaseIITrigTree->Branch("HasMRD",&fHasMRD);
    fPhaseIITrigTree->Branch("TankMRDCoinc",&fTankMRDCoinc);
    fPhaseIITrigTree->Branch("NoVeto",&fNoVeto);

    //Extended window information
    fPhaseIITrigTree->Branch("Extended",&fExtended);
    
    //Beam information
    fPhaseIITrigTree->Branch("beam_pot",&fPot);
    fPhaseIITrigTree->Branch("beam_ok",&fBeamok);

    //Event Staus Flag Information
    if(fillCleanEventsOnly){
      fPhaseIITrigTree->Branch("eventStatusApplied",&fEventStatusApplied,"eventStatusApplied/I");
      fPhaseIITrigTree->Branch("eventStatusFlagged",&fEventStatusFlagged,"eventStatusFlagged/I");
    }
    //Hit information (PMT and LAPPD)
    if(SiPMPulseInfo_fill){
      fPhaseIITrigTree->Branch("SiPM1NPulses",&fSiPM1NPulses,"SiPM1NPulses/I");
      fPhaseIITrigTree->Branch("SiPM2NPulses",&fSiPM2NPulses,"SiPM2NPulses/I");
      fPhaseIITrigTree->Branch("SiPMhitQ",&fSiPMHitQ);
      fPhaseIITrigTree->Branch("SiPMhitT",&fSiPMHitT);
      fPhaseIITrigTree->Branch("SiPMhitAmplitude",&fSiPMHitAmplitude);
      fPhaseIITrigTree->Branch("SiPMNum",&fSiPMNum);
    }

    if(TankHitInfo_fill){
      fPhaseIITrigTree->Branch("filter",&fIsFiltered);
      fPhaseIITrigTree->Branch("hitX",&fHitX);
      fPhaseIITrigTree->Branch("hitY",&fHitY);
      fPhaseIITrigTree->Branch("hitZ",&fHitZ);
      fPhaseIITrigTree->Branch("hitT",&fHitT);
      fPhaseIITrigTree->Branch("hitQ",&fHitQ);
      fPhaseIITrigTree->Branch("hitPE",&fHitPE);
      fPhaseIITrigTree->Branch("hitType", &fHitType);
      fPhaseIITrigTree->Branch("hitDetID", &fHitDetID);
    }

    if(MRDHitInfo_fill){
      fPhaseIITrigTree->Branch("MRDhitT",&fMRDHitT);
      fPhaseIITrigTree->Branch("MRDhitDetID", &fMRDHitDetID);
      fPhaseIITrigTree->Branch("MRDhitChankey", &fMRDHitChankey);
      fPhaseIITrigTree->Branch("FMVhitT",&fFMVHitT);
      fPhaseIITrigTree->Branch("FMVhitDetID",&fFMVHitDetID);
      fPhaseIITrigTree->Branch("FMVhitChankey",&fFMVHitChankey);
      fPhaseIITrigTree->Branch("vetoHit",&fVetoHit,"vetoHit/I");
    }
    if(MRDReco_fill){
      fPhaseIITrigTree->Branch("numMRDTracks",&fNumClusterTracks,"numMRDTracks/I");
    //Push back some properties
      fPhaseIITrigTree->Branch("MRDTrackAngle",&fMRDTrackAngle);
      fPhaseIITrigTree->Branch("MRDTrackAngleError",&fMRDTrackAngleError);
      fPhaseIITrigTree->Branch("MRDPenetrationDepth",&fMRDPenetrationDepth);
      fPhaseIITrigTree->Branch("MRDTrackLength",&fMRDTrackLength);
      fPhaseIITrigTree->Branch("MRDEntryPointRadius",&fMRDEntryPointRadius);
      fPhaseIITrigTree->Branch("MRDEnergyLoss",&fMRDEnergyLoss);
      fPhaseIITrigTree->Branch("MRDEnergyLossError",&fMRDEnergyLossError);
      fPhaseIITrigTree->Branch("MRDTrackStartX",&fMRDTrackStartX);
      fPhaseIITrigTree->Branch("MRDTrackStartY",&fMRDTrackStartY);
      fPhaseIITrigTree->Branch("MRDTrackStartZ",&fMRDTrackStartZ);
      fPhaseIITrigTree->Branch("MRDTrackStopX",&fMRDTrackStopX);
      fPhaseIITrigTree->Branch("MRDTrackStopY",&fMRDTrackStopY);
      fPhaseIITrigTree->Branch("MRDTrackStopZ",&fMRDTrackStopZ);
      fPhaseIITrigTree->Branch("MRDSide",&fMRDSide);
      fPhaseIITrigTree->Branch("MRDStop",&fMRDStop);
      fPhaseIITrigTree->Branch("MRDThrough",&fMRDThrough);
    }

    //Reconstructed variables after full Muon Reco Analysis
    if(TankReco_fill){
      fPhaseIITrigTree->Branch("recoVtxX",&fRecoVtxX,"recoVtxX/D");
      fPhaseIITrigTree->Branch("recoVtxY",&fRecoVtxY,"recoVtxY/D");
      fPhaseIITrigTree->Branch("recoVtxZ",&fRecoVtxZ,"recoVtxZ/D");
      fPhaseIITrigTree->Branch("recoVtxTime",&fRecoVtxTime,"recoVtxTime/D");
      fPhaseIITrigTree->Branch("recoDirX",&fRecoDirX,"recoDirX/D");
      fPhaseIITrigTree->Branch("recoDirY",&fRecoDirY,"recoDirY/D");
      fPhaseIITrigTree->Branch("recoDirZ",&fRecoDirZ,"recoDirZ/D");
      fPhaseIITrigTree->Branch("recoAngle",&fRecoAngle,"recoAngle/D");
      fPhaseIITrigTree->Branch("recoPhi",&fRecoPhi,"recoPhi/D");
      fPhaseIITrigTree->Branch("recoVtxFOM",&fRecoVtxFOM,"recoVtxFOM/D");
      fPhaseIITrigTree->Branch("recoStatus",&fRecoStatus,"recoStatus/I");
    }
  
    //MC truth information for muons
    //Output to tree when MCTruth_fill = 1 in config
    if (MCTruth_fill){
      fPhaseIITrigTree->Branch("triggerNumber",&fMCTriggerNum,"triggerNumber/I");
      fPhaseIITrigTree->Branch("mcEntryNumber",&fMCEventNum,"mcEntryNumber/I");
      fPhaseIITrigTree->Branch("trueVtxX",&fTrueVtxX,"trueVtxX/D");
      fPhaseIITrigTree->Branch("trueVtxY",&fTrueVtxY,"trueVtxY/D");
      fPhaseIITrigTree->Branch("trueVtxZ",&fTrueVtxZ,"trueVtxZ/D");
      fPhaseIITrigTree->Branch("trueVtxTime",&fTrueVtxTime,"trueVtxTime/D");
      fPhaseIITrigTree->Branch("trueDirX",&fTrueDirX,"trueDirX/D");
      fPhaseIITrigTree->Branch("trueDirY",&fTrueDirY,"trueDirY/D");
      fPhaseIITrigTree->Branch("trueDirZ",&fTrueDirZ,"trueDirZ/D");
      fPhaseIITrigTree->Branch("trueAngle",&fTrueAngle,"trueAngle/D");
      fPhaseIITrigTree->Branch("truePhi",&fTruePhi,"truePhi/D");
      fPhaseIITrigTree->Branch("trueMuonEnergy",&fTrueMuonEnergy, "trueMuonEnergy/D");
      fPhaseIITrigTree->Branch("trueTrackLengthInWater",&fTrueTrackLengthInWater,"trueTrackLengthInWater/D");
      fPhaseIITrigTree->Branch("trueTrackLengthInMRD",&fTrueTrackLengthInMRD,"trueTrackLengthInMRD/D");
      fPhaseIITrigTree->Branch("Pi0Count",&fPi0Count,"Pi0Count/I");
      fPhaseIITrigTree->Branch("PiPlusCount",&fPiPlusCount,"PiPlusCount/I");
      fPhaseIITrigTree->Branch("PiMinusCount",&fPiMinusCount,"PiMinusCount/I");
      fPhaseIITrigTree->Branch("K0Count",&fK0Count,"K0Count/I");
      fPhaseIITrigTree->Branch("KPlusCount",&fKPlusCount,"KPlusCount/I");
      fPhaseIITrigTree->Branch("KMinusCount",&fKMinusCount,"KMinusCount/I");
    }
  
    // Reconstructed variables from each step in Muon Reco Analysis
    // Currently output when RecoDebug_fill = 1 in config 
    if (RecoDebug_fill){
      fPhaseIITrigTree->Branch("seedVtxX",&fSeedVtxX); 
      fPhaseIITrigTree->Branch("seedVtxY",&fSeedVtxY); 
      fPhaseIITrigTree->Branch("seedVtxZ",&fSeedVtxZ);
      fPhaseIITrigTree->Branch("seedVtxFOM",&fSeedVtxFOM); 
      fPhaseIITrigTree->Branch("seedVtxTime",&fSeedVtxTime,"seedVtxTime/D");
    
      fPhaseIITrigTree->Branch("pointPosX",&fPointPosX,"pointPosX/D");
      fPhaseIITrigTree->Branch("pointPosY",&fPointPosY,"pointPosY/D");
      fPhaseIITrigTree->Branch("pointPosZ",&fPointPosZ,"pointPosZ/D");
      fPhaseIITrigTree->Branch("pointPosTime",&fPointPosTime,"pointPosTime/D");
      fPhaseIITrigTree->Branch("pointPosFOM",&fPointPosFOM,"pointPosFOM/D");
      fPhaseIITrigTree->Branch("pointPosStatus",&fPointPosStatus,"pointPosStatus/I");
      
      fPhaseIITrigTree->Branch("pointDirX",&fPointDirX,"pointDirX/D");
      fPhaseIITrigTree->Branch("pointDirY",&fPointDirY,"pointDirY/D");
      fPhaseIITrigTree->Branch("pointDirZ",&fPointDirZ,"pointDirZ/D");
      fPhaseIITrigTree->Branch("pointDirTime",&fPointDirTime,"pointDirTime/D");
      fPhaseIITrigTree->Branch("pointDirStatus",&fPointDirStatus,"pointDirStatus/I");
      fPhaseIITrigTree->Branch("pointDirFOM",&fPointDirFOM,"pointDirFOM/D");
    
      fPhaseIITrigTree->Branch("pointVtxPosX",&fPointVtxPosX,"pointVtxPosX/D");
      fPhaseIITrigTree->Branch("pointVtxPosY",&fPointVtxPosY,"pointVtxPosY/D");
      fPhaseIITrigTree->Branch("pointVtxPosZ",&fPointVtxPosZ,"pointVtxPosZ/D");
      fPhaseIITrigTree->Branch("pointVtxTime",&fPointVtxTime,"pointVtxTime/D");
      fPhaseIITrigTree->Branch("pointVtxDirX",&fPointVtxDirX,"pointVtxDirX/D");
      fPhaseIITrigTree->Branch("pointVtxDirY",&fPointVtxDirY,"pointVtxDirY/D");
      fPhaseIITrigTree->Branch("pointVtxDirZ",&fPointVtxDirZ,"pointVtxDirZ/D");
      fPhaseIITrigTree->Branch("pointVtxFOM",&fPointVtxFOM,"pointVtxFOM/D");
      fPhaseIITrigTree->Branch("pointVtxStatus",&fPointVtxStatus,"pointVtxStatus/I");
    } 

    // Difference in MC Truth and Muon Reconstruction Analysis
    // Output to tree when muonTruthRecoDiff_fill = 1 in config
    if (muonTruthRecoDiff_fill){
      fPhaseIITrigTree->Branch("deltaVtxX",&fDeltaVtxX,"deltaVtxX/D");
      fPhaseIITrigTree->Branch("deltaVtxY",&fDeltaVtxY,"deltaVtxY/D");
      fPhaseIITrigTree->Branch("deltaVtxZ",&fDeltaVtxZ,"deltaVtxZ/D");
      fPhaseIITrigTree->Branch("deltaVtxR",&fDeltaVtxR,"deltaVtxR/D");
      fPhaseIITrigTree->Branch("deltaVtxT",&fDeltaVtxT,"deltaVtxT/D");
      fPhaseIITrigTree->Branch("deltaParallel",&fDeltaParallel,"deltaParallel/D");
      fPhaseIITrigTree->Branch("deltaPerpendicular",&fDeltaPerpendicular,"deltaPerpendicular/D");
      fPhaseIITrigTree->Branch("deltaAzimuth",&fDeltaAzimuth,"deltaAzimuth/D");
      fPhaseIITrigTree->Branch("deltaZenith",&fDeltaZenith,"deltaZenith/D");
      fPhaseIITrigTree->Branch("deltaAngle",&fDeltaAngle,"deltaAngle/D");
    } 
  }
  return true;
}

bool PhaseIITreeMaker::Execute(){
  Log("===========================================================================================",v_debug,verbosity);
  Log("PhaseIITreeMaker Tool: Executing",v_debug,verbosity);


  // Reset variables
  this->ResetVariables();
  // Get a pointer to the ANNIEEvent Store

  //  If only clean events are built, return true for dirty events
  if(fillCleanEventsOnly){
    auto get_flagsapp = m_data->Stores.at("RecoEvent")->Get("EventFlagApplied",fEventStatusApplied);
    auto get_flags = m_data->Stores.at("RecoEvent")->Get("EventFlagged",fEventStatusFlagged); 
    if(!get_flagsapp || !get_flags) {
      Log("PhaseITreeMaker tool: No Event status applied or flagged bitmask!!", v_error, verbosity);
      return false;	
    }
    // check if event passes the cut
    if((fEventStatusFlagged) != 0) {
      Log("PhaseIITreeMaker Tool: Event was flagged with one of the active cuts.",v_debug, verbosity);
      return true;	
    }
  }



  if(TankClusterProcessing){
    Log("PhaseIITreeMaker Tool: Beginning Tank cluster processing",v_debug,verbosity);
    //bool get_clusters = m_data->Stores.at("ANNIEEvent")->Get("ClusterMap",m_all_clusters);
    bool get_clusters = m_data->CStore.Get("ClusterMap",m_all_clusters);
    if(!get_clusters){
      std::cout << "PhaseIITreeMaker tool: No clusters found!" << std::endl;
      return false;
    }
    Log("PhaseIITreeMaker Tool: Accessing pairs in all_clusters map",v_debug,verbosity);
    get_clusters = m_data->CStore.Get("ClusterMapDetkey",m_all_clusters_detkeys);
    if (!get_clusters){
      std::cout <<"PhaseIITreeMaker tool: No cluster detkeys found!" << std::endl;
      return false;
    }

    int cluster_num = 0;
    std::cout <<"PMT cluster size: "<<m_all_clusters->size()<<std::endl;
    for (std::pair<double,std::vector<Hit>>&& cluster_pair : *m_all_clusters) {
      Log("PhaseIITreeMaker Tool: Resetting variables prior to getting run level info",v_debug,verbosity);
      this->ResetVariables();
      fClusterNumber = cluster_num;

      //Standard run level information
      Log("PhaseIITreeMaker Tool: Getting run level information from ANNIEEvent",v_debug,verbosity);
      m_data->Stores.at("ANNIEEvent")->Get("RunNumber",fRunNumber);
      m_data->Stores.at("ANNIEEvent")->Get("SubrunNumber",fSubrunNumber);
      m_data->Stores.at("ANNIEEvent")->Get("RunType",fRunType);
      m_data->Stores.at("ANNIEEvent")->Get("RunStartTime",fStartTime);
  
      fStartTime_Tree = (ULong64_t) fStartTime;
      // ANNIE Event number
      m_data->Stores.at("ANNIEEvent")->Get("EventTimeTank",fEventTimeTank);
      fEventTimeTank_Tree = (ULong64_t) fEventTimeTank;
      m_data->Stores.at("ANNIEEvent")->Get("EventNumber",fEventNumber);
      uint32_t trigword_tmp;
      m_data->Stores.at("ANNIEEvent")->Get("TriggerWord",trigword_tmp);
      fTriggerword = int(trigword_tmp);
      m_data->Stores["ANNIEEvent"]->Get("TriggerExtended",fExtended);
      BeamStatus beamstat;
      m_data->Stores["ANNIEEvent"]->Get("BeamStatus",beamstat);
      if (beamstat.ok()) fBeamok = 1;
      else fBeamok = 0;
      fPot = beamstat.pot();

      bool pmtmrdcoinc, noveto;
      m_data->Stores.at("RecoEvent")->Get("PMTMRDCoinc",pmtmrdcoinc);
      m_data->Stores.at("RecoEvent")->Get("NoVeto",noveto);
      std::cout <<"Tank cluster processing: PMTMRDCoinc: "<<pmtmrdcoinc<<std::endl;

      if (pmtmrdcoinc) fTankMRDCoinc = 1;
      else fTankMRDCoinc = 0;
      if (noveto) fNoVeto = 1;
      else fNoVeto = 0;

      //Number of ADC samples
      std::map<unsigned long, std::vector<Waveform<unsigned short>>> raw_waveform_map;
      bool has_raw = m_data->Stores["ANNIEEvent"]->Get("RawADCData",raw_waveform_map);
      if (!has_raw) {
        Log("PhaseIITreeMaker tool: Did not find RawADCData in ANNIEEvent! Abort",v_error,verbosity);
        /*return false;*/
      }

      if (has_raw){
        for (auto& temp_pair : raw_waveform_map) {
          const auto& achannel_key = temp_pair.first;
          auto& araw_waveforms = temp_pair.second;
          for (unsigned int i=0; i< araw_waveforms.size(); i++){
            auto samples = araw_waveforms.at(i).GetSamples();
            int size_sample = samples->size();
            fADCWaveformChankeys.push_back(achannel_key);
            fADCWaveformSamples.push_back(size_sample);
          }
        }
      }else {
        //Try to get summary information on the ADC waveform samples in case no raw info is saved
        std::map<unsigned long, std::vector<int>> raw_acqsize_map;
        has_raw = m_data->Stores["ANNIEEvent"]->Get("RawAcqSize",raw_acqsize_map);
        if (!has_raw) {
          Log("RunValidation tool: Did not find RawAcqSize in ANNIEEvent! Abort",v_error,verbosity);
          return false;
        }
        int size_of_window = 2000;
        if (has_raw){
          for (auto& temp_pair : raw_acqsize_map) {
            const auto& achannel_key = temp_pair.first;
            auto& araw_acqsize = temp_pair.second;
            for (unsigned int i=0; i< araw_acqsize.size(); i++){
              int size_sample = araw_acqsize.at(i);
              fADCWaveformChankeys.push_back(achannel_key);
              fADCWaveformSamples.push_back(size_sample);
            }
          }
        }
      }

      std::vector<Hit> cluster_hits = cluster_pair.second;
      fClusterTime = cluster_pair.first;
      if(TankHitInfo_fill){
        Log("PhaseIITreeMaker Tool: Loading tank cluster hits into cluster tree",v_debug,verbosity);
        this->LoadTankClusterHits(cluster_hits);
      }

      bool good_class = this->LoadTankClusterClassifiers(cluster_pair.first);
      if(!good_class){
        if(verbosity>3) Log("PhaseIITreeMaker Tool: No cluster classifiers.  Continuing tree",v_debug,verbosity);
      }
      if(SiPMPulseInfo_fill) this->LoadSiPMHits();
      fPhaseIITankClusterTree->Fill();
      cluster_num += 1;
    }

  }

  if(MRDClusterProcessing){
    Log("PhaseIITreeMaker Tool: Beginning MRD cluster processing",v_debug,verbosity);

    bool get_clusters = m_data->CStore.Get("MrdTimeClusters",MrdTimeClusters);
    if(!get_clusters){
      std::cout << "PhaseIITreeMaker tool: No MRD clusters found! Did you run the TimeClustering tool?" << std::endl;
      return false;
    }
	
    int num_mrd_clusters;
    m_data->CStore.Get("NumMrdTimeClusters",num_mrd_clusters);


    if (num_mrd_clusters > 0){
      m_data->CStore.Get("MrdDigitTimes",mrddigittimesthisevent);
      m_data->CStore.Get("MrdDigitPmts",mrddigitpmtsthisevent);
      m_data->CStore.Get("MrdDigitChankeys",mrddigitchankeysthisevent);
    }
    

    //Check if there's a veto hit in the acquisition
    int TrigHasVetoHit = 0;
    bool has_tdc = m_data->Stores["ANNIEEvent"]->Get("TDCData",TDCData);  // a std::map<ChannelKey,vector<TDCHit>>
    if (!has_tdc){
      std::cout << "No TDCData store in ANNIEEvent." << std::endl;
      //return false;
    }
    
    int cluster_num = 0;
    for (int i=0; i<MrdTimeClusters.size(); i++){
      Log("PhaseIITreeMaker Tool: Resetting variables prior to getting MRD cluster info",v_debug,verbosity);
      this->ResetVariables();
      uint32_t trigword_temp;
      m_data->Stores.at("ANNIEEvent")->Get("TriggerWord",trigword_temp);
      fTriggerword = int(trigword_temp);
      m_data->Stores["ANNIEEvent"]->Get("TriggerExtended",fExtended);
      BeamStatus beamstat;
      m_data->Stores["ANNIEEvent"]->Get("BeamStatus",beamstat);
      if (beamstat.ok()) fBeamok = 1;
      else fBeamok = 0;
      fPot = beamstat.pot();
      bool pmtmrdcoinc, noveto;
      m_data->Stores.at("RecoEvent")->Get("PMTMRDCoinc",pmtmrdcoinc);
      m_data->Stores.at("RecoEvent")->Get("NoVeto",noveto);
	
      if (pmtmrdcoinc) fTankMRDCoinc = 1;
      else fTankMRDCoinc = 0;
      if (noveto) fNoVeto = 1;
      else fNoVeto = 0;
      fMRDClusterHits = 0;
      if (has_tdc){
        if(TDCData->size()>0){
          Log("PhaseIITreeMaker tool: Looping over FACC/MRD hits... looking for Veto activity",v_debug,verbosity);
          for(auto&& anmrdpmt : (*TDCData)){
            unsigned long chankey = anmrdpmt.first;
            std::vector<Hit> mrdhits = anmrdpmt.second;
            Detector* thedetector = geom->ChannelToDetector(chankey);
            unsigned long detkey = thedetector->GetDetectorID();
            if(thedetector->GetDetectorElement()=="Veto") {
              TrigHasVetoHit=1; // this is a veto hit, not an MRD hit.
              for (int j = 0; j < (int) mrdhits.size(); j++){
                fFMVHitT.push_back(mrdhits.at(j).GetTime());
                fFMVHitDetID.push_back(detkey);
                fFMVHitChankey.push_back(chankey);
              }
            }
          }
        }
      }
      fVetoHit = TrigHasVetoHit;
      std::vector<int> ThisClusterIndices = MrdTimeClusters.at(i);
      std::cout <<"Loop through cluster indices"<<std::endl;
      for (int j=0; j < (int) ThisClusterIndices.size(); j++){
        Detector *thedetector = geom->ChannelToDetector(mrddigitchankeysthisevent.at(ThisClusterIndices.at(j)));
        unsigned long detkey = thedetector->GetDetectorID();
        fMRDHitT.push_back(mrddigittimesthisevent.at(ThisClusterIndices.at(j)));
        fMRDHitDetID.push_back(detkey);
        fMRDHitChankey.push_back(mrddigitchankeysthisevent.at(ThisClusterIndices.at(j)));
        fMRDClusterHits+=1;
      }

      fMRDClusterNumber = cluster_num;
      ComputeMeanAndVariance(fMRDHitT, fMRDClusterTime, fMRDClusterTimeSigma);
      //FIXME: calculate fMRDClusterTime

      //Standard run level information
      Log("PhaseIITreeMaker Tool: Getting run level information from ANNIEEvent",v_debug,verbosity);
      m_data->Stores.at("ANNIEEvent")->Get("RunNumber",fRunNumber);
      m_data->Stores.at("ANNIEEvent")->Get("SubrunNumber",fSubrunNumber);
      m_data->Stores.at("ANNIEEvent")->Get("RunType",fRunType);
      m_data->Stores.at("ANNIEEvent")->Get("RunStartTime",fStartTime);
      fStartTime_Tree = (ULong64_t) fStartTime;
      m_data->Stores.at("ANNIEEvent")->Get("EventNumber",fEventNumber);
      m_data->Stores.at("ANNIEEvent")->Get("EventTimeTank",fEventTimeTank);
      fEventTimeTank_Tree = (ULong64_t) fEventTimeTank;
      std::cout <<"MRD EventTimeTank: "<<fEventTimeTank<<std::endl;
      m_data->Stores.at("ANNIEEvent")->Get("EventTimeMRD",fEventTimeMRD);  
      std::cout <<"MRD EventTimeMRD: "<<fEventTimeMRD.GetNs()<<std::endl;
      fEventTimeMRD_Tree = (ULong64_t) fEventTimeMRD.GetNs();

      // ANNIE Event number
      //bool got_mrdtime = m_data->Stores.at("ANNIEEvent")->Get("EventTime",mrd_timestamp);
      //if(got_mrdtime) fEventTimeMRD = static_cast<uint64_t>(mrd_timestamp->GetNs());

      if(MRDReco_fill){
        fNumClusterTracks = this->LoadMRDTrackReco(i);
        //Get the track info
      }

      fPhaseIIMRDClusterTree->Fill();
      cluster_num += 1;
    }
  }
 
  if(TriggerProcessing) {

    Log("PhaseIITreeMaker Tool: Getting run level information from ANNIEEvent",v_debug,verbosity);
    this->ResetVariables();

    m_data->Stores.at("ANNIEEvent")->Get("RunNumber",fRunNumber);
    m_data->Stores.at("ANNIEEvent")->Get("SubrunNumber",fSubrunNumber);
    m_data->Stores.at("ANNIEEvent")->Get("RunType",fRunType);
    m_data->Stores.at("ANNIEEvent")->Get("RunStartTime",fStartTime);
    fStartTime_Tree = (ULong64_t) fStartTime;  

    // ANNIE Event number
    m_data->Stores.at("ANNIEEvent")->Get("EventNumber",fEventNumber);
    m_data->Stores.at("ANNIEEvent")->Get("EventTimeTank",fEventTimeTank);
    fEventTimeTank_Tree = (ULong64_t) fEventTimeTank;

    bool pmtmrdcoinc, noveto;
    uint32_t trigword_temp;
    m_data->Stores.at("ANNIEEvent")->Get("TriggerWord",trigword_temp);
    fTriggerword = int(trigword_temp);
    m_data->Stores["ANNIEEvent"]->Get("TriggerExtended",fExtended);
    BeamStatus beamstat;
    m_data->Stores["ANNIEEvent"]->Get("BeamStatus",beamstat);
    if (beamstat.ok()) fBeamok = 1;
    else fBeamok = 0;
    fPot = beamstat.pot();
    m_data->Stores.at("ANNIEEvent")->Get("DataStreams",fDataStreams);
    m_data->Stores.at("RecoEvent")->Get("PMTMRDCoinc",pmtmrdcoinc);
    m_data->Stores.at("RecoEvent")->Get("NoVeto",noveto);

    if (pmtmrdcoinc) fTankMRDCoinc = 1;
    else fTankMRDCoinc = 0;
    if (noveto) fNoVeto = 1;
    else fNoVeto = 0;

    if (fDataStreams["Tank"]==1) fHasTank = 1;
    else fHasTank = 0;
    if (fDataStreams["MRD"]==1) fHasMRD = 1;
    else fHasMRD = 0;

    m_data->Stores.at("ANNIEEvent")->Get("EventTimeMRD",fEventTimeMRD);
    fEventTimeMRD_Tree = (ULong64_t) fEventTimeMRD.GetNs();
    // bool got_mrdtime = m_data->Stores.at("ANNIEEvent")->Get("EventTime",mrd_timestamp);
    //if(got_mrdtime) fEventTimeMRD = static_cast<uint64_t>(mrd_timestamp->GetNs());

    // Read hits and load into ntuple
    if(TankHitInfo_fill){
      this->LoadAllTankHits();
    }
    if(SiPMPulseInfo_fill) this->LoadSiPMHits();

    if(MRDHitInfo_fill) this->LoadAllMRDHits();

    if(MRDReco_fill){
      fNumClusterTracks=0;
      bool get_clusters = m_data->CStore.Get("MrdTimeClusters",MrdTimeClusters);
      if(!get_clusters){
        std::cout << "PhaseIITreeMaker tool: No MRD clusters found.  Will be no tracks." << std::endl;
        return false;
      }
      for(int i=0; i < (int) MrdTimeClusters.size(); i++) fNumClusterTracks += this->LoadMRDTrackReco(i);
    }

    bool got_reco = false;
    if(TankReco_fill) got_reco = this->FillTankRecoInfo();

    bool gotmctruth = false;
    if(MCTruth_fill)  gotmctruth = this->FillMCTruthInfo();

    if (muonTruthRecoDiff_fill) this->FillTruthRecoDiffInfo(gotmctruth,got_reco);
    if (got_reco && gotmctruth && (verbosity>4)) this->RecoSummary();

    // FIll tree with all reconstruction information
    if (RecoDebug_fill) this->FillRecoDebugInfo();

    fPhaseIITrigTree->Fill();
  }
  return true;
}

bool PhaseIITreeMaker::Finalise(){
  
  fOutput_tfile->cd();
  fPhaseIITrigTree->Write();
  fPhaseIIMRDClusterTree->Write();
  fPhaseIITankClusterTree->Write();
  fOutput_tfile->Close();
  if(verbosity>0) cout<<"PhaseIITreeMaker exitting"<<endl;

  return true;
}

void PhaseIITreeMaker::ResetVariables() {
  // tree variables
  fStartTime_Tree = 9999;
  fEventNumber = -9999;
  fEventTimeTank_Tree = 9999;
  fNHits = -9999;
  fVetoHit = -9999;
  fEventTimeMRD_Tree = 9999;
  fTriggerword = -1;
  fTankMRDCoinc = -99999;
  fNoVeto = -99999;
  fHasTank = -99999;
  fHasMRD = -99999;
  fPot = -99999;
  fBeamok = -99999;
  fExtended = -99999;

  if(SiPMPulseInfo_fill){
    fSiPM1NPulses = -9999;
    fSiPM2NPulses = -9999;
    fSiPMHitQ.clear();
    fSiPMHitT.clear();
    fSiPMHitAmplitude.clear();
    fSiPMNum.clear();
  }

  if(TankClusterProcessing){
    fClusterMaxPE = -9999;
    fClusterChargePointX = -9999;
    fClusterChargePointY = -9999;
    fClusterChargePointZ = -9999;
    fClusterChargeBalance = -9999;
    fClusterTime = -9999;
    fClusterCharge = -9999;
    fClusterNumber = -9999;
    fADCWaveformSamples.clear();
    fADCWaveformChankeys.clear();
  } 
  if(MCTruth_fill){ 
    fMCEventNum = -9999;
    fMCTriggerNum = -9999;
    fTrueVtxX = -9999;
    fTrueVtxY = -9999;
    fTrueVtxZ = -9999;
    fTrueVtxTime = -9999;
    fTrueMuonEnergy = -9999;
    fTrueTrackLengthInWater = -9999;
    fTrueTrackLengthInMRD = -9999;
    fTrueDirX = -9999;
    fTrueDirY = -9999;
    fTrueDirZ = -9999;
    fTrueAngle = -9999;
    fTruePhi = -9999;
    fPi0Count = -9999;
    fPiPlusCount = -9999;
    fPiMinusCount = -9999;
    fK0Count = -9999;
    fKPlusCount = -9999;
    fKMinusCount = -9999;
  }

  if (RecoDebug_fill){ 
    fSeedVtxX.clear();
    fSeedVtxY.clear();
    fSeedVtxZ.clear();
    fSeedVtxFOM.clear();
    fSeedVtxTime = -9999;
    fPointPosX = -9999;
    fPointPosY = -9999;
    fPointPosZ = -9999;
    fPointPosTime = -9999;
    fPointPosFOM = -9999;
    fPointPosStatus = -9999;
    fPointDirX = -9999;
    fPointDirY = -9999;
    fPointDirZ = -9999;
    fPointDirTime = -9999;
    fPointDirFOM = -9999;
    fPointDirStatus = -9999;
    fPointVtxPosX = -9999;
    fPointVtxPosY = -9999;
    fPointVtxPosZ = -9999;
    fPointVtxDirX = -9999;
    fPointVtxDirY = -9999;
    fPointVtxDirZ = -9999;
    fPointVtxTime = -9999;
    fPointVtxStatus = -9999;
    fPointVtxFOM = -9999;
  }

  if(TankReco_fill){
    fRecoVtxX = -9999;
    fRecoVtxY = -9999;
    fRecoVtxZ = -9999;
    fRecoStatus = -9999;
    fRecoVtxTime = -9999;
    fRecoVtxFOM = -9999;
    fRecoDirX = -9999;
    fRecoDirY = -9999;
    fRecoDirZ = -9999;
    fRecoAngle = -9999;
    fRecoPhi = -9999;
  }

  if(MRDClusterProcessing){
    fMRDClusterNumber = -9999;
    fMRDClusterTime = -9999;
    fMRDClusterTimeSigma = -9999;
    fMRDClusterHits = -9999;
  }
  if(MRDHitInfo_fill){
    fMRDHitT.clear();
    fMRDHitDetID.clear();
    fMRDHitChankey.clear();
    fFMVHitT.clear();
    fFMVHitDetID.clear();
    fFMVHitChankey.clear();
  }
  if(MRDReco_fill){
    fMRDTrackAngle.clear();
    fMRDTrackAngleError.clear();
    fMRDPenetrationDepth.clear();
    fMRDTrackLength.clear();
    fMRDEntryPointRadius.clear();
    fMRDEnergyLoss.clear();
    fMRDEnergyLossError.clear();
    fMRDTrackStartX.clear();
    fMRDTrackStartY.clear();
    fMRDTrackStartZ.clear();
    fMRDTrackStopX.clear();
    fMRDTrackStopY.clear();
    fMRDTrackStopZ.clear();
    fMRDSide.clear();
    fMRDStop.clear();
    fMRDThrough.clear();
  }

  if(TankHitInfo_fill){
    fIsFiltered.clear();
    fHitX.clear();
    fHitY.clear();
    fHitZ.clear();
    fHitT.clear();
    fHitQ.clear();
    fHitPE.clear();
    fHitType.clear();
    fHitDetID.clear();
    fHitChankey.clear();
  }
  
  if (muonTruthRecoDiff_fill){ 
    fDeltaVtxX = -9999;
    fDeltaVtxY = -9999;
    fDeltaVtxZ = -9999;
    fDeltaVtxR = -9999;
    fDeltaVtxT = -9999;
    fDeltaParallel = -9999;
    fDeltaPerpendicular = -9999;
    fDeltaAzimuth = -9999;
    fDeltaZenith = -9999;
    fDeltaAngle = -9999;
  }
}

bool PhaseIITreeMaker::LoadTankClusterClassifiers(double cluster_time){
  //Save classifiers to ANNIEEvent
  Log("PhaseITreeMaker tool: Getting cluster classifiers", v_debug, verbosity);

  bool got_ccp = m_data->Stores.at("ANNIEEvent")->Get("ClusterChargePoints", ClusterChargePoints);
  bool got_ccb = m_data->Stores.at("ANNIEEvent")->Get("ClusterChargeBalances", ClusterChargeBalances);
  bool got_cmpe = m_data->Stores.at("ANNIEEvent")->Get("ClusterMaxPEs", ClusterMaxPEs);
  bool good_classifiers = got_ccp && got_ccb && got_cmpe;
  if(!good_classifiers){
    Log("PhaseITreeMaker tool: One of the charge cluster classifiers is not available", v_debug, verbosity);
  } else { 
    Log("PhaseITreeMaker tool: Setting fCluster variables to classifier parameters", v_debug, verbosity);
    fClusterMaxPE = ClusterMaxPEs.at(cluster_time);
    Position ClusterChargePoint = ClusterChargePoints.at(cluster_time);
    fClusterChargePointX = ClusterChargePoint.X();
    fClusterChargePointY = ClusterChargePoint.Y();
    fClusterChargePointZ = ClusterChargePoint.Z();
    fClusterChargeBalance = ClusterChargeBalances.at(cluster_time);
  }
  return good_classifiers;
}

void PhaseIITreeMaker::LoadTankClusterHits(std::vector<Hit> cluster_hits){
  Position detector_center=geom->GetTankCentre();
  double tank_center_x = detector_center.X();
  double tank_center_y = detector_center.Y();
  double tank_center_z = detector_center.Z();

  fClusterCharge = 0;
  fClusterPE = 0;
  fClusterHits = 0;
  for (int i = 0; i<cluster_hits.size(); i++){
    int channel_key = cluster_hits.at(i).GetTubeId();
    std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(channel_key);
    if(it != ChannelKeyToSPEMap.end()){ //Charge to SPE conversion is available
      Detector* this_detector = geom->ChannelToDetector(channel_key);
      Position det_position = this_detector->GetDetectorPosition();
      double hit_charge = cluster_hits.at(i).GetCharge();
      double hit_PE = hit_charge / ChannelKeyToSPEMap.at(channel_key);
      fHitX.push_back((det_position.X()-tank_center_x));
      fHitY.push_back((det_position.Y()-tank_center_y));
      fHitZ.push_back((det_position.Z()-tank_center_z));
      fHitQ.push_back(hit_charge);
      fHitPE.push_back(hit_PE);
      fHitT.push_back(cluster_hits.at(i).GetTime());
      fHitDetID.push_back(channel_key);
      fHitChankey.push_back(channel_key);
      fHitType.push_back(RecoDigit::PMT8inch);
      fClusterCharge+=hit_charge;
      fClusterPE+=hit_PE;
      fClusterHits += 1;
    } else {
      if(verbosity>4){
        std::cout << "FOUND A HIT FOR CHANNELKEY " << channel_key << "BUT NO CONVERSION " <<
            "TO PE AVAILABLE.  SKIPPING PE." << std::endl;
      }
    }
  }
  return;
}


void PhaseIITreeMaker::LoadSiPMHits() {
  std::map<unsigned long, std::vector< std::vector<ADCPulse>> > aux_pulse_map;
  m_data->Stores.at("ANNIEEvent")->Get("RecoADCAuxHits", aux_pulse_map);
  fSiPM1NPulses = 0;
  fSiPM2NPulses = 0;
  for (const auto& temp_pair : aux_pulse_map) {
    const auto& channel_key = temp_pair.first;
    //For now, only calibrate the SiPM waveforms
    int sipm_number = -1;
    if(AuxChannelNumToTypeMap->at(channel_key) == "SiPM1"){
      sipm_number = 1;
    } else if (AuxChannelNumToTypeMap->at(channel_key) == "SiPM2"){
      sipm_number = 2;
    } else continue;

    std::vector< std::vector<ADCPulse>> sipm_minibuffers = temp_pair.second;
    size_t num_minibuffers = sipm_minibuffers.size();  //Should be size 1 in FrankDAQ mode
    for (size_t mb = 0; mb < num_minibuffers; ++mb) {
      std::vector<ADCPulse> thisbuffer_pulses = sipm_minibuffers.at(mb);
      if(sipm_number == 1) fSiPM1NPulses += thisbuffer_pulses.size();
      if(sipm_number == 2) fSiPM2NPulses += thisbuffer_pulses.size();
      for (size_t i = 0; i < thisbuffer_pulses.size(); i++){
        ADCPulse apulse = thisbuffer_pulses.at(i);
        fSiPMHitAmplitude.push_back(apulse.amplitude());
        fSiPMHitT.push_back(apulse.peak_time());
        fSiPMHitQ.push_back(apulse.charge());
        fSiPMNum.push_back(sipm_number);
      }
    }
  }
}

void PhaseIITreeMaker::LoadAllMRDHits(){
  //bool get_clusters = m_data->CStore.Get("MrdTimeClusters",MrdTimeClusters);
  //if(!get_clusters){
  //  std::cout << "BeamClusterAnalysis tool: No MRD clusters found! Did you run the TimeClustering tool?" << std::endl;
  //  return;
  //}
  
  //m_data->CStore.Get("MrdDigitTimes",mrddigittimesthisevent);
  //m_data->CStore.Get("MrdDigitPmts",mrddigitpmtsthisevent);

  fVetoHit = 0; 
  bool get_ok = m_data->Stores["ANNIEEvent"]->Get("TDCData",TDCData);  // a std::map<ChannelKey,vector<TDCHit>>
  if (!get_ok){
    std::cout << "PhaseIITreeMaker tool: No TDCData store in ANNIEEvent." << std::endl;
    return;
  }
  if(TDCData->size()>0){
    Log("PhaseIITreeMaker tool: Looping over FACC/MRD hits... looking for Veto activity",v_debug,verbosity);
    for(auto&& anmrdpmt : (*TDCData)){
      unsigned long chankey = anmrdpmt.first;
      Detector* thedetector = geom->ChannelToDetector(chankey);
      unsigned long detkey = thedetector->GetDetectorID();
      std::vector<Hit> mrdhits = anmrdpmt.second;
      if(thedetector->GetDetectorElement()=="Veto") {
        fVetoHit=1; // this is a veto hit, not an MRD hit.
        for (int j = 0; j < (int) mrdhits.size(); j++){
          fFMVHitT.push_back(mrdhits.at(j).GetTime());
          fFMVHitDetID.push_back(detkey);
          fFMVHitChankey.push_back(chankey);
        }
      } else {
        for(int j = 0; j<mrdhits.size(); j++){
          fMRDHitT.push_back(mrdhits.at(j).GetTime());
          fMRDHitDetID.push_back(detkey);
          fMRDHitChankey.push_back(chankey);
        }
      }
    }
  }
  return;
}

int PhaseIITreeMaker::LoadMRDTrackReco(int SubEventID) {
  //Check for valid track criteria
  m_data->Stores["MRDTracks"]->Get("MRDTracks",theMrdTracks);
  m_data->Stores["MRDTracks"]->Get("NumMrdTracks",numtracksinev);
  // Loop over reconstructed tracks


  Position StartVertex;
  Position StopVertex;
  double TrackLength = -9999;
  double TrackAngle = -9999;
  double TrackAngleError = -9999;
  double PenetrationDepth = -9999;
  Position MrdEntryPoint;
  double EnergyLoss = -9999; //in MeV
  double EnergyLossError = -9999;
  double EntryPointRadius = -9999;
  bool IsMrdPenetrating;
  bool IsMrdStopped;
  bool IsMrdSideExit;

  int NumClusterTracks = 0;
  for(int tracki=0; tracki<numtracksinev; tracki++){
    BoostStore* thisTrackAsBoostStore = &(theMrdTracks->at(tracki));
    int TrackEventID = -1; 
    //get track properties that are needed for the through-going muon selection
	thisTrackAsBoostStore->Get("MrdSubEventID",TrackEventID);
    if(TrackEventID!= SubEventID) continue;
 

    //If we're here, this track is associated with this cluster
    thisTrackAsBoostStore->Get("StartVertex",StartVertex);
    thisTrackAsBoostStore->Get("StopVertex",StopVertex);
    thisTrackAsBoostStore->Get("TrackAngle",TrackAngle);
    thisTrackAsBoostStore->Get("TrackAngleError",TrackAngleError);
    thisTrackAsBoostStore->Get("PenetrationDepth",PenetrationDepth);
    thisTrackAsBoostStore->Get("MrdEntryPoint",MrdEntryPoint);
    thisTrackAsBoostStore->Get("EnergyLoss",EnergyLoss);
    thisTrackAsBoostStore->Get("EnergyLossError",EnergyLossError);
    thisTrackAsBoostStore->Get("IsMrdPenetrating",IsMrdPenetrating);        // bool
    thisTrackAsBoostStore->Get("IsMrdStopped",IsMrdStopped);                // bool
    thisTrackAsBoostStore->Get("IsMrdSideExit",IsMrdSideExit);
    TrackLength = sqrt(pow((StopVertex.X()-StartVertex.X()),2)+pow(StopVertex.Y()-StartVertex.Y(),2)+pow(StopVertex.Z()-StartVertex.Z(),2)) * 100.0;
    EntryPointRadius = sqrt(pow(MrdEntryPoint.X(),2) + pow(MrdEntryPoint.Y(),2)) * 100.0; // convert to cm
    PenetrationDepth = PenetrationDepth*100.0;


    //Push back some properties
    fMRDTrackAngle.push_back(TrackAngle);
    fMRDTrackAngleError.push_back(TrackAngleError);
    fMRDTrackLength.push_back(TrackLength);
    fMRDPenetrationDepth.push_back(PenetrationDepth);
    fMRDEntryPointRadius.push_back(EntryPointRadius);
    fMRDEnergyLoss.push_back(EnergyLoss);
    fMRDEnergyLossError.push_back(EnergyLossError);
    fMRDTrackStartX.push_back(StartVertex.X());
    fMRDTrackStartY.push_back(StartVertex.Y());
    fMRDTrackStartZ.push_back(StartVertex.Z());
    fMRDTrackStopX.push_back(StopVertex.X());
    fMRDTrackStopY.push_back(StopVertex.Y());
    fMRDTrackStopZ.push_back(StopVertex.Z());
    fMRDStop.push_back(IsMrdStopped);
    fMRDSide.push_back(IsMrdSideExit);
    fMRDThrough.push_back(IsMrdPenetrating);
    NumClusterTracks+=1;
  }
  return NumClusterTracks;
}

void PhaseIITreeMaker::LoadAllTankHits() {
  std::map<unsigned long, std::vector<Hit>>* Hits = nullptr;
  bool got_hits = m_data->Stores["ANNIEEvent"]->Get("Hits", Hits);
  if (!got_hits){
    std::cout << "No Hits store in ANNIEEvent. Continuing to build tree " << std::endl;
    return;
  }
  Position detector_center=geom->GetTankCentre();
  double tank_center_x = detector_center.X();
  double tank_center_y = detector_center.Y();
  double tank_center_z = detector_center.Z();
  fNHits = 0;
  for(std::pair<unsigned long, std::vector<Hit>>&& apair : *Hits){
    unsigned long channel_key = apair.first;
    std::map<int, double>::iterator it = ChannelKeyToSPEMap.find(channel_key);
    if(it != ChannelKeyToSPEMap.end()){ //Charge to SPE conversion is available
      std::vector<Hit>& ThisPMTHits = apair.second;
      fNHits+=ThisPMTHits.size();
      for (Hit &ahit : ThisPMTHits){
        double hit_charge = ahit.GetCharge();
        double hit_PE  = hit_charge / ChannelKeyToSPEMap.at(channel_key);
        Detector* this_detector = geom->ChannelToDetector(channel_key);
        Position det_position = this_detector->GetDetectorPosition();
        fHitX.push_back((det_position.X()-tank_center_x));
        fHitY.push_back((det_position.Y()-tank_center_y));
        fHitZ.push_back((det_position.Z()-tank_center_z));
        fHitT.push_back(ahit.GetTime());
        fHitQ.push_back(hit_charge);
        fHitPE.push_back(hit_PE);
        fHitDetID.push_back(ahit.GetTubeId());
        fHitChankey.push_back(channel_key);
        fHitType.push_back(RecoDigit::PMT8inch); // 0 For PMTs
      }
    }
  }
  return;
}

bool PhaseIITreeMaker::FillTankRecoInfo() {
  bool got_reco_info = true;
  auto* reco_event = m_data->Stores["RecoEvent"];
  if (!reco_event) {
    Log("Error: The PhaseITreeMaker tool could not find the RecoEvent Store", v_error, verbosity);
    got_reco_info = false;
  }
  // Read reconstructed Vertex
  RecoVertex* recovtx = 0;
  auto get_extendedvtx = m_data->Stores.at("RecoEvent")->Get("ExtendedVertex",recovtx); 
  if(!get_extendedvtx) {
    Log("Warning: The PhaseITreeMaker tool could not find ExtendedVertex. Continuing to build tree", v_message, verbosity);
    got_reco_info = false;
  } else {
    fRecoVtxX = recovtx->GetPosition().X();
    fRecoVtxY = recovtx->GetPosition().Y();
    fRecoVtxZ = recovtx->GetPosition().Z();
    fRecoVtxTime = recovtx->GetTime();
    fRecoVtxFOM = recovtx->GetFOM();
    fRecoDirX = recovtx->GetDirection().X();
    fRecoDirY = recovtx->GetDirection().Y();
    fRecoDirZ = recovtx->GetDirection().Z();
    fRecoAngle = TMath::ACos(fRecoDirZ);
    if (fRecoDirX>0.0){
      fRecoPhi = atan(fRecoDirY/fRecoDirX);
    }
    if (fRecoDirX<0.0){
      fRecoPhi = atan(fRecoDirY/fRecoDirX);
      if( fRecoDirY>0.0) fRecoPhi += TMath::Pi();
      if( fRecoDirY<=0.0) fRecoPhi -= TMath::Pi();
    }
    if (fRecoDirX==0.0){
      if( fRecoDirY>0.0) fRecoPhi = 0.5*TMath::Pi();
      else if( fRecoDirY<0.0) fRecoPhi = -0.5*TMath::Pi();
      else fRecoPhi = 0;
    }
    fRecoStatus = recovtx->GetStatus();
  }
  return got_reco_info;
}

void PhaseIITreeMaker::FillRecoDebugInfo() {
  // Read Seed candidates   
  std::vector<RecoVertex>* seedvtxlist = 0;
  auto get_seedvtxlist = m_data->Stores.at("RecoEvent")->Get("vSeedVtxList",seedvtxlist);  ///> Get List of seeds from "RecoEvent"
  if(get_seedvtxlist){
    for( auto& seed : *seedvtxlist ){
      fSeedVtxX.push_back(seed.GetPosition().X());
      fSeedVtxY.push_back(seed.GetPosition().Y());
      fSeedVtxZ.push_back(seed.GetPosition().Z());
      fSeedVtxTime = seed.GetTime();
    }
  } else {  
	  Log("PhaseIITreeMaker  Tool: No Seed List found.  Continuing to build tree ",v_message,verbosity); 
  }
  std::vector<double>* seedfomlist = 0;
  auto get_seedfomlist = m_data->Stores.at("RecoEvent")->Get("vSeedFOMList",seedfomlist);  ///> Get List of seed FOMs from "RecoEvent" 
  if(get_seedfomlist){
    for( auto& seedFOM : *seedfomlist ){
      fSeedVtxFOM.push_back(seedFOM);
    }
  } else {  
	  Log("PhaseIITreeMaker  Tool: No Seed FOM List found.  Continuing to build tree ",v_message,verbosity); 
  }
  
  // Read PointPosition-fitted Vertex   
  RecoVertex* pointposvtx = 0;
  auto get_pointposdata = m_data->Stores.at("RecoEvent")->Get("PointPosition",pointposvtx);
  if(get_pointposdata){ 
    fPointPosX = pointposvtx->GetPosition().X();
    fPointPosY = pointposvtx->GetPosition().Y();
    fPointPosZ = pointposvtx->GetPosition().Z();
    fPointPosTime = pointposvtx->GetTime();
    fPointPosFOM = pointposvtx->GetFOM();
    fPointPosStatus = pointposvtx->GetStatus();
  } else{
    Log("PhaseIITreeMaker Tool: No PointPosition Tool data found.  Continuing to build remaining tree",v_message,verbosity);
  }

  // Read PointDirection-fitted Vertex
  RecoVertex* pointdirvtx = 0;
  auto get_pointdirdata = m_data->Stores.at("RecoEvent")->Get("PointDirection",pointdirvtx);
  if(get_pointdirdata){ 
    fPointDirX = pointdirvtx->GetDirection().X();
    fPointDirY = pointdirvtx->GetDirection().Y();
    fPointDirZ = pointdirvtx->GetDirection().Z();
    fPointDirTime = pointdirvtx->GetTime();
    fPointDirFOM = pointdirvtx->GetFOM();
    fPointDirStatus = pointdirvtx->GetStatus();
  } else{
    Log("PhaseIITreeMaker Tool: No PointDirection Tool data found.  Continuing to build remaining tree",v_message,verbosity);
  }
  
  // Read PointVertex Tool's fitted Vertex
  RecoVertex* pointvtx = 0;
  auto get_pointvtxdata = m_data->Stores.at("RecoEvent")->Get("PointVertex",pointvtx);
  if(get_pointvtxdata){ 
    fPointVtxPosX = pointvtx->GetPosition().X();
    fPointVtxPosY = pointvtx->GetPosition().Y();
    fPointVtxPosZ = pointvtx->GetPosition().Z();
    fPointVtxDirX = pointvtx->GetDirection().X();
    fPointVtxDirY = pointvtx->GetDirection().Y();
    fPointVtxDirZ = pointvtx->GetDirection().Z();
    fPointVtxTime = pointvtx->GetTime();
    fPointVtxFOM = pointvtx->GetFOM();
    fPointVtxStatus = pointvtx->GetStatus();
  } else{
    Log("PhaseIITreeMaker Tool: No PointVertex Tool data found.  Continuing to build remaining tree",v_message,verbosity);
  }
}

bool PhaseIITreeMaker::FillMCTruthInfo() {
  bool successful_load = true;
  // MC entry number
  m_data->Stores.at("ANNIEEvent")->Get("MCEventNum",fMCEventNum);  
  // MC trigger number
  m_data->Stores.at("ANNIEEvent")->Get("MCTriggernum",fMCTriggerNum); 
  std::string logmessage = "  Retriving information for MCEntry "+to_string(fMCEventNum)+
  	", MCTrigger "+ to_string(fMCTriggerNum) + ", EventNumber " + to_string(fEventNumber);
  Log(logmessage, v_message, verbosity);
 


  RecoVertex* truevtx = 0;
  auto get_muonMC = m_data->Stores.at("RecoEvent")->Get("TrueVertex",truevtx);
  auto get_muonMCEnergy = m_data->Stores.at("RecoEvent")->Get("TrueMuonEnergy",fTrueMuonEnergy);
  if(get_muonMC && get_muonMCEnergy){ 
    fTrueVtxX = truevtx->GetPosition().X();
    fTrueVtxY = truevtx->GetPosition().Y();
    fTrueVtxZ = truevtx->GetPosition().Z();
    fTrueVtxTime = truevtx->GetTime();
    fTrueDirX = truevtx->GetDirection().X();
    fTrueDirY = truevtx->GetDirection().Y();
    fTrueDirZ = truevtx->GetDirection().Z();
    double TrueAngRad = TMath::ACos(fTrueDirZ);
    fTrueAngle = TrueAngRad/(TMath::Pi()/180.0); // radians->degrees
    if (fTrueDirX>0.0){
      fTruePhi = atan(fTrueDirY/fTrueDirX);
    }
    if (fTrueDirX<0.0){
      fTruePhi = atan(fTrueDirY/fTrueDirX);
      if( fTrueDirY>0.0) fTruePhi += TMath::Pi();
      if( fTrueDirY<=0.0) fTruePhi -= TMath::Pi();
    }
    if (fTrueDirX==0.0){
      if( fTrueDirY>0.0) fTruePhi = 0.5*TMath::Pi();
      else if( fTrueDirY<0.0) fTruePhi = -0.5*TMath::Pi();
      else fTruePhi = 0;
    }
  } else {
    Log("PhaseIITreeMaker Tool: Missing MC Energy/Vertex info; is this MC?  Continuing to build remaining tree",v_message,verbosity);
    successful_load = false;
  }
  double waterT, MRDT;
  auto get_tankTrackLength = m_data->Stores.at("RecoEvent")->Get("TrueTrackLengthInWater",waterT); 
  auto get_MRDTrackLength = m_data->Stores.at("RecoEvent")->Get("TrueTrackLengthInMRD",MRDT); 
  if (get_tankTrackLength && get_MRDTrackLength){
    fTrueTrackLengthInWater = waterT;
    fTrueTrackLengthInMRD = MRDT;
  } else {
    Log("PhaseIITreeMaker Tool: True track lengths missing. Continuing to build tree",v_message,verbosity);
    successful_load = false;
  }

  int pi0count, pipcount, pimcount, K0count, Kpcount, Kmcount;
  auto get_pi0 = m_data->Stores.at("RecoEvent")->Get("MCPi0Count",pi0count);
  auto get_pim = m_data->Stores.at("RecoEvent")->Get("MCPiMinusCount",pimcount);
  auto get_pip = m_data->Stores.at("RecoEvent")->Get("MCPiPlusCount",pipcount);
  auto get_k0 = m_data->Stores.at("RecoEvent")->Get("MCK0Count",K0count);
  auto get_km = m_data->Stores.at("RecoEvent")->Get("MCKMinusCount",Kmcount);
  auto get_kp = m_data->Stores.at("RecoEvent")->Get("MCKPlusCount",Kpcount);
  if(get_pi0 && get_pim && get_pip && get_k0 && get_km && get_kp) {
    // set values in tree to thouse grabbed from the RecoEvent Store
    fPi0Count = pi0count;
    fPiPlusCount = pipcount;
    fPiMinusCount = pimcount;
    fK0Count = K0count;
    fKPlusCount = Kpcount;
    fKMinusCount = Kmcount;
  } else {
    Log("PhaseIITreeMaker Tool: Missing MC Pion/Kaon count information. Continuing to build remaining tree",v_message,verbosity);
    successful_load = false;
  }
  return successful_load;
}

void PhaseIITreeMaker::FillTruthRecoDiffInfo(bool successful_mcload,bool successful_recoload) {
  if (!successful_mcload || !successful_recoload) {
    Log("PhaseIITreeMaker Tool: Error loading True Muon Vertex or Extended Vertex information.  Continuing to build remaining tree",v_message,verbosity);
  } else {
    //Make sure MCTruth Information is loaded from store
    //Let's fill in stuff from the RecoSummary
    fDeltaVtxX = fRecoVtxX - fTrueVtxX;
    fDeltaVtxY = fRecoVtxY - fTrueVtxY;
    fDeltaVtxZ = fRecoVtxZ - fTrueVtxZ;
    fDeltaVtxT = fRecoVtxTime - fTrueVtxTime;
    fDeltaVtxR = sqrt(pow(fDeltaVtxX,2) + pow(fDeltaVtxY,2) + pow(fDeltaVtxZ,2)); 
    fDeltaParallel = fDeltaVtxX*fRecoDirX + fDeltaVtxY*fRecoDirY + fDeltaVtxZ*fRecoDirZ;
    fDeltaPerpendicular = sqrt(pow(fDeltaVtxR,2) - pow(fDeltaParallel,2));
    fDeltaAzimuth = (fRecoAngle - fTrueAngle)/(TMath::Pi()/180.0);
    fDeltaZenith = (fRecoPhi - fTruePhi)/(TMath::Pi()/180.0); 
    double cosphi = fTrueDirX*fRecoDirX+fTrueDirY*fRecoDirY+fTrueDirZ*fRecoDirZ;
    double phi = TMath::ACos(cosphi); // radians
    double TheAngle = phi/(TMath::Pi()/180.0); // radians->degrees
    fDeltaAngle = TheAngle;
  }
}

void PhaseIITreeMaker::RecoSummary() {

  // get reconstruction output
  double dx = fRecoVtxX - fTrueVtxX;
  double dy = fRecoVtxY - fTrueVtxY;
  double dz = fRecoVtxZ - fTrueVtxZ;
  double dt = fRecoVtxTime - fTrueVtxTime;
  double deltaR = sqrt(dx*dx + dy*dy + dz*dz);
  double cosphi = 0., phi = 0., DeltaAngle = 0.;
  cosphi = fTrueDirX*fRecoDirX+fTrueDirY*fRecoDirY+fTrueDirZ*fRecoDirZ;
  phi = TMath::ACos(cosphi); // radians
  DeltaAngle = phi/(TMath::Pi()/180.0); // radians->degrees
  std::cout << "============================================================================"<<std::endl;
  std::cout << " Event number " << fEventNumber << std::endl;
  std::cout << "  trueVtx=(" << fTrueVtxX << ", " << fTrueVtxY << ", " << fTrueVtxZ <<", "<< fTrueVtxTime<< std::endl
            << " TrueMuonEnergy= " << fTrueMuonEnergy << std::endl
            << "           " << fTrueDirX << ", " << fTrueDirY << ", " << fTrueDirZ << ") " << std::endl;
  std::cout << "  recoVtx=(" << fRecoVtxX << ", " << fRecoVtxY << ", " << fRecoVtxZ <<", "<< fRecoVtxTime << std::endl
            << "           " << fRecoDirX << ", " << fRecoDirY << ", " << fRecoDirZ << ") " << std::endl;
  std::cout << "  DeltaR = "<<deltaR<<"[cm]"<<"\t"<<"  DeltaAngle = "<<DeltaAngle<<" [degree]"<<std::endl;
  std::cout << "  FOM = " << fRecoVtxFOM << std::endl;
  std::cout << "  RecoStatus = " << fRecoStatus <<std::endl;
  std::cout << std::endl;
}

