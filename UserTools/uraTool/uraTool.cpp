// standard library includes
#include <ctime>

// ToolAnalysis includes
#include "uraTool.h"

uraTool::uraTool():Tool(){}


bool uraTool::Initialise(std::string configfile, DataModel &data){

  ///////////////////////// Useful header /////////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  ////////////////////////////////////////////////////////////////

  // DEFAULT VARIABLE VALUES
  verbosity = 3;
  outp_file_prefix = "test_";
  beam_evt_ctr = 0; // if defined here, increments as expected

  // Load values in config file
  m_variables.Get("verbosity", verbosity);
  m_variables.Get("OutputFile", outp_file_prefix);

  std::cout << "  Config variables loaded. Moving on.." << std::endl;

  // Load in other shit
//  m_data->CStore.Get("AuxChannelNumToTypeMap", AuxChannelNumToTypeMap);
  auto get_geo_ok = m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry", fGeo);
  if (!get_geo_ok)
  {
    std::cout << " warning! could not load AnnieGeometry!" << std::endl;
  }

  return true;
}


bool uraTool::Execute(){
  // do da real shit hurrr
  if (verbosity > 3)
  {
    std::cout << "MUCH VERBOSE! SO WORDS! MANY TALKS!" << std::endl;
    std::cout << "LET'S DO SOME QUICK MATHS!!!" << std::endl;
  }

  // quick maths
  //int a = 3, b = 9;
  //std::cout << "a + b = " << (a + b) << std::endl;
  //std::cout << "a - b = " << (a - b) << std::endl;

  // check if ANNIEEvent exists
  int annieeventexists = m_data->Stores.count("ANNIEEvent");
  if (!annieeventexists)
  {
    std::cerr << "no ANNIEEvent store!" << std::endl;
    return false;
  }

  // ptr to ANNIEEventStore
  auto* annie_event = m_data->Stores["ANNIEEvent"];
//  annie_event->Print(false);

  // geometry stuff
//  Geometry *fGeo = nullptr; // now in header
//  auto get_geometry = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry", fGeo); // now in Initialise
  bool get_ok = m_data->Stores["ANNIEEvent"]->Get("RecoADCData", pulse_map);
  if (!get_ok)
  {
    std::cout << " not ok! couldn't find RecoADCData in store!" << std::endl;
    return false;
  }
  for (const auto& temp_pair : pulse_map)
  {
    //const unsigned long& chankey = temp_pair.first;
    unsigned long chankey = temp_pair.first;
    auto& pulses = temp_pair.second;

//    std::cout << " >>>>>>>>>>>> chankey: " << chankey << std::endl;
    std::cout << " >>>>>>>>>>>> chankey, pulses size: " << chankey << ", " << pulses.size() << std::endl;

    Detector* thedetector = fGeo->ChannelToDetector(chankey); // where is chankey?
    std::cout << " >>>>>>>>>>>> chankey label: " << thedetector->GetDetectorElement() << std::endl;
  }

  // get evt numnum
  uint32_t nevt;
  m_data->Stores["ANNIEEvent"]->Get("EventNumber", nevt);
  std::cout << " >>>>>>>>>>>>>>>>>>>> evt no. : " << nevt << std::endl;
  
  // get the (mrd) trigger name/type
  std::string trigselect = "Beam";
  std::string trigname, mrdtrigtype;
  std::vector<TriggerClass>* TriggerData = nullptr;
  TriggerClass TriggerDataData;

  // note: some ppl write `Stores.at("ANNIEEvent")` instead. diff?
  m_data->Stores["ANNIEEvent"]->Get("TriggerData", TriggerDataData);
  std::cout << " >>>>>>>>>>>>>>>>>>>> Trigger Data: " << std::endl;
  TriggerDataData.Print();

  trigname = TriggerDataData.GetName();
  std::cout << " >>>>>>>>>>>>>>>>>>>> Trigger Name: " << trigname << std::endl;

  m_data->Stores["ANNIEEvent"]->Get("MRDTriggerType", mrdtrigtype);
  std::cout << " >>>>>>>>>>>>>>>>>>>> MRD Trigger Type: " << mrdtrigtype << std::endl;
  std::cout << std::endl;

  //int beam_evt_ctr = 0; // if defined here, will always print 1
  if (mrdtrigtype == trigselect)
  {
    beam_evt_ctr++;
    if (verbosity > 3) { std::cout << " FOUND A BEAM EVENT BITCHES! " << std::endl; }
    std::cout << " beam events: " << beam_evt_ctr << std::endl; 
  }

  // BEAM STUFF
  bool get_beamstatus;
  // BeamStatusData defined in header
  //BeamCondition beam_condition;
  m_data->Stores["ANNIEEvent"]->Get("BeamCondition", beam_condition);
  std::cout << "beam_condition: " << beam_condition << std::endl;
  
  get_beamstatus = m_data->Stores["ANNIEEvent"]->Get("BeamStatus", BeamStatusData);
  if (get_beamstatus)
  {
    std::cout << "OH SHIIIIIT. CHECK OUT DIS BEAM INFO YO: " << BeamStatusData.Print() << std::endl;
  }
  else
  {
    std::cout << "Oh poop. Didn't find no BeamStatus in ANNIEEvent :'(" << std::endl;
  }

  // aux pulse stuff (whatever that means..)
/*
  annie_event->Get("RecoADCAuxHits", aux_pulse_map);
  std::cout << "  about to hook the aux cord!!! " << std::endl;
  for (const auto& temp_pair : aux_pulse_map)
  {
    // not entering for loop for some reason...
    std::cout << "  made it into the aux!" << std::endl;
    const auto& channel_key = temp_pair.first;
    std::string channel_type = AuxChannelNumToTypeMap->at(channel_key);

    std::cout << "channel (key, type): (" << channel_key << ", " << channel_type << ")" << std::endl;
  }
*/

  return true;
}


bool uraTool::Finalise(){

  this->HelloWorld(12, -6);
  return true;
}

// define local functions
void uraTool::HelloWorld(int a, int b)
{
  int c = a*b;
  std::cout << "c = " << c << std::endl;
  std::cout << "Heywhatsupeverybodywelcometomyyoutubechannel!" << std::endl;
  std::cout << "IT'S PICKLE RIIIIIIICK!!! WUBBALUBBADUBDUB" << std::endl;
  std::cout << " total beam events: " << beam_evt_ctr << std::endl;
  return;
}
