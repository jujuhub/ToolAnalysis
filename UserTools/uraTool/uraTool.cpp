// standard library includes
#include <ctime>

// ToolAnalysis includes
#include "uraTool.h"

uraTool::uraTool():Tool(){}


bool uraTool::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // Default variable values
  verbosity = 3;
  outp_file_prefix = "test_";

  // Load values in config file
  m_variables.Get("verbosity", verbosity);
  m_variables.Get("OutputFile", outp_file_prefix);

  std::cout << "Config variables loaded. Moving on.." << std::endl;
  return true;
}


bool uraTool::Execute(){
  // Do things here
  if (verbosity > 3)
  {
    std::cout << "MUCH VERBOSE! SO WORDS! MANY TALKS!" << std::endl;
    std::cout << "LET'S DO SOME QUICK MATHS!!!" << std::endl;
  }

  // test
  //int a = 3, b = 9;
  //std::cout << "a + b = " << (a + b) << std::endl;
  //std::cout << "a - b = " << (a - b) << std::endl;

  // print some data variables
  int annieeventexists = m_data->Stores.count("ANNIEEvent");
  if (!annieeventexists)
  {
    std::cerr << "no ANNIEEvent store!" << std::endl;
    return false;
  }

  uint32_t nevt;
  m_data->Stores.at("ANNIEEvent")->Get("EventNumber", nevt);
  std::cout << " >>>>>>>>>>>>>>>>>>>> evt no. : " << nevt << std::endl;
  
  std::string trigname;
  std::vector<TriggerClass>* TriggerData = nullptr;
  TriggerClass TriggerDataData;
  m_data->Stores.at("ANNIEEvent")->Get("TriggerData", TriggerDataData);
  std::cout << " >>>>>>>>>>>>>>>>>>>> Trigger Data: " << std::endl;
  TriggerDataData.Print();
  std::cout << std::endl;

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
  return;
}
