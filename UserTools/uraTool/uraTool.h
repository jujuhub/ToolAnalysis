#ifndef uraTool_H
#define uraTool_H

#include <string>
#include <iostream>
#include <sstream>

#include "Tool.h"
#include "ADCPulse.h"
#include "BeamStatus.h"
#include "BeamStatusClass.h"
#include "Hit.h"
#include "Particle.h"
#include "TriggerClass.h"


/**
 * \class uraTool
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $$$  Author: J. He               $$$
* $$$  Date: 2021/03/03 14:41:00   $$$
* $$$  Contact: uratool@tools.com  $$$
*/
class uraTool: public Tool {


 public:

  uraTool(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  // my local functions
  void HelloWorld(int a, int b);

 private:
  int verbosity;
  std::string outp_file_prefix;

  int beam_evt_ctr = 0;
  BeamStatus BeamStatusData;
  BeamCondition beam_condition;

  Geometry *fGeo = nullptr;

  std::map<unsigned long, std::vector<std::vector<ADCPulse>>> pulse_map;
  std::map<unsigned long, std::vector<std::vector<ADCPulse>>> aux_pulse_map;
  std::map<int, std::string>* AuxChannelNumToTypeMap;

};


#endif
