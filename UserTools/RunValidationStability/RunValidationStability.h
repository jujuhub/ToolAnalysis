#ifndef RunValidationStability_H
#define RunValidationStability_H

#include <string>
#include <iostream>
#include <fstream>

#include "Tool.h"

#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TFitResult.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TF1.h"

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>

/**
 * \class RunValidationStability
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: M. Nieslony $
* $Date: 2021/12/13 14:40:00 $
* Contact: mnieslon@uni-mainz.de
*/
class RunValidationStability: public Tool {


 public:

  RunValidationStability(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  void ConvertDateToMSec(std::string time_str, uint64_t &time_long);
  void ConvertWeekToRuns(uint64_t timestamp_start, int &run_start, int &run_end);

 private:


   int verbosity;
   std::string input_directory;
   std::string output_file;
   int current_run = -1;
   int number_runs = 50;
   int reference_run = 3082;
   bool WeeklyMode = false;
   std::string WeekStartDate;
   std::string WeekStartDateFile;

   std::string str_week_start;
   std::string str_week_end;
   std::string string_week_start;
   std::string string_week_end;
   std::string OutputPath;
   bool SaveImages=false;

};


#endif
