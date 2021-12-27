#include "TimeDistributions.h"

TimeDistributions::TimeDistributions():Tool(){}


bool TimeDistributions::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  timedist_out = new TFile("timedistributions.root","RECREATE");

  return true;
}


bool TimeDistributions::Execute(){

  std::map<unsigned long, std::vector<Hit>> *Hits = nullptr;
  m_data->Stores["ANNIEEvent"]->Get("Hits",Hits);
  int evnum;
  m_data->Stores["ANNIEEvent"]->Get("EventNumber",evnum);
  int extended;
  m_data->Stores["ANNIEEvent"]->Get("TriggerExtended",extended);

  if (extended > 0){
  std::stringstream ss_hist;
  ss_hist << "timedist_ev"<<evnum;
  TH1F *hist = new TH1F(ss_hist.str().c_str(),ss_hist.str().c_str(),1000,0,70000);

  for (std::pair<unsigned long, std::vector<Hit>>&& apair : *Hits){
    std::vector<Hit>& thishit = apair.second;
    for (int i_hit=0; i_hit < thishit.size(); i_hit++){
      hist->Fill(thishit.at(i_hit).GetTime());
    }
  }
  timedist_out->cd();
  hist->Write();
}
  return true;
}


bool TimeDistributions::Finalise(){

  timedist_out->Close();
  delete timedist_out;

  return true;
}
