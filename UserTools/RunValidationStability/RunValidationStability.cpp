#include "RunValidationStability.h"

RunValidationStability::RunValidationStability():Tool(){}


bool RunValidationStability::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  m_variables.Get("verbosity",verbosity);
  m_variables.Get("InputDirectory",input_directory);
  m_variables.Get("OutputFile",output_file);
  m_variables.Get("CurrentRun",current_run);
  m_variables.Get("NumberOfRuns",number_runs);
  m_variables.Get("ReferenceRun",reference_run);  
  m_variables.Get("WeeklyMode",WeeklyMode);
  m_variables.Get("WeekStartDate",WeekStartDateFile);
  m_variables.Get("SaveImages",SaveImages);
  m_variables.Get("OutputPath",OutputPath);

  ifstream file_week(WeekStartDateFile);
  std::string line;
  while (std::getline(file_week,line)){
    WeekStartDate = line;
  }
  file_week.close();

  uint64_t timestamp_week;
  uint64_t timestamp_week_end;
  ConvertDateToMSec(WeekStartDate,timestamp_week);
  timestamp_week_end = timestamp_week + 7*24*60*60*1000;
  double MSEC_to_SEC = 1000.;
  double SEC_to_MIN = 60.;
  double MIN_to_HOUR = 60.;
  double HOUR_to_DAY = 24.;

  std::string epoch_start = "1970/1/1";
  boost::posix_time::ptime Epoch(boost::gregorian::from_string(epoch_start));

  boost::posix_time::ptime ptime_week;
  ptime_week = Epoch + boost::posix_time::time_duration(int(timestamp_week/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_week/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_week/MSEC_to_SEC)%60,timestamp_week%1000); 
  boost::posix_time::ptime ptime_week_end;
  ptime_week_end = Epoch + boost::posix_time::time_duration(int(timestamp_week_end/MSEC_to_SEC/SEC_to_MIN/MIN_to_HOUR),int(timestamp_week_end/MSEC_to_SEC/SEC_to_MIN)%60,int(timestamp_week_end/MSEC_to_SEC)%60,timestamp_week_end%1000); 

  struct tm converted_week_tm = boost::posix_time::to_tm(ptime_week);
  struct tm converted_week_end_tm = boost::posix_time::to_tm(ptime_week_end);
  std::stringstream ss_date_week, ss_date_week_end, ss_week, ss_week_end;
  ss_date_week << converted_week_tm.tm_year+1900 << "_" << converted_week_tm.tm_mon+1 << "-" <<converted_week_tm.tm_mday;
  ss_date_week_end << converted_week_end_tm.tm_year+1900 << "_" << converted_week_end_tm.tm_mon+1 << "-" <<converted_week_end_tm.tm_mday;
  ss_week << converted_week_tm.tm_year+1900 << "-" << converted_week_tm.tm_mon+1 << "-" <<converted_week_tm.tm_mday;
  ss_week_end << converted_week_end_tm.tm_year+1900 << "-" << converted_week_end_tm.tm_mon+1 << "-" <<converted_week_end_tm.tm_mday;
  str_week_start = ss_date_week.str();
  str_week_end = ss_date_week_end.str();
  string_week_start = ss_week.str();
  string_week_end = ss_week_end.str();

  int run_start, run_end;
  ConvertWeekToRuns(timestamp_week,run_start, run_end);
  if (run_start == -1) Log("RunValidationStability: Did not find run in the specified week of "+WeekStartDate+"!",0,verbosity);

  if (WeeklyMode){
    current_run = run_end;
    number_runs = run_end - run_start;
  }

  return true;
}


bool RunValidationStability::Execute(){

  std::vector<int> run_numbers;
  std::vector<double> rate, rate_pmt, rate_pmtmrd, rate_pmtmrdfmv, rate_pmtmrdnofmv, rate_fmv, rate_mrd;
  std::vector<double> frac, frac_pmt, frac_pmtmrd, frac_pmtmrdfmv, frac_pmtmrdnofmv, frac_fmv, frac_mrd;
  std::vector<double> align_pmt_mrd, align_pmt_fmv, align_mrd_fmv;
  std::vector<double> frac_prompt, frac_delayed, frac_beam, frac_led, frac_cosmic;
  std::vector<double> avg_mult, avg_mult_coinc, avg_mult_coinc_nofmv_cb;
  std::vector<double> rate_beam, rate_cosmic, rate_led, rate_extended;
  std::vector<double> avg_charge_prompt, avg_charge_del, avg_charge_per_pmt;
  std::vector<double> pot;

  TH1D *PMT_t_clusters_combined = new TH1D("PMT_t_clusters_combined","PMT cluster times (combined)",250,0,2000);
  TH1D *PMT_t_clusters_2pe_combined = new TH1D("PMT_t_clusters_2pe_combined","PMT cluster times (#geq 2.p.e, combined)",250,0,2000);
  TH1D *PMT_t_clusters_cosmics_combined = new TH1D("PMT_t_clusters_cosmics_combined","PMT cluster times Cosmics (combined)",250,0,2000);
  TH1D *PMT_t_clusters_2pe_cosmics_combined = new TH1D("PMT_t_clusters_2pe_cosmics_combined","PMT cluster times Cosmics (#geq 2.p.e, combined)",250,0,2000);
  TH1D *PMT_t_clusters_led_combined = new TH1D("PMT_t_clusters_led_combined","PMT cluster times LED (combined)",250,0,2000);
  TH1D *PMT_t_clusters_2pe_led_combined = new TH1D("PMT_t_clusters_2pe_led_combined","PMT cluster times LED (#geq 2.p.e, combined)",250,0,2000);
  TH1D *MRD_t_clusters_combined = new TH1D("MRD_t_clusters_combined","MRD cluster times (combined)",250,0,4000);
  TH1D *MRD_t_clusters_beam_combined = new TH1D("MRD_t_clusters_beam_combined","MRD cluster times Beam (combined)",250,0,4000);
  TH1D *MRD_t_clusters_cosmic_combined = new TH1D("MRD_t_clusters_cosmic_combined","MRD cluster times Cosmic (combined)",250,0,4000);
  TH1D *FMV_t_clusters_combined = new TH1D("FMV_t_clusters_combined","MRD cluster times (combined)",250,0,4000);
  TH1D *FMV_t_clusters_beam_combined = new TH1D("FMV_t_clusters_beam_combined","MRD cluster times Beam (combined)",250,0,4000);
  TH1D *FMV_t_clusters_cosmic_combined = new TH1D("FMV_t_clusters_cosmic_combined","MRD cluster times Cosmic (combined)",250,0,4000);
  TH1D *PMT_DelayedMult_combined = new TH1D("PMT_DelayedMult_combined","Neutron multiplicity (combined)",20,0,20);
  TH1D *PMT_DelayedMult_Coinc_combined = new TH1D("PMT_DelayedMult_Coinc_combined","Neutron multiplicity (PMT/MRD Coincidence, combined)",20,0,20);
  TH1D *PMT_DelayedMult_Coinc_NoFMV_CB_combined = new TH1D("PMT_DelayedMult_Coinc_NoFMV_CB_combined","Neutron multiplicity (PMT/MRD Coincidence, No FMV, CB < 0.4, combined)",20,0,20);
  TH1D *Triggerwords_combined = new TH1D("Triggerwords_combined","Triggerwords (combined)",60,0,60);  
  TH1D *ADCWaveform_Samples_combined = new TH1D("ADCWaveformSamples_combined","ADC Waveform Samples",5000,0,50000);
  TH1D *PMT_DelayedTime_combined = new TH1D("PMT_DelayedTime_combined","Neutron candidate time distribution (combined)",100,12000,67000);
  TH1D *PMT_DelayedTime_CB_combined = new TH1D("PMT_DelayedTime_CB_combined","Neutron candidate time distribution (CB < 0.4,combined)",100,12000,67000);
  TH1D *PMT_DelayedTime_Coinc_combined = new TH1D("PMT_DelayedTime_Coinc_combined","Neutron candidate time distribution (PMT/MRD Coinc, combined)",100,12000,67000);
  TH1D *PMT_DelayedTime_Coinc_NoFMV_combined = new TH1D("PMT_DelayedTime_Coinc_NoFMV_combined","Neutron candidate time distribution (PMT/MRD Coinc, No FMV, combined)",100,12000,67000);
  TH1D *PMT_DelayedTime_Coinc_NoFMV_CB_combined = new TH1D("PMT_DelayedTime_Coinc_NoFMV_CB_combined","Neutron candidate time distribution (PMT/MRD Coinc, No FMV, CB < 0.4, combined)",100,12000,67000);
  TH1D *PMT_prompt_charge_combined = new TH1D("PMT_prompt_charge_combined","Prompt charge distribution",200,0,5000);
  TH1D *PMT_delayed_charge_zoom_combined = new TH1D("PMT_delayed_charge_zoom_combined","Delayed charge distribution",200,0,200);
  TH1D *PMT_chargeperpmt_combined = new TH1D("PMT_chargeperpmt_combined","Charge per PMT distribution",50,0,30);
  TH2D *MRD_PMT_t_100pe_Beam_combined = new TH2D("MRD_PMT_t_100pe_Beam_combined","MRD vs. PMT Cluster Tiems Beam (>100 p.e.)",50,0,4000,50,0,2000);
  TH2D *FMV_PMT_t_100pe_Beam_combined = new TH2D("FMV_PMT_t_100pe_Beam_combined","FMV vs. PMT Cluster Times Beam (>100 p.e.)",50,0,4000,50,0,2000);
  TH1D *MRD_PMT_Deltat_100pe_Beam_combined = new TH1D("MRD_PMT_Deltat_100pe_Beam_combined","Difference MRD-PMT cluster times (>100 p.e.)",500,-2000,4000);
  TH1D *FMV_PMT_Deltat_100pe_Beam_combined = new TH1D("FMV_PMT_Deltat_100pe_Beam_combined","Difference FMV-PMT cluster times (>100 p.e.)",500,-2000,4000);

  TH1D *PMT_t_clusters_2pe_reference;
  TH1D *MRD_t_clusters_beam_reference;
  TH1D *FMV_t_clusters_beam_reference;
  TH1D *ADCWaveform_Samples_reference;
  TH1D *Triggerwords_reference;
  TH1D *FMV_PMT_Deltat_100pe_Beam_reference;
  TH1D *MRD_PMT_Deltat_100pe_Beam_reference;
  TH1D *PMT_prompt_charge_reference;
  TH1D *PMT_delayed_charge_zoom_reference;
  TH1D *PMT_chargeperpmt_reference;

  double rate_trigger_beam;
  double rate_trigger_beam_extended;
  double rate_trigger_cosmic;
  double rate_trigger_led;
  double pot_total;

  std::stringstream runvalidation_name_reference;
  runvalidation_name_reference << input_directory << "/RunValidation_R" << reference_run << "S0T3.root";
  std::cout <<"Get reference histograms from "<<runvalidation_name_reference.str().c_str()<<std::endl;
  if (!gSystem->AccessPathName(runvalidation_name_reference.str().c_str())){
    std::cout <<"File exists"<<std::endl;
    TFile *fref = new TFile(runvalidation_name_reference.str().c_str(),"READ");
    PMT_t_clusters_2pe_reference = (TH1D*) fref->Get("PMT_t_clusters_2pe");
    PMT_t_clusters_2pe_reference->SetName("PMT_t_clusters_2pe_reference");
    MRD_t_clusters_beam_reference = (TH1D*) fref->Get("MRD_t_clusters_beam");
    MRD_t_clusters_beam_reference->SetName("MRD_t_clusters_beam_reference");
    FMV_t_clusters_beam_reference = (TH1D*) fref->Get("FMV_t_clusters_beam");
    FMV_t_clusters_beam_reference->SetName("FMV_t_clusters_beam_reference");
    ADCWaveform_Samples_reference = (TH1D*) fref->Get("ADCWaveform_Samples");
    ADCWaveform_Samples_reference->SetName("ADCWaveform_Samples_reference");
    Triggerwords_reference = (TH1D*) fref->Get("Triggerwords");
    Triggerwords_reference->SetName("Triggerwords_reference");
    FMV_PMT_Deltat_100pe_Beam_reference = (TH1D*) fref->Get("FMV_PMT_Deltat_Beam_100pe");
    FMV_PMT_Deltat_100pe_Beam_reference->SetName("FMV_PMT_Deltat_100pe_Beam_reference");
    MRD_PMT_Deltat_100pe_Beam_reference = (TH1D*) fref->Get("MRD_PMT_Deltat_100pe_Beam");
    MRD_PMT_Deltat_100pe_Beam_reference->SetName("MRD_PMT_Deltat_100pe_Beam_reference");
    PMT_prompt_charge_reference = (TH1D*) fref->Get("PMT_prompt_charge");
    PMT_prompt_charge_reference->SetName("PMT_prompt_charge_reference");
    PMT_delayed_charge_zoom_reference = (TH1D*) fref->Get("PMT_delayed_charge_zoom");
    PMT_delayed_charge_zoom_reference->SetName("PMT_delayed_charge_zoom_reference");
    PMT_chargeperpmt_reference = (TH1D*) fref->Get("PMT_chargeperpmt");
    PMT_chargeperpmt_reference->SetName("PMT_chargeperpmt_reference");
    //fref->Close();
    //delete fref;
  }
  std::cout <<"Done with histogram loop"<<std::endl;
  PMT_t_clusters_2pe_reference->SetLineColor(8);
  MRD_t_clusters_beam_reference->SetLineColor(8);
  FMV_t_clusters_beam_reference->SetLineColor(8);
  ADCWaveform_Samples_reference->SetLineColor(8);
  Triggerwords_reference->SetLineColor(8);
  FMV_PMT_Deltat_100pe_Beam_reference->SetLineColor(8);
  //MRD_PMT_Deltat_100pe_Beam_reference->SetLineColor(8);
  PMT_prompt_charge_reference->SetLineColor(8);
  PMT_delayed_charge_zoom_reference->SetLineColor(8);
  PMT_chargeperpmt_reference->SetLineColor(8);

  PMT_t_clusters_2pe_reference->SetLineWidth(2);
  MRD_t_clusters_beam_reference->SetLineWidth(2);
  FMV_t_clusters_beam_reference->SetLineWidth(2);
  ADCWaveform_Samples_reference->SetLineWidth(2);
  Triggerwords_reference->SetLineWidth(2);
  FMV_PMT_Deltat_100pe_Beam_reference->SetLineWidth(2);
  //MRD_PMT_Deltat_100pe_Beam_reference->SetLineWidth(2);
  PMT_prompt_charge_reference->SetLineWidth(2);
  PMT_delayed_charge_zoom_reference->SetLineWidth(2);
  PMT_chargeperpmt_reference->SetLineWidth(2);

  for (int i_run = current_run - number_runs; i_run <= current_run; i_run++){
    std::stringstream runvalidation_name;
    runvalidation_name << input_directory << "/RunValidation_R" << i_run << "S0T3.root";
    if (!gSystem->AccessPathName(runvalidation_name.str().c_str())){
      TFile *fval = new TFile(runvalidation_name.str().c_str(),"READ");
      TH1D *PMT_t_clusters = (TH1D*) fval->Get("PMT_t_clusters");
      PMT_t_clusters_combined->Add(PMT_t_clusters);
      TH1D *PMT_t_clusters_2pe = (TH1D*) fval->Get("PMT_t_clusters_2pe");
      PMT_t_clusters_2pe_combined->Add(PMT_t_clusters_2pe);
      TH1D *PMT_t_clusters_cosmics = (TH1D*) fval->Get("PMT_t_clusters_cosmics");
      PMT_t_clusters_cosmics_combined->Add(PMT_t_clusters_cosmics);
      TH1D *PMT_t_clusters_2pe_cosmics = (TH1D*) fval->Get("PMT_t_clusters_2pe_cosmics");
      PMT_t_clusters_2pe_cosmics_combined->Add(PMT_t_clusters_2pe_cosmics);
      TH1D *PMT_t_clusters_led = (TH1D*) fval->Get("PMT_t_clusters_led");
      PMT_t_clusters_led_combined->Add(PMT_t_clusters_led);
      TH1D *PMT_t_clusters_2pe_led = (TH1D*) fval->Get("PMT_t_clusters_2pe_led");
      PMT_t_clusters_2pe_led_combined->Add(PMT_t_clusters_2pe_led);
      TH1D *MRD_t_clusters = (TH1D*) fval->Get("MRD_t_clusters");
      MRD_t_clusters_combined->Add(MRD_t_clusters);
      TH1D *MRD_t_clusters_beam = (TH1D*) fval->Get("MRD_t_clusters_beam");
      MRD_t_clusters_beam_combined->Add(MRD_t_clusters_beam);
      TH1D *MRD_t_clusters_cosmic = (TH1D*) fval->Get("MRD_t_clusters_cosmic");
      MRD_t_clusters_cosmic_combined->Add(MRD_t_clusters_cosmic);
      TH1D *FMV_t_clusters = (TH1D*) fval->Get("FMV_t_clusters");
      FMV_t_clusters_combined->Add(FMV_t_clusters);
      TH1D *FMV_t_clusters_beam = (TH1D*) fval->Get("FMV_t_clusters_beam");
      FMV_t_clusters_beam_combined->Add(FMV_t_clusters_beam);
      TH1D *FMV_t_clusters_cosmic = (TH1D*) fval->Get("FMV_t_clusters_cosmic");
      FMV_t_clusters_cosmic_combined->Add(FMV_t_clusters_cosmic);
      TH1D *MRD_PMT_Deltat_100pe = (TH1D*) fval->Get("MRD_PMT_Deltat_100pe");
      TH1D *FMV_PMT_Deltat_100pe = (TH1D*) fval->Get("FMV_PMT_Deltat_100pe");
      TH1D* MRD_FMV_Deltat = (TH1D*) fval->Get("MRD_FMV_Deltat");
      TH1D *PMT_DelayedMult = (TH1D*) fval->Get("PMT_DelayedMult");
      PMT_DelayedMult_combined->Add(PMT_DelayedMult);
      TH1D *PMT_DelayedMult_Coinc = (TH1D*) fval->Get("PMT_DelayedMult_Coinc");
      PMT_DelayedMult_Coinc_combined->Add(PMT_DelayedMult_Coinc);
      TH1D *PMT_DelayedMult_Coinc_NoFMV_CB = (TH1D*) fval->Get("PMT_DelayedMult_Coinc_NoFMV_CB");
      PMT_DelayedMult_Coinc_NoFMV_CB_combined->Add(PMT_DelayedMult_Coinc_NoFMV_CB);
      TH1D *ADCWaveform_Samples = (TH1D*) fval->Get("ADCWaveform_Samples");
      ADCWaveform_Samples_combined->Add(ADCWaveform_Samples);
      TH1D *Triggerwords = (TH1D*) fval->Get("Triggerwords");
      Triggerwords_combined->Add(Triggerwords);
      TH1D *ANNIE_rates = (TH1D*) fval->Get("ANNIE_rates");
      TH1D *ANNIE_fractions = (TH1D*) fval->Get("ANNIE_fractions");
      TH1D *PMT_DelayedTime = (TH1D*) fval->Get("PMT_DelayedTime");
      PMT_DelayedTime_combined->Add(PMT_DelayedTime);
      TH1D *PMT_DelayedTime_CB = (TH1D*) fval->Get("PMT_DelayedTime_CB");
      PMT_DelayedTime_CB_combined->Add(PMT_DelayedTime_CB);
      TH1D *PMT_DelayedTime_Coinc = (TH1D*) fval->Get("PMT_DelayedTime_Coinc");
      PMT_DelayedTime_Coinc_combined->Add(PMT_DelayedTime_Coinc);
      TH1D *PMT_DelayedTime_Coinc_NoFMV = (TH1D*) fval->Get("PMT_DelayedTime_Coinc_NoFMV");
      PMT_DelayedTime_Coinc_NoFMV_combined->Add(PMT_DelayedTime_Coinc_NoFMV);
      TH1D *PMT_DelayedTime_Coinc_NoFMV_CB = (TH1D*) fval->Get("PMT_DelayedTime_Coinc_NoFMV_CB");
      PMT_DelayedTime_Coinc_NoFMV_CB_combined->Add(PMT_DelayedTime_Coinc_NoFMV_CB);
      TH1D *PMT_prompt_charge = (TH1D*) fval->Get("PMT_prompt_charge");
      PMT_prompt_charge_combined->Add(PMT_prompt_charge);
      TH1D *PMT_delayed_charge_zoom = (TH1D*) fval->Get("PMT_delayed_charge_zoom");
      PMT_delayed_charge_zoom = (TH1D*) fval->Get("PMT_delayed_charge_zoom");
      PMT_delayed_charge_zoom_combined->Add(PMT_delayed_charge_zoom);
      TH1D *PMT_chargeperpmt = (TH1D*) fval->Get("PMT_chargeperpmt");
      PMT_chargeperpmt_combined->Add(PMT_chargeperpmt);
      TH2D *MRD_PMT_t_100pe_Beam = (TH2D*) fval->Get("MRD_PMT_t_100pe_Beam");
      MRD_PMT_t_100pe_Beam_combined->Add(MRD_PMT_t_100pe_Beam);
      TH2D *FMV_PMT_t_100pe_Beam = (TH2D*) fval->Get("FMV_PMT_t_100pe_Beam");
      FMV_PMT_t_100pe_Beam_combined->Add(FMV_PMT_t_100pe_Beam); 
      TH1D *MRD_PMT_Deltat_100pe_Beam = (TH1D*) fval->Get("MRD_PMT_Deltat_100pe_Beam");
      MRD_PMT_Deltat_100pe_Beam_combined->Add(MRD_PMT_Deltat_100pe_Beam);
      TH1D *FMV_PMT_Deltat_100pe_Beam = (TH1D*) fval->Get("FMV_PMT_Deltat_Beam_100pe");
      FMV_PMT_Deltat_100pe_Beam_combined->Add(FMV_PMT_Deltat_100pe_Beam); 

      double avg_q_prompt = PMT_prompt_charge->GetMean();
      double avg_q_del = PMT_delayed_charge_zoom->GetMean();
      double avg_qperpmt = PMT_chargeperpmt->GetMean();
      avg_charge_prompt.push_back(avg_q_prompt);
      avg_charge_del.push_back(avg_q_del);
      avg_charge_per_pmt.push_back(avg_qperpmt);

      TTree *t_stats = (TTree*) fval->Get("tree_stats");
      t_stats->SetBranchAddress("rate_trigger_beam",&rate_trigger_beam);
      t_stats->SetBranchAddress("rate_trigger_beam_extended",&rate_trigger_beam_extended);
      t_stats->SetBranchAddress("rate_trigger_cosmic",&rate_trigger_cosmic);
      t_stats->SetBranchAddress("rate_trigger_led",&rate_trigger_led);
      t_stats->SetBranchAddress("pot_total",&pot_total);
      t_stats->GetEntry(0);
      rate_beam.push_back(rate_trigger_beam);
      rate_extended.push_back(rate_trigger_beam_extended);
      rate_cosmic.push_back(rate_trigger_cosmic);
      rate_led.push_back(rate_trigger_led);
      pot.push_back(pot_total);

      std::cout <<"Getting ptr from MRD_PMT_Deltat"<<std::endl;
      TF1 *fgaus = new TF1("fgaus","gaus+pol0(3)",500,900);
      fgaus->SetParameter(0,MRD_PMT_Deltat_100pe_Beam->GetMaximum());
      fgaus->SetParLimits(0,0,2*MRD_PMT_Deltat_100pe_Beam->GetMaximum());
      fgaus->SetParameter(1,750);
      fgaus->SetParLimits(1,500,1000);
      fgaus->SetParameter(2,20);
      fgaus->SetParLimits(2,0,100);
      TFitResultPtr ptr = MRD_PMT_Deltat_100pe_Beam->Fit("fgaus","S","",500,900);
      Int_t ptrInt = ptr;
      if (ptrInt==0){
      align_pmt_mrd.push_back(ptr->Parameter(1));
      } else {
        align_pmt_mrd.push_back(MRD_PMT_Deltat_100pe_Beam->GetMean());
      }
      TFitResultPtr ptr2 = FMV_PMT_Deltat_100pe->Fit("gaus","S","",750,850);
      Int_t ptr2Int = ptr2;
      if (ptr2Int==0){
      align_pmt_fmv.push_back(ptr2->Parameter(1));
      } else align_pmt_fmv.push_back(0);
      TFitResultPtr ptr3 = MRD_FMV_Deltat->Fit("gaus","S","",0,100);
      Int_t ptr3Int = ptr3;
      if (ptr3Int==0){
      align_mrd_fmv.push_back(ptr3->Parameter(1));
      } else align_mrd_fmv.push_back(0);

      double temp_avg_mult=0;
      double temp_avg_mult_coinc=0;
      double temp_avg_mult_coinc_nofmv_cb=0;
      for (int i_bin=0; i_bin < PMT_DelayedMult->GetXaxis()->GetNbins(); i_bin++){
        temp_avg_mult+=(i_bin*PMT_DelayedMult->GetBinContent(i_bin+1));
        temp_avg_mult_coinc+=(i_bin*PMT_DelayedMult_Coinc->GetBinContent(i_bin+1));
        temp_avg_mult_coinc_nofmv_cb+=(i_bin*PMT_DelayedMult_Coinc_NoFMV_CB->GetBinContent(i_bin+1));
      }
      int nentries = PMT_DelayedMult->GetEntries();
      int nentries_coinc = PMT_DelayedMult_Coinc->GetEntries();
      int nentries_coinc_nofmv_cb = PMT_DelayedMult_Coinc_NoFMV_CB->GetEntries();
      if (nentries > 0) temp_avg_mult /= nentries;
      if (nentries_coinc > 0) temp_avg_mult_coinc /= nentries_coinc;
      if (nentries_coinc_nofmv_cb > 0) temp_avg_mult_coinc_nofmv_cb /= nentries_coinc_nofmv_cb;
      avg_mult.push_back(temp_avg_mult);
      avg_mult_coinc.push_back(temp_avg_mult_coinc);
      avg_mult_coinc_nofmv_cb.push_back(temp_avg_mult_coinc_nofmv_cb); 

      double n_beam = Triggerwords->GetBinContent(6);
      double n_led = Triggerwords->GetBinContent(32);
      double n_cosmic = Triggerwords->GetBinContent(37);
      double n_total = n_beam+n_led+n_cosmic;
      if (n_total>0){
        n_beam/=n_total;
        n_led/=n_total;
        n_cosmic/=n_total;
      }
      frac_beam.push_back(n_beam);
      frac_led.push_back(n_led);
      frac_cosmic.push_back(n_cosmic);
 
      double n_prompt = ADCWaveform_Samples->GetBinContent(100);
      double n_delayed = ADCWaveform_Samples->GetBinContent(3500);
      double n_prompt_delayed = n_prompt+n_delayed;
      if (n_prompt_delayed){
        n_prompt /= n_prompt_delayed;
        n_delayed /= n_prompt_delayed;
      }
      frac_prompt.push_back(n_prompt);
      frac_delayed.push_back(n_delayed);

      rate.push_back(ANNIE_rates->GetBinContent(1));
      rate_pmt.push_back(ANNIE_rates->GetBinContent(2));
      rate_pmtmrd.push_back(ANNIE_rates->GetBinContent(5));
      rate_pmtmrdfmv.push_back(ANNIE_rates->GetBinContent(10));
      rate_pmtmrdnofmv.push_back(ANNIE_rates->GetBinContent(6));
      rate_fmv.push_back(ANNIE_rates->GetBinContent(7));
      rate_mrd.push_back(ANNIE_rates->GetBinContent(8));
      
      frac.push_back(ANNIE_fractions->GetBinContent(1));
      frac_pmt.push_back(ANNIE_fractions->GetBinContent(2));
      frac_pmtmrd.push_back(ANNIE_fractions->GetBinContent(5));
      frac_pmtmrdfmv.push_back(ANNIE_fractions->GetBinContent(10));
      frac_pmtmrdnofmv.push_back(ANNIE_fractions->GetBinContent(6));
      frac_fmv.push_back(ANNIE_fractions->GetBinContent(7));
      frac_mrd.push_back(ANNIE_fractions->GetBinContent(8));
      run_numbers.push_back(i_run);
      fval->Close();
      delete fval;
    }
  }

  //Open output file and fill TGraphs
  std::stringstream filename_out;
  if (WeeklyMode) filename_out << "RunValidationStability_"<<str_week_start<<"-"<<str_week_end<<".root";
  else filename_out << "RunValidationStability_R" << current_run - number_runs << "-" << current_run << ".root";
  TFile *file_output = new TFile(filename_out.str().c_str(),"RECREATE");
  
  TGraph *gr_rate = new TGraph();
  TGraph *gr_rate_pmt = new TGraph();
  TGraph *gr_rate_pmtmrd = new TGraph(); 
  TGraph *gr_rate_pmtmrdfmv = new TGraph(); 
  TGraph *gr_rate_pmtmrdnofmv = new TGraph(); 
  TGraph *gr_rate_fmv = new TGraph();
  TGraph *gr_rate_mrd = new TGraph();
  TGraph *gr_frac = new TGraph();
  TGraph *gr_frac_pmt = new TGraph();
  TGraph *gr_frac_pmtmrd = new TGraph();
  TGraph *gr_frac_pmtmrdfmv = new TGraph();
  TGraph *gr_frac_pmtmrdnofmv = new TGraph();
  TGraph *gr_frac_fmv = new TGraph();
  TGraph *gr_frac_mrd = new TGraph();
  TGraph *gr_align_pmt_mrd = new TGraph();
  TGraph *gr_align_pmt_fmv = new TGraph();
  TGraph *gr_align_mrd_fmv = new TGraph();
  TGraph *gr_frac_prompt = new TGraph();
  TGraph *gr_frac_delayed = new TGraph();
  TGraph *gr_frac_beam = new TGraph();
  TGraph *gr_frac_led = new TGraph();
  TGraph *gr_frac_cosmic = new TGraph();
  TGraph *gr_avg_mult = new TGraph();
  TGraph *gr_avg_mult_coinc = new TGraph();
  TGraph *gr_avg_mult_coinc_nofmv_cb = new TGraph();
  TGraph *gr_trigger_beam = new TGraph();
  TGraph *gr_trigger_cosmic = new TGraph();
  TGraph *gr_trigger_led = new TGraph();
  TGraph *gr_trigger_extended = new TGraph();
  TGraph *gr_charge_prompt = new TGraph();
  TGraph *gr_charge_del = new TGraph();
  TGraph *gr_charge_per_pmt = new TGraph();
  TGraph *gr_pot = new TGraph();
  TGraph *gr_pot_scaled = new TGraph();

  gr_rate->SetLineStyle(2);
  gr_rate_pmt->SetLineStyle(2);
  gr_rate_pmtmrd->SetLineStyle(2);
  gr_rate_pmtmrdfmv->SetLineStyle(2);
  gr_rate_pmtmrdnofmv->SetLineStyle(2);
  gr_rate_fmv->SetLineStyle(2);
  gr_rate_mrd->SetLineStyle(2);
  gr_frac->SetLineStyle(2);
  gr_frac_pmt->SetLineStyle(2);
  gr_frac_pmtmrd->SetLineStyle(2);
  gr_frac_pmtmrdfmv->SetLineStyle(2);
  gr_frac_pmtmrdnofmv->SetLineStyle(2);
  gr_frac_fmv->SetLineStyle(2);
  gr_frac_mrd->SetLineStyle(2);
  gr_align_pmt_mrd->SetLineStyle(2);
  gr_align_pmt_fmv->SetLineStyle(2);
  gr_align_mrd_fmv->SetLineStyle(2);
  gr_frac_prompt->SetLineStyle(2);
  gr_frac_delayed->SetLineStyle(2);
  gr_frac_beam->SetLineStyle(2);
  gr_frac_led->SetLineStyle(2);
  gr_frac_cosmic->SetLineStyle(2);
  gr_avg_mult->SetLineStyle(2);
  gr_avg_mult_coinc->SetLineStyle(2);
  gr_avg_mult_coinc_nofmv_cb->SetLineStyle(2);
  gr_trigger_beam->SetLineStyle(2);
  gr_trigger_cosmic->SetLineStyle(2);
  gr_trigger_led->SetLineStyle(2);
  gr_trigger_extended->SetLineStyle(2);
  gr_charge_prompt->SetLineStyle(2);
  gr_charge_del->SetLineStyle(2);
  gr_charge_per_pmt->SetLineStyle(2);  
  gr_pot->SetLineStyle(2);
  gr_pot_scaled->SetLineStyle(2);

  gr_rate->SetLineWidth(2);
  gr_rate_pmt->SetLineWidth(2);
  gr_rate_pmtmrd->SetLineWidth(2);
  gr_rate_pmtmrdfmv->SetLineWidth(2);
  gr_rate_pmtmrdnofmv->SetLineWidth(2);
  gr_rate_fmv->SetLineWidth(2);
  gr_rate_mrd->SetLineWidth(2);
  gr_frac->SetLineWidth(2);
  gr_frac_pmt->SetLineWidth(2);
  gr_frac_pmtmrd->SetLineWidth(2);
  gr_frac_pmtmrdfmv->SetLineWidth(2);
  gr_frac_pmtmrdnofmv->SetLineWidth(2);
  gr_frac_fmv->SetLineWidth(2);
  gr_frac_mrd->SetLineWidth(2);
  gr_align_pmt_mrd->SetLineWidth(2);
  gr_align_pmt_fmv->SetLineWidth(2);
  gr_align_mrd_fmv->SetLineWidth(2);
  gr_frac_prompt->SetLineWidth(2);
  gr_frac_delayed->SetLineWidth(2);
  gr_frac_beam->SetLineWidth(2);
  gr_frac_led->SetLineWidth(2);
  gr_frac_cosmic->SetLineWidth(2);
  gr_avg_mult->SetLineWidth(2);
  gr_avg_mult_coinc->SetLineWidth(2);
  gr_avg_mult_coinc_nofmv_cb->SetLineWidth(2);
  gr_trigger_beam->SetLineWidth(2);
  gr_trigger_cosmic->SetLineWidth(2);
  gr_trigger_led->SetLineWidth(2);
  gr_trigger_extended->SetLineWidth(2);  
  gr_charge_prompt->SetLineWidth(2);
  gr_charge_del->SetLineWidth(2);
  gr_charge_per_pmt->SetLineWidth(2);
  gr_pot->SetLineWidth(2);
  gr_pot_scaled->SetLineWidth(2);

  gr_rate->SetMarkerStyle(21);
  gr_rate_pmt->SetMarkerStyle(21);
  gr_rate_pmtmrd->SetMarkerStyle(21);
  gr_rate_pmtmrdfmv->SetMarkerStyle(21);
  gr_rate_pmtmrdnofmv->SetMarkerStyle(21);
  gr_rate_fmv->SetMarkerStyle(21);
  gr_rate_mrd->SetMarkerStyle(21);
  gr_frac->SetMarkerStyle(21);
  gr_frac_pmt->SetMarkerStyle(21);
  gr_frac_pmtmrd->SetMarkerStyle(21);
  gr_frac_pmtmrdfmv->SetMarkerStyle(21);
  gr_frac_pmtmrdnofmv->SetMarkerStyle(21);
  gr_frac_fmv->SetMarkerStyle(21);
  gr_frac_mrd->SetMarkerStyle(21);
  gr_align_pmt_mrd->SetMarkerStyle(21);
  gr_align_pmt_fmv->SetMarkerStyle(21);
  gr_align_mrd_fmv->SetMarkerStyle(21);
  gr_frac_prompt->SetMarkerStyle(21);
  gr_frac_delayed->SetMarkerStyle(21);
  gr_frac_beam->SetMarkerStyle(21);
  gr_frac_led->SetMarkerStyle(21);
  gr_frac_cosmic->SetMarkerStyle(21);
  gr_avg_mult->SetMarkerStyle(21);
  gr_avg_mult_coinc->SetMarkerStyle(21);
  gr_avg_mult_coinc_nofmv_cb->SetMarkerStyle(21);
  gr_trigger_beam->SetMarkerStyle(21);
  gr_trigger_cosmic->SetMarkerStyle(21);
  gr_trigger_led->SetMarkerStyle(21);
  gr_trigger_extended->SetMarkerStyle(21);  
  gr_charge_prompt->SetMarkerStyle(21);
  gr_charge_del->SetMarkerStyle(21);
  gr_charge_per_pmt->SetMarkerStyle(21);
  gr_pot->SetMarkerStyle(21);
  gr_pot_scaled->SetMarkerStyle(21);

  gr_rate->SetMarkerSize(0.7);
  gr_rate_pmt->SetMarkerSize(0.7);
  gr_rate_pmtmrd->SetMarkerSize(0.7);
  gr_rate_pmtmrdfmv->SetMarkerSize(0.7);
  gr_rate_pmtmrdnofmv->SetMarkerSize(0.7);
  gr_rate_fmv->SetMarkerSize(0.7);
  gr_rate_mrd->SetMarkerSize(0.7);
  gr_frac->SetMarkerSize(0.7);
  gr_frac_pmt->SetMarkerSize(0.7);
  gr_frac_pmtmrd->SetMarkerSize(0.7);
  gr_frac_pmtmrdfmv->SetMarkerSize(0.7);
  gr_frac_pmtmrdnofmv->SetMarkerSize(0.7);
  gr_frac_fmv->SetMarkerSize(0.7);
  gr_frac_mrd->SetMarkerSize(0.7);
  gr_align_pmt_mrd->SetMarkerSize(0.7);
  gr_align_pmt_fmv->SetMarkerSize(0.7);
  gr_align_mrd_fmv->SetMarkerSize(0.7);
  gr_frac_prompt->SetMarkerSize(0.7);
  gr_frac_delayed->SetMarkerSize(0.7);
  gr_frac_beam->SetMarkerSize(0.7);
  gr_frac_led->SetMarkerSize(0.7);
  gr_frac_cosmic->SetMarkerSize(0.7);
  gr_avg_mult->SetMarkerSize(0.7);
  gr_avg_mult_coinc->SetMarkerSize(0.7);
  gr_avg_mult_coinc_nofmv_cb->SetMarkerSize(0.7);
  gr_trigger_beam->SetMarkerSize(0.7);
  gr_trigger_cosmic->SetMarkerSize(0.7);
  gr_trigger_led->SetMarkerSize(0.7);
  gr_trigger_extended->SetMarkerSize(0.7);
  gr_charge_prompt->SetMarkerSize(0.7);
  gr_charge_del->SetMarkerSize(0.7);
  gr_charge_per_pmt->SetMarkerSize(0.7);
  gr_pot->SetMarkerSize(0.7);
  gr_pot_scaled->SetMarkerSize(0.7);

  std::stringstream ss_title_pmt, ss_title_pmt_2pe, ss_title_pmt_cosmics, ss_title_pmt_2pe_cosmics, ss_title_pmt_led, ss_title_pmt_2pe_led, ss_title_mrd, ss_title_mrd_beam, ss_title_mrd_cosmic, ss_title_fmv, ss_title_fmv_beam, ss_title_fmv_cosmic;
  std::stringstream ss_title_mult, ss_title_mult_coinc, ss_title_mult_coinc_nofmv_cb, ss_title_trigword, ss_title_adc, ss_title_delayedt, ss_title_delayedt_coinc, ss_title_delayedt_coinc_nofmv, ss_title_delayedt_coinc_nofmv_cb, ss_title_delayedt_cb, ss_title_mult_coinc_nofmv;
  std::stringstream ss_title_charge_prompt, ss_title_charge_del, ss_title_charge_per_pmt;
  std::stringstream ss_title_mrd_pmt_100pe, ss_title_fmv_pmt_100pe;
  std::stringstream ss_title_mrd_pmt_deltat, ss_title_fmv_pmt_deltat;
  if (!WeeklyMode){
  ss_title_pmt << "PMT Cluster times (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_pmt_2pe << "PMT Cluster times (#geq 2p.e.,Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_pmt_cosmics << "PMT Cluster times Cosmics (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_pmt_2pe_cosmics << "PMT Cluster times Cosmics (#geq 2p.e., Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_pmt_led << "PMT Cluster times LED (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_pmt_2pe_led << "PMT Cluster times LED (#geq 2p.e., Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_mrd << "MRD Cluster times (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_mrd_beam << "MRD Cluster times Beam (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_mrd_cosmic << "MRD Cluster times Cosmic (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_fmv << "FMV Cluster times (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_fmv_beam << "FMV Cluster times Beam (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_fmv_cosmic << "FMV Cluster times Cosmic (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_mult << "Neutron multiplicity (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_mult_coinc << "Neutron multiplicity PMT/MRD Coincidence (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_title_mult_coinc_nofmv << "Neutron multiplicity PMT/MRD Coincidence, CB < 0.4 (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_title_mult_coinc_nofmv_cb << "Neutron multiplicity PMT/MRD Coincidence, No FMV, CB < 0.4 (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_title_trigword << "Triggerwords (Run "<<current_run - number_runs << "- Run "<<current_run <<")";
  ss_title_adc << "ADC Waveform Samples (Run "<<current_run - number_runs <<"- Run "<<current_run <<")";
  ss_title_delayedt << "Neutron candidate time distribution (Run "<<current_run - number_runs << "- Run "<<current_run <<")";
  ss_title_delayedt_cb << "Neutron candidate time distribution CB < 0.4 (Run "<<current_run - number_runs << "- Run "<<current_run <<")";
  ss_title_delayedt_coinc << "Neutron candidate time distribution PMT/MRD Coincidence (Run "<<current_run - number_runs << "- Run "<<current_run <<")";
  ss_title_delayedt_coinc_nofmv << "Neutron candidate time distribution PMT/MRD Coinc, No FMV (Run "<<current_run - number_runs << "- Run "<<current_run <<")";
  ss_title_delayedt_coinc_nofmv_cb << "Neutron candidate time distribution PMT/MRD Coinc, No FMV, CB < 0.4 (Run "<<current_run - number_runs << "- Run "<<current_run <<")";
  ss_title_charge_prompt << "Average prompt charge (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_title_charge_del << "Average delayed charge (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_title_charge_per_pmt << "Average charge per pmt (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_title_mrd_pmt_100pe << "MRD vs. PMT Cluster Times Beam (>100 p.e.) (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_title_fmv_pmt_100pe << "FMV vs. PMT Cluster Times Beam (>100 p.e.) (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_title_mrd_pmt_deltat << "Difference MRD/PMT Cluster Times (>100 p.e.) (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_title_fmv_pmt_deltat << "Difference FMV/PMT Cluster Times (>100 p.e.) (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  } else {
  ss_title_pmt << "PMT Cluster times ( "<<string_week_start<< " - "<<string_week_end<<")";
  ss_title_pmt_2pe << "PMT Cluster times (#geq 2p.e., "<<string_week_start << " - "<<string_week_end<<")";
  ss_title_pmt_cosmics << "PMT Cluster times Cosmics ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_title_pmt_2pe_cosmics << "PMT Cluster times Cosmics (#geq 2p.e., "<<string_week_start << " - "<<string_week_end<<")";
  ss_title_pmt_led << "PMT Cluster times LED ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_title_pmt_2pe_led << "PMT Cluster times LED (#geq 2p.e., "<< string_week_start << " - "<<string_week_end<<")";
  ss_title_mrd << "MRD Cluster times ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_title_mrd_beam << "MRD Cluster times Beam ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_title_mrd_cosmic << "MRD Cluster times Cosmic ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_title_fmv << "FMV Cluster times ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_title_fmv_beam << "FMV Cluster times Beam ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_title_fmv_cosmic << "FMV Cluster times Cosmic ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_title_mult << "Neutron multiplicity ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_title_mult_coinc << "Neutron multiplicity PMT/MRD Coincidence ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_title_mult_coinc_nofmv << "Neutron multiplicity PMT/MRD Coincidence, CB < 0.4 ( "<< string_week_start << " - "<<string_week_end << ")";
  ss_title_mult_coinc_nofmv_cb << "Neutron multiplicity PMT/MRD Coincidence, No FMV, CB < 0.4 ( "<<string_week_start << " - "<<string_week_end << ")";
  ss_title_trigword << "Triggerwords ( "<<string_week_start << " - "<<string_week_end <<")";
  ss_title_adc << "ADC Waveform Samples ( "<<string_week_start <<" - "<<string_week_end <<")";
  ss_title_delayedt << "Neutron candidate time distribution ( "<<string_week_start << " - "<<string_week_end <<")";
  ss_title_delayedt_cb << "Neutron candidate time distribution CB < 0.4 ( "<<string_week_start << " - "<<string_week_end <<")";
  ss_title_delayedt_coinc << "Neutron candidate time distribution PMT/MRD Coincidence ( "<<string_week_start << " - "<<string_week_end <<")";
  ss_title_delayedt_coinc_nofmv << "Neutron candidate time distribution PMT/MRD Coinc, No FMV ( "<< string_week_start << " - "<<string_week_end <<")";
  ss_title_delayedt_coinc_nofmv_cb << "Neutron candidate time distribution PMT/MRD Coinc, No FMV, CB < 0.4 ( "<<string_week_start << " - "<<string_week_end <<")";
  ss_title_charge_prompt << "Average prompt charge ( "<<string_week_start << " - "<<string_week_end << ")";
  ss_title_charge_del << "Average delayed charge ( "<<string_week_start << " - "<<string_week_end << ")";
  ss_title_charge_per_pmt << "Average charge per pmt ( "<<string_week_start << " - "<<string_week_end << ")";
  ss_title_mrd_pmt_100pe << "MRD vs. PMT Cluster Times Beam (>100 p.e.) ( "<<string_week_start << " - "<< string_week_end << ")";
  ss_title_fmv_pmt_100pe << "FMV vs. PMT Cluster Times Beam (>100 p.e.) ( "<<string_week_start << " - "<< string_week_end << ")";
  ss_title_mrd_pmt_deltat << "Difference MRD/PMT Cluster Times (>100 p.e.) ( "<<string_week_start << " - "<<string_week_end << ")";
  ss_title_fmv_pmt_deltat << "Difference FMV/PMT Cluster Times (>100 p.e.) ( "<<string_week_start << " - "<<string_week_end << ")";

  }
  PMT_t_clusters_combined->SetTitle(ss_title_pmt.str().c_str());
  PMT_t_clusters_2pe_combined->SetTitle(ss_title_pmt_2pe.str().c_str());
  PMT_t_clusters_cosmics_combined->SetTitle(ss_title_pmt_cosmics.str().c_str());
  PMT_t_clusters_2pe_cosmics_combined->SetTitle(ss_title_pmt_2pe_cosmics.str().c_str());
  PMT_t_clusters_led_combined->SetTitle(ss_title_pmt_led.str().c_str());
  PMT_t_clusters_2pe_led_combined->SetTitle(ss_title_pmt_2pe_led.str().c_str());
  MRD_t_clusters_combined->SetTitle(ss_title_mrd.str().c_str());
  MRD_t_clusters_beam_combined->SetTitle(ss_title_mrd_beam.str().c_str());
  MRD_t_clusters_cosmic_combined->SetTitle(ss_title_mrd_cosmic.str().c_str());
  FMV_t_clusters_combined->SetTitle(ss_title_fmv.str().c_str());
  FMV_t_clusters_beam_combined->SetTitle(ss_title_fmv_beam.str().c_str());
  FMV_t_clusters_cosmic_combined->SetTitle(ss_title_fmv_cosmic.str().c_str());
  PMT_DelayedMult_combined->SetTitle(ss_title_mult.str().c_str());
  PMT_DelayedMult_Coinc_combined->SetTitle(ss_title_mult_coinc.str().c_str());
  PMT_DelayedMult_Coinc_NoFMV_CB_combined->SetTitle(ss_title_mult_coinc_nofmv_cb.str().c_str());
  Triggerwords_combined->SetTitle(ss_title_trigword.str().c_str());
  PMT_DelayedTime_combined->SetTitle(ss_title_delayedt.str().c_str());
  PMT_DelayedTime_CB_combined->SetTitle(ss_title_delayedt_cb.str().c_str());
  PMT_DelayedTime_Coinc_combined->SetTitle(ss_title_delayedt_coinc.str().c_str());
  PMT_DelayedTime_Coinc_NoFMV_combined->SetTitle(ss_title_delayedt_coinc_nofmv.str().c_str());
  PMT_DelayedTime_Coinc_NoFMV_CB_combined->SetTitle(ss_title_delayedt_coinc_nofmv_cb.str().c_str());
  PMT_prompt_charge_combined->SetTitle(ss_title_charge_prompt.str().c_str());
  PMT_delayed_charge_zoom_combined->SetTitle(ss_title_charge_del.str().c_str());
  PMT_chargeperpmt_combined->SetTitle(ss_title_charge_per_pmt.str().c_str());
  MRD_PMT_t_100pe_Beam_combined->SetTitle(ss_title_mrd_pmt_100pe.str().c_str());
  FMV_PMT_t_100pe_Beam_combined->SetTitle(ss_title_fmv_pmt_100pe.str().c_str());
  MRD_PMT_Deltat_100pe_Beam_combined->SetTitle(ss_title_mrd_pmt_deltat.str().c_str());
  FMV_PMT_Deltat_100pe_Beam_combined->SetTitle(ss_title_fmv_pmt_deltat.str().c_str());

  PMT_t_clusters_combined->GetXaxis()->SetTitle("t_{PMT} [ns]");
  PMT_t_clusters_2pe_combined->GetXaxis()->SetTitle("t_{PMT} [ns]");
  PMT_t_clusters_cosmics_combined->GetXaxis()->SetTitle("t_{PMT} [ns]");
  PMT_t_clusters_2pe_cosmics_combined->GetXaxis()->SetTitle("t_{PMT} [ns]");
  PMT_t_clusters_led_combined->GetXaxis()->SetTitle("t_{PMT} [ns]");
  PMT_t_clusters_2pe_led_combined->GetXaxis()->SetTitle("t_{PMT} [ns]");
  MRD_t_clusters_combined->GetXaxis()->SetTitle("t_{MRD} [ns]");
  MRD_t_clusters_beam_combined->GetXaxis()->SetTitle("t_{MRD} [ns]");
  MRD_t_clusters_cosmic_combined->GetXaxis()->SetTitle("t_{MRD} [ns]");
  FMV_t_clusters_combined->GetXaxis()->SetTitle("t_{FMV} [ns]");
  FMV_t_clusters_beam_combined->GetXaxis()->SetTitle("t_{FMV} [ns]");
  FMV_t_clusters_cosmic_combined->GetXaxis()->SetTitle("t_{FMV} [ns]");
  PMT_DelayedMult_combined->GetXaxis()->SetTitle("N_{neutron}");
  PMT_DelayedMult_Coinc_combined->GetXaxis()->SetTitle("N_{neutron}");
  PMT_DelayedMult_Coinc_NoFMV_CB_combined->GetXaxis()->SetTitle("N_{neutron}");
  Triggerwords_combined->GetXaxis()->SetTitle("Triggerword");
  PMT_DelayedTime_combined->GetXaxis()->SetTitle("t_{cluster} [ns]"); 
  PMT_DelayedTime_CB_combined->GetXaxis()->SetTitle("t_{cluster} [ns]"); 
  PMT_DelayedTime_Coinc_combined->GetXaxis()->SetTitle("t_{cluster} [ns]"); 
  PMT_DelayedTime_Coinc_NoFMV_combined->GetXaxis()->SetTitle("t_{cluster} [ns]"); 
  PMT_DelayedTime_Coinc_NoFMV_CB_combined->GetXaxis()->SetTitle("t_{cluster} [ns]"); 
  PMT_prompt_charge_combined->GetXaxis()->SetTitle("q [p.e.]");
  PMT_delayed_charge_zoom_combined->GetXaxis()->SetTitle("q [p.e.]");
  PMT_chargeperpmt_combined->GetXaxis()->SetTitle("q/n_{PMT} [p.e.]");
  MRD_PMT_t_100pe_Beam_combined->GetXaxis()->SetTitle("t_{MRD} [ns]");
  FMV_PMT_t_100pe_Beam_combined->GetXaxis()->SetTitle("t_{FMV} [ns]");
  MRD_PMT_Deltat_100pe_Beam_combined->GetXaxis()->SetTitle("t_{MRD} - t_{PMT} [ns]");
  FMV_PMT_Deltat_100pe_Beam_combined->GetXaxis()->SetTitle("t_{FMV} - t_{PMT} [ns]");

  PMT_t_clusters_combined->GetYaxis()->SetTitle("#");
  PMT_t_clusters_2pe_combined->GetYaxis()->SetTitle("#");
  PMT_t_clusters_cosmics_combined->GetYaxis()->SetTitle("#");
  PMT_t_clusters_2pe_cosmics_combined->GetYaxis()->SetTitle("#");
  PMT_t_clusters_led_combined->GetYaxis()->SetTitle("#");
  PMT_t_clusters_2pe_led_combined->GetYaxis()->SetTitle("#");
  PMT_t_clusters_2pe_cosmics_combined->GetYaxis()->SetTitle("#");
  MRD_t_clusters_combined->GetYaxis()->SetTitle("#");
  PMT_DelayedMult_combined->GetYaxis()->SetTitle("#");
  PMT_DelayedMult_Coinc_combined->GetYaxis()->SetTitle("#");
  PMT_DelayedMult_Coinc_NoFMV_CB_combined->GetYaxis()->SetTitle("#");
  Triggerwords_combined->GetYaxis()->SetTitle("#");
  PMT_DelayedTime_combined->GetYaxis()->SetTitle("#");
  PMT_DelayedTime_CB_combined->GetYaxis()->SetTitle("#");
  PMT_DelayedTime_Coinc_combined->GetYaxis()->SetTitle("#");
  PMT_DelayedTime_Coinc_NoFMV_combined->GetYaxis()->SetTitle("#");
  PMT_DelayedTime_Coinc_NoFMV_CB_combined->GetYaxis()->SetTitle("#");
  PMT_prompt_charge_combined->GetYaxis()->SetTitle("#");
  PMT_delayed_charge_zoom_combined->GetYaxis()->SetTitle("#");
  PMT_chargeperpmt_combined->GetYaxis()->SetTitle("#");
  MRD_PMT_t_100pe_Beam_combined->GetYaxis()->SetTitle("t_{PMT} [ns]");
  FMV_PMT_t_100pe_Beam_combined->GetYaxis()->SetTitle("t_{PMT} [ns]");
  MRD_PMT_Deltat_100pe_Beam_combined->GetYaxis()->SetTitle("#");
  FMV_PMT_Deltat_100pe_Beam_combined->GetYaxis()->SetTitle("#");

  for (int i_run=0; i_run < (int) run_numbers.size(); i_run++){

    gr_rate->SetPoint(i_run,run_numbers.at(i_run),rate.at(i_run));
    gr_rate_pmt->SetPoint(i_run,run_numbers.at(i_run),rate_pmt.at(i_run));
    gr_rate_pmtmrd->SetPoint(i_run,run_numbers.at(i_run),rate_pmtmrd.at(i_run));
    gr_rate_pmtmrdfmv->SetPoint(i_run,run_numbers.at(i_run),rate_pmtmrdfmv.at(i_run));
    gr_rate_pmtmrdnofmv->SetPoint(i_run,run_numbers.at(i_run),rate_pmtmrdnofmv.at(i_run));
    gr_rate_fmv->SetPoint(i_run,run_numbers.at(i_run),rate_fmv.at(i_run));
    gr_rate_mrd->SetPoint(i_run,run_numbers.at(i_run),rate_mrd.at(i_run));
    gr_frac->SetPoint(i_run,run_numbers.at(i_run),frac.at(i_run));
    gr_frac_pmt->SetPoint(i_run,run_numbers.at(i_run),frac_pmt.at(i_run));
    gr_frac_pmtmrd->SetPoint(i_run,run_numbers.at(i_run),frac_pmtmrd.at(i_run));
    gr_frac_pmtmrdfmv->SetPoint(i_run,run_numbers.at(i_run),frac_pmtmrdfmv.at(i_run));
    gr_frac_pmtmrdnofmv->SetPoint(i_run,run_numbers.at(i_run),frac_pmtmrdnofmv.at(i_run));
    gr_frac_fmv->SetPoint(i_run,run_numbers.at(i_run),frac_fmv.at(i_run));
    gr_frac_mrd->SetPoint(i_run,run_numbers.at(i_run),frac_mrd.at(i_run));
    gr_align_pmt_mrd->SetPoint(i_run,run_numbers.at(i_run),align_pmt_mrd.at(i_run));
    gr_align_pmt_fmv->SetPoint(i_run,run_numbers.at(i_run),align_pmt_fmv.at(i_run));
    gr_align_mrd_fmv->SetPoint(i_run,run_numbers.at(i_run),align_mrd_fmv.at(i_run));
    gr_frac_prompt->SetPoint(i_run,run_numbers.at(i_run),frac_prompt.at(i_run));
    gr_frac_delayed->SetPoint(i_run,run_numbers.at(i_run),frac_delayed.at(i_run));
    gr_frac_beam->SetPoint(i_run,run_numbers.at(i_run),frac_beam.at(i_run));
    gr_frac_led->SetPoint(i_run,run_numbers.at(i_run),frac_led.at(i_run));
    gr_frac_cosmic->SetPoint(i_run,run_numbers.at(i_run),frac_cosmic.at(i_run));
    gr_avg_mult->SetPoint(i_run,run_numbers.at(i_run),avg_mult.at(i_run));
    gr_avg_mult_coinc->SetPoint(i_run,run_numbers.at(i_run),avg_mult_coinc.at(i_run));
    gr_avg_mult_coinc_nofmv_cb->SetPoint(i_run,run_numbers.at(i_run),avg_mult_coinc_nofmv_cb.at(i_run));
    gr_trigger_beam->SetPoint(i_run,run_numbers.at(i_run),rate_beam.at(i_run));
    gr_trigger_extended->SetPoint(i_run,run_numbers.at(i_run),rate_extended.at(i_run));
    gr_trigger_led->SetPoint(i_run,run_numbers.at(i_run),rate_led.at(i_run));
    gr_trigger_cosmic->SetPoint(i_run,run_numbers.at(i_run),rate_cosmic.at(i_run));
    gr_charge_prompt->SetPoint(i_run,run_numbers.at(i_run),avg_charge_prompt.at(i_run));
    gr_charge_del->SetPoint(i_run,run_numbers.at(i_run),avg_charge_del.at(i_run));
    gr_charge_per_pmt->SetPoint(i_run,run_numbers.at(i_run),avg_charge_per_pmt.at(i_run));
    gr_pot->SetPoint(i_run,run_numbers.at(i_run),pot.at(i_run));
    gr_pot_scaled->SetPoint(i_run,run_numbers.at(i_run),pot.at(i_run)/1.e14);
  }

  TMultiGraph *multi_rate = new TMultiGraph();
  TMultiGraph *multi_trigger = new TMultiGraph();
  TMultiGraph *multi_frac = new TMultiGraph();
  TMultiGraph *multi_types = new TMultiGraph();
  TMultiGraph *multi_mult = new TMultiGraph();
  TMultiGraph *multi_pot = new TMultiGraph();
  TLegend *leg_rate = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_trigger = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_frac = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_types = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_mult = new TLegend(0.5,0.6,0.88,0.88);
  TLegend *leg_ref_timing_1 = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_ref_timing_2 = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_ref_timing_3 = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_ref_charge_1 = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_ref_charge_2 = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_ref_charge_3 = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_ref_trigger_1 = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_ref_trigger_2 = new TLegend(0.7,0.7,0.88,0.88);
  TLegend *leg_pot = new TLegend(0.7,0.7,0.88,0.88);
  TCanvas *canvas_rate = new TCanvas("canvas_rate","Rate Multigraph",900,600);
  TCanvas *canvas_trigger = new TCanvas("canvas_trigger","Trigger Rate Multigraph",900,600);
  TCanvas *canvas_frac = new TCanvas("canvas_frac","Fractions Multigraph",900,600);
  TCanvas *canvas_types = new TCanvas("canvas_types","Trigger Types Multigraph",900,600);
  TCanvas *canvas_mult = new TCanvas("canvas_mult","Multiplicity Multigraph",900,600);
  TCanvas *canvas_align = new TCanvas("canvas_align","Time Alignment evolution",900,600);
  TCanvas *canvas_charge = new TCanvas("canvas_charge","Charge evolution",900,600);
  TCanvas *canvas_ref_timing = new TCanvas("canvas_ref_timing","Reference timing",900,600);
  TCanvas *canvas_ref_charge = new TCanvas("canvas_ref_charge","Reference charge",900,600);
  TCanvas *canvas_ref_trigger = new TCanvas("canvas_ref_trigger","Reference trigger",900,600);
  TCanvas *canvas_pot = new TCanvas("canvas_pot","POT Multigraph",900,600);

  std::stringstream ss_multi_rate, ss_multi_trigger, ss_multi_frac, ss_multi_types, ss_multi_mult, ss_multi_pot;
  ss_multi_rate << "Event Rates ("<<string_week_start << " - " << string_week_end << ")";
  ss_multi_trigger << "Trigger Rates ("<<string_week_start << " - " << string_week_end << ")";
  ss_multi_frac << "Event Fractions ("<<string_week_start << " - " << string_week_end << ")";
  ss_multi_types << "Event Types ("<<string_week_start << " - " << string_week_end << ")";
  ss_multi_mult << "Neutron Multiplicities ("<<string_week_start << " - " << string_week_end << ")";
  ss_multi_pot << "PMT + POT Rates ("<<string_week_start << " - " << string_week_end << ")";
 
  multi_rate->SetTitle(ss_multi_rate.str().c_str());
  multi_trigger->SetTitle(ss_multi_trigger.str().c_str());
  multi_frac->SetTitle(ss_multi_frac.str().c_str());
  multi_types->SetTitle(ss_multi_types.str().c_str());
  multi_mult->SetTitle(ss_multi_mult.str().c_str());
  multi_pot->SetTitle(ss_multi_pot.str().c_str());


  gr_rate->SetLineColor(1);
  gr_rate->SetMarkerColor(1);
  multi_rate->Add(gr_rate);
  leg_rate->AddEntry(gr_rate,"Event rate","l");
  gr_rate_pmt->SetLineColor(kOrange);
  gr_rate_pmt->SetMarkerColor(kOrange);
  multi_rate->Add(gr_rate_pmt);
  leg_rate->AddEntry(gr_rate_pmt,"PMT Cluster rate","l");
  gr_rate_pmtmrd->SetLineColor(8);
  gr_rate_pmtmrd->SetMarkerColor(8);
  multi_rate->Add(gr_rate_pmtmrd);
  leg_rate->AddEntry(gr_rate_pmtmrd,"PMT+MRD Cluster rate","l");
  gr_rate_pmtmrdfmv->SetLineColor(kRed);
  gr_rate_pmtmrdfmv->SetMarkerColor(kRed);
  multi_rate->Add(gr_rate_pmtmrdfmv);
  leg_rate->AddEntry(gr_rate_pmtmrdfmv,"PMT+MRD+FMV Cluster rate","l");
  gr_rate_pmtmrdnofmv->SetLineColor(kBlue);
  gr_rate_pmtmrdnofmv->SetMarkerColor(kBlue);
  multi_rate->Add(gr_rate_pmtmrdnofmv);
  leg_rate->AddEntry(gr_rate_pmtmrdnofmv,"PMT_MRD (No FMV) Cluster rate","l");
  canvas_rate->cd();
  multi_rate->Draw("apl");
  multi_rate->GetYaxis()->SetTitle("Rate [Hz]");
  multi_rate->GetXaxis()->SetTitle("Run number");
  leg_rate->Draw();
  multi_rate->GetYaxis()->SetRangeUser(1e-5,5);
  canvas_rate->SetLogy();

  gr_trigger_beam->SetLineColor(1);
  gr_trigger_beam->SetMarkerColor(1);
  multi_trigger->Add(gr_trigger_beam);
  leg_trigger->AddEntry(gr_trigger_beam,"Beam rate","l");
  gr_trigger_cosmic->SetLineColor(kOrange);
  gr_trigger_cosmic->SetMarkerColor(kOrange);
  multi_trigger->Add(gr_trigger_cosmic);
  leg_trigger->AddEntry(gr_trigger_cosmic,"Cosmic rate","l");
  gr_trigger_led->SetLineColor(8);
  gr_trigger_led->SetMarkerColor(8);
  multi_trigger->Add(gr_trigger_led);
  leg_trigger->AddEntry(gr_trigger_led,"LED rate","l");
  gr_trigger_extended->SetLineColor(kRed);
  gr_trigger_extended->SetMarkerColor(kRed);
  multi_trigger->Add(gr_trigger_extended);
  leg_trigger->AddEntry(gr_trigger_extended,"Extended rate","l");
  canvas_trigger->cd();
  multi_trigger->Draw("apl");
  multi_trigger->GetYaxis()->SetTitle("Trigger Rate [Hz]");
  multi_trigger->GetXaxis()->SetTitle("Run number");
  leg_trigger->Draw();
  multi_trigger->GetYaxis()->SetRangeUser(1e-3,100);
  canvas_trigger->SetLogy();

  gr_frac_pmt->SetLineColor(kOrange);
  gr_frac_pmt->SetMarkerColor(kOrange);
  multi_frac->Add(gr_frac);
  leg_frac->AddEntry(gr_frac,"Event fraction","l");
  multi_frac->Add(gr_frac_pmt);
  leg_frac->AddEntry(gr_frac_pmt,"PMT Cluster fraction","l");
  gr_frac_pmtmrd->SetLineColor(8);
  gr_frac_pmtmrd->SetMarkerColor(8);
  multi_frac->Add(gr_frac_pmtmrd);
  leg_frac->AddEntry(gr_frac_pmtmrd,"PMT+MRD Cluster fraction","l");
  gr_frac_pmtmrdfmv->SetLineColor(kRed);
  gr_frac_pmtmrdfmv->SetMarkerColor(kRed);
  multi_frac->Add(gr_frac_pmtmrdfmv);
  leg_frac->AddEntry(gr_frac_pmtmrdfmv,"PMT+MRD+FMV Cluster fraction","l");
  gr_frac_pmtmrdnofmv->SetLineColor(kBlue);
  gr_frac_pmtmrdnofmv->SetMarkerColor(kBlue);
  multi_frac->Add(gr_frac_pmtmrdnofmv);
  leg_frac->AddEntry(gr_frac_pmtmrdnofmv,"PMT_MRD (No FMV) Cluster fraction","l");
  canvas_frac->cd();
  multi_frac->Draw("apl");
  multi_frac->GetYaxis()->SetTitle("Fraction");
  multi_frac->GetXaxis()->SetTitle("Run Number");
  leg_frac->Draw();
  multi_frac->GetYaxis()->SetRangeUser(1e-4,110);
  canvas_frac->SetLogy();

  gr_frac_prompt->SetLineColor(kBlue);
  gr_frac_prompt->SetMarkerColor(kBlue);
  gr_frac_delayed->SetLineColor(kViolet);
  gr_frac_delayed->SetMarkerColor(kViolet);
  gr_frac_beam->SetLineColor(kOrange);
  gr_frac_beam->SetMarkerColor(kOrange);
  gr_frac_led->SetLineColor(kRed);
  gr_frac_led->SetMarkerColor(kRed);
  gr_frac_cosmic->SetLineColor(8);
  gr_frac_cosmic->SetMarkerColor(8);
  multi_types->Add(gr_frac_prompt);
  leg_types->AddEntry(gr_frac_prompt,"prompt window","l");
  multi_types->Add(gr_frac_delayed);
  leg_types->AddEntry(gr_frac_delayed,"delayed window","l");
  multi_types->Add(gr_frac_beam);
  leg_types->AddEntry(gr_frac_beam,"Beam Trigger","l");
  multi_types->Add(gr_frac_led);
  leg_types->AddEntry(gr_frac_led,"LED Trigger","l");
  multi_types->Add(gr_frac_cosmic);
  leg_types->AddEntry(gr_frac_cosmic,"Cosmic Trigger","l");
  canvas_types->cd();
  multi_types->Draw("apl");
  multi_types->GetYaxis()->SetTitle("Fraction");
  multi_types->GetXaxis()->SetTitle("Run Number");
  leg_types->Draw();

  gr_avg_mult->SetLineColor(1);
  gr_avg_mult->SetMarkerColor(1);
  gr_avg_mult_coinc->SetLineColor(kBlue);
  gr_avg_mult_coinc->SetMarkerColor(kBlue);
  gr_avg_mult_coinc_nofmv_cb->SetLineColor(8);
  gr_avg_mult_coinc_nofmv_cb->SetMarkerColor(8);
  multi_mult->Add(gr_avg_mult);
  leg_mult->AddEntry(gr_avg_mult,"Average neutron multiplicity","l");
  multi_mult->Add(gr_avg_mult_coinc);
  leg_mult->AddEntry(gr_avg_mult_coinc,"Average neutron multiplicity (MRD+Tank coincidence)","l");
  multi_mult->Add(gr_avg_mult_coinc_nofmv_cb);
  leg_mult->AddEntry(gr_avg_mult_coinc_nofmv_cb,"Average neutron multiplicity (MRD+Tank, no FMV)","l");
  canvas_mult->cd();
  multi_mult->Draw("apl");
  multi_mult->GetYaxis()->SetTitle("Multiplicity");
  multi_mult->GetXaxis()->SetTitle("Run Number");
  leg_mult->Draw();

  gr_pot_scaled->SetLineColor(kBlue);
  gr_pot_scaled->SetMarkerColor(kBlue);
  multi_pot->Add(gr_rate_pmt);
  leg_pot->AddEntry(gr_rate_pmt,"PMT Cluster Rate","l");
  multi_pot->Add(gr_pot_scaled);
  leg_pot->AddEntry(gr_pot_scaled,"POT Rate (scaled)","l");
  canvas_pot->cd();
  multi_pot->Draw("apl");
  multi_pot->GetYaxis()->SetTitle("Rate");
  multi_pot->GetXaxis()->SetTitle("Run Number");
  leg_pot->Draw(); 

  canvas_ref_timing->Divide(2,2);
  canvas_ref_timing->cd(1);
  PMT_t_clusters_2pe_reference->Scale(PMT_t_clusters_2pe_combined->Integral()/PMT_t_clusters_2pe_reference->Integral());
  PMT_t_clusters_2pe_combined->SetStats(0);
  PMT_t_clusters_2pe_combined->Draw();
  leg_ref_timing_1->AddEntry(PMT_t_clusters_2pe_combined,"Current week");
  PMT_t_clusters_2pe_reference->SetStats(0);
  PMT_t_clusters_2pe_reference->Draw("same");
  leg_ref_timing_1->AddEntry(PMT_t_clusters_2pe_reference,"Reference (R3082)");
  leg_ref_timing_1->Draw();
  canvas_ref_timing->cd(2);
  MRD_t_clusters_beam_reference->Scale(MRD_t_clusters_beam_combined->Integral()/MRD_t_clusters_beam_reference->Integral());
  MRD_t_clusters_beam_combined->SetStats(0);
  MRD_t_clusters_beam_combined->Draw();
  leg_ref_timing_2->AddEntry(MRD_t_clusters_beam_combined,"Current week");
  MRD_t_clusters_beam_reference->SetStats(0);
  MRD_t_clusters_beam_reference->Draw("same");
  leg_ref_timing_2->AddEntry(MRD_t_clusters_beam_reference,"Reference (R3082)");
  leg_ref_timing_2->Draw();
  canvas_ref_timing->cd(3);
  FMV_t_clusters_beam_reference->Scale(FMV_t_clusters_beam_combined->Integral()/FMV_t_clusters_beam_reference->Integral());
  FMV_t_clusters_beam_combined->SetStats(0);
  FMV_t_clusters_beam_combined->Draw();
  leg_ref_timing_3->AddEntry(FMV_t_clusters_beam_combined,"Current week");
  FMV_t_clusters_beam_reference->SetStats(0);
  FMV_t_clusters_beam_reference->Draw("same");
  leg_ref_timing_3->AddEntry(FMV_t_clusters_beam_reference,"Reference (R3082)");
  leg_ref_timing_3->Draw();
 
  canvas_ref_charge->Divide(2,2);
  TPad *p0 = (TPad*) canvas_ref_charge->cd(1);
  PMT_prompt_charge_reference->Scale(PMT_prompt_charge_combined->Integral()/PMT_prompt_charge_reference->Integral());
  PMT_prompt_charge_combined->SetStats(0);
  PMT_prompt_charge_combined->Draw();
  leg_ref_charge_1->AddEntry(PMT_prompt_charge_combined,"Current week");
  PMT_prompt_charge_reference->SetStats(0);
  PMT_prompt_charge_reference->Draw("same");
  leg_ref_charge_1->AddEntry(PMT_prompt_charge_reference,"Reference (R3082)");
  p0->SetLogy();
  leg_ref_charge_1->Draw();
  canvas_ref_charge->cd(2);
  PMT_delayed_charge_zoom_reference->Scale(PMT_delayed_charge_zoom_combined->Integral()/PMT_delayed_charge_zoom_reference->Integral());
  PMT_delayed_charge_zoom_combined->SetStats(0);
  PMT_delayed_charge_zoom_combined->Draw();
  leg_ref_charge_2->AddEntry(PMT_delayed_charge_zoom_combined,"Current week");
  PMT_delayed_charge_zoom_reference->SetStats(0);
  PMT_delayed_charge_zoom_reference->Draw("same");
  leg_ref_charge_2->AddEntry(PMT_delayed_charge_zoom_reference,"Reference (R3082)");
  leg_ref_charge_2->Draw();
  canvas_ref_charge->cd(3);
  PMT_chargeperpmt_reference->Scale(PMT_chargeperpmt_combined->Integral()/PMT_chargeperpmt_reference->Integral());
  PMT_chargeperpmt_reference->SetStats(0);
  PMT_chargeperpmt_combined->Draw();
  leg_ref_charge_3->AddEntry(PMT_chargeperpmt_combined,"Current week");
  PMT_chargeperpmt_reference->SetStats(0);
  PMT_chargeperpmt_reference->Draw("same");
  leg_ref_charge_3->AddEntry(PMT_chargeperpmt_reference,"Reference (R3082)");
  leg_ref_charge_3->Draw();

  canvas_ref_trigger->Divide(2,1);
  TPad *p1 = (TPad*) canvas_ref_trigger->cd(1);
  ADCWaveform_Samples_reference->Scale(ADCWaveform_Samples_combined->Integral()/ADCWaveform_Samples_reference->Integral());
  ADCWaveform_Samples_combined->SetStats(0);
  ADCWaveform_Samples_combined->Draw();
  leg_ref_trigger_1->AddEntry(ADCWaveform_Samples_combined,"Current week");
  ADCWaveform_Samples_reference->SetStats(0);
  ADCWaveform_Samples_reference->Draw("same");
  leg_ref_trigger_1->AddEntry(ADCWaveform_Samples_reference,"Reference (R3082)");
  p1->SetLogy();
  leg_ref_trigger_1->Draw();
  canvas_ref_trigger->cd(2);
  Triggerwords_reference->Scale(Triggerwords_combined->Integral()/Triggerwords_reference->Integral());
  Triggerwords_combined->SetStats(0);
  Triggerwords_combined->Draw();
  leg_ref_trigger_2->AddEntry(Triggerwords_combined,"Current week");
  Triggerwords_reference->SetStats(0);
  Triggerwords_reference->Draw("same");
  leg_ref_trigger_2->AddEntry(Triggerwords_reference,"Reference (R3082)");
  leg_ref_trigger_2->Draw();

  canvas_align->Divide(2,2);
  canvas_align->cd(1);
  gr_align_pmt_mrd->Draw();
  canvas_align->cd(2);
  gr_align_pmt_fmv->Draw();
  canvas_align->cd(3);
  gr_align_mrd_fmv->Draw();

  canvas_charge->Divide(2,2);
  TPad *pc1 = (TPad*) canvas_charge->cd(1);
  gr_charge_prompt->Draw();
  canvas_charge->cd(2);
  gr_charge_del->Draw();
  canvas_charge->cd(3);
  gr_charge_per_pmt->Draw();

  if (SaveImages){
    std::stringstream out_rate, out_trigger, out_frac, out_types, out_mult, out_align, out_charge, out_ref_timing, out_ref_charge, out_ref_trigger, out_pot;

  if (WeeklyMode) filename_out << "RunValidationStability_"<<str_week_start<<"-"<<str_week_end<<".root";
  else filename_out << "RunValidationStability_R" << current_run - number_runs << "-" << current_run << ".root";

  if (WeeklyMode) {
    out_rate << OutputPath << "/RunValidationStability_Rate_"<<str_week_start<<"-"<<str_week_end<<".png";
    out_trigger << OutputPath << "/RunValidationStability_Trigger_"<<str_week_start << "-" << str_week_end << ".png";
    out_frac << OutputPath << "/RunValidationStability_Fractions_"<<str_week_start << "-" << str_week_end << ".png";
    out_types << OutputPath << "/RunValidationStability_Types_"<<str_week_start << "-" << str_week_end << ".png";
    out_mult << OutputPath << "/RunValidationStability_Mult_"<<str_week_start << "-" << str_week_end << ".png";
    out_align << OutputPath << "/RunValidationStability_Align_"<<str_week_start << "-" << str_week_end << ".png";
    out_charge << OutputPath << "/RunValidationStability_Charge_"<<str_week_start << "-" << str_week_end << ".png";
    out_ref_timing << OutputPath << "/RunValidationStability_RefTiming_"<<str_week_start << "-" << str_week_end << ".png";
    out_ref_charge << OutputPath << "/RunValidationStability_RefCharge_"<<str_week_start << "-" << str_week_end << ".png";
    out_ref_trigger << OutputPath << "/RunValidationStability_RefTrigger_"<<str_week_start << "-" << str_week_end << ".png";
    out_pot << OutputPath << "/RunValidationStability_POT_"<<str_week_start << "-" << str_week_end << ".png";
  } else {
    out_rate << OutputPath << "/RunValidationStability_Rate_"<< current_run - number_runs << "-" << current_run << ".png";
    out_trigger << OutputPath << "/RunValidationStability_Trigger_"<< current_run - number_runs << "-" << current_run << ".png";
    out_frac << OutputPath << "/RunValidationStability_Fractions_"<< current_run - number_runs << "-" << current_run << ".png";
    out_types << OutputPath << "/RunValidationStability_Types_"<< current_run - number_runs << "-" << current_run << ".png";
    out_mult << OutputPath << "/RunValidationStability_Mult_"<< current_run - number_runs << "-" << current_run << ".png";
    out_align << OutputPath << "/RunValidationStability_Align_"<< current_run - number_runs << "-" << current_run << ".png";
    out_charge << OutputPath << "/RunValidationStability_Charge_"<< current_run - number_runs << "-" << current_run << ".png";
    out_ref_timing << OutputPath << "/RunValidationStability_RefTiming_"<< current_run - number_runs << "-" << current_run << ".png";
    out_ref_charge << OutputPath << "/RunValidationStability_RefCharge_"<< current_run - number_runs << "-" << current_run << ".png";
    out_ref_trigger << OutputPath << "/RunValidationStability_RefTrigger_"<< current_run - number_runs << "-" << current_run << ".png";
    out_pot << OutputPath << "/RunValidationStability_POT_"<< current_run - number_runs << "-" << current_run << ".png";
  }
  
  canvas_rate->SaveAs(out_rate.str().c_str());
  canvas_trigger->SaveAs(out_trigger.str().c_str());
  canvas_frac->SaveAs(out_frac.str().c_str());
  canvas_types->SaveAs(out_types.str().c_str());
  canvas_mult->SaveAs(out_mult.str().c_str());
  canvas_align->SaveAs(out_align.str().c_str());
  canvas_charge->SaveAs(out_charge.str().c_str());
  canvas_ref_timing->SaveAs(out_ref_timing.str().c_str());
  canvas_ref_charge->SaveAs(out_ref_charge.str().c_str());
  canvas_ref_trigger->SaveAs(out_ref_trigger.str().c_str());
  canvas_pot->SaveAs(out_pot.str().c_str());

  }

/*
  FMV_t_clusters_beam_reference->Write();
  ADCWaveform_Samples_reference->Write();
  Triggerwords_reference->Write();
  FMV_PMT_Deltat_100pe_Beam_reference->Write();
  //MRD_PMT_Deltat_100pe_Beam_reference->Write();
  PMT_prompt_charge_reference->Write();
  PMT_delayed_charge_zoom_reference->Write();
  PMT_chargeperpmt_reference->Write();
*/
  gr_rate->GetXaxis()->SetTitle("Run Number");
  gr_rate_pmt->GetXaxis()->SetTitle("Run Number");
  gr_rate_pmtmrd->GetXaxis()->SetTitle("Run Number");
  gr_rate_pmtmrdfmv->GetXaxis()->SetTitle("Run Number");
  gr_rate_pmtmrdnofmv->GetXaxis()->SetTitle("Run Number");
  gr_rate_fmv->GetXaxis()->SetTitle("Run Number");
  gr_rate_mrd->GetXaxis()->SetTitle("Run Number");
  gr_frac->GetXaxis()->SetTitle("Run Number");
  gr_frac_pmt->GetXaxis()->SetTitle("Run Number");
  gr_frac_pmtmrd->GetXaxis()->SetTitle("Run Number");
  gr_frac_pmtmrdfmv->GetXaxis()->SetTitle("Run Number");
  gr_frac_pmtmrdnofmv->GetXaxis()->SetTitle("Run Number");
  gr_frac_fmv->GetXaxis()->SetTitle("Run Number");
  gr_frac_mrd->GetXaxis()->SetTitle("Run Number");
  gr_align_pmt_mrd->GetXaxis()->SetTitle("Run Number");
  gr_align_pmt_fmv->GetXaxis()->SetTitle("Run Number");
  gr_align_mrd_fmv->GetXaxis()->SetTitle("Run Number");
  gr_frac_prompt->GetXaxis()->SetTitle("Run Number");
  gr_frac_delayed->GetXaxis()->SetTitle("Run Number");
  gr_frac_beam->GetXaxis()->SetTitle("Run Number");
  gr_frac_led->GetXaxis()->SetTitle("Run Number");
  gr_frac_cosmic->GetXaxis()->SetTitle("Run Number");
  gr_avg_mult->GetXaxis()->SetTitle("Run Number");
  gr_avg_mult_coinc->GetXaxis()->SetTitle("Run Number");
  gr_avg_mult_coinc_nofmv_cb->GetXaxis()->SetTitle("Run Number");
  gr_trigger_beam->GetXaxis()->SetTitle("Run Number");
  gr_trigger_extended->GetXaxis()->SetTitle("Run Number");
  gr_trigger_cosmic->GetXaxis()->SetTitle("Run Number");
  gr_trigger_led->GetXaxis()->SetTitle("Run Number");
  gr_charge_prompt->GetXaxis()->SetTitle("Run Number");
  gr_charge_del->GetXaxis()->SetTitle("Run Number");
  gr_charge_per_pmt->GetXaxis()->SetTitle("Run Number");
  gr_pot->GetXaxis()->SetTitle("Run Number");
  gr_pot_scaled->GetXaxis()->SetTitle("Run Number");

  gr_rate->GetYaxis()->SetTitle("R [Hz]");
  gr_rate_pmt->GetYaxis()->SetTitle("R [Hz]");
  gr_rate_pmtmrd->GetYaxis()->SetTitle("R [Hz]");
  gr_rate_pmtmrdfmv->GetYaxis()->SetTitle("R [Hz]");
  gr_rate_pmtmrdnofmv->GetYaxis()->SetTitle("R [Hz]");
  gr_rate_fmv->GetYaxis()->SetTitle("R [Hz]");
  gr_rate_mrd->GetYaxis()->SetTitle("R [Hz]");
  gr_frac->GetYaxis()->SetTitle("Fraction");
  gr_frac_pmt->GetYaxis()->SetTitle("Fraction");
  gr_frac_pmtmrd->GetYaxis()->SetTitle("Fraction");
  gr_frac_pmtmrdfmv->GetYaxis()->SetTitle("Fraction");
  gr_frac_pmtmrdnofmv->GetYaxis()->SetTitle("Fraction");
  gr_frac_fmv->GetYaxis()->SetTitle("Fraction");
  gr_frac_mrd->GetYaxis()->SetTitle("Fraction");
  gr_align_pmt_mrd->GetYaxis()->SetTitle("#Delta t_{PMT-MRD}");
  gr_align_pmt_fmv->GetYaxis()->SetTitle("#Delta t_{PMT-FMV}");
  gr_align_mrd_fmv->GetYaxis()->SetTitle("#Delta t_{MRD-FMV}");
  gr_frac_prompt->GetYaxis()->SetTitle("Fraction");
  gr_frac_delayed->GetYaxis()->SetTitle("Fraction");
  gr_frac_beam->GetYaxis()->SetTitle("Fraction");
  gr_frac_led->GetYaxis()->SetTitle("Fraction");
  gr_frac_cosmic->GetYaxis()->SetTitle("Fraction");
  gr_avg_mult->GetYaxis()->SetTitle("Neutron multiplicity");
  gr_avg_mult_coinc->GetYaxis()->SetTitle("Neutron multiplicity");
  gr_avg_mult_coinc_nofmv_cb->GetYaxis()->SetTitle("Neutron multiplicity");
  gr_trigger_beam->GetYaxis()->SetTitle("R_{beam} [Hz]");
  gr_trigger_extended->GetYaxis()->SetTitle("R_{extended} [Hz]");
  gr_trigger_led->GetYaxis()->SetTitle("R_{led} [Hz]");
  gr_trigger_cosmic->GetYaxis()->SetTitle("R_{cosmic} [Hz]");
  gr_charge_prompt->GetYaxis()->SetTitle("Q_{prompt} [p.e.]");
  gr_charge_del->GetYaxis()->SetTitle("Q_{delayed} [p.e.]");
  gr_charge_per_pmt->GetYaxis()->SetTitle("Q/n_{PMT} [p.e.]");
  gr_pot->GetYaxis()->SetTitle("POT/s");
  gr_pot_scaled->GetYaxis()->SetTitle("POT Rate [arbitrary scale]");

  std::stringstream ss_gr_rate, ss_gr_rate_pmt, ss_gr_rate_pmtmrd, ss_gr_rate_pmtmrdfmv, ss_gr_rate_pmtmrdnofmv;
  std::stringstream ss_gr_rate_fmv, ss_gr_rate_mrd, ss_gr_frac, ss_gr_frac_pmt, ss_gr_frac_pmtmrd, ss_gr_frac_pmtmrdfmv, ss_gr_frac_pmtmrdnofmv;
  std::stringstream ss_gr_frac_fmv, ss_gr_frac_mrd, ss_gr_align_pmt_mrd, ss_gr_align_pmt_fmv, ss_gr_align_mrd_fmv, ss_gr_frac_prompt, ss_gr_frac_delayed, ss_gr_frac_beam, ss_gr_frac_led, ss_gr_frac_cosmic;
  std::stringstream ss_gr_avg_mult, ss_gr_avg_mult_coinc, ss_gr_avg_mult_coinc_nofmv_cb;
  std::stringstream ss_gr_trigger_beam, ss_gr_trigger_extended, ss_gr_trigger_led, ss_gr_trigger_cosmic;  
  std::stringstream ss_gr_charge_prompt, ss_gr_charge_del, ss_gr_charge_per_pmt;
  std::stringstream ss_gr_pot, ss_gr_pot_scaled;

  if (!WeeklyMode){
  ss_gr_rate << "Rate Beam Events (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_rate_pmt << "Rate PMT Clusters [Beam Events] (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_rate_pmtmrd << "Rate PMT + MRD Clusters [Beam Events] (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_rate_pmtmrdfmv << "Rate PMT + MRD + FMV [Beam Events] (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_rate_pmtmrdnofmv << "Rate PMT + MRD + No FMV [Beam Events] (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_rate_fmv << "Rate FMV Clusters [Beam Events] (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_rate_mrd << "Rate MRD Clusters [Beam Events] (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_frac << "Fraction Beam Events (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_frac_pmt << "Fraction PMT Clusters [Beam Events] (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_frac_pmtmrd << "Fraction PMT+MRD Clusters [Beam Events] (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_frac_pmtmrdfmv << "Fraction PMT+MRD+FMV Clusters [Beam Events] (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_frac_pmtmrdnofmv << "Fraction PMT+MRD+No FMV Clusters [Beam Events] (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_frac_fmv << "Fraction FMV Clusters [Beam Events] (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_frac_mrd << "Fraction MRD Clusters [Beam Events] (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_align_pmt_mrd << "Alignment PMT/MRD timing (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_align_pmt_fmv << "Alignment PMT/FMV timing (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_align_mrd_fmv << "Alignment MRD/FMV timing (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_frac_prompt << "Fraction Prompt Acquisition [Beam Events] (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_frac_delayed << "Fraction Delayed Acquisition [Beam Events] (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_frac_beam << "Fraction Beam Events (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_frac_led << "Fraction LED Events (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_frac_cosmic << "Fraction Cosmic Events (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_avg_mult << "Average neutron multiplicity (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_avg_mult_coinc << "Average neutron multiplicity (MRD/PMT coincidence) (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_avg_mult_coinc_nofmv_cb << "Average neutron multiplicity (No FMV, CB < 0.4) (Run "<<current_run - number_runs << "- Run "<<current_run<<")";
  ss_gr_trigger_beam << "Beam Trigger Rate (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_gr_trigger_cosmic << "Cosmic Trigger Rate (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_gr_trigger_led << "LED Trigger Rate (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_gr_trigger_extended << "Extended Trigger Rate (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_gr_charge_prompt << "Average prompt charge (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_gr_charge_del << "Average delayed charge (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_gr_charge_per_pmt << "Average charge per PMT (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_gr_pot << "POT Rate (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_gr_pot_scaled << "POT Scaled Rate (Run "<<current_run - number_runs << "- Run "<<current_run << ")";
  ss_multi_rate << "Event Rates (Run "<< current_run - number_runs  << " - Run " << current_run << ")";
  ss_multi_trigger << "Trigger Rates (Run "<< current_run - number_runs << " - Run " << current_run << ")";
  ss_multi_frac << "Event Fractions (Run "<< current_run - number_runs << " - Run " << current_run << ")";
  ss_multi_types << "Event Types (Run "<< current_run - number_runs << " - Run " << current_run << ")";
  ss_multi_mult << "Neutron Multiplicities (Run "<<current_run - number_runs << " - Run " << current_run << ")";
  ss_multi_pot << "PMT + POT Rates (Run "<<current_run - number_runs << " - Run " << current_run << ")";
  } else {  
    ss_gr_rate << "Rate Beam Events ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_rate_pmt << "Rate PMT Clusters [Beam Events] ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_rate_pmtmrd << "Rate PMT + MRD Clusters [Beam Events] ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_rate_pmtmrdfmv << "Rate PMT + MRD + FMV [Beam Events] ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_rate_pmtmrdnofmv << "Rate PMT + MRD + No FMV [Beam Events] ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_rate_fmv << "Rate FMV Clusters [Beam Events] ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_rate_mrd << "Rate MRD Clusters [Beam Events] ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_frac << "Fraction Beam Events ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_frac_pmt << "Fraction PMT Clusters [Beam Events] ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_frac_pmtmrd << "Fraction PMT+MRD Clusters [Beam Events] ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_frac_pmtmrdfmv << "Fraction PMT+MRD+FMV Clusters [Beam Events] ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_frac_pmtmrdnofmv << "Fraction PMT+MRD+No FMV Clusters [Beam Events] ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_frac_fmv << "Fraction FMV Clusters [Beam Events] ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_frac_mrd << "Fraction MRD Clusters [Beam Events] ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_align_pmt_mrd << "Alignment PMT/MRD timing ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_align_pmt_fmv << "Alignment PMT/FMV timing ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_align_mrd_fmv << "Alignment MRD/FMV timing ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_frac_prompt << "Fraction Prompt Acquisition [Beam Events] ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_frac_delayed << "Fraction Delayed Acquisition [Beam Events] ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_frac_beam << "Fraction Beam Events ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_frac_led << "Fraction LED Events ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_frac_cosmic << "Fraction Cosmic Events ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_avg_mult << "Average neutron multiplicity ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_avg_mult_coinc << "Average neutron multiplicity (MRD/PMT coincidence) ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_avg_mult_coinc_nofmv_cb << "Average neutron multiplicity (No FMV, CB < 0.4) ( "<<string_week_start << " - "<<string_week_end<<")";
  ss_gr_trigger_beam << "Beam Trigger Rate ( "<<string_week_start << " - "<<string_week_end << ")";
  ss_gr_trigger_cosmic << "Cosmic Trigger Rate ( "<<string_week_start << " - "<<string_week_end << ")";
  ss_gr_trigger_led << "LED Trigger Rate ( "<<string_week_start << " - "<<string_week_end << ")";
  ss_gr_trigger_extended << "Extended Trigger Rate ( "<< string_week_start << " - "<<string_week_end << ")";
  ss_gr_charge_prompt << "Average prompt charge ( "<< string_week_start << " - "<<string_week_end << ")";
  ss_gr_charge_del << "Average delayed charge ( "<<string_week_start << " - "<<string_week_end << ")";
  ss_gr_charge_per_pmt << "Average charge per PMT ( "<<string_week_start << " - "<<string_week_end << ")";
  ss_gr_pot << "POT Rate ( "<<string_week_start << " - "<<string_week_end << ")";
  ss_gr_pot_scaled << "POT Scaled Rate ( "<<string_week_start << " - "<< string_week_end << ")";
  }

  gr_rate->SetTitle(ss_gr_rate.str().c_str());
  gr_rate_pmt->SetTitle(ss_gr_rate_pmt.str().c_str());
  gr_rate_pmtmrd->SetTitle(ss_gr_rate_pmtmrd.str().c_str());
  gr_rate_pmtmrdfmv->SetTitle(ss_gr_rate_pmtmrdfmv.str().c_str());
  gr_rate_pmtmrdnofmv->SetTitle(ss_gr_rate_pmtmrdnofmv.str().c_str());
  gr_rate_fmv->SetTitle(ss_gr_rate_fmv.str().c_str());
  gr_rate_mrd->SetTitle(ss_gr_rate_mrd.str().c_str());
  gr_frac->SetTitle(ss_gr_frac.str().c_str());
  gr_frac_pmt->SetTitle(ss_gr_frac_pmt.str().c_str());
  gr_frac_pmtmrd->SetTitle(ss_gr_frac_pmtmrd.str().c_str());
  gr_frac_pmtmrdfmv->SetTitle(ss_gr_frac_pmtmrdfmv.str().c_str());
  gr_frac_pmtmrdnofmv->SetTitle(ss_gr_frac_pmtmrdnofmv.str().c_str());
  gr_frac_fmv->SetTitle(ss_gr_frac_fmv.str().c_str());
  gr_frac_mrd->SetTitle(ss_gr_frac_mrd.str().c_str());
  gr_align_pmt_mrd->SetTitle(ss_gr_align_pmt_mrd.str().c_str());
  gr_align_pmt_fmv->SetTitle(ss_gr_align_pmt_fmv.str().c_str());
  gr_align_mrd_fmv->SetTitle(ss_gr_align_mrd_fmv.str().c_str());
  gr_frac_prompt->SetTitle(ss_gr_frac_prompt.str().c_str());
  gr_frac_delayed->SetTitle(ss_gr_frac_delayed.str().c_str());
  gr_frac_beam->SetTitle(ss_gr_frac_beam.str().c_str());
  gr_frac_led->SetTitle(ss_gr_frac_led.str().c_str());
  gr_frac_cosmic->SetTitle(ss_gr_frac_cosmic.str().c_str());
  gr_avg_mult->SetTitle(ss_gr_avg_mult.str().c_str());
  gr_avg_mult_coinc->SetTitle(ss_gr_avg_mult_coinc.str().c_str());
  gr_avg_mult_coinc_nofmv_cb->SetTitle(ss_gr_avg_mult_coinc_nofmv_cb.str().c_str());
  gr_trigger_beam->SetTitle(ss_gr_trigger_beam.str().c_str());
  gr_trigger_cosmic->SetTitle(ss_gr_trigger_cosmic.str().c_str());
  gr_trigger_led->SetTitle(ss_gr_trigger_led.str().c_str());
  gr_trigger_extended->SetTitle(ss_gr_trigger_extended.str().c_str());
  gr_charge_prompt->SetTitle(ss_gr_charge_prompt.str().c_str());
  gr_charge_del->SetTitle(ss_gr_charge_del.str().c_str());
  gr_charge_per_pmt->SetTitle(ss_gr_charge_per_pmt.str().c_str());
  gr_pot->SetTitle(ss_gr_pot.str().c_str());
  gr_pot_scaled->SetTitle(ss_gr_pot_scaled.str().c_str());
  
file_output->cd();
  TDirectory *dir_gr = file_output->mkdir("graphs");
  dir_gr->cd();
  gr_rate->Write("gr_rate");
  gr_rate_pmt->Write("gr_rate_pmt");
  gr_rate_pmtmrd->Write("gr_rate_pmtmrd");
  gr_rate_pmtmrdfmv->Write("gr_rate_pmtmrdfmv");
  gr_rate_pmtmrdnofmv->Write("gr_rate_pmtmrdnofmv");
  gr_rate_fmv->Write("gr_rate_fmv");
  gr_rate_mrd->Write("gr_rate_mrd");
  gr_frac->Write("gr_frac");
  gr_frac_pmt->Write("gr_frac_pmt");
  gr_frac_pmtmrd->Write("gr_frac_pmtmrd");
  gr_frac_pmtmrdfmv->Write("gr_frac_pmtmrdfmv");
  gr_frac_pmtmrdnofmv->Write("gr_frac_pmtmrdnofmv");
  gr_frac_fmv->Write("gr_frac_fmv");
  gr_align_pmt_mrd->Write("gr_align_pmt_mrd");
  gr_align_pmt_fmv->Write("gr_align_pmt_fmv");
  gr_align_mrd_fmv->Write("gr_align_mrd_fmv");
  gr_frac_prompt->Write("gr_frac_prompt");
  gr_frac_delayed->Write("gr_frac_delayed");
  gr_frac_beam->Write("gr_frac_beam");
  gr_frac_led->Write("gr_frac_led");
  gr_frac_cosmic->Write("gr_frac_cosmic");
  gr_avg_mult->Write("gr_avg_mult");
  gr_avg_mult_coinc->Write("gr_avg_mult_coinc");
  gr_avg_mult_coinc_nofmv_cb->Write("gr_avg_mult_coinc_nofmv_cb");
  gr_trigger_beam->Write("gr_trigger_beam");
  gr_trigger_cosmic->Write("gr_trigger_cosmic");
  gr_trigger_led->Write("gr_trigger_led");
  gr_trigger_extended->Write("gr_trigger_extended");
  gr_charge_prompt->Write("gr_charge_prompt");
  gr_charge_del->Write("gr_charge_del");
  gr_charge_per_pmt->Write("gr_charge_per_pmt");
  gr_pot->Write("gr_pot");
  gr_pot_scaled->Write("gr_pot_scaled");
  file_output->cd();
  TDirectory *dir_histograms = file_output->mkdir("histograms");
  dir_histograms->cd();
  PMT_t_clusters_combined->Write();
  PMT_t_clusters_2pe_combined->Write();
  PMT_t_clusters_cosmics_combined->Write();
  PMT_t_clusters_2pe_cosmics_combined->Write();
  PMT_t_clusters_led_combined->Write();
  PMT_t_clusters_2pe_led_combined->Write();
  MRD_t_clusters_combined->Write();
  MRD_t_clusters_beam_combined->Write();
  MRD_t_clusters_cosmic_combined->Write();
  FMV_t_clusters_combined->Write();
  FMV_t_clusters_beam_combined->Write();
  FMV_t_clusters_cosmic_combined->Write();
  MRD_PMT_t_100pe_Beam_combined->Write();
  FMV_PMT_t_100pe_Beam_combined->Write();
  MRD_PMT_Deltat_100pe_Beam_combined->Write();
  FMV_PMT_Deltat_100pe_Beam_combined->Write();
  PMT_DelayedMult_combined->Write();
  PMT_DelayedMult_Coinc_combined->Write();
  PMT_DelayedMult_Coinc_NoFMV_CB_combined->Write();
  Triggerwords_combined->Write();
  ADCWaveform_Samples_combined->Write();
  PMT_DelayedTime_combined->Write();
  PMT_DelayedTime_CB_combined->Write();
  PMT_DelayedTime_Coinc_combined->Write();
  PMT_DelayedTime_Coinc_NoFMV_combined->Write();
  PMT_DelayedTime_Coinc_NoFMV_CB_combined->Write();
  PMT_prompt_charge_combined->Write();
  PMT_delayed_charge_zoom_combined->Write();
  PMT_chargeperpmt_combined->Write();
  file_output->cd();
  TDirectory *dir_reference = file_output->mkdir("reference");
  dir_reference->cd();
  PMT_t_clusters_2pe_reference->Write();
  MRD_t_clusters_beam_reference->Write();
  FMV_t_clusters_beam_reference->Write();
  ADCWaveform_Samples_reference->Write();
  Triggerwords_reference->Write();
  FMV_PMT_Deltat_100pe_Beam_reference->Write();
  //MRD_PMT_Deltat_100pe_Beam_reference->Write();
  PMT_prompt_charge_reference->Write();
  PMT_delayed_charge_zoom_reference->Write();
  PMT_chargeperpmt_reference->Write();
  file_output->cd();
  TDirectory *dir_canvas = file_output->mkdir("canvas");
  dir_canvas->cd();
  canvas_rate->Write();
  canvas_frac->Write();
  canvas_types->Write();
  canvas_mult->Write();
  canvas_trigger->Write();
  canvas_pot->Write();
  canvas_align->Write();
  canvas_charge->Write();
  canvas_ref_timing->Write();
  canvas_ref_charge->Write();
  canvas_ref_trigger->Write();
  file_output->Close();
  delete file_output;

  return true;
}


bool RunValidationStability::Finalise(){

  return true;
}

void RunValidationStability::ConvertDateToMSec(std::string time_str, uint64_t &time_long){

  std::string epoch_start = "1970/1/1";
  boost::posix_time::ptime Epoch(boost::gregorian::from_string(epoch_start));
  boost::posix_time::ptime ptime_t(boost::posix_time::time_from_string(time_str));
  boost::posix_time::time_duration t_duration;
  t_duration = boost::posix_time::time_duration(ptime_t - Epoch);
  time_long = t_duration.total_milliseconds();//+TimeZoneShift;

}

void RunValidationStability::ConvertWeekToRuns(uint64_t timestamp_start, int &run_start, int &run_end){

   run_start = -1;
   run_end = -1;
   std::map<int,std::map<std::string,std::string>> RunInfoDB;
   bool got_runinfo_db = m_data->CStore.Get("RunInfoDB",RunInfoDB);
   uint64_t timestamp_end = timestamp_start + 7*24*60*60*1000;
   std::vector<int> vec_runs;
    if (got_runinfo_db){
      for (std::map<int,std::map<std::string,std::string>>::iterator it = RunInfoDB.begin(); it != RunInfoDB.end(); it++){
        int temp_run = it->first;
        std::map<std::string,std::string> temp_map = it->second;
        std::string temp_start = temp_map["StartTime"];
        std::string temp_end = temp_map["EndTime"];
        uint64_t time_start, time_end;
        this->ConvertDateToMSec(temp_start,time_start);
        this->ConvertDateToMSec(temp_end,time_end);
        if (time_start > timestamp_start && time_start < timestamp_end) vec_runs.push_back(temp_run);
      }
    } else {
      Log("RunValidationStability tool: Did not find RunInfoDB in CStore - Did you run the LoadRunInfo tool beforehand?",0,verbosity);
    }

  std::sort(vec_runs.begin(),vec_runs.end());

  for (int i_run = 0; i_run < (int) vec_runs.size(); i_run++){
    if (i_run == 0) run_start = vec_runs.at(i_run);
    else if (i_run == (vec_runs.size()-1)) run_end = vec_runs.at(i_run);
  } 

}
