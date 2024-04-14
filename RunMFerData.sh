#!/bin/sh

while read -u 10 RUN
do
  echo "Working on Run $RUN ..."
  filelist="/annie/app/users/jhe/MyToolAnalysis_MFer/configfiles/MuonFitter/my_inputs_R${RUN}.txt"

  sed -i "2s#.*#FileForListOfInputs ./configfiles/MuonFitter/R${RUN}_single_infile.txt#" configfiles/MuonFitter/Data/LoadANNIEEventConfig
  sed -i "18s#.*#AiEtaFile ev_ai_eta_R${RUN}.txt#" configfiles/MuonFitter/Data/MuonFitterConfig

  let i=0
  while read -r file 
  do
    echo "Processing $file"
    echo "$file" > configfiles/MuonFitter/R${RUN}_single_infile.txt

    sed -i "3s#.*#OutputFile BeamRun_R${RUN}_${i}.ntuple.root#" configfiles/MuonFitter/Data/PhaseIITreeMakerConfig
    sed -i "3s#.*#OutputFile MuonFitter_Eta_R${RUN}_${i}.root#" configfiles/MuonFitter/Data/MuonFitterConfig
    sed -i "23s#.*#OutputFile R${RUN}_EvDisplay_${i}#" configfiles/MuonFitter/Data/EventDisplayConfig

    ./Analyse MuonFitter
          #mv BeamRun_R${RUN}_${i}.ntuple.root /pnfs/annie/scratch/users/jhe/ntuples/data/
          mv MuonFitter_Eta_R${RUN}_${i}.root /pnfs/annie/scratch/users/jhe/data/
          #mv R${RUN}_EvDisplay_${i}.root /pnfs/annie/scratch/users/jhe/data/
          i=$((i+1))
  done < $filelist 

  #if [ ! -d "c_r${RUN}_new" ]; then 
  #  echo "c_r${RUN}_new does not exist. Creating one now.."
  #  mkdir c_r${RUN}_new
  #fi

  #mv c_eta_ai_r${RUN}*.png c_r${RUN}_new

done 10<runlist.txt

