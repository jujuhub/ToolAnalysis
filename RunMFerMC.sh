#!/bin/sh

filelist="/annie/app/users/jhe/MyToolAnalysis_MFer/wcsim_filelist_0-999.txt"

let i=0
while read -r file 
do
	echo "$file"
	sed -i "9s#.*#InputFile ${file}#" configfiles/MuonFitter/MC/LoadWCSimConfig
	sed -i "3s#.*#OutputFile MC_BeamRun_${i}.ntuple.root#" configfiles/MuonFitter/MC/PhaseIITreeMakerConfig
	sed -i "3s#.*#OutputFile MC_MuonFitter_${i}.root#" configfiles/MuonFitter/MC/MuonFitterConfig
	sed -i "23s#.*#OutputFile MFMC_EvDisplay_${i}#" configfiles/MuonFitter/MC/EventDisplayConfig
	#sed -i "21s#.*#OutputFile MC_LoadWCSim_${i}.root#" configfiles/MuonFitter/MC/LoadWCSimConfig
	./Analyse MuonFitter
        mv MC_MuonFitter_${i}.root /pnfs/annie/scratch/users/jhe/mc/
        #mv MC_BeamRun_${i}.ntuple.root /pnfs/annie/scratch/users/jhe/ntuples/mc/
        #mv MFMC_EvDisplay_${i}.root /pnfs/annie/scratch/users/jhe/mc/
        #mv MC_LoadWCSim_${i}.root /pnfs/annie/scratch/users/jhe/mc/
        i=$((i+1))
done < $filelist 

