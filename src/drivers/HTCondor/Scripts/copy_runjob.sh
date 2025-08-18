#!/bin/bash

export X509_USER_PROXY=x509up.cert

inputList="$1"
tag="$2"

echo "Job started on $(date)"
echo "Input list: $inputList"
echo "Tag: $tag"

echo "INFO: Run ntuplizer..."
# echo "./process_HiggsAnalysis_wCutsExtended_NoPairing_v0 --input $inputList --output output_${tag}.root --isMC false"
# echo "./process_HiggsAnalysis_wCutsExtended_NoPairing_v0 --input $inputList --output output_${tag}.root --isMC false --puDATA pileup.root"
echo "./process_HiggsAnalysis_wCutsExtended_NoPairing_v0 --input $inputList --output output_${tag}.root --isMC true --genPUMC true --puMC DY_MCPileup_${tag}.root"
# echo "./process_HiggsAnalysis_wCutsExtended_NoPairing_v0 --input $inputList --output output_${tag}.root --isMC true --genPUMC false --puMC DY_MCPileup_${tag}.root"

# ./process_HiggsAnalysis_wCutsExtended_NoPairing_v0 --input "$inputList" --output "output_${tag}.root" --isMC false
# ./process_HiggsAnalysis_wCutsExtended_NoPairing_v0 --input "$inputList" --output "output_${tag}.root" --isMC false --puDATA pileup.root
./process_HiggsAnalysis_wCutsExtended_NoPairing_v0 --input "$inputList" --output "output_${tag}.root" --isMC true --genPUMC true --puMC DY_MCPileup_${tag}.root 
# ./process_HiggsAnalysis_wCutsExtended_NoPairing_v0 --input "$inputList" --output "output_${tag}.root" --isMC true --genPUMC false --puMC DY_MCPileup_${tag}.root

echo "Job finished on $(date)"
