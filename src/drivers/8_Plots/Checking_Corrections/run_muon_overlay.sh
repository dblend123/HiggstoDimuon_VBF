#!/usr/bin/env bash

usage() {
  echo "Usage:"
  echo "  $0 -u <user> -y <year> -c <corrType> -n <configName>"
  echo
  echo "Notes:"
  echo "  <year> may be a full year (e.g. 2024) or a segmented era"
  echo "  (e.g. 12023, 22023). This value is used directly as the"
  echo "  Analyzed_SkimmedFiles subdirectory."
  echo
  echo "Directory layout assumed:"
  echo "  /eos/user/<first-letter>/<user>/Analyzed_SkimmedFiles/<year>/<corrType>/<configName>/"
  echo
  echo "Examples:"
  echo "  # 2024 data (no segmentation)"
  echo "  $0 -u memmad -y 2024  -c Rochester -n noScaleSmear_noBtagSF"
  echo
  echo "  # 2023 data, first segment"
  echo "  $0 -u memmad -y 12023 -c Rochester -n noScaleSmear_noBtagSF"
  echo
  echo "  # 2023 data, second segment with additional systematics off"
  echo "  $0 -u memmad -y 22023 -c Rochester -n noScaleSmear_noJERC_noJetVeto_noBtagSF"
  echo
  echo "  # 2023 data, stricter jet criteria"
  echo "  $0 -u memmad -y 22023 -c Rochester -n noScaleSmear_noJERC_noJetVeto_noBtagSF_noJetCriteria"
  exit 1
}
while getopts "u:y:e:c:n:" opt; do
  case $opt in
    u) USER=$OPTARG ;;
    y) YEAR=$OPTARG ;;
    e) ERA=$OPTARG ;;
    c) CORR=$OPTARG ;;
    n) CONFIG=$OPTARG ;;
    *) usage ;;
  esac
done

if [[ -z "$USER" || -z "$YEAR" || -z "$CORR" || -z "$CONFIG" ]]; then
  usage
fi

echo "[INFO] Running Muon_Overlay_AllPlots"
echo "  user        = $USER"
echo "  year        = $YEAR"
echo "  corrType    = $CORR"
echo "  configName  = $CONFIG"

root -l -b -q \
  "Muon_Overlay_AllPlots.C(\"$USER\",\"$YEAR\",\"$CORR\",\"$CONFIG\")"

