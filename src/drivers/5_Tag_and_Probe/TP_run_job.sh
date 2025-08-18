#!/bin/bash

export X509_USER_PROXY=x509up.cert

inputList="$1"
tag="$2"

echo "Job started on $(date)"
echo "Input list: $inputList"
echo "Tag: $tag"

echo "INFO: Run Tag&Probe Script..."
echo "root -l Root_Tag_and_Probe.C+"

root -l -q "Root_Tag_and_Probe.C(\"$1\", \"$2\")"

echo "Job finished on $(date)"
