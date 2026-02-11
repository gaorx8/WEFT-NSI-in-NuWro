#!/usr/bin/env bash
set -euo pipefail
 
# ASIMOV_NPOT=139000
ASIMOV_NPOT=1000000

echo "Start in $(pwd)..."

pushd atn_ana/flux_model_xsec > /dev/null
echo "Now in atn_ana/flux_model_xsec..."
echo ">>> 1. Multiply the model by xsec..."
root -l -q 'model_xsec_L.cxx()'
root -l -q 'model_xsec_R.cxx()'
root -l -b -q 'plot_event_rate_diff.cxx()'
popd > /dev/null