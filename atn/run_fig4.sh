#!/usr/bin/env bash
set -euo pipefail
 
# ASIMOV_NPOT=139000
ASIMOV_NPOT=1000000

echo "Start in $(pwd)..."

pushd atn_ana/flux_model_xsec > /dev/null
echo "Now in atn_ana/flux_model_xsec..."
echo ">>> 1. Multiply the model by xsec..."
root -l -q 'model_xsec.cxx()'
popd > /dev/null

pushd atn_ana/flux_model_xsec_res > /dev/null
echo "Now in atn_ana/flux_model_xsec_res..."
echo ">>> 2. Build resolusion model..."
root -l -b -q 'build_resolution_model_mu.cxx'
root -l -b -q 'build_resolution_model_mu_bar.cxx'

echo ">>> 3. Smear..."
root -l -b -q 'smear_model.cxx("mu","")'
root -l -b -q 'smear_model.cxx("mu_bar","")'

echo ">>> 4. Combine mu and mu_bar..."
hadd -ff atnu_mu_all_no_model_xsec_res.root atnu_mu_bar_model_xsec_res.root atnu_mu_model_xsec_res.root
popd > /dev/null

pushd atn_ana > /dev/null
echo "Now in atn_ana..."
echo ">>> 5. Generate asimov datasat..."
root -l -b -q "generate_asimov_dataset.cxx(${ASIMOV_NPOT},\"_no\")"

echo ">>> 6. Analyze..."
root -l -b -q 'analysis.cxx(0.01, 0.01, "_no")'

echo ">>> 7. Plot..."
root -l -b -q 'plot_analysis_multiple.cxx("_no")'

echo "Completed!" 
popd > /dev/null  