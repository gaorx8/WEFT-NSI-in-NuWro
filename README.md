Event-level CCQE detection NSI in NuWro and impact on oscillation-parameter extraction

This repository provides the implementation used in our paper to study CCQE detection NSI effects within the WEFT framework.

The workflow combines:
1.	NuWro (fixed commit)
2.	Oscillated atmospheric flux (nuSQUIDS + HKKM)
3.	WEFT cross-section tables
4.	ROOT-based figure generation

All figures in the paper can be reproduced following the steps below.

1. Clone repository
 
2. Get NuWro
  git clone https://github.com/NuWro/nuwro.git
  cd nuwro
  git checkout d67a4d9a9360fc642deb548328c28e23efc5bffb

3. Copy the modified files from this repository into the NuWro directory.
- Files located under `nuwro/src/` must be copied into your `nuwro/src/` directory
- Files located outside `src/` in this repository should be placed at the corresponding top-level location in the NuWro directory

4. Make NuWro
  make clean && make

Figure 1:
  python catch_sigma.py
  python plot_ratio.py

Figure 2:
  nuwro -i params_SM.txt -o SM.root
  nuwro -i params_SM_R.txt -o SM_R.root
  myroot -b -q 'plotcostheta.c()'

Figure 3:
  - move to the atn folder
  bash run_fig3.sh

Figure 4:
  bash run_fig4.sh

Figure 5: (The input text file is constructed by collecting Δχ²_bias = χ²(truth) − χ²_min from the outputs of run_fig4.sh)
  python plot_epsilon_chi2.py


