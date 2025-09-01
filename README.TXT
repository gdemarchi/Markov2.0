Markov 2.0 a.k.a Code for the 
"Nature Communications: Matters Arising" doi:123456
===============================================

Overview
--------
This repository contains MATLAB and R scripts for the analyses in Demarchi et al. (2025), including reanalysis of the original dataset 
and a pilot experiment with bias-free ordered sequences. The code supports time-generalized decoding and top-down vs. bottom-up classification.

Contents
--------
- `alajulienne/`: Scripts for the main reanalysis and hybrid Schubert/Demarchi decoding framework.
    - `om_decode_timegen_trSNDteSND_alajulienne.m`: Main decoding script using time-generalization and hybrid method.
    - `plot_TG_trSNDteSND_MVPALight_mr_alajulienne.m`: Plots time-generalization matrices and statistical results for the hybrid analysis.
    - `run_om_decode_timegen_trSNDteSND_alajulienne.m`: SLURM job submission script for batch decoding on a cluster.
- `4and8tones/`: Scripts for the pilot experiment with the 4/8 tones and tone sequence analysis.
    - `mtng_decode_timegen_trSNDteSND.m`: Decoding script for pilot data using standard MVPA as in Demarchi et al. 2019.
    - `plot_mtng_SNDtoSND.m`: Plots time-generalization results for the pilot experiment.
    - `run_mtng_decode_timegen_trSNDteSND.m`: SLURM job submission script for pilot decoding.
    - `plot_stat_sequences.R`: R script for statistical analysis of sequence regression data.


Data
----
Raw MEG data are available on Zenodo (https://doi.org/10.5281/zenodo.3268713). 
Pilot data available upon request only now, will be fully shared once published.

Prerequisites
-------------
- MATLAB R2020a or later
- FieldTrip toolbox
- MVPA-Light toolbox
- Sufficient disk space for MEG data (~20 GB)

Citing this work
----------------
If you use this code or data, please cite:
Demarchi, G., et al. (2025). Markov 2.0. Nature Communications Matters Arising xxxyyy

License
-------
MIT License

Contact
-------
For questions, comments or suggestions, please contact `gianpaolo.demarchi@plus.ac.at` (or open an issue on this repository, but YMMV).