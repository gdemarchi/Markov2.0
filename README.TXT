Markov 2.0 a.k.a Code for the 
"Nature Communications: Matters Arising" doi:123456
===============================================

Overview
--------
This repository contains the code and documentation needed to understand the analyses reported in our follow‑up to Demarchi et al. (2019) (doi:xxxx). 
The goal is to reanalyse the original dataset, address the criticisms raised by Abdoun et al. and demonstrate a new classification framework that distinguishes top‑down and bottom‑up predictive patterns. 
It also includes a pilot experiment that implements a bias‑free ordered sequence design.

Contents
--------
The repository is organised as follows:

- `reshuffled/`: MATLAB scripts for recreating the original analysis pipeline and running the decoding analyses on the reshuffled datasets. For the reshuffling it Abdoun et al. Python code has been used. Please refer to their repository for getting the necessary Python scripts.
- `alajulienne/`: MATLAB scripts implementing the top‑down vs. bottom‑up classification approach described in the paper, along with helper functions for time‑generalised decoding.
- `4and8tones/`: Code for the pilot experiment, including the MATLAB scripts to preprocess and analyse the MEG data.

Data
----
The original dataset is openly available on Zenodo (version used in Demarchi et al., 2019). Because of its size (≈20 GB), the raw MEG data are not included directly in this repository. 
To run the analyses you will need to download the data from Zenodo following the link below and extract it into a suitable` directory.

```
https://doi.org/10.5281/zenodo.3268713
```

The pilot experiment is available (by now) only upon reasonable request. We plan to release it within a follow up preregisterd study.
It comprises MEG recordings from 27 participants collected during a seminar course in the summer semester of 2025 using a bias‑free ordered tone design.

Prerequisites
-------------
- MATLAB R2020a or later.
- Sufficient disk space to store the MEG data (≈20 GB for the original dataset). About the same amount if you want to try the analysis on the reshuffled data.



Pilot experiment
----------------
The pilot study used an ordered sequence with a unique randomised frequency order for each participant to eliminate bias in the decoding due to neighbouring frequencies. 
Here you can find the code reproduce the pilot analysis, it matches 1 to 1 the analysis of Demarchi et al., 2019. 

Citing this work
----------------
If you use this code or data for your own research, please cite the following reference:

```
Demarchi, G., Hartmann, T., Schubert, J., Topalidis, P., et al. (2025). Markov 2.0. Nature Communications Matters Arising xxxyyy

```

Please also cite the Zenodo dataset with its DOI if you reuse the data.

License
-------
This repository is licensed under the MIT License.

Contact
-------
For questions, comments or suggestions, please contact `gianpaolo.demarchi@plus.ac.at` (or open an issue on this repository, but YMMV).