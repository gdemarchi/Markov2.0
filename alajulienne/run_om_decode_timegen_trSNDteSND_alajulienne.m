%% Example SLURM job submission script for SCC Salzburg Computing Cluster
% This script demonstrates how to run decoding jobs on the SCC cluster using plus-slurm-matlab.
% NOTE: Paths, job settings, and subject lists should be adapted to your own environment.
% plus-slurm-matlab is free and available via GitLab: https://gitlab.com/thht/plus-slurm-matlab
% obob_ft (FieldTrip cluster utilities) is available via GitLab: https://gitlab.com/thht/obob_ft
% YMMV for other clusters or local setups.

restoredefaultpath; matlabrc;
clear all
close all

% Add required toolboxes (update these paths for your environment)
addpath('/path/to/obob_ownft/');


cfg = [];
cfg.package.plus_slurm = true;
obob_init_ft(cfg);
addpath('/path/to/plus-slurm-matlab/');

% SLURM job configuration (adapt memory, CPUs, time, and paths as needed)
cfg = [];
cfg.mem = '64G';  % Memory per job
cfg.cpus = 4;     % Number of CPUs per job
cfg.request_time = 8*60; % Time in minutes
cfg.jobsdir      = '/path/to/jobs/'; % Change to your jobs directory
%cfg.qos          = 'high_prio'; % Uncomment if needed
cfg.container    = 'oras://ghcr.io/thht/obob-singularity-container/xfce_desktop_matlab:latest';
%cfg.exclude_nodes= 'node01,node02'; % Example: comma separated list

slurm_struct = plus_slurm.create(cfg);

% Add your analysis code directory (update path)
addpath('/path/to/matlab/');

% List of subjects to process (adapt to your dataset)
subjList = {
  'PNRK', 'KRHR', 'GBSH', 'BRHC', 'CRLE', 'ANSR', 'SSLD', 'AGSG', 'RFTM', 'SLBR',
  'GDZN', 'EEHB', 'BTKC', 'GNTA', 'SZDT', 'SBPE', 'KTAD', 'IMSH', 'ATLI', 'HLHY',
  'IGSH', 'MCSH', 'CRBC', 'GBHL', 'MNSU', 'IIQI', 'HIEC', 'KRKE', 'BRSH', 'LLZM',
  'EIFI', 'MRGU', 'IONP'
};


% Analysis configuration (adapt as needed)
cfg = [];
cfg.chanType = {'MEGMAG'}; % Channel type: MEGMAG, MEGGRAD, etc.
cfg.Fs = {'100'};           % Sampling frequency


% Add jobs to SLURM struct
slurm_struct = plus_slurm.addjob_cell(slurm_struct, 'om_decode_timegen_trSNDteSND_alajulienne', subjList, cfg);


%% Submit jobs to SLURM
plus_slurm.submit(slurm_struct)
