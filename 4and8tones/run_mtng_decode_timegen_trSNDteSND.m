%% Example SLURM job submission script for SCC Salzburg Computing Cluster
% This script demonstrates how to run decoding jobs on the SCC cluster using plus-slurm-matlab.
% NOTE: Paths, job settings, and subject lists should be adapted to your own environment.
% plus-slurm-matlab is free and available via GitLab: https://gitlab.com/thht/plus-slurm-matlab
% obob_ft (FieldTrip cluster utilities) is available via GitLab: https://gitlab.com/thht/obob_ft
% YMMV for other clusters or local setups.

restoredefaultpath; matlabrc;
clear all
close all

cfg = [];
cfg.package.plus_slurm = true;
% Add required toolboxes (update paths to your installation)
addpath('/path/to/obob_ft/');
obob_init_ft;

addpath('/path/to/plus-slurm-matlab/');

% Example: define nodes to exclude (optional, SCC specific)
baseNode = 'node%02d.scc-pilot.plus.ac.at';
nodeNumbers = [1, 3, 4, 5, 6, 7, 8, 9, 10];
nodeStrings = arrayfun(@(x) sprintf(baseNode, x), nodeNumbers, 'UniformOutput', false);
nodeStrings = string(nodeStrings);
nodeStrings = arrayfun(@(x) sprintf(baseNode, x), nodeNumbers, 'UniformOutput', false);
nodeStrings = string(nodeStrings);

% SLURM job configuration (adapt memory, CPUs, time, and paths as needed)
cfg = [];
cfg.mem = '32G';  % Memory per job
cfg.cpus = 2;     % Number of CPUs per job
cfg.request_time = 4*60; % Time in minutes
cfg.jobsdir      = '/path/to/jobs/'; % Change to your jobs directory
%cfg.qos          = 'high_prio'; % Uncomment if needed
% matlab container
cfg.container    = 'oras://ghcr.io/thht/obob-singularity-container/xfce_desktop_matlab:latest';
%cfg.exclude_nodes= strjoin(nodeStrings, ',');

% Add your analysis code directory (update path)
addpath('/path/to/matlab/scripts')

% List of subjects to process (adapt to your dataset)
subjList = {'20011126sbgi', '20031013aglv', '20040330brhr'; ...
                  '20040416hlln', '20040424nvls', '20040507rshr'; ...
                  '20040513mnbr', '20040514jmhn', '20040520tphn'; ...
                  '20040527rrhr', '20040603grsr', '20040604mrbr'; ...
                  '20040610rrhr', '20040617sbsr', '20040624rrhr'; ...
                  '20040701lvbr', '20040708rrhr', '20040715rrhr'; ...
                  '20040722rrhr', '20040729rrhr', '20040805rrhr'; ...
                  '20040812rrhr', '20040819rrhr', '20040826rrhr'; ...
                  '20040902rrhr', '20040909rrhr', '20040916rrhr'};

                  subjList = subjList(:); % Convert to column vector
% Analysis configuration (adapt as needed)
cfg = [];
cfg.chanType = {'MEGMAG'}; % Channel type: MEGMAG, MEGGRAD, etc.
cfg.Fs = '1000';             % Sampling frequency
cfg.oldmvpa = 'yes';       % Use original MVPA pipeline, or new one
%cfg.metric = 'confusion'; % Other metrics: acc, dval, auc, etc.
% Add jobs to SLURM struct
slurm_struct = plus_slurm.addjob_cell(slurm_struct, 'mtng_decode_timegen_trSNDteSND', subjList, cfg);

%% Submit jobs to SLURM
plus_slurm.submit(slurm_struct)