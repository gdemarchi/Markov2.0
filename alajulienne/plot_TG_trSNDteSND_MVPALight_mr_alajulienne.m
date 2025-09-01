%% Time Generalization Plotting Script
% This script loads decoding results for a set of subjects and plots time generalization matrices and statistics.
% Paths are set to generic placeholders for portability. Only standard FieldTrip is required.

clear; close all; clc;
addpath('/path/to/fieldtrip'); % Add FieldTrip toolbox
ft_defaults;

figDir = '/path/to/figures/'; % Directory to save figures
addpath('/path/to/decoding/functions/'); % Add decoding functions
addpath('/path/to/markov_revenge/matlab/'); % Add analysis scripts
addpath('/path/to/mattols/'); % Add custom tools, e.g smooth2a.m

% List of subjects to process
subjList = {
    'PNRK', 'KRHR', 'GBSH', 'BRHC', 'CRLE', 'ANSR', 'SSLD', 'AGSG', 'RFTM', 'SLBR', 'GDZN', 'EEHB',
    'BTKC', 'GNTA', 'SZDT', 'SBPE', 'KTAD', 'IMSH', 'ATLI', 'HLHY', 'IGSH', 'MCSH', 'CRBC', 'GBHL',
    'MNSU', 'IIQI', 'HIEC', 'KRKE', 'BRSH', 'LLZM', 'EIFI', 'MRGU', 'IONP'
};
nSubj = numel(subjList);


% Initialize variables and set parameters
clear TG* tmp* tg* ga* trainTime
Fs = 100; % Sampling frequency
fileDir = '/path/to/decoding/results/'; % Directory containing decoding results
fNamePart = '*MEGMAG_timegen_trRdSND_teSND_Fs100_alajulienne.mat'; % Filename pattern

% Loop over subjects and load decoding results
for iFile = 1:nSubj
    subjID = subjList{iFile};
    resultFiles = dir([fileDir subjID fNamePart]);
    if isempty(resultFiles)
        fprintf('\n File for subject %s missing ...', subjID);
        continue
    end
    curFile = resultFiles(1).name;
    tmptg = load([fileDir curFile]);

    % Store decoding results for each subject
    TG_acc_RD(iFile,:,:)      = tmptg.accTG_SND_RD;
    TG_acc_RD_MM(iFile,:,:)   = tmptg.accTG_SND_RD_MM;
    TG_acc_RD_MP(iFile,:,:)   = tmptg.accTG_SND_RD_MP;
    TG_acc_RD_OR(iFile,:,:)   = tmptg.accTG_SND_RD_OR;

    TG_acc_SND_fwOR_srRD(iFile,:,:) = tmptg.accTG_SND_fwOR_srRD;
    TG_acc_SND_fwOR_srMM(iFile,:,:) = tmptg.accTG_SND_fwOR_srMM;
    TG_acc_SND_fwOR_srMP(iFile,:,:) = tmptg.accTG_SND_fwOR_srMP;
    TG_acc_SND_fwOR_srOR(iFile,:,:) = tmptg.accTG_SND_fwOR_srOR;

    TG_acc_SND_fwRD_srRD(iFile,:,:) = tmptg.accTG_SND_fwRD_srRD;
    TG_acc_SND_fwRD_srMM(iFile,:,:) = tmptg.accTG_SND_fwRD_srMM;
    TG_acc_SND_fwRD_srMP(iFile,:,:) = tmptg.accTG_SND_fwRD_srMP;
    TG_acc_SND_fwRD_srOR(iFile,:,:) = tmptg.accTG_SND_fwRD_srOR;

    % Most probable decoding results
    TG_acc_SND_fwOR_srRDmp(iFile,:,:) = tmptg.accTG_SND_fwOR_srRDmp;
    TG_acc_SND_fwOR_srMMmp(iFile,:,:) = tmptg.accTG_SND_fwOR_srMMmp;
    TG_acc_SND_fwOR_srMPmp(iFile,:,:) = tmptg.accTG_SND_fwOR_srMPmp;
    TG_acc_SND_fwOR_srORmp(iFile,:,:) = tmptg.accTG_SND_fwOR_srORmp;

    TG_acc_SND_fwRD_srRDmp(iFile,:,:) = tmptg.accTG_SND_fwRD_srRDmp;
    TG_acc_SND_fwRD_srMMmp(iFile,:,:) = tmptg.accTG_SND_fwRD_srMMmp;
    TG_acc_SND_fwRD_srMPmp(iFile,:,:) = tmptg.accTG_SND_fwRD_srMPmp;
    TG_acc_SND_fwRD_srORmp(iFile,:,:) = tmptg.accTG_SND_fwRD_srORmp;

    trainTime{iFile} = tmptg.time;
    fprintf('Processing file %s for subject %s\n', curFile, subjID);
end

% Create the filename title string for plots
pattern = 'MEGMAG_timegen_trRdSND_teSND.*(?=\.mat)';
extracted_part = regexp(curFile, pattern, 'match');
title_str = extracted_part{1}; % Used for plot titles

% Compute grand averages across subjects
% (add comments for each block below as needed)
gaTG_acc_RD      = squeeze(mean(TG_acc_RD));
gaTG_acc_RD_MM   = squeeze(mean(TG_acc_RD_MM));
gaTG_acc_RD_MP   = squeeze(mean(TG_acc_RD_MP));
gaTG_acc_RD_OR   = squeeze(mean(TG_acc_RD_OR));

gaTG_acc_fwOR_srRD = squeeze(mean(TG_acc_SND_fwOR_srRD));
gaTG_acc_fwOR_srMM = squeeze(mean(TG_acc_SND_fwOR_srMM));
gaTG_acc_fwOR_srMP = squeeze(mean(TG_acc_SND_fwOR_srMP));
gaTG_acc_fwOR_srOR = squeeze(mean(TG_acc_SND_fwOR_srOR));

gaTG_acc_fwRD_srRD = squeeze(mean(TG_acc_SND_fwRD_srRD));
gaTG_acc_fwRD_srMM = squeeze(mean(TG_acc_SND_fwRD_srMM));
gaTG_acc_fwRD_srMP = squeeze(mean(TG_acc_SND_fwRD_srMP));
gaTG_acc_fwRD_srOR = squeeze(mean(TG_acc_SND_fwRD_srOR));

% Most probable grand averages
gaTG_acc_fwOR_srRDmp = squeeze(mean(TG_acc_SND_fwOR_srRDmp));
gaTG_acc_fwOR_srMMmp = squeeze(mean(TG_acc_SND_fwOR_srMMmp));
gaTG_acc_fwOR_srMPmp = squeeze(mean(TG_acc_SND_fwOR_srMPmp));
gaTG_acc_fwOR_srORmp = squeeze(mean(TG_acc_SND_fwOR_srORmp));

gaTG_acc_fwRD_srRDmp = squeeze(mean(TG_acc_SND_fwRD_srRDmp));
gaTG_acc_fwRD_srMMmp = squeeze(mean(TG_acc_SND_fwRD_srMMmp));
gaTG_acc_fwRD_srMPmp = squeeze(mean(TG_acc_SND_fwRD_srMPmp));
gaTG_acc_fwRD_srORmp = squeeze(mean(TG_acc_SND_fwRD_srORmp));



%% now plot the TG - cross decoding - one by one manually first
ft_hastoolbox('brewermap', 1);

imagesc([trainTime{1}(1) trainTime{1}(end)],[trainTime{1}(1) trainTime{1}(end)],gaTG_acc_fwRD_srRDmp);
%caxis([.20 .30]);
%caxis([-.05 .05]);
colorbar;
refline(1,.333);
refline(1,-.333);
refline(1);
xlim([-.7 .7]);
ylim([-.33 .33]);
clim([.2 .3]);
barH= colorbar;
ylabel(barH, 'accuracy');
set(gca,'YDir','normal');
title('gaTG_acc_fwOR_srORmp','Interpreter', 'none');


%% stat part, 4 entropies

nInterp = 5;

% fwRD
cond1= TG_acc_SND_fwRD_srRDmp;
cond2= TG_acc_SND_fwRD_srMMmp;
cond3= TG_acc_SND_fwRD_srMPmp;
cond4= TG_acc_SND_fwRD_srORmp;

% fwOR later
% cond1= TG_acc_SND_fwOR_srRDmp;
% cond2= TG_acc_SND_fwOR_srMMmp;
% cond3= TG_acc_SND_fwOR_srMPmp;
% cond4= TG_acc_SND_fwOR_srORmp;

for iFile=1:length(subjList)
    smoothC1(iFile,:,:) = smooth2a(squeeze(cond1(iFile,:,:)),nInterp); %
    smoothC2(iFile,:,:) = smooth2a(squeeze(cond2(iFile,:,:)),nInterp); %
    smoothC3(iFile,:,:) = smooth2a(squeeze(cond3(iFile,:,:)),nInterp); %
    smoothC4(iFile,:,:) = smooth2a(squeeze(cond4(iFile,:,:)),nInterp); %
end

% convert to pseudo TF - for fieldtrip 
clear tmpMat
tmpMat(:,1,:,:) = smoothC1;
fakeTF_cond1 = [];
fakeTF_cond1.fsample = Fs;
fakeTF_cond1.powspctrm = tmpMat;
fakeTF_cond1.label = {'accuracy'};
fakeTF_cond1.time = trainTime{1};
fakeTF_cond1.freq= trainTime{1};
fakeTF_cond1.dimord = 'rpt_chan_freq_time'; 

tmpMat(:,1,:,:) = smoothC2;
fakeTF_cond2 = [];
fakeTF_cond2.fsample = Fs;
fakeTF_cond2.powspctrm = tmpMat;
fakeTF_cond2.label = {'accuracy'};
fakeTF_cond2.time = trainTime{1};
fakeTF_cond2.freq= trainTime{1};
fakeTF_cond2.dimord = 'rpt_chan_freq_time'; 

tmpMat(:,1,:,:) = smoothC3;
fakeTF_cond3 = [];
fakeTF_cond3.fsample = Fs;
fakeTF_cond3.powspctrm = tmpMat;
fakeTF_cond3.label = {'accuracy'};
fakeTF_cond3.time = trainTime{1};
fakeTF_cond3.freq= trainTime{1};
fakeTF_cond3.dimord = 'rpt_chan_freq_time'; 


tmpMat(:,1,:,:) = smoothC4;
fakeTF_cond4 = [];
fakeTF_cond4.fsample = Fs;
fakeTF_cond4.powspctrm = tmpMat;
fakeTF_cond4.label = {'accuracy'};
fakeTF_cond4.time = trainTime{1};
fakeTF_cond4.freq= trainTime{1};
fakeTF_cond4.dimord = 'rpt_chan_freq_time';

% 4 entropies regression

% fake ordered
nRand = 1000;
cfg = [];
cfg.channel          = {'accuracy'};
cfg.latency          = [-0.3 .3]; % test time
cfg.frequency        = [-.33 0];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesregrT';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha       = 0.05;
cfg.correctm    = 'cluster';
cfg.numrandomization = nRand;
%cfg.computecritval = 'yes';

design = zeros(2,4*length(subjList));%zeros(2,4*length(subjList));
design(1,1:length(subjList)) = 1;
design(1,length(subjList)+1:2*length(subjList)) = 2;
design(1,2*length(subjList)+1:3*length(subjList)) = 3;
design(1,3*length(subjList)+1:4*length(subjList)) = 4;
design(2,:) = [1:length(subjList), 1:length(subjList), 1:length(subjList), 1:length(subjList)];

cfg.design           = design;
cfg.ivar             = 1;
cfg.uvar             = 2;

[stat] = ft_freqstatistics(cfg,fakeTF_cond1,fakeTF_cond2,fakeTF_cond3,fakeTF_cond4);

% quick check
cfg = [];
cfg.parameter      = 'stat';
cfg.maskparameter  = 'mask';
cfg.maskstyle      = 'outline';
cfg.highlight  = 'on';
ft_singleplotTFR(cfg,stat)
refline(1,.333);
refline(1,-.333);
refline(1);
xlim([stat.time(1) stat.time(end)]);
ylim([stat.freq(1) stat.freq(end)]);
barH= colorbar;
ylabel(barH, 't value');
cmax = max(abs([max(max(stat.stat)), min(min(stat.stat))]));
clim([-cmax cmax])
title(['' '' ],'Interpreter', 'none');
xlabel('Testing time /s');
ylabel('Training time /s');
title('Regression trfwRD -> tested most probable(RD-MM-MP-OR)')