%% Main script: plot time generalisation training RdSND testing SND with statistics

% Clear workspace and set up FieldTrip
clear all; close all; restoredefaultpath; matlabrc;
addpath('/path/to/fieldtrip');
ft_defaults;

%%

% Add MVPA-Light toolbox and set up paths
addpath(genpath('/path/to/MVPA-Light/'));
startup_MVPA_Light;

% Optional: add export_fig if needed
%addpath('/path/to/export_fig/');
% or where smooth2a.m resides
%addpath('/path/to/mattols/');

% Set input/output directories
fileDir = '/path/to/decoding_output/';
outDir = '/path/to/figures/';

%% List of subjects
subjList = {'19920916slsa', ...
     '19940511krln', ...
     '20010127urdr', ...
     '20010618ttbl', ...
     '20010629tnli', ...
     '20011116mrsr', ...
     '20020521bidv', ...
     '20020827ktsh', ...
     '20021005asap', ...
     '20021230urre', ...
     '20030208crsh', ...
     '20030225aehf', ...
     '20030328ptbn', ...
     '20030502rnhn', ...
     '20030516ktba', ...
     '20030717begl', ...
     '20030812argn', ...
     '20030829mroc', ...
     '20030901adjg', ...
     '20030919mrgs', ...
     '20031221mlht', ...
     '20040210gidn', ...
     '20040308pgwl', ...
     '20040330brhr', ...
     '20040416hlln', ...  
     '20031013aglv', ...  
     '20011126sbgi', ... 
};

%%

% Load decoding results for each subject
clear fake* TG_acc* data* comp* tmptg* fileToRead curFile MIA
Fs = 100;
numS = length(subjList);
oldMVPA = 'yes';

clear TG_acc* gaTF* MIA
fNamePart = ['*Fs' num2str(Fs) 'Hz_acc_oldMVPA' oldMVPA '.mat'];
for iFile=1:numS
  if ~isempty(dir([fileDir subjList{iFile} fNamePart ]))
    fileToRead  = dir([fileDir subjList{iFile} fNamePart]);
  else     
    fprintf('\n Missing file for subject %d ...', iFile)
    MIA{iFile} = subjList{iFile};
    continue
  end

  curFile = fileToRead.name;
  tmptg=load([fileDir curFile]);

  if ndims(tmptg.accTG_SND_RD8) > 2 % not accuracy, e.g. confusion matrix
    TG_confusion_RD(iFile,:,:,:,:) = tmptg.accTG_SND_RD;
  else % accuracy
    TG_acc_RD4(iFile,:,:) = tmptg.accTG_SND_RD4;
    TG_acc_RD8(iFile,:,:) = tmptg.accTG_SND_RD8;
    TG_acc_RD_OR4(iFile,:,:) = tmptg.accTG_SND_RD_OR4;
    TG_acc_RD_OR8(iFile,:,:) = tmptg.accTG_SND_RD_OR8;
  end

  trainTime{iFile} = tmptg.time;
  iFile;
  sprintf('\nDone %s,  %d %%',subjList{iFile}, round(100*iFile/length(subjList)));
end
tTime = tmptg.time;


% Remove elements with rows full of zeros (if the case)
if ndims(tmptg.accTG_SND_RD4) ==  2 % accuracy
  % Compute grand averages
  gaTG_acc_RD4 = squeeze(mean(TG_acc_RD4));
  gaTG_acc_RD8 = squeeze(mean(TG_acc_RD8));
  gaTG_acc_RD_OR4 = squeeze(mean(TG_acc_RD_OR4));
  gaTG_acc_RD_OR8 = squeeze(mean(TG_acc_RD_OR8));
end



% Plot grand averages/diagonals for visual inspection
figure;
hold all;
for iFile=1:numS
    plot(tTime, diag(squeeze(TG_acc_RD4(iFile,:,:))));
    plot(tTime, diag(squeeze(TG_acc_RD_OR4(iFile,:,:))));

end
plot(tTime, diag(gaTG_acc_RD4), 'r','LineWidth',4);
plot(tTime, diag(gaTG_acc_RD_OR4),'g', 'LineWidth',4);
title('diag RD4 (red) vs OR_RD4 (green) ', 'Interpreter','none');
xlim([-.8 .8]);
ylim([0.2 .45])


%% Check for missing subjects and update subject count
badList = [];
if exist('MIA')
  badList=MIA(~cellfun('isempty',MIA)) ;
else end

% set numSubject to the real number, without the badList
numS = numS - length(badList)
% NOTE: badList should be zero ...

%%  Smooth accuracy matrices for each subject (for stats)
clear TG_acc_smooth*;
nInterp = 5;
for iFile=1:numS
  TG_acc_smoothRD_OR8(iFile,:,:) = smooth2a(squeeze(TG_acc_RD_OR8(iFile,:,:)),nInterp); %
  TG_acc_smoothRD_OR4(iFile,:,:) = smooth2a(squeeze(TG_acc_RD_OR4(iFile,:,:)),nInterp); %
  TG_acc_smoothRD4(iFile,:,:) = smooth2a(squeeze(TG_acc_RD4(iFile,:,:)),nInterp); %
  TG_acc_smoothRD8(iFile,:,:) = smooth2a(squeeze(TG_acc_RD8(iFile,:,:)),nInterp); % 
end


%%  Plot time generalisation matrices (TG plots)
ft_hastoolbox('brewermap', 1);
     
figure;
subplot(2,2,1)
imagesc([tTime(1) tTime(end)],[tTime(1) tTime(end)],smooth2a(gaTG_acc_RD_OR8 ,5));
%caxis([.23 .27]);
colorbar;
xlim([-.4 .8]);
ylim([0 .4]);
barH= colorbar;
ylabel(barH, 'accuracy');
set(gca,'YDir','normal');
subtitle('RD to OR 8 Tones');

subplot(2,2,2)
imagesc([tTime(1) tTime(end)],[tTime(1) tTime(end)],smooth2a(gaTG_acc_RD_OR4,5));
%caxis([.23 .27]);
colorbar;
xlim([-.4 .8]);
ylim([0 .4]);;
barH= colorbar;
ylabel(barH, 'accuracy');
set(gca,'YDir','normal');
subtitle('RD to OR 4 Tones');

subplot(2,2,3)
imagesc([tTime(1) tTime(end)],[tTime(1) tTime(end)],smooth2a(gaTG_acc_RD8,5));
%caxis([.23 .27]);
colorbar;
xlim([-.4 .8]);
ylim([0 .4]);
barH= colorbar;
ylabel(barH, 'accuracy');
set(gca,'YDir','normal');
subtitle('RD to RD 8 Tones');

subplot(2,2,4)
imagesc([tTime(1) tTime(end)],[tTime(1) tTime(end)],smooth2a(gaTG_acc_RD4,5));
%caxis([.23 .27]);
colorbar;
xlim([-.4 .8]);
ylim([0 .4]);
barH= colorbar;
ylabel(barH, 'accuracy');
set(gca,'YDir','normal');
subtitle('RD to RD 4 Tones');
%titleStr=['trRdSNDteSND - single entropies - ' num2str(numP) '  Patients'];
colormap(brewermap(256, '*RdYlBu'));


%% Plot difference matrices to show 'direction' of statistical effects
diffr = (gaTG_acc_RD_OR4-gaTG_acc_RD4);
xtrm = max(abs(min(min(diffr))), abs(max(max(diffr))));
figure;
subplot(1,2,1)
imagesc([tTime(1) tTime(end)],[tTime(1) tTime(end)],smooth2a((diffr),5));
%caxis([.23 .27]);
colorbar;
xlim([-.4 .8]);
ylim([0 .4]);
clim([-xtrm/2 xtrm/2]);
barH= colorbar;
ylabel(barH, 'accuracy');
set(gca,'YDir','normal');
subtitle('RD to OR 4 minus RD to RD - 4 Tones');


diffr = (gaTG_acc_RD_OR8-gaTG_acc_RD8);
xtrm = max(abs(min(min(diffr))), abs(max(max(diffr))));
subplot(1,2,2)
imagesc([tTime(1) tTime(end)],[tTime(1) tTime(end)],smooth2a((diffr),5));
%caxis([.23 .27]);
colorbar;
xlim([-.4 .8]);
ylim([0 .4]);
clim([-xtrm/2 xtrm/2]);
barH= colorbar;
ylabel(barH, 'accuracy');
set(gca,'YDir','normal');
subtitle('RD to OR minus RD to RD - 8 Tones');


% % Prepare data for statistics (FieldTrip format)
clear tmpMat fake* stat*;
for iFile=1:numS%
  %or8
  tmpMat = TG_acc_smoothRD_OR8(iFile,:,:);
  
  fakeTF_RD_OR8{iFile} = [];
  fakeTF_RD_OR8{iFile}.fsample = Fs;
  fakeTF_RD_OR8{iFile}.powspctrm = tmpMat;
  fakeTF_RD_OR8{iFile}.label = {'accuracy'};
  fakeTF_RD_OR8{iFile}.time = tTime;
  fakeTF_RD_OR8{iFile}.freq= tTime;
  fakeTF_RD_OR8{iFile}.dimord = 'freq_time';
 
  %or4
  tmpMat = TG_acc_smoothRD_OR4(iFile,:,:); 
  
  fakeTF_RD_OR4{iFile} = [];
  fakeTF_RD_OR4{iFile}.fsample = Fs;
  fakeTF_RD_OR4{iFile}.powspctrm = tmpMat;
  fakeTF_RD_OR4{iFile}.label = {'accuracy'};
  fakeTF_RD_OR4{iFile}.time = tTime;
  fakeTF_RD_OR4{iFile}.freq= tTime;
  fakeTF_RD_OR4{iFile}.dimord = 'freq_time'; 
  
  %rd8
  tmpMat = TG_acc_smoothRD8(iFile,:,:);
  
  fakeTF_RD8{iFile} = [];
  fakeTF_RD8{iFile}.fsample = Fs;
  fakeTF_RD8{iFile}.powspctrm = tmpMat;
  fakeTF_RD8{iFile}.label = {'accuracy'};
  fakeTF_RD8{iFile}.time = tTime;
  fakeTF_RD8{iFile}.freq= tTime;
  fakeTF_RD8{iFile}.dimord = 'freq_time'; 
  
  
  %rd4
  tmpMat = TG_acc_smoothRD4(iFile,:,:);
  
  fakeTF_RD4{iFile} = [];
  fakeTF_RD4{iFile}.fsample = Fs;
  fakeTF_RD4{iFile}.powspctrm = tmpMat;
  fakeTF_RD4{iFile}.label = {'accuracy'};
  fakeTF_RD4{iFile}.time = tTime;
  fakeTF_RD4{iFile}.freq= tTime;
  fakeTF_RD4{iFile}.dimord = 'freq_time'; 
  
end


%% Run cluster-based permutation statistics (FieldTrip)
nRand = 1000; % only 1000 to speed up
cfg4 = [];
cfg4.channel          = {'accuracy'};
cfg4.latency          = [-.333 0];
cfg4.frequency        = [0 .333]; 
cfg4.method           = 'montecarlo'; % 
cfg4.statistic        = 'ft_statfun_depsamplesT'; %
cfg4.tail             = 0;
cfg4.clustertail      = 0;
cfg4.clusteralpha     = 0.05;     
cfg4.alpha       = 0.05;
cfg4.correctm    = 'cluster';
cfg4.numrandomization = nRand;

% experimental design
clear design
design(1,:) = [1:numS 1:numS];
design(2,:) = [ones(1,numS) 2*ones(1,numS)];
cfg4.design           = design;
cfg4.uvar             = 1;
cfg4.ivar             = 2;

[stat4] = ft_freqstatistics(cfg4,fakeTF_RD_OR4{:},fakeTF_RD4{:});

% Save statistics
save([outDir, 'stat4'], 'stat4');

% Plot cluster statistics for 4 tones
cfg = [];
cfg.parameter      = 'stat';
cfg.maskparameter  = 'mask';
cfg.maskstyle      = 'outline';
cfg.highlight      = 'on';
ft_singleplotTFR(cfg,stat4)

% adapt the image
colormap(brewermap(256, '*RdYlBu'));
title(['stat RD to OR - 4  tones -  smoothN=' num2str(nInterp) ''  ] , 'FontSize', 16 )
%refline(1,.333);
xlim([stat4.time(1) stat4.time(end)]);
ylim([stat4.freq(1) stat4.freq(end)]);
barH= colorbar;
ylabel(barH, 't value');
cmax = max(abs([max(max(stat4.stat)), min(min(stat4.stat))]));
clim([-cmax cmax])
title(['' '' ],'Interpreter', 'none');
xlabel('Testing time /s');
ylabel('Training time /s');
title(['stat RD to OR - 4 tones - smoothN=' num2str(nInterp) ''  ] , 'FontSize', 16 )
saveas(gcf, [outDir,'SNDtoSND_4tones'], 'fig');
saveas(gcf, [outDir,'SNDtoSND_4tones'], 'svg');
saveas(gcf, [outDir,'SNDtoSND_4tones'],'pdf');



%% Repeat statistics and plotting for 8 tones
nRand = 1000;
cfg8 = [];
cfg8.channel          = {'accuracy'};
cfg8.latency          = [-.333 0];
cfg8.frequency        = [0 .333]; 
cfg8.method           = 'montecarlo';
cfg8.statistic        = 'ft_statfun_depsamplesT'; 
cfg8.tail             = 0;
cfg8.clustertail      = 0;
cfg8.clusteralpha     = 0.05;      
cfg8.correctm    = 'cluster';
cfg8.numrandomization = nRand;
clear design
design(1,:) = [1:numS 1:numS];
design(2,:) = [ones(1,numS) 2*ones(1,numS)];

cfg8.design           = design;
cfg8.uvar             = 1;
cfg8.ivar             = 2;

[stat8] = ft_freqstatistics(cfg8,fakeTF_RD_OR8{:},fakeTF_RD8{:});
save([outDir, 'stat8'], 'stat8');


% Plot cluster statistics for 8 tones
cfg = [];
cfg.parameter      = 'stat';
cfg.maskparameter  = 'mask';
cfg.maskstyle      = 'outline'; % 'opacity', 'saturation' 'outline' 
cfg.highlight      = 'on';
ft_singleplotTFR(cfg,stat8);

colormap(brewermap(256, '*RdYlBu'));
%refline(1,.333);
xlim([stat8.time(1) stat8.time(end)]);
ylim([stat8.freq(1) stat8.freq(end)]);
barH= colorbar;
ylabel(barH, 't value');
cmax = max(abs([max(max(stat8.stat)), min(min(stat8.stat))]));
clim([-cmax cmax])
title(['' '' ],'Interpreter', 'none');
xlabel('Testing time /s');
ylabel('Training time /s');
title(['stat RD to OR - 8  tones -  smoothN=' num2str(nInterp) ''  ] , 'FontSize', 16 )
saveas(gcf, [outDir,'SNDtoSND_8tones'], 'fig');
saveas(gcf, [outDir,'SNDtoSND_8tones'], 'svg');
saveas(gcf, [outDir,'SNDtoSND_8tones'],'pdf');


% ends here