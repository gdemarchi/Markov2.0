function om_decode_timegen_trSNDteSND_alajulienne(subJ,cfg_in)
%% om_decode_timegen_trSNDteSND_alajulienne.m
%%
% This script implements two complementary decoding approaches to investigate prediction-related patterns 
% in the brain, as described more in details in Schubert et al. (2023) and 
% Topalidis et al. (2025) (https://doi.org/10.1016/j.cub.2025.03.064).
%
% In summary:
% - We train a multiclass LDA classifier on tones from the random sequence, excluding tone repetitions (bottom-up classifier).
% - We test this classifier on tone repetitions from both random and predictable sequences to detect pre-activation of stimulus-specific activity.
% - we also relabel tone repetitions in the test set to the most likely expected next tone, capturing prediction-related neural activity.
% - Separately, we train a second classifier (top-down classifier) on forward tone transitions from the predictable sequence (e.g., A/B, B/C).
% - This top-down classifier is tested on random and predictable most probable tone repetitions to identify top-down representations of expected transitions.


%clear all; close all; restoredefaultpath; matlabrc;
disp('Those are the inputs:')
cfg_in


OLDCFG = cfg_in;
% Add user's version of FieldTrip or obob_ownft (update this path as needed)
addpath('/path/to/your/fieldtrip'); % or obob_ownft
ft_defaults;

%% MVPALight
% Add user's version of MVPA-Light (update this path as needed)
addpath(genpath('/path/to/your/MVPA-Light/'));

%%% the rest
% Add user's analysis code directories (update these paths as needed)
addpath('/path/to/your/omissionMarkov/decoding');
addpath('/path/to/your/omissionMarkov/decoding/functions/');

%% out things
% Set user's data and output directories (update these paths as needed)
fileDir = '/path/to/your/data/sssoriginal/';
outDirTGSND = '/path/to/your/data/decoding/matlab/alajulienne/';

%%
clear tmpdata data % to stay on the safe side
conds={'random*','midminus*','midplus*','ordered*'};
trialinfos = [];

%%
for iFile=1:length(conds) %only on orderer if 4

  tmpFile= dir([fileDir,'*',subJ,'_block*',conds{iFile}]);
  cur_file = [tmpFile.folder,'/',tmpFile.name];

  cfg = [];
  cfg.dataset = cur_file;
  cfg.trialdef.triallength = Inf; % infinite trial length for one trial
  cfg.trialdef.ntrials = 1; % one trial
  cfg = ft_definetrial(cfg);

  cfg.channel = cfg_in.chanType; % select channel type
  cfg.hpfilter = 'yes'; % high-pass filter
  cfg.hpfreq = 0.1;
  cfg.hpinstabilityfix = 'split';
  tmpdata = ft_preprocessing(cfg);

  cfg = [];
  cfg.channel = cfg_in.chanType; 
  cfg.dataset=cur_file ;
  cfg.trialdef.prestim =  1;
  cfg.trialdef.poststim = 1;
  cfg.trialdef.eventtype = 'Trigger';
  cfg.trialdef.eventvalue = [1 2 3 4 10 20 30 40]; % only sounds needed
  cfg_wtrials = ft_definetrial(cfg);
  data{iFile}= ft_redefinetrial(cfg_wtrials, tmpdata);
  clear tmpdata;

  data{iFile}.trialinfo(:,2)= iFile;
  %%%
  % Low-pass filter at 30 Hz 
  cfg=[];
  cfg.lpfilter = 'yes'; %usually is yes here
  cfg.lpfreq = 30;
  data{iFile} = ft_preprocessing(cfg,data{iFile});

  %% fix the stimulus delay
  cfg = []; 
  cfg.offset = -24; % 24 samples at 1kHz (approx 23.5 ms delay), in the original study with the original pneumatic tubes
  data{iFile}        = ft_redefinetrial(cfg,data{iFile});

  % Downsample if needed (if fsample > 100 (files on Zenode have Fs = 100))
  if data{iFile}.fsample > 100
     cfg=[];
     if iscell(cfg_in.Fs)
       cfg.resamplefs=str2num(cfg_in.Fs{1}); 
     else
       cfg.resamplefs=str2num(cfg_in.Fs); %
     end
     data{iFile}=ft_resampledata(cfg, data{iFile});
    else end

  trialinfos=[trialinfos; data{iFile}.trialinfo];

  if ~strcmp(cfg_in.chanType,'MEGMAG') % combine planar if not MEGMAG
    cfg=[];
    data{iFile}=ft_combineplanar(cfg, data{iFile});
  else end

end
%% Remove empty blocks if it happens and append data
if max(size(data(~cellfun('isempty',data))))>1 % do not count the empty elements
    cfg = [];
    cfg.appenddim = 'rpt';
    data=ft_appenddata(cfg, data{:});
else
    foo = data(~cellfun('isempty',data)); %take out the only non empty one
    data = foo{1};
    clear foo;
end

%% Remove cfg field to save disk space
data = rmfield(data,'cfg');


%% MVPA part

clear acc* result* cfg*
%% init the magic ....
startup_MVPA_Light;

%% find the indices of the different sounds and omissions
allIdxRD = find(data.trialinfo(:,2)==1); %rd sounds and omissions
allIdxMM = find(data.trialinfo(:,2)==2); %mm sounds and omissions
allIdxMP = find(data.trialinfo(:,2)==3); %mp sounds and omissions
allIdxOR = find(data.trialinfo(:,2)==4); %or sounds and omissions

%% find the self repetitions, see the function below
tmpidx = findSelfRepetitions(data.trialinfo(allIdxRD,1));
selfRepRD = allIdxRD(tmpidx);

% same for the other entropies
tmpidx = findSelfRepetitions(data.trialinfo(allIdxMM,1));
selfRepMM = allIdxMM(tmpidx);
tmpidx = findSelfRepetitions(data.trialinfo(allIdxMP,1));
selfRepMP = allIdxMP(tmpidx);
tmpidx  = findSelfRepetitions(data.trialinfo(allIdxOR,1));
selfRepOR = allIdxOR(tmpidx);

% find the omissions and the trials immediately after
% RD
[tmpidx, tmpidy] = detectOmissions(data.trialinfo(allIdxRD,1));
omIdxRD = allIdxRD(tmpidx);  postOmIdxRD = allIdxRD(tmpidy);

% MM
[tmpidx, tmpidy] = detectOmissions(data.trialinfo(allIdxMM,1));
omIdxMM = allIdxMM(tmpidx);  postOmIdxMM = allIdxMM(tmpidy);

% MP
[tmpidx, tmpidy] = detectOmissions(data.trialinfo(allIdxMP,1));
omIdxMP = allIdxMP(tmpidx);  postOmIdxMP = allIdxMP(tmpidy);

% OR
[tmpidx, tmpidy] = detectOmissions(data.trialinfo(allIdxOR,1));
omIdxOR = allIdxOR(tmpidx);  postOmIdxOR = allIdxOR(tmpidy);


%% create the 'classic' training set, i.e. wo  omissions, and post omissions

% RD
normalIdxRD = setdiff(allIdxRD, union(omIdxRD, postOmIdxRD));
cfg=[];
cfg.keeptrials='yes';
cfg.trials = [normalIdxRD]; 
data_tl_snd_rd=ft_timelockanalysis(cfg, data);
time=data_tl_snd_rd.time;

% MM
normalIdxMM = setdiff(allIdxMM, union(omIdxMM, postOmIdxMM));
cfg=[];
cfg.keeptrials='yes';
cfg.trials = [normalIdxMM]; %
data_tl_snd_mm=ft_timelockanalysis(cfg, data);

% MP
normalIdxMP = setdiff(allIdxMP, union(omIdxMP, postOmIdxMP));
cfg=[];
cfg.keeptrials='yes';
cfg.trials = [normalIdxMP]; %
data_tl_snd_mp=ft_timelockanalysis(cfg, data);

% OR
normalIdxOR = setdiff(allIdxOR, union(omIdxOR, postOmIdxOR));
cfg=[];
cfg.keeptrials='yes';
cfg.trials = [normalIdxOR]; %
data_tl_snd_or=ft_timelockanalysis(cfg, data);


%% create the 'forward' training set, i.e. the one without self repetitions and omissions, and post omissions
% RD
forwardIdxRD = setdiff(allIdxRD, union(selfRepRD, union(omIdxRD, postOmIdxRD)));
cfg=[];
cfg.keeptrials='yes';
cfg.trials = [forwardIdxRD]; %
data_tl_rdFW=ft_timelockanalysis(cfg, data); 

% MM
forwardIdxMM = setdiff(allIdxMM, union(selfRepMM, union(omIdxMM, postOmIdxMM)));
cfg=[];
cfg.keeptrials='yes';
cfg.trials = [forwardIdxMM]; %
data_tl_mmFW=ft_timelockanalysis(cfg, data);

% MP
forwardIdxMP = setdiff(allIdxMP, union(selfRepMP, union(omIdxMP, postOmIdxMP)));
cfg=[];
cfg.keeptrials='yes';
cfg.trials = [forwardIdxMP]; %
data_tl_mpFW=ft_timelockanalysis(cfg, data);

% OR
forwardIdxOR = setdiff(allIdxOR, union(selfRepOR, union(omIdxOR, postOmIdxOR)));
cfg=[];
cfg.keeptrials='yes';
cfg.trials = [forwardIdxOR]; %
data_tl_orFW=ft_timelockanalysis(cfg, data);

%% create the self repetitions test set, i.e. the one with only self repetitions
% RD
SRidxRD = setdiff(selfRepRD, union(omIdxRD, postOmIdxRD)); %remove OMs as well
cfg=[];
cfg.keeptrials='yes';
cfg.trials =  SRidxRD;%
data_tl_rdSR=ft_timelockanalysis(cfg, data);

% MM
SRidxMM = setdiff(selfRepMM, union(omIdxMM, postOmIdxMM)); %remove OMs as well
cfg=[];
cfg.keeptrials='yes';
cfg.trials = [SRidxMM]; %
data_tl_mmSR=ft_timelockanalysis(cfg, data);

% MP
SRidxMP = setdiff(selfRepMP, union(omIdxMP, postOmIdxMP)); %remove OMs as well
cfg=[];
cfg.keeptrials='yes';
cfg.trials = [SRidxMP]; %
data_tl_mpSR=ft_timelockanalysis(cfg, data);

% OR
SRidxOR = setdiff(selfRepOR, union(omIdxOR, postOmIdxOR)); %remove OMs as well
cfg=[];
cfg.keeptrials='yes';
cfg.trials = [SRidxOR]; %
data_tl_orSR=ft_timelockanalysis(cfg, data);

%% 'classical' way a la Demarchi et al. 2019
% RD train/test
cfg =  [];
cfg.classifier = 'multiclass_lda';
cfg.metric     = 'accuracy';
cfg.preprocessing = 'undersample';
[accTG_SND_RD, result_accTG_SND_RD] = mv_classify_timextime(cfg, data_tl_snd_rd.trial, data_tl_snd_rd.trialinfo(:,1));

% Cross decoding
% Rd_SND to Mm_SND
[accTG_SND_RD_MM, result_accTG_SND_RD_MM] = mv_classify_timextime(cfg, data_tl_snd_rd.trial, data_tl_snd_rd.trialinfo(:,1),data_tl_snd_mm.trial, data_tl_snd_mm.trialinfo(:,1));

% Rd_SND to Mp_SND
[accTG_SND_RD_MP, result_accTG_SND_RD_MP] = mv_classify_timextime(cfg,data_tl_snd_rd.trial, data_tl_snd_rd.trialinfo(:,1),data_tl_snd_mp.trial, data_tl_snd_mp.trialinfo(:,1));

% Rd_SND to Or_SND 
[accTG_SND_RD_OR, result_accTG_SND_RD_OR] = mv_classify_timextime(cfg,data_tl_snd_rd.trial, data_tl_snd_rd.trialinfo(:,1),data_tl_snd_or.trial, data_tl_snd_or.trialinfo(:,1));


%%  train on ordered forward and test on self repetitions
cfg =  [];
cfg.classifier = 'multiclass_lda';
cfg.metric     = 'accuracy';
cfg.preprocessing = 'undersample';
[accTG_SND_fwOR_srRD, result_acc_fwOR_srRD] = mv_classify_timextime(cfg, data_tl_orFW.trial, data_tl_orFW.trialinfo(:,1), data_tl_rdSR.trial, data_tl_rdSR.trialinfo(:,1));
[accTG_SND_fwOR_srMM, result_acc_fwOR_srMM] = mv_classify_timextime(cfg, data_tl_orFW.trial, data_tl_orFW.trialinfo(:,1), data_tl_mmSR.trial, data_tl_mmSR.trialinfo(:,1));
[accTG_SND_fwOR_srMP, result_acc_fwOR_srMP] = mv_classify_timextime(cfg, data_tl_orFW.trial, data_tl_orFW.trialinfo(:,1), data_tl_mpSR.trial, data_tl_mpSR.trialinfo(:,1));
[accTG_SND_fwOR_srOR, result_acc_fwOR_srOR] = mv_classify_timextime(cfg, data_tl_orFW.trial, data_tl_orFW.trialinfo(:,1), data_tl_orSR.trial, data_tl_orSR.trialinfo(:,1));

% the same on the relabeled / most probable
cfg =  [];
cfg.classifier = 'multiclass_lda';
cfg.metric     = 'accuracy';
cfg.preprocessing = 'undersample';
[accTG_SND_fwOR_srRDmp, result_acc_fwOR_srRDmp] = mv_classify_timextime(cfg, data_tl_orFW.trial, data_tl_orFW.trialinfo(:,1), data_tl_rdSR.trial, makeMostProbable(data_tl_rdSR.trialinfo(:,1)));
[accTG_SND_fwOR_srMMmp, result_acc_fwOR_srMMmp] = mv_classify_timextime(cfg, data_tl_orFW.trial, data_tl_orFW.trialinfo(:,1), data_tl_mmSR.trial, makeMostProbable(data_tl_mmSR.trialinfo(:,1)));
[accTG_SND_fwOR_srMPmp, result_acc_fwOR_srMPmp] = mv_classify_timextime(cfg, data_tl_orFW.trial, data_tl_orFW.trialinfo(:,1), data_tl_mpSR.trial, makeMostProbable(data_tl_mpSR.trialinfo(:,1)));
[accTG_SND_fwOR_srORmp, result_acc_fwOR_srORmp] = mv_classify_timextime(cfg, data_tl_orFW.trial, data_tl_orFW.trialinfo(:,1), data_tl_orSR.trial, makeMostProbable(data_tl_orSR.trialinfo(:,1)));

%%  training on the random forward
cfg =  [];
cfg.classifier = 'multiclass_lda';
cfg.metric     = 'accuracy'% {'accuracy', 'confusion', 'f1', 'mae'}; %  
cfg.output_type =  'dval'; %
cfg.preprocessing = 'undersample';
[accTG_SND_fwRD_srRD, result_acc_fwRD_srRD] = mv_classify_timextime(cfg, data_tl_rdFW.trial, data_tl_rdFW.trialinfo(:,1), data_tl_rdSR.trial, data_tl_rdSR.trialinfo(:,1));
[accTG_SND_fwRD_srMM, result_acc_fwRD_srMM] = mv_classify_timextime(cfg, data_tl_rdFW.trial, data_tl_rdFW.trialinfo(:,1), data_tl_mmSR.trial, data_tl_mmSR.trialinfo(:,1));
[accTG_SND_fwRD_srMP, result_acc_fwRD_srMP] = mv_classify_timextime(cfg, data_tl_rdFW.trial, data_tl_rdFW.trialinfo(:,1), data_tl_mpSR.trial, data_tl_mpSR.trialinfo(:,1));
[accTG_SND_fwRD_srOR, result_acc_fwRD_srOR] = mv_classify_timextime(cfg, data_tl_rdFW.trial, data_tl_rdFW.trialinfo(:,1), data_tl_orSR.trial, data_tl_orSR.trialinfo(:,1));

% most probable
cfg =  [];
cfg.classifier = 'multiclass_lda';
cfg.metric     = 'accuracy'% {'accuracy', 'confusion', 'f1', 'mae'}; %  
cfg.output_type =  'dval'; %
cfg.preprocessing = 'undersample';
[accTG_SND_fwRD_srRDmp, result_acc_fwRD_srRDmp] = mv_classify_timextime(cfg, data_tl_rdFW.trial, data_tl_rdFW.trialinfo(:,1), data_tl_rdSR.trial, makeMostProbable(data_tl_rdSR.trialinfo(:,1)));
[accTG_SND_fwRD_srMMmp, result_acc_fwRD_srMMmp] = mv_classify_timextime(cfg, data_tl_rdFW.trial, data_tl_rdFW.trialinfo(:,1), data_tl_mmSR.trial, makeMostProbable(data_tl_mmSR.trialinfo(:,1)));
[accTG_SND_fwRD_srMPmp, result_acc_fwRD_srMPmp] = mv_classify_timextime(cfg, data_tl_rdFW.trial, data_tl_rdFW.trialinfo(:,1), data_tl_mpSR.trial, makeMostProbable(data_tl_mpSR.trialinfo(:,1)));
[accTG_SND_fwRD_srORmp, result_acc_fwRD_srORmp] = mv_classify_timextime(cfg, data_tl_rdFW.trial, data_tl_rdFW.trialinfo(:,1), data_tl_orSR.trial, makeMostProbable(data_tl_orSR.trialinfo(:,1)));



%% and save!

cfg_in = OLDCFG; % to keep the original cfg_in

if iscell(cfg_in.Fs)
  Fs = cfg_in.Fs{1};
else
  Fs = cfg_in.Fs;
end

outFile = [ subJ '_'  cfg_in.chanType '_timegen_trRdSND_teSND_Fs' num2str(Fs) '_alajulienne.mat' ];

if iscell(outFile)
  outFile = [outFile{:}];
else end

save(fullfile(outDirTGSND, outFile),'acc*','result*','time','cfg_in' ,'-v7.3'); 
disp(['Saved: ',fullfile(outDirTGSND, outFile)]);

end


%% local functions
function repetition_indices = findSelfRepetitions(sequence)
  % Initialize an empty array to store indices of the first self-repetitions
  repetition_indices = [];
  % Flag to track if we're in a repetition sequence
  inRepetition = false;
  % Loop through the sequence starting from the second element
  for i = 2:length(sequence)
      % Check if the current element is equal to the previous element
      if sequence(i) == sequence(i-1)
          % If not already in a repetition, this is the first repetition
          if ~inRepetition
              repetition_indices = [repetition_indices, i];
              inRepetition = true;  % Mark that we're now in a repetition
          end
      else
          % If the current element is different, reset the repetition flag
          inRepetition = false;
      end
  end
end


function [omIdx, postOmIdx] = detectOmissions(sequence)
  % Find indices where elements are 10 or more (detect omissions)
  omIdx = find(sequence >= 10); 
  
  % Initialize array to store the indices following the omissions
  postOmIdx = omIdx + 1; 
  
  % Ensure indices are within bounds
  postOmIdx(postOmIdx > length(sequence)) = []; 
  
end

function shifted_sequence = makeMostProbable(sequence)
  % mod  to perform the cyclic shift
  % 2 -> 1, 3 -> 2, 4 -> 3, 1 -> 4, ...
  shifted_sequence = mod(sequence - 2, 4) + 1;
end


