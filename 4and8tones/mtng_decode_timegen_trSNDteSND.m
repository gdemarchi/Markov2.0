function mtng_decode_timegen_trSNDteSND(subJ,cfg_in)
%% training and testing on sound

% to test, e.g.
% subJ = '20030502rnhn'
% cfg_in = [];


%%% fieldtrip
addpath('/path/to/fieldtrip');
ft_defaults;

%%% the rest if needed
addpath('/path/to/functions');

%%
if ~isfield(cfg_in,'chanType') cfg_in.chanType = 'MEGMAG'; else end
if ~isfield(cfg_in,'pre') cfg_in.pre = '1'; else end
if ~isfield(cfg_in,'post') cfg_in.post = '1'; else end
if ~isfield(cfg_in,'metric') cfg_in.metric = 'acc'; else end
if ~isfield(cfg_in,'oldmvpa') cfg_in.oldmvpa = 'yes'; else end
if ~isfield(cfg_in,'Fs') cfg_in.Fs = '100'; else end

disp(cfg_in); %show the input
%% reading the data part, always the same acrosss different types of
%%% decoding later: we do every time also to save data ...

baseDir = '/path/to/cleaned_data/';
outDir = '/path/to/decoding_output/';

%%% MVPALight
if strcmp(cfg_in.oldmvpa ,'yes') %load the raw data and apply ica later
    addpath(genpath('/path/to/oldMVPA-Light/'));
else
    addpath(genpath('/path/to/MVPA-Light/'));
end

%% Data reading
clear tmpdata data* % to stay on the safe side
trialinfos = [];


% Find all files in baseDir starting with subJ and ending with .fif
filePattern = fullfile(baseDir, [subJ, '*.fif']);
allFiles = dir(filePattern);
allFiles = allFiles(~contains({allFiles.name}, '-1')); % Exclude files containing '-1'
fprintf('Found %d files for subject %s\n', numel(allFiles), subJ);


for iFile = 1:numel(allFiles)

    fprintf('Processing file %d of %d: %s\n', iFile, numel(allFiles), allFiles(iFile).name);
    cur_file = fullfile(allFiles(iFile).folder, allFiles(iFile).name);

    % read the events and filter them
    [cond_code, trigger_codes] = file2cond(cur_file) %double check this!!


    % high pass on the continous
    cfg = [];
    cfg.dataset = cur_file;
    cfg.trialdef.triallength = Inf;
    cfg.trialdef.ntrials = 1;
    cfg = ft_definetrial(cfg);

    cfg.channel = cfg_in.chanType; %do only on magnetometers usuallt, since maxfiltered
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 0.1;
    cfg.hpinstabilityfix =  'split';
    tmpdata = ft_preprocessing(cfg);

    try
        cfg = [];
        cfg.channel = cfg_in.chanType; % mag, grads, all together - see above
        cfg.dataset = cur_file;
        cfg.trialdef.prestim = str2num(cfg_in.pre);
        cfg.trialdef.poststim = str2num(cfg_in.post);
        cfg.trialdef.eventtype = 'Trigger';
        cfg.trialdef.eventvalue = [trigger_codes]; % keep oms for later
        cfg_wtrials = ft_definetrial(cfg);
    catch ME
        fprintf('No good trials found in the file %s - %s!\n', cur_file, ME.message);
        fprintf('Moving on with the next block!\n');
        continue; % Skip to the next iteration of the loop
    end

    % moving on with the block!
    data{iFile}= ft_redefinetrial(cfg_wtrials, tmpdata);
    clear tmpdata;

    data{iFile}.trialinfo(:,2)= cond_code;

    %%% fix the stimulus delay
    cfg = [];
    cfg.offset = -16; %16 samples at 1KHz with the new tubes set
    data{iFile}        = ft_redefinetrial(cfg,data{iFile});

    % downsamp always, otherwise too big
    cfg=[];
    cfg.resamplefs=str2num(cfg_in.Fs); %change on cluster
    data{iFile}=ft_resampledata(cfg, data{iFile});

    % code the info
    trialinfos=[trialinfos; data{iFile}.trialinfo];

    % combine planar if grads or both
    if (strcmp(cfg_in.chanType,'MEGGRAD') || strcmp(cfg_in.chanType,'MEG'))
        cfg=[];
        data{iFile}=ft_combineplanar(cfg, data{iFile});
    else end

end

%% getting rid of empty blocks, if there ... shouldnt 
if max(size(data(~cellfun('isempty',data))))>1 % do not count the empty elements
    cfg = [];
    cfg.appenddim = 'rpt';
    data=ft_appenddata(cfg, data{:});
else
    foo = data(~cellfun('isempty',data)); %take out the only non empty one
    data = foo{1};
    clear foo;
end

%% a bit of useless cleaning, saves a lot of disk!
data = rmfield(data,'cfg');

% check tha in the data.trialinfo there is the all these codes are mapped and nothing is missnig
all_codes = [1:4, 65:72, 129:132, 193:200];
missing_codes = setdiff(all_codes, unique(data.trialinfo(:,1)));
if ~isempty(missing_codes)
    error('Missing trial codes in data.trialinfo: %s', num2str(missing_codes));
end


%%%%%%%%%%%%%%%%%%%%%%%% DECODING PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  time generalisation - train on RD sound test on R/OR sounds

clear acc* result*
%% init the magic ....
startup_MVPA_Light;

%% train on random 4/8 sounds, test on random and ordered 4/8 sounds

tonesRDidx4 = find(data.trialinfo(:,2)==1);
tonesORidx4 = find(data.trialinfo(:,2)==2);
tonesRDidx8 = find(data.trialinfo(:,2)==3);
tonesORidx8 = find(data.trialinfo(:,2)==4);

%% select the trials for the timelock analysis
cfg=[];
cfg.keeptrials='yes';
cfg.trials =[tonesRDidx4];
data_tl_rd4=ft_timelockanalysis(cfg, data);
cfg.trials =[tonesORidx4];
data_tl_or4=ft_timelockanalysis(cfg, data);
cfg.trials =[tonesRDidx8];
data_tl_rd8=ft_timelockanalysis(cfg, data);
cfg.trials =[tonesORidx8];
data_tl_or8=ft_timelockanalysis(cfg, data);



%% RDtoRD
cfg =  [];
cfg.classifier = 'multiclass_lda';
cfg.metric     = cfg_in.metric;
%cfg.CV = 'none'; %no cross validation, for testing if it runs

if strcmp(cfg_in.oldmvpa ,'yes') %as in Demarchi 2019
    cfg.balance = 'undersample';
else
    cfg.preprocess      = 'undersample';
end

[accTG_SND_RD4, result_accTG_SND_RD4] = mv_classify_timextime(cfg, data_tl_rd4.trial, classRelabel(data_tl_rd4.trialinfo(:,1)));
[accTG_SND_RD8, result_accTG_SND_D8] = mv_classify_timextime(cfg, data_tl_rd8.trial, classRelabel(data_tl_rd8.trialinfo(:,1)));

%% RDtoOR
cfg =  [];
cfg.classifier = 'multiclass_lda';
cfg.metric     = cfg_in.metric;

if strcmp(cfg_in.oldmvpa ,'yes') %load the raw data and apply ica later
    cfg.balance = 'undersample';
else
    cfg.preprocess      = 'undersample';
end

% HERE: it is cross decdoding so i have to relabel the trialinfo, i.e.
% 129 130 131 132 shoule become 1 2 3 4
% 193 194 195 196 197 198 199 200 should become 65 66 67 68 69 70 71 72

data_tl_or4.trialinfo(:,1) = data_tl_or4.trialinfo(:,1) - 128; % 129-128=1, 130-128=2, etc
data_tl_or8.trialinfo(:,1) = data_tl_or8.trialinfo(:,1) - 128; % 193-128=65, 194-128=66, etc

[accTG_SND_RD_OR4, result_accTG_SND_RD_OR4] = mv_classify_timextime(cfg, data_tl_rd4.trial, classRelabel(data_tl_rd4.trialinfo(:,1)), data_tl_or4.trial, classRelabel(data_tl_or4.trialinfo(:,1)));
[accTG_SND_RD_OR8, result_accTG_SND_RD_OR8] = mv_classify_timextime(cfg, data_tl_rd8.trial, classRelabel(data_tl_rd8.trialinfo(:,1)), data_tl_or8.trial, classRelabel(data_tl_or8.trialinfo(:,1)));

% save one time axis
time=data_tl_rd8.time;


%% and save!
if iscell(cfg_in.Fs)
    cfg_in.Fs = cell2mat(cfg_in.Fs);
else end

outFile = [ subJ '_markov_tng_Fs' cfg_in.Fs 'Hz_' cfg_in.metric '_oldMVPA' cfg_in.oldmvpa];

save (fullfile(outDir, outFile),'acc*','result*','time*','cfg*' ,'stuff*','-v7.3'); %weights !!!
disp('Data saved!')



%% ancillary functions
function [cond_code, trigger_codes] = file2cond(fName)

    % check the content of the events, and assign it to the right spot
    % 1 <- 4 tones, random  codes: 1, 2, 3, 4
    % 2 <- 4 tones, ordered codes: 129, 130, 131, 132
    % 3 <- 8 tones, random  codes: 65, 66, 67, 68, 69, 70, 71, 72
    % 4 <- 8 tones, ordered codes: 193, 194, 195, 196, 197, 198, 199, 200
    
    events = ft_read_event(fName, 'chanidx', find(strcmp(ft_read_header(fName).label, 'STI101'))); % read the events from the file, STI101 is the channel for triggers
    evtsnum = [events.value];
    unique_events = unique(evtsnum);
    
    % filter some out
    unique_events(unique_events == 5) = [];
    unique_events(unique_events == 64) = [];
    unique_events(unique_events == 128) = [];
    unique_events(unique_events == (128+64)) = [];

    
    % if there are events outside the range of 
    % 5, 64, 128, then 1:4, 129:132, 65:72, 193:200
    % emit an error and exit *not* gracefully ...
    
    % Define allowed event codes
    allowed_codes = [1:4, 5, 64, 65:72, 128, 129:132, 193:200];

    % Check for unexpected event codes
    if any(~ismember(unique_events, allowed_codes))
        error('Unexpected event codes found in the file: %s', fName);
    end

    if all(ismember(unique_events, 1:4))
        cond_code = 1;
    elseif all(ismember(unique_events, 129:132))
        cond_code = 2;
    elseif all(ismember(unique_events, 65:72))
        cond_code = 3;
    elseif all(ismember(unique_events, 193:200))
        cond_code = 4;
    else
        error('Unrecognized trigger pattern');
    end

    switch cond_code
        case 1
            cond_string = '4 random tones';
        case 2
            cond_string = '4 ordered tones';
        case 3
            cond_string = '8 random tones';
        case 4
            cond_string = '8 ordered tones';
        otherwise
            error('Unexpected condition code: %d', cond_code);
    end

    trigger_codes = unique_events(unique_events ~= 0 & unique_events ~= 5 & unique_events ~= 64 & unique_events ~= 128 & unique_events ~= (128+64));
    fprintf('Condition code: %d, i.e. %s - Trigger codes: %s\n', cond_code, cond_string, num2str(trigger_codes));

end


function clabel1toN = classRelabel(clabel)
% relabels a list of conditions with any number of conditions to the
% standard 1:nclasses that mvpa_light likes

classesLabel = unique(clabel);
nclasses = max(size(classesLabel));

%init output
clabel1toN=zeros(length(clabel),1);

for iCount=1:nclasses
    clabel1toN(clabel == classesLabel(iCount))= iCount;
end

