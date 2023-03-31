%%%% TMS+EEG preprocessing script %%%%
clear;
clc;
mainpath =cd;
filepath=[mainpath '/TMS+EEG_controldata/'];
%filepath =[mainpath '/TMS+EEG_MDDdata/sham/post/'];
%filepath =[mainpath '/TMS+EEG_MDDdata/active/post/'];
sub_no = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',...
    '16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'};

%% parameter settings
for ns= [72] %1:length(sub_no)
    cd([filepath 'sub' sub_no{ns}])
    filename =deblank(ls('*.vhdr'));
    disp(['loading filename: ' filename])
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadbv([filepath,'sub' sub_no{ns} '/'],filename,[],[1:63]);

    % Inputs:
    %   EEG                  - Raw continuous EEG data in EEGLAB data structure (for one condition)
    %   cfg                  - Configuration variable with hyperparameters
    cfg.EventCode           = 'S255';
    cfg.TMSEEGrootFolder    = [filepath 'sub' sub_no{ns}];
    cfg.TrialStart          = -2000; % -1000 v1
    cfg.TrialEnd            = 2000; % 1000 v1
    cfg.PulseLen            = 20;
    cfg.PulseShift          = -5;%-2 v1
    cfg.BaseLine            = [-550 -50];%[-500 -200] v1
    cfg.plottimes           = [25,45,60,100,180,280];%[30,45,60,75,100,150,200,300];
    cfg.NameProject         = 'clean'; % sham, active,clean
    cfg.NameCond            = 'control'; % pre, post, control
    cfg.NameSub             = ['sub' num2str(ns)];
    % Required fields:
    %   cfg.EventCode        - Event code for the TMS pulses in cfg.EventCode (needs to be a string or a numeric value)
    %   cfg.TMSEEGrootFolder - Path of the folder to store results
    %   Optional fields:
    %   cfg.TrialStart       - Start of each epoch in ms (default: -1000)
    %   cfg.TrialEnd         - End of each epoch in ms (default: 1500)
    %   cfg.PulseLen         - Length of the TMS pulse to be removed and interpolated in
    %                          ms (default: 10)
    %   cfg.PulseShift       - Time (ms) from zero where you want to start removing the
    %                          TMS pulse (make sure it's before the rise of the pulse. Default: -2)
    %   cfg.BaseLine         - Baseline time window in ms (default: [-300, -100])
    %   cfg.plottimes        - Times (ms) at which topo plots are shown in the butterfly plots (default: [15,25,40,60,75,100,150,200,300])
    %   cfg.NameProject      - Name of the project (default: '')
    %   cfg.NameCond         - Name of the condition (default: '')
    %   cfg.NameSub          - Name of the subject (default: '')
    
    %% step1: SET UP INPUT VARIABLES
    if isfield(cfg, 'EventCode')
        if ischar(cfg.EventCode) || isnumeric(cfg.EventCode)
            EventCode = cfg.EventCode;
        else
            error('EventCode NEEDS TO BE A STRING OR NUMBER!');
        end
    else
        error('EVENT INFO MISSING!');
    end
    if isfield(cfg, 'TrialStart')
        TrialStart = cfg.TrialStart;
    else
        TrialStart = -1000;% Default start of each epoch: -1000 ms
        cfg.TrialStart = TrialStart;
    end
    if isfield(cfg, 'TrialEnd')
        TrialEnd = cfg.TrialEnd;
    else
        TrialEnd = 1500;% Default end of each epoch: 1500 ms
        cfg.TrialEnd = TrialEnd;
    end
    T = TrialEnd - TrialStart;
    if isfield(cfg, 'PulseLen')
        PulseLen = cfg.PulseLen;
        cfg.TMSlength = cfg.PulseLen; % USED IN CLASSIFYDECAYART
    else
        PulseLen = 10;% length of the TMS pulse to be removed and interpolated (ms)
        cfg.PulseLen = PulseLen;
        cfg.TMSlength = cfg.PulseLen; % USED IN CLASSIFYDECAYART
    end
    if isfield(cfg, 'PulseShift')
        PulseShift = cfg.PulseShift;
    else
        PulseShift = -2;% time (ms) from zero where you want to start removing the TMS pulse (make sure it's before the rise of the pulse)
        cfg.PulseShift = PulseShift;
    end
    
    if isfield(cfg, 'BaseLine')
        BaseLine = cfg.BaseLine;
    else
        BaseLine(1) = -300;
        BaseLine(2) = -100;% Baseline window
        cfg.BaseLine = BaseLine;
    end
    if ~isfield(cfg, 'plottimes')
        cfg.plottimes = [15,25,40,60,75,100,150,200,300];
    end
    if ~isfield(cfg, 'TMSEEGrootFolder')
        error('NEED A DIRECTORY TO STORE DATA')
    end
    if ~isfield(cfg, 'NameProject'); cfg.NameProject = ''; end
    if ~isfield(cfg, 'NameCond'); cfg.NameCond = ''; end
    if ~isfield(cfg, 'NameSub'); cfg.NameSub = ''; end
    
    disp(['Using ' num2str(cfg.TrialStart) 'ms before TMS pulse'])
    disp(['Using ' num2str(cfg.TrialEnd) 'ms after TMS pulse'])
    disp(['Using ' num2str(cfg.BaseLine(1)) 'ms to ' num2str(cfg.BaseLine(2)) 'ms for baseline'])
    disp(['Will interpolate ' num2str(cfg.PulseLen) 'ms around the TMS pulse'])
    disp('If any of these settings are wrong, stop now and rerun');pause(2)
    
    %%  step2: SET UP FOLDERS FOR SAVING VISUALIZATION
    if ~isdir([cfg.TMSEEGrootFolder filesep cfg.NameProject]);mkdir([cfg.TMSEEGrootFolder filesep cfg.NameProject]);end;cd([cfg.TMSEEGrootFolder filesep cfg.NameProject])
    if ~isdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub]);mkdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub]);end;cd([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub])
    if ~isdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond]);mkdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond]);end;cd([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond])
    cfg.folderPath = ([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond]);
    cfg.fullCondName = [cfg.NameProject '_' cfg.NameSub '_' cfg.NameCond];
    if ~isdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond filesep 'ICA1']);mkdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond filesep 'ICA1']);end;
    if ~isdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond filesep 'ICA2']);mkdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond filesep 'ICA2']);end;
    if ~isdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond filesep 'QC']);mkdir([cfg.TMSEEGrootFolder filesep cfg.NameProject filesep cfg.NameSub filesep cfg.NameCond filesep 'QC']);end;
    
    %%  step3: REMOVE THE PULSE ARTIFACT AND DOWNSAMPLE TO 1KHZ
    disp('REMOVING THE PULSE ARTIFACT')
    EEGbeforeRemoval = EEG; % SAVE FOR PLOTTING PURPOSES
    cnt = 1;clear PulseStart PulseEnd
    RemoveLastPulse = 0;
    RemoveFirstPulse = 0;
    for kk = 1:length(EEG.event)
        if strcmp(num2str(EEG.event(kk).type), num2str(EventCode))
            TStart = ceil(EEG.event(kk).latency + TrialStart/1000*EEG.srate);
            if TStart <= 0
                warning('Insufficient preTMS data for the first trial! Remove the first event');
                RemoveFirstPulse = 1;
                FirstPulseEvent = kk;
                continue;
            end
            TEnd = ceil(EEG.event(kk).latency + TrialEnd/1000*EEG.srate - 1);
            if TEnd > EEG.pnts
                warning('Insufficient postTMS data for the last trial! Remove the last event');
                RemoveLastPulse = 1;
                break;
            end
            PulseStart(cnt) = ceil(EEG.event(kk).latency + (PulseShift/1000*EEG.srate));
            PulseEnd(cnt) = ceil(PulseStart(cnt) + floor(PulseLen/1000*EEG.srate));
            x = ceil([TStart : PulseStart(cnt)-1,PulseEnd(cnt) + 1 : TEnd]);
            y = EEG.data(:, x);
            xi = floor(PulseStart(cnt) : PulseEnd(cnt));
            EEG.data(:,xi) = interp1(x, y', xi)';
            if cnt == 1
                firstPulseMS = EEG.event(kk).latency/(EEG.srate/1000);
            end
            cnt = cnt + 1;
        end
    end
    
    RemovePoint = [];
    if RemoveFirstPulse == 1
        RemovePoint = [RemovePoint; 1 EEG.event(FirstPulseEvent).latency + 1];
    end
    
    if RemoveLastPulse == 1
        RemovePoint = [RemovePoint; EEG.event(kk).latency - 1 EEG.pnts];
    end
    
    if ~isempty(RemovePoint)
        EEG = pop_select(EEG, 'nopoint', RemovePoint);
    end
    
    disp(['FOUND:' num2str(cnt) ' EVENTS TO REMOVE PULSE'])
    cfg.srate_orig = EEG.srate;
    if EEG.srate >= 1000
        disp('Downsampling to 1KHz')
        EEG = pop_resample(EEG, 1000);
    else
        error('Sampling rate should be greater than 1000 Hz!');
    end
    disp('DOWNSAMPLING TO 1KHz')
    
    %%  !!! REMOVE BAD CHANNELS for specific subjects
    artchan = [2,39,43,52,61]; % 'C4'=24, 'F5'=36
    disp(['BAD CHANNELS: ' num2str(artchan)]); pause(1)
    EEG = eeg_interp(EEG, artchan);
    
    %% step4: PLOT THE PULSE ARTIFACT BEFORE AND AFTER PULSE REMOVAL REJECTION
    close all; figure
    plot(EEGbeforeRemoval.times-firstPulseMS,squeeze(EEGbeforeRemoval.data(10,:,:)),'b','LineWidth',1); hold all; box off;
    plot(EEG.times-firstPulseMS,squeeze(EEG.data(10,:,:)),'r','LineWidth',1); hold all; box off;
    line([PulseShift/(EEG.srate/1000) PulseShift/(EEG.srate/1000)+PulseLen],[500 500],'Color','k','LineStyle','--','LineWidth',2); hold all;
    line([PulseShift/(EEG.srate/1000) PulseShift/(EEG.srate/1000)],[-500 500],'Color','k','LineStyle','--','LineWidth',2); hold all;
    text(PulseShift/(EEG.srate/1000), 1500, [num2str(PulseShift/(EEG.srate/1000)) 'ms Pulse Shift, ' num2str(PulseLen) 'ms Pulse Duration']);
    xlim([-20 50]); box off; xlabel('Time (ms)'); ylabel('uV');
    % REMOVE UNDERLINES SO TITLE DOESN'T HAVE UNDERSCORES
    tt = cfg.fullCondName;tm = strfind(tt,'_');for ji = 1:length(tm);tt(tm(ji)) = ' ';end;title(tt);
    h = legend({'Before Pulse Removal';'After Pulse Removal'},'Location','SouthEast');legend boxoff
    cd(cfg.folderPath);cd('QC');savefig([cfg.fullCondName '_QC_PulseRemoval'],16,16,150,'',4,[10 8]);
    
    %% step5: DETREND THE CONTINUOUS DATA
    disp('DETRENDING CONTINUOUS DATA')
    stimtrial = [];
    for kk = 1:length(EEG.event)
        if strcmp(num2str(EEG.event(kk).type), num2str(EventCode))
            stimtrial = [stimtrial kk];
        end
    end
    for kk = 1:length(stimtrial)-1
        time1 = round(TrialStart/1000*EEG.srate + EEG.event(stimtrial(kk)).latency);
        time2 = round(TrialStart/1000*EEG.srate + EEG.event(stimtrial(kk + 1)).latency - 1);
        EEG.data(:, time1 : time2) = detrend(EEG.data(:, time1 : time2)','linear')';
    end
    % LAST EPOCH
    time = round([TrialStart/1000*EEG.srate:TrialEnd/1000*EEG.srate] + EEG.event(stimtrial(end)).latency);
    EEG.data(:, time) = detrend(EEG.data(:, time)','linear')';
        
    %% step6: EPOCHED DATA
    disp('EPOCHING')
    EEGepoch = pop_epoch(EEG, {EventCode}, [TrialStart/1000 TrialEnd/1000]);
    
    %% step7: 1ST STAGE - USE ICA ON THE CONTINUOUS DATA TO REMOVE BIG-AMPLITUDE DECAY ARTIFACT
    disp('ICA ROUND 1 - REMOVE LARGE AMPLITUDE DECAY ARTIFACT');pause(1)
    
    EEGepoch = pop_reref(EEGepoch,[]);
    EEG = pop_reref(EEG,[]);
    
    % INFOMAX ICA
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 0, 'pca', EEG.nbchan - 1, 'interupt', 'off');
    % COPY THE ICA PARAMETERS TO THE EPOCHED DATA
    EEGepoch.icaweights = EEG.icaweights;
    EEGepoch.icasphere = EEG.icasphere;
    % DETECT BIG-AMPLITUDE DECAY COMPONENTS FROM THE EPOCHED DATA
    artcomp = classifydecayart(EEGepoch, cfg);
    
    % PLOT ICA MAPS
    ARTIST_plotICA(EEGepoch, cfg, artcomp, 1);
    
    EEG = pop_subcomp(EEG, artcomp);
    
    %% step8: BANDPASS FILTER CONTINUOUS DATA BETWEEN 1 AND 100 Hz
    disp('BANDPASS FILTER 1:100Hz');pause(1)
    lcut = 1;hcut = 100;
    EEG = pop_eegfiltnew(EEG, lcut, 0);
    EEG = pop_eegfiltnew(EEG, 0, hcut);
    
    %% step9: NOTCH FILTERING TO REMOVE LINE NOISE
    disp('NOTCH FILTER (LINE NOISE REMOVAL) at 50 Hz');pause(1)
    EEG = pop_eegfiltnew(EEG, 48, 52, 2000*EEG.srate/1000, 1);
    
    %% step10: 2ND STAGE - IDENTIFY AND REMOVE BAD TRIALS AND CHANNELS
    disp('REMOVE BAD TRIALS AND CHANNELS');pause(1)
    EEG = pop_epoch(EEG, {EventCode}, [TrialStart/1000 TrialEnd/1000]);
    % IDENTIFY AND REMOVE BAD TRIALS
    [arttrial] = identifyarttrial(EEG, cfg);
    disp(['BAD TRIALS: ' num2str(arttrial')]); pause(1)
    EEG = pop_select(EEG, 'notrial', arttrial);
    % IDENTIFY AND REMOVE BAD CHANNELS
    % Set cfg.isransac = 1 to remove bad electrode clusters. But Use it
    % with caution as ransac tend to remove too many channels
    artchan = identifyartchan(EEG,cfg);
    disp(['BAD CHANNELS: ' num2str(artchan)]); pause(1)
    EEG = eeg_interp(EEG, artchan);
    
    %% step11: 3RD STAGE - REMOVE THE REMAINING ARTIFACTS
    disp('ICA ROUND 2 - REMOVE REMAINING ARTIFACTS');pause(1)
    % COMMON AVERAGE REFERENCE
    EEG = pop_reref(EEG, []);
    % DETERMINE OPTIMAL COMPONENT NUMBER NUMBER USING PCA
    CovM = EEG.data(:, :)*EEG.data(:, :)'/size(EEG.data(:, :), 2);
    [~, D] = eig(CovM);
    d = sort(diag(D), 'descend');
    dd = zeros(1, length(d));
    for l = 1:length(d)
        dd(l) = sum(d(1:l));
    end
    cntNumCompForICA = sum(dd/dd(end) < .999);
    % INFOMAX ICA
    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 0, 'pca', cntNumCompForICA, 'interupt', 'off');
    % LOAD TRAINED CLASSIFIER
    load classifierweight.mat;
    cfg.w = classifier.w;
    cfg.b = classifier.b;
    artcomp = classifyartcomp(cfg, EEG);
    
    % PLOT ICA MAPS
    ARTIST_plotICA(EEG, cfg, artcomp, 2)
    
    EEG = pop_subcomp(EEG, artcomp);
    %% step12: COMMON AVERAGE REFERENCE
    disp('COMMON AVERAGE REFERENCING'); pause(1)
    EEG = pop_reref(EEG, []);
    
    %% step13: BASELINE CORRECTION w.r.t -550 ~ -50 ms
    disp('BASELINE CORRECTION'); pause(1)
    EEG = pop_rmbase(EEG, [BaseLine(1), BaseLine(2)]);
    % save data
    cd(filepath)
    eeg_all{1,ns}=EEG;
    save('control_v3.mat','eeg_all','cfg','-v7.3')
    %save('sham_post_v3.mat','eeg_all','cfg','-v7.3')
    %save('active_post_v3.mat','eeg_all','cfg','-v7.3')

    %% step14: PLOT TEP AND TOPOGRAPHY
    close all;figure;xlimm = [-500 500];colormap(jet);fnameTitle = cfg.fullCondName;
    fnameTitletmp = strfind(fnameTitle,'_');for bb = 1:length(fnameTitletmp);fnameTitle(fnameTitletmp(bb)) = ' ';end
    tit = sprintf([fnameTitle ' \n Total trials: ' num2str(size(EEG.data,3)) ', Srate: ' num2str(round(cfg.srate_orig)) 'Hz, \n '...
        num2str(length(artchan)) ' Bad Channels, ' num2str(length(arttrial)) ' Bad Trials, '...
        num2str(length(artcomp)) ' Bad Components']);
    title(tit,'FontSize',10);
    dat = squeeze(mean(EEG.data(:,EEG.times>xlimm(1) & EEG.times<xlimm(2),:),3));
    timtopo(dat,EEG.chanlocs,[xlimm(1),xlimm(2),-1.5*max(max(dat)),1.5*max(max(dat))],cfg.plottimes,'',0,0,'shrink','on');box off
    cd(cfg.folderPath);cd('QC');
    savefig([cfg.fullCondName '_TEP'],16,16,150,'',4,[10 8]);
    disp('SAVING FILE...')
end

%% PLOT AVERAGED TEPS
k=1;
eeglab;
%data = pre;
%data = pre_EEG;
data = eeg_all;
for i =  1:length(data)
    if ~isempty(data{1,i})
   temp(:,:,k)=squeeze(mean(data{1,i}.data,3));
   k=k+1;
    end
end
EEG.data=temp;
EEG.times = data{1, 1}.times;
EEG.chanlocs = data{1, 1}.chanlocs;
cfg.plottimes           = [30,45,60,100,200];
figure;xlimm = [-100 500];colormap(jet);
dat = squeeze(mean(EEG.data(:,EEG.times>xlimm(1) & EEG.times<xlimm(2),:),3));
%timtopo(dat,EEG.chanlocs,[xlimm(1),xlimm(2),-1.5*max(max(dat)),1.5*max(max(dat))],cfg.plottimes,'',0,0,'shrink','on');box off
timtopo(dat,EEG.chanlocs,[xlimm(1),xlimm(2),-8,8],cfg.plottimes,'',0,0,'shrink','on');box off
    
    
