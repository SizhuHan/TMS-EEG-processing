%% step 1: compute forward model
clear;
clc;
mainpath =cd;
cd([mainpath '/TMS+EEGstudy/TMS+EEG_controldata/'])
load(['reformat_control_v3.mat']);
cfg=[];
data=[];
data=ft_timelockanalysis(cfg,data);

% load template headmodel/volume conduction model
load('standard_bem.mat','vol'); % dipoli

% compute source model - dipoli source models
cfg             = [];
cfg.headmodel   = vol; % used to estimate extent of grid
cfg.resolution  = 8; % a source per 8 mm
cfg.inwardshift = 0; % moving sources 0mm inwards from the skull, ...
% since BEM models may be unstable here
sourcemodel = ft_prepare_sourcemodel(cfg);

% load elec position
elec = ft_read_sens('standard_1020.elc');
temp = elec.label;

% realignment step 1: coregister the electrodes to scalp
cfg = [];
cfg.method    = 'project'; % onto scalp surface
cfg.elec      = elec;
cfg.headshape = vol.bnd(1); % scalp surface
elec = ft_electroderealign(cfg);


% realignment step 2: manual modulation (optional)
cfg = [];
cfg.method    = 'interactive';
cfg.elec      = elec;
cfg.headshape = vol.bnd(1);% scalp
elec = ft_electroderealign(cfg);
elec.label = temp;

% visual check
inside = sourcemodel;
outside = sourcemodel;

inside.pos = sourcemodel.pos(sourcemodel.inside, :);
outside.pos = sourcemodel.pos(~sourcemodel.inside, :);

figure;hold on
ft_plot_mesh(inside, 'vertexsize', 20, 'vertexcolor', 'red');
%ft_plot_mesh(outside, 'vertexsize', 20)
ft_plot_sens(elec, 'elecsize', 40, 'color','blue');
%ft_plot_headshape(sourcemodel);
ft_plot_vol(vol, 'facealpha', 0.5);
ft_plot_mesh(vol.bnd(3),'edgecolor',[1,1,0],'facecolor','none'); %brain
%view(125, 10)
view(90, 0)

% computing forward model
cfg             = [];
cfg.elec        = elec;   % sensor information
cfg.channel     = data_clean{1, 1}.label;  % the used channels
cfg.grid        = sourcemodel;   % source points
cfg.headmodel   = vol;   % volume conduction model
cfg.lcmv.reducerank = 3; % 3 for eeg, 2 for meg
cfg.lcmv.batchsize = 5000;
cfg.normalize       = 'yes';
leadfield   = ft_prepare_leadfield(cfg);

cd([mainpath '/TMS+EEGstudy/'])
save('Files4source.mat','elec','vol','leadfield');

%% step 2: source estimation: ERP-based source activity
sub_no = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14',...
    '15','16','17','18','19','20','21','22','23','24','25','26','27','28','29',...
    '30','31','32','33','34','35','36','37','38','39','40','41','42','43'};
subpath = {'active/pre/','active/post/','sham/pre/','sham/post/'};
cond = {'active_pre','active_post','sham_pre','sham_post','HC'};
time = {[0.025 0.035],[0.040 0.070],[0.080 0.120],[0.175 0.225]};

tic;
for nc = 5
    cd([mainpath '/TMS+EEGstudy/TMS+EEG_controldata/'])
    load(['reformat_control_v3.mat']);
    %cd([mainpath '/TMS+EEGstudy/TMS+EEG_MDDdata/' subpath{1,nc}])
    %load(['reformat_' cond{1,nc} '_v3.mat']);
    mkdir(['Source_' cond{1,nc}]);
    cd(['Source_' cond{1,nc}])
    k=1;
    for ns=[1:64] % control
        %for ns=[1:28] % active
        %for ns=[1:25] % sham
        fprintf('Loading for for subject %s\n', sub_no{ns});
        %% redefine trials and compute covariance
        cfg=[];
        cfg.toilim = [-0.50 0.50];
        data_clean{1,ns} = ft_redefinetrial(cfg, HC_clean{1,ns});
        
        cfg=[];
        cfg.covariance = 'yes';
        cfg.covariancwindow = [-0.50 0.50];
        avgElec =ft_timelockanalysis(cfg, HC_clean{1,ns}); % timelock
        %%  source analysis
        cfg                 = [];
        cfg.method          = 'lcmv';
        cfg.elec            = elec;
        cfg.grid            = leadfield; % forward model
        cfg.headmodel       = vol;
        cfg.rawtrial        = 'no';
        cfg.lcmv.keepfilter ='yes';
        cfg.lcmv.lambda     = '5%';
        cfg.lcmv.fixedori   = 'yes';
        AvgFullSource           = ft_sourceanalysis(cfg,avgElec); % compute sptial filter
        
        %% compute ERP-based source activity
        % average across time the dipole moments within a component latency range
        tmpmom = AvgFullSource.avg.mom(AvgFullSource.inside);
        source_tw=[];
        for t = 1:length(time)
            ind    = find(AvgFullSource.time>=time{1,t}(1) & AvgFullSource.time<=time{1,t}(2));
            mom    = AvgFullSource.avg.pow(AvgFullSource.inside);
            for ii = 1:length(tmpmom)
                mom(ii) = mean(abs(tmpmom{ii}(ind)));
            end
            
            % insert the component amplitude in the 'pow' field and save to disk, the
            % original pow contains the mean amplitude-squared across the
            % time-window used for the channel-level covariance computation
            source_tw{t} = AvgFullSource;
            source_tw{t}.avg.pow(source_tw{t}.inside) = abs(mom);
            source_tw{t}.cfg = rmfield(source_tw{t}.cfg, {'headmodel' 'callinfo'}); % this is removed because it takes up a lot of memory
            
            TempData = source_tw{t};
            TempData.cfg=[];% to save memory
            switch t
                case 1
                    all_p30{k}  = TempData;
                case 2
                    all_p60{k}  = TempData;
                case 3
                    all_n100{k} = TempData;
                case 4
                    all_p200{k} = TempData;
            end
        end
        k=k+1;
        save(['AvgSource_' cond{1,nc} '_sub' sub_no{1,ns} '.mat'],'AvgFullSource','source_tw')
    end
    save(['AvgSource_ERP_' cond{1,nc} '.mat'],'all_*','-v7.3')
end
%clear HC*
toc;
%%
% grand average
cfg           = [];
cfg.parameter = 'avg.pow';
source_allsub = ft_sourcegrandaverage(cfg, all_p200{:});
% interpolate onto MRI
mri              = ft_read_mri('single_subj_T1.nii');
mri = ft_volumereslice([], mri);
cfg              = [];
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int       = ft_sourceinterpolate(cfg, source_allsub, mri);
% plot
%aal               = ft_read_atlas('/Users/tristahan/Lib_Matlab/fieldtrip-20161231/template/atlas/aal/ROI_MNI_V4.nii');
cfg               = [];
cfg.method        = 'ortho';
%cfg.method       = 'surface';
cfg.funparameter  = 'pow';
%cfg.surffile       = 'surface_pial_left.mat';
%cfg.maskparameter = 'mask';
%cfg.location = [-42 -18 67];
%cfg.atlas         = aal;
cfg.opacitymap    = 'rampup';
%cfg.opacitylim    = [0 2e3];
cfg.funcolorlim = [0 2e4];
ft_sourceplot(cfg,source_int);
%view(150,10)

%% step 3: source estimation: compute binary matrix and SCD
sub_no = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15'...
    '16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'};
subpath = {'active/pre/','active/post/','sham/pre/','sham/post/'};
cond = {'active_pre','active_post','sham_pre','sham_post','HC'};

tic;
for nc = 5%1:5
    if nc==5
        cd([mainpath '/TMS+EEGstudy/TMS+EEG_controldata/'])
        load(['reformat_control_v3.mat']);
        data_clean = HC_clean;
        clear HC_clean
        mkdir(['Source2_' cond{1,nc}])
        cd(['Source2_' cond{1,nc}])
    else
        cd([mainpath '/TMS+EEGstudy/TMS+EEG_MDDdata/' subpath{1,nc}])
        load(['reformat_' cond{1,nc} '_v3.mat']);
        mkdir(['Source2_' cond{1,nc}])
        cd(['Source2_' cond{1,nc}])
    end
    
    switch nc
        case {1,2} % active
            range = [1:28];
        case {3,4} % sham
            range = [1:25];
        case 5   % control
            range = [1:64];
    end
    for ns= [1:64]%range
        fprintf('Loading for for subject %s\n', sub_no{ns});
        %% redefine trials and compute covariance
%         cfg=[];
%         cfg.toilim = [-0.50 -0.001];
%         pre = ft_redefinetrial(cfg, data_clean{1,ns});
%         cfg=[];
%         cfg.covariance = 'yes';
%         cfg.keeptrials = 'no';
%         AvgPreElec =ft_timelockanalysis(cfg,pre);
%         [~,s,~] = svd(AvgPreElec.cov);
%         s=diag(s);
%         d = -diff(log10(s)); d = d./std(d);
%         lambda=s(find(d>5,1,'first'));
%         
%         clear pre AvgPreElec
        
        cfg=[];
        cfg.toilim = [0.001 0.50];
        post = ft_redefinetrial(cfg, data_clean{1,ns});
        cfg=[];
        cfg.covariance = 'yes';
        cfg.keeptrials = 'no';
        AvgPostElec =ft_timelockanalysis(cfg,post);
        %
        cfg=[];
        cfg.toilim = [-0.50 0.50];
        full = ft_redefinetrial(cfg, data_clean{1,ns});
        cfg=[];
        cfg.covariance = 'yes';
        cfg.keeptrials = 'no';
        cfg.covariancwindow = [-0.50 0.50];
        AvgFullElec =ft_timelockanalysis(cfg, full); % timelock
        cfg.keeptrials = 'yes';
        AllFullElec =ft_timelockanalysis(cfg, full);
        %%  source analysis
        cfg               = [];
        cfg.method        = 'lcmv';
        cfg.elec          = elec;
        cfg.grid          = leadfield; % forward model
        cfg.headmodel     = vol;
        cfg.rawtrial      = 'no';
        cfg.lcmv.lambda        = '5%';%'0.00005%';
        cfg.lcmv.weightnorm    = 'nai';
        cfg.lcmv.projectnoise  = 'yes';
        cfg.lcmv.keepfilter    = 'yes';
        cfg.lcmv.fixedori      = 'yes';
        AvgFullSource          = ft_sourceanalysis(cfg,AvgFullElec); % compute sptial filter
        commonfilter = cell2mat(AvgFullSource.avg.filter);
        
%         cfg               = [];
%         cfg.method        = 'lcmv';
%         cfg.elec          = elec;
%         cfg.grid          = leadfield; % forward model
%         cfg.headmodel     = vol;
%         cfg.rawtrial      = 'no';
%         cfg.lcmv.lambda        = lambda;
%         cfg.lcmv.weightnorm    = 'nai';
%         cfg.lcmv.projectnoise  = 'yes';
%         cfg.lcmv.fixedori = 'yes';
        cfg.lcmv.filter = AvgFullSource.avg.filter;%commonfilter;% important!!!
        
        AvgPostSource     = ft_sourceanalysis(cfg,AvgPostElec); % average post-tms
        AvgPostSource.cfg = rmfield(AvgPostSource.cfg, {'headmodel' 'callinfo'}); % this is removed because it takes up a lot of memory
        
                
        alldata=[];
        for i = 1:size(AllFullElec.trial,1)
            alldata(:,:,i) = commonfilter * squeeze(AllFullElec.trial(i,:,:));% voxcel * time * trials
        end
        avgdata = mean(alldata,3);
        
        %avgdata = cell2mat(AvgFullSource.avg.mom);
        
        clear AvgPostElec AvgFullSource               
        %% compute binary matrix
        
        pre_origin = avgdata(:,1:500);
        post_origin = avgdata(:,502:1001);
        % normalization for original dataset
        origin_mean = mean(pre_origin,2);
        origin_var = std(pre_origin,0,2);
        pre_norm = (pre_origin - origin_mean)./origin_var; % voxcel * time
        post_norm = (post_origin - origin_mean)./origin_var; % voxcel * time
        
        mergedata = [pre_origin,post_origin];
        c = size(mergedata,2); % number of time points
        max_p =[];thresh=[];
        % permutation
        for nt = 1:1000
            disp(['running perm set no.: ', num2str(nt)]);
            d=randperm(c);
            pre_perm = mergedata(:,d(1:c/2));
            post_perm = mergedata(:,d(c/2+1:c));
            perm_mean = mean(pre_perm,2);
            perm_var = std(pre_perm,0,2);
            perm_norm = (post_perm-perm_mean)./perm_var; % normalization for perm dataset
            max_p(nt,:)=max(abs(perm_norm),[],1);
        end
        
        Mask = zeros(size(post_origin,1),size(post_origin,2)); % voxcel * time
        for t = 1:size(post_origin,2) % time points
            temp = sort(max_p(:,t));
            thresh(t) = temp(990); % alpha = 0.01
            Mask(find(abs(post_norm(:,t)) > thresh(t)),t)=1;
        end
        close all;
        figure;plot(1:500,max(post_norm,[],1),1:500,max(pre_norm,[],1),1:500,thresh);
        title([cond{1,nc} ' Sub' sub_no{1,ns}]);
        
        save(['PostSource_' cond{1,nc} '_sub' sub_no{1,ns} '.mat'],'AvgPostSource','Mask','post_origin')
        
    end
end
%clear HC*
toc;

%% parameter setting
sub_no = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14',...
    '15','16','17','18','19','20','21','22','23','24','25','26','27','28','29',...
    '30','31','32','33','34','35','36','37','38','39','40','41','42','43',...
    '44','45','46','47','48','49','50','51','52','53','54','55','56','57','58',...
    '59','60','61','62','63','64','65','66','67','68','69','70','71','72'};
subpath = {'active/pre/','active/post/','sham/pre/','sham/post/'};
cond = {'active_pre','active_post','sham_pre','sham_post','HC'};
%ROI_Name = {'Frontal_Mid_L','Frontal_Sup_L'}; %ROI_MNI_V4
ROI_Name = {{'Frontal_Mid_2_L'},{'Amygdala_L','Amygdala_R'},...
    {'ACC_sub_L','ACC_sub_R'},{'Hippocampus_L','Hippocampus_R'}...
    {'OFCant_L','OFCant_R'},{'Parietal_Sup_L','Parietal_Sup_R'}...
    {'Precentral_L','Precentral_R'},{'ACC_sup_L','ACC_sup_R'}}; %ROI_MNI_V7

aal   = ft_read_atlas('/Users/tristahan/Lib_Matlab/fieldtrip-20161231/template/atlas/aal/ROI_MNI_V7.nii');
aal.coordsys = 'mni';
aal.tissue = aal.brick0;
aal.tissuelabel = aal.brick0label;
load('/Users/tristahan/Lib_Matlab/spm12/toolbox/AAL3/ROI_MNI_V7_List.mat');% label name
for i=1:length(ROI)
    aal.tissuelabel{ROI(i).ID} = ROI(i).Nom_L;
end
aal = rmfield(aal,{'brick0label','brick0'});

roi_id = 1;
% find coordinates of ROIs
%aal               = ft_read_atlas('/Users/tristahan/Lib_Matlab/fieldtrip-20161231/template/atlas/aal/ROI_MNI_V4.nii');
cfg=[];
cfg.atlas = aal;
cfg.inputcoord = 'mni';
cfg.roi = ROI_Name{roi_id};
roi_mask = ft_volumelookup(cfg,AvgPostSource);
roi_mask =logical(reshape(roi_mask,[18*23*18,1]));
roi_loc = AvgPostSource.pos(roi_mask,:);
target_loc = repmat([-36,26,42],size(roi_loc,1),1);% seed at lDLPFC
distance = sqrt(sum((roi_loc-target_loc).^2,2));
%% compute SCD and SCS for ROIs
for nc=1:5
    ROI_CD=[];ROI_SCD=[];ROI_SCS=[];
    switch nc
        case {1,2}
            cd(['/Volumes/SizhuFiles/TMS+EEGstudy/TMS+EEG_MDDdata/', subpath{1,nc} 'Source2_' cond{1,nc}])
            range =[1:28];   % active        
        case {3,4}
            cd(['/Volumes/SizhuFiles/TMS+EEGstudy/TMS+EEG_MDDdata/', subpath{1,nc} 'Source2_' cond{1,nc}])
            range = [1:25]; % sham
        case 5
            cd(['/Volumes/SizhuFiles/TMS+EEGstudy/TMS+EEG_controldata/','Source2_' cond{1,nc}])
            range = [1:64]; % control
    end
    k=1;
    for ns = range
        load(['PostSource_' cond{1,nc} '_sub' sub_no{1,ns} '.mat'],'AvgPostSource','Mask','post_origin');
        TimeCourse = abs(post_origin).* Mask;  % voxcel * time
        for t = 1:length(AvgPostSource.time)
            [ns, t]
            AvgPostSource.avg.pow(AvgPostSource.inside) = abs(post_origin(:,t));
            ROI_CD(k,t) = nansum(AvgPostSource.avg.pow(roi_mask));
            
            AvgPostSource.avg.pow(AvgPostSource.inside) = TimeCourse(:,t);            
            ROI_SCD(k,t) = nansum(AvgPostSource.avg.pow(roi_mask));
                                   
            AvgPostSource.avg.pow(AvgPostSource.inside) = Mask(:,t);
            ROI_SCS(k,t) = nansum(AvgPostSource.avg.pow(roi_mask).* distance);

        end
        k=k+1;
    end
    
    switch nc
        case 1
            Active_pre_CD = ROI_CD;
            Active_pre_SCD = ROI_SCD;
            Active_pre_SCS = ROI_SCS;
        case 2
            Active_post_CD = ROI_CD;
            Active_post_SCD = ROI_SCD;
            Active_post_SCS = ROI_SCS;
        case 3
            Sham_pre_CD = ROI_CD;
            Sham_pre_SCD = ROI_SCD;
            Sham_pre_SCS = ROI_SCS;
        case 4
            Sham_post_CD = ROI_CD;
            Sham_post_SCD = ROI_SCD;
            Sham_post_SCS = ROI_SCS;
        case 5
            HC_CD = ROI_CD;
            HC_SCD = ROI_SCD;
            HC_SCS = ROI_SCS;
    end
end
cd /Volumes/SizhuFiles/TMS+EEGstudy/Source2_Data
save(['Frontal_Mid_2_L.mat'],'Active*','Sham*','HC*'); % 73
%save(['OFCant_LR.mat'],'Active*','Sham*','HC*'); % 19
%save(['Hippocampus_LR.mat'],'Active*','Sham*','HC*'); % 30
%save(['Parietal_Sup_LR.mat'],'Active*','Sham*','HC*'); % 66

%% SCD dynamic plot: plotting
%load(['Frontal_Mid_2_L.mat']);
clear A B
col = {[0.89,0.15,0],[0.13,0.13,0.13],[0.06,0.56,0.25],[0.52,0.52,0.52]};
Pre_SCD = [Active_pre_SCD;Sham_pre_SCD];
n=8;
for i=1:500-9
    A(:,i) = mean(Pre_SCD(:,i:i+9)/n,2);
    B(:,i) = mean(HC_SCD(:,i:i+9)/n,2);
    %A(:,i) = mean(Active_post_SCD(:,i:i+9)/n,2)-mean(Active_pre_SCD(:,i:i+9)/n,2);
    %B(:,i) = mean(Sham_post_SCD(:,i:i+9)/n,2)-mean(Sham_pre_SCD(:,i:i+9)/n,2);
end
y1=mean(A,1);
errBar1=std(A)/sqrt(size(A,1));
uE1=y1+errBar1;
lE1=y1-errBar1;
yP1=[lE1,fliplr(uE1)];

%B=HC_SCD/n;
y2=mean(B,1);
errBar2=std(B)/sqrt(size(B,1));
uE2=y2+errBar2;
lE2=y2-errBar2;
yP2=[lE2,fliplr(uE2)];

x = 5:495;
xP= [x,fliplr(x)];
y=[y1;y2];
yy=NaN(length(x),1);
ind1 = x>150 & x<180;
%yy(ind1)=5e-4;

figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1,...
    'Position',[0.128257839721254 0.186679841897233 0.775 0.785750988142292]);
hold(axes1,'on');
area(x,yy,...
    'FaceColor',[0.9 0.9 0.9],...
    'EdgeColor',[1 1 1]);
plot(x,y1,'LineWidth',1,'Color',col{1,3}); hold on
patch('XData',xP,'YData',yP1,'FaceAlpha',0.2,'EdgeColor','none','FaceColor',col{1,3})
plot(x,y2,'LineWidth',1,'Color',col{1,4});
patch('XData',xP,'YData',yP2,'FaceAlpha',0.2,'EdgeColor','none','FaceColor',col{1,4})
xlabel({'Time (ms)'});
ylabel({'Current Density'});
%ylabel({'SCD'});
%ylabel({'SCS, mm'});
annotation(figure1,'textbox',...
    [0.71702787456446 0.784273166447079 0.078397212543554 0.0968379446640316],...
    'Color',col{1,3},...
    'String',{'Active'},...
    'FontSize',12,...
    'EdgeColor','none');
annotation(figure1,'textbox',...
    [0.71702787456446 0.724273166447079 0.078397212543554 0.0968379446640316],...
    'Color',col{1,4},...
    'String',{'Sham'},...
    'FontSize',12,...
    'EdgeColor','none');
%axis([0,500,0,5e-4])
%axis([0,500,0,50])
%% %%% extract mean SCD %%%%
clear A B C D
Pre_Data = [Active_pre_SCD;Sham_pre_SCD];
HC_Data = HC_SCD;
Active_Data = [Active_pre_SCD;Active_post_SCD];
Sham_Data = [Sham_pre_SCD;Sham_post_SCD];
n=66;
for i=1:500-9
    A(:,i) = mean(Pre_Data(:,i:i+9)/n,2);
    B(:,i) = mean(HC_Data(:,i:i+9)/n,2);
    C(:,i) = mean(Active_Data(:,i:i+9)/n,2);
    D(:,i) = mean(Sham_Data(:,i:i+9)/n,2);
end
x=5:495;
%x=1:500;
toi_1 = find(x >= 164 & x <= 215); 
toi_2 = find(x >= 150 & x <= 185); 
% Pre_avg = mean(Pre_Data(:,toi),2);
% HC_avg = mean(HC_Data(:,toi),2);
Pre_avg = mean(A(:,toi_1),2);
HC_avg = mean(B(:,toi_1),2);
HCvsPats = [HC_avg;Pre_avg];

% Active_avg = mean(Active_Data(:,toi),2);
% Sham_avg = mean(Sham_Data(:,toi),2);
Active_avg = mean(C(:,toi_2),2);
Sham_avg = mean(D(:,toi_2),2);
PrevsPost = [Active_avg;Sham_avg];

active_subtract = Active_avg(29:56) - Active_avg(1:28);
sham_subtract = Sham_avg(26:50) - Sham_avg(1:25);
%subtract_data = [active_subtract;sham_subtract]';




%% SCD/SCS source plot: load data
for nc=1:5
    switch nc
        case {1,2}
            cd(['/Volumes/SizhuFiles/TMS+EEGstudy/TMS+EEG_MDDdata/', subpath{1,nc} 'Source_' cond{1,nc}])
            range =[1:28];   % active        
        case {3,4}
            cd(['/Volumes/SizhuFiles/TMS+EEGstudy/TMS+EEG_MDDdata/', subpath{1,nc} 'Source_' cond{1,nc}])
            range = [1:25]; % sham
        case 5
            cd(['/Volumes/SizhuFiles/TMS+EEGstudy/TMS+EEG_controldata/','Source_' cond{1,nc}])
            range = [1:64]; % control
    end
    k=1;AllPow=[];AllDist=[];
    for ns=range
        fprintf('Loading for subject %s\n', sub_no{ns});
        load(['PostSource_' cond{1,nc} '_sub' sub_no{1,ns} '.mat'],'AvgPostSource','Mask','post_origin');
        AllPow(k,:,:) = abs(post_origin).* Mask;  % nsubj * voxcel * time
        
        Target_loc = repmat([-36,26,42],size(AvgPostSource.pos(AvgPostSource.inside),1),1);% as an example
        Distance = sqrt(sum((AvgPostSource.pos(AvgPostSource.inside)-Target_loc).^2,2));
        AllDist(k,:,:) = repmat(Distance,1,size(Mask,2)).*Mask;
        k=k+1;
    end
    switch nc
        case 1
            Active_pre_pow = AllPow;
            Active_pre_dist = AllDist;
        case 2
            Active_post_pow = AllPow;
            Active_post_dist = AllDist;
        case 3
            Sham_pre_pow = AllPow;
            Sham_pre_dist = AllDist;
        case 4
            Sham_post_pow = AllPow;
            Sham_post_dist = AllDist;
        case 5
            HC_pow = AllPow;
            HC_dist = AllDist;
    end
end
cd('/Volumes/SizhuFiles/TMS+EEGstudy')
save('AllSubSourcePow.mat','*_pow','*_dist')
%% SCD/SCS source plot: plotting

toi = find(AvgPostSource.time >= 0.204 & AvgPostSource.time <= 0.212);
AllPow = Sham_post_pow;
%AllPow = cat(1,Active_pre_pow,Sham_pre_pow);
AvgPow = squeeze(mean(AllPow,1));
AvgPostSource.avg.pow(AvgPostSource.inside) = mean(AvgPow(:,toi),2);
%AllDist = Sham_post_dist;
%AllDist = cat(1,Active_pre_dist,Sham_pre_dist);
%AvgDist = squeeze(mean(AllDist,1));
%AvgPostSource.avg.pow(AvgPostSource.inside) = mean(AvgDist(:,toi),2);

mri                 = ft_read_mri('single_subj_T1.nii');
cfg              = [];
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int  = ft_sourceinterpolate(cfg, AvgPostSource, mri);

% find coordinates of ROIs
%aal               = ft_read_atlas('/Users/tristahan/Lib_Matlab/fieldtrip-20161231/template/atlas/aal/ROI_MNI_V4.nii');
aal         = ft_read_atlas('/Users/tristahan/Lib_Matlab/fieldtrip-20161231/template/atlas/aal/ROI_MNI_V7.nii');
aal.coordsys = 'mni';
aal.tissue = aal.brick0;
aal.tissuelabel = aal.brick0label;
load('/Users/tristahan/Lib_Matlab/spm12/toolbox/AAL3/ROI_MNI_V7_List.mat');% label name
for i=1:length(ROI)
    aal.tissuelabel{ROI(i).ID} = ROI(i).Nom_L;
end
aal = rmfield(aal,{'brick0label','brick0'});

cfg=[];
cfg.atlas = aal;
cfg.inputcoord = 'mni';
%cfg.roi = {'Frontal_Mid_2_L'};
%cfg.roi = {'ACC_sub_L','ACC_sub_R'};
cfg.roi = {'Hippocampus_R'};
cor_mask = ft_volumelookup(cfg,source_int);
source_int.mask = logical(reshape(cor_mask,[91*109*91,1]));

% plot source
cfg               = [];
cfg.method        = 'ortho';
%cfg.method       = 'surface';
cfg.funparameter  = 'pow';
cfg.surffile       = 'surface_pial_left.mat';
%cfg.maskparameter = 'mask';
%cfg.location = [-42 -18 67];
%cfg.atlas         = aal;
%cfg.opacitymap    = 'rampup';
%cfg.opacitylim    = [0 1];
%cfg.funcolorlim = [0 5e-4];
ft_sourceplot(cfg,source_int);
view(90,0)


%% step 3: add the source data into two cell arrays
subpath = {'active/pre/','active/post/','sham/pre/','sham/post/'};
cond = {'active_pre','active_post','sham_pre','sham_post','HC'};

source_hc_all          = [];
k=1;
for sub =[1:64] % control
    fprintf('Loading for for HC group %s\n', sub_no{sub});
    path ='/Volumes/SizhuFiles/TMS+EEGstudy/TMS+EEG_controldata/Source_HC/';
    cd(path);
    % Load hc data
    load(['Source_hc_sub' sub_no{sub} '.mat']);
    % Put these into arrays outside the loop
    source_hc_all{k}     = source;
    k=k+1;
end

source_active_pre_all  = [];
k=1;
for nc=1
    for sub = [1:28] % active
        fprintf('Loading for for active group %s\n', sub_no{sub});
        path= ['/Volumes/SizhuFiles/TMS+EEGstudy/TMS+EEG_MDDdata/' subpath{1,nc}];
        cd([path 'Source_' cond{1,nc}]);
        % Load hc data
        load(['Source_' cond{1,nc} '_sub' sub_no{sub} '.mat']);
        source_active_pre_all{k}     = source;
        k=k+1;
    end
end

source_sham_pre_all  = [];
k=1;
for nc=3
    for sub = [1:25] % sham
        fprintf('Loading for for sham group %s\n', sub_no{sub});
        path= ['/Volumes/SizhuFiles/TMS+EEGstudy/TMS+EEG_MDDdata/' subpath{1,nc}];
        cd([path 'Source_' cond{1,nc}]);
        % Load hc data
        load(['Source_' cond{1,nc} '_sub' sub_no{sub} '.mat']);
        source_sham_pre_all{k}     = source;
        k=k+1;
    end
end

source_pre_all = [source_active_pre_all,source_sham_pre_all];


%% grandaverage source plot
cfg=[];
cfg.latency = [0.08 0.12];
cfg.parameter    = 'avg.pow';
%[grandavg] = ft_sourcegrandaverage(cfg,source_hc_all{:});
[grandavg] = ft_sourcegrandaverage(cfg,source_pre_all{:});

mri                 = ft_read_mri('single_subj_T1.nii');

cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'mom';
cfg.interpmethod = 'nearest';
source_int  = ft_sourceinterpolate(cfg, grandavg, mri);

source_int.mask = source_int.pow > median(source_int.pow(:)); % 30 % of maximum
cfg               = [];
cfg.method        = 'ortho';
%cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
%cfg.location = [-42 -18 67];
cfg.funcolormap = 'jet';
cfg.funcolorlim = [-50 50];
ft_sourceplot(cfg,source_int);


%% Perform Statistical Analysis

cfg                     = [];
cfg.dim                 = source.dim;
%cfg.time             = [0.15 0.21];%'all';
%cfg.avgovertime         = 'yes';
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_indepsamplesT';
cfg.parameter           = 'pow';
cfg.correctm            = 'cluster';
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic    = 'maxsum';
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.025;
cfg.numrandomization    = 1000;

% Design Matrix
Ncond_1 = length(source_pre_all);
Ncond_2 = length(source_hc_all);
cfg.design(1,:)  = [ones(1,Ncond_1) 2*ones(1,Ncond_2)];
cfg.ivar         = 1;

% Perform statistical analysis
[stat]  = ft_sourcestatistics(cfg,source_pre_all{:}, source_hc_all{:});


%% Sourceinterpolate
mri                 = ft_read_mri('single_subj_T1.nii');

cfg                 = [];
cfg.voxelcoord      = 'no';
cfg.parameter       = 'stat';
cfg.interpmethod    = 'nearest';
statint             = ft_sourceinterpolate(cfg, stat, mri);
%sourceint             = ft_sourceinterpolate(cfg, source, mri);

cfg.parameter       = 'mask';
maskint             = ft_sourceinterpolate(cfg, stat,mri);
statint.mask        = maskint.mask;

%% Show raw source level statistics (3D plot)

cfg                = [];
cfg.method         = 'surface';
%cfg.maskparameter  = 'mask';
cfg.funparameter   = 'stat';%'pow'
cfg.projmethod     = 'nearest';
%cfg.surffile        = 'surface_pial_left.mat';
cfg.surfinflated   = 'surface_pial_both.mat';
cfg.camlight       = 'no';
ft_sourceplot(cfg, statint);
%ft_sourceplot(cfg, sourceint);
%colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
view([0, 90]);

light ('Position',[0 0 50])
light ('Position',[0 -50 0])
material dull;
% drawnow;
view([0, 90]);
set(gca,'FontSize',14);

%%
atlas        = ft_read_atlas('ROI_MNI_V4.nii');
atlas.anatomy   = mri.anatomy;
cfg                     = [];
cfg.method              = 'surface';
cfg.projmethod          = 'project';
cfg.camlight            = 'yes';
cfg.surffile            ='surface_pial_left.mat';
% uncomment to project half brain
cfg.locationcoordinates = 'voxel';
%cfg.funcolormap         = cfg.cmap;
cfg.funparameter        = 'tissue';
cfg.atlas               = 'ROI_MNI_V4.nii';
ft_sourceplot(cfg, atlas)
view([90 0]);

%% ortho-plot
cfg                  = [];
cfg.method           = 'ortho';
%cfg.zlim             = 'maxabs';
cfg.funparameter     = 'stat';
cfg.latency          = [0.15,0.2];
%cfg.maskparameter   = 'mask';
cfg.location         = 'max';
cfg.funcolorlim   = [0.0 5];
cfg.opacitylim    = [0.0 5];

ft_sourceplot(cfg,statint);
%ft_sourceplot(cfg,sourceint);
%colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

%% slice plot

maxval = max(statint.stat);
cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'stat';
cfg.nslices       = 50;
%cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 5];
cfg.opacitylim    = [0.0 5];
%cfg.opacitymap    = 'rampup';
ft_sourceplot(cfg, statint);

%% Export to nifti

% Run this if you want to export the clusters rather than the raw stats
statint.stat(isnan(statint.stat)) = 0;
statint.stat                      = (statint.stat(:).*statint.mask(:));

% Use ft_sourcewrite to export to nifti
cfg                         = [];
cfg.filetype                = 'nifti';
cfg.filename                = 'group_clustered';
cfg.parameter               = 'stat';
ft_sourcewrite(cfg,statint);

