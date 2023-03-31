
clear;
clc;
mainpath =cd;
filepath = [mainpath '/TMS+EEG_controldata/'];
%filepath = [mainpath '/TMS+EEGstudy/TMS+EEG_MDDdata/'];
sub_no = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',...
    '16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'};
subpath = {'active/pre/','active/post/','sham/pre/','sham/post/'};
cond = {'active_pre','active_post','sham_pre','sham_post','control'};

%% step 1:transfer from eeglab to fieldtrip format and save data
for nc = 1:length(cond)-1
    data_clean=[];
    load([filepath subpath{1,nc} cond{1,nc} '_v3.mat']);
    %load([filepath cond{1,nc} '_v3.mat']); % nc=5
    for ns = 1:length(eeg_all)
        cd([filepath subpath{1,nc} 'sub' sub_no{ns}])
        %cd([filepath 'sub' sub_no{ns}]) % nc=5
        filename = deblank(ls('*.vhdr'));
        % load data
        cfg                 = [];
        cfg.dataset         = filename;
        cfg.baselinewindow  = [];
        data = ft_preprocessing(cfg);
        
        % epoch
        cfg                     = [];
        cfg.dataset             = filename;
        cfg.trialfun            = 'ft_trialfun_general'; % this is the default
        cfg.trialdef.eventtype  = 'Stimulus';
        cfg.trialdef.eventvalue = {'S255'};
        cfg.trialdef.prestim    = 2.0;
        cfg.trialdef.poststim   = 2.0;
        cfg                     = ft_definetrial(cfg);
        data_seg                = ft_redefinetrial(cfg,data);
        clear data
        
        % resample
        cfg                = [];
        cfg.resamplefs     = 1000;
        data_rsp           = ft_resampledata(cfg,data_seg);
        clear data_seg
        
        % replace data
        temp=eeg_all{1,ns};
        data_rsp.trial=[];
        data_rsp.time=[];
        for i = 1:size(temp.data,3)
            data_rsp.trial{1,i} = squeeze(temp.data(:,:,i));
            data_rsp.time{1,i} = temp.times/1000;
        end
        
        data_rsp.trialinfo = repmat(255,[size(temp.data,3) 1]);
        
        
        data_clean{1,ns} = data_rsp;
        % HC_clean{1,ns} = data_rsp; % nc=5
        clear data_rsp temp
        % save data
        cd([filepath subpath{1,nc}])
        save(['reformat_' cond{1,nc} '_v3.mat'],'data_clean','-v7.3')
        %cd(filepath)
        %save(['reformat_' cond{1,nc} '_v3.mat'],'HC_clean','-v7.3')        
    end
end

%% step 2: calculate grand average for each condition
%load([mainpath '/TMS+EEG_controldata/reformat_control_v3.mat'])
% load([mainpath '/TMS+EEG_MDDdata/preMDD/reformat_preMDD.mat'])
hc_avg=[];
k=1;
for i=[1:22 24:29 32:45 47:61 63:72] % sub 30 46 62 is missing, sub23 unfinished, sub 31 is bad
    cfg =[];
    hc_avg{1,k} =ft_timelockanalysis(cfg, HC_clean{1,i});
    k=k+1;
end
clear HC_clean

filepath = ['/TMS+EEGstudy/TMS+EEG_MDDdata/'];
pre_active=[];
pre_sham=[];
for nc=[1,3]
    load([filepath subpath{1,nc} 'reformat_' cond{1,nc} '_v3.mat'],'data_clean');
    k=1;
    if nc==1
        for i=[1:9 11:18 20:33] % sub 10,19 missing
            pre_active{1,k}= data_clean{1,i};
            k=k+1;
        end
    end
    
    if nc==3
        for i=[1:25 29 34:35] % sub 26:28, 30:33 no marker,36 is missing
            pre_sham{1,k}= data_clean{1,i};
            k=k+1;
        end
    end
end

pre_clean=[pre_active,pre_sham];
for i=1:length(pre_clean)
    cfg=[];
    pre_avg{1,i} =ft_timelockanalysis(cfg, pre_clean{1,i});
end



cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
grand_avgHC  = ft_timelockgrandaverage(cfg, hc_avg{:});
grand_avgPre   = ft_timelockgrandaverage(cfg, pre_avg{:});

active_avg=[];
sham_avg=[];
grand_avg=[];
for nc = 1:4
    load([filepath subpath{1,nc} 'reformat_' cond{1,nc} '_v3.mat'],'data_clean');
    k=1;
    if nc<3
        for i=[1:9 11:18 20:33] % sub 10,19 missing
            cfg=[];
            active_avg{nc,k} =ft_timelockanalysis(cfg, data_clean{1,i});
            k=k+1;
        end
        cfg = [];
        cfg.channel   = 'all';
        cfg.latency   = 'all';
        cfg.parameter = 'avg';
        grand_avg{1,nc}  = ft_timelockgrandaverage(cfg, active_avg{nc,:});
    else
        for i=[1:25 29 34:35] % sub 26:28, 30:33 no marker,36 is missing
            cfg=[];
            sham_avg{nc-2,k} =ft_timelockanalysis(cfg, data_clean{1,i});
            k=k+1;
        end
        cfg = [];
        cfg.channel   = 'all';
        cfg.latency   = 'all';
        cfg.parameter = 'avg';
        grand_avg{1,nc}  = ft_timelockgrandaverage(cfg, sham_avg{nc-2,:});
    end
end

clear data_clean

%% step 2.1 remove 6 patients & 3 healthy subjects
hc_ind = [1:22 24:29 32:45 47:61 63:72];
ex_hc = find(hc_ind ~= 15 & hc_ind ~= 19 & hc_ind ~= 38);% sub 15,19,38 HDRS>7
hc_avg = hc_avg(:,ex_hc);
active_ind = [1:9 11:18 20:33];
ex_ac = find(active_ind ~= 4 & active_ind~=25 & active_ind~=26); % sub 4,25,26 HDRS<20
active_avg = active_avg(:,ex_ac);
sham_ind = [1:25 29 34:35];
ex_sh = find(sham_ind ~= 5 & sham_ind~=8 & sham_ind~=12); % sub 5,8,12 HDRS<20
sham_avg = sham_avg(:,ex_sh);
pre_avg = [active_avg(1,:),sham_avg(1,:)];
grand_avg=[];
for nc =1:4
    if nc<3
      cfg = [];
      cfg.channel   = 'all';
      cfg.latency   = 'all';
      cfg.parameter = 'avg';
      grand_avg{1,nc}  = ft_timelockgrandaverage(cfg, active_avg{nc,:});
    else
      cfg = [];
      cfg.channel   = 'all';
      cfg.latency   = 'all';
      cfg.parameter = 'avg';
      grand_avg{1,nc}  = ft_timelockgrandaverage(cfg, sham_avg{nc-2,:});
    end
end

save('TEP_data_rmsubj.mat','*avg')

%% step 3:  grand average plot

time = {[0.025 0.035],[0.040 0.070],[0.080 0.120],[0.160 0.200]};
tle = {'P30','P60','N100','P180'};
z={[-2,2],[-2,2],[-2,2],[-2,2]};
cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';

figure;
for t=1:length(time)
    cfg.xlim = time{1,t};
    cfg.zlim = z{t};
     subplot(2,length(time),t)
    ft_topoplotER(cfg,grand_avgHC);
    %title(tle{t});
    subplot(2,length(time),t+length(time))
    ft_topoplotER(cfg,grand_avgPre);
    %ft_topoplotER(cfg,grand_avg{1,1});% active_pre
    %subplot(2,length(time),t+length(time))
    %ft_topoplotER(cfg,grand_avg{1,2});% active_post
    %subplot(2,length(time),t+length(time))
     %ft_topoplotER(cfg,grand_avg{1,3});% sham_pre
     %subplot(2,length(time),t+length(time))
     %ft_topoplotER(cfg,grand_avg{1,4});% sham_post
end
%colormap(flipud(brewermap(64,'RdBu')))

% multichannel plot
cfg = [];
cfg.xlim = [-0.1,0.5];
cfg.showlabels  = 'yes';
cfg.layout      = 'acticap-64ch-standard2.mat';
figure;
%ft_multiplotER(cfg, grand_avgHC, grand_avgPre);
ft_multiplotER(cfg, grand_avg{1,1}, grand_avg{1,3}); % pre: active and sham
figure;
ft_multiplotER(cfg, grand_avg{1,1}, grand_avg{1,2}); % active: pre and post
figure;
ft_multiplotER(cfg, grand_avg{1,3}, grand_avg{1,4}); % sham: pre and post

% single channel plot
cfg = [];
cfg.xlim = [-0.1,0.5];
cfg.channel = {'F3','F5','F1','FC3','AF3'};
figure;
subplot(2,1,1)
ft_singleplotER(cfg, grand_avgHC, grand_avgPre,grand_avg{1,2});
legend({'hc','pre','active_post'})
subplot(2,1,2)
ft_singleplotER(cfg, grand_avg{1,1}, grand_avg{1,2},grand_avg{1,3}, grand_avg{1,4}); % active: pre and post
legend({'active_pre','active_post','sham_pre','sham_post'})


for nc=1:4
    figure;
    for t=1:length(time)
        cfg.xlim = time{1,t};
        cfg.zlim = [-3,3];
        subplot(1,length(time),t)
        ft_topoplotER(cfg,grand_avg{1,nc});
        title(tle{t});
    end
end

%% step 4: GMFP/TEPs difference and t-value plot
cfg = [];
cfg.method      = 'template';                         % try 'distance' as well              % specify type of template
cfg.layout      = 'acticap-64ch-standard2.mat';                % specify layout of channels
%cfg.feedback    = 'yes';                              % show a neighbour plot
neighbours      = ft_prepare_neighbours(cfg, pre_avg{1,1}); % define neighbouring channels

cfg = [];
cfg.channel     = {'all'};
cfg.neighbours  = neighbours; % defined as above
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.latency     = [-0.1 0.5];%'all';
%cfg.latency     = time{1,2};
%cfg.avgoverchan = 'yes'; % average over channel
%cfg.avgovertime = 'yes'; % average over time
% cfg.method      = 'montecarlo';
cfg.statistic   = 'indepsamplesT';
% %cfg.statistic   = 'depsamplesT';
% cfg.correctm    = 'cluster';
% cfg.correcttail = 'prob';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.minnbchan        = 2;
% cfg.tail             = 0;
% cfg.clustertail      = 0;
% cfg.numrandomization = 1000;
cfg.alpha            = 0.05;

%Ncond_1 = length(hc_avg);
%Ncond_2 = length(pre_avg);
Ncond_1 = length(active_avg(2,:));
Ncond_2 = length(sham_avg(2,:));
cfg.design(1,:)  = [ones(1,Ncond_1) 2*ones(1,Ncond_2)];
cfg.ivar         = 1; % the 1st row in cfg.design contains the independent variable

% subj=length(sham_avg);%length(active_TFR);
% design = zeros(2,2*subj);
% design(1,1:subj)=1;
% design(1,(subj+1):2*subj)=2;
% for i=1:subj
%     design(2,i)=i;
%     design(2,subj+i)=i;
% end
% cfg.ivar=1;
% cfg.uvar=2;
% cfg.design=design;


%stat = ft_timelockstatistics(cfg, HC_gmfp_sub{:}, pre_gmfp_sub{:});
%stat = ft_timelockstatistics(cfg, hc_avg{:}, pre_avg{:});
%stat = ft_timelockstatistics(cfg, active_avg{1,:}, active_avg{2,:});
stat = ft_timelockstatistics(cfg, active_avg{2,:}, sham_avg{2,:});
%stat = ft_timelockstatistics(cfg, active_gmfp_sub{1,:}, active_gmfp_sub{2,:});
%stat = ft_timelockstatistics(cfg, sham_gmfp_sub{1,:}, sham_gmfp_sub{2,:});
%stat = ft_timelockstatistics(cfg, active_contrast_avg{:}, sham_contrast_avg{:});
%stat = ft_timelockstatistics(cfg, active_contrast_gmfp{:}, sham_contrast_gmfp{:});
% 
cfg=[];
cfg.alpha = 0.05;
cfg.parameter = 'stat';
%cfg.channel = {'F3'};
cfg.xlim= [-0.1,0.5];
cfg.maskparameter = 'mask';
cfg.maskstyle ='box';
figure;
ft_singleplotER(cfg,stat);

% figure;
% cfg=[];
% cfg.parameter = 'stat';
% cfg.layout = 'acticap-64ch-standard2.mat';
% %cfg.xlim= [0.135,0.149];
% %cfg.zlim = [-2 2];
% cfg.showlabels   = 'yes';
% cfg.maskparameter = 'mask';
% cfg.maskstyle ='outline';
% %ft_topoplotER(cfg,stat); colorbar
% ft_clusterplot(cfg,stat); colorbar

%% step 5.1: GMFA/TEP analysis
% GMFA calculation for each subject
x=-0.1:0.001:0.5;
ind=1901:2501;
k=1;
HC_tmp=[];pre_tmp=[];
for i=1:length(hc_avg)
    cfg = [];
    cfg.method =  'amplitude';
    cfg.channel = {'F3','F5','F1','FC3','AF3'};
    HC_gmfp_sub{k} = ft_globalmeanfield(cfg, hc_avg{i});
    HC_tmp(k,:) = HC_gmfp_sub{k}.avg(ind);
    %HC_tmp(k,:) = squeeze(mean(hc_avg{1,k}.avg([2 33 35 36 38],ind),1));
    k=k+1;
end
k=1;
pre_tmp=[];
for i=1:length(pre_avg)
    cfg = [];
    cfg.method = 'amplitude';
    cfg.channel = {'F3','F5','F1','FC3','AF3'};
    pre_gmfp_sub{k} = ft_globalmeanfield(cfg, pre_avg{i});
    pre_tmp(k,:)=pre_gmfp_sub{k}.avg(ind);
    %pre_tmp(k,:) = squeeze(mean(pre_avg{1,k}.avg([2 33 35 36 38],ind),1));
    k=k+1;
end

%%%%%%% pre vs. post %%%%%
x=-0.1:0.001:0.5;
ind=1901:2501;
active_tmp=[];
sham_tmp=[];
for nc =1:2
    k=1;
    for i=1:length(active_avg)
        cfg = [];
        cfg.method = 'amplitude';
        cfg.channel = {'F3','F5','F1','FC3','AF3'};
        active_gmfp_sub{nc,k} = ft_globalmeanfield(cfg, active_avg{nc,i});
        active_tmp{nc}(k,:)=active_gmfp_sub{nc,k}.avg(ind);
        %active_tmp{nc}(k,:) = squeeze(mean(active_avg{nc,i}.avg([2 33 35 36 38],ind),1));
        k=k+1;
    end
    k=1;
    for i=1:length(sham_avg)
        cfg = [];
        cfg.method = 'amplitude';
        cfg.channel = {'F3','F5','F1','FC3','AF3'};
        sham_gmfp_sub{nc,k} = ft_globalmeanfield(cfg, sham_avg{nc,i});
        sham_tmp{nc}(k,:)=sham_gmfp_sub{nc,k}.avg(ind);
        %sham_tmp{nc}(k,:) = squeeze(mean(sham_avg{nc,i}.avg([2 33 35 36 38],ind),1));
        k=k+1;
    end
    
end

%% step 5.2: GMFA/TEP analysis and plot
% plot GMFP/TEPs with errorbar
A=active_tmp{2}-active_tmp{1};
A = A - mean(A(:,[1:100]),2);
%A=pre_tmp - mean(pre_tmp(:,[1:100]),2);
y1=mean(A,1);
errBar1=std(A)/sqrt(size(A,1));
uE1=y1+errBar1;
lE1=y1-errBar1;
yP1=[lE1,fliplr(uE1)];


B=sham_tmp{2}-sham_tmp{1};
B = B - mean(B(:,[1:100]),2);
%B=HC_tmp - mean(HC_tmp(:,[1:100]),2);
y2=mean(B,1);
errBar2=std(B)/sqrt(size(B,1));
uE2=y2+errBar2;
lE2=y2-errBar2;
yP2=[lE2,fliplr(uE2)];

x = -100:500;
xP= [x,fliplr(x)];
y=[y1;y2];
yy=NaN(length(x),1);
h1=ttest(A);h2=ttest2(A,B);
h3=h1+h2;
%ind1 = find(h3==2);
ind1 = find(x>=150 & x<=185);
%ind1 = find(x>=164 & x<=215);
yy(ind1)=4;
%yy(ind1)=y2(ind1);

col = {[0.75,0.08,0.13],[0,0,0],[0.95,0.49,0.44],[0.96,0.75,0.48]};
figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1,...
    'Position',[0.128257839721254 0.186679841897233 0.775 0.785750988142292]);
hold(axes1,'on');
area(x,yy,...
    'FaceColor',[0.9 0.9 0.9],...
    'EdgeColor',[1 1 1]);
patch('XData',xP,'YData',yP1,'FaceAlpha',0.2,'EdgeColor','none','FaceColor',col{1,3});hold on
plot(x,y1,'LineWidth',1,'Color',col{1,3});
patch('XData',xP,'YData',yP2,'FaceAlpha',0.2,'EdgeColor','none','FaceColor',col{1,4})
plot(x,y2,'LineWidth',1,'Color',col{1,4});
xlabel({'Time (ms)'});
%ylabel({'TEPs (uV)'});
ylabel({'LMFA (uV)'});
annotation(figure1,'textbox',...
    [0.71702787456446 0.784273166447079 0.078397212543554 0.0968379446640316],...
    'Color',col{1,3},...
    'String',{'Active'},...
    'FontSize',14,...
    'EdgeColor','none');
annotation(figure1,'textbox',...
    [0.71702787456446 0.724273166447079 0.078397212543554 0.0968379446640316],...
    'Color',col{1,4},...
    'String',{'Sham'},...
    'FontSize',14,...
    'EdgeColor','none');
%axis([-100,500,-4,4])
%axis([-100,500,-2,6])
axis([-100,500,-2,4])

%% step 5.3: time window selection - cluster analysis

%%% calculate cluster p values %%%
AA=HC_tmp;
BB=pre_tmp;
All = [AA;BB];
[H,P,CI,STATS]=ttest2(AA,BB);
x = -100:500;
tw_ind = x>=164 & x<=215; % identify time window where H ~=0

% AA=active_tmp{2};
% BB=active_tmp{1};
% All = [AA;BB];
% [H,P,CI,STATS]=ttest(AA,BB);
% tw_ind = x>=150 & x<=185;

real = sum(STATS.tstat(tw_ind));

shuff=[];
for i=1:1000
    ord = randperm(117);
    [H,P,CI,STATS]=ttest2(All(ord(1:64),:),All(ord(65:117),:));
    %ord = randperm(56);
    %[H,P,CI,STATS]=ttest(All(ord(1:28),:),All(ord(29:56),:));
    shuff(i)= sum(STATS.tstat(tw_ind));
end

% shuff=[];
% for i=1:1000
%     %ord = randperm(117);
%     %[H,P,CI,STATS]=ttest2(All(ord(1:53),:),All(ord(54:117),:));
%     ord = randperm(56);
%     [H,P,CI,STATS]=ttest(All(ord(1:28),:),All(ord(29:56),:));
%     %shuff(i)= sum(STATS.tstat(tw_ind));
%     ind = find(H~=0);
%     if length(ind)<=3 
%         shuff(i)=max(STATS.tstat);
%     else
%     c1 = 1;
%     arrset = cell(0,0);
%     while (c1<numel(ind))
%         c2 = 0;
%         while (c1+c2+1 <= numel(ind) && ind(c1)+c2+1==ind(c1+c2+1))
%             c2=c2+1;
%         end
%         if(c2>=1)
%             arrset = [arrset;(ind(c1:1:c1+c2))]; % 
%         end
%         c1=c1+c2+1;
%     end
%     
%     tsum=[];
%     for k=1:length(arrset)
%         temp = cell2num(arrset(k,1));
%         tsum(k)=sum(STATS.tstat(temp));
%     end
%     shuff(i)= max(tsum);
%     end
%             
% end


shuff =sort(shuff);
figure;
hist(shuff,15); hold on
plot(repmat(real,100),1:100)

%% step 6: extract output: LMFA - AUC,TEPs
%%%%% TEPS %%%%%%%%%%
x=-0.1:0.001:0.5;
time = {[0.025 0.035],[0.035,0.055],[0.040 0.070],[0.080 0.120],[0.16 0.20]};
TEP_mean=[];TEP_trt=[];
for k = 1:length(time)
    ind = x>=time{1,k}(1) & x<time{1,k}(2);
    TEP_mean.data(:,k) = [mean(HC_tmp(:,ind),2);mean(pre_tmp(:,ind),2)]; % HC+ MDD   
%     TEP_trt.data(:,k)=[mean(active_tmp{1,1}(:,ind),2);mean(active_tmp{1,2}(:,ind),2);...
%         mean(sham_tmp{1,1}(:,ind),2);mean(sham_tmp{1,2}(:,ind),2)];
    
end

%%%%% LMFA %%%%%
x=-0.1:0.001:0.5;
time ={[0.164 0.215]};
HC_AUC=[];pre_AUC=[];
for k = 1:length(time)
    index = x>=time{1,k}(1) & x<time{1,k}(2);
    for i = 1:size(HC_tmp,1)
        HC_AUC(i,k) = max(cumtrapz(x(index),HC_tmp(i,index))); % later
    end
    for i=1:size(pre_tmp,1)
        pre_AUC(i,k) = max(cumtrapz(x(index),pre_tmp(i,index))); % later
    end
    
end
[HC_AUC;pre_AUC];
[H0,P0,CI0] = ttest2(HC_AUC,pre_AUC)

x=-0.1:0.001:0.5;
active_AUC=[];sham_AUC=[];
time ={[0.150,0.185]};
for k = 1:length(time)
    index = x>=time{1,k}(1) & x<time{1,k}(2);
    for i = 1:size(active_tmp{1,1},1)
    active_AUC{1}(i,k) = max(cumtrapz(x(index),active_tmp{1}(i,index))); % pre
    active_AUC{2}(i,k) = max(cumtrapz(x(index),active_tmp{2}(i,index))); % post
    end
    for i=1:size(sham_tmp{1,1},1)
    sham_AUC{1}(i,k) = max(cumtrapz(x(index),sham_tmp{1}(i,index))); % pre
    sham_AUC{2}(i,k) = max(cumtrapz(x(index),sham_tmp{2}(i,index))); % post
    end
    
end

%[H0,P0,CI0] = ttest2(active_AUC{2}(:,2),sham_AUC{2}(:,2))
[active_AUC{1}(:,2);active_AUC{2}(:,2);sham_AUC{1}(:,2);sham_AUC{2}(:,2)];
