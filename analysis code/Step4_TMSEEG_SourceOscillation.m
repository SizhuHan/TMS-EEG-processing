clear;
clc;
sub_no = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15'...
    '16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'};
subpath = {'active/pre/','active/post/','sham/pre/','sham/post/'};
cond = {'active_pre','active_post','sham_pre','sham_post','HC'};
%ROI_Name = {'Frontal_Mid_L','Frontal_Sup_L'}; %ROI_MNI_V4
ROI_Name = {{'Frontal_Mid_2_L'},{'Amygdala_L','Amygdala_R'},...
    {'ACC_sub_L','ACC_sub_R'},{'Hippocampus_L','Hippocampus_R'}...
    {'OFCant_L','OFCant_R'},{'Parietal_Sup_L','Parietal_Sup_R'},...
    {'Precentral_L','Precentral_R'}}; %ROI_MNI_V7

roi_id = 6;

aal   = ft_read_atlas('/fieldtrip-20161231/template/atlas/aal/ROI_MNI_V7.nii');
aal.coordsys = 'mni';
aal.tissue = aal.brick0;
aal.tissuelabel = aal.brick0label;
load('/spm12/toolbox/AAL3/ROI_MNI_V7_List.mat');% label name
for i=1:length(ROI)
    aal.tissuelabel{ROI(i).ID} = ROI(i).Nom_L;
end
aal = rmfield(aal,{'brick0label','brick0'});

% find coordinates of ROIs
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
band = {[1,4],[4,8],[8,13],[13,30],[30,50]};
fc = 1000; % sample rate = 1000Hz
N = 500;   
n = 0:N-1; %
f = n*fc/N; % frequency sequence
for nc=1:4
    ROI_freqSCD=[];ROI_freqSCS=[];
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
        
        for nf=1:5
           Wn = band{nf}*2/fc; %??
           [m,l] = butter(2,Wn);%4?IIR???
           result = filtfilt(m,l,post_origin');
           tmp = hilbert(result);
           trans_data= abs(tmp);
           TimeCourse = trans_data'.* Mask;  % voxcel * time
        
        for t = 1:length(AvgPostSource.time)
            [ns, t]
           % AvgPostSource.avg.pow(AvgPostSource.inside) = TimeCourse(:,t);            
           % ROI_freqSCD{nf}(k,t) = nansum(AvgPostSource.avg.pow(roi_mask));
            
            AvgPostSource.avg.pow(AvgPostSource.inside) = Mask(:,t);
            ROI_freqSCS{nf}(k,t) = nansum(AvgPostSource.avg.pow(roi_mask).* distance);
                                   
        end
        end
        k=k+1;
    end
    
    switch nc
        case 1
           % Active_pre_freqSCD = ROI_freqSCD;
            Active_pre_freqSCS = ROI_freqSCS;
        case 2
           % Active_post_freqSCD = ROI_freqSCD;
            Active_post_freqSCS = ROI_freqSCS;
        case 3
           % Sham_pre_freqSCD = ROI_freqSCD;
            Sham_pre_freqSCS = ROI_freqSCS;
        case 4
           % Sham_post_freqSCD = ROI_freqSCD;
            Sham_post_freqSCS = ROI_freqSCS;
        case 5
           % HC_freqSCD = ROI_freqSCD;
            HC_freqSCS = ROI_freqSCS;
    end
end
cd /Volumes/SizhuFiles/TMS+EEGstudy/Source2_Data
%save(['Frontal_Mid_2_L_freq.mat'],'Active*','Sham*','HC*'); % 73
%save(['OFCant_LR_freq.mat'],'Active*','Sham*','HC*'); % 19
%save(['Hippocampus_LR_freq.mat'],'Active*','Sham*','HC*'); % 30
save(['Parietal_Sup_LR_freq.mat'],'Active*','Sham*','HC*'); % 66


%% SCD/SCS dynamic plot
clear A B
col = {[0.89,0.15,0],[0.13,0.13,0.13],[0.06,0.56,0.25],[0.52,0.52,0.52]};
n=113;
for i=1:500-9
    A(:,i) = mean(Sham_pre_freqSCD{1}(:,i:i+9)/n,2);
    B(:,i) = mean(Sham_post_freqSCD{1}(:,i:i+9)/n,2);
end

y1=mean(A,1);
errBar1=std(A)/sqrt(size(A,1));
uE1=y1+errBar1;
lE1=y1-errBar1;
yP1=[lE1,fliplr(uE1)];

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
plot(x,y1,'LineWidth',1,'Color',col{1,1}); hold on
patch('XData',xP,'YData',yP1,'FaceAlpha',0.2,'EdgeColor','none','FaceColor',col{1,1})
plot(x,y2,'LineWidth',1,'Color',col{1,2});
patch('XData',xP,'YData',yP2,'FaceAlpha',0.2,'EdgeColor','none','FaceColor',col{1,2})
xlabel({'Time (ms)'});
%ylabel({'Current Density, uA/mm^2'});
ylabel({'SCD, uA/mm^2'});
%ylabel({'SCS, mm'});
annotation(figure1,'textbox',...
    [0.71702787456446 0.784273166447079 0.078397212543554 0.0968379446640316],...
    'Color',col{1,1},...
    'String',{'MDD'},...
    'FontSize',12,...
    'EdgeColor','none');
annotation(figure1,'textbox',...
    [0.71702787456446 0.724273166447079 0.078397212543554 0.0968379446640316],...
    'Color',col{1,2},...
    'String',{'HC'},...
    'FontSize',12,...
    'EdgeColor','none');
%axis([0,500,0,5e-4])
%axis([0,500,0,50])

%% extract data
n=6;
for foi = 1:5
    Pre_Data = [Active_pre_freqSCD{foi};Sham_pre_freqSCD{foi}];
    HC_Data = HC_freqSCD{foi};
    Active_Data = [Active_pre_freqSCD{foi};Active_post_freqSCD{foi}];
    Sham_Data = [Sham_pre_freqSCD{foi};Sham_post_freqSCD{foi}];
for i=1:500-9
    A(:,i) = mean(Pre_Data(:,i:i+9)/n,2);
    B(:,i) = mean(HC_Data(:,i:i+9)/n,2);
    C(:,i) = mean(Active_Data(:,i:i+9)/n,2);
    D(:,i) = mean(Sham_Data(:,i:i+9)/n,2);
end
x=5:495;
toi_1 = find(x >= 164 & x <= 215); %
toi_2 = find(x >= 150 & x <= 185); %  

Pre_Avg = mean(A(:,toi_1),2);
HC_Avg = mean(B(:,toi_1),2);
Active_Avg = mean(C(:,toi_2),2);
Sham_Avg = mean(D(:,toi_2),2);

HCvsPats(:,foi) = [HC_Avg;Pre_Avg];
PrevsPost(:,foi) = [Active_Avg;Sham_Avg];
end

active_subtract = PrevsPost([29:56],:) - PrevsPost([1:28],:);
sham_subtract = PrevsPost([82:106],:) - PrevsPost([57:81],:);

