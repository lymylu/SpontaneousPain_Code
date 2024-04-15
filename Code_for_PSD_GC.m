% Analysis LFP from ExtractData/ using NeuroView 1.3.0
% analysis three type of activities: Spontaneouspain, Groom, Move..
% no active means NoSpontaneouspain, NoGroom, NoMove..
Basepath='ExtractData';
eventtype={'Spontaneousnopain','Spontaneouspain'};
filelist=dir(Basepath);
PSDparams=PowerSpectralDensity.getParams();
GCparams=getGrangerCausualityParams();
for i=3:length(filelist)
    tmpmat=matfile(fullfile(filelist(i).folder,filelist(i).name));
    tmpmatsave=matfile(fullfile('PSDresult',filelist(i).name),"Writable",true);
    try 
        tmpdata=eval(['tmpmat.',eventtype{2}]);
    for j=1:2
        tmpdata=eval(['tmpmat.',eventtype{j},';']);
        tmpdata=NeuroResult(tmpdata);
        tmpdata=PowerSpectralDensity.recal(PSDparams,tmpdata,'PSD');  
        tmpdata=GrangerCausaulity(GCparams,tmpdata,'GC');
        eval(['tmpmatsave.',eventtype{j},'=tmpdata;']);
    end
    end
end
%% summarize the PSD and GC for Figure 1 and 
Resultpath='PSDresult';
filelist=dir(Resultpath);
for j=1:2
    c=1;
    for i=3:length(filelist)
        if contains(filelist(i).name,{'_1d_','_3d_','_7d_','_9d_','_14d_'})
            tmpmat=matfile(fullfile('PSDresult',filelist(i).name),"Writable",true);
        try
        if j==2
            Active(c)=NeuroResult(eval(['tmpmat.',eventtype{j},';']));
            c=c+1;
        else
            NoActive(c)=NeuroResult(eval(['tmpmat.',eventtype{j},';']));
            c=c+1;
        end
        end
        end
    end
end
Brainname=unique(Active(1).LFPinfo.channeldescription);
%%
Brainname={'PL','IL','NACcore','NACshell','S1','MD','dCA1','BLA','CeA','vCA1'};
channeldescription=Active(1).LFPinfo.channeldescription;

for j=1:length(Active)
    ActivePSD(:,:,j)=Active(j).PSD.S{:};
    NoActivePSD(:,:,j)=NoActive(j).PSD.S{:};
    ActiveGC(:,:,j)=Active(j).GC.Granger;
    NoActiveGC(:,:,j)=NoActive(j).GC.Granger;
end
NoActivePSDnew=[];ActivePSDnew=[];
figure; % Figure S1B
for i=1:length(Brainname)
        NoActivePSDnew(:,i,:)=mean(NoActivePSD(:,ismember(channeldescription,Brainname{i}),:),2);
        ActivePSDnew(:,i,:)=mean(ActivePSD(:,ismember(channeldescription,Brainname{i}),:),2);
end
for i=1:length(Brainname)
    valid=true(size(ActivePSDnew,3),1);
[~,invalidindex1{i}]=deleteoutliers(squeeze(mean(ActivePSDnew(:,i,:),1)),0.01);
[~,invalidindex2{i}]=deleteoutliers(squeeze(mean(NoActivePSDnew(:,i,:),1)),0.01);
subplot(2,5,i); 
plot(Active(1).PSD.f_lfp,squeeze(mean(ActivePSDnew(:,i,valid),3)),'blue');
hold on;
plot(Active(1).PSD.f_lfp,squeeze(mean(NoActivePSDnew(:,i,valid),3)),'red');
shadebar(Active(1).PSD.f_lfp,squeeze(ActivePSDnew(:,i,valid)),'blue');
shadebar(Active(1).PSD.f_lfp,squeeze(NoActivePSDnew(:,i,valid)),'red');
set(gca,'YScale','log');
title(Brainname{i}); 
end
% specific plot for PL Figure 1F
figure;hold on;
imagesc(Active(1).PSD.f_lfp,1:10,repmat(PSD_pval_correct(:,1)',[10,1]));caxis([0,0.05]);colormap gray;
% plot(Active(1).PSD.f_lfp,squeeze(mean(ActivePSDnew(:,7,valid),3)),'blue');
% plot(Active(1).PSD.f_lfp,squeeze(mean(NoActivePSDnew(:,7,valid),3)),'red');
shadebar(Active(1).PSD.f_lfp,squeeze(ActivePSDnew(:,1,valid)),'blue');
shadebar(Active(1).PSD.f_lfp,squeeze(NoActivePSDnew(:,1,valid)),'red');
set(gca,'YScale','log');axis tight;

for i=1:size(ActivePSDnew,1)
    for j=1:size(ActivePSDnew,2)
        [~,PSD_pval(i,j),~,stat]=ttest(squeeze(ActivePSDnew(i,j,:)),squeeze(NoActivePSDnew(i,j,:)));
        PSD_Tstat(i,j)=stat.tstat;
    end
end
PSD_pval_correct=reshape(fdr_BH(PSD_pval,0.05),164,10);
% Figure 1E
figure; subplot(1,2,1); 
imagesc(Active(1).PSD.f_lfp,1:10,PSD_pval_correct'); caxis([0,0.05]); yticklabels(Brainname); axis xy;
subplot(1,2,2); % 
imagesc(Active(1).PSD.f_lfp,1:10,PSD_Tstat');yticklabels(Brainname);axis xy;
%%
Brainname=unique(Brainname);
Brainname={'PL','IL','NACcore','NACshell','S1','MD','dCA1','BLA','CeA','vCA1'};
for i=1:length(Brainname)
    for j=1:length(Brainname)
        NoActiveGCnew(i,j,:)=mean(mean(NoActiveGC(ismember(channeldescription,Brainname{i}),ismember(channeldescription,Brainname{j}),:),1),2);
        ActiveGCnew(i,j,:)=mean(mean(ActiveGC(ismember(channeldescription,Brainname{i}),ismember(channeldescription,Brainname{j}),:),1),2);
    end
end
for i=1:length(Brainname)
    for j=1:length(Brainname)
        valid=true(size(ActiveGC,3),1);
        try
        [~,T_pval(i,j),~,tstat]=ttest(ActiveGCnew(i,j,valid),NoActiveGCnew(i,j,valid));
        Tstat(i,j)=tstat.tstat;
        catch
            T_pval(i,j)=nan;Tstat(i,j)=[];
        end
    end
end
T_pval_correct= reshape(fdr_BH(T_pval,0.05),10,10);
% Figure 1C
figure; subplot(1,2,2);plot_pw(T_pval_correct);xticks(1:10);yticks(1:10);xticklabels(Brainname);yticklabels(Brainname); clim([0,0.05]);colormap gray;colorbar;
subplot(1,2,1); plot_pw(Tstat);xticks(1:10);yticks(1:10);xticklabels(Brainname);yticklabels(Brainname);clim([-3,3]); colormap bone;colorbar;
function tmpdata=GrangerCausaulity(GCparams,tmpdata,GCname)
    % cal the GC value from the NeuroResult data format
    tmpdata=tmpdata.Split2Splice();
    for i=1:size(tmpdata.LFPdata{1},2)
        tmpdata.LFPdata{1}(:,i)=detrend(tmpdata.LFPdata{1}(:,i));
        tmpdata.LFPdata{1}(:,i)=zscore(tmpdata.LFPdata{1}(:,i));
    end
    [epochtime]=windowepoched(tmpdata.LFPdata{1},[GCparams.epochsize,GCparams.epochsize],1,inf,tmpdata.LFPinfo.Fs);
    tmplfp=[];
    for i=1:size(epochtime,2)
        tmplfp(:,:,i)=downsample(tmpdata.LFPdata{1}(epochtime(:,i),:),5);
    end 
    tmplfp=basecorrect(tmplfp,1:size(tmplfp,1),1,size(tmplfp,1),'Zscore');
    tmplfp=permute(tmplfp,[2,1,3]);
    channelname=unique(tmpdata.LFPinfo.channeldescription);
    [AIC,BIC] = tsdata_to_infocrit(tmplfp,GCparams.maxorder,GCparams.inforegress);
    [~,bmo_AIC] = min(AIC);
    [F,A,SIG] = GCCA_tsdata_to_pwcgc(tmplfp,bmo_AIC,GCparams.varregress);
    tmpdata.addprop(GCname);
    eval(['tmpdata.',GCname,'.Granger=F;']);
    eval(['tmpdata.',GCname,'.sources=channelname;']);
end
function params=getGrangerCausualityParams()
 prompt={'max model order','information criteria regression mode','VAR model estimation regression mode','epoched window size'};
 title='params';
lines=4;
def={'20','LWR','OLS','1'};
x=inputdlg(prompt,title,lines,def,'on');
params.maxorder=str2num(x{1});
params.inforegress=x{2};
params.varregress=x{3};
params.epochsize=str2num(x{4});
end