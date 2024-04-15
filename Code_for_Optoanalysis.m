% opto tagged neurons modulated the stimulus-induced neurons
BaseDir='OptoinducedFiles';
filelist=dir(BaseDir);
Eventtype={'VonfreyInON','VonfreyUnON','InjuredON','UninjuredON','VonfreyInOFF','VonfreyUnOFF','InjuredOFF','UninjuredOFF','OptoTag','airpuffON','airpuffOFF'};
for i=3:length(filelist)
      tmpmat=matfile(fullfile(BaseDir,filelist(i).name),'Writable',true);
      tmpmat=cal_binspikes(tmpmat,50,'Combined',Eventtype); % binspike_short as 50Hz bin; binspike 10Hz bin
      tmpmat=cal_optotag(tmpmat);
      tmpmat=cal_response(tmpmat,Eventtype);
      tmpmat=cal_spectrogram(tmpmat);
      tmpmat=cal_binspikes(tmpmat,1,'Spon','SPKdata'); % long time stimulus: 1Hz bin
end
%% plot the stimulus induced responses of different event types summarized from the opto tagged and opto untagged neurons.
optospike=cell(1,11); unoptospike=cell(1,11); optopsth=cell(1,11); unoptopsth=cell(1,11); t=linspace(-2,4,61);
binspikesummary=cell(1,11);
eventresponse=cell(1,11);
psthspikesummary=cell(1,11);
optoindex_all=cell(1,11);
neuroduration=cell(1,11);
binspikesummary_short=cell(1,11);
neurofiringrate=cell(1,11);
spikename=cell(1,1);
 Blacklist=matfile('optoinduced/Blacklist.mat');% 
 for j=1:11
     for i=3:length(filelist)
    tmpmat=matfile(fullfile(BaseDir,filelist(i).name)); 
    tmpmat=tmpmat.Combined;
     eventinvalid=eval(['Blacklist.',filelist(i).name(1:end-4)]);
     eventinvalid=cellfun(@(x) str2num(x),eventinvalid.Eventindex,'UniformOutput',1);
     eventselect=tmpmat.EVTinfo.eventselect;
     eventinvalid=ismember(eventselect,eventinvalid); 
    tmpoptoindex=tmpmat.Optoresponse==1;
    eventindex=ismember(tmpmat.EVTinfo.eventdescription,Eventtype{j});
    tmpbinspike=cell2mat(cellfun(@(x) nanmean(x(:,eventindex),2),tmpmat.binspike,'UniformOutput',0));
    tmpbinspike_short=cell2mat(cellfun(@(x) nanmean(x(:,eventindex),2),tmpmat.binspike_short,'UniformOutput',0));
    tmppsth=eval(['tmpmat.',Eventtype{j},'psth']);
    tmpneurontype=tmpmat.SPKinfo.troughToPeak;
    tmpneuronfiringrate=tmpmat.SPKinfo.firingRate;
    if j==1
    tmpspikename=cellfun(@(x) [tmpmat.Subjectname,x],tmpmat.SPKinfo.name,'UniformOutput',0);
    spikename=cat(2,spikename,tmpspikename);
    end
    tmppsthspike=eval(['tmpmat.',Eventtype{j},'psth']);
    binspikesummary{j}=cat(2,binspikesummary{j},tmpbinspike);
    psthspikesummary{j}=cat(2,psthspikesummary{j},tmppsth);
    tmpeventresponse=eval(['tmpmat.',Eventtype{j},'response']);  
    eventresponse{j}=cat(2,eventresponse{j},tmpeventresponse);
    optoindex_all{j}=cat(2,optoindex_all{j},tmpoptoindex);
    neuroduration{j}=cat(2,neuroduration{j},tmpneurontype);
    neurofiringrate{j}=cat(2,neurofiringrate{j},tmpneuronfiringrate);
    binspikesummary_short{j}=cat(2,binspikesummary_short{j},tmpbinspike_short);
     end
 end

%% delet the invalid neurons.
optoindex=optoindex_all{1};neuroduration=neuroduration{1};neurofiringrate=neurofiringrate{1};
neuroinvalid=neuroduration>0.5&neurofiringrate>0.5&neurofiringrate<10; % choose the pyramidal neuron to calculate.
%%
binspikesummary_spon=[];
for i=3:length(filelist)
    tmpmat=matfile(fullfile(SponDir,filelist(i).name));
    tmpmat=tmpmat.Spon;
    tmpbinspike=cell2mat(tmpmat.binspike);
    binspikesummary_spon=cat(2,binspikesummary_spon,tmpbinspike);
end
%%
t_short=linspace(-2,4,301);
% calculate negative (irrevalent group)
negative=negativeresponse(binspikesummary_short_correct{9},t_short);
negative(optoindex==1&negative==1)=0;
t_spon=linspace(-3,6,541);
%%
binspikesummary_correct=cellfun(@(x) basecorrect(x,t,-1,0,'zscore'),binspikesummary,'UniformOutput',0);
binspikesummary_short_correct=cellfun(@(x) basecorrect(x,t_short,-1,0,'zscore'),binspikesummary_short,'UniformOutput',0);
binspikesummary_spon_correct=basecorrect(binspikesummary_spon,t_spon,-3,0,'Zscore');
invalid=cell2mat(cellfun(@(x) isnan(x(1,:))',binspikesummary_correct,'UniformOutput',0));
invalid=logical(mean(invalid,2));
neuroinvalid=neuroinvalid&~invalid';
%% plot the opto neurons with three groups 
close all;
% load Optoanalysisworkspace.mat
for i=1:3
    figure; 
    if i==3
        sgtitle('Optotagged neurons');
        tmpindex=optoindex==1&neuroinvalid==1;
        tmpindex=find(tmpindex==1);
       [~,sortindex]=sort(mean(binspikesummary_short_correct{9}(t_short>=0&t_short<0.04,tmpindex),1),2);
       tmp=binspikesummary_short_correct{9}(:,tmpindex);
       imagesc(t,1:size(tmp,2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight; xlim([-0.5,1]);caxis([-1,2]);
       typeindex{i}=tmpindex(sortindex);
    elseif i==2
        sgtitle('Opto Modulated neurons');
        tmpindex=optoindex==0&negative==1&neuroinvalid==1;
        tmpindex=find(tmpindex==1);
        [~,sortindex]=sort(mean(binspikesummary_short_correct{9}(t_short>=0&t_short<0.1,tmpindex),1),2);
       tmp=binspikesummary_short_correct{9}(:,tmpindex);
       imagesc(t,1:size(tmp,2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight; xlim([-0.5,1]);caxis([-1,2]);
       typeindex{i}=tmpindex(sortindex);

    elseif i==1
        sgtitle('Opto Unmodulated neurons');
        tmpindex=optoindex==0&negative==0&neuroinvalid==1;
        tmpindex=find(tmpindex==1);
        [~,sortindex]=sort(mean(binspikesummary_short_correct{9}(t_short>=0&t_short<0.1,tmpindex),1),2);
       tmp=binspikesummary_short_correct{9}(:,tmpindex);
       imagesc(t,1:size(tmp,2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight; xlim([-0.5,1]);caxis([-1,2]);
        typeindex{i}=tmpindex(sortindex);
    end
end
%%
% something wrong in the Salt? the opto tag neurons should be modified
tmp_typeindex1=cat(2,typeindex{3}(7:end),typeindex{1}(end-4:end));
tmp_typeindex2=cat(2,typeindex{3}(1:6),typeindex{1}(1:end-5));
typeindex{3}=unique(tmp_typeindex1);
typeindex{1}=unique(tmp_typeindex2);
typeindex{2}=unique(typeindex{2});
%% plot the corrected classfication of optoneurons
% Figure 4G
for i=1:3
    figure;
    if i==3
        sgtitle('Optotagged neurons');
       [~,sortindex]=sort(mean(binspikesummary_short_correct{9}(t_short>=0&t_short<0.04,typeindex{3}),1),2);
       tmp=binspikesummary_short_correct{9}(:,typeindex{i});
%        valid=[1:6,8:11,13:18,20:size(tmp,2)];
%        valid=1:size(tmp,2);
       imagesc(t,1:size(tmp(:,valid),2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight; xlim([-0.5,1]);caxis([-1,2]);
       typeindex_sort{i}=typeindex{i}(sortindex(valid));

    elseif i==2
        sgtitle('Opto Modulated neurons');
        [~,sortindex]=sort(mean(binspikesummary_short_correct{9}(t_short>=0&t_short<0.1,typeindex{2}),1),2);
       tmp=binspikesummary_short_correct{9}(:,typeindex{i});
       imagesc(t,1:size(tmp,2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight; xlim([-0.5,1]);caxis([-1,2]);
       typeindex_sort{i}=typeindex{i}(sortindex);

    elseif i==1
        sgtitle('Opto Unmodulated neurons');
        [~,sortindex]=sort(mean(binspikesummary_short_correct{9}(t_short>=0&t_short<0.1,typeindex{1}),1),2);
       tmp=binspikesummary_short_correct{9}(:,typeindex{i});
       imagesc(t,1:size(tmp,2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight; xlim([-0.5,1]);caxis([-1,2]);
        typeindex_sort{i}=typeindex{i}(sortindex);
    end
end
% Figure 4F
figure;hold on;
linecolor={'r','g','b'};
for i=1:3
    shadebar(t_short,binspikesummary_short_correct{9}(:,typeindex{i}),linecolor{i});
end
xlim([-0.2,0.5]);
%% for review, calculate the peak latency for each group of spike
for i=2:3
    tmpspike=smoothdata(binspikesummary_short_correct{9}(t_short>-0.02&t_short<0.1,typeindex{i}),1,'gaussian',5);
    if i==2
    [~,modulatelat]=min(tmpspike);
    modulatelat=modulatelat/50;
    elseif i==3
        [~,taglat]=max(tmpspike);
        taglat=taglat/50;
    end
end

%% Figure 4J
t=linspace(-2,4,61);
figtitle={'Opto Unmodulated neurons','Opto Modulated neurons','Opto tagged neurons'};
for i=1:3
    figure;
sgtitle(figtitle{i});
 %    tmpindex=invalidindex_modified{i};
subplot(2,4,1);
plot(t,smoothdata(nanmean(binspikesummary_correct{3}(:,typeindex_sort{i}),2),1,'gaussian',5));title('Light ON Injured side');
xlim([-0.5,1]); ylim([-1,3]);
%plot(pstht,unoptopsth{3});title('Light ON Injured side');xlim([-0.5,1]); %ylim([-0.5,3.5]);
subplot(2,4,2);
plot(t,smoothdata(nanmean(binspikesummary_correct{7}(:,typeindex_sort{i}),2),1,'gaussian',5)); title('Light OFF Injured side'); xlim([-0.5,1]); ylim([-1,3]);
%plot(pstht,unoptopsth{7});title('Light OFF Injured side');xlim([-0.5,1]); %ylim([-0.5,3.5]);
subplot(2,4,3);
plot(t,smoothdata(nanmean(binspikesummary_correct{4}(:,typeindex_sort{i}),2),1,'gaussian',5)); title('Light ON UnInjured side');xlim([-0.5,1]); ylim([-1,3]);
%plot(pstht,unoptopsth{4});title('Light ON UnInjured side');xlim([-0.5,1]); %ylim([-0.5,3.5]);
subplot(2,4,4);
plot(t,smoothdata(nanmean(binspikesummary_correct{8}(:,typeindex_sort{i}),2),1,'gaussian',5)); title('Light OFF UnInjured side');xlim([-0.5,1]); ylim([-1,3]);
%plot(pstht,unoptopsth{8});title('Light OFF UnInjured side');xlim([-0.5,1]); %ylim([-0.5,3.5]);
subplot(2,4,5);
[~,sortindex]=sort(mean(binspikesummary_correct{3}(t>0&t<1,typeindex_sort{i}),1),2);
tmp=binspikesummary_correct{3}(:,typeindex_sort{i});
imagesc(t,1:size(tmp,2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight; xlim([-1,2]);caxis([-1,3]);
subplot(2,4,6);
[~,sortindex]=sort(mean(binspikesummary_correct{7}(t>0&t<1,typeindex_sort{i}),1),2);
tmp=binspikesummary_correct{7}(:,typeindex_sort{i});
imagesc(t,1:size(tmp,2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight; xlim([-1,2]);caxis([-1,3]);
subplot(2,4,7);
[~,sortindex]=sort(mean(binspikesummary_correct{4}(t>0&t<1,typeindex_sort{i}),1),2);
tmp=binspikesummary_correct{4}(:,typeindex_sort{i});
imagesc(t,1:size(tmp,2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight; xlim([-1,2]);caxis([-1,3]);
subplot(2,4,8);
[~,sortindex]=sort(mean(binspikesummary_correct{8}(t>0&t<1,typeindex_sort{i}),1),2);
tmp=binspikesummary_correct{8}(:,typeindex_sort{i});
imagesc(t,1:size(tmp,2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight; xlim([-1,2]);caxis([-1,3]);
end
%%
optotaggedlongfiring_base=([nanmean(binspikesummary_spon(t_spon<0,typeindex_sort{3}),1);nanmean(binspikesummary_spon(t_spon>0&t_spon<3,typeindex_sort{3}),1);nanmean(binspikesummary_spon(t_spon>3,typeindex_sort{3}),1)])';
optomodulatedlongfiring_base=([nanmean(binspikesummary_spon(t_spon<0,typeindex_sort{2}),1);nanmean(binspikesummary_spon(t_spon>0&t_spon<3,typeindex_sort{2}),1);nanmean(binspikesummary_spon(t_spon>3,typeindex_sort{2}),1)])';
optounmodulatedlongfiring_base=([nanmean(binspikesummary_spon(t_spon<0,typeindex_sort{1}),1);nanmean(binspikesummary_spon(t_spon>0&t_spon<3,typeindex_sort{1}),1);nanmean(binspikesummary_spon(t_spon>3,typeindex_sort{1}),1)])';


%% bar plot the injured unoptospike and uninjured unoptospike in ON and OFF
optotaggedfiring=cellfun(@(x) nanmean(x(t>0&t<0.5,typeindex_sort{3}),1)/0.5,binspikesummary_correct([1:8,10,11]),'UniformOutput',0);
optomodulatedfiring=cellfun(@(x) nanmean(x(t>0&t<0.5,typeindex_sort{2}),1)/0.5,binspikesummary_correct([1:8,10,11]),'UniformOutput',0);
optounmodulatedfiring=cellfun(@(x) nanmean(x(t>0&t<0.5,typeindex_sort{1}),1)/0.5,binspikesummary_correct([1:8,10,11]),'UniformOutput',0);
optotaggedfiring=cell2mat(optotaggedfiring([8,4,7,3])')';
optomodulatedfiring=cell2mat(optomodulatedfiring([8,4,7,3])')';
optounmodulatedfiring=cell2mat(optounmodulatedfiring([8,4,7,3])')';
optoinducedfiring_base=cellfun(@(x) nansum(x(t>-1&t<0,typeindex_sort{3}),1)/1,binspikesummary([1:8,10,11]),'UniformOutput',0);
optomodulatedfiring_base=cellfun(@(x) nansum(x(t>-1&t<0,typeindex_sort{2}),1)/1,binspikesummary([1:8,10,11]),'UniformOutput',0);
optounmodulatedfiring_base=cellfun(@(x) nansum(x(t>-1&t<0,typeindex_sort{1}),1)/1,binspikesummary([1:8,10,11]),'UniformOutput',0);
optoinducedfiring_base=cell2mat(optoinducedfiring_base')';
optomodulatedfiring_base=cell2mat(optomodulatedfiring_base')';
optounmodulatedfiring_base=cell2mat(optounmodulatedfiring_base')';

%%
% Figure G-H
figure;t_spon=linspace(-3,6,541);
subplot(2,3,1);
tmp=binspikesummary_short_correct{9}(:,typeindex_sort{3});
[~,sortindex]=sort(mean(tmp(t_short>=0&t_short<=0.2,:),1),2);
%invalidindex{3}=invalidindex{3}(sortindex); % wrong classified the first 11; 
imagesc(t_short,1:size(typeindex_sort{3},2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight; xlim([-0.5,1]);caxis([-1,2]);title('Optotagged');
subplot(2,3,4);
tmp=binspikesummary_spon_correct(:,typeindex_sort{3});
imagesc(t_spon,1:size(typeindex_sort{3},2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight;title('Optotagged long');caxis([-1,2]);
subplot(2,3,2);
tmp=binspikesummary_short_correct{9}(:,typeindex_sort{2});
[~,sortindex]=sort(mean(tmp(t_short>=0&t_short<=0.2,:),1),2);
%invalidindex{2}=invalidindex{2}(sortindex);% wrong classified the last 6
imagesc(t_short,1:size(typeindex_sort{2},2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight; xlim([-0.5,1]);caxis([-1,2]);title('Optomodulated');
subplot(2,3,5);
tmp=binspikesummary_spon_correct(:,typeindex_sort{2});
imagesc(t_spon,1:size(typeindex_sort{2},2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight;title('Optomodulated long');caxis([-1,2]);
subplot(2,3,3);
tmp=binspikesummary_short_correct{9}(:,typeindex_sort{1});
[~,sortindex]=sort(mean(tmp(t_short>=0&t_short<=0.2,:),1),2);
%invalidindex{1}=invalidindex{1}(sortindex); % why? wrong classified the last 5 (opto) and first 14 (modulate)
imagesc(t_short,1:size(typeindex_sort{1},2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight; xlim([-0.5,1]);caxis([-1,2]);title('Optounmodulated');
subplot(2,3,6);
tmp=binspikesummary_spon_correct(:,typeindex_sort{1});
imagesc(t_spon,1:size(typeindex_sort{1},2),smoothdata(tmp(:,sortindex),1,'gaussian',5)');axis tight;title('Optounmodulated long');caxis([-1,2]);
%%
%% show the spectrogram during the laser stimulus and opto stimulus
lfp=cell(1,9);lfpt=linspace(-2,4,7501);
for i=1:9
    for j=3:length(filelist)
        tmpmat=matfile(fullfile(BaseDir,filelist(j).name));
        tmpmat=tmpmat.Combined;
        %eventinvalid=eval(['Blacklist.',filelist(j).name(1:end-4)]);
        %eventinvalid=cellfun(@(x) str2num(x),eventinvalid.Eventindex,'UniformOutput',1);
        %eventselect=tmpmat.EVTinfo.eventselect;
        %eventinvalid=ismember(eventselect,eventinvalid); 
        %eventindex=ismember(tmpmat.EVTinfo.eventdescription,Eventtype{i})&eventinvalid==0;
        lfp{i}=cat(3,lfp{i},nanmean(tmpmat.Spectrom_STFT(:,:,:),3));
    end
end
lfp_correct=cellfun(@(x) basecorrect(x,lfpt,-1,0,'Changepercent'),lfp,'UniformOutput',0);
% Figure L-M
figure;
    subplot(1,4,1);
    imagesc(lfpt,1:100,nanmean(lfp_correct{3},3)');title('Light ON Injured side');xlim([-0.5,1]); caxis([-0.2,1]); axis xy 
    subplot(1,4,2);
    imagesc(lfpt,1:100,nanmean(lfp_correct{7},3)');title('Light OFF Injured side');xlim([-0.5,1]);caxis([-0.2,1]);axis xy
    subplot(1,4,3);
    imagesc(lfpt,1:100,nanmean(lfp_correct{4},3)');title('Light ON Uninjured side');xlim([-0.5,1]);caxis([-0.2,1]);axis xy
    subplot(1,4,4);
    imagesc(lfpt,1:100,nanmean(lfp_correct{8},3)');title('Light OFF Uninjured side');xlim([-0.5,1]);caxis([-0.2,1]);axis xy
    
    figure;
   
    imagesc(lfpt,1:100,nanmean(lfp{9},3)');title('Opto induced spectrogram');xlim([-0.5,1]); axis xy
   gamma=cellfun(@(x) squeeze(nanmean(nanmean(x(:,30:100,:),2),3)),lfp_correct,'UniformOutput',0);
   figure;hold on;
   plot(lfpt,gamma{3});xlim([-0.5,1]); ylim([-0.5,1.5]);
   plot(lfpt,gamma{7});xlim([-0.5,1]); ylim([-0.5,1.5]);
   plot(lfpt,gamma{4});xlim([-0.5,1]); ylim([-0.5,1.5]);
   plot(lfpt,gamma{8}); xlim([-0.5,1]); ylim([-0.5,1.5]);legend({'Light ON injured','Light OFF injured','Light ON Uninjured','Light OFF Unjured'});
    gammma_summary=cell2mat(cellfun(@(x) squeeze(nanmean(nanmean(x(lfpt>0.2&lfpt<0.5,35:100,:),2),1)),lfp_correct,'UniformOutput',0));
    
%%
function response=negativeresponse(untaggedbinspikes,t)
for i=1:size(untaggedbinspikes,2)
x=mean(untaggedbinspikes(t>=0&t<=0.1,i),1);
if x<mean(untaggedbinspikes(t>-0.5&t<0,i),1)-0.5*std(untaggedbinspikes(t>-0.5&t<0,i))
    response(i)=1;
else
    response(i)=0;
end
end
end
function datamat=cal_binspikes(datamat,Fs,fieldname,eventtype) 
tmpmat=eval(['datamat.',fieldname,';']);
try
spiketime=arrayfun(@(x) x{:}+tmpmat.EVTinfo.timerange(1),tmpmat.SPKdata,'UniformOutput',0);
catch
    spiketime=arrayfun(@(x) x{:}-tmpmat.EVTinfo.timestart(1),tmpmat.SPKdata,'UniformOutput',0);
end
for k=1:size(spiketime,1)
    for i=1:size(spiketime,2)
    st.spiketime=spiketime{k,i};
    try
    binspike{k}(:,i)=binspikes(st.spiketime,Fs,tmpmat.EVTinfo.timerange);
    catch
         binspike{k}(:,i)=binspikes(st.spiketime,Fs,[0,round(tmpmat.EVTinfo.timestop-tmpmat.EVTinfo.timestart)]);
    end
    end
end
% for i=1:length(eventtype)
%     index=ismember(tmpmat.EVTinfo.eventdescription,eventtype{i});
%     for k=1:size(spiketime)
%         tmpspiketime=spiketime(k,index);
%         for j=1:length(tmpspiketime)
%             st(j).spiketime=tmpspiketime{j};
%         end
%         [tmppsth,tmpmat.pstht]=psth(st,0.05,'n',[-2,4]);
%         try
%         eval(['tmpmat.',eventtype{i},'psth(:,k)=tmppsth;']);
%         catch ME
%             tmppsth=zeros(1,length(tmpmat.pstht));
%             eval(['tmpmat.',eventtype{i},'psth(:,k)=tmppsth;']);
%         end
%     end
% end
    tmpmat.binspike=binspike;
    eval(['datamat.',fieldname,'=tmpmat;']);
end% using the psth from chronux toolbox ().
function datamat=cal_optotag(datamat)
    tmpmat=datamat.Combined;
    index=ismember(tmpmat.EVTinfo.eventdescription,'OptoTag');
    tmpspiketime=arrayfun(@(x) x{:}+tmpmat.EVTinfo.timerange(1),tmpmat.SPKdata(:,index),'UniformOutput',0);
    for k=1:size(tmpspiketime)
        for i=1:size(tmpspiketime,2)
            st.spiketime=tmpspiketime{k,i};
            tmpbinspike{k}(:,i)=binspikes(st.spiketime,1000,tmpmat.EVTinfo.timerange);
        end
    end
    for i=1:length(tmpbinspike)
    t=linspace(-2,4,size(tmpbinspike{i},1));    
    p(i)=salt(tmpbinspike{i}(t<0&t>-0.5,:)',tmpbinspike{i}(t>0,:)',0.001,0.005);
    if p(i)<0.05
        tmpmat.Optoresponse(i)=1;
    else
        firingbefore=mean(tmpbinspike{i}(t<0&t>-0.5,:),1);
        firingafter=mean(tmpbinspike{i}(t>0&t<0.5,:),1);
        [~,p_tmp]=ttest(firingbefore,firingafter);
        if p_tmp<0.05 && mean(firingbefore)>mean(firingafter)
            tmpmat.Optoresponse(i)=-1;
        else
            tmpmat.Optoresponse(i)=0;
        end
    end
    end
    datamat.Combined=tmpmat;
end
function datamat=cal_response(datamat,eventtype)
timerange=[-1,1];
tmpmat=datamat.Combined;
for i=1:length(eventtype)
index=ismember(tmpmat.EVTinfo.eventdescription,eventtype{i});
spiketime=arrayfun(@(x) x{:}+tmpmat.EVTinfo.timerange(1),tmpmat.SPKdata(:,index),'UniformOutput',0);
spikebefore=cellfun(@(x) length(find(x<0&x>timerange(1)))/abs(timerange(1)),spiketime,'UniformOutput',1);
spikeafter=cellfun(@(x) length(find(x>0&x<timerange(2)))/abs(timerange(2)),spiketime,'UniformOutput',1); 

    [~,p]=ttest(spikebefore',spikeafter');
    for k=1:length(p)
    if mean(spikebefore(k,:))>mean(spikeafter(k,:)) && p(k)<0.05
        response(k)=-1;
    elseif mean(spikebefore(k,:))<mean(spikeafter(k,:)) && p(k)<0.05
        response(k)=1;
    else
        response(k)=0;
    end
    end
    eval(['tmpmat.',eventtype{i},'response=response;']);
end
    datamat.Combined=tmpmat;
end
function datamat=cal_spectrogram(datamat)
tmpmat=datamat.Combined;
for i=1:length(tmpmat.LFPdata)
    if size(tmpmat.LFPdata{i},1)>7501
        tmpmat.LFPdata{i}(end,:)=[];
    end
    %tmpmat.Spectrom(:,:,i)=abs(awt_freqlist(mean(tmpmat.LFPdata{i},2),tmpmat.LFPinfo.Fs,1:100));
    [~,tmp]=sub_stft(mean(tmpmat.LFPdata{i},2),linspace(-2,4,7501),linspace(-2,4,7501),1:100,tmpmat.LFPinfo.Fs,0.2);
    tmpmat.Spectrom_STFT(:,:,i)=tmp';
end
datamat.Combined=tmpmat;
end
    