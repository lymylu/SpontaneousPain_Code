%% cal the spike-field coherence among the different brain regions in both spontaneous pain epochs and no spotaneous pain epochs
% final version 2023/7/2
% In fact the file contains all type of events (include 6 type stimuli, and
% 3 type duration) is only 7 file
% the file contains 6 type stimuli and spontaneous pain is for
% multimodalties analysis
BaseDir='ExtractData'; 
filelist=dir(BaseDir);
for i=3:length(filelist)
    try
        tmpmat=matfile(fullfile(BaseDir,filelist(i).name),'Writable',true);

    tmpmat=cal_durationresponse(tmpmat,'Spontaneouspain','Spontaneousnopain');
    end
    try
    tmpmat=cal_durationresponse(tmpmat,'Move','NoMove');
    end
    try
    tmpmat=cal_durationresponse(tmpmat,'Groom','NoGroom');
    end
    try
    tmpmat=cal_stimulusresponse(tmpmat,{'AirPuff','Brush','Laser','VonfreyHigh','VonfreyLow','Pin','Sound'});
%       tmpmat=cal_PSD(tmpmat,{'Spontaneouspain','Spontaneousnopain'});
%        tmpmat=cal_SFC(tmpmat,{'Spontaneousnopain','Spontaneouspain'});
    catch ME
    disp(ME);
    end
end
%% for review, summarized the neurons responses to spontaneous pain in famale rats
for i=3:length(filelist)
    tmpmat=matfile(fullfile(BaseDir,filelist(i).name));
    Spontaneous(i-2)=NeuroResult(tmpmat.Spontaneouspain);
    NoSpontaneous(i-2)=NeuroResult(tmpmat.Spontaneousnopain);
end
Spontaneouspain_all=Spontaneous.CollectSpikeVariables({'SPKinfo.spikename','SPKinfo.SPKchanneldescription','SPKinfo.putativeCellType','SPKinfo.firingRate','SPKinfo.troughToPeak','durationfrequency','permutationfreq','response'},[2,2,2,2,2,2,2,2]);
Spontaneousnopain_all=NoSpontaneous.CollectSpikeVariables({'SPKinfo.spikename','SPKinfo.SPKchanneldescription','SPKinfo.putativeCellType','SPKinfo.firingRate','SPKinfo.troughToPeak','durationfrequency','timeall'},[2,2,2,2,2,2,2]);


%% summary Spontaneous responses  used after run Spikeresponse_cal.m
% Group 1 is sham group (not used in this manuscript.
% Only data with time from 1d to 14d were used in this manuscript.
Group={'(\<RAT6)|(\<RAT11)|(\<RAT14)|(\<RAT16)','(\<RAT4)|(\<RAT7)|(\<RAT8)|(\<RAT15)|(\<RAT17)|(\<RAT19)|(\<RAT21)|(\<RAT22)|(\<RAT24)'};
Time={'(_baseline)|(_Baseline)|(_BASELINE)','(_1d_)|(_3d_)|(_7d_)|(_9d_)|(_14d_)'}; % baseline ,early stage or late stage
namelist=struct2table(filelist);
namelist=namelist.name;
index=cellfun(@(x) ~isempty(regexpi(x,Time{2},'match')),namelist,'UniformOutput',1)&cellfun(@(x) ~isempty(regexpi(x,Group{2},'match')),namelist,'UniformOutput',1);
filelist=filelist(index);
%%
n=1;m=1;k=1;n1=1;n2=1;n3=1;n4=1;n5=1;n6=1;n7=1;
% Spontaneous=[];NoSpontaneous=[];Groom=[];NoGroom=[];Move=[];NoMove=[];AirPuff=[];Brush=[];Pin=[];VonfreyHigh=[];VonfreyLow=[];Sound=[];
modalities={'Pin','VonfreyHigh','VonfreyLow','Brush','Laser','AirPuff','Sound','Spontaneouspain','Groom','Move','Spontaneousnopain','NoGroom', 'NoMove'};
modindex={'n1','n2','n3','n4','n5','n6','n7'};
m1=1;m2=1;m3=1;
for i=1:length(filelist)
        tmpmat=matfile(fullfile(BaseDir,filelist(i).name),'Writable',true); 
        invalid=contains(modalities,fieldnames(tmpmat));
         if sum(invalid(1:8))==8 %%  contains all stimulus type
            for j=1:length(modalities)-6
                eval([modalities{j},'(',modindex{j},')=NeuroResult(tmpmat.',modalities{j},');']);
                eval([modindex{j},'=',modindex{j},'+1;']);
            end
         end
        if sum(invalid(9))==1 % contains Groom
             Groom(m1)=NeuroResult(tmpmat.Groom);
             NoGroom(m1)=NeuroResult(tmpmat.NoGroom);
             m1=m1+1;
        end
        if sum(invalid(10)==1)% contains Move
             Move(m2)=NeuroResult(tmpmat.Move);
             NoMove(m2)=NeuroResult(tmpmat.NoMove);
             m2=m2+1;
        end
        if sum(invalid(8)==1) % contains Spontaneous pain
            Spontaneous(m3)=NeuroResult(tmpmat.Spontaneouspain);
            NoSpontaneous(m3)=NeuroResult(tmpmat.Spontaneousnopain);  
             m3=m3+1;
        end
end
%% cat all Spike data from the NeuroResult objects
Spontaneouspain_all=Spontaneous.CollectSpikeVariables({'SPKinfo.spikename','SPKinfo.SPKchanneldescription','SPKinfo.putativeCellType','SPKinfo.firingRate','SPKinfo.troughToPeak','durationfrequency','permutationfreq','response'},[2,2,2,2,2,2,2,2]);
Spontaneousnopain_all=NoSpontaneous.CollectSpikeVariables({'SPKinfo.spikename','SPKinfo.SPKchanneldescription','SPKinfo.putativeCellType','SPKinfo.firingRate','SPKinfo.troughToPeak','durationfrequency','timeall'},[2,2,2,2,2,2,2]);
Groom_all=Groom.CollectSpikeVariables({'SPKinfo.spikename','SPKinfo.SPKchanneldescription','SPKinfo.putativeCellType','SPKinfo.firingRate','SPKinfo.troughToPeak','durationfrequency','permutationfreq','response'},[2,2,2,2,2,2,2,2]);
NoGroom_all=NoGroom.CollectSpikeVariables({'SPKinfo.spikename','SPKinfo.SPKchanneldescription','SPKinfo.putativeCellType','SPKinfo.firingRate','SPKinfo.troughToPeak','durationfrequency','timeall'},[2,2,2,2,2,2,2]);
Move_all=Move.CollectSpikeVariables({'SPKinfo.spikename','SPKinfo.SPKchanneldescription','SPKinfo.putativeCellType','SPKinfo.firingRate','SPKinfo.troughToPeak','durationfrequency','permutationfreq','response'},[2,2,2,2,2,2,2,2]);
NoMove_all=NoMove.CollectSpikeVariables({'SPKinfo.spikename','SPKinfo.SPKchanneldescription','SPKinfo.putativeCellType','SPKinfo.firingRate','SPKinfo.troughToPeak','durationfrequency','timeall'},[2,2,2,2,2,2,2]);
for i=1:length(modalities)-6
    eval([modalities{i},'_all=',modalities{i},'.CollectSpikeVariables({''response'',''SPKinfo.SPKchanneldescription'',''meanbinspike'',},[2,2,1]);']);
    tmp=eval([modalities{i}]);
    for j=1:length(tmp)
        trialnum(i,j)=size(tmp(j).binspike,3);
    end
end
% generate the spike matrix for naive bayes classifier. 
% min trial number of all modalities is 15.
category=[];spikematrix=[];
for i=1:length(modalities)-6
    tmp=eval([modalities{i}]);
    modalityspike=[];
    for j=1:length(tmp)
        modalityspike=cat(1,modalityspike,tmp(j).binspike(:,:,randperm(size(tmp(j).binspike,3),15)));
    end
    spikematrix=cat(3,spikematrix,modalityspike);
    category=cat(1,category,repmat(modalities(i),[15,1]));
end

%% hypergeometric analysis for the multimodalities overlap.
index=ismember(Spontaneouspain_all.SPKinfo.SPKchanneldescription,'PL')&Spontaneouspain_all.SPKinfo.troughToPeak>0.5;
valid=true(length(index),1);
binspiket=linspace(-2,2,41);
for i=1:length(modalities)-6
    tmp1=eval([modalities{i},'_all.meanbinspike;']);
    tmp1=basecorrect(tmp1',binspiket,-2,0,'zscore')';
    tmp1(isinf(tmp1))=nan;
    valid=valid&~isnan(tmp1(:,1));
end
index=index&valid';
spikeindex=index;
% plot throughToPeak of valid neurons.
troughToPeak=Spontaneouspain_all.SPKinfo.troughToPeak(ismember(Spontaneouspain_all.SPKinfo.SPKchanneldescription,'PL')&valid');
figure;histogram(troughToPeak,linspace(0.1,1.8,18));xlim([0.1,1.5]);
for i=1:length(modalities)-5
    for j=1:length(modalities)-5
        tmp1=eval([modalities{i},'_all.response(index);']);
        try tmp1=cell2mat(tmp1); end
        tmp1(isnan(tmp1))=0;
        tmp2=eval([modalities{j},'_all.response(index);']);
        try tmp2=cell2mat(tmp2);end
        tmp2(isnan(tmp2))=0;
        overlap(i,j)=sum((tmp1==1&tmp2==1)|(tmp1==-1&tmp2==-1));
        %overlap(i,j)=sum((tmp1==1&tmp2==1));
        p(i,j)=hygepdf(overlap(i,j),sum(index),sum(tmp1==1|tmp1==-1),sum(tmp2==1|tmp2==-1));
        if overlap(i,j)<=hygestat(sum(index),sum(tmp1==1|tmp1==-1),sum(tmp2==1|tmp2==-1))&& p(i,j)<0.05
          maxtrix(i,j)=-1;
        elseif overlap(i,j)>hygestat(sum(index),sum(tmp1==1|tmp1==-1),sum(tmp2==1|tmp2==-1))&& p(i,j)<0.05
            matrix(i,j)=1;
        else
            matrix(i,j)=nan;
        end
    end
end
% Figure 4D
figure;subplot(1,3,1);heatmap(modalities(1:8),modalities(1:8),overlap);
subplot(1,3,2);imagesc(p);caxis([0,0.05]);
subplot(1,3,3);heatmap(modalities(1:8),modalities(1:8),matrix);

%% Naive Bayes classifier
index=ismember(Spontaneouspain_all.SPKinfo.SPKchanneldescription,'PL')&Spontaneouspain_all.SPKinfo.troughToPeak>0.5;
valid=true(length(index),1);
binspiket=linspace(-2,2,41);
for i=1:length(modalities)-6
    tmp1=eval([modalities{i},'_all.meanbinspike;']);
    tmp1=basecorrect(tmp1',binspiket,-2,0,'zscore')';
    tmp1(isinf(tmp1))=nan;
    valid=valid&~isnan(tmp1(:,1));
end
index=index&valid';
spikeindex=index;
%%
rng("default");
% for replicable
positiveresponse=cellfun(@(x) x==1,Spontaneouspain_all.response,'UniformOutput',1);
negativeresponse=cellfun(@(x) x==-1,Spontaneouspain_all.response,'UniformOutput',1);
spikematrix1=spikematrix(spikeindex,:,:);
dataset=[];
for i=1:size(spikematrix1,3)
    for j=1:size(spikematrix1,1)
        dataset(j,i)=svd(spikematrix1(j,21:30,i));
    end
end
[cmatrix,pmatrix,cmatrix_shuffle]=naiveBayers(dataset,category,modalities(1:7));
%Figure 4E
figure; subplot(1,2,1); imagesc(1:7,1:7,cmatrix'./15);
xticklabels(modalities(1:7)); yticklabels(modalities(1:7));
subplot(1,2,2);imagesc(1:7,1:7,pmatrix');clim([0,0.05]);
xticklabels(modalities(1:7)); yticklabels(modalities(1:7));
%%
rng(1);
spikematrix2=spikematrix(spikeindex&~positiveresponse&~negativeresponse,:,:);
dataset=[];
for i=1:size(spikematrix2,3)
    for j=1:size(spikematrix2,1)
        dataset(j,i)=svd(spikematrix2(j,21:30,i));
    end
end
[cmatrix_exclude,pmatrix_exclude,cmatrix_shuffle_exclude]=naiveBayers(dataset,category,modalities(1:7));
% Figure 4F
figure; subplot(1,2,1); imagesc(1:7,1:7,cmatrix_exclude'./15);
xticklabels(modalities(1:7)); yticklabels(modalities(1:7));
subplot(1,2,2);imagesc(1:7,1:7,pmatrix_exclude');clim([0,0.05]);
xticklabels(modalities(1:7)); yticklabels(modalities(1:7));
%% 
rng(1);
spikematrix2=spikematrix(spikeindex&~positiveresponse&~negativeresponse,:,:); % excluded spon responses amount
dataset=[];
for i=1:size(spikematrix2,3)
    for j=1:size(spikematrix2,1)
        dataset(j,i)=svd(spikematrix2(j,21:30,i));
    end
end
accurarcy=naiveBayers2(dataset,category,modalities(1:7));
%
for i=1:50
    spikematrix3=spikematrix(spikeindex,:,:); % random sampled amount
    randindex=randperm(349,230); 
    spikematrix3=spikematrix3(randindex,:,:);
    dataset=[];
    for a=1:size(spikematrix3,3)
    for b=1:size(spikematrix3,1)
        dataset(b,a)=svd(spikematrix3(b,21:30,a));
    end
    end
    accurarcy_permute(:,i)=naiveBayers2(dataset,category,modalities(1:7));
end
%
rng("default");
spikematrix3=spikematrix(spikeindex,:,:); %origin amount
dataset=[];
for i=1:size(spikematrix3,3)
    for j=1:size(spikematrix3,1)
        dataset(j,i)=svd(spikematrix3(j,21:30,i));
    end
end
accurarcy_origin=naiveBayers2(dataset,category,modalities(1:7));
figure; % kernel fit to the distribution of NB accuracy
% Figure S2H
for i=1:7
    subplot(1,7,i)
    histfit(accurarcy_permute(i,:),10,'kernel');
    line([accurarcy(i),accurarcy(i)],[0,20],'color','r');
    line([accurarcy_origin(i),accurarcy_origin(i)],[0,20],'color','b');title(modalities(i));
    xlim([0,1]);
end
%% plot spike show for each stimulus modalities 
%index=ismember(Spontaneouspain_all.SPKinfo.SPKchanneldescription,'PL')&Spontaneouspain_all.SPKinfo.troughToPeak>0.5;
painresponse=Spontaneouspain_all.response(index);
[painresponse,index2]=sort(cell2mat(painresponse),'descend');
binspiket=linspace(-2,2,41);
figure;% Figure 3B
for i=1:length(modalities)-6
    tmp1=eval([modalities{i},'_all.meanbinspike(index,:);']);
    tmp1=basecorrect(tmp1',binspiket,-2,0,'zscore');
    tmp1(isinf(tmp1))=nan;
    tmp1=tmp1(:,index2);
    if i==1
    [~,indexa]=sort(mean(tmp1(21:end,painresponse==1),1));
    [~,indexb]=sort(mean(tmp1(21:end,painresponse==0),1));
    [~,indexc]=sort(mean(tmp1(21:end,painresponse==-1),1));
    end
    tmpa=tmp1(:,painresponse==1);
    tmpb=tmp1(:,painresponse==0);
    tmpc=tmp1(:,painresponse==-1);
    tmp1=cat(2,tmpa(:,indexa),tmpb(:,indexb),tmpc(:,indexc));
    tmp1=smoothdata(tmp1,'movmean',5);
    subplot(1,7,i)
    imagesc(binspiket,1:size(tmp1,2),tmp1'); caxis([-2,2]); title(modalities{i});
    colormap jet;
    if i==7
        hold on;
        fill([1,1,2,2],[0,sum(painresponse==1),sum(painresponse==1),0],'r');
        fill([1,1,2,2],[sum(painresponse==1),sum(painresponse==1)+sum(painresponse==0),sum(painresponse==1)+sum(painresponse==0),sum(painresponse==1)],'black');
        fill([1,1,2,2],[sum(painresponse==1)+sum(painresponse==0),length(painresponse),length(painresponse),sum(painresponse==0)+sum(painresponse==1)],'blue');
    end
end
%% plot the averaged spike firing of each reflexive modality
% not showed in the manuscript.
figure;
for i=1:length(modalities)-6
    tmp1=eval([modalities{i},'_all.meanbinspike(index,:);']);
    tmp1=basecorrect(tmp1',binspiket,-2,0,'zscore');
    tmp1(isinf(tmp1))=nan;
    tmp1=tmp1(:,index2);
    if i==1
    [~,indexa]=sort(mean(tmp1(21:end,painresponse==1),1));
    [~,indexb]=sort(mean(tmp1(21:end,painresponse==0),1));
    [~,indexc]=sort(mean(tmp1(21:end,painresponse==-1),1));
    end
    tmpa=tmp1(:,painresponse==1);
    tmpb=tmp1(:,painresponse==0);
    tmpc=tmp1(:,painresponse==-1);
    firingchange{1}(:,i)=mean(tmpa(21:30,:),1)';
    firingchange{2}(:,i)=mean(tmpb(21:30,:),1)';
    firingchange{3}(:,i)=mean(tmpc(21:30,:),1)';
%     tmpa=smoothdata(tmpa,'movmean',5);
%     tmpb=smoothdata(tmpb,'movmean',5);
%     tmpc=smoothdata(tmpc,'movmean',5);
    subplot(1,7,i) 
    hold on;
    shadebar(binspiket(1:40),tmpa(1:40,:),'r'); ylim([-0.5,1.5]);
    shadebar(binspiket(1:40),tmpb(1:40,:),'black'); ylim([-0.5,1.5]);
    shadebar(binspiket(1:40),tmpc(1:40,:),'blue'); ylim([-0.5,1.5]);title(modalities{i});
end
%%
tmpvariable=Spontaneouspain_all;
notmpvariable=Spontaneousnopain_all;
earlyindex=cellfun(@(x) ~isempty(regexpi(x,'(_1d_)','match')),tmpvariable.Subjectname,'UniformOutput',1);
lateindex=cellfun(@(x) isempty(regexpi(x,'(_1d_)','match')),tmpvariable.Subjectname,'UniformOutput',1);
PL_pyrindex=ismember(tmpvariable.SPKinfo.SPKchanneldescription,'PL')&tmpvariable.SPKinfo.troughToPeak>0.5&~isnan(cell2mat(tmpvariable.response));
IL_pyrindex=ismember(tmpvariable.SPKinfo.SPKchanneldescription,'IL')&tmpvariable.SPKinfo.troughToPeak>0.5&~isnan(cell2mat(tmpvariable.response));
PL_intindex=ismember(tmpvariable.SPKinfo.SPKchanneldescription,'PL')&tmpvariable.SPKinfo.troughToPeak<0.5&~isnan(cell2mat(tmpvariable.response));
IL_intindex=ismember(tmpvariable.SPKinfo.SPKchanneldescription,'IL')&tmpvariable.SPKinfo.troughToPeak<0.5&~isnan(cell2mat(tmpvariable.response));
PL_pyr=tabulate(cell2mat(tmpvariable.response(PL_pyrindex)));
IL_pyr=tabulate(cell2mat(tmpvariable.response(IL_pyrindex)));
PL_int=tabulate(cell2mat(tmpvariable.response(PL_intindex)));
IL_int=tabulate(cell2mat(tmpvariable.response(IL_intindex)));
a={cell2mat(tmpvariable.durationfrequency(PL_pyrindex))',cell2mat(notmpvariable.durationfrequency(PL_pyrindex))',cell2mat(tmpvariable.durationfrequency(IL_pyrindex))',cell2mat(notmpvariable.durationfrequency(IL_pyrindex))'};
b={cell2mat(tmpvariable.durationfrequency(PL_intindex))',cell2mat(notmpvariable.durationfrequency(PL_intindex))',cell2mat(tmpvariable.durationfrequency(IL_intindex))',cell2mat(notmpvariable.durationfrequency(IL_intindex))'};
figure; % Figure 2E and Figure S2C :total amount of PL and IL
subplot(2,2,1)
pie(PL_pyr(:,2),num2str(PL_pyr(:,2))); title('PL PYR');
subplot(2,2,2);
pie(IL_pyr(:,2),num2str(IL_pyr(:,2))); title('IL PYR');
subplot(2,2,3)
pie(PL_int(:,2),num2str(PL_int(:,2))); title('PL Int');
subplot(2,2,4);
pie(IL_int(:,2),num2str(IL_int(:,2))); title('IL Int');
%% show the overlap between groom, spontaneous and move
MoveSubject=cellfun(@(x,y) [x,y],Move_all.Subjectname,Move_all.SPKinfo.spikename,'UniformOutput',0);
PainSubject=cellfun(@(x,y) [x,y],Spontaneouspain_all.Subjectname,Spontaneouspain_all.SPKinfo.spikename,'UniformOutput',0);
GroomSubject=cellfun(@(x,y) [x,y],Groom_all.Subjectname,Groom_all.SPKinfo.spikename,'UniformOutput',0);
Subjectname=intersect(intersect(MoveSubject,GroomSubject),PainSubject);
Movevalid=ismember(MoveSubject,Subjectname);
Groomvalid=ismember(GroomSubject,Subjectname);
Painvalid=ismember(PainSubject,Subjectname);
List=[];n=1;
Response=[1,0,-1];
for i=1:length(Response)
    index1=ismember(Move_all.SPKinfo.SPKchanneldescription,'PL')&Move_all.SPKinfo.troughToPeak>0.5&(cell2mat(Move_all.response)==Response(i));
    index1=index1(Movevalid);
    for j=1:length(Response)
        index2=ismember(Spontaneouspain_all.SPKinfo.SPKchanneldescription,'PL')&Spontaneouspain_all.SPKinfo.troughToPeak>0.5&(cell2mat(Spontaneouspain_all.response)==Response(j));
        index2=index2(Painvalid);
        List=cat(1,List,[{['Move ',num2str(Response(i))]},{sum(index1&index2)},{['Pain ',num2str(Response(j))]}]);
        n=n+1;
    end
end
% Venn plot from Move and Pain
Subjectname=intersect(MoveSubject,PainSubject);
Movevalid=ismember(MoveSubject,Subjectname);
Painvalid=ismember(PainSubject,Subjectname);
indexpain=ismember(Spontaneouspain_all.SPKinfo.SPKchanneldescription,'PL')&Spontaneouspain_all.SPKinfo.troughToPeak>0.5&Painvalid;
Responsepain=cell2mat(Spontaneouspain_all.response(indexpain));
responsepain=sum(Responsepain==1|Responsepain==-1);
indexmove=ismember(Move_all.SPKinfo.SPKchanneldescription,'PL')&Move_all.SPKinfo.troughToPeak>0.5&Movevalid;
Responsemove=cell2mat(Move_all.response(indexmove));
responsemove=sum(Responsemove==1|Responsemove==-1);
overlay=sum(Responsemove==1&Responsepain==1);
% Figure 3C
figure;subplot(1,2,1);venn([responsepain,responsemove],overlay);
p=hygepdf(1:min([responsemove,responsepain]),length(Responsepain),responsemove,responsepain);subplot(1,2,2);plot(1:1:min([responsemove,responsepain]),p); hold on; line([overlay,overlay],[overlay,0]); ylim([0,0.15]);

% Venn plot from Groom and Pain
Subjectname=intersect(GroomSubject,PainSubject);
Groomvalid=ismember(GroomSubject,Subjectname);
Painvalid=ismember(PainSubject,Subjectname);
indexpain=ismember(Spontaneouspain_all.SPKinfo.SPKchanneldescription,'PL')&Spontaneouspain_all.SPKinfo.troughToPeak>0.5&Painvalid;
Responsepain=cell2mat(Spontaneouspain_all.response(indexpain));
responsepain=sum(Responsepain==1|Responsepain==-1);
indexgroom=ismember(Groom_all.SPKinfo.SPKchanneldescription,'PL')&Groom_all.SPKinfo.troughToPeak>0.5&Groomvalid;
Responsegroom=cell2mat(Groom_all.response(indexgroom));
responsegroom=sum(Responsegroom==1|Responsepain==-1);
overlay=sum(Responsegroom==1&Responsepain==1);
% Figure 3C
figure;subplot(1,2,1);venn([responsepain,responsegroom],overlay);
p=hygepdf(1:min([responsegroom,responsepain]),length(Responsepain),responsegroom,responsepain);subplot(1,2,2);plot(1:1:min([responsegroom,responsepain]),p); hold on; line([overlay,overlay],[overlay,0]); ylim([0,0.15]);
for i=1:length(Response)
    index1=ismember(Spontaneouspain_all.SPKinfo.SPKchanneldescription,'PL')&Spontaneouspain_all.SPKinfo.troughToPeak>0.5&(cell2mat(Spontaneouspain_all.response)==Response(i));
    index1=index1(Painvalid);
    for j=1:length(Response)
        index2=ismember(Groom_all.SPKinfo.SPKchanneldescription,'PL')&Groom_all.SPKinfo.troughToPeak>0.5&(cell2mat(Groom_all.response)==Response(j));
        index2=index2(Groomvalid);
        List=cat(1,List,[{['Pain ',num2str(Response(i))]},{sum(index1&index2)},{['Groom ',num2str(Response(j))]}]);
        n=n+1;
    end
end  
%% response neuron show (positive=1136 negative=1131 no response=1149 for spontaneous pain,
% Figure 2D
tmpvariable=Spontaneouspain_all;
notmpvariable=Spontaneousnopain_all;
tmpvariablenew=Spontaneous;
for c=1136 % or 1131 or 1149
index=c;
if tmpvariable.response{index}==-1
spiketime=Spontaneousnopain_all.timeall{index};
spikefrequency=tmpvariable.durationfrequency{index};
permutationfreq=tmpvariable.permutationfreq{index};
Subjectname=tmpvariable.Subjectname{index};
for i=1:length(tmpvariablenew)
index_sub(i)=strcmp(tmpvariablenew(i).Subjectname,Subjectname);
end
yestimerange=tmpvariablenew(index_sub).SPKinfo.timerange;
Subjectname=[tmpvariable.Subjectname{index},tmpvariable.SPKinfo.spikename{index}];
Spontaneousshow(spiketime,spikefrequency,Subjectname,permutationfreq,yestimerange);
end
end
function tmpdata=cal_durationresponse(tmpdata,yesdataname,nodataname)
    % cal the spike responses to the Spontaneous behavior
    % the classification of responses using the comparision of the true
    % firingrate of goal segements among that of the surrogate segments.
    nodata=eval(['tmpdata.',nodataname,';']);
    yesdata=eval(['tmpdata.',yesdataname,';']);
%     frequency=sum(cellfun(@(x) length(x),yesdata.SPKdata,'UniformOutput',1))/sum(yesdata.timerange(:,2)-yesdata.timerange(:,1));
%     nofrequency=sum(cellfun(@(x) length(x),nodata.spiketime,'UniformOutput',1))/sum(nodata.werange(:,2)-nodata.timerange(:,1));
    % permutation test?
    % nopain data and pain data
    nodata=NeuroResult(nodata);
    yesdata=NeuroResult(yesdata);
    nodata=nodata.Splice2Split;
    yesdata=yesdata.Splice2Split;
    nodata=nodata.NeuroResult2Struct;
    yesdata=yesdata.NeuroResult2Struct;
    notimerange=[min(reshape(nodata.SPKinfo.timerange,[],1)),max(reshape(nodata.SPKinfo.timerange,[],1))];
    yestimerange=[min(reshape(yesdata.SPKinfo.timerange,[],1)),max(reshape(yesdata.SPKinfo.timerange,[],1))];
    timerangeall=[min(min(notimerange),min(yestimerange)),max(max(notimerange),max(yestimerange))];
    yestimeall=sum(yesdata.SPKinfo.timerange(:,2)-yesdata.SPKinfo.timerange(:,1)); % the total time of spontaneouspain 
    notimeall=timerangeall(2)-timerangeall(1)-yestimeall; % the total time of spontaneousnopain
    notime=cell(1,size(nodata.SPKdata,1));
    yestime=cell(1,size(yesdata.SPKdata,1));
    for i=1:size(nodata.SPKdata,1)
        for j=1:size(nodata.SPKdata,2)
            notime{i}=cat(1,notime{i},nodata.SPKdata{i,j});
        end
        nofrequency{i}=length(notime{i})/notimeall;
    end
    for i=1:size(yesdata.SPKdata,1)
        for j=1:size(yesdata.SPKdata,2)
            yestime{i}=cat(1,yestime{i},yesdata.SPKdata{i,j});
        end
        yesfrequency{i}=length(yestime{i})/yestimeall;
    end
    spiketime=cellfun(@(x,y) sort(cat(1,x,y)),yestime,notime,'UniformOutput',0); % the spike time list from the whole time duration
    permutationpaintimerange=getPermutationtimerange(yesdata.SPKinfo.timerange,nodata.SPKinfo.timerange,notimeall,500);
    [response,permutationfreq]=cellfun(@(x,y) cal_permutationresponse(x,permutationpaintimerange,yestimeall,y),spiketime,yesfrequency,'UniformOutput',0);
    nodata.durationfrequency=nofrequency;
    yesdata.durationfrequency=yesfrequency;
    yesdata.response=response;
    yesdata.permutationfreq=permutationfreq;
    nodata.timeall=spiketime;
    eval(['tmpdata.',nodataname,'=nodata;']);
    eval(['tmpdata.',yesdataname,'=yesdata;']);
end
function permutationpaintimerange=getPermutationtimerange(yestimerange,notimerange,notimeall,nperm)
        rng("default");
        for i=1:nperm
        seqofyes=randperm(size(yestimerange,1)); % randperm the sequence of pain time segment
%         notimeseq=diff([sort(randsample(1:notimeall,length(yestimerange)+1)),notimeall(end)]); % randperm the time sequences between permutation spontanouespain.
        notimeseq=diff(sort([notimeall*rand(size(yestimerange,1)+2,1)]));
        yesduration=yestimerange(seqofyes,2)-yestimerange(seqofyes,1); % each time duration of the permutation spontaneouspain.
        timerangeall=[min(min(notimerange),min(yestimerange)),max(max(notimerange),max(yestimerange))];
        for j=1:length(notimeseq)
            if j==2
            permutationpaintimerange(j-1,1,i)=sum(notimeseq(1:j-1));
            permutationpaintimerange(j-1,2,i)=permutationpaintimerange(j-1,1,i)+yesduration(j-1);
            elseif j>2
                permutationpaintimerange(j-1,1,i)=sum(notimeseq(1:j-1))+sum(yesduration(1:j-2));
                permutationpaintimerange(j-1,2,i)=permutationpaintimerange(j-1,1,i)+yesduration(j-1);
            end
        end
        permutationpaintimerange(:,:,i)=permutationpaintimerange(:,:,i)+timerangeall(1);
        if permutationpaintimerange(j-1,2,i)>timerangeall(end)
            a=1;
        end
    end
end
function [response,permutationfreq]=cal_permutationresponse(spiketime,permutationpaintimerange,yesduration,yesfrequency)
        permutationfreq=[];
        for i=1:size(permutationpaintimerange,3) 
            permutationspike=[];
        for j=1:size(permutationpaintimerange,1)
            permutationspike=cat(1,permutationspike,length(find(spiketime>permutationpaintimerange(j,1,i)&spiketime<permutationpaintimerange(j,2,i))));
        end
            permutationfreq(i)=sum(permutationspike)/sum(yesduration);
        end
    p=length(find(permutationfreq>yesfrequency))/size(permutationpaintimerange,3);
    if ~isempty(spiketime)
        if p>=0.95
        response=-1;
        elseif p<=0.05
        response=1;
        else
        response=0;
        end
    else
        response=0;
    end
end
function Spontaneousshow(spiketime,spikefrequency,spikename,permutationfreq,paintimerange)
figure;
subplot(1,3,[2,3]);
binsize=10;
binrange=linspace(min(spiketime),max(spiketime),(max(spiketime)-min(spiketime))/binsize+1);
a=hist(spiketime,binrange);
hist(spiketime,binrange);
tmptimerange=paintimerange;
xlim([1000,1700]);
hold on
for i=1:size(tmptimerange,1)
    fill([tmptimerange(i,1),tmptimerange(i,1),tmptimerange(i,2),tmptimerange(i,2)],[0,max(a),max(a),0],'r','FaceAlpha',.3,'EdgeAlpha',0);
end
title(spikename);
line([min(spiketime),max(spiketime)],[mean(permutationfreq)*binsize,mean(permutationfreq)*binsize]);
ylabel('Spike counts per bin');
xlabel('Time (s)'); axis tight;
subplot(1,3,1);
a=hist(permutationfreq,20);
hist(permutationfreq,20);
hold on;
line([spikefrequency,spikefrequency],[0,max(a)]);
axis tight;
ylabel('Permutation counts');
xlabel('Spike frequency');
end
function structdata=cal_binspikefrom_LFPspike(structdata)
spiketime=cellfun(@(x) x+structdata.EVTinfo.timerange(1),structdata.SPKdata,'UniformOutput',0);
for i=1:size(spiketime,1)
    st=[];
    try
    for j=1:size(spiketime,2)
        st(j).spiketime=spiketime{i,j}; 
        if ~isempty(st(j).spiketime)
            [gaussspike{i}(:,j),t]=psth(st,0.05,'n',structdata.EVTinfo.timerange);
        end
    end
    catch ME
        disp(ME);
    end
      % gaussian smooth the stimulus data;
end
try
  structdata.Result.gaussspike=gaussspike;
  structdata.Result.gausst=t;
catch
    structdata.Result.gaussspike=[];
    structdata.Result.gausst=t;
end
end % using the psth from chronux toolbox ().
function tmpdata=cal_stimulusresponse(tmpdata,eventtype)
for i=1:length(eventtype)
    try 
        data=eval(['tmpdata.',eventtype{i}]);
        data=cal_binspikefrom_LFPspike(data);
        data=cal_responsefrom_LFPspike(data);
        eval(['tmpdata.',eventtype{i},'=data;']);
    catch ME
        disp(ME)
    end
end
end
function structdata=cal_responsefrom_LFPspike(structdata)
% timeinitial is the begin time of each trial (-1)
% using the ttest to evaluate the responses to the given stimulus
spiketime=cellfun(@(x) x+structdata.EVTinfo.timerange(1),structdata.SPKdata,'UniformOutput',0);
spikebefore=cellfun(@(x) length(find(x<0&x>structdata.EVTinfo.timerange(1)))/abs(structdata.EVTinfo.timerange(1)),spiketime,'UniformOutput',1);
spikeafter=cellfun(@(x) length(find(x>0&x<structdata.EVTinfo.timerange(2)))/abs(structdata.EVTinfo.timerange(2)),spiketime,'UniformOutput',1);
for i=1:size(spikebefore,1)
[~,p]=ttest(spikebefore(i,:),spikeafter(i,:));
if mean(spikebefore(i,:))>mean(spikeafter(i,:)) && p<0.05
    response(i)=-1;
elseif mean(spikebefore(i,:))<mean(spikeafter(i,:)) && p<0.05
    response(i)=1;
else
    response(i)=0;
end
end
structdata.response=response;
end
function [cmatrix1,pmatrix,cmatrix2]=naiveBayers(dataset,category,modality)
% % % Bayers Classifier using the Spikeinformation and categroy
% Spike is a cell matrix with n cells (predictors), 
% each cell contains the trials * time bins, the categroy is the class information about each trials
% using the trials in each bins to classifier the predictors to which
% stimulus.
% error1 is the true predict, error2 is the shuffled predict

% limited the trials among the cells


%  for i=21:60
%     dataset=cell2mat(cellfun(@(x) x(:,i),binspike,'UniformOutput',0));
    mdl1=fitcnb(dataset',category,'ClassNames',modality,'DistributionNames','kernel','CrossVal','on','KFold',5);
    label1=kfoldPredict(mdl1);
    cmatrix1=confusionmat(category,label1);
    figure;
    confusionchart(cmatrix1,modality); 
    %cmatrix1=cmatrix1./sum(cmatrix1,2);
    parfor s=1:200
        mdl2=fitcnb(dataset',datasample(category,length(category),'replace',false),'ClassNames',modality,'DistributionNames','kernel','CrossVal','on','KFold',5);
        label2=kfoldPredict(mdl2);
        cmatrix2(:,:,s)=confusionmat(category,label2);
        % cmatrix2(:,:,s)=cmatrix./sum(cmatrix,2);
    end
    for i=1:size(cmatrix1,1)
       for j=1:size(cmatrix1,2)
           pmatrix(i,j)=sum(squeeze(cmatrix1(i,j)<cmatrix2(i,j,:))==1)/200;
       end
    end
%  end
end
function [accurate]=naiveBayers2(dataset,category,modality)
    mdl1=fitcnb(dataset',category,'ClassNames',modality,'DistributionNames','kernel','CrossVal','on','KFold',5);
    label1=kfoldPredict(mdl1);
    cmatrix1=confusionmat(category,label1);
    accurate=diag(cmatrix1)/15;
end