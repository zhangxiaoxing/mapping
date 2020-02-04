%% laser-off then laser on
%% sum_nphr_welltrained.m

cd('K:\Mapping\mapping');

tasks={'DPA-ED','DPA-LD','DPA-DM'};

value_arr=[];

for t=tasks
    perTask=all_branches(strcmp(all_branches(:,1),t),2:3);
    effSizeMean=zeros(size(perTask,1),1);
    effSizeSEM=zeros(size(perTask,1),1);
    effSize=cell(size(perTask,1),1);
    for i=1:size(perTask,1)
        [effSizeMean(i),effSizeSEM(i),effSize{i}]=cohen_s_d(perTask{i,2}(:,1),perTask{i,2}(:,2));
    end
    value_arr=vertcat(value_arr,effSizeMean');
    [effSizeMean,effIdx]=sort(effSizeMean);
    fh=figure('Color','w');
    hold on;
    bar(effSizeMean,'FaceColor','w','EdgeColor','k');
    errorbar(effSizeMean,effSizeSEM(effIdx),'k.','LineWidth',1);
    effSize=effSize(effIdx);
    for i=1:length(effSize)
        plot(i,effSize{i},'b.')
    end
    set(gca(),'XTick',1:size(perTask,1),'XTickLabel',perTask(effIdx,1),'XTickLabelRotation',90)
    title(t);
    %% hit and rejection scatter
    
    rejectsEffMean=zeros(size(perTask,1),1);
    hitsEffMean=zeros(size(perTask,1),1);
    rejects=cell(size(perTask,1),1);
    hits=cell(size(perTask,1),1);
    %[correct_off,correct_on,hit_off,hit_on,reject_off,reject_on];
    for i=1:size(perTask,1)
        [rejectsEffMean(i),~,rejects{i}]=cohen_s_d(perTask{i,2}(:,5),perTask{i,2}(:,6));
        [hitsEffMean(i),~,hits{i}]=cohen_s_d(perTask{i,2}(:,3),perTask{i,2}(:,4));
    end    
    value_arr=vertcat(value_arr,hitsEffMean',rejectsEffMean'); %%
    colors={'r','g','b','c','m','k'};
    markers={'+','o','x','s','d','^','p'}; 
    
    figure('Color','w')
    hold on
    for i=1:length(rejects)
        colorIdx=rem(i,6)+1;
        markIdx=fix(i/6)+1;
        markStr=sprintf('%s%s',colors{colorIdx},markers{markIdx});
        plot(hits{i},rejects{i},markStr,'MarkerSize',10,'MarkerFaceColor',colors{colorIdx},'LineWidth',2)
        
    end
    legend(perTask(:,1),'Location','eastoutside')
    xlabel('effect on hits');
    ylabel('effect on rejection');
    title(t)
  
end


cd('K:\Mapping\mapping');
value_labels={'ED-correct-rate','ED-hits','ED-rejects','LD-correct-rate','LD-hits','LD-rejects','DM-correct-rate','DM-hits','DM-rejects'};
region_list=perTask(:,1);
save('nphr_mat.mat','value_arr','value_labels','region_list');

function [effOpto_mean,sem,effSize]=cohen_s_d(optoOff,optoOn)
nOpto=numel(optoOff)-1;
sOpto=sqrt((nOpto*var(optoOff)+nOpto*var(optoOn))/(2*nOpto));
effOpto=(optoOn-optoOff)/sOpto;
effOpto_mean=mean(effOpto);
sem=std(effOpto)/sqrt(length(effOpto));
effSize=effOpto;
end