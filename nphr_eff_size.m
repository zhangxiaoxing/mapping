%% laser-off then laser on
%% sum_nphr_welltrained.m

cd('D:\Mapping\mapping');
if ~exist('all_branches','var')
    load('nphr_perf.mat','all_branches')
end
tasks={'DPA-ED','DPA-LD','DPA-DM'};
load('nphr_mat.mat','nphr_lut');
mean_arr=[];
median_arr=[];

for t=tasks
    perTask=all_branches(strcmp(all_branches(:,1),t),2:3);
    effSizeMean=zeros(size(perTask,1),1);
    effSizeMedian=zeros(size(perTask,1),1);
    effSizeSEM=zeros(size(perTask,1),1);
    effSize=cell(size(perTask,1),1);
    for i=1:size(perTask,1)
        [effSizeMean(i),effSizeMedian(i),effSizeSEM(i),effSize{i}]=cohen_s_d(perTask{i,2}(:,1),perTask{i,2}(:,2));
    end
    mean_arr=vertcat(mean_arr,effSizeMean');
    median_arr=vertcat(median_arr,effSizeMedian');
    [effSizeMean,effIdx]=sort(effSizeMean);
    effSizeMedian=effSizeMedian(effIdx);
    effSize=effSize(effIdx);

    fh=figure('Color','w','Position',[50,50,720,480]);
    hold on;
%     bar(effSizeMean,'FaceColor','w','EdgeColor','k');
    yline(0,'k-');    
    for i=1:numel(effSizeMean)
        plot([i-0.4,i+0.4],effSizeMean(i)*[1,1],'b-','LineWidth',2);
        plot([i-0.4,i+0.4],effSizeMedian(i)*[1,1],'r-','LineWidth',2);
        plot(i,effSize{i},'k.')
    end
%     errorbar(effSizeMean,effSizeSEM(effIdx),'.','LineWidth',1,'Color',[0.5,0.5,0.5]);
    
    set(gca(),'XTick',1:size(perTask,1),'XTickLabel',perTask(effIdx,1),'XTickLabelRotation',90);
    ylabel('Effect size')
    print(sprintf('EFFSIZE_%s_%s',char(t),'NPHR'),'-dpng','-r300')
    title([t, '-NPHR']);
    %% hit and rejection scatter
    
    rejectsEffMean=zeros(size(perTask,1),1);
    hitsEffMean=zeros(size(perTask,1),1);
    rejectsEffMedian=zeros(size(perTask,1),1);
    hitsEffMedian=zeros(size(perTask,1),1);
    rejects=cell(size(perTask,1),1);
    hits=cell(size(perTask,1),1);
    %[correct_off,correct_on,hit_off,hit_on,reject_off,reject_on];
    for i=1:size(perTask,1)
        [rejectsEffMean(i),rejectsEffMedian(i),~,rejects{i}]=cohen_s_d(perTask{i,2}(:,5),perTask{i,2}(:,6));
        [hitsEffMean(i),hitsEffMedian(i),~,hits{i}]=cohen_s_d(perTask{i,2}(:,3),perTask{i,2}(:,4));
    end    
    mean_arr=vertcat(mean_arr,hitsEffMean',rejectsEffMean'); %%
    median_arr=vertcat(median_arr,hitsEffMedian',rejectsEffMedian'); %%
    colors={'r','g','b','c','m','y','k'};
    markers={'+','o','x','s','d','^','p','h'}; 
    
    figure('Color','w')
    hold on
    for i=1:length(rejects)
        colorIdx=rem(i,7)+1;
        markIdx=fix(i/7)+1;
        markStr=sprintf('%s%s',colors{colorIdx},markers{markIdx});
        plot(hits{i},rejects{i},markStr,'MarkerSize',5,'MarkerFaceColor',colors{colorIdx},'LineWidth',2)
        
    end
    legend(perTask(:,1),'Location','eastoutside')
    xlabel('effect on hits');
    ylabel('effect on rejection');
    title(t)
  
end


cd('D:\Mapping\mapping');
value_labels={'ED-correct-rate','ED-hits','ED-rejects','LD-correct-rate','LD-hits','LD-rejects','DM-correct-rate','DM-hits','DM-rejects'};
region_list=perTask(:,1);
save('nphr_mat.mat','mean_arr','median_arr','value_labels','region_list','nphr_lut');

function [effOpto_mean,effOpto_median,sem,effSize]=cohen_s_d(optoOff,optoOn)
nOpto=numel(optoOff)-1;
sOpto=sqrt((nOpto*var(optoOff)+nOpto*var(optoOn))/(2*nOpto));
effOpto=(optoOn-optoOff)/sOpto;
effOpto_mean=mean(effOpto);
effOpto_median=median(effOpto);
sem=std(effOpto)/sqrt(length(effOpto));
effSize=effOpto;
end