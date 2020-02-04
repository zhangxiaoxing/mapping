%% laser-off then laser on
cd('K:\Mapping\codes');

%% Per Task Vgat
tasks={'DPA-ED','DPA-LD','DPA-DM'};
load('WT_LN.mat','wtRows','wtRowsHitMiss','wtRowsCRFalse');
wtRows=wtRows(strcmpi(wtRows(:,5),'simple_dpa'),:);
regionSet=unique(wtRows(:,2))';
vgat_branches=cell(0,3);
for t=tasks
    for r=regionSet
        perf_correct=cell2mat(wtRows(strcmpi(wtRows(:,1),t) & strcmpi(wtRows(:,2),r),6:7));
        hits=cell2mat(wtRowsHitMiss(strcmpi(wtRowsHitMiss(:,1),t) & strcmpi(wtRowsHitMiss(:,2),r),6:7));
        rejects=cell2mat(wtRowsCRFalse(strcmpi(wtRowsCRFalse(:,1),t) & strcmpi(wtRowsCRFalse(:,2),r),6:7));
        if numel(unique([size(perf_correct,1),size(hits,1),size(rejects,1)]))==1
            vgat_branches(end+1,:)={char(t),char(r),horzcat(perf_correct,hits,rejects)};
        else
            fprintf('inconsistency, %s, %s\n',char(t),char(r));
        end
    end
end


%% Per Task nphr

load('nphr_perf.mat','all_branches');
for t=tasks
    perTask=all_branches(strcmp(all_branches(:,1),t),2:3);
    vgatTask=vgat_branches(strcmp(vgat_branches(:,1),t),2:3);
    labels=vertcat(cellfun(@(x) [x,'-nphr'],perTask(:,1),'UniformOutput',false),...
                   cellfun(@(x) [x,'-vgat'],vgatTask(:,1),'UniformOutput',false)) ;
    perTask=horzcat(labels,vertcat(perTask(:,2),vgatTask(:,2)));
    dir_sel=zeros(size(perTask,1),1);
    for i=1:size(perTask,1)
        dir_sel(i)=zx_s_r(perTask{i,2}(:,1),perTask{i,2}(:,2));
    end
    [dir_sel,effIdx]=sort(dir_sel);
    fh=figure('Color','w');
    hold on;
    bar(dir_sel,'FaceColor','w','EdgeColor','k');
    set(gca(),'XTick',1:size(labels),'XTickLabel',labels(effIdx),'XTickLabelRotation',90)
    title(t);
    %% hit and rejection scatter
    
    rejects=zeros(size(perTask,1),1);
    hits=zeros(size(perTask,1),1);
    for i=1:size(perTask,1)
        hits(i)=zx_s_r(perTask{i,2}(:,3),perTask{i,2}(:,4));
        rejects(i)=zx_s_r(perTask{i,2}(:,5),perTask{i,2}(:,6));
    end    
    
    colors={'r','g','b','c','m','k'};
    markers={'+','o','x','s','d','^','p'}; 
    
    figure('Color','w')
    hold on
    for i=1:length(rejects)
%         colorIdx=rem(i,6)+1;
%         markIdx=fix(i/6)+1;
%         markStr=sprintf('%s%s',colors{colorIdx},markers{markIdx});
        if strcmpi(perTask(i,1),'ctrl')
            color='r';
        else
            color='k';
        end
        text(hits(i),rejects(i),perTask(i,1),'Color',color)
        
    end
    xlim([-1.2,1.2]);
    ylim([-1.2,1.2]);
%     legend(perTask(:,1),'Location','eastoutside')
    xlabel('effect on hits');
    ylabel('effect on rejection');
    title(t)
    
    
    
    
    
    
end
function dir_sel=zx_s_r(optoOff,optoOn)
prob_increased=nnz((optoOff ==1 & optoOn ==1) | (optoOn>optoOff));
prob_decreased=nnz((optoOff ==1 & optoOn ==1) | (optoOn<optoOff));
dir_sel=(prob_increased-prob_decreased)/(prob_increased+prob_decreased);
end