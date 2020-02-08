%% laser-off then laser on
% sum_nphr_welltrained.m
cd('K:\Mapping\mapping');

%% Per Task Vgat
tasks={'DPA-ED','DPA-LD','DPA-DM'};
load('vgat_WT_LN.mat','wtRows','wtRowsHitMiss','wtRowsCRFalse');
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
[vgat_arr,vgat_regions]=oneOpto(vgat_branches,'-VGAT',tasks);
[nphr_arr,nphr_regions]=oneOpto(all_branches,'-NPHR',tasks);
value_labels={'ED-correct-rate','ED-hits','ED-rejects','LD-correct-rate','LD-hits','LD-rejects','DM-correct-rate','DM-hits','DM-rejects'};
save('zx_specificity_mat.mat','vgat_arr','nphr_arr','vgat_regions','nphr_regions','value_labels');


function [value_arr,regions]=oneOpto(all_tasks,suffix,tasks)


value_arr=[];
regions=cell(0);
for t=tasks
    perTask=all_tasks(strcmp(all_tasks(:,1),t),2:3);
    dir_sel=zeros(size(perTask,1),1);
    for i=1:size(perTask,1)
        dir_sel(i)=zx_s_r(perTask{i,2}(:,1),perTask{i,2}(:,2));
    end
    value_arr=vertcat(value_arr,dir_sel');
    regions{end+1}=perTask(:,1);
    [dir_sel,effIdx]=sort(dir_sel);
    fh=figure('Color','w','Position',[50,50,720,480]);
    hold on;
    bar(dir_sel,'FaceColor','w','EdgeColor','k');
    set(gca(),'XTick',1:size(dir_sel),'XTickLabel',perTask(effIdx,1),'XTickLabelRotation',90)
    ylabel('OpGen specificity')
    title([t,suffix]);
    print(sprintf('SPEC_%s_%s',char(t),char(suffix)),'-dpng','-r300')
    %% hit and rejection scatter
    
    rejects=zeros(size(perTask,1),1);
    hits=zeros(size(perTask,1),1);
    for i=1:size(perTask,1)
        hits(i)=zx_s_r(perTask{i,2}(:,3),perTask{i,2}(:,4));
        rejects(i)=zx_s_r(perTask{i,2}(:,5),perTask{i,2}(:,6));
    end
    value_arr=vertcat(value_arr,hits',rejects');
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
    title([t,suffix])
end

end
function dir_sel=zx_s_r(optoOff,optoOn)
prob_increased=nnz((optoOff ==1 & optoOn ==1) | (optoOn>optoOff));
prob_decreased=nnz((optoOff ==1 & optoOn ==1) | (optoOn<optoOff));
dir_sel=(prob_increased-prob_decreased)/(prob_increased+prob_decreased);
end