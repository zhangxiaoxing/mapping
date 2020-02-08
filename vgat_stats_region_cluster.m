load('vgat_heat_mat.mat','vgat_heat_mat','vgat_SEM','vgat_individual_eff');

allStats=cell(0,3);
idx=1;
for tId=[13,16,19]
    col_sel=~cellfun('isempty',vgat_individual_eff(tId,2:end));
    [impaired,improved]=statsOne(cell2mat(vgat_heat_mat(tId,[false,col_sel])),...
        cell2mat(vgat_SEM(tId,[false,col_sel])),...
        vgat_individual_eff(tId,[false,col_sel]),...
        vgat_heat_mat(1,[false,col_sel]),vgat_heat_mat{tId,1});
        allStats(end+1,:)={impaired,improved,vgat_heat_mat{tId,1}};
%     allStats{end,2}=tId;
end
save('cluster_by_opgen.mat','allStats');

function [statsImpair,statsImprove]=statsOne(effSizeMean,effSizeSEM,effSize,region_list,task)
%     if contains(task,'Simple_DPA_')
%         task=replace(task,'Simple_DPA_','');
%     end
    [effSizeMean,effIdx]=sort(effSizeMean);
%     fh=figure('Color','w','Position',[50,50,720,480]);
%     hold on;
%     bar(effSizeMean,'FaceColor','w','EdgeColor','k');
%     errorbar(effSizeMean,effSizeSEM(effIdx),'.','LineWidth',1,'Color',[0.5,0.5,0.5]);

    effSize=effSize(effIdx);
    region_list=region_list(effIdx);
    ctrlIdx=contains(region_list,'ctrl','IgnoreCase',true);
    
    ctrlStats=cell2mat(effSize(ctrlIdx));
    
    statsImpair=cell(0,4);
    
    for i=1:find(ctrlIdx)
        region_cluster=effSize{i};
        if mean(region_cluster)<mean(ctrlStats)
           p=ranksum(region_cluster,ctrlStats);
           if p<0.2
               statsImpair(end+1,:)={region_list{i},find(ctrlIdx),mean(region_cluster),p};
           end
        else
            break
        end
    end
    
    statsImprove=cell(0,4);
    for i=numel(effSizeMean):-1:(find(ctrlIdx)+1)
        region_cluster=effSize{i};
        if mean(region_cluster)>mean(ctrlStats)
           p=ranksum(region_cluster,ctrlStats);
           if p<0.2
               statsImprove(end+1,:)={region_list{i},find(ctrlIdx),mean(region_cluster),p};
           end
        else
            break
        end
    end
end