
cd('K:\Mapping\mapping')

load('vgat_WT_LN.mat','wtRows','lnRows','wtRowsHitMiss','wtRowsCRFalse');
load('vgatLUT.mat','vgatLUT')
%% WT

for i=1:size(wtRows,1)
    if ~isempty(wtRows{i,8})
        wtRows(end+1,1:4)=wtRows(i,1:4);
        wtRows(end,5:7)=wtRows(i,8:10);
    end
    if ~isempty(wtRows{i,11})
        wtRows(end+1,1:4)=wtRows(i,1:4);
        wtRows(end,5:7)=wtRows(i,11:13);
    end
end
wtRows(:,8:end)=[];

%WT hit-miss and CR-false
wtRowsHitMiss(cellfun('isempty',wtRowsHitMiss(:,5)),:)=[];
wtRowsHitMiss(:,1)=cellfun(@(x) sprintf('%s-Hit',x),wtRowsHitMiss(:,1),'UniformOutput',false);

wtRowsCRFalse(cellfun('isempty',wtRowsCRFalse(:,5)),:)=[];
wtRowsCRFalse(:,1)=cellfun(@(x) sprintf('%s-CR',x),wtRowsCRFalse(:,1),'UniformOutput',false);

wtRows=[wtRows;wtRowsHitMiss;wtRowsCRFalse];

%% Learning

for i=1:size(lnRows,1)
lnRows{i,1}=sprintf('%s-%s',lnRows{i,1},lnRows{i,5});
end
lnRows(:,5)=[];
for i=1:size(lnRows,1)
    if ~isempty(lnRows{i,8})
        lnRows(end+1,1:4)=lnRows(i,1:4);
        lnRows(end,5:7)=lnRows(i,8:10);
    end
    if ~isempty(lnRows{i,11})
        lnRows(end+1,1:4)=lnRows(i,1:4);
        lnRows(end,5:7)=lnRows(i,11:13);
    end
end
lnRows(:,8:end)=[];

wtRows=[wtRows;lnRows];
wtRows(cellfun('isempty',wtRows(:,6)),:)=[];

%% calculate learning baseline
allbase=wtRows(contains(wtRows(:,1),'DPA-Le') & contains(wtRows(:,2),'Ct'),:);
dpaCtrlDayAvg=arrayfun(@(i) nanmean([allbase{contains(allbase(:,1),num2str(i)), 7}]),1:5);
for i=1:5
    wtRows(contains(wtRows(:,1),'DPA-Le') & contains(wtRows(:,1),num2str(i)) & isnan([wtRows{:,6}])',6)={dpaCtrlDayAvg(i)};
end 
%% continue    

wtRows(isnan(cell2mat(wtRows(:,7))),:)=[];
wtRows(isnan(cell2mat(wtRows(:,6))),:)=[];

wtRows=translateName(wtRows,vgatLUT);

regions=unique(wtRows(:,2));
tasks=unique(wtRows(:,1));
parts=unique(wtRows(:,5));
vgat_heat_mat=cell(0);
vgat_SEM=cell(0);
vgat_individual_eff=cell(0);
tRowIdx=2;
skipped=false;

for tIdx=1:length(tasks)
    for pIdx=1:length(parts)
        for rIdx=1:length(regions)
            if tRowIdx==2
                vgat_heat_mat{1,rIdx+1}=regions{rIdx};
            end
            %%%%%col 5
%             if contains(tasks(tIdx),'DualTask-learning-day3-catch') && pIdx==2 && rIdx==1
%                 disp('day3cat');
%             end
            perfT=cell2mat(wtRows(strcmp(wtRows(:,1),tasks(tIdx))...
                & strcmp(wtRows(:,2),regions(rIdx))...
                & strcmp(wtRows(:,5),parts(pIdx)),6:7));
            if isempty(perfT) || all(isnan(perfT(:)))
                 vgat_heat_mat{tRowIdx,rIdx+1}=[];
            else
                skipped=false;
                [effMean,effSEM,effAll]=cohen_s_d(perfT(:,1),perfT(:,2));%% laser-on minus laser off
                vgat_heat_mat{tRowIdx,1}=sprintf('%s_%s',parts{pIdx},tasks{tIdx});
                vgat_heat_mat{tRowIdx,rIdx+1}=effMean;
                vgat_SEM{tRowIdx,rIdx+1}=effSEM;
                vgat_individual_eff{tRowIdx,rIdx+1}=effAll;
            end
        end
        if ~skipped && any(cellfun(@(x) ~isempty(x) && ~isnan(x),vgat_heat_mat(tRowIdx,2:end)))
            tRowIdx=tRowIdx+1;
        end
    end
end


if ~any(cellfun(@(x) ~isempty(x) && ~isnan(x),vgat_heat_mat(end,2:end)))
    vgat_heat_mat(end,:)=[];
    vgat_SEM(end,:)=[];
    vgat_individual_eff(end,:)=[];
end

%% unnecessarily complex sort of tasks
[~,simpleIdx]=sort(vgat_heat_mat(2:end,1));
vgat_heat_mat(2:end,:)=vgat_heat_mat(simpleIdx+1,:);
vgat_SEM(2:end,:)=vgat_SEM(simpleIdx+1,:);
vgat_individual_eff(2:end,:)=vgat_individual_eff(simpleIdx+1,:);


sortIdx=1:(size(vgat_heat_mat,1)-1);
sortTags=vgat_heat_mat(2:end,1);
sortIdx(~contains(sortTags,'ear'))=sortIdx(~contains(sortTags,'ear'))+1000;
sortIdx(contains(sortTags,'DPA_S'))=sortIdx(contains(sortTags,'DPA_S'))+100;
sortIdx(contains(sortTags,'Dual_DPA_D'))=sortIdx(contains(sortTags,'Dual_DPA_D'))+200;
sortIdx(contains(sortTags,'Dual_DRT'))=sortIdx(contains(sortTags,'Dual_DRT'))+300;
[~,tIdx]=sort(sortIdx);
vgat_heat_mat(2:end,:)=vgat_heat_mat(tIdx+1,:);
vgat_SEM(2:end,:)=vgat_SEM(tIdx+1,:);
vgat_individual_eff(2:end,:)=vgat_individual_eff(tIdx+1,:);
%%


% vgat_heat_matBak=vgat_heat_mat;
% vgat_heat_mat(:, find(contains(vgat_heat_mat(1,2:end),'Ctrl'))+1)=[];
vgat_heat_mat(cellfun('isempty',vgat_heat_mat))={nan};
regions=vgat_heat_mat(1,2:end);
imgMat=cell2mat(vgat_heat_mat(2:end,2:end));
% [~,rIdx]=sort(nanmean(imgMat([1:11,13,14],:)));
[~,rIdx]=sort(imgMat(12,:));
imgMat=imgMat(:,rIdx);

close all;
figure('Color','w','Position',[1441,-93,1280,917]);
hi=imagesc(imgMat,[-1.5,1.5]);
set(hi,'AlphaData',~isnan(imgMat));
set(gca,'Color',[0.5,0.5,0.5])
colormap('jet');
colorbar();
set(gca,'XTick',1:size(imgMat,2),'XTickLabel',regions(rIdx),'XTickLabelRotation',90,...
'YTick',1:size(imgMat,1),'YTickLabel',vgat_heat_mat(2:end,1),'TickLabelInterpreter','none','FontSize',14);
cd('K:\Mapping\mapping')
save('vgat_heat_mat.mat','vgat_heat_mat','vgat_SEM','vgat_individual_eff');

for tId=[13,16,19]
    col_sel=~cellfun('isempty',vgat_individual_eff(tId,2:end));
    plotOne(cell2mat(vgat_heat_mat(tId,[false,col_sel])),...
        cell2mat(vgat_SEM(tId,[false,col_sel])),...
        vgat_individual_eff(tId,[false,col_sel]),...
        vgat_heat_mat(1,[false,col_sel]),vgat_heat_mat{tId,1});
end





function [effOpto_mean,sem,effSize]=cohen_s_d(optoOff,optoOn)
nOpto=numel(optoOff)-1;
sOpto=sqrt((nOpto*var(optoOff)+nOpto*var(optoOn))/(2*nOpto));
effOpto=(optoOn-optoOff)/sOpto;
effOpto_mean=mean(effOpto);
sem=std(effOpto)/sqrt(length(effOpto));
effSize=effOpto;
end


function out=translateName(wtRows,vgatLUT)
    for i=1:size(vgatLUT,1)
        if ~isempty(vgatLUT{i,2}) && ~strcmp(vgatLUT{i,1},vgatLUT{i,2})
            wtRows(strcmp(wtRows(:,2),vgatLUT{i,1}),2)={vgatLUT{i,2}};
        end
    end
    out=wtRows;
end

function plotOne(effSizeMean,effSizeSEM,effSize,region_list,task)
    if contains(task,'Simple_DPA_')
        task=replace(task,'Simple_DPA_','');
    end
    [effSizeMean,effIdx]=sort(effSizeMean);
    fh=figure('Color','w','Position',[50,50,720,480]);
    hold on;
    bar(effSizeMean,'FaceColor','w','EdgeColor','k');
    errorbar(effSizeMean,effSizeSEM(effIdx),'.','LineWidth',1,'Color',[0.5,0.5,0.5]);
    effSize=effSize(effIdx);
    for i=1:length(effSize)
        plot(i,effSize{i},'b.')
    end
    set(gca(),'XTick',1:numel(region_list),'XTickLabel',region_list(effIdx),'XTickLabelRotation',90);
    ylabel('Effect size')
    title([task, '-VGAT']);
    print(sprintf('EFFSIZE_%s_%s',task,'VGAT'),'-dpng','-r300')
    
end