
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
                eff=cohensD(perfT(:,1),perfT(:,2));%% laser-on minus laser off
                vgat_heat_mat{tRowIdx,1}=sprintf('%s_%s',parts{pIdx},tasks{tIdx});
                vgat_heat_mat{tRowIdx,rIdx+1}=eff;
            end
        end
        if ~skipped && any(cellfun(@(x) ~isempty(x) && ~isnan(x),vgat_heat_mat(tRowIdx,2:end)))
            tRowIdx=tRowIdx+1;
        end
    end
end


if ~any(cellfun(@(x) ~isempty(x) && ~isnan(x),vgat_heat_mat(end,2:end)))
    vgat_heat_mat(end,:)=[];
end

%% unnecessarily complex sort of tasks
[~,simpleIdx]=sort(vgat_heat_mat(2:end,1));
vgat_heat_mat(2:end,:)=vgat_heat_mat(simpleIdx+1,:);
sortIdx=1:(size(vgat_heat_mat,1)-1);
sortTags=vgat_heat_mat(2:end,1);
sortIdx(~contains(sortTags,'ear'))=sortIdx(~contains(sortTags,'ear'))+1000;
sortIdx(contains(sortTags,'DPA_S'))=sortIdx(contains(sortTags,'DPA_S'))+100;
sortIdx(contains(sortTags,'Dual_DPA_D'))=sortIdx(contains(sortTags,'Dual_DPA_D'))+200;
sortIdx(contains(sortTags,'Dual_DRT'))=sortIdx(contains(sortTags,'Dual_DRT'))+300;
[~,tIdx]=sort(sortIdx);
vgat_heat_mat(2:end,:)=vgat_heat_mat(tIdx+1,:);

%%


% vgat_heat_matBak=vgat_heat_mat;
vgat_heat_mat(:, find(contains(vgat_heat_mat(1,2:end),'Ctrl'))+1)=[];
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
save('vgat_heat_mat.mat','vgat_heat_mat');

function out=cohensD(cond1,cond2)
nOpto=sum(~isnan(cond1))-1;
stdv=sqrt((nOpto*nanvar(cond1)+nOpto*nanvar(cond2))/(2*nOpto));
out=(nanmean(cond2)-nanmean(cond1))/stdv;
end


function out=translateName(wtRows,vgatLUT)
    for i=1:size(vgatLUT,1)
        if ~isempty(vgatLUT{i,2}) && ~strcmp(vgatLUT{i,1},vgatLUT{i,2})
            wtRows(strcmp(wtRows(:,2),vgatLUT{i,1}),2)={vgatLUT{i,2}};
        end
    end
    out=wtRows;
end
