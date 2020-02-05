path='k:\code\delayStats.csv';
delayStats=importdata(path,',',0);

path='k:\code\testStats.csv';
testStats=importdata(path,',',0);

path='K:\code\tca\Batch30B\tca_su_factors.csv';
tcaStats=importdata(path,',',0);
tcaStats.raw_data=tcaStats.data;

load('vgatLUT.mat','vgatLUT');
load('opgene_ephys_corr.mat','translate')
nphrLUT=translate;

vfstr=load('vgat_heat_mat.mat','vgat_heat_mat');
nfstr=load('nphr_mat.mat','value_arr','value_labels','region_list');

figIdx=20020600;

%% nphr, vgat tca components supp

tcaStats.data=tcaStats.raw_data(:,[3,3]);
plotOneCorr(nfstr.value_arr(7,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR decision-making correct-rate effect-size',...
    'TCA component 3 conefficient') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(9,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR decision-making reject effect-size',...
    'TCA component 3 conefficient') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(8,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR decision-making hit effect-size',...
    'TCA component 3 conefficient') %nphr_eff_size.m

plotOneCorr(vgat_arr(7,:),vgatLUT,vgat_regions{3},tcaStats,...
    'VGAT-ChR2 decision-making correct-rate OpGen specificity',...
    'TCA component 3 conefficient') %nphr_model.m

plotOneCorr(vgat_arr(8,:),vgatLUT,vgat_regions{3},tcaStats,...
    'VGAT-ChR2 decision-making hit OpGen specificity',...
    'TCA component 3 conefficient') %nphr_model.m

plotOneCorr(vgat_arr(9,:),vgatLUT,vgat_regions{3},tcaStats,...
    'VGAT-ChR2 decision-making rejection OpGen specificity',...
    'TCA component 3 conefficient') %nphr_model.m


%% component 4
tcaStats.data=tcaStats.raw_data(:,[4,4]);
plotOneCorr(nfstr.value_arr(7,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR decision-making correct-rate effect-size',...
    'TCA component 4 conefficient') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(9,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR decision-making reject effect-size',...
    'TCA component 4 conefficient') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(8,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR decision-making hit effect-size',...
    'TCA component 4 conefficient') %nphr_eff_size.m


plotOneCorr(vgat_arr(7,:),vgatLUT,vgat_regions{3},tcaStats,...
    'VGAT-ChR2 decision-making correct-rate OpGen specificity',...
    'TCA component 4 conefficient') %nphr_model.m

plotOneCorr(vgat_arr(8,:),vgatLUT,vgat_regions{3},tcaStats,...
    'VGAT-ChR2 decision-making hit OpGen specificity',...
    'TCA component 4 conefficient') %nphr_model.m

plotOneCorr(vgat_arr(9,:),vgatLUT,vgat_regions{3},tcaStats,...
    'VGAT-ChR2 decision-making rejection OpGen specificity',...
    'TCA component 4 conefficient') %nphr_model.m

%% component 6
tcaStats.data=tcaStats.raw_data(:,[6,6]);
plotOneCorr(nfstr.value_arr(7,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR decision-making correct-rate effect-size',...
    'TCA component 6 conefficient') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(9,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR decision-making reject effect-size',...
    'TCA component 6 conefficient') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(8,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR decision-making hit effect-size',...
    'TCA component 6 conefficient') %nphr_eff_size.m


plotOneCorr(vgat_arr(7,:),vgatLUT,vgat_regions{3},tcaStats,...
    'VGAT-ChR2 decision-making correct-rate OpGen specificity',...
    'TCA component 6 conefficient') %nphr_model.m

plotOneCorr(vgat_arr(8,:),vgatLUT,vgat_regions{3},tcaStats,...
    'VGAT-ChR2 decision-making hit OpGen specificity',...
    'TCA component 6 conefficient') %nphr_model.m

plotOneCorr(vgat_arr(9,:),vgatLUT,vgat_regions{3},tcaStats,...
    'VGAT-ChR2 decision-making rejection OpGen specificity',...
    'TCA component 6 conefficient') %nphr_model.m


%% supp end
return


%% nphr, vgat-tca components 3

tcaStats.data=tcaStats.raw_data(:,[3,3]);
plotOneCorr(nphr_arr(7,:),nphrLUT,nphr_regions{3},tcaStats,...
    'NPHR decision-making correct-rate OpGen specificity',...
    'TCA component 3 conefficient') %nphr_model.m

plotOneCorr(nphr_arr(8,:),nphrLUT,nphr_regions{3},tcaStats,...
    'NPHR decision-making hit OpGen specificity',...
    'TCA component 3 conefficient') %nphr_model.m

plotOneCorr(nphr_arr(9,:),nphrLUT,nphr_regions{3},tcaStats,...
    'NPHR decision-making rejection OpGen specificity',...
    'TCA component 3 conefficient') %nphr_model.m


plotOneCorr(cell2mat(vfstr.vgat_heat_mat(13,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 decision-making correct-rate effect-size',...
    'TCA component 3 conefficient') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(14,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 decision-making correct-rejection effect-size',...
    'TCA component 3 conefficient') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(15,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 decision-making hit effect-size',...
    'TCA component 3 conefficient') %vgat_heatmap.m


%% component 4
tcaStats.data=tcaStats.raw_data(:,[4,4]);
plotOneCorr(nphr_arr(7,:),nphrLUT,nphr_regions{3},tcaStats,...
    'NPHR decision-making correct-rate OpGen specificity',...
    'TCA component 4 conefficient') %nphr_model.m

plotOneCorr(nphr_arr(8,:),nphrLUT,nphr_regions{3},tcaStats,...
    'NPHR decision-making hit OpGen specificity',...
    'TCA component 4 conefficient') %nphr_model.m

plotOneCorr(nphr_arr(9,:),nphrLUT,nphr_regions{3},tcaStats,...
    'NPHR decision-making rejection OpGen specificity',...
    'TCA component 4 conefficient') %nphr_model.m


plotOneCorr(cell2mat(vfstr.vgat_heat_mat(13,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 decision-making correct-rate effect-size',...
    'TCA component 4 conefficient') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(14,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 decision-making correct-rejection effect-size',...
    'TCA component 4 conefficient') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(15,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 decision-making hit effect-size',...
    'TCA component 4 conefficient') %vgat_heatmap.m

%% component 6
tcaStats.data=tcaStats.raw_data(:,[6,6]);
plotOneCorr(nphr_arr(7,:),nphrLUT,nphr_regions{3},tcaStats,...
    'NPHR decision-making correct-rate OpGen specificity',...
    'TCA component 6 conefficient') %nphr_model.m

plotOneCorr(nphr_arr(8,:),nphrLUT,nphr_regions{3},tcaStats,...
    'NPHR decision-making hit OpGen specificity',...
    'TCA component 6 conefficient') %nphr_model.m

plotOneCorr(nphr_arr(9,:),nphrLUT,nphr_regions{3},tcaStats,...
    'NPHR decision-making rejection OpGen specificity',...
    'TCA component 6 conefficient') %nphr_model.m


plotOneCorr(cell2mat(vfstr.vgat_heat_mat(13,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 decision-making correct-rate effect-size',...
    'TCA component 6 conefficient') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(14,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 decision-making correct-rejection effect-size',...
    'TCA component 6 conefficient') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(15,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 decision-making hit effect-size',...
    'TCA component 6 conefficient') %vgat_heatmap.m



%% nphr vgat early delay tca component 7

tcaStats.data=tcaStats.raw_data(:,[7,7]);
plotOneCorr(nfstr.value_arr(1,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR early-delay correct-rate effect-size',...
    'TCA component 7 conefficient') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(3,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR early-delay reject effect-size',...
    'TCA component 7 conefficient') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(2,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR early-delay hits effect-size',...
    'TCA component 7 conefficient') %nphr_eff_size.m


plotOneCorr(cell2mat(vfstr.vgat_heat_mat(16,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 early-delay correct-rate effect-size',...
    'TCA component 7 conefficient') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(17,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 early-delay correct-rejection effect-size',...
    'TCA component 7 conefficient') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(18,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 early-delay hit effect-size',...
    'TCA component 7 conefficient') %vgat_heatmap.m


%% nphr vgat early delay tca component 8
tcaStats.data=tcaStats.raw_data(:,[8,8]);

plotOneCorr(nfstr.value_arr(4,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR late-delay correct-rate effect-size',...
    'TCA component 8 conefficient') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(6,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR late-delay reject effect-size',...
    'TCA component 8 conefficient') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(5,:),nphrLUT,nfstr.region_list,tcaStats,...
    'NPHR late-delay hit effect-size',...
    'TCA component 8 conefficient') %nphr_eff_size.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(19,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 late-delay correct-rate effect-size',...
    'TCA component 8 conefficient') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(20,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 late-delay correct-rejection effect-size',...
    'TCA component 8 conefficient') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(21,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),tcaStats,...
    'VGAT-ChR2 late-delay hit effect-size',...
    'TCA component 8 conefficient') %vgat_heatmap.m



%% nphr - EFF Size
plotOneCorr(nfstr.value_arr(1,:),nphrLUT,nfstr.region_list,delayStats,...
    'NPHR early-delay correct-rate effect-size',...
    'Sample-selective fraction during delay') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(3,:),nphrLUT,nfstr.region_list,delayStats,...
    'NPHR early-delay reject effect-size',...
    'Sample-selective fraction during delay') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(2,:),nphrLUT,nfstr.region_list,delayStats,...
    'NPHR early-delay hits effect-size',...
    'Sample-selective fraction during delay') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(4,:),nphrLUT,nfstr.region_list,delayStats,...
    'NPHR late-delay correct-rate effect-size',...
    'Sample-selective fraction during delay') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(6,:),nphrLUT,nfstr.region_list,delayStats,...
    'NPHR late-delay reject effect-size',...
    'Sample-selective fraction during delay') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(5,:),nphrLUT,nfstr.region_list,delayStats,...
    'NPHR late-delay hit effect-size',...
    'Sample-selective fraction during delay') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(7,:),nphrLUT,nfstr.region_list,testStats,...
    'NPHR decision-making correct-rate effect-size',...
    'Pair-selective fraction during decision-making') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(9,:),nphrLUT,nfstr.region_list,testStats,...
    'NPHR decision-making reject effect-size',...
    'Pair-selective fraction during decision-making') %nphr_eff_size.m

plotOneCorr(nfstr.value_arr(8,:),nphrLUT,nfstr.region_list,testStats,...
    'NPHR decision-making hit effect-size',...
    'Pair-selective fraction during decision-making') %nphr_eff_size.m

%% VGAT - EFF Size
plotOneCorr(cell2mat(vfstr.vgat_heat_mat(13,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),testStats,...
    'VGAT-ChR2 decision-making correct-rate effect-size',...
    'Pair-selective fraction during decision-making') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(14,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),testStats,...
    'VGAT-ChR2 decision-making correct-rejection effect-size',...
    'Pair-selective fraction during decision-making') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(15,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),testStats,...
    'VGAT-ChR2 decision-making hit effect-size',...
    'Pair-selective fraction during decision-making') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(16,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),delayStats,...
    'VGAT-ChR2 early-delay correct-rate effect-size',...
    'Sample-selective fraction during delay') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(17,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),delayStats,...
    'VGAT-ChR2 early-delay correct-rejection effect-size',...
    'Sample-selective fraction during delay') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(18,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),delayStats,...
    'VGAT-ChR2 early-delay hit effect-size',...
    'Sample-selective fraction during delay') %vgat_heatmap.m


plotOneCorr(cell2mat(vfstr.vgat_heat_mat(19,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),delayStats,...
    'VGAT-ChR2 late-delay correct-rate effect-size',...
    'Sample-selective fraction during delay') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(20,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),delayStats,...
    'VGAT-ChR2 late-delay correct-rejection effect-size',...
    'Sample-selective fraction during delay') %vgat_heatmap.m

plotOneCorr(cell2mat(vfstr.vgat_heat_mat(21,2:end)),vgatLUT,vfstr.vgat_heat_mat(1,2:end),delayStats,...
    'VGAT-ChR2 late-delay hit effect-size',...
    'Sample-selective fraction during delay') %vgat_heatmap.m

%% VGAT - specificity

load('zx_specificity_mat.mat','vgat_arr','nphr_arr','vgat_regions','nphr_regions','value_labels');

plotOneCorr(vgat_arr(1,:),vgatLUT,vgat_regions{1},delayStats,...
    'VGAT-ChR2 early-delay correct-rate OpGen specificity',...
    'Sample-selective fraction during delay') %nphr_model.m
plotOneCorr(vgat_arr(2,:),vgatLUT,vgat_regions{1},delayStats,...
    'VGAT-ChR2 early-delay hit OpGen specificity',...
    'Sample-selective fraction during delay') %nphr_model.m
plotOneCorr(vgat_arr(3,:),vgatLUT,vgat_regions{1},delayStats,...
    'VGAT-ChR2 early-delay rejection OpGen specificity',...
    'Sample-selective fraction during delay') %nphr_model.m

plotOneCorr(vgat_arr(4,:),vgatLUT,vgat_regions{2},delayStats,...
    'VGAT-ChR2 late-delay correct-rate OpGen specificity',...
    'Sample-selective fraction during delay') %nphr_model.m

plotOneCorr(vgat_arr(5,:),vgatLUT,vgat_regions{2},delayStats,...
    'VGAT-ChR2 late-delay hit OpGen specificity',...
    'Sample-selective fraction during delay') %nphr_model.m

plotOneCorr(vgat_arr(6,:),vgatLUT,vgat_regions{2},delayStats,...
    'VGAT-ChR2 late-delay rejection OpGen specificity',...
    'Sample-selective fraction during delay') %nphr_model.m

plotOneCorr(vgat_arr(7,:),vgatLUT,vgat_regions{3},testStats,...
    'VGAT-ChR2 decision-making correct-rate OpGen specificity',...
    'Pair-selective fraction during decision making') %nphr_model.m

plotOneCorr(vgat_arr(8,:),vgatLUT,vgat_regions{3},testStats,...
    'VGAT-ChR2 decision-making hit OpGen specificity',...
    'Pair-selective fraction during decision making') %nphr_model.m

plotOneCorr(vgat_arr(9,:),vgatLUT,vgat_regions{3},testStats,...
    'VGAT-ChR2 decision-making rejection OpGen specificity',...
    'Pair-selective fraction during decision making') %nphr_model.m


%% NPHR - specificity
plotOneCorr(nphr_arr(1,:),nphrLUT,nphr_regions{1},delayStats,...
    'NPHR early-delay correct-rate OpGen specificity',...
    'Sample-selective fraction during delay') %nphr_model.m

plotOneCorr(nphr_arr(2,:),nphrLUT,nphr_regions{1},delayStats,...
    'NPHR early-delay hit OpGen specificity',...
    'Sample-selective fraction during delay') %nphr_model.m

plotOneCorr(nphr_arr(3,:),nphrLUT,nphr_regions{1},delayStats,...
    'NPHR early-delay rejection OpGen specificity',...
    'Sample-selective fraction during delay') %nphr_model.m

plotOneCorr(nphr_arr(4,:),nphrLUT,nphr_regions{2},delayStats,...
    'NPHR late-delay correct-rate OpGen specificity',...
    'Sample-selective fraction during delay') %nphr_model.m
plotOneCorr(nphr_arr(5,:),nphrLUT,nphr_regions{2},delayStats,...
    'NPHR late-delay hit OpGen specificity',...
    'Sample-selective fraction during delay') %nphr_model.m
plotOneCorr(nphr_arr(6,:),nphrLUT,nphr_regions{2},delayStats,...
    'NPHR late-delay rejection OpGen specificity',...
    'Sample-selective fraction during delay') %nphr_model.m


plotOneCorr(nphr_arr(7,:),nphrLUT,nphr_regions{3},testStats,...
    'NPHR decision-making correct-rate OpGen specificity',...
    'Pair-selective fraction during decision making') %nphr_model.m

plotOneCorr(nphr_arr(8,:),nphrLUT,nphr_regions{3},testStats,...
    'NPHR decision-making hit OpGen specificity',...
    'Pair-selective fraction during decision making') %nphr_model.m

plotOneCorr(nphr_arr(9,:),nphrLUT,nphr_regions{3},testStats,...
    'NPHR decision-making rejection OpGen specificity',...
    'Pair-selective fraction during decision making') %nphr_model.m




function plotOneCorr(opgen_vec,LUT,region_label,ephysStats,xlbl,ylbl)
fWT=figure('Color','w','Position',[50,50,1250,1250]);
hold on;
corrMatWT=[];
for i=1:length(opgen_vec)
    opRegion=region_label{i};
    opEffWT=opgen_vec(i);
    if any(strcmp(opRegion,LUT(:,1))) && ~isempty(LUT{strcmp(opRegion,LUT(:,1)),2})
        ephysRegion=LUT(strcmp(opRegion,LUT(:,1)),2);
    else
        ephysRegion=opRegion;
    end
    if any(strcmp(ephysRegion,ephysStats.rowheaders))
        selRatWT=ephysStats.data(strcmp(ephysRegion,ephysStats.rowheaders),2);
        yspan=max(ephysStats.data(:,2))-min(ephysStats.data(:,2));
        dy=yspan*0.005;
        figure(fWT);
        LWT=plot(opEffWT,selRatWT,'ro','MarkerFaceColor','r');
        text(opEffWT,selRatWT-dy,opRegion,'HorizontalAlignment','center','FontSize',12);
        corrMatWT=[corrMatWT;opEffWT,selRatWT];
    end
    
end


[rWT,pWT]=corrcoef(corrMatWT(:,1),corrMatWT(:,2));


fprintf('corr coef WT r = %0.3f, p = %0.3f\n',rWT(1,2),pWT(1,2));

figure(fWT);
% legend([LWT],{'Welltrained'});
xlabel(xlbl,'FontSize',14);
ylabel(ylbl,'FontSize',14);
text(min(xlim())+0.1*diff(xlim()),min(ylim)+0.1*diff(ylim()),sprintf('r = %0.3f, p = %0.3f',rWT(2),pWT(2)),'FontSize',16);
title('Ephys-optogene correlation, DPA, Welltrained')
figIdx=evalin('base','figIdx');
assignin('base','figIdx',figIdx+1);
print(sprintf('F%08d_r%0.3f_p%0.3f.png',figIdx,rWT(1,2),pWT(1,2)),'-dpng','-r300');
end
