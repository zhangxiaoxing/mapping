close all
fAll=figure('Color','w','Position',[1500,-50,1080,900]);
hold on;
fLN=figure('Color','w','Position',[1500,-50,1080,900]);
hold on;
fWT=figure('Color','w','Position',[1500,-50,1080,900]);
hold on;
fDelta=figure('Color','w','Position',[1500,-50,1080,900]);
hold on;
corrMat=[];
corrMatLN=[];
corrMatWT=[];
corrMatDelta=[];
for i=2:size(nphrHeatMat,2)
    opRegion=nphrHeatMat{1,i};
%     opEffLN=mean(cell2mat(nphrHeatMat(2:4,i)));
%     opEffWT=mean(cell2mat(nphrHeatMat(14:15,i)));
    opEffLN=nanmean(cell2mat(nphrHeatMat(2:4,i)));
    opEffWT=cell2mat(nphrHeatMat(15,i));
    if any(strcmp(opRegion,translate(:,1))) && ~isempty(translate{strcmp(opRegion,translate(:,1)),2})
        ephysRegion=translate(strcmp(opRegion,translate(:,1)),2);
        if any(strcmp(ephysRegion,delayStats.rowheaders))
            selRatLN=delayStats.data(strcmp(ephysRegion,delayStats.rowheaders),1);
            selRatWT=delayStats.data(strcmp(ephysRegion,delayStats.rowheaders),2);
            figure(fAll);
            LALN=plot(opEffLN,selRatLN,'bo','MarkerFaceColor','b');
            text(opEffLN,selRatLN-0.006,opRegion,'HorizontalAlignment','center');
            LAWT=plot(opEffWT,selRatWT,'ro','MarkerFaceColor','r');
            plot([opEffLN;opEffWT],[selRatLN;selRatWT],':k');
            figure(fDelta);
            LDelta=plot(opEffWT-opEffLN,selRatWT-selRatLN,'ko','MarkerFaceColor','k');
            text(opEffWT-opEffLN-0.005,selRatWT-selRatLN-0.005,opRegion,'HorizontalAlignment','center');
           
            figure(fLN);
            LLN=plot(opEffLN,selRatLN,'bo','MarkerFaceColor','b');
            text(opEffLN,selRatLN-0.006,opRegion,'HorizontalAlignment','center');
            
            figure(fWT);
            LWT=plot(opEffWT,selRatWT,'ro','MarkerFaceColor','r');
            text(opEffWT,selRatWT-0.006,opRegion,'HorizontalAlignment','center');
            
            
            corrMat=[corrMat;opEffLN,selRatLN;opEffWT,selRatWT];
            corrMatLN=[corrMatLN;opEffLN,selRatLN];
            corrMatWT=[corrMatWT;opEffWT,selRatWT];
            corrMatDelta=[corrMatDelta;opEffWT-opEffLN,selRatWT-selRatLN];
        end
    end
end

[rall,pall]=corrcoef(corrMat(:,1),corrMat(:,2));
[rLN,pLN]=corrcoef(corrMatLN(:,1),corrMatLN(:,2));
[rWT,pWT]=corrcoef(corrMatWT(:,1),corrMatWT(:,2));
[rDelta,pDelta]=corrcoef(corrMatDelta(:,1),corrMatDelta(:,2));

fprintf('corr coef all r = %0.3f, p = %0.3f',rall,pall);
fprintf('corr coef LN r = %0.3f, p = %0.3f',rLN,pLN);
fprintf('corr coef WT r = %0.3f, p = %0.3f',rWT,pWT);
fprintf('corr coef Delta r = %0.3f, p = %0.3f',rDelta,pDelta);


figure(fAll);
legend([LALN,LAWT],{'Learning','Welltrained'});
xlabel('DPA delay effect size');
ylabel('Sample selective fraction during delay');
text(min(xlim())+0.1*diff(xlim()),min(ylim)+0.1*diff(ylim()),sprintf('r = %0.3f, p = %0.3f',rall(2),pall(2)));
% title('Ephys-optogene correlation, Dualtask-simple trials, LN+WT')a
title('Ephys-optogene correlation, DPA, LN+WT')

figure(fLN);
legend([LLN],{'Learning'});
xlabel('DPA delay effect size');
ylabel('Sample-selective fraction during delay');
text(min(xlim())+0.1*diff(xlim()),min(ylim)+0.1*diff(ylim()),sprintf('r = %0.3f, p = %0.3f',rLN(2),pLN(2)));
title('Ephys-optogene correlation, DPA, Learning')

figure(fWT);
legend([LWT],{'Welltrained'});
xlabel('DPA delay effect size');
ylabel('Sample-selective fraction during delay');
text(min(xlim())+0.1*diff(xlim()),min(ylim)+0.1*diff(ylim()),sprintf('r = %0.3f, p = %0.3f',rWT(2),pWT(2)));
title('Ephys-optogene correlation, DPA, Welltrained')

figure(fDelta);
legend([LDelta],{'Change in Select frac., Effect size'});
xlabel('Change in DPA delay effect size');
ylabel('Change in Sample-selective fraction during delay');
text(min(xlim())+0.1*diff(xlim()),min(ylim)+0.1*diff(ylim()),sprintf('r = %0.3f, p = %0.3f',rDelta(2),pDelta(2)));
title('Ephys-optogene correlation, Dualtask-simple trials, Learning-welltrained delta')


return

