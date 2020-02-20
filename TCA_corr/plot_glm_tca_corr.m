fl=ls('GLM_TCA*.mat');
for i=1:size(fl,1)
    load(strtrim(fl(i,:)));
    fh=figure('Color','w','Position',[50,50,960,960]);
    subplot(2,2,1);
    [~,Iaic]=min(cell2mat(int_result(:,2)));
    mdl=int_result{Iaic,1};
    plot(mdl.Fitted{:,1},mdl.Variables{:,end},'r.','MarkerSize',10);
    text(min(xlim())+0.15*diff(xlim()),min(ylim())+0.85*diff(ylim()),sprintf('rsq = %.3f, p = %.3f',int_result{Iaic,3},int_result{Iaic,4}),'HorizontalAlignment','left')
    xlabel('GLM regression');
    ylabel('Experimental data')
    subplot(2,2,2);
    plot(cv_results(:,2),cv_results(:,1),'b.','MarkerSize',10);
    text(min(xlim())+0.15*diff(xlim()),min(ylim())+0.85*diff(ylim()),sprintf('rsq = %.3f, p = %.3f',r*r,p),'HorizontalAlignment','left')
    xlabel('Leave-one-region-out prediction');
    ylabel('Experimental data')
    
    
    uit=uitable(fh);
    uit.Data=table2array(mdl.Coefficients);
    uit.Position=[32,32,896,384];
    uit.ColumnName=mdl.Coefficients.Properties.VariableNames;
    if mdl.NumPredictors==4 && strcmp(int_result{Iaic,5},'interactions')
        RowNames={'Intercept',...
            sprintf('Comp_%02d', int_result{Iaic,6}(1)),...
            sprintf('Comp_%02d', int_result{Iaic,6}(2)),...
            sprintf('Comp_%02d', int_result{Iaic,6}(3)),...
            sprintf('Comp_%02d', int_result{Iaic,6}(4)),...
            sprintf('C%02d:C%02d interact', int_result{Iaic,6}(1),int_result{Iaic,6}(2)),...
            sprintf('C%02d:C%02d interact', int_result{Iaic,6}(1),int_result{Iaic,6}(3)),...
            sprintf('C%02d:C%02d interact', int_result{Iaic,6}(1),int_result{Iaic,6}(4)),...
            sprintf('C%02d:C%02d interact', int_result{Iaic,6}(2),int_result{Iaic,6}(3)),...
            sprintf('C%02d:C%02d interact', int_result{Iaic,6}(2),int_result{Iaic,6}(4)),...
            sprintf('C%02d:C%02d interact', int_result{Iaic,6}(3),int_result{Iaic,6}(4))};
    end
    
    uit.RowName=RowNames;
    uit.FontSize=12;
    uit.ColumnWidth={160,160,160,160};
    tiStr=num2str(i);
    if contains(fl(i,:),'vgat')
        tiStr=strcat(tiStr, 'VGAT');
    else
        tiStr=strcat(tiStr,'NPHR');
    end
    if contains(fl(i,:),'EFF')
        tiStr=strcat(tiStr, ' Effect Size');
    else
        tiStr=strcat(tiStr,' Opgen Selectivity');
    end
    tiStr=strcat(tiStr,'-ephys correlation-',regexp(fl(i,:),'(DM|ED|LD)-\w{2,8}','match'));
    sgtitle(tiStr)
    print([char(replace(tiStr,' ','_')),'.png'],'-dpng','-r150');
    close(fh)
end
