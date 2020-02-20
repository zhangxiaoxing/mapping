use_eff=false;
path='PCA_Comp.csv';
load('vgatLUT.mat')
if use_eff
    load('vgat_heat_mat.mat')
    vgat_arr=cell2mat(vgat_heat_mat(13:21,2:end));
    vgat_regions=vgat_heat_mat(1,2:end);
    value_labels=vgat_heat_mat(13:21,1);
else
    load('zx_specificity_mat.mat','vgat_arr','vgat_regions','value_labels');
end
tca_coeff=importdata(path,',',0);
for task_idx=1:size(vgat_arr,1)
    glm_mat=[];
    regions=cell(0,2);
    for i=1:size(vgatLUT,1)
        if ~isempty(vgatLUT{i,2})
            opgen_reg=vgatLUT{i,1};
            ephys_reg=vgatLUT{i,2};
            if use_eff
                opgen_idx=strcmp(ephys_reg,vgat_regions);
            else
                opgen_idx=strcmp(opgen_reg,vgat_regions{1});
            end
            ephys_idx=strcmp(ephys_reg,tca_coeff.rowheaders);
            
            if any(opgen_idx) && any(ephys_idx)
                opgen=vgat_arr(task_idx,opgen_idx);
                ephys=tca_coeff.data(ephys_idx,:);
                glm_mat(end+1,:)=[opgen,ephys];
                regions(end+1,:)={opgen_reg,ephys_reg};
            else
                continue
            end
        end
    end
    
    int_result=cell(0,6);
    
    for i=2:12
        mdl=fitglm(glm_mat(:,i),glm_mat(:,1),'linear');
        int_result{end+1,1}=mdl;
        int_result{end,2}=mdl.ModelCriterion.AIC;
        int_result{end,3}=mdl.Rsquared.Ordinary;
        int_result{end,4}=coefTest(mdl);
        int_result{end,5}='linear';
        int_result{end,6}=i;
    end
    
    
    for n=2:11
        C=nchoosek(2:12,n);
        for i=1:size(C,1)
            mdl=fitglm(glm_mat(:,C(i,:)),glm_mat(:,1),'linear');
            int_result{end+1,1}=mdl;
            int_result{end,2}=mdl.ModelCriterion.AIC;
            int_result{end,3}=mdl.Rsquared.Ordinary;
            int_result{end,4}=coefTest(mdl);
            int_result{end,5}='linear';
            int_result{end,6}=C(i,:);
        end
    end
    
    for n=2:4
        C=nchoosek(2:12,n);
        for i=1:size(C,1)
            mdl=fitglm(glm_mat(:,C(i,:)),glm_mat(:,1),'interactions');
            int_result{end+1,1}=mdl;
            int_result{end,2}=mdl.ModelCriterion.AIC;
            int_result{end,3}=mdl.Rsquared.Ordinary;
            int_result{end,4}=coefTest(mdl);
            int_result{end,5}='interactions';
            int_result{end,6}=C(i,:);
        end
    end
    
    
    [~,Imin_aic]=min([int_result{:,2}]);
    
    %% cross validation
    
    
    cv_results=nan(size(glm_mat,1),2);
    for i=1:size(glm_mat,1)
        cv_mat=glm_mat;
        cv_mat(i,:)=[];
        mdl=fitglm(cv_mat(:,int_result{Imin_aic,6}),cv_mat(:,1),int_result{Imin_aic,5});
        pred=mdl.predict(glm_mat(i,int_result{Imin_aic,6}));
        cv_results(i,:)=[glm_mat(i,1),pred];
    end
    % figure()
    % plot(cv_results(:,1),cv_results(:,2),'k.')
    % ylim([-1,1])
    % xlim([-1,1])
    [r,p]=corrcoef(cv_results(:,1),cv_results(:,2));
    r=r(1,2);
    p=p(1,2);
    suffix=value_labels{task_idx};
    if use_eff
        suffix=['EFF_',suffix];
    end
    save(sprintf('GLM_PCA_vgat_%s.mat',suffix),'int_result','cv_results','r','p')
end
