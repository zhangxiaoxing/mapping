use_eff=false;
path='K:\code\tca\Batch30B\tca_su_factors.csv';
features=importdata(path,',',0);
load('opgene_ephys_corr.mat','translate')
if use_eff
    load('nphr_mat.mat','value_arr','value_labels','region_list')
else
    load('zx_specificity_mat.mat','nphr_arr','nphr_regions','value_labels');
    region_list=nphr_regions{1};
end
for task_idx=1:size(value_arr,1)
    glm_mat=[];
    regions=cell(0,2);
    for i=1:size(translate,1)
        if ~isempty(translate{i,2})
            opgen_reg=translate{i,1};
            ephys_reg=translate{i,2};
            opgen_idx=strcmp(opgen_reg,region_list);
            ephys_idx=strcmp(ephys_reg,features.rowheaders);
            
            if any(opgen_idx) && any(ephys_idx)
                opgen=value_arr(task_idx,opgen_idx);
                ephys=features.data(ephys_idx,:);
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
        cv_results(i,:)=[glm_mat(i,1),pred]; % target prediction
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
    save(sprintf('GLM_selec_nphr_%s.mat',suffix),'int_result','cv_results','r','p')
end