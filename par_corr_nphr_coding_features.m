function par_corr_nphr_coding_features(effs, tasks,subset)
if ~exist('subset','var')
    subset=false;
end

p = gcp();
futures=parallel.FevalFuture.empty(0);
for use_eff=effs
    for task_idx=tasks
        futures(end+1) = parfeval(p,@one_corr,0,use_eff,task_idx, subset); % Square size determined by idx
    end
end
% Collect the results as they become available.
for idx = 1:numel(futures)
    % fetchNext blocks until next results are available.
    completedIdx = fetchNext(futures);
    fprintf('Got result with index: %d.\n', completedIdx);
end

end


function one_corr(use_eff,task_idx,subset)
path='glm_coding_features.csv';
features=importdata(path,',',0);
load('opgene_ephys_corr.mat','translate')
if use_eff
    load('nphr_mat.mat','value_arr','value_labels','region_list')
else
    load('zx_specificity_mat.mat','nphr_arr','nphr_regions','value_labels');
    value_arr=nphr_arr;
    region_list=nphr_regions{1};
end

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

if subset
    sub_sel=ismember(regions(:,1),{'AId','AON','M2','SS','RSPd','CA1v','DG'});
    %TODO finish subset list
end


curr_aic=realmax;
last_mdl_idx=0;
int_result=cell(0,6);



for n=1:(size(glm_mat,2)-1)
    C=nchoosek(2:size(glm_mat,2),n);
    for i=1:size(C,1)
        mdl=fitglm(glm_mat(:,C(i,:)),glm_mat(:,1),'linear');
        new_idx=size(int_result,1)+1;
        int_result{new_idx,2}=mdl.ModelCriterion.AIC;
        int_result{new_idx,3}=mdl.Rsquared.Ordinary;
        int_result{new_idx,4}=coefTest(mdl);
        int_result{new_idx,5}='linear';
        int_result{new_idx,6}=C(i,:);
        if(mdl.ModelCriterion.AIC)<curr_aic
            if last_mdl_idx>0
                int_result{last_mdl_idx,1}=NaN;
            end
            int_result{new_idx,1}=mdl;
            last_mdl_idx=new_idx;
            curr_aic=mdl.ModelCriterion.AIC;
        end
    end
end

for n=2:4
    C=nchoosek(2:size(glm_mat,2),n);
    for i=1:size(C,1)
        mdl=fitglm(glm_mat(:,C(i,:)),glm_mat(:,1),'interactions');
        new_idx=size(int_result,1)+1;
        int_result{new_idx,2}=mdl.ModelCriterion.AIC;
        int_result{new_idx,3}=mdl.Rsquared.Ordinary;
        int_result{new_idx,4}=coefTest(mdl);
        int_result{new_idx,5}='interactions';
        int_result{new_idx,6}=C(i,:);
        
        if(mdl.ModelCriterion.AIC)<curr_aic
            if last_mdl_idx>0
                int_result{last_mdl_idx,1}=NaN;
            end
            int_result{new_idx,1}=mdl;
            last_mdl_idx=new_idx;
            curr_aic=mdl.ModelCriterion.AIC;
        end
    end
end

%% cross validation

cv_results=nan(size(glm_mat,1),2);
for i=1:size(glm_mat,1)
    cv_mat=glm_mat;
    cv_mat(i,:)=[];
    mdl=fitglm(cv_mat(:,int_result{last_mdl_idx,6}),cv_mat(:,1),int_result{last_mdl_idx,5});
    pred=mdl.predict(glm_mat(i,int_result{last_mdl_idx,6}));
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
save(sprintf('GLM_select_feat_nphr_%s.mat',suffix),'int_result','cv_results','r','p')

end
