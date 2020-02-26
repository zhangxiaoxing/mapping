function par_ridge_nphr_selectivity(effs, tasks,lambda,use_mean,subset)
addpath('FastRidge');
if ~exist('subset','var')
    subset=false;
end
if ~exist('lambda','var')
   lambda=0.1; 
end
if ~exist('use_mean','var')
    use_mean=true;
end

p = gcp();
futures=parallel.FevalFuture.empty(0);
for use_eff=effs
    for task_idx=tasks
        futures(end+1) = parfeval(p,@one_corr,0,use_eff,task_idx, lambda, use_mean, subset); % Square size determined by idx
    end
end
% Collect the results as they become available.
for idx = 1:numel(futures)
    % fetchNext blocks until next results are available.
    completedIdx = fetchNext(futures);
    fprintf('Got result with index: %d.\n', completedIdx);
end

end


function one_corr(use_eff,task_idx,lambda,use_mean,subset)
path='glm_coding_features.csv';
features=importdata(path,',',0);
if use_eff
    load('nphr_mat.mat','mean_arr','median_arr','value_labels','region_list','nphr_lut')
    if use_mean
        value_arr=mean_arr;
    else
        value_arr=median_arr;
    end
else
    load('zx_specificity_mat.mat','nphr_arr','nphr_regions','value_labels');
    load('nphr_mat.mat','nphr_lut')
    value_arr=nphr_arr;
    region_list=nphr_regions{1};
end

glm_mat=[];
regions=cell(0,2);
for i=1:size(nphr_lut,1)
    if ~isempty(nphr_lut{i,2})
        opgen_reg=nphr_lut{i,1};
        ephys_reg=nphr_lut{i,2};
        opgen_idx=strcmp(opgen_reg,region_list);
        ephys_idx=strcmp(ephys_reg,features.rowheaders);
        
        if any(opgen_idx) && any(ephys_idx) && features.data(ephys_idx,1)>=50
            opgen=value_arr(task_idx,opgen_idx);
            ephys=features.data(ephys_idx,2:end);
            glm_mat(end+1,:)=[opgen,ephys];
            regions(end+1,:)={opgen_reg,ephys_reg};
        else
            continue
        end
    end
end

% if subset
%     sub_sel=ismember(regions(:,1),{'AId','AON','M2','SS','RSPd','CA1v','DG'});
%     %TODO finish subset list
% end

% normalize X, center Y
% it's actually implemented in the downloaded function
% glm_mat(:,2:end)=normalize(glm_mat(:,2:end));
% glm_mat(:,1)=glm_mat(:,1)-mean(glm_mat(:,1));

int_result=cell(0,7);



for n=1:(size(glm_mat,2)-1)
    C=nchoosek(2:size(glm_mat,2),n);
    for i=1:size(C,1)
        [beta, b0, tau2, DOF, lambda, score] = fastridge(glm_mat(:,C(i,:)), glm_mat(:,1), 'criterion', 'aicc', 'lambda',lambda);
%         mdl=fitglm(glm_mat(:,C(i,:)),glm_mat(:,1),'linear');
        [rr,pp]=corrcoef(glm_mat(:,1),glm_mat(:,C(i,:))*beta);
        new_idx=size(int_result,1)+1;
        int_result{new_idx,1}=beta;
        int_result{new_idx,2}=score;
        int_result{new_idx,3}=rr(1,2)^2;
        int_result{new_idx,4}=pp(1,2);
        int_result{new_idx,5}='linear';
        int_result{new_idx,6}=C(i,:);
        int_result{new_idx,7}=[glm_mat(:,1),glm_mat(:,C(i,:))*beta];
    end
end

for n=2:4
    C=nchoosek(2:size(glm_mat,2),n);
    for i=1:size(C,1)
        X_from=glm_mat(:,C(i,:));
        X=x2fx(X_from,'interaction');
        X(:,1)=[];
        [beta, b0, tau2, DOF, lambda, score] = fastridge(X, glm_mat(:,1), 'criterion', 'aicc', 'lambda',lambda);
%         mdl=fitglm(glm_mat(:,C(i,:)),glm_mat(:,1),'linear');
        [rr,pp]=corrcoef(glm_mat(:,1),X*beta);
        new_idx=size(int_result,1)+1;
        int_result{new_idx,1}=beta;
        int_result{new_idx,2}=score;
        int_result{new_idx,3}=rr(1,2)^2;
        int_result{new_idx,4}=pp(1,2);
        int_result{new_idx,5}='interact';
        int_result{new_idx,6}=C(i,:);
        int_result{new_idx,7}=[glm_mat(:,1),X*beta];
    end
end

%% cross validation
[~,low_aic_idx]=min(cell2mat(int_result(:,2)));
cv_results=nan(size(glm_mat,1),2);
for i=1:size(glm_mat,1)
    cv_mat=glm_mat;
    cv_mat(i,:)=[];
    if strcmp(int_result{low_aic_idx,5},'linear')
        X=cv_mat(:,int_result{low_aic_idx,6});
        [beta, b0, tau2, DOF, lambda, score] = fastridge(X, cv_mat(:,1), 'criterion', 'aicc', 'lambda',lambda);
        X=glm_mat(i,int_result{low_aic_idx,6});
        pred=X*beta+b0;
        cv_results(i,:)=[glm_mat(i,1),pred]; % target prediction
    else
        X_from=cv_mat(:,int_result{low_aic_idx,6});
        X=x2fx(X_from,'interaction');
        X(:,1)=[];
        [beta, b0, tau2, DOF, lambda, score] = fastridge(X, cv_mat(:,1), 'criterion', 'aicc', 'lambda',lambda);
        X_from=glm_mat(i,int_result{low_aic_idx,6});
        X=x2fx(X_from,'interaction');
        pred=X(2:end)*beta+b0;
        cv_results(i,:)=[glm_mat(i,1),pred]; % target prediction
    end
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
if use_mean
    suffix=['Mean_',suffix];
else
    suffix=['Median_',suffix];
end
save(sprintf('Ridge_select_feat_nphr_lambda_%0.3f_%s.mat',lambda,suffix),'int_result','cv_results','r','p','lambda','regions')

end
