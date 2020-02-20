function par_corr(effs, tasks)
p = gcp();
% To request multiple evaluations, use a loop.
futures=parallel.FevalFuture.empty(0);
for use_eff=effs
    for task_idx=tasks
        futures(end+1) = parfeval(p,@one_corr,0,use_eff,task_idx); % Square size determined by idx
    end
end
% Collect the results as they become available.
for idx = 1:18
    % fetchNext blocks until next results are available.
    completedIdx = fetchNext(futures);
    fprintf('Got result with index: %d.\n', completedIdx);
end

end

function one_corr(use_eff,task_idx)

path='coding_features.csv';
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
%% up to 4 comp interaction
currAIC=realmax;
int_result=cell(1,6);

for n=2:4
    C=nchoosek(2:23,n);
    for i=1:size(C,1)
        mdl=fitglm(glm_mat(:,C(i,:)),glm_mat(:,1),'interactions');
        if mdl.ModelCriterion.AIC<currAIC
            currAIC=mdl.ModelCriterion.AIC;
            int_result{1}=mdl;
            int_result{2}=mdl.ModelCriterion.AIC;
            int_result{3}=mdl.Rsquared.Ordinary;
            int_result{4}=coefTest(mdl);
            int_result{5}='interactions';
            int_result{6}=C(i,:);
        end
        disp(size(int_result,1));
    end
end

%% cross validation

cv_results=nan(size(glm_mat,1),2);
for i=1:size(glm_mat,1)
    cv_mat=glm_mat;
    cv_mat(i,:)=[];
    mdl=fitglm(cv_mat(:,int_result{6}),cv_mat(:,1),int_result{5});
    pred=mdl.predict(glm_mat(i,int_result{6}));
    cv_results(i,:)=[glm_mat(i,1),pred];
end

[r,p]=corrcoef(cv_results(:,1),cv_results(:,2));
r=r(1,2);
p=p(1,2);
suffix=value_labels{task_idx};
if use_eff
    suffix=['EFF_',suffix];
end
save(sprintf('GLM_INT_SEL_PCA_HYB_vgat_%s.mat',suffix),'int_result','cv_results','r','p')


%% 22 comp linear combination

currAIC=realmax;
int_result=cell(1,6);

for n=2:22
    C=nchoosek(2:23,n);
    for i=1:size(C,1)
        mdl=fitglm(glm_mat(:,C(i,:)),glm_mat(:,1),'linear');
        if mdl.ModelCriterion.AIC<currAIC
            currAIC=mdl.ModelCriterion.AIC;
            int_result{1}=mdl;
            int_result{2}=mdl.ModelCriterion.AIC;
            int_result{3}=mdl.Rsquared.Ordinary;
            int_result{4}=coefTest(mdl);
            int_result{5}='linear';
            int_result{6}=C(i,:);
        end
        disp(size(int_result,1));
    end
end

%% cross validation

cv_results=nan(size(glm_mat,1),2);
for i=1:size(glm_mat,1)
    cv_mat=glm_mat;
    cv_mat(i,:)=[];
    mdl=fitglm(cv_mat(:,int_result{6}),cv_mat(:,1),int_result{5});
    pred=mdl.predict(glm_mat(i,int_result{6}));
    cv_results(i,:)=[glm_mat(i,1),pred];
end

[r,p]=corrcoef(cv_results(:,1),cv_results(:,2));
r=r(1,2);
p=p(1,2);
suffix=value_labels{task_idx};
if use_eff
    suffix=['EFF_',suffix];
end
save(sprintf('GLM_LNR_SEL_PCA_HYB_vgat_%s.mat',suffix),'int_result','cv_results','r','p')

end


