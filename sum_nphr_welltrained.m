%% laser-off then laser on
cd('D:\Mapping\mapping');
%% list regions
rootpath='D:\Mapping\mapping NpHR';
folders=dir(rootpath);
folders(startsWith({folders.name},'.') | ~[folders.isdir] |startsWith({folders.name},'effect'))=[];
%% into subgroups of regions
all_branches=cell(0,3);
for i=1:length(folders)
    cd(fullfile(rootpath,folders(i).name));
    subfolders=dir();
    if any(strcmpi({subfolders.name},'anadata'))
        cd('anadata')
    end
    subfolders=dir();
    subfolders(startsWith({subfolders.name},'.') | ~[subfolders.isdir])=[];
    if any(strcmp({subfolders.name},'DPA-DM')) && ...
            any(strcmp({subfolders.name},'DPA-ED')) && ...
            any(strcmp({subfolders.name},'DPA-LD'))
    
    %% process 3 well-trained tasks
        one_branch_DM=processOneBranch(pwd(),'DPA-DM');
        one_branch_ED=processOneBranch(pwd(),'DPA-ED');
        one_branch_LD=processOneBranch(pwd(),'DPA-LD');
        all_branches=vertcat(all_branches,one_branch_ED,one_branch_LD,one_branch_DM);
    end    
end
all_branches=combineDuplicate(all_branches);

cd('D:\Mapping\mapping');
save('nphr_perf.mat','all_branches');



%% regional pref
function out=processOneBranch(rootpath,task)
%% one region group, one task

branch=fullfile(rootpath,task);
%files=dir('*.mat'); % not going to use ready made per trial data due to
%performance window
folders=dir(branch);
folders(startsWith({folders.name},'.') | ~[folders.isdir] |startsWith({folders.name},'effect'))=[];
all_region=cell(0,3);
for i=1:length(folders)
    curr_region=mergeRegions(folders(i).name);
    one_region=processOneRegion(fullfile(folders(i).folder,folders(i).name));
    all_region(end+1,:)={task,curr_region,one_region};
end
out=all_region;
end
  

function out=processOneRegion(path)
files=dir(fullfile(path,'*.mat'));
one_region=[];

for f=files'
    fstr=load(fullfile(f.folder,f.name),'DPATrial');
    cleared_trial=clearBadPerf(fstr.DPATrial);
    if size(cleared_trial,1)<40
        continue
    end
    correct_on=nnz(ismember(cleared_trial(:,6),[5,7]) & cleared_trial(:,4)==1)/...
        nnz(cleared_trial(:,4)==1);
    correct_off=nnz(ismember(cleared_trial(:,6),[5,7]) & cleared_trial(:,4)==0)/...
        nnz(cleared_trial(:,4)==0);
    
    hit_on=nnz(cleared_trial(:,6)==7 & cleared_trial(:,4)==1)/...
        nnz(ismember(cleared_trial(:,6),[6,7]) & cleared_trial(:,4)==1);
    hit_off=nnz(cleared_trial(:,6)==7 & cleared_trial(:,4)==0)/...
        nnz(ismember(cleared_trial(:,6),[6,7]) & cleared_trial(:,4)==0);
    
    
    reject_on=nnz(cleared_trial(:,6)==5 & cleared_trial(:,4)==1)/...
        nnz(ismember(cleared_trial(:,6),[4,5]) & cleared_trial(:,4)==1);
    reject_off=nnz(cleared_trial(:,6)==5 & cleared_trial(:,4)==0)/...
        nnz(ismember(cleared_trial(:,6),[4,5]) & cleared_trial(:,4)==0);
    one_region=[one_region;correct_off,correct_on,hit_off,hit_on,reject_off,reject_on];
end
out=one_region;
end


function out=clearBadPerf(facSeq)

    if length(facSeq)>=80
        facSeq(:,8)=0;
        i=80;
        while i<=length(facSeq)
            goodOff=nnz((facSeq(i-79:i,6)==5 | facSeq(i-79:i,6)==7)& facSeq(i-79:i,4)==0);
            if goodOff>=nnz(facSeq(i-79:i,3)==0)*75/100
                facSeq(i-79:i,8)=1;
            end
            i=i+1;
        end
        out=facSeq(facSeq(:,8)==1,:);
    else
        out=[];
    end
end

function out=mergeRegions(in)
if ismember(in,{'CA1d','CA1v','CA3d','CA3v','DGd','DGv','SUBd','SUBv','M2p'})
    out=in(1:end-1);
elseif ismember(in,{'APC','PPC'})
    out='PIR';
elseif ismember(in,{'AONup','AONv'})
    out='AON';
else
    out=in;
end
end


function out=combineDuplicate(in)
tasks=unique(in(:,1));
regs=unique(in(:,2));
out=cell(0,3);
for t=tasks'
    for r=regs'
        out(end+1,:)={char(t),char(r),cell2mat(in(strcmp(in(:,1),t) & strcmp(in(:,2),r),3))};
    end
end
end


