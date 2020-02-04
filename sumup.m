%% laser-off then laser on
cd('K:\Mapping\codes');
%% list regions
rootpath='K:\Mapping\mapping VGAT';
folders=dir(rootpath);

folders(startsWith({folders.name},'.') | ~[folders.isdir])=[];
regions=cell(0);
keep=[];
for i=1:length(folders)
    cd(fullfile(rootpath,folders(i).name));
    subfolders=dir();
    subfolders(startsWith({subfolders.name},'.') | ~[subfolders.isdir])=[];
    for j=1:length(subfolders)
        regions{end+1}=subfolders(j).name;
    end
    cd('..');
    
end
regions=unique(regions);
%% regional pref
wtRows=cell(0);
lnRows=cell(0);

wtRowsHitMiss=cell(0);

wtRowsCRFalse=cell(0);

for ridx=1:length(regions)
    for fidx=1:length(folders)
        if ~contains(folders(fidx).name,'lear','IgnoreCase',true)
            %% well-trained
            cd(fullfile(rootpath,folders(fidx).name));
            if exist(regions{ridx},'dir')
                cd(regions{ridx});
                files=dir('*.mat');
                for fiidx=1:length(files)
                    miceId=regexp(files(fiidx).name,'^\d+(?=-)','match');
                    prog=regexp(files(fiidx).name,'(?<=^\d+-)\d+(?=-)','match');
                    fs=load(files(fiidx).name);
                    fsn=fieldnames(fs);
                    ftrialIdx=find(endsWith(fsn,'Trial'));
                    if numel(ftrialIdx)>=1
                        wtRows(end+1,1:4)={folders(fidx).name,regions{ridx},miceId{1},prog{1}};
                        wtRowsHitMiss(end+1,1:4)={folders(fidx).name,regions{ridx},miceId{1},prog{1}};
                        wtRowsCRFalse(end+1,1:4)={folders(fidx).name,regions{ridx},miceId{1},prog{1}};
                        dual_diff=0;
                        for taskIdx=1:numel(ftrialIdx)
                            perfT=fs.(fsn{ftrialIdx(taskIdx)});

                            
                            expType=strrep(fsn{ftrialIdx(taskIdx)},'Trial','');
                            if strcmp(expType,'DPA') && (strcmp(prog{1},'4351') || strcmp(prog{1},'4353') || strcmp(prog{1},'4341'))
                                perfT=clearBadPerf(perfT);
                                if isempty(perfT)
                                    continue;
                                end
                                
                                wtRows{end,4+taskIdx*3-2}='Simple_DPA';
                                wtRowsHitMiss{end,4+taskIdx*3-2}='Simple_DPA';
                                wtRowsCRFalse{end,4+taskIdx*3-2}='Simple_DPA';
                                perfOff=sum(ismember(perfT(perfT(:,4)==0,6),[5,7]))./sum(perfT(:,4)==0);
                                perfOn=sum(ismember(perfT(perfT(:,4)==1,6),[5,7]))./sum(perfT(:,4)==1);
                                wtRows{end,4+taskIdx*3-1}=perfOff;
                                wtRows{end,4+taskIdx*3}=perfOn;
                                HitOff=sum(perfT(perfT(:,4)==0,6)==7)./sum(perfT(:,4)==0 & ismember(perfT(:,6),[6,7]));
                                HitOn=sum(perfT(perfT(:,4)==1,6)==7)./sum(perfT(:,4)==1 & ismember(perfT(:,6),[6,7]));
                                wtRowsHitMiss{end,4+taskIdx*3-1}=HitOff;
                                wtRowsHitMiss{end,4+taskIdx*3}=HitOn;
                                CROff=sum(perfT(perfT(:,4)==0,6)==5)./sum(perfT(:,4)==0 & ismember(perfT(:,6),[4,5]));
                                CROn=sum(perfT(perfT(:,4)==1,6)==5)./sum(perfT(:,4)==1 & ismember(perfT(:,6),[4,5]));
                                wtRowsCRFalse{end,4+taskIdx*3-1}=CROff;
                                wtRowsCRFalse{end,4+taskIdx*3}=CROn;
                            elseif strcmp(expType,'DPA') && (strcmp(prog{1},'4364') || strcmp(prog{1},'4363'))
                                wtRows{end,4+taskIdx*3-2}='Dual_DPA_S';
                                perfOffS=sum(ismember(perfT(perfT(:,4)==0 & perfT(:,8)==0,6),[5,7]))./sum(perfT(:,4)==0 & perfT(:,8)==0);
                                perfOnS=sum(ismember(perfT(perfT(:,4)==1 & perfT(:,8)==0,6),[5,7]))./sum(perfT(:,4)==1 & perfT(:,8)==0);
                                wtRows{end,4+taskIdx*3-1}=perfOffS;
                                wtRows{end,4+taskIdx*3}=perfOnS;
                                dual_diff=1;
                                taskIdx=taskIdx+dual_diff;
                                wtRows{end,4+taskIdx*3-2}='Dual_DPA_D';
                                perfOffD=sum(ismember(perfT(perfT(:,4)==0 & perfT(:,8)==1,6),[5,7]))./sum(perfT(:,4)==0 & perfT(:,8)==1);
                                perfOnD=sum(ismember(perfT(perfT(:,4)==1 & perfT(:,8)==1,6),[5,7]))./sum(perfT(:,4)==1 & perfT(:,8)==1);
                                wtRows{end,4+taskIdx*3-1}=perfOffD;
                                wtRows{end,4+taskIdx*3}=perfOnD;
                            elseif contains(expType,'DRT') && strcmp(prog{1},'4341')
                                wtRows{end,4+taskIdx*3-2}='Simple_DRT';
                                perfOff=sum(ismember(perfT(perfT(:,4)==0,6),[5,7]))./sum(perfT(:,4)==0);
                                perfOn=sum(ismember(perfT(perfT(:,4)==1,6),[5,7]))./sum(perfT(:,4)==1);
                                wtRows{end,4+taskIdx*3-1}=perfOff;
                                wtRows{end,4+taskIdx*3}=perfOn;
                            elseif contains(expType,'DRT') && (strcmp(prog{1},'4364') || strcmp(prog{1},'4363'))
                                taskIdx=taskIdx+dual_diff;
                                wtRows{end,4+taskIdx*3-2}='Dual_DRT';
                                perfOff=sum(ismember(perfT(perfT(:,7)==0,5),[5,7]))./sum(perfT(:,7)==0);
                                perfOn=sum(ismember(perfT(perfT(:,7)==1,5),[5,7]))./sum(perfT(:,7)==1);
                                wtRows{end,4+taskIdx*3-1}=perfOff;
                                wtRows{end,4+taskIdx*3}=perfOn;
                            else
                                disp('Unexpected design');
                                perf=-1;
                            end
                        end
                    end
                end
                cd('..');
            end
            cd('..');
        elseif contains(folders(fidx).name,'lear','IgnoreCase',true)
            %% learning
            cd(folders(fidx).name);
            if exist(regions{ridx},'dir')
                cd(regions{ridx});
                dayFolders=dir();
                dayFolders(startsWith({dayFolders.name},'.') | ~[dayFolders.isdir])=[];
                for dayIdx=1:length(dayFolders)
                    cd(dayFolders(dayIdx).name);
                    files=dir('*.mat');
                    for fiidx=1:length(files)
                        miceId=regexp(files(fiidx).name,'^\d+(?=-)','match');
                        prog=regexp(files(fiidx).name,'(?<=^\d+-)\d+(?=-)','match');
                        day=str2double(regexp(files(fiidx).name,'(?<=day)\d+(?=-)','match','once'));
                        fs=load(files(fiidx).name);
                        fsn=fieldnames(fs);
                        ftrialIdx=find(endsWith(fsn,'Trial'));
                        if numel(ftrialIdx)>=1
                            lnRows(end+1,1:5)={folders(fidx).name,regions{ridx},miceId{1},prog{1},sprintf('day%d',day)};
                            dual_diff=0;
                            for taskIdx=1:numel(ftrialIdx)
                                perfT=fs.(fsn{ftrialIdx(taskIdx)});
                                expType=strrep(fsn{ftrialIdx(taskIdx)},'Trial','');
                                if strcmp(expType,'DPA') && (strcmp(prog{1},'4351') || strcmp(prog{1},'4353') || strcmp(prog{1},'4341'))
                                    lnRows{end,5+taskIdx*3-2}='Simple_DPA';
                                    perfOff=sum(ismember(perfT(perfT(:,4)==0,6),[5,7]))./sum(perfT(:,4)==0);
                                    perfOn=sum(ismember(perfT(perfT(:,4)==1,6),[5,7]))./sum(perfT(:,4)==1);
                                    lnRows{end,5+taskIdx*3-1}=perfOff;
                                    lnRows{end,5+taskIdx*3}=perfOn;
                                elseif strcmp(expType,'DPA') && (strcmp(prog{1},'4364') || strcmp(prog{1},'4363'))
                                    lnRows{end,5+taskIdx*3-2}='Dual_DPA_S';
                                    perfOffS=sum(ismember(perfT(perfT(:,4)==0 & perfT(:,8)==0,6),[5,7]))./sum(perfT(:,4)==0 & perfT(:,8)==0);
                                    perfOnS=sum(ismember(perfT(perfT(:,4)==1 & perfT(:,8)==0,6),[5,7]))./sum(perfT(:,4)==1 & perfT(:,8)==0);
                                    lnRows{end,5+taskIdx*3-1}=perfOffS;
                                    lnRows{end,5+taskIdx*3}=perfOnS;
                                    dual_diff=1;
                                    taskIdx=taskIdx+dual_diff;
                                    lnRows{end,5+taskIdx*3-2}='Dual_DPA_D';
                                    perfOffD=sum(ismember(perfT(perfT(:,4)==0 & perfT(:,8)==1,6),[5,7]))./sum(perfT(:,4)==0 & perfT(:,8)==1);
                                    perfOnD=sum(ismember(perfT(perfT(:,4)==1 & perfT(:,8)==1,6),[5,7]))./sum(perfT(:,4)==1 & perfT(:,8)==1);
                                    lnRows{end,5+taskIdx*3-1}=perfOffD;
                                    lnRows{end,5+taskIdx*3}=perfOnD;
                                elseif contains(expType,'DRT') && strcmp(prog{1},'4341')
                                    lnRows{end,5+taskIdx*3-2}='Simple_DRT';
                                    perfOff=sum(ismember(perfT(perfT(:,4)==0,6),[5,7]))./sum(perfT(:,4)==0);
                                    perfOn=sum(ismember(perfT(perfT(:,4)==1,6),[5,7]))./sum(perfT(:,4)==1);
                                    lnRows{end,5+taskIdx*3-1}=perfOff;
                                    lnRows{end,5+taskIdx*3}=perfOn;
                                elseif contains(expType,'DRT') && (strcmp(prog{1},'4364') || strcmp(prog{1},'4363'))
                                    taskIdx=taskIdx+dual_diff;
                                    lnRows{end,5+taskIdx*3-2}='Dual_DRT';
                                    perfOff=sum(ismember(perfT(perfT(:,7)==0,5),[5,7]))./sum(perfT(:,7)==0);
                                    perfOn=sum(ismember(perfT(perfT(:,7)==1,5),[5,7]))./sum(perfT(:,7)==1);
                                    lnRows{end,5+taskIdx*3-1}=perfOff;
                                    lnRows{end,5+taskIdx*3}=perfOn;
                                else
                                    disp('Unexpected design');
                                    perf=-1;
                                end
                                
                            end
                        end
                        
                    end
                    cd('..');
                end
                cd('..');
            end
        cd('..');    
        end
        
    end
end
save('WT_LN.mat','wtRows','lnRows','wtRowsHitMiss','wtRowsCRFalse');




function out=clearBadPerf(facSeq)

    if length(facSeq)>=80
        facSeq(:,8)=0;
        i=80;
        while i<=length(facSeq)
            goodOff=nnz((facSeq(i-79:i,6)==5 | facSeq(i-79:i,6)==7)& facSeq(i-79:i,4)==0);
            if goodOff>=nnz(facSeq(i-79:i,3)==0)*80/100
                facSeq(i-79:i,8)=1;
            end
            i=i+1;
        end
        out=facSeq(facSeq(:,8)==1,:);
    else
        out=[];
    end
end


