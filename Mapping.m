classdef Mapping < handle
    methods (Static)
        function out=brainRegionList()
            fs=ls('.\**\*.mat');
        end
        
        function mCount=APC()
            fs=dir('.\**\*.mat');
            sel=find(contains({fs.name},'APC.mat','IgnoreCase',true));
            mCount=sel.'*[0 0];
            
            for i=1:numel(sel)
                fstr=load(fullfile(fs(sel(i)).folder,fs(sel(i)).name));
                mCount(i,1)=sel(i);
                if isfield(fstr,'GeneInfoCol')
                    mCount(i,2)=size(fstr.GeneInfoCol{1},1);
                elseif isfield(fstr,'DPAGeneInfoCol')
                    mCount(i,2)=size(fstr.DPAGeneInfoCol{1},1);
                end
            end
        end
    end
end