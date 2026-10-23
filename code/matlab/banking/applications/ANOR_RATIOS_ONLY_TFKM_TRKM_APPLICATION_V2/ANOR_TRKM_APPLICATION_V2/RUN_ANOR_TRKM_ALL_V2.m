function OUT=RUN_ANOR_TRKM_ALL_V2(varargin)
% Convenience serial wrapper.
args=cell(1,4);args(:)={""};for k=1:min(numel(varargin),4),args{k}=varargin{k};end
OUT=RUN_ANOR_TRKM_APPLICATION_V2("ALL",1,1,args{:});
end
