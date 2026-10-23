function OUT = RUN_ECB_TRKM_ALL_V1(ecb_mat_file,names364_file)
%RUN_ECB_TRKM_ALL_V1 Convenience serial workflow.
if nargin<1,ecb_mat_file="";end
if nargin<2,names364_file="";end
OUT=RUN_ECB_TRKM_APPLICATION_V1("ALL",1,1,ecb_mat_file,names364_file);
end
