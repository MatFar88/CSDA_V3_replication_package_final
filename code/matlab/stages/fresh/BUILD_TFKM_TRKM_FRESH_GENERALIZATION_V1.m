function OUT = BUILD_TFKM_TRKM_FRESH_GENERALIZATION_V1(documents_root)
%BUILD_TFKM_TRKM_FRESH_GENERALIZATION_V1
% Build the immutable fresh H2/generalization study specification.
%
% Design frozen before outcomes:
%   n = 120. Main geometry uses r0 = 3, g0 = 4.
%   alpha0 = 0 (clean), 0.10 (bridge/leverage)
%   MAIN: p/n in {0.5,1,2}, rho in {0.25,1,2}, delta=10,
%         contamination in {clean,bridge,leverage}: 27 cells.
%   WEAK-SEPARATION: leverage only, delta=7, at
%         (p/n,rho)={(1,1),(1,2),(2,.25),(2,1),(2,2)}: 5 cells.
%   STRUCTURE: p/n=1, rho=1, delta=10, all contamination types,
%         truth pairs (r0,g0) in {(2,3),(3,5),(4,5)}: 9 cells.
%   Total = 41 cells x 30 pre-specified replicates = 1230 units.
%
% This builder creates NO method outcomes and opens NO truth performance.
% It only writes study.mat, STUDY_CATALOG.csv, UNIT_MANIFEST.csv, PROTOCOL.txt.
%
% Usage:
%   S = BUILD_TFKM_TRKM_FRESH_GENERALIZATION_V1;
%   S = BUILD_TFKM_TRKM_FRESH_GENERALIZATION_V1("C:\Users\name\Documents");
%   BUILD_TFKM_TRKM_FRESH_GENERALIZATION_V1("SELFTEST");

if nargin==1 && strcmpi(string(documents_root),"SELFTEST")
    [catalog,manifest,cfg]=make_spec_local();
    assert(height(catalog)==41);
    assert(height(manifest)==1230);
    assert(sum(catalog.panel=="MAIN" & catalog.delta==10)==27);
    assert(sum(catalog.delta==7)==5);
    assert(all(catalog.n==120));
    assert(all(ismember(catalog.p,[60 120 240])));
    assert(all(catalog.g0==catalog.c0));
    assert(sum(catalog.panel=="STRUCTURE")==9);
    Cstruct=catalog(catalog.panel=="STRUCTURE",:);
    Pstruct=unique([double(Cstruct.r0),double(Cstruct.g0)],"rows");
    assert(isequal(sortrows(Pstruct),[2 3;3 5;4 5]));
    assert(numel(unique(manifest.data_seed))==1230);
    fprintf('SELFTEST PASS: 41-cell fresh generalization design + 1230 fixed seeds + c0/g0 DGP contract.\n');
    OUT=struct('selftest_pass',true,'config',cfg);
    return
end
if nargin<1, documents_root=""; end
[root,template]=resolve_root_and_template_local(documents_root);
S=load(template);
assert(isfield(S,'options'),'Template study.mat lacks options.');
options=S.options;
validate_options_local(options);
[catalog,manifest,CFG]=make_spec_local();

outroot=fullfile(root,'TFKM_TRKM_FRESH_GENERALIZATION_V1_output');
if ~isfolder(outroot),mkdir(outroot);end
writetable(catalog,fullfile(outroot,'STUDY_CATALOG.csv'));
writetable(manifest,fullfile(outroot,'UNIT_MANIFEST.csv'));
save(fullfile(outroot,'study.mat'),'options','catalog','manifest','CFG','-v7.3');
write_protocol_local(fullfile(outroot,'PROTOCOL.txt'),CFG,catalog,manifest,template);

fprintf('\n============================================================\n');
fprintf('TFKM/TRKM FRESH GENERALIZATION V1 - DESIGN BUILT\n');
fprintf('41 cells x 30 fixed replicates = 1230 units\n');
fprintf('n=120 | p={60,120,240} | p/n={.5,1,2}\n');
fprintf('rho={.25,1,2} | delta=10 main + delta=7 leverage stress + structure panel\n');
fprintf('No method outcomes generated.\n');
fprintf('Output: %s\n',outroot);
fprintf('============================================================\n');

OUT=struct('catalog',catalog,'manifest',manifest,'config',CFG, ...
    'output_root',string(outroot),'template_study',string(template));
end

function [catalog,manifest,CFG]=make_spec_local()
CFG=struct();
CFG.study_name="TFKM_TRKM_FRESH_GENERALIZATION_V1";
CFG.n=120;CFG.r0=3;CFG.g0=4;CFG.alpha_leverage=.10;CFG.alpha_bridge=.10;
CFG.structure_pairs=[2 3;3 5;4 5];CFG.structure_ratio=1;CFG.structure_rho=1;CFG.structure_delta=10;
CFG.ratios=[.5 1 2];CFG.rhos=[.25 1 2];CFG.delta_main=10;CFG.delta_weak=7;
CFG.alpha_grid=[0 .025 .05 .075 .10 .125 .15 .20 .25 .30];
CFG.n_replicates=30;
CFG.seed_base=314159265; % independent namespace, safely below uint32 bound
CFG.tfkm_pair_offset=110000;
CFG.trkm_pair_offset=220000;
CFG.primary_batch=1:10;
CFG.extension_batch=11:30;
CFG.truth_firewall="Do not compute/print truth performance during FIT phase.";

rows=cell(0,11);k=0;
conts=["clean","bridge","leverage"];
for ratio=CFG.ratios
    p=round(CFG.n*ratio);
    for rho=CFG.rhos
        for cont=conts
            k=k+1;sid=51000+k;
            a0=0;if cont~="clean",a0=.10;end
            rows(end+1,:)={sid,"MAIN",CFG.n,p,ratio,CFG.r0,CFG.g0,CFG.delta_main,rho,cont,a0}; %#ok<AGROW>
        end
    end
end
weak=[1 1;1 2;2 .25;2 1;2 2];
for j=1:size(weak,1)
    ratio=weak(j,1);rho=weak(j,2);p=round(CFG.n*ratio);sid=51100+j;
    rows(end+1,:)={sid,"WEAK_DELTA7",CFG.n,p,ratio,CFG.r0,CFG.g0,CFG.delta_weak,rho,"leverage",.10}; %#ok<AGROW>
end
for j=1:size(CFG.structure_pairs,1)
    r0=CFG.structure_pairs(j,1);g0=CFG.structure_pairs(j,2);ratio=CFG.structure_ratio;rho=CFG.structure_rho;p=round(CFG.n*ratio);
    for cont=conts
        k=k+1;sid=51200+k; a0=0;if cont~="clean",a0=.10;end
        rows(end+1,:)={sid,"STRUCTURE",CFG.n,p,ratio,r0,g0,CFG.structure_delta,rho,cont,a0}; %#ok<AGROW>
    end
end
catalog=cell2table(rows,'VariableNames',{'scenario_id','panel','n','p','p_over_n','r0','g0','delta','rho','contamination','alpha0'});
catalog.scenario_id=double(catalog.scenario_id);catalog.n=double(catalog.n);catalog.p=double(catalog.p);
catalog.p_over_n=double(catalog.p_over_n);catalog.r0=double(catalog.r0);catalog.g0=double(catalog.g0);
catalog.c0=catalog.g0; % canonical DGP generator contract; g0 is retained as reporting alias
catalog.delta=double(catalog.delta);catalog.rho=double(catalog.rho);catalog.alpha0=double(catalog.alpha0);
catalog.panel=string(catalog.panel);catalog.contamination=string(catalog.contamination);
assert(height(catalog)==41);

manifest=table();
for i=1:height(catalog)
    sid=catalog.scenario_id(i);
    for rep=1:CFG.n_replicates
        seed=CFG.seed_base + sid*100 + rep;
        assert(seed+CFG.trkm_pair_offset+6000 < 2^32-1,'Seed overflow.');
        batch="EXTENSION";if rep<=10,batch="PRIMARY10";end
        manifest=[manifest;table(sid,rep,seed,seed+CFG.tfkm_pair_offset,seed+CFG.trkm_pair_offset,string(batch), ...
            'VariableNames',{'scenario_id','replicate','data_seed','tfkm_pair_seed','trkm_pair_seed','batch'})]; %#ok<AGROW>
    end
end
manifest=sortrows(manifest,{'scenario_id','replicate'});
end

function validate_options_local(options)
req={'hd_init_names','hd_tfkm_allocation','hd_trkm_allocation','trim_stability_weight'};
for j=1:numel(req),assert(isfield(options,req{j}),'options missing %s.',req{j});end
assert(isequal(string(options.hd_init_names(:)'),["WINSOR_PCA_TKM","SIGN_PCA_TKM","TKM_ORIGINAL","RP_TKM"]),'Initialization family drift.');
assert(isequal(double(options.hd_tfkm_allocation(:)'),[50 35 10 5]),'TFKM allocation drift.');
assert(isequal(double(options.hd_trkm_allocation(:)'),[40 30 15 15]),'TRKM allocation drift.');
assert(abs(double(options.trim_stability_weight)-.2)<=1e-12,'trim_stability_weight drift.');
end

function [root,template]=resolve_root_and_template_local(arg)
if strlength(string(arg))>0,cands=string(arg);else,cands=[string(fullfile(getenv('USERPROFILE'),'Documents'));string(pwd)];end
root="";template="";
for c=cands(:)'
    t1=fullfile(c,'CSDA6TER_PN_PATH_V1_output','study.mat');
    t2=fullfile(c,'TFKM_TRKM_UNSEEN_FINAL_v1_output','study.mat');
    if isfile(t1),root=c;template=string(t1);break;end
    if isfile(t2),root=c;template=string(t2);break;end
end
assert(strlength(root)>0,'Cannot locate template study.mat. Pass Documents explicitly.');
end

function write_protocol_local(f,CFG,catalog,manifest,template)
fid=fopen(f,'w');assert(fid>0);c=onCleanup(@()fclose(fid));
fprintf(fid,'TFKM_TRKM_FRESH_GENERALIZATION_V1\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'Status: FRESH GENERALIZATION / H2 STRESS PANEL.\n');
fprintf(fid,'Procedure must remain frozen. No truth-driven retuning.\n\n');
fprintf(fid,'DGP main truth: n=%d, r0=%d, g0=%d.\n',CFG.n,CFG.r0,CFG.g0);
fprintf(fid,'Main: p/n={0.5,1,2}, rho={0.25,1,2}, delta=10, clean/bridge/leverage.\n');
fprintf(fid,'Weak-separation: delta=7, leverage only, at (1,1),(1,2),(2,.25),(2,1),(2,2).\n');
fprintf(fid,'Structure panel: p/n=1, rho=1, delta=10, clean/bridge/leverage, truth pairs (2,3),(3,5),(4,5).\n');
fprintf(fid,'alpha0=0 clean; alpha0=.10 bridge/leverage.\n');
fprintf(fid,'Replicates: 30 per cell, all seeds fixed now. First 10 and remaining 20 are workload batches only.\n');
fprintf(fid,'Candidate alpha grid: 0,.025,.05,.075,.10,.125,.15,.20,.25,.30.\n');
fprintf(fid,'Candidate structural grid: all 1<=r<g<=6 (15 pairs).\n');
fprintf(fid,'TFKM pair offset=%d; TRKM pair offset=%d.\n',CFG.tfkm_pair_offset,CFG.trkm_pair_offset);
fprintf(fid,'Template options: %s\n',template);
fprintf(fid,'Cells=%d; units=%d.\n',height(catalog),height(manifest));
fprintf(fid,'Truth firewall: %s\n',CFG.truth_firewall);
end
