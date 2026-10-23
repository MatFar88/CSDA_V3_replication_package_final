function out=select_final_pipeline_v1(U,method)
%SELECT_FINAL_PIPELINE_V1 Frozen final structural pipeline, PATCH1.
% Stage1: TFKM CURRENT + GMAX_R_E guard; TRKM CURRENT.
% Stage2: g strict W Pareto, then r strict R2/native Pareto.
% PATCH1 TRKM only: one non-iterated native-010 g safeguard.
% Stage3: frozen alpha selector along the post-structural pair.
method=string(method);
assert(height(U)>0,'Empty candidate cube unit.');
[s1r,s1g,s1a,s1row,guarded,anchor_r]=stage1_start(U,method);
[g_r,g_g,g_a,g_row,gchg]=g_strict_w(U,s1r,s1g,s1a);
[gr_r,gr_g,gr_a,gr_row,rchg]=r_strict_nr2(U,method,g_r,g_g,g_a);

% PATCH1. The safeguard was pre-specified and accepted on the consumed
% Confirm4 cube only. It is deliberately ONE STEP: no cumulative .01 loss.
[sg_r,sg_g,sg_a,sg_row,sg_trigger,sg_dR2,sg_dCV]=trkm_native010_one_step(U,method,gr_r,gr_g,gr_a);

P=sortrows(U(double(U.r)==sg_r & double(U.g)==sg_g,:),'alpha');
assert(~isempty(P),'Final pair absent from cube.');
as=select_alpha_performance_frozen_v1(P,method); final_row=P(as.row_index,:);
out=struct('method',method,'r',sg_r,'g',sg_g,'alpha',double(final_row.alpha),'selected_row',final_row, ...
 'stage1_r',s1r,'stage1_g',s1g,'stage1_alpha',s1a,'stage1_row',s1row, ...
 'stage1_guarded',guarded,'stage1_anchor_r',anchor_r, ...
 'stage2g_r',g_r,'stage2g_g',g_g,'stage2g_alpha',g_a,'stage2g_row',g_row,'g_changed',gchg, ...
 'stage2gr_r',gr_r,'stage2gr_g',gr_g,'stage2gr_alpha',gr_a,'stage2gr_row',gr_row,'r_changed',rchg, ...
 'native010_one_step_triggered',sg_trigger,'native010_r',sg_r,'native010_g',sg_g, ...
 'native010_alpha_ref',sg_a,'native010_row',sg_row,'native010_delta_R2',sg_dR2,'native010_delta_cv',sg_dCV, ...
 'alpha_selector',as);
end

function [rr,gg,aa,row,guarded,anchor_r]=stage1_start(U,method)
s=select_joint_rg_vote_then_alpha_v1(U,method); rr0=double(s.r); gg=double(s.g); guarded=false; anchor_r=rr0;
if method=="REG_TFKM_HD"
    gmax=max(double(U.g(logical(U.valid)))); rr=rr0;
    if gg==gmax
        guarded=true; av=sort(unique(double(U.alpha))); AR=zeros(numel(av),1);
        for i=1:numel(av)
            Q=U(abs(double(U.alpha)-av(i))<1e-12 & double(U.g)==gg,:);
            valid=logical(Q.valid)&isfinite(double(Q.W_ratio))&isfinite(double(Q.R2)); assert(any(valid),'No valid rank-anchor candidates.');
            x=double(Q.W_ratio).*double(Q.R2); x(~valid)=-Inf; mx=max(x);
            cand=find(valid & abs(x-mx)<=1e-12); AR(i)=min(double(Q.r(cand)));
        end
        anchor_r=modal_scalar_parsimonious(AR); rr=max(rr0,anchor_r);
    end
    P=sortrows(U(double(U.r)==rr & double(U.g)==gg,:),'alpha'); assert(~isempty(P),'Guarded Stage1 pair absent.');
    as=select_alpha_performance_frozen_v1(P,method); row=P(as.row_index,:); aa=double(row.alpha);
else
    rr=rr0; aa=double(s.alpha); row=s.selected_row;
end
end

function [rr,gg,aa,row,changed]=g_strict_w(U,sr,sg,sa)
ref=cube_row(U,sr,sg,sa); rr=sr; gg=sg; aa=sa; gmin=min(double(U.g));
for gc=(sg-1):-1:gmin
    rcap=min(sr,gc-1);
    Q=U(double(U.g)==gc & abs(double(U.alpha)-sa)<1e-12 & double(U.r)<=rcap & logical(U.valid),:);
    Q=Q(isfinite(double(Q.W_raw)) & isfinite(double(Q.W_ratio)),:); if isempty(Q),continue;end
    rc=max(double(Q.r)); cand=cube_row(U,rc,gc,sa);
    if double(cand.W_raw)<=double(ref.W_raw)+1e-12 && double(cand.W_ratio)>=double(ref.W_ratio)-1e-12
        rr=rc; gg=gc;
    end
end
row=cube_row(U,rr,gg,aa); changed=(rr~=sr)||(gg~=sg);
end

function [rr,gg,aa,row,changed]=r_strict_nr2(U,method,sr,sg,sa)
ref=cube_row(U,sr,sg,sa); rr=sr; gg=sg; aa=sa;
if method=="REG_TFKM_HD", n0=double(ref.selector_score); else, n0=double(ref.selector_cv_score); end
rmin=min(double(U.r));
for rc=(sr-1):-1:rmin
    if rc>=sg,continue;end
    cand=cube_row(U,rc,sg,sa); if ~logical(cand.valid),continue;end
    if method=="REG_TFKM_HD", nc=double(cand.selector_score); else, nc=double(cand.selector_cv_score); end
    if isfinite(double(cand.R2)) && isfinite(nc) && isfinite(n0) && double(cand.R2)>=double(ref.R2)-1e-12 && nc>=n0-1e-12
        rr=rc;
    end
end
row=cube_row(U,rr,gg,aa); changed=(rr~=sr);
end

function [rr,gg,aa,row,triggered,dR2,dCV]=trkm_native010_one_step(U,method,sr,sg,sa)
% Exactly the ONE_STEP rule that passed the consumed Confirm4 audit:
%   candidate g = g-1, r = min(r,g_candidate-1)
%   accept iff R2 loss <= .01 and RV_CV does not decrease.
% It is never iterated.
rr=sr; gg=sg; aa=sa; row=cube_row(U,sr,sg,sa); triggered=false; dR2=NaN; dCV=NaN;
if method~="REG_TRKM_HD", return; end
gmin=min(double(U.g)); if sg<=gmin, return; end
gc=sg-1; rc=min(sr,gc-1);
ref=row; cand=cube_row(U,rc,gc,sa);
valid=logical(cand.valid) && ...
    isfinite(double(ref.R2)) && isfinite(double(cand.R2)) && ...
    isfinite(double(ref.selector_cv_score)) && isfinite(double(cand.selector_cv_score));
if ~valid, return; end
dR2=double(cand.R2)-double(ref.R2);
dCV=double(cand.selector_cv_score)-double(ref.selector_cv_score);
pass=double(cand.R2)>=double(ref.R2)-0.01-1e-12 && ...
     double(cand.selector_cv_score)>=double(ref.selector_cv_score)-1e-12;
if pass
    rr=rc; gg=gc; row=cand; triggered=true;
end
end

function row=cube_row(U,r,g,a)
idx=double(U.r)==r & double(U.g)==g & abs(double(U.alpha)-a)<1e-12;
assert(sum(idx)==1,'Cube row not unique/absent for (%d,%d,%.6f).',r,g,a); row=U(idx,:);
end
function a=modal_scalar_parsimonious(x)
u=unique(x); cnt=zeros(numel(u),1); for i=1:numel(u),cnt(i)=sum(x==u(i));end
mx=max(cnt); a=min(u(cnt==mx));
end
