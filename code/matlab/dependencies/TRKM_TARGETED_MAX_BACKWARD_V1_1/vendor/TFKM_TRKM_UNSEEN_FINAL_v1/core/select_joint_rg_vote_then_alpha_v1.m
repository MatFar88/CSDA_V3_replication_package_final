function sel=select_joint_rg_vote_then_alpha_v1(C,method)
%SELECT_JOINT_RG_VOTE_THEN_ALPHA_V1 Frozen fresh-confirmation composition.
% 1) At each alpha in the full grid, apply the frozen method-specific RG selector.
% 2) Choose the modal (r,g) across alpha. Frequency ties use the already used
%    parsimony convention: min(r+g), then smaller g, then smaller r.
% 3) On the selected fixed (r,g) alpha path, apply the frozen alpha selector.
% No truth quantities are read and no new numerical tuning parameter is introduced.
method=string(method);C=C(string(C.method)==method,:);assert(~isempty(C),'No candidates for method %s.',method);
alpha_grid=sort(unique(double(C.alpha)));pairs=unique(C(:,{'r','g'}),'rows');
assert(numel(alpha_grid)>=3,'Need an alpha grid.');assert(height(pairs)>=3,'Need RG candidates.');
R=table();
for i=1:numel(alpha_grid)
 a=alpha_grid(i);Q=C(abs(double(C.alpha)-a)<1e-12,:);rs=select_rg_frozen_v1(Q,method);row=Q(rs.row_index,:);
 R=[R;table(a,double(row.r),double(row.g),string(rs.rule),'VariableNames',{'alpha','r','g','rg_rule'})]; %#ok<AGROW>
end
[rr,gg,votes,vote_fraction]=modal_rg(R);
Q=sortrows(C(C.r==rr & C.g==gg,:),'alpha');as=select_alpha_performance_frozen_v1(Q,method);row=Q(as.row_index,:);
% Diagnostics only: whether final solution is also RG-selected at selected alpha.
k=find(abs(R.alpha-double(row.alpha))<1e-12,1);rg_consistent=~isempty(k)&&R.r(k)==rr&&R.g(k)==gg;
sel=struct('method',method,'r',rr,'g',gg,'alpha',double(row.alpha),'rg_rule',string(R.rg_rule(1)), ...
 'alpha_rule',string(as.rule),'selected_row',row,'rg_by_alpha',R,'vote_count',votes,'vote_fraction',vote_fraction, ...
 'rg_consistent_at_selected_alpha',rg_consistent,'composition_rule',"RG_VOTE_THEN_ALPHA");
end

function [r,g,votes,vote_fraction]=modal_rg(R)
key=string(R.r)+"|"+string(R.g);u=unique(key);cnt=zeros(numel(u),1);
for i=1:numel(u),cnt(i)=sum(key==u(i));end
mx=max(cnt);cand=find(cnt==mx);vals=zeros(numel(cand),3);
for j=1:numel(cand)
 k=find(key==u(cand(j)),1);rr=R.r(k);gg=R.g(k);vals(j,:)=[rr+gg gg rr];
end
[~,o]=sortrows(vals,[1 2 3]);chosen=u(cand(o(1)));k=find(key==chosen,1);r=R.r(k);g=R.g(k);
votes=mx;vote_fraction=mx/height(R);
end
