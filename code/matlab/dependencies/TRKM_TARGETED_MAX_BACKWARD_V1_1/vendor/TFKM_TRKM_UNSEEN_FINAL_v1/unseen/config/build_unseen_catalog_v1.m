function catalog=build_unseen_catalog_v1(options)
%BUILD_UNSEEN_CATALOG_V1 Fixed four-cell unseen confirmation catalog.
% Same four DGP definitions as development, but NEW scenario IDs and seeds.
% These values are frozen prospectively before any unseen result is inspected.

full=build_generic_fixed_catalog_v1(options);
rows=table();
spec={ ...
 33001,.5,2,3,"TFKM_favorable","moderate","clean",    0,   "GENERAL_PRIMARY"; ...
 33004,.5,3,4,"TRKM_favorable","moderate","bridge",  .10, "GENERAL_PRIMARY"; ...
 33005,2,  4,5,"TRKM_favorable","moderate","leverage",.10,"GENERAL_PRIMARY"; ...
 33007,.5,3,4,"TFKM_favorable","perfect", "leverage",.05, "SENTINEL_LEVEL"};

for j=1:size(spec,1)
    sid=spec{j,1}; pn=spec{j,2}; r=spec{j,3}; g=spec{j,4};
    arm=string(spec{j,5}); diff=string(spec{j,6}); cont=string(spec{j,7});
    a=spec{j,8}; role=string(spec{j,9});
    b=pick_base(full,pn,r,g,arm,diff,cont);
    q=b;
    q.base_scenario_id=b.scenario_id;
    q.scenario_id=sid;
    q.alpha0=a;
    q.joint_role=role;
    q.source_confirm_id=sid-10000; % development analogue: 23001/4/5/7
    q.dgp_recalibrated=(sid==33001);
    if sid==33001
        q.delta=10;
        q.rho=2;
    end
    rows=[rows;q]; %#ok<AGROW>
end

catalog=rows;
catalog.primary_tfkm=string(catalog.joint_role)=="GENERAL_PRIMARY";
catalog.primary_trkm=string(catalog.joint_role)=="GENERAL_PRIMARY";
catalog.sentinel_tfkm=string(catalog.joint_role)=="SENTINEL_LEVEL";
catalog.sentinel_trkm=string(catalog.joint_role)=="SENTINEL_LEVEL";
assert(height(catalog)==4 && numel(unique(catalog.scenario_id))==4, ...
    'Unseen catalog construction failed.');
end

function b=pick_base(full,pn,r,g,arm,diff,cont)
q=abs(double(full.p_over_n)-pn)<1e-12 & double(full.r0)==r & double(full.c0)==g & ...
  string(full.arm)==arm & string(full.difficulty)==diff & string(full.contamination)==cont;
b=full(q,:);
assert(height(b)==1,'Expected one frozen registry cell.');
end
