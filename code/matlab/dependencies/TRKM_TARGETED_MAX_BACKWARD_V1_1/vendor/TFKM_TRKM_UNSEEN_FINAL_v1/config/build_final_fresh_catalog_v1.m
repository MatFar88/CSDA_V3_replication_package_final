function catalog=build_final_fresh_catalog_v1(options)
%BUILD_FINAL_FRESH_CATALOG_V1 Fresh four-cell validation catalog.
% New scenario IDs deliberately break every old seed path.
% The clean TFKM-favourable moderate cell is the pre-specified DGP repair:
% delta=10 unchanged, rho=2 instead of the old rho=8.
full=build_generic_fixed_catalog_v1(options);
rows=table();
spec={ ...
 23001,20001,.5,2,3,"TFKM_favorable","moderate","clean",    0,   "GENERAL_PRIMARY"; ...
 23004,20004,.5,3,4,"TRKM_favorable","moderate","bridge",  .10, "GENERAL_PRIMARY"; ...
 23005,20005,2,  4,5,"TRKM_favorable","moderate","leverage",.10,"GENERAL_PRIMARY"; ...
 23007,20007,.5,3,4,"TFKM_favorable","perfect", "leverage",.05, "SENTINEL_LEVEL"};
for j=1:size(spec,1)
    sid=spec{j,1}; src=spec{j,2}; pn=spec{j,3}; r=spec{j,4}; g=spec{j,5};
    arm=string(spec{j,6}); diff=string(spec{j,7}); cont=string(spec{j,8});
    a=spec{j,9}; role=string(spec{j,10});
    b=pick_base(full,pn,r,g,arm,diff,cont);
    q=b; q.base_scenario_id=b.scenario_id; q.scenario_id=sid; q.alpha0=a;
    q.joint_role=role; q.source_confirm_id=src; q.dgp_recalibrated=(sid==23001);
    if sid==23001
        q.delta=10; q.rho=2;
    end
    rows=[rows;q]; %#ok<AGROW>
end
catalog=rows;
catalog.primary_tfkm=string(catalog.joint_role)=="GENERAL_PRIMARY";
catalog.primary_trkm=string(catalog.joint_role)=="GENERAL_PRIMARY";
catalog.sentinel_tfkm=string(catalog.joint_role)=="SENTINEL_LEVEL";
catalog.sentinel_trkm=string(catalog.joint_role)=="SENTINEL_LEVEL";
assert(height(catalog)==4 && numel(unique(catalog.scenario_id))==4,'Fresh catalog construction failed.');
end

function b=pick_base(full,pn,r,g,arm,diff,cont)
q=abs(double(full.p_over_n)-pn)<1e-12 & double(full.r0)==r & double(full.c0)==g & ...
  string(full.arm)==arm & string(full.difficulty)==diff & string(full.contamination)==cont;
b=full(q,:); assert(height(b)==1,'Expected one frozen registry cell.');
end
