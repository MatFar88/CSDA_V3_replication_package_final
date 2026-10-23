function catalog = build_generic_fixed_catalog_v1(options)
rows = cell(0,12); id=0;
for pn=options.p_over_n
 for s=1:size(options.structures,1)
  r=options.structures(s,1); c=options.structures(s,2);
  for arm=options.arms
   for diff=options.difficulty
    [delta,rho]=fixed_separation_registry_v1(pn,r,c,arm,diff);
    for q=1:numel(options.contamination)
     id=id+1; cont=options.contamination(q); alpha=options.alpha_by_contamination(q);
     p=round(options.n*pn);
     rows(end+1,:)={id,options.n,p,pn,r,c,arm,diff,delta,rho,cont,alpha}; %#ok<AGROW>
    end
   end
  end
 end
end
catalog=cell2table(rows,'VariableNames',{'scenario_id','n','p','p_over_n','r0','c0','arm','difficulty','delta','rho','contamination','alpha0'});
end
