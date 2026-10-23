function out = score_candidates_robust_views_v3(pool,X,c,r)
%SCORE_CANDIDATES_ROBUST_VIEWS_V3 Truth-blind candidate validation scores.
%
% cv_score: 5-fold centroid self-prediction coherence, averaged across
% SIGN-PCA and WINSOR-PCA robust views.
% geom_score: between/total trace fraction, averaged across the same views.
%
% Candidate truth labels, outlier indicators and candidate ARI are never read.

views = build_robust_views_v3(X,r,c);
valid = logical(pool.valid(:));
N = numel(valid); n=size(X,1);
cv = nan(N,1); geom=nan(N,1);
view_cv = nan(N,numel(views)); view_geom=nan(N,numel(views));
fold = mod((1:n)'-1,5)+1; % deterministic and data-independent

for i=find(valid)'
    lab = double(pool.labels(:,i)); tr = logical(pool.trim(:,i)); kept=~tr;
    cvs=nan(numel(views),1); gs=nan(numel(views),1);
    for v=1:numel(views)
        Z=views{v};
        % Geometry: fraction of kept-data variation explained by clustering.
        Zk=Z(kept,:); lk=lab(kept);
        if size(Zk,1)<=c || numel(unique(lk))<c
            gs(v)=-Inf; cvs(v)=-Inf; continue;
        end
        mu=mean(Zk,1); tss=sum((Zk-mu).^2,'all'); wss=0;
        for k=1:c
            q=(lk==k);
            if ~any(q), wss=Inf; break; end
            mk=mean(Zk(q,:),1); wss=wss+sum((Zk(q,:)-mk).^2,'all');
        end
        if isfinite(wss) && tss>eps, gs(v)=1-wss/tss; else, gs(v)=-Inf; end

        % Cross-fold coherence: centroid model is learned without the test fold.
        fold_scores=nan(5,1);
        for f=1:5
            train=kept & fold~=f; test=kept & fold==f;
            if sum(test)<2 || sum(train)<=c, continue; end
            centers=zeros(c,size(Z,2)); ok=true;
            for k=1:c
                q=train & lab==k;
                if ~any(q), ok=false; break; end
                centers(k,:)=mean(Z(q,:),1);
            end
            if ~ok, fold_scores(f)=0; continue; end
            D=zeros(sum(test),c); Zt=Z(test,:);
            for k=1:c
                R=Zt-centers(k,:); D(:,k)=sum(R.^2,2);
            end
            [~,pred]=min(D,[],2);
            % Direct agreement is valid because centroids inherit candidate labels.
            fold_scores(f)=mean(pred==lab(test));
        end
        cvs(v)=mean(fold_scores,'omitnan');
    end
    view_cv(i,:)=cvs; view_geom(i,:)=gs;
    cv(i)=mean(cvs(isfinite(cvs)),'omitnan');
    geom(i)=mean(gs(isfinite(gs)),'omitnan');
end
out=struct('cv_score',cv,'geom_score',geom,'view_cv',view_cv,'view_geom',view_geom);
end
