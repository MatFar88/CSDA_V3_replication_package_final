function fit = regularized_trimmed_fit_corrected_v1(X,init,r,c,alpha,method,maxit,tauW,tauT)
%REGULARIZED_TRIMMED_FIT_CORRECTED_V1 Common TFKM/TRKM alternating fit.
%
% Setting tauW=tauT=0 yields the ORIGINAL method. Positive values yield the
% HD regularized counterpart. No pseudoinverse, null-space projection or
% hidden fallback is used for the original method.

method = string(method);

n = size(X,1);
p = size(X,2);
g = init.labels(:);
trim = logical(init.trim(:));
U = [];
labels_changed_first = NaN;
trim_changed_first = NaN;
objective = NaN;
converged = false;

for it = 1:maxit
    oldg = g;
    oldtrim = trim;
    in = ~trim;

    if sum(in) <= max(c,r)
        error('genericFixed:TooFewInliers', ...
            'Only %d current inliers for c=%d and r=%d.',sum(in),c,r);
    end

    mu = mean(X(in,:),1);
    Xc = X - mu;

    rows = find(in);
    G = sparse(rows,g(in),1,n,c);
    nk = full(sum(G,1));
    if any(nk == 0)
        Yin = Xc(in,:);
        [~,farord] = sort(sum(Yin.^2,2),'descend');
        missing = find(nk == 0);
        for j = 1:numel(missing)
            g(rows(farord(j))) = missing(j);
        end
        G = sparse(rows,g(in),1,n,c);
        nk = full(sum(G,1));
    end

    P = G * spdiags(1./max(nk(:),1),0,c,c) * G';
    Din = spdiags(double(in),0,n,n);

    SW = Xc' * (Din - P) * Xc;
    ST = Xc' * Din * Xc;
    SB = ST - SW;
    SW = (SW + SW')/2;
    ST = (ST + ST')/2;
    SB = (SB + SB')/2;

    rankW = max(rank(full(SW)),1);
    rankT = max(rank(full(ST)),1);
    lamW = tauW * trace(SW) / rankW;
    lamT = tauT * trace(ST) / rankT;

    if method == "REG_TFKM"
        A = SW + lamW*eye(p);
        B = ST + lamT*eye(p);
        [V,D] = eig(A,B,'vector');
        [U,selected_values] = select_generalized_vectors(V,D,r,"ascend");
    elseif method == "REG_TRKM"
        A = SB;
        B = ST + lamT*eye(p);
        [V,D] = eig(A,B,'vector');
        [U,selected_values] = select_generalized_vectors(V,D,r,"descend");
    else
        error('genericFixed:UnknownMethod','Unknown method: %s',method);
    end

    [U,~] = qr(real(U),0);
    if size(U,2) < r || any(~isfinite(U(:)))
        error('genericFixed:InvalidSubspace','Unable to recover %d finite directions.',r);
    end

    Y = Xc*U;
    centers = zeros(c,r);
    for k = 1:c
        q = (g==k & in);
        if any(q)
            centers(k,:) = mean(Y(q,:),1);
        else
            centers(k,:) = Y(rows(randi(numel(rows))),:);
        end
    end

    DD = zeros(n,c);
    for k = 1:c
        Z = Y - centers(k,:);
        DD(:,k) = sum(Z.^2,2);
    end
    [dist,gnew] = min(DD,[],2);

    m = floor(alpha*n);
    trimnew = false(n,1);
    if m > 0
        [~,ix] = maxk(dist,m);
        trimnew(ix) = true;
    end

    if it == 1
        labels_changed_first = sum(gnew ~= oldg);
        trim_changed_first = sum(trimnew ~= oldtrim);
    end

    g = gnew;
    trim = trimnew;
    objective = sum(dist(~trim));

    if isequal(g,oldg) && isequal(trim,oldtrim)
        converged = true;
        break
    end
end

fit = struct('labels',g,'trim',trim,'U',U,'iterations',it, ...
    'labels_changed_first',labels_changed_first, ...
    'trim_changed_first',trim_changed_first, ...
    'objective',objective,'converged',converged, ...
    'selected_eigenvalues',selected_values(:)');
end

function [U,vals] = select_generalized_vectors(V,D,r,direction)
D = real(D(:));
finite_cols = isfinite(D) & all(isfinite(real(V)),1)';
idx = find(finite_cols);
if numel(idx) < r
    error('genericFixed:SingularGeneralizedEigenproblem', ...
        'Only %d finite generalized eigenvectors are available; %d required.', ...
        numel(idx),r);
end
if direction == "ascend"
    [~,ord] = sort(D(idx),'ascend');
else
    [~,ord] = sort(D(idx),'descend');
end
sel = idx(ord(1:r));
U = real(V(:,sel));
vals = D(sel);
end
