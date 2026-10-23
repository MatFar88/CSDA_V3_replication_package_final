function fit = normalize_fit_contract_v1(fit,n,c,alpha,method)
% Normalize labels and trimming to one strict contract for every method.
% labels: n-by-1 integer in 0:c; zero is permitted only for trimmed units.
% trim:   n-by-1 logical; exactly floor(alpha*n) true values.
method = string(method);
assert(isfield(fit,'labels'),'Missing labels for %s.',method);
labels = fit.labels(:);
assert(numel(labels)==n && isnumeric(labels) && all(isfinite(labels)), ...
    'Invalid labels returned by %s.',method);
if isfield(fit,'trim') && ~isempty(fit.trim)
    trim = logical(fit.trim(:));
elseif isfield(fit,'outlier_indices')
    trim = false(n,1);
    idx = fit.outlier_indices(:);
    if islogical(idx) && numel(idx)==n
        trim = idx;
    elseif ~isempty(idx)
        assert(all(isfinite(idx)) && all(idx==floor(idx)) && ...
            all(idx>=1) && all(idx<=n), ...
            'Invalid outlier indices returned by %s.',method);
        trim(idx)=true;
    end
else
    trim = labels==0;
end
assert(numel(trim)==n,'Invalid trim length returned by %s.',method);
expected = floor(alpha*n);
assert(sum(trim)==expected, ...
    '%s returned %d trimmed observations; expected %d.', ...
    method,sum(trim),expected);
% Preserve a full membership for metric computation. For methods that
% expose zero labels only on trimmed observations, zero is retained.
assert(all(labels==floor(labels)) && all(labels>=0) && all(labels<=c), ...
    'Invalid label values returned by %s.',method);
assert(all(labels(~trim)>=1), ...
    '%s assigned zero to a retained observation.',method);
fit.labels = labels;
fit.trim = trim;
fit.trim_count = sum(trim);
fit.trim_contract = "LOGICAL_MASK_V1";
fit.method = method;
end
