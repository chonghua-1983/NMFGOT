% preprocess microbial text and taxa data
% X: row represents samples, columns represent features
function proData = filter_feas(X,smps,feas, minCounts, minSmps, flag)
addpath('external/')
% filter features
% X = tfidf2(X);  % - word count vectors (one column = one document)
if flag == 1
    X = log(X+1);
end

dataTemp = X;
dataTemp(X > 0) = 1;
msum = sum(dataTemp,1);
idx1 = find(msum <= minCounts);  % filter features
X(:,idx1) = [];  feas(idx1) = [];

% filter samples
% minSmps = 50;
dataTemp = X;
dataTemp(X > 0) = 1;
msum = sum(dataTemp,2);
idx2 = find(msum <= minSmps);
X(idx2,:) = []; smps(idx2,:) = [];

proData.data = X;
proData.Features = feas;
proData.samples = smps;

end

