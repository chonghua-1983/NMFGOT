function [WMD_D, WMD_S] = computing_WMD(X)
% compting word moving distance
% input: 
% X, data matrix, row represents featrures and column represents samples
% output:
% WMD_D: word moving distance between two samples
% WMD_S: similarity between two samples, 1-WMD_D for simplicity

[id,mu,F_scaled] = HVGs(X,0.05,3.5,0.05,1);
[tmp, I] = sort(F_scaled, 'descend');
if length(I) > 100
    X = X(I(1:100),:);
end
% if length(id)>100
%     X = X(id(1:100),:);
% end

[m, n] = size(X);
WMD_D = zeros(n,n);

X = X./ repmat(sum(X),m,1);
MI = zeros(m,m); 
% computing mutual information
k = 10;
for i = 1:m
    for j = 1:m
        MI(i,j) = mi_cont_cont(X(i,:), X(j,:), k);
    end
end
MI = MI./max(MI(:));
MI = 1 - MI; D = MI - diag(diag(MI));

% X1_norm = X1./ repmat(sum(X1),m,1);
X_norm = X;
for i = 1:n
    Ei = zeros(1,n);
    % x1 = X1(:,i)./sum(X1(:,i));
    x1 = X_norm(:,i);
    for j = (i+1):n
        % x2 = X1(:,j)./sum(X1(:,j));
        x2 = X_norm(:,j);
        [emd,flow] = emd_mex(x1',x2',D);  % MAX_SIG_SIZE 100
        Ei(j) = emd;
    end
    WMD_D(i,:) = Ei;
end 
WMD_D = WMD_D + WMD_D';
WMD_S = 1 - WMD_D;

% save(save_file,'WMD_D')