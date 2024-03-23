% parameter selection for OrthNMF algorithm
function params = parameter_selection_mfop(X1, X2, WD1, WD2, num_factor)
% WD1,WD2: WMD similairty, equals to be :1-WMD
[W1,H1] = nndsvd(X1,num_factor,0); 
[W2,H2] = nndsvd(X2,num_factor,0);
D1 = dist2(X1',X1'); A1 = affinityMatrix(D1,20); 
D2 = dist2(X2',X2'); A2 = affinityMatrix(D2,20); 
S = SNF({A1,A2},20);

n = size(X1,2);
L1 = diag(sum(WD1,1)) - WD1; 
L2 = diag(sum(WD2,1)) - WD2; 

H1tH1 = H1'*H1; H1tH1 = H1tH1./max(H1tH1(:));
H2tH2 = H2'*H2; H2tH2 = H2tH2./max(H2tH2(:));
error_1 = norm(X1-W1*H1,'fro')^2+norm(X2-W2*H2,'fro')^2;
error_2 = norm(S-H1'*H1,'fro')^2+norm(S-H2'*H2,'fro')^2;
error_3 = trace(H1tH1 *(1 - H2tH2));
error_4 = trace(H1*L1*H1'+H2*L2*H2');
error_5 = norm(S*ones(n,1)-ones(n,1));

alpha = error_1/error_2;
% beta = abs(error_1/error_3);
gamma = abs(error_1/error_4);
phi = error_1/error_5;

% demo_1
params.alpha = alpha/5;
%params.beta = beta/100;
params.phi = phi/1000;
params.gamma = gamma/1000;

end
