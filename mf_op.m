% Dissect host gene and microbe association by optimal transport and matrix  factorization
% inputs:
% X1: gene expression matrix
% X2: microbiome abundance matrix
% W1,W2,H1,H2: the basic matirx and coefficient matrix initialized from X1 and X2
% 
% outputs:
% ||Xi-W_i*H_i|| + alpha*||s-H_i'H_i|| + beta*H1'H1*(1-H2'H2) + gamma*tr(Hi*Li*H_i')
% st. S*1 = 1

function [S,W1,W2, H1,H2,obj] = mf_op(X1,X2,A1,A2,Inits,parameters)

alpha = parameters.alpha; 
%beta = parameters.beta; 
gamma = parameters.gamma;
phi = parameters.phi;
% alpha = 0.1; gamma = 1; phi = 2;
beta = 0.1;
% ablation experiments
% beta = 0;
% gamma =0;

W1 = Inits.W1; W2 = Inits.W2; 
H1 = Inits.H1; H2 = Inits.H2; 
S = Inits.S; 
D1 = diag(sum(A1,1)); D2 = diag(sum(A2,1));

[m,n] = size(X1);
obj_old = 1; 
stop_rule = 2; 
clear Inits

disp('iteration starts!');
Maxiter = 2000; objs = zeros(Maxiter,1); 

for iter = 1:Maxiter
    % update W
    H1H1t = H1*H1'; H2H2t = H2*H2';
    W1 = W1.* (X1*H1')./(W1*H1H1t);
    W2 = W2.* (X2*H2')./(W2*H2H2t);
    % updata H
    H1tH1 = H1'*H1; H2tH2 = H2'*H2;
    H1 = H1.*(W1'*X1 + gamma*H1*A1 + alpha*(2*H1*S') + beta*H1*H2tH2)./(W1'*W1*H1 + gamma*H1*D1 + 2*alpha*H1*H1tH1); %+ beta*H1*ones(n)
    H2 = H2.*(W2'*X2 + gamma*H2*A2 + alpha*(2*H2*S') + beta*H2*H1tH1)./(W2'*W2*H2 + gamma*H2*D2 + 2*alpha*H2*H2tH2);
    % update S
    H1tH1 = H1'*H1; H1tH1 = H1tH1./max(H1tH1(:));
    H2tH2 = H2'*H2; H2tH2 = H2tH2./max(H2tH2(:));
    %[W1, H1] = Normalize(W1, H1');
    %[W2, H2] = Normalize(W2, H2');
    S = S.*(phi*ones(n) + alpha*(H1tH1 + H2tH2))./(phi*S*ones(n) + 2*alpha*S);
         
    if stop_rule == 2
        obj = compute_obj(X1,X2,W1,W2,S,H1,H2,A1,A2,D1,D2,alpha,beta,gamma,phi);
        objs(iter,1) = obj;
        if (abs(obj_old-obj)/obj_old < 10^(-6) && iter > 1) || iter == Maxiter
            disp('converged!');
            break;
        end
        obj_old = obj;
    end
    if mod(iter,50) == 0
        disp(['number of iteration:',num2str(iter),'  obj:',num2str(obj)]);
    end
    
end
S = (S+S')/2;
end


function obj = compute_obj(X1,X2,W1,W2,S,H1,H2,A1,A2,D1,D2,alpha,beta,gamma,phi)
  L1 = D1 - A1; 
  L2 = D2 - A2;
  n = size(X1,2);
  error_1 = norm(X1-W1*H1,'fro')^2+norm(X2-W2*H2,'fro')^2;
  error_2 = alpha*(norm(S-H1'*H1,'fro')^2+norm(S-H2'*H2,'fro')^2);
  error_3 = beta*trace(H1'*H1 *(1-H2'*H2));
  error_4 = gamma*(trace(H1*L1*H1'+H2*L2*H2'));
  error_5 = phi*norm(S*ones(n,1)-ones(n,1));
  obj = error_1 + error_2 - error_3 + error_4 + error_5;
end

function obj = compute_obj_shareH(X1,X2,C,P,W1,W2,H,mu,lambda,gamma)
%   D1 = diag(sum(A1,1)); D2 = diag(sum(A2,1));
%   L1 = D1 - A1; L2 = D2 - A2;
  sigma = diag(1./sum(P)); S = P*sigma*P'; Deg = diag(sum(S)); L = Deg - S;
  error_1 = norm(X1-W1*H,'fro')^2+norm(X2-W2*H,'fro')^2;
  error_2 = mu*sum(sum(C.*P));
  error_3 = 1/lambda * sum(sum((P.*log(P))));
  error_4 = gamma*(trace(H*L*H'));
  
  obj = error_1 + error_2 - error_3 + error_4;
end
