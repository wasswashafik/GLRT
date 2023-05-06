function InvCov= func_pseudoInv(imgVec,OptPara)
% Authored: Wasswa Shafik
% Email: wasswashafik@gmail.com
% Date: 2022-02-23
%% Function Usage
% Input:
%     imgVec -- The input matrix
%     OptPara --  Optinal parameters,the default value is 1
% Output:
%     InvCov -- The preudo inverse matrix of imgVec 
%% Main Function

if nargin==1
    OptPara = 1;
elseif OptPara~=fix(OptPara)
    error('MATLAB:func_pseudoInv:WrongNValue',...
        'The value of opt can only be a strictly positive integer 1 or 2!')
end
    

if OptPara == 1
    %% Method 1: Principal component analysis method to find pseudo-inverse (same effect as method 2)
    % This is equivalent to PCA transformation, taking the first PCs eigenvalues and eigenvectors
    [eig_XL,eig_Z]=eig(imgVec);
    [Deig_Z,ind]=sort(diag(eig_Z),'descend');
    D_eigXL=eig_XL(:,ind');

    %% Automatically determine the number of principal components to choose
    rate = 0.9999;% This parameter is adjustable  
    Sumva1 = rate * sum(Deig_Z); % Select the eigenvalues according to the ratio of the sum to 0.99999
    T0=cumsum(Deig_Z);           % Cumsum is an accumulative function, accumulating downward  
    ki=find(T0>Sumva1);   
    PCs=ki(1);

    InvCov=D_eigXL(:,1:PCs)*inv(diag(Deig_Z(1:PCs)))*D_eigXL(:,1:PCs)';

  
elseif OptPara == 2
    %% Method 2: Singular value decomposition to find pseudo-inverse (verified, feasible)
    % Here is equivalent to PCA transformation
    [U_tmp,S_tmp,V_tmp]=svd(imgVec);    
    % PCs=size(S_tmp,2);              % S_tmp (As a square matrix, only meaningful principal component scores are displayed)  
    rate = 0.9999;  % This parameter is adjustable
    Sumva1 = rate * sum(diag(S_tmp)); % Select the eigenvalues according to the ratio of the sum to 0.99999
    T0=cumsum(diag(S_tmp));           % cumsum (an accumulative function, accumulating downward)
    ki=find(T0>Sumva1);
    
    PCs=ki(1);
    inv_S_tmp=S_tmp(1:PCs,1:PCs);        % num_Pcs x  num_Pcs                
    inv_S_tmp=diag(ones(PCs,1)./diag(inv_S_tmp));   % num_Pcs x  num_Pcs 
    InvCov=V_tmp(:,1:PCs)*(inv_S_tmp)*U_tmp(:,1:PCs)'; % num_Pcs x  num_Pcs
end

end

