function Dict=ConstructionD_lilu(data,K,P)
%% Function Usage
% Input:
%     data -- (data£ºbands x pixels£©
%        K --  K is the number of clusters
%        P --  P number of atoms reserved for each class
% K=15,P=20, The setting of these values refers to the LRASR algorithm
% Output:
%    Dict --dictionary, each column is an atom, data£ºbands x nums
%% Main Function

IDX=kmeans(data',K,'start','plus'); % 'plus'Default initial centroid position
Dict=[];
for i=1:K
    pos=find(IDX==i);
    D_temp=data(:,pos);
    if(size(D_temp,2)<P)
        continue; % continue the function is to end this loop, skip the code after continue, and continue to the next loop operation
    end
    %% The following code is for RX detection
    mu=mean(D_temp,2); % mean of all vectors, a column vector
    COV_inv=pinv(cov(D_temp')); % find the covariance of this class
    D_temp_C=D_temp-repmat(mu,[1,size(D_temp,2)]);
    Dis=zeros(1,size(D_temp,2));
    for j=1:size(D_temp,2)
        Dis(j)=D_temp_C(:,j)'*COV_inv*D_temp_C(:,j); % RX detection algorithm
    end
    %% Sort the Mahalanobis distance from small to large
    [~,Ind]=sort(Dis);  
    Dict=[Dict,D_temp(:,Ind(1:P))]; % Take the P with the smallest distance as the dictionary of this class
end
