function [ sample ] = spectral_clustering( W,k )
% Executes spectral clustering algorithm using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}
%   Input: W --- Adjacency matrix;
%          k --- number of clusters to look for
%   Output: groups --- N-dimensional vector containing the memberships of the N points 
%                        to the n groups obtained by spectral clustering

N = size(W,1);
%MaxIter = 1000;
%Replicates = 20;

D = diag(1./sqrt(sum(W)+eps));
L = eye(N)-D*W*D;

[~,~,V] = svd(L);
sample = V(:,N-k+1:N);

for i = 1:N
   sample(i,:) =  sample(i,:)./norm(sample(i,:)+eps);
end

%groups = kmeans(sample,k,'maxiter',MaxIter,'replicates',Replicates,'emptyaction','singleton');

end

