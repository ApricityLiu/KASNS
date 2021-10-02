function [cluster_lables, center_idxs] = cluster_dp(dist, para, k)

%% Input and Output
% INPUT :
% dist : A nCases*nCases matrix. each dist(i, j) represents the distance
%        between the i-th datapoint and j-th datapoint.
% para : options
%        percent : The parameter to determine dc. 1.0 to 2.0 can often yield good performance.
%        method  : alternative ways to compute density. 'gaussian' or
%                  'cut_off'. 'gaussian' often yields better performance.
% OUTPUT :
% cluster_labels : a nCases vector contains the cluster labls. Lable equals to 0 means it's in the halo region
% center_idxs    : a nCluster vector contains the center indexes.

%% Estimate dc
% disp('Estimating dc...');
percent = para.percent;
N = size(dist,1);
position = round(N*(N-1)*percent/100);
tri_u = triu(dist,1);
sda = sort(tri_u(tri_u~=0));
dc = sda(position);
clear sda; clear tri_u;
% dc = 1.5
%% Compute rho(density)
% fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);
switch para.method
    % Gaussian kernel
    case 'gaussian'
        rho = sum(exp(-(dist./dc).^2),2)-1;
        % "Cut off" kernel
    case 'cut_off'
        rho = sum((dist-dc)<0, 2);
end
[~,ordrho]=sort(rho,'descend');

%% Compute delta
% disp('Computing delta...');
delta = zeros(size(rho));
nneigh = zeros(size(rho));

delta(ordrho(1)) = -1;
nneigh(ordrho(1)) = 0;
for i = 2:size(dist,1)
    range = ordrho(1:i-1);
    [delta(ordrho(i)), tmp_idx] = min(dist(ordrho(i),range));
    nneigh(ordrho(i)) = range(tmp_idx); 
end
delta(ordrho(1)) = max(delta(:));

%% Decision graph, choose min rho and min delta
rho = (rho-min(rho))/(max(rho)-min(rho));
delta = (delta-min(delta))/(max(delta)-min(delta));
% rs = (rho-1).^2+(delta-1).^2;
rs = rho.*delta;
[idxs,~] = argsort(rs,'descend');
center_idxs = idxs(1:k);
% disp([num2str(length(center_idxs)),' cluster centers found...']);
% figure(10000);
% plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
% hold on;
% plot(rho(center_idxs),delta(center_idxs),'o','MarkerSize',5,'MarkerFaceColor','r','MarkerEdgeColor','r');
% title ('Decision Graph','FontSize',15.0)
% xlabel ('\rho')
% ylabel ('\delta')
% hold off;
%% Assignment
% raw assignment
% disp('Assigning data-points into clusters...');
cluster_lables = -1*ones(size(dist,1),1);
for i = 1:length(center_idxs)
    cluster_lables(center_idxs(i)) = i;
end
for i=1:length(cluster_lables)
    if (cluster_lables(ordrho(i))==-1)
        cluster_lables(ordrho(i)) = cluster_lables(nneigh(ordrho(i)));
    end
end
raw_cluster_lables = cluster_lables;

% % find and cut off halo region
% disp('Cut off halo regions...');
% for i = 1:length(center_idxs)
%     tmp_idx = find(raw_cluster_lables==i);
%     tmp_dist = dist(tmp_idx,:);
%     tmp_dist(:,tmp_idx) = max(dist(:));
%     tmp_rho = rho(tmp_idx);
%     tmp_lables = raw_cluster_lables(tmp_idx);
%     tmp_border = find(sum(tmp_dist<dc,2)>0);
%     if ~isempty(tmp_border)
%         rho_b = max(tmp_rho(tmp_border));
%         halo_idx = rho(tmp_idx) < rho_b;
%         tmp_lables(halo_idx) = 0;
%         % lable equals to 0 means it's in the halo region
%         cluster_lables(tmp_idx) = tmp_lables;
%     end
% end

end

function [sorted_indices, sorted_M] = argsort(M, mode)   
    if (~exist('mode', 'var'))
        mode = 'ascend';
    end    
    
    nRows = size(M,1);
    indices = (1 : nRows)';
    M_with_indices = horzcat(M, indices);
    [sorted_M, sorted_indices] = sort(M_with_indices(:,1), mode);
end