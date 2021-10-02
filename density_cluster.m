function [acc,nmi] = density_cluster(x,y,N)
    N = 1;
    dist = pdist2(x,x);
    para.method = 'gaussian';
    n = 0;
    for pp = 1:0.1:2
        para.percent = pp;
        [n_class,~] = size(unique(y));
        [cluster_lables, ~] = cluster_dp(dist, para, n_class);
        acc1 = Accuracy(cluster_lables,y);
        n = n+1;
        accs(n) = acc1;
    end
    acc = max(accs);
    nmi = compute_nmi(y,cluster_lables);
end