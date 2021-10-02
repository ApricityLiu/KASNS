function [W1,W2,W_all] = gen_general(X,Y)

[W1,~]=LDA(X,Y);

W1 = W1*diag(1./sqrt(sum(W1.*W1,1)));

W2 = null(W1');

W2 = W2*diag(1./sqrt(sum(W2.*W2,1)));

W_all = [W1,W2];

end