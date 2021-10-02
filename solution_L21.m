function [S] = solution_L21(A,lambda)
%SOLUTION_L21 此处显示有关此函数的摘要
%   min_{S} 0.5 * || S - A ||^2_F + lambda * || S ||_{2,1}
% || S ||_{2,1} = \sum_{i=1}|| S_{:,i}||_2
% S(n*m) A(n*m)

[n,m] = size(A);
S = zeros(n,m);

for i = 1:m
   normA = norm(A(:,i));
   if(normA > lambda)
       Si = normA-lambda;
       Si = Si/normA;
       Si = Si * A(:,i);
   else
       Si = zeros(n,1);
   end
   S(:,i) = Si;
end


end

