function [ ZB,P,C,S,Y2,i,flag,Obj,errorRe,Accden] = solution_KASNS( B,XB,YB,dz,gamma,lambda,alpha,maxiter )
%SOLVE20210530 Summary of this function goes here
%   Detailed explanation goes here

% min gamma*||BS-P||^2_F+lambda*||S^T||_{2,1}+alpha*||P*ZB-XB||^2_F
% st.t.,C=S,P^T*P=I
% L = gamma*||BC-P||^2_F+lambda*||S^T||_{2,1}+alpha*||P*ZB-XB||^2_F+0.5*mu*||C-S+Y2*(1/mu)||^2_F
% st.t., P^T*P=I
% Input:  B(d*d) XB(d*n) dz(1*1) lambda alpha 
% Output:  S(m*dz) C(m*dz) P(d*dz) ZB(dz*n)  Y2(m*dz)

flag = 0;
[d,n] = size(XB);
m = size(B,2);

mu = 1e-3;
mu_max = 1e6;
rho = 1.9;


ZB = zeros(dz,n);
Y2 = zeros(m,dz);
P = zeros(d,dz);
S = zeros(m,dz);
C = S;
oldobj = 0;

for i = 1:maxiter
    
    S = (solution_L21((C+Y2*(1/mu))',lambda/mu))';
    for sr = 1:size(S,1)
       if (norm(S(sr,:))<0.00001)
          S(sr,:) = zeros(size(S(sr,:))); 
       end
    end
    
    temp1 = S-Y2*(1/mu);
    C = ((2*gamma*(B'*B))+mu*eye(size(B'*B)))\(B'*P+mu*temp1);
    
     temp2 = 2*gamma*B*C+2*alpha*XB*ZB';
     [U,~,V] = svd(temp2);
     P = (V*eye(dz,d)*U')';
    
    ZB_old = ZB;
    ZB = ((P'*P))\(P'*XB);
    
    Y2 = Y2 + mu*(C-S);
    mu = min(rho*mu,mu_max);

    Snorm = 0;
    for j = 1:size(S,1)
       Snorm =  Snorm + norm(S(j,:));
    end
    obj = gamma*(norm(B*S-P))^2+lambda*Snorm+alpha*((norm(P*ZB-XB))^2);
    errRe = norm(P*ZB-XB);
    errZB = norm(ZB-ZB_old,'inf');
    errCS = norm(C-S,'inf');
    errobj = abs(obj-oldobj);
    Accden = density_cluster(ZB',YB);
    fprintf('iteration=%d, obj=%f, errRe=%f, errZB=%f, errCS=%f, ACCden=%f\n',...
        i,obj,errRe,errZB,errCS,Accden);
    Obj(i) = obj;
    errorRe(i) = errRe;
    
    if(i<=maxiter&&i>1)
        if(errCS<1e-4)
            fprintf('Converge!!!\n');
            flag = 1;
            break;
        end
    end
    
    oldobj = obj;

end

end

