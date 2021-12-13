
clear;
clc;
close all;
addpath(genpath(cd));

dbs1 = 'C';
dbs2 = 'V';

name = ['W_',dbs1];
load(name);
[tx,ty] = load_data(dbs2);

dz = 70;
gamma = 0.5;
lambda = 1e-3;
alpha = 1e-4;
maxiter = 100;

[ ZB,P,C,S,Y2,iter,flag,Obj,errorRe,Accden] = solution_KASNS( W,tx',ty,dz,gamma,lambda,alpha,maxiter );
 



