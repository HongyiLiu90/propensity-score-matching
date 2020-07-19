%   Propensity to score matching
%
%   Hongyi Liu
%   email:hongyi.liu@wustl.edu 
%   02/19/17
clear;
%% Data input

% intial values

h = 0.05; % bandwith
nMC = 200; % number of Monte Carlo replications
nobs = 2000; % number of simulated obs
N = nMC*nobs;
matchATT = zeros(nMC,1);
matchATNT = zeros(nMC,1);

% type: matrix
% rows: N
% columns: 10
% mcrep indid theta z x y0 y1 e s d

A =csvread('/Users/MCdata_corr.csv', 1, 0);

    
    mcrep = A(:, 1);
    indid = A(:, 2); %this is individual id
    theta = A(:, 3);
    z = A(:, 4);
    x = A(:, 5);
    y0 = A(:, 6);
    y1 = A(:, 7);
    e = A(:, 8); % e stands for effort
    s = A(:, 9); % s stands for test score
    d = A(:, 10); % stands for schooling decision
    true_ATE = sum(log(y1) - log(y0))/N;  % true param.
    true_ATT = sum((log(y1) - log(y0)).* d(:,1))/sum(d(:,1)); % true param.
    true_ATNT = sum((log(y1) - log(y0)).* (1 - d(:,1)))/sum(1 - d(:,1)); 
    % true param.
    
%% Propensity Score Estimation

 conditioning = [x];
% conditioning = [x z s];
% conditioning = [x theta];
% conditioning = [z s theta];
% conditioning = [x z s theta];
% conditioning = [x s theta];

[b,dev,stats] = glmfit(conditioning,d,'binomial','link','probit');
ps = glmval(b, conditioning ,'probit');

for i = 1:nMC,
    
    lb = (i-1)*nobs+1;
    up = i*nobs;
    p = ps(lb:up,1);
    y = log(y1(lb:up,1)).*d(lb:up,1) + log(y0(lb:up,1)).*(1-d(lb:up,1));
        
    % Epanechnikov Kernel Matching    
    group = d(lb:up,1);
    [ey0 ey1]= epan_kw(h,y,group,p);
    matchATT(i,1) = mean(ey1)- mean(ey0);

    group = 1 - d(lb:up,1);
    [ey1 ey0]= epan_kw(h,y,group,p);
    matchATNT(i,1) = mean(ey1)- mean(ey0);

end;
%% Output

matchATT = mean(matchATT)
matchATNT = mean(matchATNT)

bias_ATT=(mean(matchATT)-true_ATT)/true_ATT
bias_ATNT=(mean(matchATNT)-true_ATNT)/true_ATNT
