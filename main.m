%% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Mortensen-Pissarides 
%
%   Michael G. Zdinak
%   Washington University in St. Louis
%
%   Problem: Find labor market tightness (Theta) for each shock value.
%
%   1) Set up all of the environment, 
%   2) Define a function that takes values for log THETA 
%                  and pumps out the N-dimensional residual. 
%   3) Newton? Root finder minimize this residual.
%
%   Note: - Very persistent process, so use a lot of grid points (30-50).
%         - Solving for log theta ensures values 0 < theta < inf.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ENVIRONMENT
clc;
root = 'C:\Users\micha\OneDrive\Documents\Classes\Numerical';
cd([root '\HW\HW1']); addpath('func');
clear;

disp('Mortensen-Pissarides');
disp('');

delta = 0.0081;     %      Separation rate
alpha = 0.72;       %      Elasticity of matching
A     = 0.158;      %      Matching effciency
                    %      Matching function
rho   = 0.9895;     %      Autocorrelation of weekly productivity
sigma = 0.0034;     %      Standard Deviation for innovations
mu    = 0.72;       %      Bargaining weight for workers
b     = 0.4;        %      Unemployment utility
kappa = 0.6;        %      Posting cost is about 3 days of output
beta  = 0.999;      %      Weekly discount rate
n = 50;             %      Grid size, recommend 30-50

%% MAIN
tic
% Discreetize AR Process
[z,P,~] = rouwenhorst(rho,sigma,n);
    
% Make a starting guess at the solution (zero) & call solver 
[theta,fval,exitflag] = MP(kappa,beta,A,alpha,P,mu,z',b,delta,zeros(n,1));

% Simulate a stochastic process for series of aggregate shocks
[z_star,ug] = simulateMC(z,P,n*100);

% Check simulation accuracy by estimating AR(1) process on output
[rho_star, s2, ~] = burg(z_star,1);          % yuwaest(z_star,1); 
sigma_star = sqrt(s2);

% Find corresponding theta_star/theta for each z_star/z
[~,ind] = ismember(z_star',z');
theta_star = theta(ind);

% Realizations of productivity & time-series of the endogenous variables.

toc

%% Output
fprintf('\n\tRho\t\tRho*\tSigma\tSigma*\n');
fprintf('\t---\t\t----\t-----\t-----\n\t');
fprintf('%0.3f\t',[rho rho_star sigma sigma_star]);
fprintf('\n\n');

plot(theta_star);
% plot(z_star);

% Computational limits, repeating grid bootstrap moments
