%% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Alternative Specification :: Mortensen-Pissarides
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

fprintf('\n\tHagedorn-Manovskii :: Mortensen-Pissarides\n');

delta = 0.0081;     %      Separation rate
alpha = 0.72;       %      Elasticity of matching
A     = 0.158;      %      Matching effciency
                    %      Matching function
rho   = 0.9895;     %      Autocorrelation of weekly productivity
sigma = 0.0034;     %      Standard Deviation for innovations
mu    = 0.05;       %      NEW! Bargaining weight for workers
b     = 0.95;       %      NEW! Unemployment utility
kappa = 0.6;        %      Posting cost is about 3 days of output
beta  = 0.999;      %      Weekly discount rate
n = 50;             %      Grid size, recommend 30-50
t = n*100;          %      Number of simulations

%% MAIN
tic
% Discreetize AR Process
[z,P,~] = rouwenhorst(rho,sigma,1.1,n);
    
% Make a starting guess at the solution (zero) & call solver 
[theta,fval,exitflag] = MP(kappa,beta,A,alpha,P,mu,z',b,delta,ones(n,1));

% Impose eower bound on theta
theta( le(theta,(1/A).^(-1./alpha)) ) = (1/A).^(-1./alpha);

% Simulate a stochastic process for series of aggregate shocks
[z_star,ug] = simulateMC(z,P,t);

% Check simulation accuracy by estimating AR(1) process on output
[rho_star, s2, ~] = burg(z_star,1);          % yuwaest(z_star,1); 
sigma_star = sqrt(s2);

% Find corresponding theta_star/theta for each z_star/z
[~,ind] = ismember(z_star',z');
theta_star = theta(ind);

% Realizations of productivity & time-series of the endogenous variables.
[w_star, u_star, d_star, q_star] = eMP(kappa,beta,A,alpha,P,mu,z_star, ...
                                          b,delta,theta_star);
toc

%% Output

fprintf('\n\t\tRho\t\tRho*\tSigma\tSigma*\n');
fprintf('\t\t---\t\t----\t-----\t-----\n\t');
fprintf('ar\t');
fprintf('%0.3f\t',[rho rho_star sigma sigma_star]);
fprintf('\n');

fprintf('\n\t\tMin\t\tMax\t\tMean\tVar\n');
fprintf('\t\t---\t\t----\t-----\t-----\n\t');
fprintf('z*\t');
fprintf('%0.3f\t',[min(z_star) max(z_star) ...
                  mean(z_star) var(z_star)]);
fprintf('\n\tt*\t');
fprintf('%0.3f\t',[min(theta_star) max(theta_star) ...
                  mean(theta_star) var(theta_star)]);
fprintf('\n\tw*\t');
fprintf('%0.3f\t',[min(w_star) max(w_star) ...
                  mean(w_star) var(w_star)]);
fprintf('\n\tq*\t');
fprintf('%0.3f\t',[min(q_star) max(q_star) ...
                  mean(q_star) var(q_star)]);
fprintf('\n\t\t---\t\t----\t-----\t-----\n\t');
fprintf('u\t');
fprintf('%0.3f\t',[min(u_star) max(u_star) ...
                  mean(u_star) var(u_star)]);
fprintf('\n\td\t');
fprintf('%0.3f\t',[min(d_star) max(d_star) ...
                  mean(d_star) var(d_star)]); 
fprintf('\n\t\t---\t\t----\t-----\t-----\n\n');
 
fig = figure('Name','alt','NumberTitle','off');
subplot(3,2,1)
plot(theta_star);
    hline = refline([0 mean(theta_star)]);
    hline.Color = 'r';
    hline.LineStyle = '--';
    hline.LineWidth = 1;
    title('theta')
    axis tight
    xlim([0 t])

subplot(3,2,2)
plot(q_star);
    hline = refline([0 mean(q_star)]);
    hline.Color = 'r';
    hline.LineStyle = '--';
    hline.LineWidth = 1;
    title('q')
    axis tight
    xlim([0 t])

subplot(3,2,3)
plot(z_star);
    hline = refline([0 mean(z_star)]);
    hline.Color = 'r';
    hline.LineStyle = '--';
    hline.LineWidth = 1;
    title('z')
    axis tight
     xlim([0 t])
    set(gca,'XTickLabel',[]);
       
subplot(3,2,4)
plot(w_star);
    hline = refline([0 mean(w_star)]);
    hline.Color = 'r';
    hline.LineStyle = '--';
    hline.LineWidth = 1;
    title('wage')
    axis tight
    xlim([0 t])
    set(gca,'XTickLabel',[]);
    
subplot(3,2,5)
plot(u_star*100);
    hline = refline([0 mean(u_star*100)]);
    hline.Color = 'r';
    hline.LineStyle = '--';
    hline.LineWidth = 1;
    title('unemployment %')
    axis tight
    xlim([0 t])
    set(gca,'XTickLabel',[]);
    
subplot(3,2,6)
plot(d_star);
    hline = refline([0 mean(d_star)]);
    hline.Color = 'r';
    hline.LineStyle = '--';
    hline.LineWidth = 1;
    title('duration')
    axis tight
    xlim([0 t])
    set(gca,'XTickLabel',[]);
saveas(fig,'.\out\alt.pdf');     
fprintf('\tfin!\n\n');
