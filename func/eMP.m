function [w,u,d] = eMP(kappa,~,A,alpha,~,mu,z,b,delta,theta)
%   
%   Solve for endogenous wage, unemployment rate, & duration of unemploy-
%   ment in a Mortensen-Pissarides model with aggregate fluctuations.
%    
% INPUTS:   kappa  - Posting cost is about 3 days of output
%           beta   - Weekly discount rate
%           A      - Matching effciency
%           alpha  - Elasticity of matching
%           PI     - Transition matrix of discreetized AR process
%           mu     - Bargaining weight for workers
%           z*     - Simulated Productivity shocks
%           b      - Unemployment utility
%           delta  - Separation rate
%           theta* - theta corresponding to z* shock 
%
% OUTPUT:   w  - worker's wages
%           u  - unemployment rate
%           d  - duration of unemployment
%
%  Authors:  Michael Zdinak
%            Department of Economics
%            Washington University in St. Louis
%            zdinakmg@wustl.edu
%
%  Author:   Juan Ignacio Vizcaino (the Nacho)
%            Department of Economics
%            Washington University in St. Louis
%            jivizcaino@go.wustl.edu
%            https://github.com/jivizcaino/Numerical-Methods/
%

    w = mu.*z' + (1-mu)*b + mu.*kappa.*theta;
    
        q = A.*(theta.^-alpha);
        p = q.*theta;
    
    u = delta./(delta+p);
    d = 1./q;
 
end