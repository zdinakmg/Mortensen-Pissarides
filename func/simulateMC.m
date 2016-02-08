function [theta,ug] = simulateMC(y,P,n)
%   
%   Simulates a stochastic process following a discrete time Markov chain.
%   Starting point drawn from stationary distribution of transition matrix. 
%
% INPUTS:   y  - 1 x n grid points of Markov chain, evenly spaced states
%           P  - transition matrix for the Markov chain
%           n  - number of periods to simulate
%
% OUTPUT:   theta   - pseudo-random numbers, Markov chain for simulation
%           ug      - pseudo-random numbers, uniform distribution
%
% Reference: Yigit Saglam - Intro to Probab Theory
%            Chapter 5 - Stochastic Processes
%            http://home.uchicago.edu/hickmanbr/uploads/chapter5_1.pdf
%
%  Author:  Michael Zdinak
%           Department of Economics
%           Washington University in St. Louis
%           zdinakmg@wustl.edu
%
                               % Starting from steady state      (1)
P_bar = invdist(P);            % Unconditional distribution
CDF_bar = cumsum(P_bar,2);     % CDF of P_bar
                               % Index  initial steady state                                
i = find( CDF_bar >= ... 
          rand(1),1,'First');
                               % Starting from arbitrary state   (2)
% i = 1;                       % Index initial arbitrary shock
ind = zeros(1,n);
ind(1) = i;
CDF = cumsum(P,2);             % CDF of P

    for t = 1:n;                
        ug(t) = rand(1);        
        j = find(CDF(i,:)>=ug(t),1,'First');
        ind(t) = j;
        i = j; 
    end
    
    theta = y(ind);             % If using arbitrary starting point
                                % Remember to burn-in simulation
end