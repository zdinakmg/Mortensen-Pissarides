function [y,P,s]=rouwenhorst(rho,sigma,mu,n)
%
% Rowenhurst method to approximate univariate process by Markov chain
%      y(t) = rho*y(t-1)+ e(t)
%
% Where y(t) is AR(1) process
%       e(t) is iid random variable with mean 0 and std.dev. sigma
%
% INPUTS: rho    - autocorrelation coefficient
%         sigma  - standard deviation of e(t)
%         mu     - mean of y(t)
%         n      - number of states in Markov chain
%
% OUTPUT: y - 1 x n grid  points of Markov chain, evenly spaced states
%         P - n x n transition matrix 
%         s - n x 1 normalized pdf
%
% Reference:
% "Finite State Markov-Chain Approximations to Highly Persistent Processes."
%         by Karen A. Kopecky and Richard M. H. Suen.
%         Review of Economic Dynamics 13 (2010), pp. 701-714.
%         URL: http://www.karenkopecky.net/RouwenhorstPaperFinal.pdf
%
%  Author: Iskander Karibzhanov
%          Department of Economics
%          University of Minnesota
%          karib003@umn.edu
%          http://karibzhanov.com/
%

ybar=sqrt((n-1)/(1-rho^2))*sigma;   
y=linspace(mu-ybar,mu+ybar,n);      % Grid size, evenly space around mean
p=(1+rho)/2; q=p;                   % Defined to hit autocorrelation      
P=rh(n);                            % Call rh function n times
s=zeros(n,1);
for j=1:n
    s(j)=nchoosek(n-1,j-1);         % Binomial coefficient
end
s=s/2^(n-1);                        % Normalized pdf

    function P=rh(h)
        if h==2
            P=[p 1-p; 1-q q];       % 2 x 2 transition matrix
        else
            P1=rh(h-1);             % Call rh h-1 times 
            z=zeros(1,h);
            z1=zeros(h-1,1);
            P=[p*P1 z1; z]+[z1 (1-p)*P1; z]+... 
              [z; (1-q)*P1 z1]+[z; z1 q*P1];    % h x h transition matrix
            P(2:h-1,:)=P(2:h-1,:)/2;        % Normalize stochastic matrix
        end
    end
end