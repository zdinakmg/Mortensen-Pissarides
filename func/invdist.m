function[Q] = invdist(P)
% 
%   Computes the invariant, limiting, or sationary distribution given a 
%           transition matrix for a discrete time Markov chain.
%
% INPUTS:   P - n x n   Transition matrix, conditional prob given state.
%           
% OUTPUT:   Q - n x 1   Invariant row matrix, unconditional prob inde
%
% Reference: Yigit Saglam - Intro to Probab Theory
%            Chapter 5 - Stochastic Processes
%            http://home.uchicago.edu/hickmanbr/uploads/chapter5_1.pdf
%
%            Wiki "Time-homogenous Markov chain w/ finite state space."
%            If MC is irreducible & aperiodic, then unique stationary 
%                     distribution Q = lim as k->inf of P^k exists.
%
%  Author:  Michael Zdinak
%           Department of Economics
%           Washington University in St. Louis
%           zdinakmg@wustl.edu
%

    n = size(P,2);           %   QP = Q          By definition, P is n x n
    B = (P-eye(n));          %   Q(P-I) = O      Subtract Q both sides
                             %   B               n+1 equations, n unknowns
    B(:,n) = ones(n,1);      %                   Replace column n w/ 1                            
    O = zeros(1,n);          %   O               Vector of zeros
    O(n) = 1;                %                   Replace column n w/ 1
    Q = O*B^(-1);            %   Q = O(P-I)^(-1) Multiply O by inverse B
                            
end