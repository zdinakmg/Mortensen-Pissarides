function [sol,fval,exitflag] = MP(kappa,beta,A,alpha,PI,mu,z,b,delta,x0)
%
%   Mortensen-Pissarides MP model with aggregate fluctuations.
%
%       Minimize the residual of N+1 equations to solve for N+1 values 
%       of theta, here labor market tightness, v/u?
%    
% INPUTS:   kappa - Posting cost is about 3 days of output
%           beta  - Weekly discount rate
%           A     - Matching effciency
%           alpha - Elasticity of matching
%           PI    - Transition matrix of discreetized AR process
%           mu    - Bargaining weight for workers
%           z     - Productivity shocks
%           b     - Unemployment utility
%           delta - Separation rate
%           x0    - theta0 = zeros(N,1)
% 
% OUTPUT:   sol   - theta n x 1
%           fval  -  value of the equations FUN at x
%           exitflag  - 1 if fsolve converged to a root.
%                       2  Change in x too small.
%
% Reference:
%
%  Author: Juan Ignacio Vizcaino (the Nacho)
%          Department of Economics
%          Washington University in St. Louis
%          jivizcaino@go.wustl.edu
%          https://github.com/jivizcaino/Numerical-Methods/
%

    % TO DO ...
    % optimoptions(@fsolve);
    % options.TolX        = 0.001;
    % options.TolFun      = 0.001;
    % options.MaxIter     = 2000;
    % Make these arguments, not sure they are used?
    % Set option to display information after each iteration
    % options=optimset('Display','iter');

    function sol = residual(theta)
        % 
        % Minimize residual 
        % Solving for log theta ensures 0 < theta < inf
        % Matching function
        %       m(u,v) = min[1,A*1/theta^-alpha] 
        %
            sol = (kappa/(beta*A)).*(theta.^alpha)-PI*((1-mu).* ...
                   (z-b)-(kappa*mu).*(theta)+(((1-delta)*kappa)/A).* ... 
                   (theta.^alpha));
             %q = min(1,A.*theta.^alpha);  
             %sol = q.^-alpha/A;
             sol = log(sol);     
    end

    [sol,fval,exitflag] = fsolve (@residual,x0);
    
end
