function  [fi,s2,C] = yuwaest(y,p);
% function  [fi,s2,C] = yuwaest(y,p);
% estimates AR(p) model parameter by solving the
% Yule-Walker equations
% y time series
% fi   autoregressive model parameter  
% s2   estimated WN varians 
% C estimated covariance matrice of the estimated fi
% see Brockwell p 130
% http://www.math.kth.se/matstat/gru/5b1545/homeworks_2004.html

g=acvf(y(:));
g=g(1:p+1); 
g=g(:); %column
p=length(g)-1;

if (p>195),     % use higher overhead, but asymptotically faster algorithm
	%       Levinson-Durbin Algorithm
	fi=zeros(p,1);
	fi=g(1);     % Inititialization
	for K=1:p-1,
		fi(K+1)=(g(K+1)-fi(1:K)*g(K:-1:1))/(1-fi(1:K)*g(1:K));
		fi(1:K)=PHI(1:K)-fi(K+1)*fi(K:-1:1); 
	end;
else   % use the good old \ command
	fi = toeplitz(g(1:p))\g(2:p+1);
end;


s2=g(:)'*[1; -fi];

fi=fi';

C=s2*inv(toeplitz(g(1:p)))/length(y);
