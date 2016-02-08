function [fi, s2, C]=burg(x,p);
% [fi, s2, C]=burg(x,p), estimation of AR(p) parameters 
% using Burg's algorithm, Brockwell p 146
% x time series, p AR order
% fi estimated parameter vector, 
% s2 estimated variance
% C estimetad covariance matrice of the estimated fi1:s
% http://www.math.kth.se/matstat/gru/5b1545/homeworks_2004.html

x=x(:);							% column
n=length(x);
d=x'*x-0.5*x(1)^2-0.5*x(n)^2;
u=x(n:-1:1); 
v=u;

for i=1:p
        fi1(i)=v(2:n-i+1)'*u(1:n-i)/d;
        u1=u;
        u=u(1:n-i)-fi1(i)*v(2:n-i+1);
        v=v(2:n-i+1)-fi1(i)*u1(1:n-i);
        s2=(1-fi1(i)^2)*d/((n-i));
        d=(1-fi1(i)^2)*d-0.5*(v(1)^2+u(n-i)^2);
end

fi=fi1(1);
 
for k=2:p
        fi=fi-fi1(k)*fi(length(fi):-1:1);
        fi=[fi fi1(k)];
end
     
g=acvf(x);
g=g(1:p);
C=s2*inv(toeplitz(g))/n;

