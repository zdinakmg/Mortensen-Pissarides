function y=acvf(Z)
% acvf(Z) sample autocovarians; Brockwell page 17
 Z=Z(:)'-mean(Z);
N = length(Z);
y = filter(Z(N:-1:1),1,Z);
y = y(N:-1:1)/N;
