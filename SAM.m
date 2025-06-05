function [theta]=SAM(A,B)
% Calculate the spectral angle distance (SAD) between spectrum A and spectrum B
AB = A'*B;
Anorm = norm(A,2);
Bnorm = norm(B,2);
V = AB/(Anorm*Bnorm);
if V > 1
    theta = 0;
else
    theta = acos(V);
end
