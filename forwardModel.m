function [rs] = forwardModel(theta_w,theta_v,a,b_b,H,rho)
% The forward model of HOPE proposed by Lee et al.
u = b_b./(a+b_b);
rs_deep = (0.084+0.17*u).*u;  % Lee et al.,1998
rs = rs_deep.*(1-exp(-(1/cos(theta_w)+1.03*sqrt(1+2.4*u)/cos(theta_v)).*(a+b_b)*H))+rho/pi.*exp(-(1/cos(theta_w)+1.04*sqrt(1+5.4*u)/cos(theta_v)).*(a+b_b)*H);
end