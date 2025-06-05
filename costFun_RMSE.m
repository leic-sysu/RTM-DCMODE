function [obj] = costFun_RMSE(x0,R_rs_true,theta_w,theta_v,a_w,b_w,a0,a1,endmembers,SRF)
% Calculate the RMSE objective function
% Input:
%   x0:         Values of the unknown variables (i.e., X,P,G,H,B1,B2,B3)
%   R_rs_true:  Observed spectrum
%   theta_w:    Solar zenith angle
%   theta_v:    View zenith angle
%   a_w:        Absorption coefficient of pure water from 420 nm to 720 nm with an interval of 10 nm
%   b_w:        Backscattering coefficients of pure water from 420 nm to 720 nm with an interval of 10 nm
%   a0:         Empirical parameter for the calculation of chlorophyll-a absorption coefficient
%   a1:         Empirical parameter for the calculation of chlorophyll-a absorption coefficient
%   endmembers: Endmember spectra of sand, sea grass and coral
%   SRF:        Spectral response function
% Output:
%   obj:        RMSE value between the observed spectrum and the simulated spectrum

% Read X and calculate the total backscattering coefficient
X = x0(1);
lambda = 420:10:720;
lambda = lambda';
b_bp = X*(400./lambda).^0.681;  % Particles backscattering coefficient
b_b = b_w + b_bp;

% Read P, G and calculate the total absorption coefficient
P = x0(2);
G = x0(3);
a_phy = (a0+a1*log(P))*P;           % Chlorophyll-a absorption coefficient
a_g = G*exp(-0.0166*(lambda-440));  % CDOM chlorophyll-a absorption coefficient
a = a_w+a_phy+a_g;

% Read endmember abundance and calculate the bottom reflectance
B1 = x0(5);  % Sand
B2 = x0(6);  % Algae
B3 = x0(7);  % Coral
rho = B1*endmembers(:,1)+B2*endmembers(:,2)+B3*endmembers(:,3);

% Read water depth
H = x0(4);

% Generate the simulated spectrum
rs = forwardModel(theta_w,theta_v,a,b_b,H,rho);  % The forward model proposed by Lee et al.
R_rs = 0.5*rs./(1-1.5*rs);                       % Water-leaving reflectance
R_rs_multi = R_rs'*SRF./sum(SRF);                % Match the spectral range of Landsat-8 OLI image
R_rs_model = R_rs_multi';

obj = sqrt(mean((R_rs_model-R_rs_true).^2))/mean(R_rs_true); % RMSE
end