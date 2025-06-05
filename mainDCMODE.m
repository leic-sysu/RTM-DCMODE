clear;clc;

% Use this code to execute DCMODE
% Input:
%   img:              3-D matrix of the image data after atomospheric correction (row×col×L)
%   endmembers：      Endmember spectra (Sand, algae, and coral) 
%   endmembers_band:  Wavelength of endmember spectra
%   theta_w:          Solar zenith angle
%   theta_v:          View zenith angle
%   a_w:              Absorption coefficient of pure water from 420 nm to 720 nm with an interval of 10 nm
%   b_w:              Backscattering coefficients of pure water from 420 nm to 720 nm with an interval of 10 nm
%   a0:               Empirical parameter for the calculation of chlorophyll-a absorption coefficient
%   a1:               Empirical parameter for the calculation of chlorophyll-a absorption coefficient
%   SRF:              Spectral response function
%   bound:            3-D classification map of optically shallow water (1 = shallow water, 0 = others)

%Output:
%   H_esti3d:         3-D Bathymetry map

%% Read image data
[img,R]=readgeoraster('testdata\ZhaoshuIsland.tif');
info = geotiffinfo('testdata\ZhaoshuIsland.tif');
img = double(img)/pi;
[row,col,L] = size(img);
N = row*col;
img2d = reshape(img,N,L)';

%% Read parameters
% Read endmember spectra
endmember_data = xlsread("testdata\Parameters.xlsx",'ρ','A2:D129');
endmembers = endmember_data(:,2:4);
endmembers_band = endmember_data(:,1);
range1 = round(endmembers_band);  % Original wavelength range of endmember spectra
range2 = 420:10:720;              % Resample to 420 nm ~ 720 nm
range2 = range2';
for i = 1:size(endmembers,2)
    endmembers_intercp(:,i) = interp1(range1,endmembers(:,i),range2,'linear');  % Resample to 420 nm ~ 720 nm with an interval of 10 nm
end

% Read geometric parameters
theta_w = 53.594987/180*pi;  % Solar zenith angle
theta_v = 0;                 % View zenith angle

% Read absorption coefficient of pure water
a_w = xlsread("testdata\Parameters.xlsx",'aw','B18:B138');  % Pope et al. (1997)
range = 1:4:121;
a_w = a_w(range);

% Read backscattering coefficients of pure water
b_w = xlsread("testdata\Parameters.xlsx",'bw','B24:B54')/2;

% Read a0 and a1 (Lee et al., 1998)
a0 = xlsread("testdata\Parameters.xlsx",'a0a1','B5:B35');
a1 = xlsread("testdata\Parameters.xlsx",'a0a1','C5:C35');

% Read the spectral response function (SRF) of Landsat-8 OLI image
SRF_data = xlsread("testdata\L8SRF.xlsx",1,'A2:F302');
range = 1:10:301;
SRF = SRF_data(range,1:4);

% Read the bound of optically shallow water area
bound = double(imread("testdata\bound_ZhaoshuIsland.tif"));
bound_2d = reshape(bound,1,[]);

%% Initialization
X = 0.06;
P = 0.03;
G = 0.01;
B1 = 0.2;
B2 = 0.2;
B3 = 0.2;
H = 3;

%% RTM_DCMODE
% Image masking
NDWI = (img2d(3,:)-img2d(5,:))./(img2d(3,:)+img2d(5,:)); % Landsat8
mask = logical((NDWI>0).*(bound_2d==1));

% Optimization
H_esti = nan(1,N);   % Save the bathymetry estimation result
for i = 1:N
    R_rs_true = img2d(1:4,i);
    if mask(i) == 0
        continue;
    end
    [x,f_x] = DCMODE(R_rs_true,theta_w,theta_v,a_w,b_w,a0,a1,endmembers_intercp,SRF);
    
    % Select the optimal solution
    f_max = max(f_x);
    f_min = min(f_x);
    f_norm = (f_x-repmat(f_min,size(f_x,1),1))./repmat(f_max-f_min,size(f_x,1),1);
    f_sum = sum(f_norm,2);
    [~,I] = sort(f_sum);
    x = x(I(1),:);
    H_esti(i) = x(4);
    fprintf('End of the %dth run\n',i);
end
H_esti3d = reshape(H_esti',row,col);

% Output bathymetry map
geotiffwrite('H_ZhaoshuIsland',H_esti3d,R,'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);

