function [x,f_x] = DCMODE(R_rs_true,theta_w,theta_v,a_w,b_w,a0,a1,endmembers,SRF)
% Decomposition-based constrained multi-objective differential evolution (DCMODE) algorithm for shallow water bathymetry retrieval
% Input:
%   R_rs_true:  Observed spectrum (band×1)
%   theta_w:    Solar zenith angle
%   theta_v:    View zenith angle
%   a_w:        Absorption coefficient of pure water from 420 nm to 720 nm with an interval of 10 nm
%   b_w:        Backscattering coefficients of pure water from 420 nm to 720 nm with an interval of 10 nm
%   a0:         Empirical parameter for the calculation of chlorophyll-a absorption coefficient
%   a1:         Empirical parameter for the calculation of chlorophyll-a absorption coefficient
%   endmembers: Endmember spectra of sand, sea grass and coral (band×3)
%   SRF:        Spectral response function
% Output:
%   x:          Final solutions
%   f_x:        Objective function values of the final solutions

%% Initialization
pop_num = 20;
maxiter = 300;
lb = [0.000775,0.005567,0.003926,0,0.002,0.002,0.002];
ub = [0.1,0.03,0.15,40,1.2,1.2,1.2];
[W,pop_num] = UniformPoint(pop_num,2);
x = lhsdesign(pop_num,7);  % LHS sampling
x = x.*repmat((ub-lb)+lb,pop_num,1);
f_x = ones(pop_num,2);

x0 = [0.06 0.03 0.01 3 0.2 0.2 0.2];
options = optimoptions('fmincon','Algorithm','active-set');
options.TolCon = 1e-6;
A = [0 0 0 0 1 1 1;0 0 0 0 -1 -1 -1];
b = [1.2;-0.1];

stop = 0;
repeat = 0;
while ~stop
    [x1, fval, exitflag, output] = fmincon(@costFun_RMSE,x0,A,b,[],[],lb,ub,[],options,R_rs_true,theta_w,theta_v,a_w,b_w,a0,a1,endmembers,SRF);
    if x1(4)>=39
        x0(4) = x0(4)+4;
        repeat = repeat+1;
    end
    if x1(4)<39 || repeat>=3
        stop = 1;
    end
end

stop = 0;
repeat = 0;
x0 = [0.06 0.03 0.01 3 0.2 0.2 0.2];
while ~stop
    [x2, fval, exitflag, output] = fmincon(@costFun_SAD,x0,A,b,[],[],lb,ub,[],options,R_rs_true,theta_w,theta_v,a_w,b_w,a0,a1,endmembers,SRF);
    if x2(4)>=39
        x0(4) = x0(4)+2;
        repeat = repeat+1;
    end
    if x2(4)<39 || repeat>=3
        stop = 1;
    end
end

x(1,:) = x1;
x(2,:) = x2;

% Calculate objective function
for i = 1:pop_num
    f_x(i,1) = costFun_RMSE(x(i,:),R_rs_true,theta_w,theta_v,a_w,b_w,a0,a1,endmembers,SRF);
    f_x(i,2) = costFun_SAD(x(i,:),R_rs_true,theta_w,theta_v,a_w,b_w,a0,a1,endmembers,SRF);
end
cons = [max([zeros(pop_num,1) sum(x(:,5:7),2)-1.2],[],2) max([zeros(pop_num,1) -sum(x(:,5:7),2)+0.1],[],2)];
[z,znad]      = deal(min(f_x),max(f_x));
[z_c,znad_c]  = deal(min(cons),max(cons));

%% Optimization
iter = 0;
while iter < maxiter
    % Genarate offspring by DE/current-to-best/1
    MatingPool = randi(pop_num,1,pop_num);
    Offspring  = OperatorDE(x,lb,ub,f_x,cons);

    % Calculate objective function
    pop_num_off = size(Offspring,1);
    objs = zeros(pop_num_off,2);
    for i = 1:pop_num_off
        objs(i,1) = costFun_RMSE(Offspring(i,:),R_rs_true,theta_w,theta_v,a_w,b_w,a0,a1,endmembers,SRF);
        objs(i,2) = costFun_SAD(Offspring(i,:),R_rs_true,theta_w,theta_v,a_w,b_w,a0,a1,endmembers,SRF);
    end
    f_x = [f_x;objs];
    cons = [cons;max([zeros(pop_num_off,1) sum(Offspring(:,5:7),2)-1.2],[],2) max([zeros(pop_num_off,1) -sum(Offspring(:,5:7),2)+0.1],[],2)];
    
    % Environmental selection
    [x,z,znad,z_c,znad_c,f_x,cons] = EnvironmentalSelection([x;Offspring],W,pop_num,z,znad,z_c,znad_c,f_x,cons);
    iter = iter + 1;
end
end