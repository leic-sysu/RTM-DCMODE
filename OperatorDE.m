function Offspring = OperatorDE(Parent,lb,ub,f_x,cons)
%OperatorDE - Crossover and mutation operators of DE algorithm.
% Input:
%   Parent:     Parent solutions
%   lb:         Lower bound of decision variables
%   ub:         Upper bound of decision variables
%   f_x:        Objective function values of parent solutions
%   b_w:        Backscattering coefficients of pure water from 420 nm to 720 nm with an interval of 10 nm
%   cons:       Normalized CV values of parent solutions
% Output:
%   Offspring:  Offspring solutions

[CR,F,proM,disM] = deal(1,0.5,1,20);
[N,D] = size(Parent);
[FrontNo,MaxFNo] = NDSort(f_x,cons,N);
Gbest = Parent(find(FrontNo == 1),:);
pop_num = size(Parent,1);
MatingPool = randi(pop_num,1,2*pop_num);
MatingPool_Gbest = randi(size(Gbest,1),1,pop_num);
Parent1 = Parent;
Parent2 = Parent(MatingPool(1:end/2),:);
Parent3 = Parent(MatingPool(end/2+1:end),:);
Parent_Gbest = Gbest(MatingPool_Gbest,:);

% %% DE/rand/1
% Site = rand(N,D) < CR;
% Offspring       = Parent1;
% Offspring(Site) = Offspring(Site) + F*(Parent2(Site)-Parent3(Site));
% % Polynomial mutation
% Lower = repmat(lb,N,1);
% Upper = repmat(ub,N,1);
% Site  = rand(N,D) < proM/D;
% mu    = rand(N,D);
% temp  = Site & mu<=0.5;
% Offspring       = min(max(Offspring,Lower),Upper);
% Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
%     (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
% temp = Site & mu>0.5;
% Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
%     (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));

%% DE/rand-to-Gbest/1
Site = rand(N,D) < CR;
Offspring2       = Parent1;
Offspring2(Site) = Offspring2(Site) + rand()*(Parent_Gbest(Site)-Offspring2(Site)) + F*(Parent2(Site)-Parent3(Site));
% Polynomial mutation
Lower = repmat(lb,N,1);
Upper = repmat(ub,N,1);
Site  = rand(N,D) < proM/D;
mu    = rand(N,D);
temp  = Site & mu<=0.5;
Offspring2       = min(max(Offspring2,Lower),Upper);
Offspring2(temp) = Offspring2(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
    (1-(Offspring2(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
temp = Site & mu>0.5;
Offspring2(temp) = Offspring2(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
    (1-(Upper(temp)-Offspring2(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));

Offspring = Offspring2;
end