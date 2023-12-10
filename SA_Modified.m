% Data of the Multi-Objective Optimization Model
% Number of components and periods
N = 10;
T = 24;
J = 24;
L = T/J;
% Specification of the components
% Parameters of the Failure function
Lambda = [0.00013 0.00086 0.00027 0.00033 0.00057 0.00043 0.00023 0.00024 0.00049 0.00066];
Beta = [2.55 2.02 2.26 2.35 2.39 2.49 2.68 2.64 2.66 2.86];
% Improvement factor (Age reduction coefficient)
Alpha = [0.71 0.67 0.92 0.51 0.71 0.51 0.75 0.93 0.54 0.88];
% Failure cost
Failure_Cost = [350 260 250 310 300 250 280 310 300 270];
% Maintenance cost
M_Cost = [35 64 60 34 65 53 79 36 33 60];
% Replacement cost
R_Cost = [200 240 190 240 280 280 200 220 200 200];
% Fixed cost
Fixed_Cost = 800; 

% Engineering economics parameters
% Inflation rates of failure cost, maintenance cost, replacement cost and fixed cost
Inf_Failure = 0; %0.01/12;
Inf_M = 0; %0.015/12;
Inf_R = 0; %0.02/12;
Inf_Fix = 0; %0.01/12;
% Interest rate
Int_Rate = 0; %0.03/12;

% Parameters of the multi-objective optimization model
% Weights of the objective functions in weighted method, W1+W2 = 1
%W1 = 0.0; W2 = 1.0;
%W1 = 0.1; W2 = 0.9;
%W1 = 0.2; W2 = 0.8;
%W1 = 0.3; W2 = 0.7;
%W1 = 0.4; W2 = 0.6;
%W1 = 0.5; W2 = 0.5;
%W1 = 0.6; W2 = 0.4;
%W1 = 0.7; W2 = 0.3;
W1 = 0.8; W2 = 0.2;
%W1 = 0.9; W2 = 0.1;
%W1 = 1.0; W2 = 0.0;
% Design goals for the objective functions in goal attainment method
% Given budget
GB = 5000;
% Required reliability
RR = 0.50;

%------------------------SA Algorithm--------------------------------------

% Simulated Annealing Algorithm
% Simulated annealing algorithm parameters
% Initial temperature: 1000000
% Final temperature: 0.01
% Decreasing rate: 0.98
t_initial = 1000000;
t_final = 0.01;
t_rate = 0.98;
min = 0;
max = 2;
%Data;
N = 10;
T = 24;
J = 24;
L = T/J;
% Initial solution
a = zeros(1,T*N);
for j = 1:1:T*N
 a(j) = fix((max-min+1)*rand+min);
end
[cost,reliability,fit1,fit2,fit3] = Fitness(a);
initial_solution(1,1:N*T) = a ;
initial_solution(1,N*T+1:N*T+5) = [cost,reliability,fit1,fit2,fit3];
x = initial_solution;


t_current = t_initial;
i = 1;
while t_final <= t_current
    
    % Transition procedure
    y = Transition(x);
    [cost,reliability,fit1,fit2,fit3] = Fitness(y);
    y(1,N*T+1:N*T+5) = [cost,reliability,fit1,fit2,fit3];

    % Acceptation procedure
    if y(1,N*T+2) > 0.95
        if y(1,N*T+1) < x(1,N*T+1)
            x = y;
        elseif y(1,N*T+1) >= x(1,N*T+1)
            if rand <= exp(-(y(1,N*T+1)-x(1,N*T+1))/t_current)
                x = y;
            end
        end
    end
    
    
    solution_improvement(i,1:N*T+5) = x;
    t_current = t_rate*t_current;
    i = i+1;
end


% This section changes the final solution(1,N*T) to PMR_Schedule(N,T)
ss = sortrows(solution_improvement,N*T+1);
final_solution = ss(1:1,:);
PMR_Schedule = zeros(N,T);
for i = 1:1:N
    for j = 1:1:T
    PMR_Schedule(i,j) = final_solution(1,(i-1)*T+j);
    end
end



%---------------Functions Used in SA algorithm-----------------------------
% Fitness Functions of the Cost and Reliability Functions
function [cost,reliability,fit1,fit2,fit3] = Fitness(a)
%Data;
% This section changes a(1,N*T) to A(N,T)
N = 10;
T = 24;
J = 24;
L = T/J;
% Specification of the components
% Parameters of the Failure function
Lambda = [0.00013 0.00086 0.00027 0.00033 0.00057 0.00043 0.00023 0.00024 0.00049 0.00066];
Beta = [2.55 2.02 2.26 2.35 2.39 2.49 2.68 2.64 2.66 2.86];
% Improvement factor (Age reduction coefficient)
Alpha = [0.71 0.67 0.92 0.51 0.71 0.51 0.75 0.93 0.54 0.88];
% Failure cost
Failure_Cost = [350 260 250 310 300 250 280 310 300 270];
% Maintenance cost
M_Cost = [35 64 60 34 65 53 79 36 33 60];
% Replacement cost
R_Cost = [200 240 190 240 280 280 200 220 200 200];
% Fixed cost
Fixed_Cost = 800; 

% Engineering economics parameters
% Inflation rates of failure cost, maintenance cost, replacement cost and fixed cost
Inf_Failure = 0; %0.01/12;
Inf_M = 0; %0.015/12;
Inf_R = 0; %0.02/12;
Inf_Fix = 0; %0.01/12;
% Interest rate
Int_Rate = 0; %0.03/12;
% Weights of the objective functions in weighted method, W1+W2 = 1
%W1 = 0.0; W2 = 1.0;
%W1 = 0.1; W2 = 0.9;
%W1 = 0.2; W2 = 0.8;
%W1 = 0.3; W2 = 0.7;
%W1 = 0.4; W2 = 0.6;
%W1 = 0.5; W2 = 0.5;
%W1 = 0.6; W2 = 0.4;
%W1 = 0.7; W2 = 0.3;
W1 = 0.8; W2 = 0.2;
%W1 = 0.9; W2 = 0.1;
%W1 = 1.0; W2 = 0.0;
% Design goals for the objective functions in goal attainment method
% Given budget
GB = 5000;
% Required reliability
RR = 0.90;


for i = 1:1:N
    for j = 1:1:T
    A(i,j) = a(1,(i-1)*T+j);
    end
end

% This section calculates the x(starting effective age) and xp(ending effective age)
x = zeros(N,T);
for i = 1:1:N
    for j = 1:1:T-1
        if A(i,j) == 0
            x(i,j+1) = x(i,j)+L;
        elseif A(i,j) == 1
            x(i,j+1) = Alpha(i)*(x(i,j)+L);
        elseif A(i,j) == 2
            x(i,j+1) = 0;
        end
    end
end
xp = x+L;
% This section calculates the cost and reliability functions for series system of components
cost = 0;
max_cost = 0;
xx = zeros(N,T);
xxp = xx+L;
reliability = 1;
reliability_by_comp = zeros(N,T);
for j = 1:1:T
    counter = 0;
    for i = 1:1:N
        if A(i,j) == 0
            cost = cost+((Failure_Cost(i)*Lambda(i)*((xp(i,j)^Beta(i))-(x(i,j)^Beta(i)))*(1+Inf_Failure)^j));
        elseif A(i,j) == 1
            cost = cost+((Failure_Cost(i)*Lambda(i)*((xp(i,j)^Beta(i))- (x(i,j)^Beta(i)))*(1+Inf_Failure)^j)+(M_Cost(i)*(1+Inf_M)^j));
        elseif A(i,j) == 2
        cost = cost+((Failure_Cost(i)*Lambda(i)*((xp(i,j)^Beta(i))-(x(i,j)^Beta(i)))*(1+Inf_Failure)^j)+(R_Cost(i)*(1+Inf_R)^j));
        end
        if A(i,j) == 1 || A(i,j) == 2
            counter = 1;
        end
        max_cost = max_cost+((Failure_Cost(i)*Lambda(i)*((xxp(i,j)^Beta(i))-(xx(i,j)^Beta(i)))*(1+Inf_Failure)^j)+(R_Cost(i)*(1+Inf_R)^j));
        reliability = reliability*exp(-Lambda(i)*((xp(i,j)^Beta(i))-(x(i,j)^Beta(i))));
        reliability_by_comp(i,j) = exp(-Lambda(i)*((xp(i,j)^Beta(i))-(x(i,j)^Beta(i))));
    end
 
    if counter == 1
        cost = cost+(Fixed_Cost*(1+Inf_Fix)^j);
    end
    cost = cost*(1+Int_Rate)^(-j);
    max_cost = max_cost+(Fixed_Cost*(1+Inf_Fix)^j);
    max_cost = max_cost+(1+Int_Rate)^(-j);
end

reliability = 1;
for j = 1:1:T
    reliability = reliability * reliability_by_comp(1,j) * (1 - ((1 - reliability_by_comp(2,j))*(1-reliability_by_comp(3,j))*(1-reliability_by_comp(4,j)))) * ...
            reliability_by_comp(5,j)*reliability_by_comp(6,j)*(1 - ((1-reliability_by_comp(7,j))*(1-reliability_by_comp(8,j)))) *(1 - ((1-reliability_by_comp(9,j))*(1-reliability_by_comp(10,j))));
end


% The fitness functions,
fit1 = W1*(cost/max_cost)+W2*(-reliability);
fit2 = -reliability+(1/max_cost)*abs(GB-cost);
fit3 = (cost/max_cost)+abs(RR-reliability);
 %end function
end

% Transition Function
function [x] = Transition(x)
%Data;
N = 10;
T = 24;
J = 24;
L = T/J;
transition_point = fix(N*T*rand+1);
if x(:,transition_point) == 0
    if (rand < 0.5)
        for k = 1:1:N
            if mod(transition_point,T) == 0
                x(:,(mod(transition_point,T)+k*T)) = 1;
            else
                x(:,(mod(transition_point,T)+(k-1)*T)) = 1;
            end
        end
    elseif (rand >= 0.5)
        for k = 1:1:N
            if mod(transition_point,T) == 0
                x(:,(mod(transition_point,T)+k*T)) = 2;
            else
                x(:,(mod(transition_point,T)+(k-1)*T)) = 2;
            end
        end
    end
elseif x(:,transition_point) == 1 || x(:,transition_point) == 2
    for k = 1:1:N
        if mod(transition_point,T) == 0
            x(:,(mod(transition_point,T)+k*T)) = 0;
        else
            x(:,(mod(transition_point,T)+(k-1)*T)) = 0;
        end
    end
end
end %end function

