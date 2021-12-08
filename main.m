%% Rosenbrock
disp('Rosenbrock')
nonlinearmin(@rosenbrock,[200;200],'DFP',1e-6,1)

%% Penalty subroutine DFP
disp('Min using DFP on exp(x) function with constraint and penalty method')
F = @(x) exp(x(1)*x(2)*x(3)*x(4)*x(5));
h_1 = @(x) x(1).^2+x(2).^2+x(3).^2+x(4).^2+x(5).^2-10;
h_2 = @(x) (x(2)*x(3))-(5*x(4)*x(5));
h_3 = @(x) x(1).^3+x(3).^3+1;
h_k = {h_1;h_2;h_3};
g_k = {@(x) 0};
mu=[0.01;0.1;1;10;100];
start = [0;2;2;-1;-1];

for i=1:numel(mu)
    F_aux = penalty(F, g_k, h_k, mu(i));
    start = nonlinearmin(F_aux,start, 'DFP',1e-6,1);
end

%% Penalty subroutine BFGS
disp('Min using BFGS on exp(x) function with constraint and penalty method')
F = @(x) exp(x(1)*x(2)*x(3)*x(4)*x(5));
h_1 = @(x) x(1).^2+x(2).^2+x(3).^2+x(4).^2+x(5).^2-10;
h_2 = @(x) (x(2)*x(3))-(5*x(4)*x(5));
h_3 = @(x) x(1).^3+x(3).^3+1;
h_k = {h_1;h_2;h_3};
g_k = {@(x) 0};
mu=[0.01;0.1;1;10;100];
start = [0;2;2;-1;-1];

for i=1:numel(mu)
    F_aux = penalty(F, g_k, h_k, mu(i));
    start = nonlinearmin(F_aux,start, 'BFGS',1e-6,1);
end

%% Exercise 9.3
disp('Exercise 9.3 min using penalty')
F = @(x) exp(x(1))+x(1).^2+x(1)*x(2);
h_k = {@(x) 0.5*x(1)+x(2)-1};
g_k = {@(x) 0};
mu = 4;
start=[0;0];




disp('DFP method')


for i=1:5
    disp('Mu = ');
    disp(mu);
    F_aux = penalty(F,g_k,h_k,mu);
    start=nonlinearmin(F_aux,start,'DFP',1e-6,1);
    mu = mu*10;
end

start = [0;0];
mu = 4;
disp('BFGS method')
for i=1:5
    disp('Mu = ');
    disp(mu);
    F_aux = penalty(F,g_k,h_k,mu);
    start=nonlinearmin(F_aux,start,'BFGS',1e-6,1);
    mu = mu*10;
end



%% Exercise 9.5
disp('Exercise 9.5 min using barrier')
F = @(x) (x(1)-5).^2+(x(2)-3).^2;
g_1 = @(x) x(1)+x(2)-3;
g_2 = @(x) -x(1)+2*x(2)-4;
g_k = {g_1;g_2};
epsilon = 1;
start=[0;0];



disp('DFP method')
for i=1:11
    disp('Epsilon = ');
    disp(epsilon);
    F_aux = barrier(F,g_k,epsilon);
    start=nonlinearmin(F_aux,start,'DFP',1e-6,1);
    epsilon = epsilon/10;
end

disp('BFGS method')
epsilon = 1;
start = [0;0];
for i=1:11
    disp('Epsilon = ');
    disp(epsilon);
    F_aux = barrier(F,g_k,epsilon);
    start=nonlinearmin(F_aux,start,'BFGS',1e-6,1);
    epsilon = epsilon/10;
end


