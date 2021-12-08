%% Temp
F = @(x) exp(x(1))+x(1).^2+x(1)*x(2);
h_k = {@(x) 0.5*x(1)+x(2)-1};
g_k = {0};
mu = 4;
start=[0;0];

F_aux = penalty(F,g_k,h_k,mu);
nonlinearmin(F_aux,start,'BFGS',1e-6,1)

%% Linesearch testing
func = @(x) (1e58*x(1)-1)^100+(1e-58*x(2)-1)^2-2;

display('Test 1');
x=[0;0];
d=[0;1];
[lambda,no_it]=linesearch(func,x,d);
lambda
no_it
func_value = func(x+lambda.*d)

display('Test 2');
x=[0;0];
d=[1;0];
[lambda,no_it]=linesearch(func,x,d);
lambda
no_it
func_value = func(x+lambda.*d)

display('Test 3');
a=2;
func=@(x)(1-10^a*x)^2;
x=0;
d=1;
[lambda,no_it]=linesearch(func,x,d);
lambda
no_it
func_value=func(x+lambda.*d)

display('Test 4');
a=-2;
func=@(x)(1-10^a*x)^2;
x=0;
d=1;
[lambda,no_it]=linesearch(func,x,d);
lambda
no_it
func_value=func(x+lambda.*d)

display('Test 5');
a=5;
func=@(x)(1-10^a*x)^2;
x=0;
d=1;
[lambda,no_it]=linesearch(func,x,d);
lambda
no_it
func_value=func(x+lambda.*d)

display('Test 6');
a=-5;
func=@(x)(1-10^a*x)^2;
x=0;
d=1;
[lambda,no_it]=linesearch(func,x,d);
lambda
no_it
func_value=func(x+lambda.*d)

display('Test 7');
a=10;
func=@(x)(1-10^a*x)^2;
x=0;
d=1;
[lambda,no_it]=linesearch(func,x,d);
lambda
no_it
func_value=func(x+lambda.*d)

display('Test 8');
a=-10;
func=@(x)(1-10^a*x)^2;
x=0;
d=1;
[lambda,no_it]=linesearch(func,x,d);
lambda
no_it
func_value=func(x+lambda.*d)

%% Rosenbrock testing
nonlinearmin(@rosenbrock,[200;200],'DFP',1e-6,1)

