function barrier = barrier(F, g_k,epsilon)
%F is the function to be minimized, g_k is the matrix of inequality
%constraints and epsilon is the coefficient of beta function
barrier = @(x) F(x)+epsilon*beta(g_k,x);
end

function beta = beta(g_k, x)
g_func_values = g_values(g_k, x); %array of all functional values of contraint function at point x

if any(g_func_values>=0)%checks feasability of point x, if outside barrier beta = inf
    beta = exp(1000);
else
    g_func_values = 1./g_func_values;
    beta = -sum(g_func_values);
end
end

function g_values = g_values(g_k, x)%adds all the functional values of the contraint function (g_k) at point x
func_values(1)=0;

for i=1:numel(g_k)%iterate through function g_k or h_k and add the functional value at point x to an array
    temp = g_k{i};
    func_values(i) = temp(x);
end

g_values = func_values;
end
