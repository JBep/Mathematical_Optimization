function penalty = penalty(F, g_k, h_k, mu)
%F is the function to be minimized, g_k are the inequality contraints, h_k
%are the equality constraints and mu is the coefficient for the penalty
penalty = @(x) F(x)+mu*alpha(g_k,h_k,x);
end

function alpha = alpha(g_k,h_k,x)
alpha_g = max(func_values(g_k,x),0).^2; %array of penalty for inequalities at point x
alpha_h = func_values(h_k,x).^2; %array of penalty for equalities at point x

alpha = sum(alpha_g) + sum(alpha_h); %sum of both penalties
end

function func_values = func_values(w_k,x) %array of the functional value at point x

values(1)=0;

for i=1:numel(w_k)%iterate through function g_k or h_k and add the functional value at point x to an array
    temp = w_k{i};
    values(i) = temp(x);
end

func_values = values;
end