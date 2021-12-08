function [lambda, no_of_iterations]=linesearch(func,x,d)

F = @(lambda) func(x+lambda.*d); %Function F to be minimized
[lambda, no_of_iterations] = armijo(F); %minimizes funciton F in respect to lambda


if isnan(func(x+lambda.*d)) || func(x+lambda.*d)>func(x)
       error('Bad job of the line search!')
end

lambda;
no_of_iterations;

end


function [lambda, no_of_iterations]=armijo(F)%actual linesearch using Armijo's rule

epsillon=0.1;
alpha=2;
max_iterations = 500;
F0=F(0);
[lambda,iterations]=lambda_brack(F,F0,1,max_iterations); %initial suitable lambda and how many iterations it took

F_prime = derivative(F,lambda,0);%%Defining the derivative of function F

T = @(lambda) F0+epsillon*lambda*F_prime;

iterate = F_prime==0;

while iterations<max_iterations && ~iterate
    
    %%condition to check if both armijo statements are satisfied
    if F(lambda)<=T(lambda) && F(lambda*alpha)>=T(lambda*alpha) 
        break;
    end
    
    if ~(F(lambda)<=T(lambda))%condition for graph to be under overshot tangent line, i.e. step length is not too big
        lambda = lambda/alpha;
    elseif ~(F(alpha*lambda)>=T(alpha*lambda))%condition for step length to be big enough
        lambda = lambda*alpha;
    else
        error('Something went wrong in Armijo statements')
    end
    
    iterations=iterations+1;
end

no_of_iterations = iterations;
if(F_prime == 0)
    lambda=0;
end
end

function [lambda_suitable,no_iterations]=lambda_brack(F,F0,lambda,max_iterations)
F0=F0+max(1e-6*abs(F0), 1e-6); %machine tolerance compensation

no_iterations = 0;

while ((isfinite(F(lambda))==0 || F(lambda)>F0) && no_iterations < max_iterations) %%to check the max value for lambda so it is still calculable
    lambda = lambda/2;
    no_iterations = no_iterations+1;
end


while(F(lambda)<F0 && no_iterations < max_iterations) %to check the lambda is not too small/i.e. there is a difference between F(lambda) and F(0)
    lambda = lambda*2;
    no_iterations = no_iterations+1;
end

if(no_iterations==max_iterations)
    display('Too many itarations!')
end

lambda_suitable=lambda;
end

function [diff]=derivative(f,lambda,x)%Finds the numerical derivative of function F at point x
h=lambda*1e-5;%h is put to relevance to the suitable lambda value (big enough, or small enough)
diff = (f(x+h)-f(x))/(h);

if(diff>0)%minimum condition
    diff=0;
end
end


