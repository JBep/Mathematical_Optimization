function [x_opt,opt_value] = nonlinearmin(f,start,method,tol,printout)

if strcmp(method,'DFP') == 1

    n = numel(start);
    yplus1 = start;
    D = eye(n);
    k = 0;
    max_iterations = 500;
    x_k = 1.1.*start;
    x_k1 = start;
    if printout == 1
      disp("iteration           x          step size       f(x)      norm(grad)     ls iters      lambda");
    end
    while k<max_iterations%(norm(x_k1-x_k) > tol && norm(f(x_k1)-f(x_k)) > tol) 
        j=0;
        lambda=1;
    while j<n+1 && ~lambda==0
        y = yplus1;
        direction = -D*grad(f,y);
        [linesearchres1,linesearchres2] = linesearch(f,y,direction); 
        lambda = linesearchres1;
        yplus1 = y + lambda.*direction;
        p = lambda.*direction;
        q = grad(f,yplus1) - grad(f,y);
        D = D + (p*p')/(p'*q) - (D*(q*q')*D)/(q'*D*q);
        j=j+1;
    end
    k = k+1;
    x_k = x_k1;
    x_k1 = yplus1;
    D = eye(n);
    
    if printout == 1
    fprintf('%0.0f %25f %12.4f %12.4f %12.4f %12.4f %12.4f\n',k,x_k1(1),norm(x_k - x_k1),f(x_k1),norm(grad(f,x_k1)),linesearchres2,linesearchres1);
            for i = 2:n
        fprintf('%28f\n',x_k1(i));
            end
           
    end

    if ~((norm(x_k1-x_k) > tol && norm(f(x_k1)-f(x_k)) > tol))%termination criteria
        break;
    end
    
    x_opt = x_k1;
    opt_value = f(x_opt);
    end
    
    
    
elseif strcmp(method,'BFGS')
    n = numel(start);
    yplus1 = start;
    D = eye(n);
    max_iterations = 500;
    k = 0;
    x_k = 1000.*start;
    x_k1 = start;
    if printout == 1
      
      disp("iteration           x          step size       f(x)      norm(grad)     ls iters      lambda");
    end
    while k<max_iterations%(norm(x_k1-x_k) > tol && norm(f(x_k1)-f(x_k)) > tol) 
        j=0;
        lambda=1;
    while j<n+1 && ~lambda==0
        y = yplus1;
        direction = -D*grad(f,y);
        [linesearchres1,linesearchres2] = linesearch(f,y,direction); 
        lambda = linesearchres1;
        yplus1 = y + lambda.*direction;
        p = lambda.*direction;
        q = grad(f,yplus1) - grad(f,y);
        D = D + (1+(q'*D*q)/(p'*q))*(p*p')/(p'*q)-(1/(p'*q))*(p*q'*D+D*q*p');
        j=j+1;
    end
    k = k+1;
    x_k = x_k1;
    x_k1 = yplus1;
     D = eye(n);
    
    if printout == 1
    fprintf('%0.0f %25f %12.4f %12.4f %12.4f %12.4f %12.4f\n',k,x_k1(1),norm(x_k - x_k1),f(x_k1),norm(grad(f,x_k1)),linesearchres2,linesearchres1);
            for i = 2:n
        fprintf('%28f\n',x_k1(i));
            end
           
    end
    
    if ~(norm(x_k1-x_k) > tol && norm(f(x_k1)-f(x_k)) > tol)%termination criteria
        break;
    end

    x_opt = x_k1;
    opt_value = f(x_opt);
    end    
   
           
       else
           error('Wrong method input!')
           
end
end
       