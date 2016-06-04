function [t,x] = rk4(ode_fun,t,x0,forward,varargin)

%RK4 integrator
%INPUT:  ode_fun = function of the ODE
%        t = vector with integration time
%        x0 = vector with initial integration point
%        forward = variable that indicates if to integrate forward (1) or
%        backward (2)
%        varargin = in case more variables are needed to the ode_fun
%
%OUTPUT: t = vector with the evaluation times (all the points)
%        x = matrix [m*n], where m is the number of variables (the
%        dimension of x0, and n the number of time evaluation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate parameters and variables
%Calculate the length (for the cycle)
N = length(t);

%Write the initial guess in the output matrix (as column vector)
x = x0(:);

%Write x_old for the cycle (and transpose to obtain row vector)
x_old = x';
 
%% Integration 

if (forward==1)

%Integration forward
%Cycle for all the evaluation times
for i=2:N
    
    % Calculate the time step
    h = t(i) - t(i-1);
    
    %Define the k_i functions of RK4
    k1 = ode_fun(t(i-1),x_old,varargin{:});
    k2 = ode_fun(t(i-1) + h/2 , x_old + h*k1/2 ,varargin{:});
    k3 = ode_fun(t(i-1) + h/2 , x_old + h*k2/2 ,varargin{:});
    k4 = ode_fun(t(i-1) + h , x_old + h*k3 ,varargin{:});
   
    %Calculate new vector x_new and add to the matrix x
    x_new = x_old + h*  (k1 + 2*k2 + 2*k3 + k4)/6;
    x = [x x_new'];
    x_old = x_new;
    
end

else
    
%Integration backwards
%Cycle for all the evaluation times
for i=2:N
    
    % Calculate the time step
    h = t(N-i+1) - t(N-i+2);
    
    %Define the k_i functions of RK4
    k1 = ode_fun(t(N-i+2),x_old,varargin{:});
    k2 = ode_fun(t(N-i+2) + h/2 , x_old + h*k1/2 ,varargin{:});
    k3 = ode_fun(t(N-i+2) + h/2 , x_old + h*k2/2 ,varargin{:});
    k4 = ode_fun(t(N-i+2) + h , x_old + h*k3 ,varargin{:});
   
    %Calculate new vector x_new and add to the matrix x
    x_new = x_old + h*  (k1 + 2*k2 + 2*k3 + k4)/6;
    x = [x x_new'];
    x_old = x_new;
    
end
    
    % Now invert direction of elements of x (since they are from time N to
    % 0 and they shall be from time 0 to N) so its direction is the same of
    % timeVec
    
    x = fliplr(x);
end

end
