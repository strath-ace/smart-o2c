function M = ascent_state_equation(x,u,t)

g = 0.0016; 
T = 4e-3;

M = [x(2); 
    T*cos((u(1))); 
    x(4); 
    T*sin((u(1)))-g]; 

end