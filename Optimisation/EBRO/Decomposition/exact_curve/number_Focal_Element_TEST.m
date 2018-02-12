function [num_FE] = number_Focal_Element_TEST(problem)

num_FE = 1;
for i = 1:problem.dim_u
    
    num_FE = num_FE*length(problem.lb_u{1}{i});
    
end

end