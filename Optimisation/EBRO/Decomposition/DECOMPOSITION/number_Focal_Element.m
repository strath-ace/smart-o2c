function [num_FE] = number_Focal_Element(problem, problem_inner)

num_FE = 1;
for i = 1:problem.dim_u
    
    num_FE = num_FE*length(problem_inner.bpa{i});  % num_FE = num_FE*length(problem_inner.bpa{1}{i});
    
end

end