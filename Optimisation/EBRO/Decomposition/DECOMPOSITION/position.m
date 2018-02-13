function [P] = position(FE, problem_inner,problem_max)

%FE=9;


for j=1:problem_inner.dim
%     n_int(j) = problem_inner.par_objfun.map_u_info{1, 1}.n_int{j,1};
%n_int(j) = length(problem_max.lb_u{1,1}{j,1});
     n_int(j) = length(problem_max.lb_u{j,1});
end
prodotto = prod(n_int);   

dim=problem_inner.dim;
%n_int = 3*ones(1,dim);
%prodotto = prod(n_int);

A =1;
for j = dim:-1:2
    A = A*n_int(j);
    P(j) = ceil(FE*A/prodotto);
    FE = FE-prodotto/A*(P(j)-1);
end
P(1) = FE;

end