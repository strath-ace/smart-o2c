function [problem_minmax] = init_tc(varargin)

if (nargin == 0)
    id = 1;
else
    id = varargin{1};
end


switch (id)
    case 1
        %% initialise problem TC1

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 2;
        problem_minmax.objfun = {@mv1,@mv3};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct,struct};

        % design variables
        dim_d = 2;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(1,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 2;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[-5,-3,-1]},[dim_u,1]),repmat({[-5,-3,-1]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[-4,0,3]},[dim_u,1]),repmat({[-4,0,3]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 2.0e5;

    case 2
        %% initialise problem TC2

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 2;
        problem_minmax.objfun = {@mv2,@mv8};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct,struct};

        % design variables
        dim_d = 8;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(0,[dim_d,1]);
        problem_minmax.ub_d = repmat(3,[dim_d,1]);

        % uncertain variables
        dim_u = 8;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[-5,-3,-1]},[dim_u,1]),repmat({[0,2,3]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[-4,0,3]},[dim_u,1]),repmat({[1,4,2*pi]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1.0e6;

    case 3
        %% initialise problem TC1

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 2;
        problem_minmax.objfun = {@mv2,@em1};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct,struct};

        % design variables
        dim_d = 8;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(1,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 8;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[-5,-3,-1]},[dim_u,1]),repmat({[0,7,12]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[-4,0,3]},[dim_u,1]),repmat({[5,14,20]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1.0e6;

    case 4
        %% initialise problem TC4

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 2;
        problem_minmax.objfun = {@mv8,@mv9};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct,struct};

        % design variables
        dim_d = 2;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(1,[dim_d,1]);
        problem_minmax.ub_d = repmat(3,[dim_d,1]);

        % uncertain variables
        dim_u = 2;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[0,2,3]},[dim_u,1]),repmat({[-pi/2,0,3*pi/4]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[1,4,2*pi]},[dim_u,1]),repmat({[-pi/6,pi,3*pi/2]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 4.0e5;

    case 5
        %% initialise problem TC4

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 2;
        problem_minmax.objfun = {@mv8,@em1};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct,struct};

        % design variables
        dim_d = 4;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(1,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 4;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[0,2,3]},[dim_u,1]),repmat({[0,7,12]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[1,4,2*pi]},[dim_u,1]),repmat({[5,14,20]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1.0e6;

    case 6
        %% initialise problem TC4

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 2;
        problem_minmax.objfun = {@mv10,@mv9};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct,struct};

        % design variables
        dim_d = 1;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(-4,[dim_d,1]);
        problem_minmax.ub_d = repmat(2*pi,[dim_d,1]);

        % uncertain variables
        dim_u = 1;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[pi,5,5.5]},[dim_u,1]),repmat({[-pi/2,0,3*pi/4]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[4,6,2*pi]},[dim_u,1]),repmat({[-pi/6,pi,3*pi/2]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1.0e5;

    case 7
        %% initialise problem TC4

        % type of problem
        problem_minmax.sign_inner = 1;          % -1 will run minmin

        % objectives
        problem_minmax.n_obj = 3;
        problem_minmax.objfun = {@mv2,@mv8,@em1};            %each function is in the form [f_i] = objfun(i)(d,u,par)
        problem_minmax.par_objfun = {struct,struct,struct};

        % design variables
        dim_d = 4;
        problem_minmax.dim_d = dim_d;
        problem_minmax.lb_d = repmat(1,[dim_d,1]);
        problem_minmax.ub_d = repmat(5,[dim_d,1]);

        % uncertain variables
        dim_u = 4;
        problem_minmax.dim_u = dim_u;
        problem_minmax.lb_u = {repmat({[-5,-3,-1]},[dim_u,1]),repmat({[0,2,3]},[dim_u,1]),repmat({[0,7,12]},[dim_u,1])};
        problem_minmax.ub_u = {repmat({[-4,0,3]},[dim_u,1]),repmat({[1,4,2*pi]},[dim_u,1]),repmat({[5,14,20]},[dim_u,1])};

        % maxnfeval: optional to give it here or in the algorithm
        problem_minmax.maxnfeval = 1.0e6;

    otherwise
        error('init_testcase: trying to initialise testcase with unknown id')
end
return