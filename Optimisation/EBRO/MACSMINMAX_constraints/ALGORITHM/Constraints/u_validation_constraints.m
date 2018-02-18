function [f,u,nfeval,all_f, violation] = u_validation_constraints(problem_fix_d,d_record,u_record,local_search_flag,objectives)

violation   = [];
violation_C = [];

minmaxsign = problem_fix_d.par_objfun.sign;
archsize = size(d_record, 1);
u = cell(1,max(objectives));
f = nan(archsize,max(objectives));
for obj = objectives
    u{obj} = nan(archsize,problem_fix_d.dim);
end
nfeval = 0;
func = problem_fix_d.objfun;
all_f=[];
if nargout > 3
    all_f = cell(1,max(objectives)); % all evaluations stored like all_f{obj}(idx_d_record,idx_u_record)
    % careful when local search, we are not storing u_f only u_0 and it's f(d,u_f)
end

for idx_d = 1:archsize
    % fix d and scale
    problem_fix_d.par_objfun.d = d_record(idx_d,:);
    for obj = objectives
        problem_fix_d.par_objfun.objective = obj;
        f_du = [];
        u_du=[];
        f_d=nan;
        u_d=nan(1,problem_fix_d.dim);
        for idx_u = 1:size(u_record{obj},1)
            u_0=u_record{obj}(idx_u,:);
            if (local_search_flag)
                stop = 0;
                iter = 0;
                crowding = 0.1;
                lb_local = u_0 - crowding/2;
                ub_local = u_0 + crowding/2;
                lb_local(lb_local < 0) = 0;
                ub_local(ub_local > 1) = 1;
                options = optimset('Display','none','MaxFunEvals',50*problem_fix_d.dim,'TolFun',1e-8,...%'LargeScale','off',...
                    'Algorithm','sqp'); % add a converged stop condition. in the original there was one but wrongly implemented
                [u_du,f_du,~,output] = fmincon(func,u_0,[],[],[],[],lb_local,ub_local,[],options,problem_fix_d.par_objfun); %unconstrained
                
                
                % c = output.constrviolation;
                nfeval = nfeval + output.funcCount;
            else
                u_du = u_0;
                [f_du] = func (u_du, problem_fix_d.par_objfun);
                
                %----------------------------------------------------------
                % CONSTRAINTS
                if ~isempty(problem_fix_d.fitnessfcn.constr)                 % problem_fix_d.par_objfun.constraints_flag  == 1
                    
                    
                    constraint(1,:) = feval(problem_fix_d.par_objfun.mask_constraints, u_du, problem_fix_d.par_objfun); % @mask_constraints_macsminmax_inner %[g,h]
                    con_row = constraint(1,:);
                    con_sum(1) = sum(con_row(con_row>0));
                    
                    if  (con_sum(1) > 0)
                        
                        %%%%%%%%%%%
%                         stop = 0;
%                         iter = 0;
%                         crowding = 0.1;
%                         lb_local = problem_fix_d.lb;
%                         ub_local = problem_fix_d.ub;
%                         lb_local(lb_local < 0) = 0;
%                         ub_local(ub_local > 1) = 1;
%                         options = optimset('Display','none','MaxFunEvals',problem_fix_d.dim,'TolFun',1e-8,...%'LargeScale','off',...
%                             'Algorithm','sqp'); % add a converged stop condition. in the original there was one but wrongly implemented
%                         [u_du,f_du,~,output] = fmincon(func, u_0, [],[],[],[], lb_local, ub_local, problem_fix_d.par_objfun.mask_constraints, options, problem_fix_d.par_objfun); %unconstrained
%                         
%                         % c = output.constrviolation;
%                         nfeval = nfeval + output.funcCount;
                        %%%%%%%%
%                         f_du = abs(f_d) + 1;
                        violation_C(idx_d, idx_u) = constraint;
                    else
                        violation_C(idx_d, idx_u) = 0;
                    end
                end
                %----------------------------------------------------------
                
                
                nfeval = nfeval + 1;
            end
            if (idx_u==1 || f_du < f_d)
                f_d = f_du;
                u_d = u_du;
            end
            if (nargout > 3)
                all_f{obj}(idx_d,idx_u) = -minmaxsign*f_du;
            end
        end
        f(idx_d,obj) = -minmaxsign*f_d;
        u{obj}(idx_d,:) = u_d;
    end
end

if ~isempty(problem_fix_d.fitnessfcn.constr) 
    violation = max(max(violation_C));
else
    violation = 0;
end


return