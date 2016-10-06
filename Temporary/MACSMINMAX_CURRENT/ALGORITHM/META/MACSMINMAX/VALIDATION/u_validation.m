function [f,u,nfeval,all_f] = u_validation(problem_fix_d,d_record,u_record,local_search_flag,objectives);

minmaxsign = problem_fix_d.par_objfun.sign;
archsize = size(d_record, 1);
u = cell(1,max(objectives));
f = zeros(archsize,max(objectives));
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
        f_du = [];
        u_du=[];
        f_d=[];
        u_d=[];
        problem_fix_d.par_objfun.objective = obj;
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
                    options = optimset('Display','none','MaxFunEvals',50*problem_fix_d.dim,...%'LargeScale','off',...
                    'Algorithm','sqp'); % add a converged stop condition. in the original there was one but wrongly implemented
                    [u_du,f_du,~,output] = fmincon(func,u_0,[],[],[],[],lb_local,ub_local,[],options,problem_fix_d.par_objfun); %unconstrained
                    % c = output.constrviolation;
                    nfeval = nfeval + output.funcCount;
            else     
                u_du = u_0;
                [f_du] = func (u_du, problem_fix_d.par_objfun);
                nfeval = nfeval + 1;
            end
            if (isempty(f_d) || f_du < f_d)
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

return