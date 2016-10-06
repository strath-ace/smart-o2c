function [f] = func(d,problem_minmax)
    n_obj = problem_minmax.n_obj;
    f = [];
    for obj = 1:n_obj
        dim =0;
        fobj = 0;
        for di=d
            dim = dim +1;
            objfun = problem_minmax.objfun{obj};
            par_objfun = problem_minmax.par_objfun{obj};
            lbu = problem_minmax.lb_u{obj}{dim};
            ubu = problem_minmax.ub_u{obj}{dim};
            n_intervals = length(lbu);
            fi = -inf;
            for int=1:n_intervals
                u0s = linspace(lbu(int),ubu(int),10);
                fint = -inf;
                for u0 = u0s
                    options = optimset('Display','none','MaxFunEvals',500,'Algorithm','sqp','TolFun',1e-8);
                    [~,fu0,~,~] = fmincon(@(u)mask(u,di,objfun,par_objfun),u0,[],[],[],[],lbu(int),ubu(int),[],options);
                    fu0 = -fu0;
                    fint = max([fint,fu0]);
                end
                fi = max([fint,fi]);
            end
            fobj = fobj + fi;
        end
        f(obj) = fobj;
    end

return