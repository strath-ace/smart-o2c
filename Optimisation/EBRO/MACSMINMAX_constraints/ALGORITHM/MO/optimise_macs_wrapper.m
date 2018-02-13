function [ x, fval, exitflag, output ] = optimise_macs_wrapper(problem,par)

[ x, fval, exitflag, output ] = optimise_macs(problem.objfun, problem.lb, problem.ub, par, problem.par_objfun);

if (~isempty (par.max_arch_out) && ~isempty(par.max_arch) && par.max_arch_out ~= par.max_arch && size(x,1)>par.max_arch_out)

    oldarch = output.memory;

    %rch = [x fval];
    %archout = archive(1:par.max_arch_out,:);

    [memories,dd,energy,ener2,mins,maxs]=arch_shrk6([],[],oldarch(1:par.max_arch_out,:),0,[],min(fval(1:par.max_arch_out,:)),max(fval(1:par.max_arch_out,:)),size(x,2),size(fval,2),par.max_arch_out);

    [memories,dd,energy,ener2,mins,maxs]=arch_shrk6(memories,dd,oldarch(par.max_arch_out+1:end,:),energy,ener2,min(fval(par.max_arch_out+1:end,:)),max(fval(par.max_arch_out+1:end,:)),size(x,2),size(fval,2),par.max_arch_out);

    x = memories(:,1:size(x,2));
    fval = memories(:,size(x,2)+1:end-2);
    output.memory = memories;
    
end

return
