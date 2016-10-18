% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

for runid=1:20
    for tc = 5:7
        disp(strcat(num2str(tc),'_',num2str(runid)))

        global nfevalglobal;
        nfevalglobal = 0;

        %% initialise problem
        [ problem_minmax ] = init_tc(tc);

        
        par_macs.maxnfeval = 1e3*problem_minmax.dim_d;      % Max function evaluations
        par_macs.popsize = 10;                              % Population size
        par_macs.rhoini = 1;                                % Initial local hypercube size
        par_macs.F = 1;                                     % F
        par_macs.CR = 0.1;                                  % CR
        par_macs.p_social = 0.5;                            % Ratio between elite and total population
        par_macs.max_arch = 20;                             % Output archive size
        par_macs.coord_ratio = 1;                           % Quota of coordinates examined for each set of individualistic actions
        par_macs.contr_ratio = 0.5;                         % Contraction ratio
        par_macs.draw_flag = 0;                             % Print itarations status  
        par_macs.cp = 0;                                    % constraint flag
        par_macs.MBHflag = 0;                               % number of MBH steps

        [dmin,fminmax,exitflag,output] = optimise_macs(@(d)func(d,problem_minmax),problem_minmax.lb_d',problem_minmax.ub_d',par_macs);

        %create results directory
        savefolder = strcat('RESULTS/REF/');
        mkdir(savefolder);
        save(strcat(savefolder,'testcase_results_TC_REF_',num2str(tc),'_',num2str(runid)));
    end
end