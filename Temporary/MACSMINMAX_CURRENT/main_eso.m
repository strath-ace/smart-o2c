% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%% initialisation
init = str2func('init_ideaminmax_s');

%%
for runid=1
    for tc = 8:13
        disp(strcat(num2str(tc),'_',num2str(runid)))

        global nfevalglobal;
        nfevalglobal = 0;

        %% initialise problem
        [ problem_minmax ] = init_tc_so(tc);

        %% initialise algorithm minmax (META-ALGO)
        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

        %create results directory
        savefolder = strcat('RESULTS/ESO/');
        mkdir(savefolder);
        save(strcat(savefolder,'testcase_results_TC_eso_',num2str(tc),'_',num2str(runid)));
    end
end