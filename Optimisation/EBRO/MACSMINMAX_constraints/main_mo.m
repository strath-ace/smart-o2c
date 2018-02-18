% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%% initialisation
init = str2func(strcat('init_algo_mo_esteco'));
savefolder = strcat('RESULTS/MO/');

%%
for runid = 1
    for tc = 1:6

        disp(strcat('MO_',num2str(tc),'_',num2str(runid)))
        global nfevalglobal;
        nfevalglobal = 0;

        %% initialise problem 
        [ problem_minmax ] = init_tc(tc);
        
        % ONLY FOR BENCHMARK REPEATABILITY PURPOSES! {{{
%             init = str2func(strcat('init_algo_mo_ausa_tc',num2str(tc)));
        % }}}

        %% initialise algorithm minmax (META-ALGO).
        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);
        
        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);
        
        %create results directory
%         mkdir(savefolder);
%         save(strcat(savefolder,'testcase_results_TC_MO_',num2str(tc),'_',num2str(runid)));
    end
end
