% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%% initialisation
init = str2func(strcat('init_ideaminmax_s'));

%%
for runid=1
    for tc = 1:7
        disp(strcat(num2str(tc),'_',num2str(runid)))       
        
        figure()

        global nfevalglobal;
        nfevalglobal = 0;

        %% initialise problem 
        [ problem_minmax ] = init_tc(tc);

        %% initialise algorithm minmax (META-ALGO). different settings for tc2
        % if (tc==2)
        %     [ algo_minmax, algo_outer, algo_inner ] = init_macsminmax_tc2(problem_minmax);
        % else
            [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);
        % end


        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

        %create results directory
        savefolder = strcat('RESULTS/MO/');
        mkdir(savefolder);
        save(strcat(savefolder,'testcase_results_TC_MO_',num2str(tc),'_',num2str(runid)));
    end
end

% error_runs