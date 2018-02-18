% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_algo_so'));
savefolder = strcat('RESULTS/SO/');

% ONLY FOR BENCHMARK REPEATABILITY PURPOSES! {{{
    nfevals_inner_outer = [100,500; 100,200; 100,100; 100,500; 200,500; 500,500; 500,500; 100,100; 50,100; 50,100; 50,100; 100,500; 500,1000];
% }}}

%%
for runid= 1:20     % 20 runs
    for tc = 1:13% 13 testcases

        disp(strcat('SO_',num2str(tc),'_',num2str(runid)))
        global nfevalglobal;
        nfevalglobal = 0;

        %% initialise problem
        [ problem_minmax ] = init_tc_so(tc);

        %% initialise algorithm minmax (META-ALGO)
        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

        % ONLY FOR BENCHMARK REPEATABILITY PURPOSES! {{{
            algo_outer.par.nFeValMax = nfevals_inner_outer(tc,1);
            algo_inner.par.nFeValMax = nfevals_inner_outer(tc,2);
        % }}}

        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

        %create results directory
%         mkdir(savefolder);
%         save(strcat(savefolder,'testcase_results_TC_SO_',num2str(tc),'_',num2str(runid)));

    end
end