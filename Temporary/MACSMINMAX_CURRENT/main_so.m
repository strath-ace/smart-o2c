% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

%%initialisation
init = str2func(strcat('init_macsminmax_evolve'));
% for benchmarking purposes
nfevals_io = [100,500; 100,200; 100,100; 100,500; 200,500; 500,500; 500,500; 100,100; 50,100; 50,100; 50,100; 100,500; 500,1000];

%%
for runid= 1
    for tc = 1:13

        % tic
        global nfevalglobal;
        nfevalglobal = 0;

        %% initialise problem TC1
        [ problem_minmax ] = init_tc_so(tc);

        %% initialise algorithm minmax (META-ALGO)
        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

        %% modify nfeval (for benchmarking purposes)
            algo_outer.par.options(1) = nfevals_io(tc,1);
            algo_inner.par.options(1) = nfevals_io(tc,2);

        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

        %create results directory
        savefolder = strcat('RESULTS/SO/');
        mkdir(savefolder);
        save(strcat(savefolder,'testcase_results_TC_so_',num2str(tc),'_',num2str(runid)));
    end
end