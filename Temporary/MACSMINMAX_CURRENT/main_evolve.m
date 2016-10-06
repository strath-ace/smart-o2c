% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

% algorithm = 'macsminmax';
algorithm = 'macsminmax';

init = str2func(strcat('init_',algorithm,'_evolve'));

nfevals_io = [100,500; 100,200; 100,100; 100,500; 200,500; 500,500; 500,500; 100,100; 50,100; 50,100; 50,100; 100,500; 500,1000];
nfevals_io_minmarek = [50, 50 ; 50, 100; 50, 100 ; 400, 400; 400, 400; 400, 400];

for runid= 1:100
    for tc = 1:13
%         figure(tc)
%         hold on

        % tic
        global nfevalglobal;
        nfevalglobal = 0;

        %% initialise problem TC1
        [ problem_minmax ] = init_tc_so(tc);

        %% initialise algorithm minmax (META-ALGO)
        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

        % if strcmp(algorithm,'macsminmax')
            % set nfeval_d and nfeval_u
            algo_outer.par.options(1) = nfevals_io(tc,1);
            algo_inner.par.options(1) = nfevals_io(tc,2);

        % elseif strcmp(algorithm,'macsminmax')
        %     % set nfeval_d and nfeval_u
        %     algo_outer.par.options(1) = nfevals_io_minmarek(tc-7,1);
        %     algo_inner.par.options(1) = nfevals_io_minmarek(tc-7,2);
        % end
        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);

        %create results directory
        savefolder = strcat('RESULTS/EVOLVE/');
        mkdir(savefolder);
        save(strcat(savefolder,'testcase_results_TC_evolve_',num2str(tc),'_',num2str(runid)));
    end
end

% error_runs