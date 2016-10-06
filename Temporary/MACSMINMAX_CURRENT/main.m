% Reset random numbers generator
s = RandStream('mt19937ar','Seed',(sum(100*clock)));
RandStream.setGlobalStream(s);

% algorithm = 'macsminmax';
algorithm = 'ideaminmax_s';

init = str2func(strcat('init_',algorithm));

% error_runs= zeros(2,1);

% parpool(10);
% mycluster = parcluster('local');
% delete(mycluster.Jobs);
% parfor(runid=1:10,10)
for runid=1:100
    for tc = 12
        disp(strcat(num2str(tc),'_',num2str(runid)))
%         figure(tc)
%         hold on

        % tic
        global nfevalglobal;
        nfevalglobal = 0;

        %% initialise problem TC1
        [ problem_minmax ] = init_tc_so(tc);

        %% initialise algorithm minmax (META-ALGO)
        [ algo_minmax, algo_outer, algo_inner ] = init(problem_minmax);

% try
        %% optimise
        [ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax);
% catch
%     warning(strcat('tc',num2str(tc),' run',num2str(runid),' failed'));
%     error_runs(:,end+1)=[tc;runid];
% end
        % nfevalglobal
        % toc
        
        %create results directory
        savefolder = strcat('RESULTS/',algorithm,'/');
        mkdir(savefolder);
        save(strcat(savefolder,'testcase_results_TC_ideaminmax_',num2str(tc),'_',num2str(runid)));
    end
end

% error_runs