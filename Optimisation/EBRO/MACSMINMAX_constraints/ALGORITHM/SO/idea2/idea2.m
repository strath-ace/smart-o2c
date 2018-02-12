function [memories,nfeval,flag_restart]=idea2(func,cfunc,y,vlb0,vub0,options,varargin)
%%
flag_restart=0;
lx=length(vlb0);
ng             = options(1);  % Number of f evaluations
npop           = options(4);  % size of the population
delta          = options(13);
verbose        = options(15); % verbousity of the output
idimp(1)       = options(19);
idimp(2)       = options(20);
cflag          = options(21);
crowding       = options(22);  % crowding factor for global restart
convtol        = options(27);  %  of max distance for local convergence
F              = options(28);
CR             = options(29);
max_local_restart=options(31); % number of local restarts
strategy       = options(32);

if length(options)==34
    thrash=options(34);
    delta=thrash;
else
    thrash=[];
end
mfit= 1;
alphaini=crowding;
MAXmaxd=0;
local_restart=0;
agentgoal=zeros(npop,2);
popold    = zeros(size(y));     % toggle population
fitness   = zeros(1,npop);          % create and reset the "cost array"
bestmem   = zeros(1,lx);           % best population member ever
bestmemit = zeros(1,lx);           % best population member in iteration
nfeval    = 0;                    % number of function evaluations

%% KEY PARAMETERS
%CR=0.1;
%F=0.5;
%%

nfevalmem=0;
memories=[];
vub=vub0;
vlb=vlb0;
fworst=-1e36;

big     = 1e260;
Cf      = big*ones(1,npop);
fworst  =-big;
mu=0;
behav_com=3;
agentgoal  = -ones(npop,2);
nfeval=0;
% migration radius initialisation
alpha=ones(1,npop)*alphaini;
fid=1;
maxd=big;
%==================================================================================
if verbose>=0
    fprintf(fid,'Start IDEA\n')
    fprintf(fid,'N. of agents=%d\n',npop);
    fprintf(fid,'Deploying Agents\n');
end

%% COM SECTION
ibest   = 1;                      % start with first population membe

fitness(1) = func(y(ibest,:),varargin{:});                 % best objective function value so far
bestval=fitness(1);
nfeval  = nfeval + 1;
for i=2:npop                        % check the remaining members
    fitness(i) = func(y(i,:),varargin{:});
    nfeval  = nfeval + 1;
    if (fitness(i) < bestval)           % if member is better
        ibest   = i;                 % save its location
        bestval = fitness(i);
    end
end

bestmemit = y(ibest,:);         % best member of current iteration
bestvalit = bestval;              % best value of current iteration
bestmem = bestmemit;              % best member ever
pm1 = zeros(npop,lx);              % initialize population matrix 1
pm2 = zeros(npop,lx);              % initialize population matrix 2
pm3 = zeros(npop,lx);              % initialize population matrix 3
pm4 = zeros(npop,lx);              % initialize population matrix 4
pm5 = zeros(npop,lx);              % initialize population matrix 5
bm  = zeros(npop,lx);              % initialize bestmember  matrix
ui  = zeros(npop,lx);              % intermediate population of perturbed vectors
mui = zeros(npop,lx);              % mask for intermediate population
mpo = zeros(npop,lx);              % mask for old population
rot = (0:1:npop-1);               % rotating index array (size npop)
rotd= (0:1:lx-1);                % rotating index array (size lx)
rt  = zeros(npop);                % another rotating index array
rtd = zeros(lx);                 % rotating index array for exponential crossover
a1  = zeros(npop);                % index array
a2  = zeros(npop);                % index array
a3  = zeros(npop);                % index array
a4  = zeros(npop);                % index array
a5  = zeros(npop);                % index array
ind = zeros(4);

iter = 1;
done=0;
mem_fmin=[];
jplus=0;
nostop=1;
j=0;
while nostop
    
    %% COM SECTION
    %%
    popold = y;                   % save the old population
    
    ind = randperm(4);              % index pointer array
    
    a1  = randperm(npop);             % shuffle locations of vectors
    rt = rem(rot+ind(1),npop);        % rotate indices by ind(1) positions
    a2  = a1(rt+1);                 % rotate vector locations
    rt = rem(rot+ind(2),npop);
    a3  = a2(rt+1);
    rt = rem(rot+ind(3),npop);
    a4  = a3(rt+1);
    rt = rem(rot+ind(4),npop);
    a5  = a4(rt+1);
    
    pm1 = popold(a1,:);             % shuffled population 1
    pm2 = popold(a2,:);             % shuffled population 2
    pm3 = popold(a3,:);             % shuffled population 3
    pm4 = popold(a4,:);             % shuffled population 4
    pm5 = popold(a5,:);             % shuffled population 5
    
    for i=1:npop                      % population filled with the best member
        bm(i,:) = bestmemit;          % of the last iteration
    end
    
    mui = rand(npop,lx) < CR;          % all random numbers < CR are 1, 0 otherwise
    
    if (strategy > 3)
        st = strategy-3;		  % binomial crossover
    else
        st = strategy;		  % exponential crossover
        mui=sort(mui');	          % transpose, collect 1's in each column
        for i=1:npop
            n=floor(rand*lx);
            if n > 0
                rtd = rem(rotd+n,lx);
                mui(:,i) = mui(rtd+1,i); %rotate column i by n
            end
        end
        mui = mui';			  % transpose back
    end
    
    mpo = mui < 0.5;                % inverse mask to mui
    
    if (st == 1)                      % DE/best/1
        ui = bm + F*(pm1 - pm2);        % differential variation
        ui = popold.*mpo + ui.*mui;     % crossover
    elseif (st == 2)                  % DE/rand/1
        ui = pm3 + F*(pm1 - pm2);       % differential variation
        ui = popold.*mpo + ui.*mui;     % crossover
    elseif (st == 3)                  % DE/rand-to-best/1
        ui = popold + F*(bm-popold) + F*(pm1 - pm2);
        ui = popold.*mpo + ui.*mui;     % crossover
    end
    
    %-----Select which vectors are allowed to enter the new population------------
    for i=1:npop
        
        duil=ui(i,:)-vlb0;
        ioutl=find(duil<0);
        lioutl=length(ioutl);
        duir=vub0-ui(i,:);
        ioutr=find(duir<0);
        lioutr=length(ioutr);
        
        if lioutl>0
            ui(i,ioutl) = vlb0(ioutl) + rand(1,lioutl).*(vub0(ioutl) - vlb0(ioutl));
        end
        if lioutr>0
            ui(i,ioutr) = vlb0(ioutr) + rand(1,lioutr).*(vub0(ioutr) - vlb0(ioutr));
        end
        
        x=ui(i,:);
        
        if norm(x-y(i,:))>0
            
            tempval = func(x,varargin{:});% check cost of competitor
            
            nfeval  = nfeval + 1;
        else
            tempval=fitness(i);
        end
        
        impr=0;
        
        if (tempval < fitness(i))  % if competitor is better than value in "cost array"
            y(i,:) = ui(i,:);  % replace old vector with new one (for new iteration)
            fitness(i)   = tempval;  % save value in "cost array"
            %----we update bestval only in case of success to save time-----------
            if (tempval < bestval)     % if competitor better than the best one ever
                bestval = tempval;      % new best value
                bestmem = ui(i,:);      % new best parameter vector ever
                ibest=i;
            end
            impr=1;
        end
        
        %        [y,fitness,impr,bestval,bestmem,ibest]=select(fitness,y,ui,bestval,tempval,i,impr,bestmem,ibest);
        
        % Introduces pattern search perturbation step if DE fails
        if impr==100
            x=y(i,:);
            ids=randperm(lx);
            dy=x-y(i,:);
            Delta=(vub0-vlb0)/2;
            miny=min(y(:,1:lx));
            maxy=max(y(:,1:lx));
            rho=rand*norm(maxy-miny);
            j=0;
            xmem=x;
            while j<1&&impr==0
                j=j+1;
                x(ids(j))=xmem(ids(j))+(2*rand-1)*rho;
                if x(ids(j))<vlb0(ids(j))
                    x(ids(j))=vlb0(ids(j));
                elseif x(ids(j))>vub0(ids(j))
                    x(ids(j))=vub0(ids(j));
                end
                % evaluates step
                if norm(x-y(i,:))>0
                    tempval = func(x,varargin{:});
                    nfeval  = nfeval + 1;
                else
                    tempval=fitness(i);
                end
                [y,fitness,impr,bestval,bestmem,ibest]=select(fitness,y,ui,bestval,tempval,i,impr,bestmem,ibest);
            end
        end
        
    end %---end for imember=1:npop
    
    
    
    bestmemit = bestmem;       % freeze the best member of this iteration for the coming
    % iteration. This is needed for some of the strategies.
    iter = iter + 1;
    %%  END COM SECTION
    
    %% Evaluate distance
    meany=mean(y(:,1:lx));
    maxd=0;
    for i=1:npop
        maxd=max([maxd norm(y(i,1:lx)-meany)]);
    end
    %     miny=min(y(:,1:lx));
    %     maxy=max(y(:,1:lx));
    %     maxd=norm(maxy-miny);
    MAXmaxd=max([MAXmaxd maxd]);
    %% Restart mechanisms
    if maxd<convtol*MAXmaxd
        flag_restart=flag_restart+1;
        nfev=min([300*lx ng-nfeval]);
        nfev=round(max([nfev 0]));
        arch=[y fitness'];
        foptionsNLP=optimset('Display','none','MaxFunEvals',nfev,...%'LargeScale','off',...
            'Algorithm','sqp');
        [bestf,idbestf]=min(arch(:,lx+1));
        
        [xgrad,fvalgrad,exitflag,output]=fmincon(func,arch(idbestf,1:lx),[],[],[],[],vlb0,vub0,cfunc,foptionsNLP,varargin{:});
        nfeval=nfeval+output.funcCount;
        
        memories=[memories; xgrad fvalgrad fvalgrad 0 ];
        for k=1:npop
            if agentgoal(k,2)<=0
                memories=[memories; y(k,:) fitness(k) agentgoal(k,:) ];
            end
        end
        if local_restart<max_local_restart
            %             disp('local_restart')
            [fbest,idfbest]=min(memories(:,lx+1));
            local_restart=local_restart+1;
            vlb=memories(idfbest,1:lx)-crowding*(vub0-vlb0)/2;
            vub=memories(idfbest,1:lx)+crowding*(vub0-vlb0)/2;
            
            vlb=max([vlb; vlb0]);
            vub=min([vub; vub0]);
            y(1:npop,:)=lhsu(vlb,vub,npop);
            for k=1:npop
                fitness(k,:)=func(y(k,:),varargin{:});
                C=cfunc(y(k,:),varargin{:});
                nfeval=nfeval+1;
                fworst=max([fworst fitness(k,1:mfit)]);
                Cf(k)=fcon_red(fworst,C);
                agentgoal(k,2)=Cf(k);
                agentgoal(k,1)=fitness(k,1);
                
            end
            
            vlb=vlb0;
            vub=vub0;
        else
            maxd=1e36;
            MAXd=maxd;
            local_restart=0;
            vlb=vlb0;
            vub=vub0;
            % clustering of the current archive %
            if lx>1
                [clustCent,data2cluster,cluster2dataCell] = MeanShiftCluster(memories(:,1:lx)',crowding*sqrt(sum((vub-vlb).^2)));
                usc=size(clustCent,2);
                meanclu=clustCent';
            else
                usc=size(memories,1);
                meanclu=memories(:,1:lx);
            end
            filpop=0;
            nattempts=0;
            while filpop<npop&&nattempts<10
                nattempts=nattempts+1;
                xrestart=lhsu(vlb,vub,npop*2);
                for k=1:npop*2
                    yes=1;
                    for kk=1:usc %length(memories(:,1))
                        if norm(xrestart(k,:)-meanclu(kk,:))<crowding*sqrt(sum((vub-vlb).^2)) %memories(kk,1:lx))<0.1
                            yes=0;
                        end
                    end
                    if yes&&filpop<npop
                        filpop=filpop+1;
                        y(filpop,:)=xrestart(k,:);
                        alpha(filpop)=alphaini;
                        fitness(filpop,:)=func(y(filpop,:),varargin{:});
                        C=cfunc(y(filpop,:),varargin{:});
                        
                        nfeval=nfeval+1;
                        fworst=max([fworst fitness(filpop,1:mfit)]);
                        Cf(filpop)=fcon_red(fworst,C);
                        agentgoal(filpop,2)=Cf(filpop);
                        agentgoal(filpop,1)=fitness(filpop,1);
                        
                    end
                end
                
                
            end
            if filpop<npop
                y(filpop+1:npop,:)=lhsu(vlb,vub,npop-filpop);
                for k=1:npop-filpop
                    alpha(filpop+k)=alphaini;
                    fitness(filpop+k,:)=func(y(filpop+k,:),varargin{:});
                    C=cfunc(y(filpop+k,:),varargin{:});
                    
                    nfeval=nfeval+1;
                    fworst=max([fworst fitness(filpop+k,1:mfit)]);
                    Cf(filpop+k)=fcon_red(fworst,C);
                    agentgoal(filpop+k,2)=Cf(filpop+k);
                    agentgoal(filpop+k,1)=fitness(filpop+k,1);
                    
                end
                
                
            end
        end
    end
    
    if verbose==0
        fprintf(fid,'.')
        if(floor(j/100)==j/100)
            fprintf(fid,'%d\n',nfeval)
        end
    elseif verbose==1
        if ~isempty(memories)
            [nfeval/ng maxd min(memories(:,lx+1)) min(fitness(:,1))]
        else
            [nfeval/ng maxd min(fitness(:,1)) min(fitness(:,1))]
        end
        
    end
    %%
    nostop=(nfeval<ng)&&(jplus<20);
    if nfeval==nfevalmem
        jplus=jplus+1;
    else
        jplus=0;
    end
    nfevalmem=nfeval;
end
if verbose>=0
    fprintf(fid,'\nSearch Ended: %d function evaluations\n',nfeval);
end
%-------------------------------------------------------
options(8)=nfeval;
for i=1:length(y(:,1))
    if agentgoal(i,2)<=0
        memories=[memories; y(i,:) fitness(i) agentgoal(i,:)];
    end
end
memories=sortrows(memories,lx+1);

return
%%
function  [y,fitness,impr,bestval,bestmem,ibest]=...
    select(fitness,y,ui,bestval,tempval,i,impr,bestmem,ibest)


if (tempval < fitness(i))     % if competitor is better than value in "cost array"
    y(i,:) = ui(i,:);          % replace old vector with new one (for new iteration)
    fitness(i)   = tempval;    % save value in "cost array"
    impr=1;
    if (tempval < bestval)     % if competitor better than the best one ever
        bestval = tempval;      % new best value
        bestmem = ui(i,:);      % new best parameter vector ever
        ibest=i;
    end
end

return
%%

function [f,Cmax,x_u0,nfeval]=superproblem(x_d,x_u0,fitness,y,upop,Cmin,memories,func,cfunc,vlb_d,vub_d,vlb_u,vub_u,options,peval,flag,varargin)
%
%  [f,Cmax,x_u0,nfeval]=superproblem(x_d,x_u0,fitness,y,upop,Cmin,memories,func,cfunc,vlb_d,vub_d,vlb_u,vub_u,options,peval,flag,varargin)
%  Cost function for local search
%
%
%
%
%

% Massimiliano Vasile 2014
lu=length(vlb_u);
i=0;
lnfeval=0;
while i<length(memories);
    i=i+1;
    [ftest(i),C(i),u(i,:),lfeval]=subproblem(x_d,memories(i).u(1,1:lu),func,cfunc,vlb_d,vub_d,vlb_u,vub_u,options,0,1,varargin{:});
    lnfeval=lnfeval+lfeval;
end
[ftest(i+1),C(i+1),u(i+1,:),lfeval]=subproblem(x_d,x_u0(1,:),func,cfunc,vlb_d,vub_d,vlb_u,vub_u,options,0,1,varargin{:});
lnfeval=lnfeval+lfeval;
[f,id]=max(ftest);
Cmax=C(id);
x_u0=u(id,:);
nfeval=lnfeval;
return
