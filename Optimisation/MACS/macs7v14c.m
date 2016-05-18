function [memories,nfeval]=macs7v14c(func,vlb,vub,options,filename,fileload,varargin)
%
%  memories=macs7v14(func,vlb,vub,options)
%
%  INPUT
%         func    : function handle
%         vlb,vub : boundaries of the search domain
%         options : vector of optimisation parameters
%                   maxnfeval=options(1);
%                   popsize=options(2);
%                   rhoini=options(3);
%                   tolconv=options(4);
%                   F=options(5);
%                   CR=options(6);
%                   mfit=options(7);
%                   p_local=options(8);
%                   max_arch=options(9);
%                   tollocconv=options(10);
%                   supbr_upd_freq=options(11);
%                   T=options(12);
%                   p_arch_vs_pop=options(13);
%                   pigr_thres=options(14);
%                   rest_freq=options(15);
%                   coord_ratio=options(16);
%                   draw_flag=options(18);
%
%  OUTPUT
%         memories : matrix containing the archived solutions and the
%                    associated value of the cost function
%                    memories=[x f]

% CHANGE LOG
% Created by: Federico Zuiani 2011

if length(options)<18
    draw_flag=0;
    if length(options)<16
        options(16)=0.25;
        if length(options)<15
            options(15)=20;
        end
    end
else
    draw_flag=options(18);
end

if draw_flag>=1
    tmp_str=func2str(func);
    if draw_flag>1
        aviobj = avifile('example.avi','compression','None','quality',25);
    end
    fig=figure;
    tmp_igd_ax=1;
end
%%  PARAMETER SETTINGS
%   Essential
maxnfeval      = options(1);  % Max function evaluations
popsize        = options(2);  % Population size
p_local        = options(8);  % Ratio between elite and total population
cpat           = options(21); % pattern to DE vs pattern to local
%supbr_upd_freq = options(11); % Number of iterations between updates of the elite composition
contr_ratio    = options(17); % Contraction ratio
MBHflag        = options(20); % number of MBH steps  !! currently only for unconstrained !!
% Problem definition
mfit           = options(7);  % Number of objectives
max_arch       = options(9);  % Output archive size
cp             = options(19); % constraint flag
%  Desired
rhoini         = options(3);  % Initial local hypercube size
F              = options(5);  % F
CR             = options(6);  % CR
if rhoini==0
    rhoini=1;
end
if F==0
    F=0.9;
end
if CR==0
    CR=0.9;
end
% NOT USED
tolconv        = options(4);  % Convergence tolerance (to reset local hypercube)
n_local        = round(p_local*popsize);
supbr_upd_freq = n_local;
tollocconv     = options(10); % Local convergence tolerance   !! NOT USED !!
T              = n_local;     % options(12); % Size of the individual's neighbourhood subset on which to perform DE steps.
p_arch_vs_pop  = [];          % options(13); % Probability of performing the DE step on the archive or on the population
rest_freq      = options(15); % Frequency of local restart !! NOT USED !!!
coord_ratio    = options(16); % Quota of coordinates examined for each set of individualistic actions !! NOT USED !!!
pigr_thres     = options(14); % Threshold for utility update  !!NOT USED
%% Initialise
lx=length(vlb);
x=lhsu(vlb,vub,popsize);
maxC=zeros(popsize,1);
% Checks number of objectives 
if cp==0
    [f_dummy]=func(x(1,:),varargin{:});
    mfit=length(f_dummy);
    ncon=1;
else
    [f_dummy,c_dummy]=func(x(1,:),varargin{:});
    mfit=length(f_dummy);
    ncon=length(c_dummy);
end
%%         Prototype adaptiviy for CR -> use the utility function to adapt in the main loop
%           CR=ones(1,2);
%                 idCR=dprob(CR/max(CR));
%                 c=idCR/100;
%              CR(idCR)=CR(idCR)+2*(impr-0.5);
%            CR=CR-min(CR);
%%
switch mfit
    case 1
        % Single Objective        
        n_lambda=2;
        tmp_max_arch=round(max_arch*1.5);        
        lambda=[1;1];        
    case 2
        % Bi-objective
        n_lambda=mfit*100;
        tmp_max_arch=round(min([n_lambda max_arch])*1.5);
        alfas=linspace(0,pi/2,n_lambda)';
        alfas=alfas(2:end-1);
        alfas=alfas(randperm(n_lambda-2));
        lambda=[eye(mfit); sin(alfas) cos(alfas)];
    case 3
        % Tri-objective
        n_lambda=mfit*100;
        tmp_max_arch=round(min([n_lambda max_arch])*1.5);
        n_alpha=round(sqrt(n_lambda));
        n_beta=round(n_lambda/n_alpha);
        alphas=linspace(0,1/4,n_alpha)'*2*pi;
        betas=acos(linspace(-0.5,0.5,n_beta)'-0.5)-pi/2;
        lambda=eye(3);
        for i=1:n_alpha
            for j=1:n_beta
                if (i*j~=1)&&~(i==n_alpha&&j==1)&&(i*j~=n_alpha*n_beta)
                    
                    lambda=[lambda;cos(alphas(i))*cos(betas(j)),sin(alphas(i))*cos(betas(j)),sin(betas(j))];
                end
            end
        end
        n_lambda=n_alpha*n_beta;
        lambda(4:end,:)=lambda(3+randperm(n_lambda-3),:);
    otherwise
        % Many objective
        n_lambda=mfit*100;
        tmp_max_arch=round(min([n_lambda max_arch])*1.5);
        lambda=[eye(mfit); lhsu(zeros(1,mfit),ones(1,mfit),n_lambda-mfit)];
end

for i=1:n_lambda
    lambda(i,lambda(i,:)==0)=0.0001;
    lambda(i,:)=lambda(i,:)/norm(lambda(i,:));
end
pigr=ones(1,n_lambda);

v=zeros(popsize,lx);
rho(:,1)=rhoini*ones(popsize,1);
rho(:,2)=zeros(popsize,1);

if nargin>6&&~isempty(varargin)
    fevalstr='func(y,varargin{:})';
else
    fevalstr='func(y)';
end

nfeval=0;
if nargin>5&&exist(fileload,'file')
    load(fileload,'memories','x','f')
    if exist('f','var')&&(size(f,1)>=popsize)
        x=x(1:popsize,:);
        f=f(1:popsize,:);
    else
        ids=randperm(size(memories,1));
        ids=ids(1:popsize);
        x=memories(ids,1:lx);
        f=memories(ids,lx+1:lx+mfit);
    end
    n_mem=size(memories,1);
else
    f=zeros(popsize,mfit);
    
    for i=1:popsize
        y=x(i,:);
        if cp==0            
            f(i,:)=func(y,varargin{:});
            maxC(i)=0;
        else
            [f(i,:),c(i,:)]=func(y,varargin{:});
            maxC(i)=max(c(i,:));
        end
        nfeval=nfeval+1;
        
    end
    
    dom=dominance(f,0);
    cid=maxC;
    j=0;
    for i=1:popsize
        %if cid(i)<=0
        j=j+1;
        memories(j,:)=[ x(i,:) f(i,:) dom(i) cid(i)];
        %end
    end
    n_mem=size(memories,1);
end
delta=max(memories(:,lx+1:lx+mfit),[],1)-min(memories(:,lx+1:lx+mfit),[],1);
z=min(memories(:,lx+1:lx+mfit),[],1);

archsize=length(memories(:,1));

% Checks which is the best approximation of the i-th subproblem
lambda_f_best=zeros(1,mfit);
lambda_x_best=zeros(1,lx);
for i=1:n_lambda
    g_tmp=zeros(1,n_mem);
    for j=1:n_mem
        g_tmp(j)=g_fun(memories(j,lx+1:lx+mfit),lambda(i,:),z);
    end
    [tmp,tmp_id]=min(g_tmp);
    lambda_f_best(i,:)=memories(tmp_id,lx+1:lx+mfit);
    lambda_x_best(i,:)=memories(tmp_id,1:lx);
end
 
% Selects the initial population of subproblems to be treated
if n_local>=mfit
    act_subpr=randperm(n_lambda-mfit)+mfit;
    act_subpr=[1:mfit act_subpr(1:n_local-mfit)];
    if archsize>1
        id_pop_act_subpr=zeros(1,n_local);
        id_social=1:popsize;
        for i=1:n_local
            tmp=Inf;
            for j=id_social
                if isnan(norm((lambda_f_best(i,:)-f(j,:))./delta))
                    error('The objective vector contains NaN or has no variation')
                end
                if norm((lambda_f_best(i,:)-f(j,:))./delta)<tmp
                    tmp=norm((lambda_f_best(i,:)-f(j,:))./delta);
                    id_pop_act_subpr(i)=j;
                end
            end
            x(id_pop_act_subpr(i),:)=lambda_x_best(i,:);
            f(id_pop_act_subpr(i),:)=lambda_f_best(i,:);
            id_social=id_social(id_social~=id_pop_act_subpr(i));
        end
    else
        id_pop_act_subpr=randperm(n_local);
        id_social=setdiff(1:popsize,id_pop_act_subpr);
    end
else
    id_pop_act_subpr=[];
    id_social=1:popsize;
    n_local=0;
end

iter=0;
kkk=1:n_lambda;

%%  MAIN LOOP
while nfeval<maxnfeval
    iter=iter+1;

    all_pop=1;
    if all_pop
        compop=1:popsize;
    else
        compop=id_social;
    end
    
    % Individualistic Moves
    [xtrial,vtrial,ftrial,maxC,nfeval,discarded,rho]=explore(compop,memories,x,v,f,cid,...
        id_pop_act_subpr,act_subpr,lambda,z,nfeval,lx,mfit,F,CR,rho,vlb,vub,coord_ratio,contr_ratio,tolconv,rhoini,func,cp,pigr,MBHflag,cpat,ncon,varargin);
    
    % Individualistic Selection & Archiving
    ztmp=ftrial;
    for i=compop
        ztmp=[ztmp;discarded(i).f];
    end
     
    z=min([z;ztmp],[],1);
    
    for i=compop
                
        if (maxC(i)<=0)
            
            if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z))
                v(i,:)=vtrial(i,:);
                x(i,:)=xtrial(i,:);
                f(i,:)=ftrial(i,:);
                cid(i)=maxC(i);
            end
            
            f_tmp=[ftrial(i,:);discarded(i).f];
            x_tmp=[xtrial(i,:);discarded(i).x];
            c_tmp=[maxC(i);discarded(i).c];
            domtmp=dominance(f_tmp,0);
            x_tmp=x_tmp(domtmp==0,:);
            f_tmp=f_tmp(domtmp==0,:);
            c_tmp=c_tmp(domtmp==0,:);
            
            [domA,domB]=dominant(memories(:,lx+1:lx+mfit),f_tmp);
            
            if sum(domB==0)>0
                memories=[memories(domA==0,:);x_tmp(domB==0,:) f_tmp(domB==0,:) domB(domB==0) c_tmp(domB==0)];
            else
                memories=[memories(domA==0,:)];
            end
                        
        elseif (maxC(i)>0)&&(maxC(i)<cid(i))
            v(i,:)=vtrial(i,:);
            x(i,:)=xtrial(i,:);
            f(i,:)=ftrial(i,:);
            cid(i)=maxC(i);
            f_tmp=[ftrial(i,:);discarded(i).f];
            x_tmp=[xtrial(i,:);discarded(i).x];
            c_tmp=[maxC(i);discarded(i).c];
            domtmp=dominance(c_tmp,0);
            x_tmp=x_tmp(domtmp==0,:);
            f_tmp=f_tmp(domtmp==0,:);
            c_tmp=c_tmp(domtmp==0,:);

            [domA,domB]=dominant(memories(:,lx+mfit+2),c_tmp);
            if sum(domB==0)>0
                memories=[memories(domA==0,:);x_tmp(domB==0,:) f_tmp(domB==0,:) zeros(sum(domB==0),1) c_tmp(domB==0)];
            else
                memories=[memories(domA==0,:)];                
            end
        end
    end
    
    % Social Moves
    delta=max(memories(:,lx+1:lx+mfit),[],1)-min(memories(:,lx+1:lx+mfit),[],1);
    for i=1:n_local
        [xsamp,fsamp,maxCsamp,nfeval]=social(x(id_pop_act_subpr(i),:),f(id_pop_act_subpr(i),:),memories,x,f,cid(id_pop_act_subpr(i)),lambda(act_subpr(i),:),z,delta,rho(id_pop_act_subpr(i)),vlb,vub,nfeval,lx,mfit,T,F,CR,p_arch_vs_pop,func,cp,varargin);
        z=min([z;fsamp],[],1);
        if (maxCsamp<=0)
            if ((g_fun(fsamp,lambda(act_subpr(i),:),z)<g_fun(f(id_pop_act_subpr(i),:),lambda(act_subpr(i),:),z)))
                % v(id_pop_act_subpr(i),:)=xsamp-x(id_pop_act_subpr(i),:);
                f(id_pop_act_subpr(i),:)=fsamp;
                x(id_pop_act_subpr(i),:)=xsamp;
                cid(id_pop_act_subpr(i))=maxCsamp;
            end
            % Adds sample to global archive
            [domA,domB]=dominant(memories(:,lx+1:lx+mfit),fsamp);
            
            if any(domB==0)
                    memories=[memories(domA==0,:);xsamp fsamp domB maxCsamp];
                    archsize=sum(domA==0)+1;
            end
            
        elseif (maxCsamp>0)&&(maxCsamp<cid(id_pop_act_subpr(i)))
            f(id_pop_act_subpr(i),:)=fsamp;
            x(id_pop_act_subpr(i),:)=xsamp;
            cid(id_pop_act_subpr(i))=maxCsamp;
            % Adds sample to global archive
            [domA,domB]=dominant(memories(:,lx+mfit+2),maxCsamp);
            if any(domB==0)
                memories=[memories(domA==0,:);xsamp fsamp domB maxCsamp];
                archsize=sum(domA==0)+1;
            end
            
        end
        
    end
    
    sel=memories(:,lx+mfit+1)==0;
    memories=memories(sel,:);
    archsize=sum(sel);
    %     if (mod(iter,rest_freq)==0)
    %         % Local restart
    %         [f,x,rho,id_pop_act_subpr,nfeval]=local_restart(x,f,id_pop_act_subpr,act_subpr,lambda,z,lx,func,vlb,vub,rho,tollocconv,tolconv,rhoini,nfeval,varargin);
    %         id_social=setdiff(1:popsize,id_pop_act_subpr);
    
    % Archive shrinking    
    if (archsize>tmp_max_arch)
        memories=arch_shrk(memories,lx,mfit,archsize,tmp_max_arch);
        archsize=tmp_max_arch;
    end
    % SUbproblems update
    if n_local&&(mod(iter,supbr_upd_freq)==0)
        [act_subpr,id_pop_act_subpr,pigr,lambda_f_best,x,f]=upd_act_subpr(pigr,lambda,lambda_f_best,x,f,memories,delta,z,lx,mfit,n_local,pigr_thres);
        id_social=setdiff(1:popsize,id_pop_act_subpr);
    end
    
    %% OUTPUT ON SCREEN
    if draw_flag>=1
        subplot(2,2,2)
        plot([0:1],CR/max(CR))
        drawnow
        subplot(2,2,3)
        
        [pp,id_pigr]=sort(-pigr);
        switch mfit
            case 2
                plot(lambda(:,1),lambda(:,2),'.')
                hold on
                plot(lambda(act_subpr,1),lambda(act_subpr,2),'ro')
                plot(lambda(id_pigr(1:n_local),1),lambda(id_pigr(1:n_local),2),'ko')
            case 3
                plot3(lambda(:,1),lambda(:,2),lambda(:,3),'.')
                hold on
                plot3(lambda(act_subpr,1),lambda(act_subpr,2),lambda(act_subpr,3),'r*')
                plot3(lambda(id_pigr(1:n_local),1),lambda(id_pigr(1:n_local),2),lambda(id_pigr(1:n_local),3),'ko')
        end
        drawnow
        %         figure(1)
        subplot(2,2,4)
        [ss,kkk]=sort(lambda(:,1));
        hold off
        switch mfit
            case 2
                %                     hold(AX(1),'off')
                [AX,H1,H2]=plotyy(1:n_lambda,pigr(kkk),1:popsize,rho,'plot','semilogy');
                set(H1,'LineStyle','-.','color','b')
                set(H2,'LineStyle','-.','color','g')
                %                     hold(AX(1),'on')
                set(AX(1),'YLim',[0 1])
                set(AX(2),'YLim',[tolconv rhoini],'YGrid','on')
                axis(AX(2),[0 popsize tolconv 0.5])
                %                     hold(AX(1),'on')
                %                     drawnow
            case 3
                betat=asin(lambda(:,3));
                alphat=atan(lambda(:,2)./lambda(:,1));
                plot3(alphat,betat,pigr,'.')
                axis([0 pi/2 0 pi/2 0 1])
                hold on
        end
        %          figure(1)
        subplot(2,2,1)
        switch mfit
            case 2
                plot(memories(:,lx+1),memories(:,lx+2),'k.',f(:,1),f(:,2),'ro')
                hold on
                plot(f(id_pop_act_subpr(1:2),1),f(id_pop_act_subpr(1:2),2),'gd')
                plot(f(id_pop_act_subpr(3:end),1),f(id_pop_act_subpr(3:end),2),'r*')
            case 3
                plot3(memories(:,lx+1),memories(:,lx+2),memories(:,lx+3),'k.')
                hold on
                plot3(f(id_pop_act_subpr(1:3),1),f(id_pop_act_subpr(1:3),2),f(id_pop_act_subpr(1:3),3),'ro')
                plot3(f(id_pop_act_subpr(4:end),1),f(id_pop_act_subpr(4:end),2),f(id_pop_act_subpr(4:end),3),'r*')
        end
        drawnow
        hold off
        if draw_flag>1
            FF = getframe(fig);
            aviobj = addframe(aviobj,FF);
        end
    end
    
    if (mod(iter,50)==0)&&~isempty(filename)
        save([filename '_tmp'],'memories','x','f','iter')
    end
    
    if draw_flag>=0
        fprintf('Iter: %d - feval: %d/%d.',iter,nfeval,maxnfeval)
        fprintf('\tArch. size: %d.',size(memories,1))
        fprintf('\tNon-dominated pop.: %d/%d.',sum(dominance(f,0)==0),popsize)
        fprintf('\n')
    end
end

% Archive shrinking before output
if archsize>max_arch
    memories=arch_shrk(memories,lx,mfit,archsize,max_arch);
end

if draw_flag>=1
%    figure(1)
    subplot(2,2,1)
    hold off
    switch mfit
        case 2
            plot(memories(:,lx+1),memories(:,lx+2),'k.')
            hold on
        case 3
            plot3(memories(:,lx+1),memories(:,lx+2),memories(:,lx+3),'k.')
            hold on
    end
    drawnow
    if draw_flag>1
        FF = getframe(fig);
        aviobj = addframe(aviobj,FF);
        aviobj=close(aviobj);
    end
end

return
%%
function [xsamp,fsamp,maxC,nfeval]=social(x,f,memories,x_DE,f_DE,cid,lambda,z,delta,rho,vlb,vub,nfeval,lx,mfit,T,F,CR,p_arch_vs_pop,func,cp,arg)
%
% [xsamp,fsamp,bindom,nfeval]=social(x,f,lambda,z,vlb,vub,nfeval,lx,mfit,func,varargin)
%
%  INPUT
%           x       : positions
%           v       : velocities
%           f       : matrix with objective values associated to x
%           rho     : local search space radius
%           i       : individual number
%           vlb,vub : boundaries of the search domain
%           nfeval  : number of function evaluations
%           dact    : action probability distribution
%           lx      : problem dimensions
%           mfit    : number of objective functions
%           func    : function handle to objective functions
%           arg     : argument of func
%
%  OUTPUT
%           xsamp   : sampled values
%           fsamp   : objective functions associated to xsamp
%           nfeval  : number of function evaluations
%           dact    : action probability distribution
%

% CHANGE LOG
% Created by: Massimiliano Vasile 2010

%% DE

n_mem=size(memories,1);
pop_DE=size(x_DE,1);
p_arch_vs_pop=1-exp(-n_mem/pop_DE);
T=max([T 4]);
xsamp=x;
fsamp=f;
csamp=zeros(size(x,1),1);
maxC=csamp;
arch_vs_pop=rand<p_arch_vs_pop;
if arch_vs_pop&&n_mem>4
    tmp_dist=zeros(1,n_mem);
    for j=1:n_mem
        tmp_dist(j)=norm((memories(j,1:lx)-xsamp)./(vub-vlb));
    end
    [tmp_dist,P]=sort(tmp_dist);
    %     P=P(randperm(length(tmp_dist)));
    P=P(randperm(min([T length(tmp_dist)])));
else
    tmp_dist=Inf*ones(1,pop_DE);
    for j=1:pop_DE
        tmp_dist(j)=norm((x_DE(j,:)-xsamp)./(vub-vlb));
        if tmp_dist(j)==0
            tmp_dist(j)=Inf;
        end
    end
    [tmp_dist,id_pop]=sort(tmp_dist);
    id_pop=id_pop(randperm(min([T length(tmp_dist)])));
end

DE_type=2;
switch DE_type
    case 1 % DE/rand/1
        
        e=rand(1,lx)<CR;
        if arch_vs_pop&&n_mem>4&&length(P)>3
            xsamp=xsamp.*(1-e)+e.*(memories(P(3),1:lx)+F*(memories(P(1),1:lx)-memories(P(2),1:lx)));
        else
            xsamp=xsamp.*(1-e)+e.*(x_DE(id_pop(3),:)+F*(x_DE(id_pop(1),:)-x_DE(id_pop(2),:)));
        end
        
    case 2 % DE/current-to-rand/1
        K=rand;
        count=0;
        if arch_vs_pop&&n_mem>4
            while norm(memories(P(3),1:lx)-xsamp)==0||norm(memories(P(1),1:lx)-memories(P(2),1:lx))==0
                count=count+1;
                if count>n_mem
                    break
                end
                P=[P(2:end) P(1)];
            end
            xsamp=xsamp+K*(memories(P(3),1:lx)-xsamp)+K*F*(memories(P(1),1:lx)-memories(P(2),1:lx));
        else
            while norm(x_DE(id_pop(3))-xsamp)==0||norm(x_DE(id_pop(1))-x_DE(id_pop(2)))==0
                count=count+1;
                if count>pop_DE
                    break
                end
                id_pop=[id_pop(2:end) id_pop(1)];
            end
            xsamp=xsamp+K*(x_DE(id_pop(3))-xsamp)+K*F*(x_DE(id_pop(1))-x_DE(id_pop(2)));
        end
end
for j=1:lx
    if (xsamp(j)<vlb(j))
        xsamp(j)=vlb(j)+rand*(x(j)-vlb(j));
    elseif (xsamp(j)>vub(j))
        xsamp(j)=vub(j)-rand*(vub(j)-x(j));
    end
end
if norm(xsamp-x)>0
    if cp==0
        fsamp=func(xsamp,arg{:});
        maxC=0;
    else
        [fsamp,csamp]=func(xsamp,arg{:});
        maxC=max(csamp);
        if cid<=0&&maxC>0
            fsamp=f+maxC;
        end
    end
    nfeval=nfeval+1;
end
return
%%
function [xtrial,vtrial,ftrial,maxC,nfeval,discarded,rho]=explore(compop,memories,x,v,f,cid,...
         id_pop_act_subpr,act_subpr,lambda,z,nfeval,lx,mfit,F,CR,rho,vlb,vub,coord_ratio,contr_ratio,tolconv,rhoini,func,cp,pigr,MBHflag,cpat,ncon,arg)
%
%  [x,v,f,nfeval]=explore(memories,box,idsamp_number,x,v,f,nfeval,lx,F,CR,archsize,popsize,vlb,vub,func,arg)
%
%
%  INPUT
%           memories       : archive
%           box            : boxes
%           idsamp_number  : selected box id
%           x              : positions
%           v              : velocities
%           f              : objective values
%           nfeval         : number of function evaluations
%           lx             : problem dimensions
%           F              : perturbation factor
%           CR             :
%           archsize       : archive size
%           popsize        : population size
%           vlb,vub        : boundaries of the search space
%           func           : function handle to objective function
%           arg            : arguments of func
%
%  OUTPUT
%           x       : updated x
%           v       : updated velocity
%           f       : updated function value
%           nfeval  : updated number of function evaluations

% CHANGE LOG
% Created by: Massimiliano Vasile 2010


vtrial=zeros(size(x));
ctrial=zeros(size(x,1),ncon);
maxC=zeros(size(x,1),1);
xtrial=x;
ftrial=f;
discarded.f=[];
discarded.x=[];
discarded.c=[];
Delta=(vub-vlb)/2;
[npop,nf]=size(f);
for i=compop
    discarded(i).f=[];
    discarded(i).x=[];
    discarded(i).c=[];
    
    n_ids=ceil(lx*coord_ratio);
    impr=0;
    ids=randperm(lx);
    r=-1+2*rand(1,n_ids);
    %%
    %          INERTIA
    %          continue along the direction of the previous improvement
    if norm(v(i,:))~=0
        alph=rand;
        % Computes contraction of the step length
        for jj=1:lx
            if((v(i,jj)+x(i,jj))>vub(jj))&&(v(i,jj)~=0)
                alph=min([alph,abs((vub(jj)-x(i,jj))/v(i,jj))]);
            elseif ((v(i,jj)+x(i,jj))<vlb(jj))&&(v(i,jj)~=0)
                alph=min([alph,abs((x(i,jj)-vlb(jj))/v(i,jj))]);
            end
        end
        
        xtrial(i,:)=alph*v(i,:)+x(i,:);
        
        % Double check feasibility error
        for jj=1:lx
            if((xtrial(i,jj)-vlb(jj))<0)
                xtrial(i,jj)=vlb(jj);
            elseif((xtrial(i,jj)-vub(jj))>0)
                xtrial(i,jj)=vub(jj);
            end
        end
        % evaluates move
        if norm(xtrial(i,:)-x(i,:))>0
            
            if cp==0
                ftrial(i,:)=func(xtrial(i,:),arg{:});
                maxC(i)=0;
                if ~all(ftrial(i,:)>=f(i,:))
                    discarded(i).f=[discarded(i).f; ftrial(i,:)];
                    discarded(i).x=[discarded(i).x; xtrial(i,:)];
                    discarded(i).c=[discarded(i).c; maxC(i)];
                    impr=1;
                end
            else
                [ftrial(i,:),ctrial(i,:)]=func(xtrial(i,:),arg{:});
                maxC(i)=max(ctrial(i,:));
                if cid(i)<=0&&maxC(i)>0
                    ftrial(i,:)=f(i,:)+maxC(i);
                end
                if (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                    discarded(i).f=[discarded(i).f; ftrial(i,:)];
                    discarded(i).x=[discarded(i).x; xtrial(i,:)];
                    discarded(i).c=[discarded(i).c; maxC(i)];
                    impr=1;
                end
            end
            
            nfeval=nfeval+1;
            
        end
    end
     
    %%
    %      DIFFERENTIAL
    if impr==0
        %    take the direction suggested by other agents
        ind = randperm(4);              % index pointer array
        rot = (0:1:npop-1);               % rotating index array (size NP)
        
        a1  = randperm(npop);             % shuffle locations of vectors
        rt = rem(rot+ind(1),npop);        % rotate indices by ind(1) positions
        a2  = a1(rt+1);                 % rotate vector locations
        rt = rem(rot+ind(2),npop);
        a3  = a2(rt+1);
        
        e = rand(1,lx) < CR;          % all random numbers < 0.5 are 1, 0 otherwise
        
        strategy='rand';%'best';
        switch strategy
            case 'best'
                af=randperm(nf);
                for j=1:nf
%                    [dummy,idf(j)]=min(f(:,j));
                    [dummy,idf(j)]=min(memories(:,lx+j));                                        
                end                
                %xbest=x(idf(af(1)),:)                
                xbest=memories(idf(af(1)),1:lx);                
                dx =((xbest-x(i,:)) + F*(x(a1(i),:) - x(a2(i),:)));       % differential variation
            case 'rand'
                dx =(x(a3(i),:)-x(i,:) + F*(x(a1(i),:) - x(a2(i),:)));       % differential variation
        end
        xtrial(i,:)=e.*dx+x(i,:);
        
        
        % Double check feasibility error
        for jj=1:lx
            if((xtrial(i,jj)-vlb(jj))<0||(xtrial(i,jj)-vub(jj))>0)
                xtrial(i,jj)=rand*(vub(jj)-vlb(jj))+vlb(jj);
            end
        end
        
        % evaluates move
        if norm(xtrial(i,:)-x(i,:))>0           
            if cp==0
                ftrial(i,:)=func(xtrial(i,:),arg{:});
                maxC(i)=0;
                if ~all(ftrial(i,:)>=f(i,:))
                    discarded(i).f=[discarded(i).f; ftrial(i,:)];
                    discarded(i).x=[discarded(i).x; xtrial(i,:)];
                    discarded(i).c=[discarded(i).c; maxC(i)];
                    %---------------------------------------------------------------                        
                    %  CAREFUL THIS HAS TO BE TESTED
                    if cpat
                    vtrial(i,:)=xtrial(i,:)-x(i,:);
                    end
                    %---------------------------------------------------------------                        
                    impr=1;
                end
            else
                
                [ftrial(i,:),ctrial(i,:)]=func(xtrial(i,:),arg{:});
                maxC(i)=max(ctrial(i,:));
                if cid(i)<=0&&maxC(i)>0
                    ftrial(i,:)=f(i,:)+maxC(i);
                end
                
                if (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                    discarded(i).f=[discarded(i).f; ftrial(i,:)];
                    discarded(i).x=[discarded(i).x; xtrial(i,:)];
                    discarded(i).c=[discarded(i).c; maxC(i)];
                    %---------------------------------------------------------------                        
                    %  CAREFUL THIS HAS TO BE TESTED
                    if cpat
                    vtrial(i,:)=xtrial(i,:)-x(i,:);
                    end
                    %---------------------------------------------------------------                        
                    
                    impr=1;
                end
            end
            
            nfeval=nfeval+1;
            
        else
            ftrial(i,:)=f(i,:);
        end
    end
    %%
    %   PATTERN SEARCH
    if impr==0
 
         xtrial(i,:)=cpat*xtrial(i,:)+(1-cpat)*x(i,:);
%         ftrial(i,:)=f(i,:);
        xsample=[];
        fsample=[];
        
        for j=1:n_ids
            xtrial(i,ids(j))=xtrial(i,ids(j))+r(j)*rho(i,1)*Delta(ids(j));
            if xtrial(i,ids(j))<vlb(ids(j))
                xtrial(i,ids(j))=vlb(ids(j));
            elseif xtrial(i,ids(j))>vub(ids(j))
                xtrial(i,ids(j))=vub(ids(j));
            end
            
            if norm(xtrial(i,:)-x(i,:))>0
                
                if cp==0
                    ftrial(i,:)=func(xtrial(i,:),arg{:});
                    maxC(i)=0;
                    if ~all(ftrial(i,:)>=f(i,:))
                        discarded(i).f=[discarded(i).f; ftrial(i,:)];
                        discarded(i).x=[discarded(i).x; xtrial(i,:)];
                        discarded(i).c=[discarded(i).c; maxC(i)];
                        impr=1;
                    end
                else
                    [ftrial(i,:),ctrial(i,:)]=func(xtrial(i,:),arg{:});
                    maxC(i)=max(ctrial(i,:));
                    if cid(i)<=0&&maxC(i)>0
                        ftrial(i,:)=f(i,:)+maxC(i);
                    end
                    %                     xtrial(i,:)
                    %                     ftrial(i,:)
                    if (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                        discarded(i).f=[discarded(i).f; ftrial(i,:)];
                        discarded(i).x=[discarded(i).x; xtrial(i,:)];
                        discarded(i).c=[discarded(i).c; maxC(i)];
                        impr=1;
                    end
                end
                nfeval=nfeval+1;
                
            end
            
            
            if (all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(maxC(i)>0&&maxC(i)>cid(i))
                
                xsample=[xsample;xtrial(i,:)];
                fsample=[fsample;ftrial(i,:)];
                
                if r(j)>0
                    rr=rand;
                else
                    rr=-rand;
                end
                xtrial(i,ids(j))=x(i,ids(j))-rr*rho(i,1)*Delta(ids(j));
                if xtrial(i,ids(j))<vlb(ids(j))
                    xtrial(i,ids(j))=vlb(ids(j));
                elseif xtrial(i,ids(j))>vub(ids(j))
                    xtrial(i,ids(j))=vub(ids(j));
                end
                
                if norm(xtrial(i,:)-x(i,:))>0
                    
                    if cp==0
                        ftrial(i,:)=func(xtrial(i,:),arg{:});
                        maxC(i)=0;
                        if ~all(ftrial(i,:)>=f(i,:))
                            discarded(i).f=[discarded(i).f; ftrial(i,:)];
                            discarded(i).x=[discarded(i).x; xtrial(i,:)];
                            discarded(i).c=[discarded(i).c; maxC(i)];
                            impr=1;
                        end
                    else
                        [ftrial(i,:),ctrial(i,:)]=func(xtrial(i,:),arg{:});
                        maxC(i)=max(ctrial(i,:));
                        if cid(i)<=0&&maxC(i)>0
                            ftrial(i,:)=f(i,:)+maxC(i);
                        end
                        %                         xtrial(i,:)
                        %                         ftrial(i,:)
                        if (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                            discarded(i).f=[discarded(i).f; ftrial(i,:)];
                            discarded(i).x=[discarded(i).x; xtrial(i,:)];
                            discarded(i).c=[discarded(i).c; maxC(i)];
                            impr=1;
                        end
                    end
                    
                    nfeval=nfeval+1;
                    
                    if (all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(maxC(i)>0&&maxC(i)>cid(i))
                        
                        xsample=[xsample;xtrial(i,:)];
                        fsample=[fsample;ftrial(i,:)];
                        xtrial(i,ids(j))=x(i,ids(j));
                        ftrial(i,:)=f(i,:);
                    elseif (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                        impr=1;
                        vtrial(i,:)=xtrial(i,:)-x(i,:);
                    end
                else
                    xtrial(i,ids(j))=x(i,ids(j));
                    ftrial(i,:)=f(i,:);
                end
                
            elseif (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                impr=1;
                vtrial(i,:)=xtrial(i,:)-x(i,:);
            end
            if impr==1           
                break
            end
        end
        
    end
%     impr
%     pause

    if impr==0 
%         kkk=0;
%         for kk=1:2:length(fsample(:,1))-1
%             kkk=kkk+1;
%             for l=1:mfit
%                 dxsamp(kkk,:,l)=quad1D([x(i,:);xsample(kk,:);xsample(kk+1,:)],[f(i,l);fsample(kk,l);fsample(kk+1,l)],0);
%             end
%         end

%         [i pigr(act_subpr(id_pop_act_subpr==i))]
%         pause
        if  any(i==id_pop_act_subpr)&&MBHflag>0&&pigr(act_subpr(id_pop_act_subpr==i))<0.8%<- TO BE REVISED WITH A PROPPER PARAMETER TO BE SET BY THE USER
            niter=0;
            fu=1e37;
            fcurrent=g_MBHfun(x(i,:),func,lambda(act_subpr(id_pop_act_subpr==i),:),x(i,:),arg);
            while fu>=fcurrent&&niter<MBHflag                
                niter=niter+1
                pause
                foptionsNLP=optimset('Display','none','MaxFunEvals',100*lx,'LargeScale','off','Algorithm','sqp','TolFun',1e-8);
                rvlb=x(i,:)-rho(i,1)*Delta;
                rvub=x(i,:)+rho(i,1)*Delta;
                rvlb=max([rvlb;vlb]);
                rvub=min([rvub;vub]);
                xtrial(i,:)=x(i,:)+(2*rand-1)*rho(i,1)*Delta;
                xtrial(i,:)=max([xtrial(i,:);rvlb]);
                xtrial(i,:)=min([xtrial(i,:);rvub]);
                Z=x(i,:);
                [u,fu,exitflag,output]=fmincon(@(x)g_MBHfun(x,func,lambda(act_subpr(id_pop_act_subpr==i),:),Z,arg),xtrial(i,:),[],[],[],[],rvlb,rvub,@(x)c_MBHfun(x,func,lambda(act_subpr(id_pop_act_subpr==i),:),Z,arg),foptionsNLP);
                nfeval=nfeval+output.funcCount;
                if fu<fcurrent
                    xtrial(i,:)=u;
                    ftrial(i,:)=func(xtrial(i,:),arg{:});                
                    nfeval=nfeval+1;                    
                    impr=1;
                end
            end
        end
      
    
    end
    if impr==0 
        rho(i,1)=rho(i,1)*contr_ratio;
 %       rho
 %       pause
        rho(i,2)=rho(i,2)+1;
        if rho(i,2)>3
            rho(i,1)=rhoini;
            rho(i,2)=0;
        end
    end
end

return
%%
% function [memories,x,f,archsize]=upd_arch(memories,archsize,x,f,popsize,id,xsamp,fsamp,nsamp,lx,mfit)
% 
% for j=1:nsamp
%     sempdom=archsize*ones(archsize,1);
%     for k=1:archsize
%         sempdom(k)=sum(fsamp(j,:)>= memories(k,lx+1:lx+mfit));
%         if sempdom(k)==mfit
%             break
%         elseif sempdom(k)==0
%             memories(k,lx+mfit+1)=memories(k,lx+mfit+1)+1;
%         end
%     end
%     if max(sempdom)<mfit
%         archsize=archsize+1;
%         memories(archsize,1:lx)=xsamp(j,1:lx);
%         memories(archsize,lx+1:lx+mfit)=fsamp(j,:);
%         memories(archsize,1+lx+mfit)=0;
%     end
%     lock=0;
%     for k=1:popsize
%         if k~=id
%             sempdom=sum(fsamp(j,:)>= f(k,:));
%             if sempdom==0&&lock==0
%                 lock=1;
%                 f(k,:)=fsamp(j,:);
%                 x(k,:)=xsamp(j,:);
%                 break
%             end
%         end
%     end
% end
% 
% return
%%
function memories=arch_shrk(memories,lx,mfit,archsize,max_arch)

% Extreme points
[maximums,id_max]=max(memories(:,lx+1:lx+mfit),[],1);

% Archive shrinking
dd=Inf*ones(archsize);
for ip=1:(archsize-1)
    for jp=(ip+1):archsize
        dmem=((memories(ip,lx+1:lx+mfit)-memories(jp,lx+1:lx+mfit)));
        dd(ip,jp)=dmem*dmem';
        dd(jp,ip)=dd(ip,jp);
    end
end

selected=zeros(max_arch,1);
selected(1:mfit)=id_max;
notselected=1:archsize;
for i=1:mfit
    notselected=notselected(notselected~=selected(i));
end
for i=mfit+1:max_arch
    [tmp,id2]=max(min(dd(selected(1:i-1),notselected)));
    selected(i)=notselected(id2);
    notselected=notselected(notselected~=selected(i));
end

memories=memories(selected,:);

return
%% SCALARISATION
function g=g_fun(f,lambda,z)

g=max(lambda.*abs(f-z));

return
%%  Cost Function for MBH Step
function g=g_MBHfun(x,func,lambda,z,arg)

f=func(x,arg{:});
g=-(f-z)*(f-z)';

return
%%  Constraint Function for MBH Step
function [c,ceq]=c_MBHfun(x,func,lambda,z,arg)

f=func(x,arg{:});
c(1)=0;
ceq(1)=0;
lf=length(f);
index=radnperm(lf);
for i=1:lf-1
    c(i)=lambda(index(i+1))*(f(index(i))-z(index(i)))-lambda(index(i))*(f(index(i+1))-z(index(i+1)));
end
ceq=c;
return
%%
function [act_subpr,id_pop_act_subpr,pigr,lambda_f_best,x,f]=upd_act_subpr(pigr,lambda,lambda_f_best_old,x,f,memories,delta,z,lx,mfit,n_local,pigr_thres)

n_lambda=size(lambda,1);
n_mem=size(memories,1);
if n_mem==1
    delta=ones(1,mfit);
end
% Checks which is the best approximation of the i-th subproblem and
% computes the new value for pigr
lambda_x_best=zeros(n_lambda,lx);
lambda_f_best=zeros(n_lambda,mfit);
for i=1:n_lambda
    g_tmp=zeros(1,n_mem);
    for j=1:n_mem
        g_tmp(j)=g_fun(memories(j,lx+1:lx+mfit),lambda(i,:),z);
    end
    [tmp,tmp_id]=min(g_tmp);
    pigr_thres=max(g_tmp);    
    lambda_x_best(i,:)=memories(tmp_id,1:lx);
    lambda_f_best(i,:)=memories(tmp_id,lx+1:lx+mfit);
    old=g_fun(lambda_f_best_old(i,:),lambda(i,:),z);
    new=g_fun(lambda_f_best(i,:),lambda(i,:),z);
    Delta=(old-new);%/old;
    %     if Delta<0
    %         Delta=-pigr_thres;
    %         [lambda_g_best_old(i) lambda_g_best(i) g_fun(lambda_f_best_old(i,:),lambda(i,:),z)]
    %         [lambda_f_best_old(i,:);lambda_f_best(i,:)]
    %         pause
    %     end
    if Delta>pigr_thres
        pigr(i)=1;
    else
        pigr(i)=(0.95+0.05*Delta/pigr_thres)*pigr(i);
    end
end
% Selects the new population of subproblems to be treated
if mfit>1
tourn_size=round(n_lambda/60);
else
tourn_size=round(n_lambda/2);    
end
% tourn_size=n_lambda;
tmp_id=(mfit+1):n_lambda;
act_subpr=[1:mfit zeros(1,n_local-mfit)];
for i=mfit+1:n_local
    sel=randperm(length(tmp_id));
    sel=sel(1:min([tourn_size length(sel)]));
    [tmp,id]=max(pigr(tmp_id(sel)));
    act_subpr(i)=tmp_id(sel(id));
    tmp_id=tmp_id(tmp_id~=act_subpr(i));
end

% if n_mem>1 % This if can actually be eliminated
id_pop_act_subpr=zeros(1,n_local);
list=1:size(f,1);
for i=1:n_local
    g_tmp=Inf;
    for j=list
        tmp=g_fun(f(j,:),lambda(act_subpr(i),:),z);
        if tmp<g_tmp
            g_tmp=tmp;
            id_pop_act_subpr(i)=j;
        end
    end
    list=list(list~=id_pop_act_subpr(i));
    %         if rand<0.1
    %             f(id_pop_act_subpr(i),:)=lambda_f_best(act_subpr(i),:);
    %             x(id_pop_act_subpr(i),:)=lambda_x_best(act_subpr(i),:);
    %         end
end
% else
%     id_pop_act_subpr=randperm(n_local);
% end
% n_local
return
%%

