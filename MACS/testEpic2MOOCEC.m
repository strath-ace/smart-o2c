clear all
sigma=0;
        weights=ones(1,2);

%func='twoimpmask';
func='moo_robust';
test_function='UF1';
func=@(name)cec09(name);
switch test_function
    case 'schwefel'
        problem='soo';
        np=10;
        fp(1)=0;
        xp=zeros(1,np);
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        weights=ones(1,2);        
    case 'rastrigin'
        problem='soo';
        np=10;
        fp(1)=0;
        xp=zeros(1,np);
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        weights=ones(1,2);                
    case 'zdt6'
        problem='zdt';        
        np=10;
        fp(:,1)=0.388:(1-0.388)/500:1;
        fp(:,2)=1-fp(:,1).^2;
        xp=fp;
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        weights=ones(1,2);

    case 'zdt2'
        problem='zdt';                
        np=30;
        fp(:,1)=0.0:1/500:1;
        fp(:,2)=1-fp(:,1).^2;
        xp=fp;

        vlb0=zeros(1,np);
        vub0=ones(1,np);
        weights=ones(1,2);

    case 'zdt3'
        problem='zdt';                
        np=30;
        fp(:,1)=0.0:1/500:1;
        fp(:,2)=(1-sqrt(fp(:,1))-fp(:,1).*sin(10*pi*fp(:,1)));
        xp=fp;
 
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        weights=ones(1,2);
    case 'zdt4'
        problem='zdt';                
        np=10;
        fp(:,1)=0.0:1/500:1;    
        fp(:,2)=1-sqrt(fp(:,1));
        xp=fp;
        vlb0(1)=0;
        vub0(1)=1;

        vlb0(2:np)=-5*ones(1,np-1);
        vub0(2:np)=5*ones(1,np-1);
        weights=ones(1,2);

    case 'kur1'
        load refkur1.mat
        np=3;
        vlb0=-ones(1,np)*5;
        vub0=ones(1,np)*5;
        weights=ones(1,2);
        fp=refkur1(:,4:5);
        xp=refkur1(:,1:3);
    case 'UF1'
        problem='cec09';                
        np=30;
        name='UF1';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF1.dat
        fp=UF1;
        xp=fp;
    case 'UF2'
        problem='cec09';                        
        np=30;
        name='UF2';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF2.dat        
        fp=UF2;
        xp=fp;        
    case 'UF3'
        problem='cec09';                        
        np=30;
        name='UF3';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF3.dat        
        fp=UF3;
        xp=fp;        
    case 'UF4'
        problem='cec09';                        
        np=30;
        name='UF4';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF4.dat        
        fp=UF4;
        xp=fp;        
    case 'UF5'
        problem='cec09';                        
        np=30;
        name='UF5';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF5.dat        
        fp=UF5;
        xp=fp;        
    case 'UF6'
        problem='cec09';                       
        np=30;
        name='UF6';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF6.dat        
        fp=UF6;
        xp=fp;        
    case 'UF7'
        problem='cec09';                        
        np=30;
        name='UF7';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF7.dat        
        fp=UF7;
        xp=fp;        
    case 'UF8'
        problem='cec09';                        
        np=30;
        name='UF8';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF8.dat        
        fp=UF8;
        xp=fp;        
    case 'UF9'
        problem='cec09';                        
        np=30;
        name='UF9';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF9.dat        
        fp=UF9;
        xp=fp;        
    case 'UF10'
        problem='cec09';                        
        np=30;
        name='UF10';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF10.dat
        fp=UF10;
        xp=fp;     
    case 'CF1'
        problem='cec09';                        
        np=30;
        name='CF1';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/CF1.dat
        fp=CF1;
        xp=fp;     
    case 'CF2'
        problem='cec09';                        
        np=30;
        name='CF2';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/CF2.dat
        fp=CF2;
        xp=fp;     
    case 'CF3'
        problem='cec09';                        
        np=30;
        name='CF3';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/CF3.dat
        fp=CF3;
        xp=fp;     
    case 'CF4'
        problem='cec09';                        
        np=30;
        name='CF4';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/CF4.dat
        fp=CF4;
        xp=fp;    
    case 'Cassini'
        problem='cassini';                        
        np=6;
        vlb0=zeros(1,np);
        vub0=ones(1,np);
    case '3imp'
        problem='triimp';                        
        load Test3impPareto_Global09092010noMOEAD.txt
        fp=Test3impPareto_Global09092010noMOEAD;
        xp=fp;
        np=5;
        vlb0=zeros(1,np);
        vub0=ones(1,np);
    case 'multi_const_test'
        problem='multi';                        
        np=2;
        vlb0=-ones(1,np);
        vub0=ones(1,np);  
    case 'moo_robust'
        problem='robust';                        
        np=2;
        vlb0=zeros(1,np);
        vub0=ones(1,np)*10;          
end
test_function
test_opt=2;
switch test_opt
    case 1
% Options for EPIC
%% EPIC MOO
for j=3:3
    for k=3:3
        options(1) = 25000;             % number of function evaluations
        options(4) = j;               % number of individuals
        options(7) = 1e-5;             % convergence tollerance
        options(11)= 1.0;
        options(12) = k-1; %round(options(4)/2.5);      % evolution filter dimension (half to 2/3 of the number of individuals)
        options(22) = 1e-4;            % crowding
        options(26) = 200;             % nparetoset
        options(29) = 2;               % number of objectives
        options(2) =0;                 % depth of branching (1=one layer)
        options(28)=0;                 %  max subdomains
        options(18)=2;                 % max number of branches
        options(21)=0;                % inequality constraints <0
        %options(11)=1;                % alphaini
        options(15)=-1;
        [j,k]
        fprintf('Iter: ');
        for i=1:200
            fprintf('%d,',i)
            x  = lhsu(vlb0,vub0,options(4));
            lx = length(x(1,:));

            schema = [1];
            mode = 'free';

            %% EPIC call
            [memories,memorieslist,block]=epic2(func,x,vlb0,vub0,options,schema,mode,test_function); % Epic 3.2
            [M1(i),M2(i),M3(i),M4(i)]=ZDTmetrics(memories(:,1:np),memories(:,np+1:np+2),xp,fp,sigma,weights);
            %  cum(i).mem=memories;
            plot(memories(:,31),memories(:,32),'.',fp(:,1),fp(:,2))
            
          [i M1(i) M4(i)]
        pause
        end
        memM(j,k).M1=M1;
        memM(j,k).M4=M4;
       
        % cummem=cum(1).mem;
        % for i=2:10
        %     cummem=[cummem; cum(i).mem];
        % end
        mM1(j,k)=mean(M1);
        sM1(j,k)=std(M1);
        mM4(j,k)=mean(M4)/100;
        sM4(j,k)=std(M4)/100;
    end
end
    case 2
%% test MACS2
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
%contr_ratio    = options(17); % Contraction ratio
%cp             = options(19); % constraint flag
%MBHflag        = options(20); % number of MBH steps  !! currently only for unconstrained !!

opt(1)=300000 ; % nfeval 
opt(2)=150;     % popsize
opt(8)=0.2;      %popratio
%opt(11)=30;    % utility update frequency
opt(17)=0.5;   % contraction speed
opt(21)=0;     % pattern to DE
opt(20)=0;     % number of MBH steps

% easy to set
opt(9)=100;    % archive size
opt(18)=1;   %draw flag
opt(19)=0;   % constraints yes/no 
opt(7)=2;    % number of objective functions

% not used
opt(5)=0.9;   %F
opt(6)=0.9;   %CR
opt(4)=1e-3;  % contraction limit
opt(3)=1;  % constant
opt(10)=1e-6; %not used 
opt(12)=4; %not used
opt(13)=0.5; % not used
opt(15)=1;  % not used
opt(16)=1;  % not used
opt(14)=1;    % threshold for utility update
% opt(1)=30000;  % Number of f evaluations
% opt(4)=10;  % size of the population
% opt(13)=0.2;
% opt(14)=0.2;
% opt(15)=0; % verbousity of the output
% opt(21)=0;
% opt(27)=0.25;  % % of max distance for local convergence
% opt(31)=1; % number of local restarts
% opt(32)=2;  % communication heuristics
% opt(22)=0.2;



x=lhsu(vlb0,vub0,opt(2));
clear M1 M4
for i=1:30
    
    switch problem
        case 'soo'
            mem(i).memories_min7=macs7v15c(@simpleSOO,vlb0,vub0,opt,[],[],test_function);
        case 'zdt'
            mem(i).memories_min7=macs7v14c(@ZDtestfun,vlb0,vub0,opt,[],[],test_function);
            [M1(i),M2(i),M3(i),M4(i)]=ZDTmetrics(mem(i).memories_min7(:,1:np),mem(i).memories_min7(:,np+1:np+2),xp,fp,sigma,weights);
            [i M1(i) M4(i) sum(M4<1e-3)*100/i]

            %  cum(i).mem=memories;
            %plot(memories_min7(:,11),memories_min7(:,12),'.',fp(:,1),fp(:,2))
        case 'cec09'
            mem(i).memories_min7=macs7v14c(@test_func_cec09,vlb0,vub0,opt,[],[],name,np);
            memories_min7=mem(i).memories_min7;
        %    memories_min7=macs7v13COAPversion(@test_func_cec09,vlb0,vub0,opt,[],[],name,np);
            %            plot(memories_min7(:,np+1),memories_min7(:,np+2),'.')
              [M1(i),M2(i),M3(i),M4(i)]=ZDTmetrics(memories_min7(:,1:np),memories_min7(:,np+1:np+2),xp,fp,sigma,weights);
            %
            [i M1(i) M4(i) sum(M4<1e-3)*100/i]
        case 'cassini'
            %    mem(i).memories_min7=macs7v14c(@mexspaceartCassiniNOdsmA_MO,vlb0,vub0,opt,[],[]);
            %     memo=macs_red('triimplamb_mod','const0',x,vlb0,vub0,opt,[]);
        case 'triimp'
            %    mem(i).memories_min7=macs7v14c(@triimplamb_mod,vlb0,vub0,opt,[],[]);
            %    mem(i).memories_min7=macs7v13(@triimplamb_mod,vlb0,vub0,opt,[],[]);            
        case 'multi'            
            %    mem(i).memories_min7=macs7v14c(@multi_const_test,vlb0,vub0,opt,[],[]);
        case 'robust'            
            %    mem(i).memories_min7=macs7v14c(@moo_robust_test,vlb0,vub0,opt,[],[],[],1);            
    end
    
end
%save test_macs7v14c_3imp_correct1.mat 

%save test_macs7v14c_Cassini_correct3.mat 
%save test_macs7v14c_zdt6_15k_20pop_0DE_05con_inertiaDE_correct4.mat
end
%%
% for j=1:5
%     for k=1:j+1
%         psCon(j,k)=sum(memM(j,k).M1<0.01)*100/200;
%         psSpr(j,k)=sum(memM(j,k).M4/100<0.00065)*100/200;
%     end
% end
%save multitest_zdt2_run8_43.mat
% save (['test_',test_function,'_',int2str(options(1)),'_',int2str(j),int2str(k),'_noatt'])
  