function [risul,fobiet,best,funvin]=ottimizza_WT(numindp0,numind0,numop,D)
numindp=numindp0;
nmax=70;
aes0=1;
pscmassi=3;
pam=1.1;
%%limitfile='lim_F_WindTurbine.txt';
Fname='calcola_F_Glen';
limiti=[zeros(D,1) ones(D,1)]; %load(limitfile);
numgen=size(limiti,1);
liminf=repmat(limiti(:,1)',numindp,1);
limsup=repmat(limiti(:,2)',numindp,1);
%load risul.dat;
%y=risul;
%load yi.dat;
%y=[y;yi];

%load y_corr1.dat;
%y=y_corr1(:,2:end);

%y=[y(1:50,:);rand(numindp,numgen).*(limsup-liminf)+liminf];
Y.div=[3e4      % max Altitude Phase 1
    -2.0e3      % min Altitude Phase 1
    5.4         % max Mach Phase 1
    -400        % min Velocity Phase 1
    799.9       % max Mass Flow Phase 1
    -1.1        % min Mass Flow Phase 1
    335.0       % EAS %98e9        % max Press. Dyn. Phase 1 % to modify
    .82e5       % max Altitude Phase 2
    -0.8e5      % min Altitude Phase 2
    7950-460    % max Velocity Phase 2
    -(7900-460) % min Velocity Phase 2
    deg2rad(2)  % max Path Angle Phase 2
    deg2rad(2)  % min Path Angle Phase 2
    -1e4        % min Final Mass Phase 2
    30          % max Axial Acceleration Phase 2
    20          % max normal Acceleration Phase 1
    1000        % ----- Ind
    1000        % ----- Gen
    1e5]';      % ----- Ind*gen
div=Y.div;%[9.6e5, 2e7, 12, 1.0001, 1000, 1000];
divisore=ones(1,size(div,2));
pesi=ones(1,size(div,2));
%y=rand(numindp,numgen).*(limsup-liminf)+liminf;
%y=lhsu(limiti(:,1)',limiti(:,2)',numindp);
y = lhsdesign(numindp,D);
%y=load('y_corr73.dat');
x=[];
igen=0;
Ncl=10;
[fobietf,gcf]=calcola(Fname,y,igen);
[fobietf,funvinf,funv,cdis]=constraints(pesi,div,fobietf,gcf);
[fobiet,x,fit,psc,funvin,cdis,GC]=classifica_nsga2(y,fobietf,funvinf,cdis,gcf,Ncl,pam);

save('starting','fobiet','x','fit','psc','funvin','cdis','GC')

gccurr=GC;
y=x;
x=[];
%[fit]=order(y,fobiet,psc);
pscmass=90
nbest=0;

iord=zeros(size(fit,1),1);
iord=(size(fit,1)-1)/(2-2*pam).*(fit-pam)+1;
best=[];
fobest=[];
pscbest=[];
fitbest=[];
funvbest=[];
for ipcl=1:size(y,1)
    if psc(ipcl)==1
        best=[best;y(ipcl,:)];
        fobest=[fobest;fobiet(ipcl,:)];
        pscbest=[pscbest;psc(ipcl,1)];
        fitbest=[fitbest;fit(ipcl,1)];
        funvbest=[funvbest,funvin(1,ipcl)];
    end
end
nbest=size(best,1);
if nbest==1
   for ipcl=1:size(y,1)
        if iord(ipcl)==2
            best=[best;y(ipcl,:)];
            fobest=[fobest;fobiet(ipcl,:)];
            pscbest=[pscbest;psc(ipcl,1)];
            fitbest=[fitbest;fit(ipcl,1)];
            funvbest=[funvbest,funvin(1,ipcl)];
        end
    end
end
genvec=[];
gamvec=[];
ellecvec=[];
elledif=[];
cevec=[];
valvec=[];
val=size(y,1);
gammam=100;
funvec=[];
vinvec=[];
cdisvec=[];
% while div(1,1)>10000
%     div=div./divisore;
%     aes0=1;
%     [fobietf,funvinf,funv,cdis]=calcola_3funT3(div,y);
%     [fobiet,x,fit,psc,funvin]=classifica2(y,fobietf,funvinf);
% while igen<numop & sum(funvin==0)<size(y,1)/3
%     iterazione
% end
% end
aes0=1;
genvec=[];
gamvec=[];
ellecvec=[];
elledif=[];
cevec=[];
valvec=[];
val=size(y,1);
gammam=100;
funvec=[];
vinvec=[];
cdisvec=[];

while igen<numop
    iterazione
end
val
figure(3)
plot(genvec(:,1),vinvec(:,1),genvec(:,1),vinvec(:,2),genvec(:,1),vinvec(:,3)),legend('Vmin','Vmax','Vmed');
figure(4)
plot(genvec(:,1),cdisvec(:,1),genvec(:,1),cdisvec(:,2)),legend('Cmed1','Cmed2');
figure(5)
plot(genvec(:,1),cdisvec(:,3),genvec(:,1),cdisvec(:,4)),legend('Cmed3','Cmed4');
figure(6)
plot(genvec(:,1),cdisvec(:,5)),legend('Cmed5');

risul=y;
save risul.dat -double -ascii risul;
save fobiet.dat -double -ascii fobiet;
A=[fobiet,risul];
save risultati.out -double -ascii A;

function seed = get_seed ( dummy )

%% GET_SEED returns a random seed for the random number generator.
%
%  Modified:
%
%    16 November 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer DUMMY, a dummy input value, since MATLAB
%    will not allow functions with no arguments.
%
%    Output, integer SEED, a random seed value.
%
time_array = clock;

hour = time_array(4);
minute = time_array(5);
second = time_array(6);

seed = second + 60 * ( minute + 60 * hour );
%
%  We want values in [1,43200], not [0,43199].
%
seed = seed + 1;
%
%  Remap SEED from [1,43200] to [1,I_MAX].
%
seed = 2147483647 * ( seed  / ( 60.0 * 60.0 * 24.0 ) );

seed = floor ( seed );
%
%  Never use a seed of 0.
%
if ( seed == 0 )
    seed = 1;
end

