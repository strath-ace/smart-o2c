function [fobiet,gc]=calcola_F_Glen(ya,igen)
global population;

itorque=0;
[m,n]=size(ya);

ll=[deg2rad(ones(1,17).*(0))  ones(1,9)*0.5   80   80
    deg2rad(ones(1,17)*30)    ones(1,9)     1800 1800];
if m>=1 && itorque==0
    for i=1:m
        x=(ll(2,:)-ll(1,:)).*ya(i,:)+ll(1,:);
        %try
            [Jsimul]=mySkylon(x);
            
            
            fobiet(i,1)=(1e6-Jsimul(19))./1e6;
            fobiet(i,2)=(1e6-Jsimul(19))./1e6;
            %                                                ||  
%            gc(i,:)=[J(5) J(2) -J(2) J(13) -J(14) J(15)  J(16) -J(16) J(17) -J(17) J(18) -J(18) -J(19) J(28) i igen igen*m+i];
            gc(i,:)=[Jsimul(5) -Jsimul(9) Jsimul(end-1) -Jsimul(10) Jsimul(13) -Jsimul(14) Jsimul(15)  Jsimul(16) -Jsimul(16) Jsimul(17) -Jsimul(17) Jsimul(18) -Jsimul(18) -Jsimul(19) Jsimul(28) Jsimul(end) i igen igen*m+i];
%         catch
%             fobiet(i,1)=rand(1,1)*1e6+1e6;
%             fobiet(i,2)=rand(1,1)*1e6+1e6;
%             
%             gc(i,:)=[rand(1,1)*1e6+1e6 -(rand(1,1)*1e6+1e6) rand(1,1)*1e6+1e6 -(rand(1,1)*1e6+1e6) rand(1,1)*1e6+1e6 -(rand(1,1)*1e6+1e6) ...
%                 rand(1,1)*1e6+1e6 -(rand(1,1)*1e6+1e6) rand(1,1)*1e6+1e6 -(rand(1,1)*1e6+1e6) rand(1,1)*1e6+1e6 -(rand(1,1)*1e6+1e6) ...
%                 rand(1,1)*1e6+1e6 -(rand(1,1)*1e6+1e6) i igen igen*m+i];
%             
%             population(igen*m+i).fobiet=fobiet(i,:);
%             population(igen*m+i).gc=gc(i,:);
%         end
    end
elseif m>1 && itorque==1
    unix('rm dataM_*');
    unix('rm resulM_*');
    unix('rm Run_WT_*');
    unix('rm TRAJ*.e*');
    unix('rm TRAJ*.o*');
    for i=1:m
        muX=ya(i,1:ng);
        stdX= [0.0125 0.0125 0.0125 0.02 0.025 0.0333 0.05 0.05 0.05 0.1 0.2 0.2 0.1];%0.03*ones(1,ng);
        skeX=zeros(1,ng);
        kurX=3*ones(1,ng);
        xg=ya(i,ng+1:n);
        save(['dataM_' int2str(i)],'fob','muX','stdX','skeX','kurX','xg','i')
        fdata=['dataM_' int2str(i)];
        fresu=['resulM_' int2str(i)];
        s1=['#!/bin/sh'];
        s2=['#PBS -l nodes=1:ppn=1'];
        s3=['#PBS -q massive4'];
        s4=['#PBS -N TRAJ'];
        s5=['cd $PBS_O_WORKDIR'];
        s6=['/SOFTWARE/matlab2007/bin/matlab -nodisplay -nojvm -nodesktop -r "load(''',fdata,''');'];
        s7=['[muF,varF]=ComputeStat_MC_gaussian(fob,muX,stdX,skeX,kurX,xg,' int2str(i) '); save(''',fresu,''',''muF'',''varF'');exit"'];
        fid=fopen(['Run_WT_' int2str(i) '.sh'],'w');
        fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s%s\n', ...
            s1,s2,s3,s4,s5,s6,s7);
        fclose(fid);
        
        %        s1=['#!/bin/sh'];
        %        s2=['#$ -S /bin/bash'];
        %        s3=['#$ -N TRAJ'];
        %        s4=['#$ -q all.q'];
        %        s5=['#$ -V'];
        %        s6=['#$ -cwd'];
        %        s7=['/home/eminisci/Documents/MATLAB/bin/matlab -nodisplay -nojvm -nodesktop -r "load(''',fdata,''');'];
        %        s8=['[muF,varF]=ComputeStat_URQ(fob,muX,stdX,skeX,kurX,xg); save(''',fresu,''',''muF'',''varF'');exit"'];
        %
        %        fid=fopen(['Run_WT_' int2str(i) '.sh'],'w');
        %        fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n%s%s\n', ...
        %            s1,s2,s3,s4,s5,s6,s7,s8);
        %        fclose(fid);
        unix('qstat -Q massive4 >outmassive.out');
        fid = fopen('outmassive.out');
        Sline1 = fscanf(fid, '%s', 12);
        Sline2 = fscanf(fid, '%s', 12);
        Sline3 = fscanf(fid, '%s', 2);
        Sline3 = fscanf(fid, '%d', 2);
        Sline3 = fscanf(fid, '%s', 2);
        vec = fscanf(fid, '%d', 6);
        Sline3 = fscanf(fid, '%s', 1);
        fclose(fid);
        
        while vec(2,1)>=31
            pause(3)
            unix('qstat -Q massive4 >outmassive.out');
            fid = fopen('outmassive.out');
            Sline1 = fscanf(fid, '%s', 12);
            Sline2 = fscanf(fid, '%s', 12);
            Sline3 = fscanf(fid, '%s', 2);
            Sline3 = fscanf(fid, '%d', 2);
            Sline3 = fscanf(fid, '%s', 2);
            vec = fscanf(fid, '%d', 6);
            Sline3 = fscanf(fid, '%s', 1);
            fclose(fid);
            
        end
        unix(['qsub Run_WT_' int2str(i) '.sh']);
        pause(2)
    end
    icrun=10;
    while icrun>0
        pause(20)
        try
            fprintf(1,'\b%d',icrun);
            for i=1:m
                fprintf(1,'\b%d',i);
                load(['resulM_' int2str(i)])
                if isnan(muF)
                    fobiet(i,1)=1e9+rand(1,1)*1e4;
                    fobiet(i,2)=1e9+rand(1,1)*1e4;
                    
                    gc(i,:)=[fobiet(i,1) fobiet(i,2) fobiet(i,2) fobiet(i,2) i igen];
                    %gc(i,:)=[fobiet(i,1) fobiet(i,2)];
                else
                    fobiet(i,1)=1e6-muF(1,1);
                    fobiet(i,2)=varF(1,1);
                    
                    gc(i,:)=[1e6-muF(1,1)  varF(1,1) muF(1,2) muF(1,3) i igen];
                    %gc(i,:)=[1e6-muF(1,1)  muF(1,2)];
                end
            end
            icrun=0;
        catch
            icrun=10;
        end
    end
end
