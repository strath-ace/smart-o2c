    div
    sum(funvin==0)
    numind=numind0;
    x=[];
            aes0=aes0+.1;
            aes=aes0;
            %aes0=min(aes0,200);
            %aes=round(rand*aes0)
            ft=fit';
            ft=ft.^aes;
            %ft=ft./sum(ft)*size(y,1);
            numsamp=numind;
            val=val+numsamp;
            if fix(igen/5)==igen/5
                    x=parzenself_k(ft,y,fobiet,numsamp,'t1',0);
            else
                    x=parzenself_k(ft,y,fobiet,numsamp,'norm',0);
            end 
            x=VerifyBounds(x,limiti);
    %save (strcat(int2str(igen+1),'.dat'),'x','-ascii','-double')
    [fobietl,gcl]=calcola(Fname,x,igen);
    [fobietx,funvinx,funvx,cdisx]=constraints(pesi,div,fobietl,gcl);
    [fobiety,funviny,funvy,cdisy]=constraints(pesi,div,fobiet,gccurr);
    
    %save (strcat(int2str(igen+1),'c.dat'),'cdisx','-ascii','-double')
    y=[y;x];
%    y=x;
    fobietf=[fobiety;fobietx];
    funvinf=[funviny,funvinx];
    cdis=[cdisy;cdisx];
    GC=[gccurr; gcl];
    x=[];
    [fobiet,x,fit,psc,funvin,cdis,GC]=classifica_nsga2(y,fobietf,funvinf,cdis,GC,Ncl,pam);
    
    y=x;
    x=[];
    
    iord=(1:1:size(y,1));

    
%   
%   ***************************************************************************************
%    for ind=1:size(y,1)
%        if iord(ind)<size(y,1)/2.
%            for jc=1:(size(cdis,2))
%                if cdis(ind,jc)>0
%                    pesi(jc)=pesi(jc)+(1./size(y,1)/12.);
%                end
%            end
%        end
%    end
    pesi
%   ***************************************************************************************
    igen=igen+1;
    genvec=[genvec;igen];
    s=sum(psc==1);
       if s>numindp
            nbest=[];
            nfobest=[];
            npscbest=[];
            nfitbest=[];
            nfunvbest=[];
            ncbest=[];
            ngcbest=[];
            npsc=0;
            for ipcl=1:size(y,1)
                if iord(ipcl)<numindp+1
                    nbest=[nbest;y(ipcl,:)];
                    nfobest=[nfobest;fobiet(ipcl,:)];
                    npscbest=[npscbest;psc(ipcl,1)];
                    nfitbest=[nfitbest;fit(ipcl,1)];
                    nfunvbest=[nfunvbest,funvin(1,ipcl)];
                    ncbest=[ncbest;cdis(ipcl,:)];
                    ngcbest=[ngcbest;GC(ipcl,:)];
                end
            end
            y=nbest;
            fobiet=nfobest;
            psc=npscbest;
            fit=nfitbest;
            funvin=nfunvbest;
            cdis=ncbest;
            gccurr=ngcbest;
        else
            if size(y,1)>numindp
                nbest=[];
            nfobest=[];
            npscbest=[];
            nfitbest=[];
            nfunvbest=[];
            ncbest=[];
            ngcbest=[];
            npsc=0;
            for ipcl=1:size(y,1)
                if iord(ipcl)<numindp+1
                    nbest=[nbest;y(ipcl,:)];
                    nfobest=[nfobest;fobiet(ipcl,:)];
                    npscbest=[npscbest;psc(ipcl,1)];
                    nfitbest=[nfitbest;fit(ipcl,1)];
                    nfunvbest=[nfunvbest,funvin(1,ipcl)];
                    ncbest=[ncbest;cdis(ipcl,:)];
                    ngcbest=[ngcbest;GC(ipcl,:)];
                end
            end
            y=nbest;
            fobiet=nfobest;
            psc=npscbest;
            fit=nfitbest;
            funvin=nfunvbest;
            cdis=ncbest;
            gccurr=ngcbest;
            end
        end
    fobietf=fobiet;
    funvinf=funvin;
    [fobiet,x,fit,psc,funvin,cdis,gccurr]=classifica_nsga2(y,fobietf,funvinf,cdis,gccurr,Ncl,pam);
    y=x;
    x=[];
    %[fit]=order(y,fobiet,psc);

    iord=(1:1:size(y,1));
    % save(['current_pop_',int2str(igen)],'y','fobiet','cdis','iord','igen')
%% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
%    if fix(igen/5)==igen/5
%        iupg=iupg+1;
%        upgrade_DB_level1(y,fobiet,cdis,iupg,Ncl,pam)
%        upgrade_DB_level0_CFD(iupg)
%        unix('date')
%        for iti=1:195
%            pause(300)
%            unix('date')
%        end
%        unix('python del_all_CFD.py')
%        pause(20)
%        upgrade_DB_level0_CFD_read(iupg)
%        create_ann_HTV2_iupg_CL(iupg,1)
%
%        [fobietl,gcl]=calcola(Fname,y,igen);
%        [fobietf,funvinf,funvf,cdis]=constraints(pesi,div,fobietl,gcl);
%        [fobiet,x,fit,psc,funvin,cdis,gcl]=classifica_nsga2(y,fobietf,funvinf,cdis,gcl,Ncl,pam);
%        
%        y=x;
%        x=[]; 
%        iord=(1:1:size(y,1));
%    end
%% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
    best=[];
    fobest=[];
    pscbest=[];
    fitbest=[];
    funvbest=[];
    for ipcl=1:size(y,1)
            if funvin(ipcl)==0
            best=[best;y(ipcl,:)];
            fobest=[fobest;fobiet(ipcl,:)];
            pscbest=[pscbest;psc(ipcl,1)];
            fitbest=[fitbest;fit(ipcl,1)];
            funvbest=[funvbest,funvin(1,ipcl)];
        end
    end
    numbest=size(best,1);
    if numbest<=1
        for ipcl=1:size(y,1)
            if iord(ipcl)<=2
                best=[best;y(ipcl,:)];
                fobest=[fobest;fobiet(ipcl,:)];
                pscbest=[pscbest;psc(ipcl,1)];
                fitbest=[fitbest;fit(ipcl,1)];
                funvbest=[funvbest,funvin(1,ipcl)];
            end
        end
    end
    valvec=[valvec;val];
    max(psc);
    funvec=[funvec;min(fobiet),mean(fobiet)];
    vinvec=[vinvec;min(funvin'),max(funvin'),mean(funvin')];
    vinvec(end,:)
    cdisvec=[cdisvec;mean(cdis)];
    cdisvec(end,:)
    mean(cdis);
%     save(['y_corr' num2str(igen) '.dat'], '-double', '-ascii', 'y');
%     save(['f_corr' num2str(igen) '.dat'], '-double', '-ascii', 'fobiet');
%     save(['c_corr' num2str(igen) '.dat'], '-double', '-ascii', 'cdis');
%     save(['g_corr' num2str(igen) '.dat'], '-double', '-ascii', 'gccurr');
    disp([igen min(fobiet(:,1)) min(fobiet(:,2))])
%     if fix(igen/20)==igen/20
%         load('../data_all_caso_1')
%         close all
%         figure(1)
%         plot(data_front(:,1),data_front(:,2),'g*')
%         hold on
%         
%         if size(fobiet,2)<=2
%             plot(fobiet(:,1),fobiet(:,2),'*',fobest(:,1),fobest(:,2),'*r'),title(strcat('Pareto Front after',' ',int2str(igen),' generations',' Population',int2str(size(y,1)),' Best',int2str(size(best,1)))),xlabel('F1'),ylabel('F2');
%         elseif size(fobiet,2)==3
%             plot3(fobiet(:,1),fobiet(:,2),fobiet(:,3),'*',fobest(:,1),fobest(:,2),fobest(:,2),'*r'),title(strcat('Pareto Front after',' ',int2str(igen),' generations',' Population',int2str(size(y,1)),' Best',int2str(size(best,1)))),xlabel('F1'),ylabel('F2'),zlabel('F3');
%         end
%     end
%     %plot(fobiet(:,1),ones(length(fobiet),1),'*',fobest(:,1),ones(length(fobest),1),'*r'),title(strcat('Pareto Front after',' ',int2str(igen),' generations',' Population',int2str(size(y,1)),' Best',int2str(size(best,1)))),xlabel('F1'),ylabel('F2');
%     figure(2)
%     %subplot(2,3,1)
%     plot(genvec(:,1),vinvec(:,1),genvec(:,1),vinvec(:,2),genvec(:,1),vinvec(:,3)),legend('Vmin','Vmax','Vmed');
%     subplot(2,3,2)
%     plot(genvec(:,1),cdisvec(:,1),genvec(:,1),cdisvec(:,2)),legend('Cmed1','Cmed2');
%     subplot(2,3,3)
%     plot(genvec(:,1),cdisvec(:,3),genvec(:,1),cdisvec(:,4)),legend('Cmed3','Cmed4');
%     subplot(2,3,4)
%     plot(genvec(:,1),cdisvec(:,5),genvec(:,1),cdisvec(:,6)),legend('Cmed5','Cmed6');
%     subplot(2,3,5)
%     plot(genvec(:,1),cdisvec(:,7),genvec(:,1),cdisvec(:,8)),legend('Cmed7','Cmed8');
%     subplot(2,3,6)
%     plot(genvec(:,1),cdisvec(:,9),genvec(:,1),cdisvec(:,10)),legend('Cmed9','Cmed10');
   drawnow;
   pause(3)
