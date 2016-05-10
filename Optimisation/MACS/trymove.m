function [discarded,impr] = trymove (xtrial,f,discarded,cp,cid,func,arg)

impr = 0;

if cp==0
    
    ftrial=func(xtrial,arg{:});
    maxC=0;
    
    if any(ftrial<f)                                  % if any component of the new trial is better than the same component of the previous x, we have an improvement!
        
        discarded.f=[discarded.f; ftrial];           % put the values computed with x_trial in the discarded structure (why???)
        discarded.x=[discarded.x; xtrial];
        discarded.c=[discarded.c; maxC];
        impr=1;                                                 % and say an improvement has been made, terinating the cycle
        
        % IDEA, BETTER TO CONTINE THE LOOP THAN CHECK THIS VAR ONWARDS, SAVES TIE AND READABILITY
        
    end
    
else
    
    [ftrial,ctrial]=func(xtrial,arg{:});
    maxC=max(ctrial);
    
    if cid<=0&&maxC>0                                     % need to understand what cid and maxC are
        
        ftrial=f+maxC;
        
    end
    
    if (any(ftrial<f)&&maxC<=0)||(cid>0&&maxC<cid)	% as above: if at least 1 obective improves and respects all constraints OR at least one constraint improves from the previous iteration, keep the agent (i.e, discard it...)
        
        discarded.f=[discarded.f; ftrial];
        discarded.x=[discarded.x; xtrial];
        discarded.c=[discarded.c; maxC];
        impr=1;
        
    end
    
end