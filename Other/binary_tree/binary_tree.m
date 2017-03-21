nlevels=30;
max_cost=0.7;
tree(1).branch=[];
tree(1).cost=[];

nbranch=1;
for i=1:nlevels
    k=0;
    for j=1:nbranch
        k=k+1;                
        newtree(k).branch=[tree(j).branch;1];
        newtree(k).cost=tree(j).cost;        
        k=k+1;
        newtree(k).branch=[tree(j).branch;0];        
        newtree(k).cost=tree(j).cost;                
    end
    newnbranch=k;
    m=0;
    clear tree
    tree(1).branch=[];
    tree(1).cost=[];
    %newnbranch
    for j=1:newnbranch
        newtree(j).cost=rand; % replace with actual cost
       % [newtree(j).cost max_cost]
        if newtree(j).cost<max_cost
            m=m+1;
            tree(m).branch=newtree(j).branch;
            tree(m).cost=newtree(j).cost;
        end
      %  [j m]
      %  pause
    end
    nbranch=m;
end