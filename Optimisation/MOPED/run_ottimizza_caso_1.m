function run_ottimizza_caso_1(iind)
addpath('../')
numindp0=100;
numind0=100;
numop=600;
D=25;
format short g

[risul,fobiet,best,funvin]=ottimizza_Caso(numindp0,numind0,numop,D);

save(['risultati_ottimizzazione_',ind2str(iind)],'risul','fobiet','best','funvin')
    