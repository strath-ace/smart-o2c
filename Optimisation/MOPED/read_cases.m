figure(101)
load('../data_all_caso_1')
plot(data_front(:,1),data_front(:,2),'go')
hold on
for irun=1:100
    load(['./Case_',int2str(irun),'/risultati_ottimizzazione_',int2str(irun)])
    ff=[];
    for iind=1:size(fobiet,1)
        if funvin(iind)<=0
            ff=[ff; fobiet(iind,:)];
        end
    end
    plot(ff(:,1),ff(:,2),'r*')
end

        