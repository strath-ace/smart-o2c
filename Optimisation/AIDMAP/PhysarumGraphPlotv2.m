function [G,f] = PhysarumGraphPlotv2(Inputs, ListNodes, History, outputfile)
%This function plots the tree corresponding to the physarum graph
%
% Inputs:
% * ListNodes       : Structure containing the graph
%
% Outputs: 
% * 
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Retireve the node names
nodenames = fieldnames(ListNodes);



%Loop over all the nodes
for i = 2:length(nodenames)
    %Find each node's parent
    treeplotnodes(i) = strmatch(char(ListNodes.(char(nodenames(i))).parent), nodenames, 'exact');
    radius(i-1) = ListNodes.(char(nodenames(i))).radius;
end

for i = 2:length(treeplotnodes);
    s(i-1) = treeplotnodes(i);
    t(i-1) = i;
end

figure('units','normalized','outerposition',[0 0 1 1])
G = graph(s,t);
nodeclr = ones(height(G.Edges)+1,3);
edgeclr = ones(height(G.Edges),3);
h = plot(G,'Layout','force','EdgeColor',edgeclr,'NodeColor',nodeclr);
set(h,'NodeLabel',[])

%Prepare radiushistory for plotting
RadiusHistoryPad = History.radius;
for i = 1:length(History.radius)
    RadiusHistoryPad{i} = History.radius{i}(2:end);
    arrlength(i) = length(RadiusHistoryPad{i});
    maxradius(i) = max(RadiusHistoryPad{i});
end
maxlength = max(arrlength);

for i = 1:length(RadiusHistoryPad)
    RadiusHistoryPad{i}(end+1:maxlength) = 1e-50;
end

BestSolutionHist = History.BestSolution;
for i = 1:length(BestSolutionHist);
    for j = 1:length(BestSolutionHist{i})
    [temp, bestindextrack{i}(j)] = ismember(BestSolutionHist{i}(j),nodenames);
    end
end


%Sanity check
if (length(fieldnames(ListNodes))==1)
    disp('No chains have been found');
    return
end

vidObj = VideoWriter(strcat(outputfile,'.mp4'),'MPEG-4');
%vidObj.Quality = 100;
vidObj.FrameRate = 5;
open(vidObj);
set(gcf,'Renderer','zbuffer');

[generations, agents] = size(History.AgentMovement);
AgentMovementVec = [];
dep = {}; arr = {}; depvec = []; arrvec = [];
for i= 1:generations
    for j = 1:agents
        for k = 1:length(History.AgentMovement{i,j})
            [AgentMovement{i,j}(k), ~] = find(strcmp(nodenames,History.AgentMovement{i,j}(k)));
            AgentMovementVec(end+1) = AgentMovement{i,j}(k);
        end
        dep{i,j} = AgentMovement{i,j}(1:end-1);
        arr{i,j} = AgentMovement{i,j}(2:end);        
        
        for k = 1:length(dep{i,j})            
            [~,  agentmovementindex{i,j}(k)]=ismember([dep{i,j}(k) arr{i,j}(k)],G.Edges.EndNodes,'rows');
        end
           
        
        depvec(end+1:end+length(dep{i,j})) = dep{i,j};
        arrvec(end+1:end+length(arr{i,j})) = arr{i,j};
    end    
end



for i = 1:length(depvec)
[~, movementindex(i)]=ismember([depvec(i) arrvec(i)],G.Edges.EndNodes,'rows');
end

nodeclr = ones(height(G.Edges)+1,3);
edgeclr = ones(height(G.Edges),3);
h = plot(G,'Layout','force','EdgeColor',edgeclr,'NodeColor',nodeclr);
set(h,'NodeLabel',[])
ylim=get(gca,'ylim');
xlim=get(gca,'xlim')
txt3 = text(xlim(2),ylim(2),{char(strcat('Generation',{' '},num2str(1))),char(strcat('Agent',{' '},num2str(1)))},'HorizontalAlignment','right','VerticalAlignment','top','FontSize',16);                   


for i = 1:length(movementindex)
%Plotting code here. Loop over agent movement -> set color etc.
end

%OR 
tracker = 1;
%txt3 = text(xlim(2),ylim(2),{char(strcat('Generation ',num2str(1))),'Test2'},'HorizontalAlignment','right','VerticalAlignment','top','FontSize',16);
       hold on
for i = 1:Inputs.Generations
    for j = 1:Inputs.NumberOfAgents
            highlight(h,dep{i,j}(1),'NodeColor','r','EdgeColor','r','MarkerSize',3);      
            delete(txt3)
            txt3 = text(xlim(2),ylim(2),{char(strcat('Generation',{' '},num2str(i))),char(strcat('Agent',{' '},num2str(j)))},'HorizontalAlignment','right','VerticalAlignment','top','FontSize',16);        
             pause(0.2)
        writeVideo(vidObj,getframe);
        set(h,'MarkerSize',2)
        for k = 1:length(dep{i,j})            
            delete(txt3);    
            txt3 = text(xlim(2),ylim(2),{char(strcat('Generation',{' '},num2str(i))),char(strcat('Agent',{' '},num2str(j)))},'HorizontalAlignment','right','VerticalAlignment','top','FontSize',16);        
            
            index = agentmovementindex{i,j}(k);
            edgeclr(index,:) = zeros(1,3);
            nodeclr(dep{i,j}(k),:) = zeros(1,3);        
            nodeclr(arr{i,j}(k),:) = zeros(1,3);
            set(h,'EdgeColor',edgeclr,'NodeColor',nodeclr);            
            highlight(h,arr{i,j}(k),'NodeColor','r','EdgeColor','r','MarkerSize',3);
            %drawnow update
        pause(0.2)
        writeVideo(vidObj,getframe);
        set(h,'MarkerSize',2,'NodeColor',nodeclr,'EdgeColor',edgeclr)
            %Show movement -> plot
        end
        delete(txt3)
         txt3 = text(xlim(2),ylim(2),{char(strcat('Generation',{' '},num2str(i))),char(strcat('Agent',{' '},num2str(j)))},'HorizontalAlignment','right','VerticalAlignment','top','FontSize',16);
       
           G = graph(s,t,RadiusHistoryPad{tracker});   
            G.Edges.LWidths = 1*G.Edges.Weight/max(maxradius);
            h.LineWidth = G.Edges.LWidths;
            tracker = tracker+1;
            %drawnow update     
            pause(0.2)
            writeVideo(vidObj,getframe);
            delete(txt3)
        %Show dilation -> plot
    end
    delete(txt3)
    txt3 = text(xlim(2),ylim(2),{char(strcat('Generation',{' '},num2str(i))),char(strcat('Agent',{' '},num2str(j)))},'HorizontalAlignment','right','VerticalAlignment','top','FontSize',16);
    
    highlight(h,bestindextrack{i},'NodeColor','g','EdgeColor','g')
    pause(0.2)
       
            writeVideo(vidObj,getframe);
            delete(txt3)
    %Show GF & evaporation -> plot
end

close(vidObj)
hold off

end

