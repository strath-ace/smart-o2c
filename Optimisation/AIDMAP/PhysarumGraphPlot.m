function [G,f] = PhysarumGraphPlot(Inputs, ListNodes, History, outputfile)
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
    lastnode(i) = bestindextrack{i}(end);
end


%Sanity check
if (length(fieldnames(ListNodes))==1)
    disp('No chains have been found');
    return
end

%Loop over all the nodes
for i = 2:length(nodenames)
    %Find each node's parent
    treeplotnodes(i) = strmatch(char(ListNodes.(char(nodenames(i))).parent), nodenames, 'exact');
    radius(i-1) = ListNodes.(char(nodenames(i))).radius;
end

AdjMat = zeros(length(treeplotnodes));
for i = 2:length(treeplotnodes);
    s(i-1) = treeplotnodes(i);
    t(i-1) = i;
end


figure('units','normalized','outerposition',[0 0 1 1])
G = graph(s,t,radius);
vidObj = VideoWriter(strcat(outputfile,'.mp4'),'MPEG-4');
vidObj.Quality = 100;
vidObj.FrameRate = 5;
open(vidObj);
set(gcf,'Renderer','zbuffer');

nodeclr = ones(height(G.Edges)+1,3);
edgeclr = ones(height(G.Edges),3);
h = plot(G,'Layout','force','EdgeColor',edgeclr,'NodeColor',nodeclr);
set(h,'NodeLabel',[])
for i = 1:length(G.Edges.EndNodes);  
     [~, loc]=ismember([s(i) t(i)],G.Edges.EndNodes,'rows');
     edgeclr(loc,:) = zeros(1,3);
     nodeclr(s(i),:) = zeros(1,3);        
     nodeclr(t(i),:) = zeros(1,3);
    
    [indexcheck, indexloc] = ismember(i,arrlength);
    [bestindexcheck, bestindexloc] = ismember(i,arrlength(Inputs.NumberOfAgents:Inputs.NumberOfAgents:end));
    if indexcheck
        G = graph(s,t,RadiusHistoryPad{indexloc});   
        G.Edges.LWidths = 6*G.Edges.Weight/max(maxradius);
        h.LineWidth = G.Edges.LWidths;
    end
    
    set(h,'EdgeColor',edgeclr,'NodeColor',nodeclr);
    if bestindexcheck
        highlight(h,bestindextrack{bestindexloc},'NodeColor','g','EdgeColor','g')
    end
        
    drawnow
    writeVideo(vidObj,getframe);
end
close(vidObj)
%Plot the tree
%treeplot(treeplotnodes);

% %Find the x & y locations of each node
% [x,y] = treelayout(treeplotnodes);
% x = x';
% y = y';
% 
% %Obtain the current axis size
% ax = axis;
% 
% %Set the node_ID as label for each node
% txt = text(x(:,1), y(:,1), nodenames, 'FontSize', 5/(ax(4)-ax(3)), 'VerticalAlignment', 'bottom',...
%     'HorizontalAlignment', 'right', 'Interpreter', 'none');
% 
% %Remove labels & axis ticks
% set(gca, 'XTick', [])
% set(gca, 'YTick', [])
% xlabel('')
% ylabel('')
% 
% h = zoom; % get handle to zoom utility
% set(h,'ActionPostCallback',@zoomCallBack);
% set(h,'Enable','on');
% %Every time ones zooms in, this function is executed
% function zoomCallBack(~, evd)      
%     % Since it's expect that ax(4)-ax(3) gets smaller when one zooms in, the
%     % fontsize gets bigger.
%     ax = axis(evd.Axes); % get axis size
%     % change font size accordingly      
%     set(txt,'FontSize',5/(ax(4)-ax(3))); 
% end

end

