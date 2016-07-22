function [] = PhysarumGraphPlot(Inputs, ListNodes, History)
%% PhysarumGraphPlot: This function creates an animation that shows the evolution of the graph
%
%% Inputs:
% * InitialisedInputs : Structure containing the PhysarumSolver inputs
% * ListNodes         : Structure containing the graph
% * History           : The vein radii and agent movement throughout the generations
%
%% Outputs: 
% * 
%
%% Author: Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

%Retrieve the node names
nodenames = fieldnames(ListNodes);

%Loop over all the nodes
for i = 2:length(nodenames)
    
    %Find each node's parent
    treeplotnodes(i) = strmatch(char(ListNodes.(char(nodenames(i))).parent), nodenames, 'exact');
    radius(i-1) = ListNodes.(char(nodenames(i))).radius;
end

%Prepare the variables required by MATLAB's graph() function, being the
%departure and arrival nodes
for i = 2:length(treeplotnodes);
    s(i-1) = treeplotnodes(i);
    t(i-1) = i;
end

%Initialise a figure
figure('units','normalized','outerposition',[0 0 1 1])

%Create the graph
G = graph(s,t);

%Set the node and edge color to black & remove the node labels
nodeclr = ones(height(G.Edges)+1,3);
edgeclr = ones(height(G.Edges),3);
h = plot(G,'Layout','force','EdgeColor',edgeclr,'NodeColor',nodeclr);
set(h,'NodeLabel',[])

%Prepare radiushistory for plotting
RadiusHistoryPad = History.radius;

%Loop over all radius snapshots
for i = 1:length(History.radius)
    
    %Save the radii of each vein in this generation
    RadiusHistoryPad{i} = History.radius{i}(2:end);
    
    %Determine the number of veins this agent has walked
    arrlength(i) = length(RadiusHistoryPad{i});
    
    %Find the maximum radius in this generation
    maxradius(i) = max(RadiusHistoryPad{i});
end

%Find the maximum number of nodes in total
maxlength = max(arrlength);

%Ensure that all cells within the RadiusHistoryPad array have the same
%number of elements
for i = 1:length(RadiusHistoryPad)
    RadiusHistoryPad{i}(end+1:maxlength) = 1e-50;
end

%Retrieve the best chain throughout the generations
BestSolutionHist = History.BestSolution;
for i = 1:length(BestSolutionHist);
    for j = 1:length(BestSolutionHist{i})
        
        %Find the indices of the nodes that are part of the best chain
        [~, bestindextrack{i}(j)] = ismember(BestSolutionHist{i}(j),nodenames);
    end
end


%Sanity check
if (length(fieldnames(ListNodes))==1)
    disp('No chains have been found');
    return
end

%Initialise the video writing and set the number of frames
vidObj = VideoWriter(strcat(char(Inputs.GraphPlotFileName),'.mp4'),'MPEG-4');
vidObj.FrameRate = 10;
open(vidObj);
set(gcf,'Renderer','zbuffer');

%Retrieve the number of generations and agents
[generations, agents] = size(History.AgentMovement);

%Initialise the vectors that will hold the agent movement, departure nodes
%and arrival nodes
AgentMovementVec = [];
depvec = []; arrvec = [];

%Loop over all generations and agents
for i= 1:generations
    for j = 1:agents
	
		%Loop over the movement of the current agent
        for k = 1:length(History.AgentMovement{i,j})
		
			%Retrieve the index in the full list of all nodes, of the node
            %that the agent moved to
            [AgentMovement{i,j}(k), ~] = find(strcmp(nodenames,History.AgentMovement{i,j}(k)));
			
			%Save the movement of the current agent
            AgentMovementVec(end+1) = AgentMovement{i,j}(k);
        end
		
		%Save the indices of the departure nodes
        dep{i,j} = AgentMovement{i,j}(1:end-1);
		
		%Save the indices of the arrival nodes
        arr{i,j} = AgentMovement{i,j}(2:end);        
		
        %Loop over all departure nodes (agent moves)
        for k = 1:length(dep{i,j})  

			%Find the index of each move of the agent within the G.Edges.Endnodes array
            [~,  agentmovementindex{i,j}(k)]=ismember([dep{i,j}(k) arr{i,j}(k)],G.Edges.EndNodes,'rows');
        end

		%Generate the vectors that hold the departure and arrival nodes
        depvec(end+1:end+length(dep{i,j})) = dep{i,j};
        arrvec(end+1:end+length(arr{i,j})) = arr{i,j};
    end    
end

%Set all the node and edge colors to be white
nodeclr = ones(height(G.Edges)+1,3);
edgeclr = ones(height(G.Edges),3);

%Plot the graph & remove the node labels
h = plot(G,'Layout','force','EdgeColor',edgeclr,'NodeColor',nodeclr);
set(h,'NodeLabel',[])
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');

%Add the generation and agent counter
txt3 = text(xlim(1),ylim(2),{char(strcat({' '},'Generation',{' '},num2str(1),{' '})),char(strcat({' '},'Agent',{' '},num2str(1),{' '}))},'HorizontalAlignment','left','VerticalAlignment','top','FontSize',16);                  
 
%Initialise a tracker to track the current radius & movement snapshot being plotted
tracker = 1;
hold on

%Loop over all generations and agents
for i = 1:Inputs.Generations
    for j = 1:Inputs.NumberOfAgents
	
		%Highlight the initial starting point (being dep{i,j}(1)) 
        highlight(h,dep{i,j}(1),'NodeColor','r','EdgeColor','r','MarkerSize',3);   

		%Update the current generation and agent counter
        delete(txt3)
        txt3 = text(xlim(1),ylim(2),{char(strcat({' '},'Generation',{' '},num2str(i),{' '})),char(strcat({' '},'Agent',{' '},num2str(j),{' '}))},'HorizontalAlignment','left','VerticalAlignment','top','FontSize',16);
        
		%Save the frame
		pause(1/vidObj.FrameRate)
        writeVideo(vidObj,getframe);
		
		%Reset the marker size
        set(h,'MarkerSize',2)
		
		%Loop over all the moves performed by the agent
        for k = 1:length(dep{i,j}) 

			%Update the current generation and agent counter 
            delete(txt3);    
            txt3 = text(xlim(1),ylim(2),{char(strcat({' '},'Generation',{' '},num2str(i),{' '})),char(strcat({' '},'Agent',{' '},num2str(j),{' '}))},'HorizontalAlignment','left','VerticalAlignment','top','FontSize',16);      

			%Set the current departure and arrival node of the agent from white to black
            index = agentmovementindex{i,j}(k);
            edgeclr(index,:) = zeros(1,3);
            nodeclr(dep{i,j}(k),:) = zeros(1,3);        
            nodeclr(arr{i,j}(k),:) = zeros(1,3);
            set(h,'EdgeColor',edgeclr,'NodeColor',nodeclr);     

			%Highly the current arrival node
            highlight(h,arr{i,j}(k),'NodeColor','r','EdgeColor','r','MarkerSize',3);
            
			%Save the frame
            pause(1/vidObj.FrameRate)
            writeVideo(vidObj,getframe);
			
			%Update marker sizes and color
            set(h,'MarkerSize',2,'NodeColor',nodeclr,'EdgeColor',edgeclr)
            
        end
        
		%Update the current generation and agent counter 
        delete(txt3)
        txt3 = text(xlim(1),ylim(2),{char(strcat({' '},'Generation',{' '},num2str(i),{' '})),char(strcat({' '},'Agent',{' '},num2str(j),{' '}))},'HorizontalAlignment','left','VerticalAlignment','top','FontSize',16);

		%Once agent has completed moving, show update of vein radii 
        G = graph(s,t,RadiusHistoryPad{tracker});   
        G.Edges.LWidths = G.Edges.Weight/max(maxradius);
        h.LineWidth = G.Edges.LWidths;
        
		%Increment tracker to track the current history snapshot to be shown
		tracker = tracker+1;
        
		%Display frame
        pause(1/vidObj.FrameRate)
        writeVideo(vidObj,getframe);
        delete(txt3)
    end
	
	%Update the current generation and agent counter once all agents of a generation have moved 
    delete(txt3)
    txt3 = text(xlim(1),ylim(2),{char(strcat({' '},'Generation',{' '},num2str(i),{' '})),char(strcat({' '},'Agent',{' '},num2str(j),{' '}))},'HorizontalAlignment','left','VerticalAlignment','top','FontSize',16);
	
	%Highly best chain of the current generation and show growth and evaporation
    highlight(h,bestindextrack{i},'NodeColor','g','EdgeColor','g')
    h.LineWidth = G.Edges.LWidths;            
    
	%Save the frame
	pause(1/vidObj.FrameRate)
    writeVideo(vidObj,getframe);
	
	%Remove the current generation and agent counter to prevent overlapping text 
    delete(txt3)
end

%Update the current generation and agent counter once all generations have been shown
txt3 = text(xlim(1),ylim(2),{char(strcat({' '},'Generation',{' '},num2str(i))),char(strcat({' '},'Agent',{' '},num2str(j)))},'HorizontalAlignment','left','VerticalAlignment','top','FontSize',16);

%Close the video file
close(vidObj)
hold off

end

