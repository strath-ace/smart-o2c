function [r] = PlotTrajectories(Solutions, ListNodes, vvecsize, sizefont, dt)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
%% PlotTrajectories: This function is used to plot the trajectories found with the AIDMAP algorithm
% 
%% Inputs:
% * Solutions     : Array containing the unique IDs of the nodes to be
%                   plotted. Can be cell array of multiple solutions
% * ListNodes     : Structure containing (at least) all the nodes to be
%                   plotted
% * vvecsize      : The scaling factor for the thrust vector size
% * sizefont      : Font size
% * dt            : Time-step for the plotting
% 
%% Outputs: 
% * r             : The Cartesian elements over time [cell array]
% 
%% Author(s): Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

if iscellstr(Solutions)
    Solutions = {Solutions};
end

% Loop over all the solutions to be plotted
for i = 1:length(Solutions)
        
    % Initialise the Cartesian coordinates vector
    r{i} = [];
    totalcost = 0;
    
    % Create a figure and set the axis and grid
    figure('Color', [1 1 1])
    axis equal; grid on; hold on
    
    
    
    % Loop over all the nodes in the solution
    for j = 2:length(Solutions{i})
        
        % For ease of reading, define the name of the arrival & departure
        % node
        arrnode = ListNodes.(char(Solutions{i}(j)));
        depnode = ListNodes.(arrnode.parent);
        
        % Retrieve the transfer orbit
        kep_trans = depnode.attributes.kep_trans;
        
        % Generate the time vector
        tvec = depnode.attributes.t_dep:dt:arrnode.attributes.t_dep;
        
        % Retrieve the dV vector
        vvec = (arrnode.attributes.dV_dep.*vvecsize);
        
        % Obtain the current size of the r vector
        [a, ~] = size(r{i});
        
        % add the Cartesian elements to the end of the vector
        for k = (1+a):(length(tvec)+a)
            r{i}(k, :) = StardustTool.CartesianElementsAt(kep_trans, tvec(k-a));
        end
        
        % Plot the arrival node
        h = plot3(depnode.attributes.r_arr(1), depnode.attributes.r_arr(2), depnode.attributes.r_arr(3), 'ko');
        set(h, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k')
        
        % Plot the dV vector
        quiver3(arrnode.attributes.r_dep(1), arrnode.attributes.r_dep(2), arrnode.attributes.r_dep(3), vvec(1), vvec(2), vvec(3), 'g', 'linewidth', 1.5, 'MaxHeadSize', 1);
        
        % Retrieve the asteroid/planet name
        temp = strsplit(depnode.node_ID, '___');
        temp = strsplit(char(temp(end)), '__');
        depname = temp(1);
        text(depnode.attributes.r_arr(1), depnode.attributes.r_arr(2), depnode.attributes.r_arr(3), depname, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'FontSize', sizefont);
        
        % Generate the text for the dV vector
        currentcost = ListNodes.(char(Solutions{i}(j))).length;
        totalcost = totalcost+currentcost;
        dVtext = strcat(num2str(currentcost), ' km/s');
        text(arrnode.attributes.r_dep(1), arrnode.attributes.r_dep(2), arrnode.attributes.r_dep(3), dVtext, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'FontSize', sizefont);
        
        % Check if the last node is being evaluated. 
        % If so, perform a number of additional calculations    
        if j == length(Solutions{i})
            
            % Retrieve the transfer orbit & time vector
            kep_trans = arrnode.attributes.kep_trans;
            tvec = arrnode.attributes.t_dep:dt:arrnode.attributes.t_arr;
            
            % Obtain the current size of the r vector and add the last arc
            % to it
            [a, ~] = size(r{i});
            for k = (1+a):(length(tvec)+a)
                r{i}(k, :) = StardustTool.CartesianElementsAt(kep_trans, tvec(k-a));
            end
            
            % Plot the final asteroid/planet
            h = plot3(arrnode.attributes.r_arr(1), arrnode.attributes.r_arr(2), arrnode.attributes.r_arr(3), 'ko');
            set(h, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
            
            % Retrieve the final asteroid/planet's name
            temp = strsplit(arrnode.node_ID, '___');
            temp = strsplit(char(temp(end)), '__');
            arrname = temp(1);
             
            text(arrnode.attributes.r_arr(1), arrnode.attributes.r_arr(2), arrnode.attributes.r_arr(3), arrname, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'FontSize', sizefont);
        end
            
    
    end
    
    % Set the font size
    set(gca, 'FontSize', sizefont);
    
    % Once the entire trajectory has been calculated, plot the r vector
    plot3(r{i}(:, 1), r{i}(:, 2), r{i}(:, 3), 'b')
    
    % Add axis labels
    xlabel('X [km]')
    ylabel('Y [km]')
    
    % Add the total dV & number of asteroids to plot
    ylim=get(gca, 'ylim');
    xlim=get(gca, 'xlim');
    asteroidnum = length(Solutions{i})-1; % -1 for root
    text(xlim(2), ylim(2), {char(strcat('Number of Asteroids:', {' '}, num2str(asteroidnum))), char(strcat('Total \DeltaV:', {' '}, num2str(totalcost), ' km/s'))}, 'HorizontalAlignment', 'right', 'FontSize', 16);
    

    
end

end

