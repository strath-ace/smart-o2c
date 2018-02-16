function plot_solution_vs_time(x,u,x_0,x_b,t_0,t_f,els,structure,varargin)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% This function plots the states and the controls over time, on two
% parallel plots. 

% Internal nodals values are marked with a '*', while the value of the 
% interpolation/approximation polynomial is given in a solid line. The
% polynomial approximation/interpolation evaluated at the boundaries of
% each element is marked with a 'o'

if ~isempty(varargin)
    
    plot_nmbr = varargin{1};
    
else
    
    plot_nmbr = 1;
    
end

t_plot = linspace(-1,1,100);    % points per elements, used just to plot

real_els = t_0+els*(t_f-t_0);%t_0+structure.els*(t_f-t_0);

tns = structure.integr_nodes;%structure.s_nodes;
tnc = structure.integr_nodes;%structure.c_nodes;

if isnan(tns)
    
    tns = 0;
    
end

if isnan(tnc)
    
    tnc = 0;
    
end


nplots = 2;%1+(sum(structure.num_controls)>0);
figure(plot_nmbr)
subplot(nplots,1,1)
hold on
ColorSet=varycolor(structure.num_eqs);

%% Solution

if structure.state_order>0
    
    for i = 1:structure.num_elems
        
        % eval polynomial inside element, 100 points
        xx = (real_els(i,1)+real_els(i,2))/2+(real_els(i,2)-real_els(i,1))/2*t_plot;
        xns = (real_els(i,1)+real_els(i,2))/2+(real_els(i,2)-real_els(i,1))/2*tns;
        
        for q = 1:structure.num_eqs
            
            % compute function handle
            p = coeffs_to_handle(structure.PC(:,1+(i-1)*(structure.state_order+1):(structure.state_order+1)+(i-1)*(structure.state_order+1),q)');
            
            plot(xx,p(t_plot)'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'Color', ColorSet(q,:));
            
            % extremal nodal values
            if i>1
                
                plot(xx(1),p(t_plot(1))'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'o','Color', ColorSet(q,:));
                
            end
            
            if i<structure.num_elems
                
                plot(xx(end),p(t_plot(end))'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'o','Color', ColorSet(q,:));
                
            end
            
            % boundary values
            if structure.DFET==0
                
                plot(xx(1),p(t_plot(1))'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'o','Color', ColorSet(q,:));
                plot(xx(end),p(t_plot(end))'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'o','Color', ColorSet(q,:));
                
            else
                
                if i==1
                    
                    plot(xx(1),x_0(q)*(i==1),'o','Color', ColorSet(q,:));
                    
                end
                
                if i==structure.num_elems
                    
                    plot(xx(end),x_b(q)*(i==structure.num_elems),'o','Color', ColorSet(q,:));
                    
                end
                
            end
            
            % internal nodal values
            plot(xns,p(tns)'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'*','Color', ColorSet(q,:));
            
        end
        
    end
    
else
    
    for i = 1:structure.num_elems
        
        % eval polynomial inside element, 100 points
        xx = (real_els(i,1)+real_els(i,2))/2+(real_els(i,2)-real_els(i,1))/2*t_plot;
        xns = (real_els(i,1)+real_els(i,2))/2+(real_els(i,2)-real_els(i,1))/2*tns;
        
        for q = 1:structure.num_eqs
            
            % compute function handle
            p = coeffs_to_handle(structure.PC(:,1+(i-1)*(structure.state_order+1):(structure.state_order+1)+(i-1)*(structure.state_order+1),q)');
            
            plot(xx,repmat(p(t_plot),length(xx),1)*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'Color', ColorSet(q,:));
            
            % extremal nodal values
            if structure.DFET==0
                
                plot(xx(1),p(t_plot(1))'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'o','Color', ColorSet(q,:));
                plot(xx(end),p(t_plot(end))'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'o','Color', ColorSet(q,:));
                
            else
                
                if i==1
                    
                    plot(xx(1),x_0(q)*(i==1),'o','Color', ColorSet(q,:));
                    
                end
                
                if i==structure.num_elems
                    
                    plot(xx(end),x_b(q)*(i==structure.num_elems),'o','Color', ColorSet(q,:));
                    
                end
                
            end
            
            % internal nodal values
            plot(xns,p(tns)'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'*','Color', ColorSet(q,:));
            
        end
        
    end
    
end

%% Control variables
figure(plot_nmbr)
subplot(nplots,1,nplots)
hold on
ColorSet=varycolor(structure.num_controls);

if isempty(structure.control_bounds)
    
    plot([real_els(1,1) real_els(end,2)],[nan nan]);
    
else
    
    if structure.control_order>0
        
        for i = 1:structure.num_elems
            
            % eval polynomial inside element, 100 points
            xx = (real_els(i,1)+real_els(i,2))/2+(real_els(i,2)-real_els(i,1))/2*t_plot;
            xnc = (real_els(i,1)+real_els(i,2))/2+(real_els(i,2)-real_els(i,1))/2*tnc;
            
            for q = 1:structure.num_controls
                
                % compute function handle
                p = coeffs_to_handle(structure.PCU(:,1+(i-1)*(structure.control_order+1):(structure.control_order+1)+(i-1)*(structure.control_order+1),q)');
                
                plot(xx,p(t_plot)'*u(1+(structure.control_order+1)*(q-1):(structure.control_order+1)*q,i),'Color', ColorSet(q,:));
                
                % extremal nodal values
                plot(xx(1),p(t_plot(1))'*u(1+(structure.control_order+1)*(q-1):(structure.control_order+1)*q,i),'o','Color', ColorSet(q,:));
                plot(xx(end),p(t_plot(end))'*u(1+(structure.control_order+1)*(q-1):(structure.control_order+1)*q,i),'o','Color', ColorSet(q,:));
                
                % internal nodal values
                plot(xnc,p(tnc)'*u(1+(structure.control_order+1)*(q-1):(structure.control_order+1)*q,i),'*','Color', ColorSet(q,:));
                
            end
            
        end
        
    else
        
        for i = 1:structure.num_elems
            
            % eval polynomial inside element, 100 points
            xx = (real_els(i,1)+real_els(i,2))/2+(real_els(i,2)-real_els(i,1))/2*t_plot;
            xnc = (real_els(i,1)+real_els(i,2))/2+(real_els(i,2)-real_els(i,1))/2*tnc;
            
            for q = 1:structure.num_controls
                
                % compute function handle
                p = coeffs_to_handle(structure.PCU(:,1+(i-1)*(structure.control_order+1):(structure.control_order+1)+(i-1)*(structure.control_order+1),q)');
                
                plot(xx,repmat(p(t_plot),length(xx),1)*u(1+(structure.control_order+1)*(q-1):(structure.control_order+1)*q,i),'Color', ColorSet(q,:));
                
                % extremal nodal values
                plot(xx(1),p(t_plot(1))'*u(1+(structure.control_order+1)*(q-1):(structure.control_order+1)*q,i),'o','Color', ColorSet(q,:));
                plot(xx(end),p(t_plot(end))'*u(1+(structure.control_order+1)*(q-1):(structure.control_order+1)*q,i),'o','Color', ColorSet(q,:));
                
                % internal nodal values
                plot(xnc,p(tnc)'*u(1+(structure.control_order+1)*(q-1):(structure.control_order+1)*q,i),'*','Color', ColorSet(q,:));
                
                
            end
            
        end
        
    end
    
end

end