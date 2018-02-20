function plot_solution_vs_time(x,u,x_0,x_b,t_0,t_f,els,structure,varargin)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%

if ~isempty(varargin)
   
    plot_nmbr = varargin{1};
    
else
   
    plot_nmbr = 1;
    
end

t_plot = linspace(-1,1,100);    % points per elements, used just to plot

real_els = t_0+els*(t_f-t_0);%t_0+structure.els*(t_f-t_0);

tns = structure.s_nodes;
tnc = structure.c_nodes;

if isnan(tns)
    
    tns = 0;
    
end

if isnan(tnc)
    
    tnc = 0;
    
end


nplots = 1+(sum(structure.num_controls)>0);
figure(plot_nmbr)
subplot(nplots,1,1)
hold on
ColorSet=varycolor(structure.num_eqs);

%% Solution

if strcmp(structure.state_distrib,'linear')
    
    if structure.state_order>0
        
        for i = 1:structure.num_elems
            
            % eval polynomial inside element, 100 points
            xx = (real_els(i,1)+real_els(i,2))/2+(real_els(i,2)-real_els(i,1))/2*t_plot;
            xns = (structure.in_nodes_state(i,1)+structure.in_nodes_state(i,end))/2+(structure.in_nodes_state(i,end)-structure.in_nodes_state(i,1))/2*tns;
            
            for q = 1:structure.num_eqs
                
                % compute function handle
                p = coeffs_to_handle(structure.PC(:,1+(i-1)*(structure.state_order+1):(structure.state_order+1)+(i-1)*(structure.state_order+1),q)');
                
                plot(xx,p(t_plot)'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'Color', ColorSet(q,:));
                
                % extremal nodal values
                
                plot(xns(1),p(tns(1))'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'o','Color', ColorSet(q,:));
                plot(xns(end),p(tns(end))'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'o','Color', ColorSet(q,:));
                
                % internal nodal values
                if length(xns)>1
                    
                    plot(xns(2:end-1),p(tns(2:end-1))'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'*','Color', ColorSet(q,:));
                    
                end
                
            end
            
        end
        
    else
        
        for i = 1:structure.num_elems
            
            % eval polynomial inside element, 100 points
            xx = (real_els(i,1)+real_els(i,2))/2+(real_els(i,2)-real_els(i,1))/2*t_plot;
            xns = (structure.in_nodes_state(i,1)+structure.in_nodes_state(i,end))/2+(structure.in_nodes_state(i,end)-structure.in_nodes_state(i,1))/2*tns;
            
            for q = 1:structure.num_eqs
                
                % nodal solutions
                %plot(els(i,:),sols(i,:,q),'o','Color',ColorSet(q,:));
                
                % compute function handle
                p = coeffs_to_handle(structure.PC(:,1+(i-1)*(structure.state_order+1):(structure.state_order+1)+(i-1)*(structure.state_order+1),q)');
                
                plot(xx,repmat(p(t_plot),length(xx),1)*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'Color', ColorSet(q,:));
                
                % extremal nodal values
                plot(xns(1),p(tns(1))'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'o','Color', ColorSet(q,:));
                plot(xns(end),p(tns(end))'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'o','Color', ColorSet(q,:));
                
                % internal nodal values
                if length(xns)>1
                    
                    plot(xns(2:end-1),p(tns(2:end-1))'*x(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,i),'*','Color', ColorSet(q,:));
                    
                end
                
            end
            
        end
        
    end
    
else
    
    if strcmp(structure.state_distrib,'Cheby')||strcmp(structure.state_distrib,'Lobatto')||strcmp(structure.state_distrib,'Legendre')
        
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
        
    end
    
end

%% Control variables
figure(plot_nmbr)
subplot(nplots,1,nplots)
hold on
ColorSet=varycolor(structure.num_controls);

if strcmp(structure.control_distrib,'linear')
    
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
                plot(xnc(1),p(tnc(1))'*u(1+(structure.control_order+1)*(q-1):(structure.control_order+1)*q,i),'o','Color', ColorSet(q,:));
                plot(xnc(end),p(tnc(end))'*u(1+(structure.control_order+1)*(q-1):(structure.control_order+1)*q,i),'o','Color', ColorSet(q,:));
                
                % internal nodal values
                if length(xnc)>1
                    
                    plot(xnc(2:end-1),p(tnc(2:end-1))'*u(1+(structure.control_order+1)*(q-1):(structure.control_order+1)*q,i),'*','Color', ColorSet(q,:));
                    
                end
                
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
                plot(xnc(1),p(tnc(1))'*u(1+(structure.control_order+1)*(q-1):(structure.control_order+1)*q,i),'o','Color', ColorSet(q,:));
                plot(xnc(end),p(tnc(end))'*u(1+(structure.control_order+1)*(q-1):(structure.control_order+1)*q,i),'o','Color', ColorSet(q,:));
                
                % internal nodal values
                if length(xnc)>1
                    
                    plot(xnc(2:end-1),p(tnc(2:end-1))'*u(1+(structure.control_order+1)*(q-1):(structure.control_order+1)*q,i),'*','Color', ColorSet(q,:));
                    
                end
                
            end
            
        end
        
    end
    
else
    
    if strcmp(structure.control_distrib,'Cheby')||strcmp(structure.control_distrib,'Legendre')||strcmp(structure.control_distrib,'Lobatto')
        
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

end
