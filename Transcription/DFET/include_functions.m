function structure = include_functions(structure,f,dfx,dfu,g,weights,dgu0,dgxf,dguf,dgxi,dgui,c,dcx,dcu,e,dex,deu,h,wh,dhu0,dhxf,dhuf,dhxi,dhui,q,wq,dqu0,dqxf,dquf,dqxi,dqui,constants)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Includes dynamics, objective functions, constraints, and weights. Where
% derivatives are not supplied, it will generate a function handle to its 
% automatically generated finite difference approximation.

%% Step size is hardcoded here

stepsize = 1e-7;

%% commonly used strings
one_point_args = '(x,u,time,static,scales,constants)';
two_points_args = '(x0,u0,t0,xf,uf,tf,x,u,time,static,scales,constants)';

%% Dynamics

if ~isempty(f)
    
    if ~ischar(f)
        
        error('f must be a string containing the name of the function for the dynamics');
        
    else
        
        fhand = str2func(['@' one_point_args f one_point_args]);
        
    end
    
    % generate jacobian wrt x, if not explicitly given
    
    if isempty(dfx)
        
        dfx = @(x,u,t,static,scales,constants)generate_fd1_deriv_x(x,u,t,static,scales,constants,fhand,stepsize);
        
    end
    
    % generate jacobian wrt u, if not explicitly given
    
    if isempty(dfu)
        
        dfu = @(x,u,t,static,scales,constants)generate_fd1_deriv_u(x,u,t,static,scales,constants,fhand,stepsize);
        
    end
    
else
    
    error('Dynamics not defined, include a function name');
    
end

%% Objective function

if ~isempty(g)
    
    if ~ischar(g)
        
        error('g must be a string containing the name of the function for the dynamics');
        
    else
        
        funname = g;
        ghand = str2func(['@' two_points_args g two_points_args]);
        
    end
    
    % generate jacobian wrt u0, if not explicitly given
    
    if isempty(dgu0)
        
        dgu0 = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_u0(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,ghand,stepsize);
        
    end
    
    % generate jacobian wrt xf, if not explicitly given
    
    if isempty(dgxf)
        
        dgxf = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_xf(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,ghand,stepsize);
        
    end
    
    % generate jacobian wrt uf, if not explicitly given
    
    if isempty(dguf)
        
        dguf = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_uf(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,ghand,stepsize);
        
    end
    
    % generate jacobian wrt x, if not explicitly given
    
    if isempty(dgxi)
        
        dgxi = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_xi(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,ghand,stepsize);
        
    end
    
    % generate jacobian wrt u, if not explicitly given
    
    if isempty(dgui)
        
        dgui = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_ui(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,ghand,stepsize);
        
    end
    
else
    
    ghand = [];
    
end

%% Inequality path constraints

if ~isempty(c)
    
    if ~ischar(c)
        
        error('c must be a string containing the name of the function for the dynamics');
        
    else
        
        chand = str2func(['@' one_point_args c one_point_args]);
        
    end
    
    % generate jacobian wrt x, if not explicitly given
    
    if isempty(dcx)
        
        dcx = @(x,u,t,static,scales,constants)generate_fd1_deriv_x(x,u,t,static,scales,constants,chand,stepsize);
        
    end
    
    % generate jacobian wrt u, if not explicitly given
    
    if isempty(dcu)
        
        dcu = @(x,u,t,static,scales,constants)generate_fd1_deriv_u(x,u,t,static,scales,constants,chand,stepsize);
        
    end
    
else
    
    chand = [];
    
end

%% Equality path constraints

if ~isempty(e)
    
    if ~ischar(e)
        
        error('e must be a string containing the name of the function for the dynamics');
        
    else
        
        ehand = str2func(['@' one_point_args e one_point_args]);
        
    end
    
    % generate jacobian wrt x, if not explicitly given
    
    if isempty(dex)
        
        dex = @(x,u,t,static,scales,constants)generate_fd1_deriv_x(x,u,t,static,scales,constants,ehand,stepsize);
        
    end
    
    % generate jacobian wrt u, if not explicitly given
    
    if isempty(deu)
        
        deu = @(x,u,t,static,scales,constants)generate_fd1_deriv_u(x,u,t,static,scales,constants,ehand,stepsize);
        
    end
    
else
    
    ehand = [];
    
end

%% Boundary and integral inequality constraints

if ~isempty(h)
    
    if ~ischar(h)
        
        error('h must be a string containing the name of the function for the dynamics');
        
    else
        
        hhand = str2func(['@' two_points_args h two_points_args]);
        
    end
    
    % generate jacobian wrt u0, if not explicitly given
    
    if isempty(dhu0)
        
        dhu0 = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_u0(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,hhand,stepsize);
        
    end
    
    % generate jacobian wrt xf, if not explicitly given
    
    if isempty(dhxf)
        
        dhxf = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_xf(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,hhand,stepsize);
        
    end
    
    % generate jacobian wrt uf, if not explicitly given
    
    if isempty(dhuf)
        
        dhuf = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_uf(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,hhand,stepsize);
        
    end
    
    % generate jacobian wrt x, if not explicitly given
    
    if isempty(dhxi)
        
        dhxi = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_xi(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,hhand,stepsize);
        
    end
    
    % generate jacobian wrt u, if not explicitly given
    
    if isempty(dhui)
        
        dhui = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_ui(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,hhand,stepsize);
        
    end
    
else
    
    hhand = [];
    
end

%% Boundary and integral equality constraints

if ~isempty(q)
    
    if ~ischar(q)
        
        error('q must be a string containing the name of the function for the dynamics');
        
    else
        
        qhand = str2func(['@' two_points_args q two_points_args]);
        
    end
    
    % generate jacobian wrt u0, if not explicitly given
    
    if isempty(dqu0)
        
        dqu0 = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_u0(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,qhand,stepsize);
        
    end
    
    % generate jacobian wrt xf, if not explicitly given
    
    if isempty(dqxf)
        
        dqxf = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_xf(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,qhand,stepsize);
        
    end
    
    % generate jacobian wrt uf, if not explicitly given
    
    if isempty(dquf)
        
        dquf = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_uf(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,qhand,stepsize);
        
    end
    
    % generate jacobian wrt x, if not explicitly given
    
    if isempty(dqxi)
        
        dqxi = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_xi(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,qhand,stepsize);
        
    end
    
    % generate jacobian wrt u, if not explicitly given
    
    if isempty(dqui)
        
        dqui = @(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants)generate_fd1_deriv_ui(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,qhand,stepsize);
        
    end
    
else
   
    qhand = [];
   
end

%% include functions in structure

% Dynamics

structure.f = fhand;
structure.dfx = dfx;
structure.dfu = dfu;

% Objectives

structure.g = ghand;
structure.dgu0 = dgu0;
structure.dgxf = dgxf;
structure.dguf = dguf;
structure.dgxi = dgxi;
structure.dgui = dgui;
structure.weights = weights;

% Path constraints

structure.c = chand;
structure.dcx = dcx;
structure.dcu = dcu;

structure.e = ehand;
structure.dex = dex;
structure.deu = deu;

% Final state and integral inequality constraints

structure.h = hhand;
structure.dhu0 = dhu0;
structure.dhxf = dhxf;
structure.dhuf = dhuf;
structure.dhxi = dhxi;
structure.dhui = dhui;
structure.wh = wh;

% Final state and integral equality constraints

structure.q = qhand;
structure.dqu0 = dqu0;
structure.dqxf = dqxf;
structure.dquf = dquf;
structure.dqxi = dqxi;
structure.dqui = dqui;
structure.wq = wq;

end