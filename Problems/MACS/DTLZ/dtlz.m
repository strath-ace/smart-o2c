%
% DTLZ.M
%
% Scalable test suite.
%
% f = dtlz(x, M, k, testNo)
%
% Inputs: x      - candidate solutions
%         M      - no. of objectives
%         k      - no. of variables related to G(M) 
%         testNo - test problem no.
%
% Output: f      - results
%
% Ref: Deb, K., Thiele, L., Laumanns, M., and Zitzler, E., 2001,
%      'Scalable Test Problems for Evolutionary Multi-Objective
%      Optimization', TIK-Technical Report No. 112, ETH Zurich,
%      Switzerland.
%
% Created by Robin Purshouse, 17-Oct-2001
% Add inverted problems
% testNo == 11  % Inverted DTLZ1
% testNo == 12  % Inverted DTLZ2
% testNo == 22  % DTLZ2 + 5
% testNo == 24  % DTLZ2 + 5
% Updatd by Rui, 17-Oct-2014
%% Add constrained problems
% Ref: H. Jain and K. Deb, "An evolutionary many-objective optimization algorithm 
% using reference-point based nondominated sorting approach, part II: handling
% constraints and extending to an adaptive approach," Evolutionary Computation,
% IEEE Transactions on, vol. 18, pp. 602-622, 2014.

% testNo == 91  % Type-1 C-DTLZ1  barrier to approach the PF
% testNo == 92  % Type-1 C-DTLZ3  barrier to approach the PF
% testNo == 93  % Type-2 C-DTLZ2  discontinuities
% testNo == 94  % Type-2 C-DTLZ2  discontinuities
% testNo == 95  % Type-3 C-DTLZ1  PF = portion of contraints
% testNo == 96  % Type-2 C-DTLZ2  PF = portion of contraints
% By Rui, 19-Apr-2016

%
function f = dtlz(x, M,k,testNo)

% Check for correct number of inputs.
if nargin ~= 4
   error('Four inputs are required.');
end

% Get total number of decision variables and no. of candidates.
[noSols, n] = size(x);

% Check for consistency.
if( n ~= (M + k -1) )
   error('Input data is inconsistent.');
end

% Could also check that input data is in range [0, 1].

% Set-up the output matrix.
f = NaN * ones(noSols, M);

% Select the correct test problem.
if testNo == 1
  % Compute 'g' functional.
  if k > 0
    g = 100 * (k + sum( (x(:,M:n) - 0.5).^2 - cos(20*pi*(x(:,M:n)-0.5)), 2) );
  else
    g = 0;
  end
  % Compute the objectives.
  f(:,1) = 0.5 * prod( x(:,1:M-1), 2) .* (1 + g); 
  for(fNo = 2:M-1)
    f(:,fNo) = 0.5 * prod( x(:,1:M-fNo), 2) .* (1 - x(:,M-fNo+1)) .* (1 + g);
  end
  f(:,M) = 0.5 * (1 - x(:,1)) .* (1+g);
elseif testNo == 2
  if k > 0
    g = sum( (x(:,M:n) - 0.5).^2 ,2);
  else
    g = 0;
  end
  f(:,1) = prod( cos(x(:,1:M-1)*pi/2) ,2) .* (1 + g);
  for(fNo = 2:M-1)
    f(:,fNo) = prod( cos(x(:,1:M-fNo)*pi/2) ,2) .* sin( x(:,M-fNo+1)*pi/2 ) .* (1 + g);
  end
  f(:,M) = sin( x(:,1)*pi/2 ) .* (1 + g);
elseif testNo == 3
  if k > 0
    g = 100 * (k + sum( (x(:,M:n) - 0.5).^2 - cos(20*pi*(x(:,M:n)-0.5)), 2) );
  else
    g = 0;
  end
  f(:,1) = prod( cos(x(:,1:M-1)*pi/2) ,2) .* (1 + g);
  for(fNo = 2:M-1)
    f(:,fNo) = prod( cos(x(:,1:M-fNo)*pi/2) ,2) .* sin( x(:,M-fNo+1)*pi/2 ) .* (1 + g);
  end
  f(:,M) = sin( x(:,1)*pi/2 ) .* (1 + g);
elseif testNo == 4
  if k > 0
    g = sum( (x(:,M:n) - 0.5).^2 ,2);
  else
    g = 0;
  end
  % Compute meta-variables.
  xm = x.^100;
  f(:,1) = prod( cos(xm(:,1:M-1)*pi/2) ,2) .* (1 + g);
  for(fNo = 2:M-1)
    f(:,fNo) = prod( cos(xm(:,1:M-fNo)*pi/2) ,2) .* sin( xm(:,M-fNo+1)*pi/2 ) .* (1 + g);
  end
  f(:,M) = sin( xm(:,1)*pi/2 ) .* (1 + g);
elseif testNo == 5
  if k > 0
    g = sum( (x(:,M:n) - 0.5).^2 ,2);
  else
    g = 0;
  end
  % Compute meta-variables.
  xm = [x(:,1)*pi/2, rep(pi ./ (4 * (1+g)), [1, M-2]) .* (1 + 2*rep(g, [1, M-2]).*x(:,2:M-1))];
  f(:,1) = prod( cos(xm(:,1:M-1)) ,2) .* (1 + g);
  for(fNo = 2:M-1)
    f(:,fNo) = prod( cos(xm(:,1:M-fNo)) ,2) .* sin(xm(:,M-fNo+1)) .* (1 + g);
  end
  f(:,M) = sin( xm(:,1)) .* (1 + g);

elseif testNo == 6
  if k > 0
    g = sum( x(:,M:n).^(0.1) ,2);
  else
    g = 0;
  end
  % Compute meta-variables.
  xm = [x(:,1)*pi/2, rep(pi ./ (4 * (1+g)), [1, M-2]) .* (1 + 2*rep(g, [1, M-2]).*x(:,2:M-1))];
  f(:,1) = prod( cos(xm(:,1:M-1)) ,2) .* (1 + g);
  for(fNo = 2:M-1)
    f(:,fNo) = prod( cos(xm(:,1:M-fNo)) ,2) .* sin(xm(:,M-fNo+1)) .* (1 + g);
  end
  f(:,M) = sin( xm(:,1)) .* (1 + g);
  
  
elseif testNo == 7
  if k > 0 
    g = 1 + (9 / k) * sum( x(:,M:n), 2);
  else
    g = 0;
  end
  f(:,1:M-1) = x(:,1:M-1);
  % Compute 'h' functional.
  h = M - sum( ( f(:,1:M-1) ./ (1 + rep(g, [1, M-1])) ) .* (1 + sin(3*pi*f(:,1:M-1)))  ,2);
  f(:,M) = (1+g) .* h;
  
elseif testNo == 11  % Inverted DTLZ1
  % Compute 'g' functional.
  if k > 0
    g = 100 * (k + sum( (x(:,M:n) - 0.5).^2 - cos(20*pi*(x(:,M:n)-0.5)), 2) );
  else
    g = 0;
  end
  % Compute the objectives.
  f(:,1) = 0.5 * prod( x(:,1:M-1), 2) .* (1 + g); 
  for(fNo = 2:M-1)
    f(:,fNo) = 0.5 * prod( x(:,1:M-fNo), 2) .* (1 - x(:,M-fNo+1)) .* (1 + g);
  end
  f(:,M) = 0.5 * (1 - x(:,1)) .* (1+g);
  % transformation
  f = rep(0.5*(1+g),[1,M])-f;
  
elseif testNo == 12
    if k > 0
        g = sum( (x(:,M:n) - 0.5).^2 ,2);
    else
        g = 0;
    end
    f(:,1) = prod( cos(x(:,1:M-1)*pi/2) ,2) .* (1 + g);
    for(fNo = 2:M-1)
        f(:,fNo) = prod( cos(x(:,1:M-fNo)*pi/2) ,2) .* sin( x(:,M-fNo+1)*pi/2 ) .* (1 + g);
    end
    f(:,M) = sin( x(:,1)*pi/2 ) .* (1 + g);
    f = rep((1+g),[1,M])-f;
    % transformation
    %f(:,1:M-1) = rep((1+g),[1,M-1]) - f(:,1:M-1).^8;
    %f(:,M) = (1+g) - f(:,M).^4;
elseif testNo == 22 % Equal to  f(DTLZ2) + 5
    if k > 0
        g = sum( (x(:,M:n) - 0.5).^2 ,2);
    else
        g = 0;
    end
    f(:,1) = prod( cos(x(:,1:M-1)*pi/2) ,2) .* (1 + g);
    for(fNo = 2:M-1)
        f(:,fNo) = prod( cos(x(:,1:M-fNo)*pi/2) ,2) .* sin( x(:,M-fNo+1)*pi/2 ) .* (1 + g);
    end
    f(:,M) = sin( x(:,1)*pi/2 ) .* (1 + g);
    f = f+5;
elseif testNo == 24   % Equal to  f(DTLZ4) + 5
    if k > 0
        g = sum( (x(:,M:n) - 0.5).^2 ,2);
    else
        g = 0;
    end
    % Compute meta-variables.
    xm = x.^100;
    f(:,1) = prod( cos(xm(:,1:M-1)*pi/2) ,2) .* (1 + g);
    for(fNo = 2:M-1)
        f(:,fNo) = prod( cos(xm(:,1:M-fNo)*pi/2) ,2) .* sin( xm(:,M-fNo+1)*pi/2 ) .* (1 + g);
    end
    f(:,M) = sin( xm(:,1)*pi/2 ) .* (1 + g);
    f = f+5;
elseif testNo == 91
    %% DTLZ 1 part
    % Compute 'g' functional.
    if k > 0
        g = 100 * (k + sum( (x(:,M:n) - 0.5).^2 - cos(20*pi*(x(:,M:n)-0.5)), 2) );
    else
        g = 0;
    end
    % Compute the objectives.
    f(:,1) = 0.5 * prod( x(:,1:M-1), 2) .* (1 + g);
    for(fNo = 2:M-1)
        f(:,fNo) = 0.5 * prod( x(:,1:M-fNo), 2) .* (1 - x(:,M-fNo+1)) .* (1 + g);
    end
    f(:,M) = 0.5 * (1 - x(:,1)) .* (1+g);
    %% Add a constraint  C>0
    C = 1-f(:,M)/0.6-sum(f(:,1:M-1),2)/0.5;
    f = [f,-1*C];
elseif testNo == 92
    %% DTLZ3 part
    if k > 0
        g = 100 * (k + sum( (x(:,M:n) - 0.5).^2 - cos(20*pi*(x(:,M:n)-0.5)), 2) );
    else
        g = 0;
    end
    f(:,1) = prod( cos(x(:,1:M-1)*pi/2) ,2) .* (1 + g);
    for(fNo = 2:M-1)
        f(:,fNo) = prod( cos(x(:,1:M-fNo)*pi/2) ,2) .* sin( x(:,M-fNo+1)*pi/2 ) .* (1 + g);
    end
    f(:,M) = sin( x(:,1)*pi/2 ) .* (1 + g);
    %% Add a constraint
    temp1 = sum(f.^2,2)-16;
    if M==3; r = 9;
    elseif M==5; r = 12.5;
    elseif M==8; r = 12.5;
    elseif M==10; r = 15;
    elseif M==15; r = 15;
    else  r = 9;  
    end
    tempr = sum(f.^2,2)-r^2;
    C = temp1.*tempr;
    f = [f,-1*C];
elseif testNo == 93
    %% DTLZ2 part
    if k > 0
        g = sum( (x(:,M:n) - 0.5).^2 ,2);
    else
        g = 0;
    end
    f(:,1) = prod( cos(x(:,1:M-1)*pi/2) ,2) .* (1 + g);
    for(fNo = 2:M-1)
        f(:,fNo) = prod( cos(x(:,1:M-fNo)*pi/2) ,2) .* sin( x(:,M-fNo+1)*pi/2 ) .* (1 + g);
    end
    f(:,M) = sin( x(:,1)*pi/2 ) .* (1 + g);
   %% Add a constaint
   if M==3; r = 0.4;
   else  r = 0.2; end
   
   temp = sum(f.^2,2);
   for ci = 1:M
       C(:,ci) = (f(:,ci)-1).^2+temp-f(:,ci).^2-r.^2;
   end
   P1 = min(C,[],2);
   P2 = sum((f-1/sqrt(M)).^2,2)-r.^2;
   C = -1*min([P1,P2],[],2);
   f = [f,-1*C];
elseif testNo == 94
    %% convex DTLZ2 part
    if k > 0
        g = sum( (x(:,M:n) - 0.5).^2 ,2);
    else
        g = 0;
    end
    f(:,1) = prod( cos(x(:,1:M-1)*pi/2) ,2) .* (1 + g);
    for(fNo = 2:M-1)
        f(:,fNo) = prod( cos(x(:,1:M-fNo)*pi/2) ,2) .* sin( x(:,M-fNo+1)*pi/2 ) .* (1 + g);
    end
    f(:,M) = sin( x(:,1)*pi/2 ) .* (1 + g);
    f(:,1:M-1) = f(:,1:M-1).^4;
    f(:,M) = f(:,M).^2;
   %% Add a constaint
   lambdaC = sum(f,2)/M;
   if M==3; r = 0.225;
   elseif M==5; r = 0.225;
   elseif M==8; r = 0.26;
   elseif M==10; r = 0.26;
   elseif M==15; r = 0.27;
   else  r = 0.2;
   end
   C = sum((f-rep(lambdaC,[1,M])).^2,2)-r.^2;
   f = [f,-1*C];
elseif testNo == 95 
  %% DTLZ 1 part
    % Compute 'g' functional.
  if k > 0
    g = 100 * (k + sum( (x(:,M:n) - 0.5).^2 - cos(20*pi*(x(:,M:n)-0.5)), 2) );
  else
    g = 0;
  end
  % Compute the objectives.
  f(:,1) = 0.5 * prod( x(:,1:M-1), 2) .* (1 + g); 
  for(fNo = 2:M-1)
    f(:,fNo) = 0.5 * prod( x(:,1:M-fNo), 2) .* (1 - x(:,M-fNo+1)) .* (1 + g);
  end
  f(:,M) = 0.5 * (1 - x(:,1)) .* (1+g);
  %% Add a constraint  M constraints Ci>0
  temp = sum(f,2);
  for ci = 1:M
      C(:,ci) = temp-f(:,ci)+f(:,ci)/0.5-1;
  end
  f = [f,-1*C];
elseif testNo == 96
    %% DTLZ 4 part
    if k > 0
        g = sum( (x(:,M:n) - 0.5).^2 ,2);
    else
        g = 0;
    end
    % Compute meta-variables.
    xm = x.^100+0.00001;
    f(:,1) = prod( cos(xm(:,1:M-1)*pi/2) ,2) .* (1 + g);
    for(fNo = 2:M-1)
        f(:,fNo) = prod( cos(xm(:,1:M-fNo)*pi/2) ,2) .* sin( xm(:,M-fNo+1)*pi/2 ) .* (1 + g);
    end
    f(:,M) = sin( xm(:,1)*pi/2 ) .* (1 + g);
    %% Add a constraint  M constraints Ci>0
    temp = sum(f.^2,2);
    for ci = 1:M
        C(:,ci) = (f(:,ci).^2)/4+(temp-f(:,ci).^2)-1;
    end
    f = 0.5*[f,-1*C];
end
    
end