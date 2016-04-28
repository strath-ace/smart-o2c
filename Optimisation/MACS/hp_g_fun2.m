function g=hp_g_fun2(f,lambda,z,maxs)

% High performance g_fun. Identical to original one but is vectorized,
% sparing a whole for loop. Works on a single subproblem (lambda) for now
%
% INPUT
%       f       :       vector of objective function values
%       lambda  :       weigths vector of current problem
%       z       :       minimas of each objective function
%
% OUTPUT
%       g       :       vector of g_fun values

n_a = size(f,1);                                                            % number of agents

delta = maxs-z;
delta = delta.*(delta~=zeros(size(delta)))+ones(size(delta)).*(delta==zeros(size(delta))); % if there's only one value, max and min are the same, so this avoids a useless division by zero

f = f-repmat(z,size(f,1),1);
f = f./repmat(delta,size(f,1),1);

g = max(repmat(lambda,n_a,1).*abs(f-repmat(z,n_a,1)),[],2)';                % efficient evaluation of g_function

return