function [b_list,D] = MADS_generate_spanning_set (b_list,l,n,deltam,type)

% following exactly Mesh Adaptive Direct Search Algorithms for Constrained
% Optimisation

if l>size(b_list,2)
   
    b_list = [b_list MADS_generate_new_poll_direction(l,n)];
    
end

b = b_list(:,l);
    
% type 1 = minimal (n+1 directions), 2 = maximal (2n directions)

index = 1:length(b);
ihat = index(abs(b)==2^l);

L = zeros(n-1,n-1);

set = -2^l+1:2^l-1;

% generation of lower triangular matrix L

for i = 1:n-1
    
    plus_or_min = randi(2);
    L(i,i) = (-1)^plus_or_min*2^l;
    
    % just a quick way to generate the elements below the diagonal for
    % matrix L
    
    pick = randi(length(set),n-1-i,1);
    L(i+1:end,i) = set(pick);    
    
end

% permutation of lines of matrix L and completion to a basis B in R^n

N_without_ihat = 1:n;
N_without_ihat(ihat) = [];

index = randperm(length(N_without_ihat));
Perms = N_without_ihat(index);

B = zeros(n,n);

for i = 1:n-1
   
    B(Perms(i),1:end-1) = L(i,:);
    
end

B(ihat,1:end-1) = zeros(1,n-1); % should not be needed, as ihat is removed from indexes in Perms...
B(:,end) = b;

% permutation of columns of matrix B

Perms = randperm(n);

Bprime = zeros(n,n);

for i = 1:n-1
   
    Bprime(:,Perms(i)) = B(:,i);
    
end

% completion to a positive basis

if type ==1 
   
    D = [Bprime sum(Bprime,2)];
    %deltap = n*(deltam^0.5);
    
else
    
    D = [Bprime -Bprime];
    %deltap = deltam^0.5;
    
end