function [y, dy, update] = rbf_surrogate(x, surrogate, ybest, subobj, minmax)

model = surrogate.model;
X = model.X;
basis = model.basis;
s = model.s;

N = size(X,1);
[n,m] = size(x);

r = zeros(N,n);
D = zeros(N,m);

rbf_func = strcat('rbf_',lower(basis));
rbf_func = str2func(rbf_func);

for i = 1:N
    for j = 1:n
        D(i,:) = x(j,:) - X(i,:);
        r(i,j) = norm(D(i,:));
    end
end
[psi,dpsidr] = rbf_func(r,s);

drdx = D./repmat(r,1,m);
dpsi = repmat(dpsidr,1,m).*drdx;

w = model.w;
y = w.'*psi;
y = y.';

dy = w.'*dpsi;
% if any(isnan(dy))
%     keyboard
% end

update = rbf_update(x, X, surrogate.update_threshold(1));

end
