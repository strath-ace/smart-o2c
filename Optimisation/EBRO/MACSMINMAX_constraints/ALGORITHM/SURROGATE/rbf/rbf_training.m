function model = rbf_training(x, y, surrogate)

basis = surrogate.corrfun;
s = surrogate.s;

N = size(x,1);
R = zeros(N);

rbf_func = strcat('rbf_',lower(basis));
rbf_func = str2func(rbf_func);

for i = 1:N
    for j = 1:N
        R(i,j) = norm(x(i,:) - x(j,:));
    end
end
PSI = rbf_func(R,s);

w = PSI\y;

model.X = x;
model.Y = y;
model.PSI = PSI;
model.w = w;
model.basis = basis;
model.s = s;

end
