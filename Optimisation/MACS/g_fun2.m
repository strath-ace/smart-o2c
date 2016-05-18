function g=g_fun2(f,lambda,z,maxs)

delta = maxs-z;
delta = delta.*(delta~=zeros(size(delta)))+ones(size(delta)).*(delta==zeros(size(delta))); % if there's only one value, max and min are the same, so this avoids a useless division by zero

f = f-repmat(z,size(f,1),1);
f = f./repmat(delta,size(f,1),1);

g=max(lambda.*abs(f-z));

return