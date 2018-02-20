% Downloaded from:
% http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC%202011%2d%20RWP&FolderCTID=&View=%7bDAF31868%2d97D8%2d4779%2dAE49%2d9CEC4DC3F310%7d

% clear all
% mex cec14_func.cpp -DWINDOWS
func_num=1;
D=10;
Xmin=-100;
Xmax=100;
pop_size=100;
iter_max=5000;
runs=1;
fhd=str2func('cec14_func');
for i=24:24
    func_num=i;
    for j=1:runs
        i,j,
        [gbest,gbestval,FES]= PSO_func(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
        xbest(j,:)=gbest;
        fbest(i,j)=gbestval;
        fbest(i,j)
    end
    f_mean(i)=mean(fbest(i,:));
end



% for i=1:30
% eval(['load input_data/shift_data_' num2str(i) '.txt']);
% eval(['O=shift_data_' num2str(i) '(1:10);']);
% f(i)=cec14_func(O',i);i,f(i)
% end
