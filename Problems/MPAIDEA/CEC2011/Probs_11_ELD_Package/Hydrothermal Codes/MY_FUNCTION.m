% Downloaded from:
% http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC%202011%2d%20RWP&FolderCTID=&View=%7bDAF31868%2d97D8%2d4779%2dAE49%2d9CEC4DC3F310%7d

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HYDROTHERMAL SCHEDULING (Fitness Function)
%% Guided by : Dr. P.K. Rout, SOA University
%% Coded by  : Krishnanand K.R., Santanu Kumar Nayak
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fitness Function to be called from the optimization algorithm
%% Evaluates a population of row vectors
function [y Count] = MY_FUNCTION(input_array)
siz = size(input_array,1);
Count = siz ;
y = zeros(siz,1);
for i =1:siz
    y(i,1) =  Fn_Eval(input_array(i,:));
end
end
%% Evaluates a single row vector
function y = Fn_Eval(x)
x=round(x*10000)/10000; %% For fixing the 4 digit precision
y = fn_HT_ELD_Case_1(x); %% Change Here to Call a Different System
end
%%