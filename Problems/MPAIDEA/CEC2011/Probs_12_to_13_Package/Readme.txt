% Downloaded from:
% http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC%202011%2d%20RWP&FolderCTID=&View=%7bDAF31868%2d97D8%2d4779%2dAE49%2d9CEC4DC3F310%7d



While Using the problem formulation of Trajectory Problems i.e. Cassini2 and messengerfull.mat


1)First of all the messengerfull.mat or Cassini2.mat file should be loaded.

2)Boundary values for each dimension of messenger problem is given in the getlimit_messenger file.
Here in the .rar folder getlimit_cassini2.m file is prepared for Cassini2 problem that can be called while assigning the boundary values for each dimension.

3) At the time of function evaluation each population ('1'X'd')['d'- Is the dimension of the population.] should be passed as an argument along with another argument i.e. 'MGADSMproblem'.

Example- 

for i = 1 : NP %['NP'=> Population Size]
    calcutaed_fitness_value(i,1) = cassini2(x(i,:),MGADSMproblem); %['x'=> The variable denoting %  population]
end

                                                  OR

for i = 1 : NP %['NP'=> Population Size]
    calcutaed_fitness_value(i,1) = messengerfull(x(i,:),MGADSMproblem); %['x'=> The variable denoting  %  population]
end