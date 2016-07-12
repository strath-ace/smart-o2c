function [] = SaveTrajectorySolution(Sequence,ListNodes,filename)
% This function handles the ramification to new nodes. It does so by
% generating a preset number of random nodes and making a probabilistic
% selection based on the cost function.
%
% Inputs:
% * Sequence    : Structure containing the graph
% * Inputs      : Structure containing the PhysarumSolver inputs
% * Agents      : The structure containing the agents
% * agent       : A string containing the name of the current agent
%
% Outputs: 
% * ListNodes       : ListNodes structure where the radii have been updated with
%                 the dilation and evaporation
% * Agents      : The updated structure containing the agents
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

fid = fopen(filename,'w');
fprintf(fid,'-----------------------------------------------------------------------------------------------\n');
fprintf(fid,'\tAsteroid\tWaiting Time\tLambert T_dep\tLambert ToF\t\tLambert T_arr\t\tdV\n');
fprintf(fid,'-----------------------------------------------------------------------------------------------\n');
formatspec = '%12s\t%12.4f\t%12s\t%12.4f\t%12s\t%8.2f\n';
%formatspec = '%s';
Asteroids = Sequence';
dV_tot = 0;

for i = 2:length(Asteroids)
%Find asteroid names
asteroidsplit = strsplit(Asteroids{i},'___');
asteroidnamesplit = strsplit(char(asteroidsplit(2)),'__');
asteoroidname = asteroidnamesplit{1};
AsteroidsNames{i-1}= char(asteoroidname);
firstchar(i-1) = asteoroidname(1);
end

if sum(firstchar=='a')==length(firstchar)
    for i = 1:length(AsteroidsNames)
        AsteroidsNames{i} = AsteroidsNames{i}(2:end);
    end
end
    
for i = 2:length(Asteroids)
%Waiting time
asteroidparent = ListNodes.(Sequence{i}).parent;
waitingtime = ListNodes.(Sequence{i}).attributes.t_dep - ListNodes.(asteroidparent).attributes.t_arr;

%Departure Date Lambert Arc
lambertdepdate2000 = ListNodes.(Sequence{i}).attributes.t_dep;
lambertdepdate = datestr(mjd20002date(lambertdepdate2000),'yyyy/mm/dd');

%ToF Lambert Arc [Days]
lamberttof = ListNodes.(Sequence{i}).attributes.tof;

%Arrival Date at Asteroid Node
lambertarrdate2000 = ListNodes.(Sequence{i}).attributes.t_arr;
lambertarrdate = datestr(mjd20002date(lambertarrdate2000),'yyyy/mm/dd');

%dV [km/s]
dV = ListNodes.(Sequence{i}).attributes.dV_sum;
dV_tot = dV_tot+dV;
fprintf(fid,formatspec,AsteroidsNames{i-1},waitingtime,lambertdepdate,lamberttof,lambertarrdate,dV);
end

fprintf(fid,'-----------------------------------------------------------------------------------------------\n');
fprintf(fid,'%12s\t%12.4f\t%12s\t%12.4f\t%12s\t\t\t\t\t\t\t%8.2f\n','Total:',[],[],[],[],dV_tot);
% fprintf(fid,formatspec,char(AsteroidsNames),lambertdepdate,lamberttof,lamberttof,lambertarrdate,dV)
fid = fclose('all');
end

