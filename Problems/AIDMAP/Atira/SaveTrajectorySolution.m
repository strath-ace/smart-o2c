function [] = SaveTrajectorySolution(Sequence, ListNodes, filename)
%% SaveTrajectorySolution: This function can be used to save a trajectory to a text file
%
%% Inputs:
% * Sequence    : The string array with the node IDs of the nodes to be
%                 saved
% * ListNodes   : Structure containing (at least) all the information on
%                 the nodes in the solution to be saved
% * filename    : The name of the file to be saved [string]
%
%% Outputs: 
% * 
%
%% Author: Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

%Open the file
fid = fopen(filename, 'w');

%Create the header
fprintf(fid, '-----------------------------------------------------------------------------------------------\n');
fprintf(fid, '\tAsteroid\tWaiting Time\tLambert T_dep\tLambert ToF\t\tLambert T_arr\t\tdV\n');
fprintf(fid, '-----------------------------------------------------------------------------------------------\n');
formatspec = '%12s\t%12.4f\t%12s\t%12.4f\t%12s\t%8.2f\n';

%Load the sequence
Sequence = Sequence';

%Initialise the total dV so far
dV_tot = 0;

%Loop over all the planets/asteroids in the sequence
for i = 2:length(Sequence)
    
%Find names of the planets/asteroids
bodysplit = strsplit(Sequence{i}, '___');
bodynamesplit = strsplit(char(bodysplit(2)), '__');
bodyname = bodynamesplit{1};
BodyNames{i-1}= char(bodyname);

firstchar(i-1) = bodyname(1);
end

%If all the first characters are "a", AIDMAP has added this first character
%to prevent structure names from starting with a number. Thus, this initial
%character should be removed to obtain the original asteroid/planet name
if sum(firstchar=='a')==length(firstchar)
    for i = 1:length(BodyNames)
        BodyNames{i} = BodyNames{i}(2:end);
    end
end

%Loop over all the planets/asteroids in the sequence
for i = 2:length(Sequence)
    
%Find the waiting time
asteroidparent = ListNodes.(Sequence{i}).parent;
waitingtime = ListNodes.(Sequence{i}).attributes.t_dep - ListNodes.(asteroidparent).attributes.t_arr;

%Obtain the departure date of the Lambert arc
lambertdepdate2000 = ListNodes.(Sequence{i}).attributes.t_dep;
lambertdepdate = datestr(mjd20002date(lambertdepdate2000), 'yyyy/mm/dd');

%Get the time of flight of the Lambert arc [Days]
lamberttof = ListNodes.(Sequence{i}).attributes.tof;

%Find the arrival date of the Lambert arc
lambertarrdate2000 = ListNodes.(Sequence{i}).attributes.t_arr;
lambertarrdate = datestr(mjd20002date(lambertarrdate2000), 'yyyy/mm/dd');

%Get the dV [km/s]
dV = ListNodes.(Sequence{i}).attributes.dV_sum;
dV_tot = dV_tot+dV;

%Add line to the file
fprintf(fid, formatspec, BodyNames{i-1}, waitingtime, lambertdepdate, lamberttof, lambertarrdate, dV);
end

%Add summary at the bottom of the file
fprintf(fid, '-----------------------------------------------------------------------------------------------\n');
fprintf(fid, '%12s\t%12.4f\t%12s\t%12.4f\t%12s\t\t\t\t\t\t\t%8.2f\n', 'Total:', [], [], [], [], dV_tot);

%Close the file
fid = fclose('all');
end

