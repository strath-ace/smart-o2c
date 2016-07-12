function [] = SaveTrajectorySolution(Sequence,ListNodes,filename)
%This function can be used to save the solution to a .mat file
%
% Inputs:
% * Sequence     : List of unique IDs to be saved
% * ListNodes    : Structure containing (at least) the nodes to be saved
% * filename     : Name of the .mat file
%
% Outputs: 
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Initialize the writing to the file
fid = fopen(filename,'w');
fprintf(fid,'-----------------------------------------------------------------------------------------------\n');
fprintf(fid,'\tAsteroid\tWaiting Time\tLambert T_dep\tLambert ToF\t\tLambert T_arr\t\tdV\n');
fprintf(fid,'-----------------------------------------------------------------------------------------------\n');
formatspec = '%12s\t%12.4f\t%12s\t%12.4f\t%12s\t%8.2f\n';

%Retrieve the sequence & initialize the total dV so far
Targets = Sequence';
dV_tot = 0;

%Loop over all the targets
for i = 2:length(Targets)
    
    %Find the target names
    targetsplit = strsplit(Targets{i},'___');
    targetnamesplit = strsplit(char(targetsplit(2)),'__');
    targetname = targetnamesplit{1};
    TargetsNames{i-1}= char(targetname);

    %Obtain the first character of each target
    firstchar(i-1) = targetname(1);
end

%If all the first characters are "a", the AIDMAP algorith added "a" in
%front of each target, to ensure the field names of the structure are
%valid. As such, this first character is removed in this caes
if sum(firstchar=='a')==length(firstchar)
    for i = 1:length(TargetsNames)
        TargetsNames{i} = TargetsNames{i}(2:end);
    end
end

%Loop over all targets
for i = 2:length(Targets)
    
    %Retrieve the waiting time
    targetparent = ListNodes.(Sequence{i}).parent;
    waitingtime = ListNodes.(Sequence{i}).attributes.t_dep - ListNodes.(targetparent).attributes.t_arr;

    %Find the departure date of the Lambert arc
    lambertdepdate2000 = ListNodes.(Sequence{i}).attributes.t_dep;
    lambertdepdate = datestr(mjd20002date(lambertdepdate2000),'yyyy/mm/dd');

    %Obtain the time of flight of the Lambert arc
    lamberttof = ListNodes.(Sequence{i}).attributes.tof;

    %Retrieve the arrival date at target
    lambertarrdate2000 = ListNodes.(Sequence{i}).attributes.t_arr;
    lambertarrdate = datestr(mjd20002date(lambertarrdate2000),'yyyy/mm/dd');

    %Obtain the dV for this arc and the total dV so far
    dV = ListNodes.(Sequence{i}).attributes.dV_sum;
    dV_tot = dV_tot+dV;

    %Print to file
    fprintf(fid,formatspec,TargetsNames{i-1},waitingtime,lambertdepdate,lamberttof,lambertarrdate,dV);
end

%Add summary to the bottom of the file
fprintf(fid,'-----------------------------------------------------------------------------------------------\n');
fprintf(fid,'%12s\t%12.4f\t%12s\t%12.4f\t%12s\t\t\t\t\t\t\t%8.2f\n','Total:',[],[],[],[],dV_tot);

%Close the file
fid = fclose('all');
end

