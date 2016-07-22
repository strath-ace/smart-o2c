function [Asteroids] = EvaluateAsteroidsXLS(AsteroidsFileName,MatFileName,NameFile)
%% MyEndConditions: This function evaluates an XLS file containing the orbital elements of asteroids, and outputs a structure with the CelestialBody objects
%
%% Inputs:
% * AsteroidsFileName    : The filename of the XLS file containing the orbital elements[string]
% * MatFileName          : The name of the mat file as which the output
%                              structure should be saved [string]
% * NameFile             : The name of the file that contains all the
%                              asteroid names (needed for the options.Cities input of AIDMAP)
%
%% Outputs: 
% * Asteroids  : The structure containing the CelestialBody objects of the
%                   asteroids
%
%% Author: Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

%Read the XLS file
[ephemerides, txt] = xlsread(AsteroidsFileName);

%Initialise the structure
Asteroids.initialise = [];

%Convert MJD to MJD2000
for i = 1:length(ephemerides)
    ephemerides(i,7) = mjd2mjd2000(ephemerides(i,7));
end

%Initialise the cell array that will contain the correct names
correctnames(1) = txt(1,end);

%Initialise a flag to see if one of the names start with a number (AIDMAP
%will otherwise not work)
containsnum = 0;

%Check if any asteroid names start with a number
for i = 1:length(txt)
    firstchar = txt{i}(1);
    firstcharnum = str2num(firstchar);
    if (~isempty(firstcharnum))
        containsnum = 1;
    end
end

%Make sure asteroid names have the correct format (letter first, only
%numbers and letters in the name)
for i = 2:length(txt)
    currentname = char(txt(i,end));
    if containsnum == 1
        correctnames{i} = strcat('a',strrep(strrep(currentname,' ',''),'-',''));
    else
        correctnames{i} = strrep(strrep(currentname(1:end-1),' ',''),'-','');
    end
end

%Save the names to a secondary variable
names = correctnames';
    
%Loop over all the asteroids    
for k = 1 : length(names)
    
    %For ease of reading, define the name currently evaluated
    name = names{k};
    
    %Create a  temporary CelestialBody object
    temporary = CelestialBody(names{k}, ...
                              ephemerides(k,1), ...
                              ephemerides(k,2), ...
                              ephemerides(k,3), ...
                              ephemerides(k,4), ...
                              ephemerides(k,5), ...
                              ephemerides(k,6), ...
                              ephemerides(k,7));

    %Add this object to the Asteroids structure
    Asteroids.(char(name)) = temporary;
    
    
end

%Remove the field that was needed for the initialisation
Asteroids = rmfield(Asteroids,'initialise');

%Save the structure
save(MatFileName,'Asteroids');

%Save the names of all the asteroids
fileID = fopen(NameFile,'w');
formatspec = '%12s\n';
for i = 2:length(names)
    fprintf(fileID,formatspec,char(names{i}));
end

end

