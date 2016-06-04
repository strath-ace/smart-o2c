function [Asteroids] = EvaluateAsteroidsXLS(AsteroidsFileName,MatFileName,NameFile)

[ephemerides, txt] = xlsread(AsteroidsFileName);
%Note: all asteroids names must start with a letter
Asteroids.initialize = [];

%Convert MJD to MJD2000
for i = 1:length(ephemerides)
    ephemerides(i,7) = mjd2mjd2000(ephemerides(i,7));
end

correctnames(1) = txt(1,end);

containsnum = 0;
%Check if any target names start with a number
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
    if isempty(strfind(currentname,'('))
        splittedname = strsplit(currentname,' ');
        if containsnum == 1
            correctnames{i} = strcat('a',splittedname(end));
        else
            correctnames{i} = splittedname(end);
        end
    else
        splittedname = strsplit(currentname,'(');
        splittedname2 = char(splittedname(end));
        if containsnum == 1
            correctnames{i} = strcat('a',strrep(strrep(splittedname2(1:end-1),' ',''),'-',''));
        else
            correctnames{i} = strrep(strrep(splittedname2(1:end-1),' ',''),'-','');
        end
    end
end

names = correctnames';
    
    


for k = 1 : length(names)
    
    name = names{k};
    
    temporary = CelestialBody(names{k}, ...
                              ephemerides(k,1), ...
                              ephemerides(k,2), ...
                              ephemerides(k,3), ...
                              ephemerides(k,4), ...
                              ephemerides(k,5), ...
                              ephemerides(k,6), ...
                              ephemerides(k,7));

    %assignin('base', names{k}, temporary);
    
    Asteroids.(char(name)) = temporary;
    
%     Asteroids = mergetest(Asteroids,name,temporary);
    
end

Asteroids = rmfield(Asteroids,'initialize');
%clearvars -except Asteroids AsteroidsFileName MatFileName

save(MatFileName,'Asteroids');

fileID = fopen(NameFile,'w');
formatspec = '%12s\n';
for i = 2:length(names)
    fprintf(fileID,formatspec,char(names{i}));
end
end

