[ephemerides, names] = xlsread('DiameterGreater10kmM90.xlsx');
%Note: all asteroids names must start with a letter
Asteroids.initialize = [];
                     
for k = 1 : length(names)
    
    name = names(k);
    
    temporary = CelestialBody(names{k}, ...
                              ephemerides(k,1), ...
                              ephemerides(k,2), ...
                              ephemerides(k,3), ...
                              ephemerides(k,4), ...
                              ephemerides(k,5), ...
                              ephemerides(k,6), ...
                              ephemerides(k,7));

    assignin('base', names{k}, temporary);
    
    Asteroids.(char(name)) = temporary;
    
%     Asteroids = mergetest(Asteroids,name,temporary);
    
end

Asteroids = rmfield(Asteroids,'initialize');
clearvars -except Asteroids
save('AsteroidsMainBelt90.mat')


