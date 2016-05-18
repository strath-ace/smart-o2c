%Find epochs of all asteroids
close all; clear all; clc
addpath('astro_tool')

AsteroidsMainBelt();

epoch_start = '2020/01/01';
epoch_end = '2030/01/01';
decimalplaces = 2;


epoch_start = date2mjd2000([datevec(epoch_start,'yyyy/mm/dd')]);
epoch_end = date2mjd2000([datevec(epoch_end,'yyyy/mm/dd')]);

asteroidnames = fieldnames(Asteroids);
for i = 2:length(asteroidnames)
asteroid = char(asteroidnames(i));

mu = AstroConstants.Sun_Planetary_Const;


satorbit = struct('a',1,...
                  'e',0.0167,...
                  'i',0,...
                  'OM',-11.26064,...
                  'W',102.94719,...
                  'M',0,... %Current mean anomaly
                  't',0); %Current  time
    
              
satellite = CelestialBody('satellite',...
        satorbit.a,...
        satorbit.e,...
        satorbit.i,...
        satorbit.OM,...
        satorbit.W,...
        satorbit.M,...
        satorbit.t);


asteroidorbit = Asteroids.(asteroid).getKeplerianElements();
satelliteorbit = satellite.getKeplerianElements();

%Compute Nodal points
[Mnode1, Mnode2, error_status, theta1, E1, theta2, E2 ] = computeNodalPoints_M0(satellite,Asteroids.(asteroid),mu);

epochsnode1 = computeNodalPassesEpochs(asteroidorbit,Mnode1,epoch_start, epoch_end, mu);
idnode1 = ones(1,length(epochsnode1));
epochsnode2 = computeNodalPassesEpochs(asteroidorbit,Mnode2,epoch_start, epoch_end, mu);
idnode2 = 2*ones(1,length(epochsnode2));

epochsnode{i,1} = [epochsnode1 epochsnode2];
identifiers{i,1} = [idnode1 idnode2];
for j = 1:length(epochsnode{i,1})
    cartnode{i,1}(j,:) = StardustTool.CartesianElementsAt(Asteroids.(asteroid),epochsnode{i,1}(j));
    xnode{i,1}(j) = cartnode{i,1}(j,1);
    ynode{i,1}(j) = cartnode{i,1}(j,2);
    znode{i,1}(j) = cartnode{i,1}(j,3);
    
end
orbitchars(i,:) = [asteroidorbit.a asteroidorbit.e asteroidorbit.i asteroidorbit.OM asteroidorbit.W asteroidorbit.M0 asteroidorbit.t0];

end

for i = 1:length(epochsnode)
    for j = 1:length(epochsnode{i})
        epochsnode{i,:} = round(epochsnode{i,:},decimalplaces);
    end
end

save('epochsnode','epochsnode')
save('asteroidorbitchars','orbitchars')


