% This script calculates the passing epochs of all the asteroids through
% their nodal points, given a time-frame
%
% Inputs:
% * epoch_start : The epoch from which the nodal points are to be
%                 calculated
% * epoch_end   : The epoch up to which the nodal points should be
%                 calculated
% * Asteroids   : File containing the characteristics of the asteroids
%
% Outputs: 
% * epochsnode  : Epochs at which the asteroid passes its nodal points
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Initialize
close all; clear all; clc
addpath(genpath('astro_tool'));

SaveDir = 'IO_Dir/';
addpath(SaveDir);

%Input start and end epochs
epoch_start = 7304.5;
epoch_end = 10957.5;

%Load the file containing the asteroids
Asteroids = Asteroids();

%Retrieve the Sun's standard gravitational parameter
mu = AstroConstants.Sun_Planetary_Const;

%Obtain the asteroid names
asteroidnames = fieldnames(Asteroids);

%Generate a virtual orbit that the asteroid's orbit will be intersected with
virtualorbit = struct('a',1,...
                  'e',0.035999947927132,...
                  'i',0,...
                  'OM',-11.26064,...
                  'W',102.94719,...
                  'M',0,... 
                  't',0);

%Save the virtual orbit as a CelestialBody object         
virtorbit = CelestialBody('virtualorbit',...
        virtualorbit.a,...
        virtualorbit.e,...
        virtualorbit.i,...
        virtualorbit.OM,...
        virtualorbit.W,...
        virtualorbit.M,...
        virtualorbit.t);

%Loop over all the asteroids
for i = 1:length(asteroidnames)
    
%Convert the name of the asteorid currently being evaluated to characters
asteroid = char(asteroidnames(i));


%Retrieve the Keplerian elements
asteroidorbit = Asteroids.(asteroid).getKeplerianElements();

%Compute the location of Nodal points (the intersection between the two
%orbits)
[Mnode1, Mnode2, error_status, theta1, E1, theta2, E2 ] = computeNodalPoints_M0(virtorbit,Asteroids.(asteroid),mu);

%Find the passing epochs
epochsnode1 = computeNodalPassesEpochs(asteroidorbit,Mnode1,epoch_start, epoch_end, mu);
epochsnode2 = computeNodalPassesEpochs(asteroidorbit,Mnode2,epoch_start, epoch_end, mu);

%Save the passing epochs into a single array
epochsnode{i,1} = sort([epochsnode1 epochsnode2]);

end

%Output the passing epoch into a file
save(strcat(SaveDir,'epochsnode'),'epochsnode')


