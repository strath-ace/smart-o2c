function [toNodeAttributes, veinlength] = MyCostFunctionMainBelt(Inputs, fromNode, toNodeAttributes)
% This function calculates the cost of a certain connection.
% It can be altered such that it is applicable to the problem at hand.
%
% Inputs:
% * fromNode  : The node from which the cost is calculated (structure)
% * toNode    : The node to which the cost is calculated (structure)
%
% Outputs: 
% * newNode   : The node currently being evaluated (structure)
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%attributenames = fieldnames(Attributes);

%Calculate cost to change orbital elements with orbit characterstics such
%as ToF set. Currently simple formula to test functionality.
%toNode.length = ((Attributes.t_arr - Attributes.tof) - fromNode.attributes.t_arr)^5;

curr_orbit = fromNode.attributes.kep_trans;
mu = AstroConstants.Sun_Planetary_Const;

[departure_r, departure_v] = StardustTool.CartesianElementsAt(curr_orbit,toNodeAttributes.t_dep);

%Save the departure coordinates
toNodeAttributes.r_dep = departure_r;
toNodeAttributes.v_dep = departure_v;

%Retrieve the arrival coordinates
arrival_r = toNodeAttributes.r_arr;
ToF = toNodeAttributes.tof;

if Inputs.AdditionalInputs{2} == 0
[~,~,~,err,vel_initial,vel_final,~,~] = lambertMR(departure_r,      ... % Initial Vector Position
                                                  arrival_r,        ... % Final position vector
                                                  ToF*86400, ... % Time of flight [seconds]
                                                  mu,               ... % Planetary constant of the planet (mu = mass * G) [L^3/T^2]
                                                  0,                ... % Logical variable defining whether transfer is
                                                                                ... %   0: direct transfer from R1 to R2 (counterclockwise)
                                                                                ... %   1: retrograde transfer from R1 to R2 (clockwise)
                                                  0,                ... % Number of revolutions (Nrev):
                                                                                ... %   if Nrev = 0 ZERO-REVOLUTION transfer is calculated
                                                                                ... %   if Nrev > 0 two transfers are possible. Ncase should be
                                                                                ... %   defined to select one of the two.
                                                  0,                ... % Logical variable defining the small-a or large-a option in
                                                                                ... % case of Nrev>0:
                                                                                ... %   0: small-a option
                                                                                ... %   1: large-a option
                                                  0);                   % LambertMR options:
                                                                                    %   optionsLMR(1) = display options:
                                                                                    %     - 0: no display
                                                                                    %     - 1: warnings are displayed only when the algorithm does not converge
                                                                                    %     - 2: full warnings displayed
if (err~=0)
    veinlength = Inf;
    return
end
%Save initial & final lambert velocity
toNodeAttributes.lambertV_ini = vel_initial;
toNodeAttributes.lambertV_final = vel_final;
                                                                                    
% Compute the Total DeltaV - ignore arrival dV due to flyby
dv1(1) = vel_initial(1) - departure_v(1);
dv1(2) = vel_initial(2) - departure_v(2);
dv1(3) = vel_initial(3) - departure_v(3);

% dv2(1) = arrival_v(1) - vel_final(1);
% dv2(2) = arrival_v(2) - vel_final(2);
% dv2(3) = arrival_v(3) - vel_final(3);

deltaV_Departure =  abs(norm(dv1));
%deltaV_Arrival   =  abs(norm(dv2));
%deltaV_Total     = (deltaV_Departure+ deltaV_Arrival);
deltaV_Total     = deltaV_Departure;

%Save found dV in node attributes
toNodeAttributes.dV_dep = dv1;
toNodeAttributes.dV_sum = deltaV_Total;
if isempty(fromNode.attributes.dV_tot)
    toNodeAttributes.dV_tot = deltaV_Total;
else
    toNodeAttributes.dV_tot = fromNode.attributes.dV_tot + deltaV_Total; 
end
            
% ==========================================================================================                                                                        
% MARILENA
% Compute the Keplerian elements of the transfer orbit at the departure position (a in km and angles in rad)               
kep_transfer_orbit = cart2kep([departure_r, vel_initial],mu);

% Eccentricity of transfer orbit
ecc = kep_transfer_orbit(2);
 
if ecc >= 1
    veinlength = Inf;
    return
end

% True anomaly at departure position over transfer orbit [rad]
theta = kep_transfer_orbit(6);

% Compute the eccentric anomaly at the departure position [rad]
cos_E = ( cos(theta) + ecc ) / (1 + ecc * cos(theta) );
sin_E = ( sin(theta) * sqrt(1 - ecc^2) ) / (1 + ecc * cos(theta) );

E = atan2(sin_E, cos_E);
E = mod(E, 2*pi);

% Compute the mean anomaly at the departure position [rad]
M = E - ecc*sin(E); 

% Mean anomaly at the departure position [deg]
M = M * 180/pi;

% Creat object to use as 1st input in SearchAndCreateNewArc
% (angles in degrees and distance in AU!)
au2km = AstroConstants.Astronomical_Unit;
curr_departure_orbit = CelestialBody('transfer_orbit',             ... % Name of the CelestialBody (This case current trajectory
                                     kep_transfer_orbit(1)/au2km,  ... % Semimajor axis [AU]
                                     kep_transfer_orbit(2),        ... % Eccentricity
                                     kep_transfer_orbit(3)*180/pi, ... % Inclination [deg]
                                     kep_transfer_orbit(4)*180/pi, ... % Asc. Node/Raan [deg]
                                     kep_transfer_orbit(5)*180/pi, ... % Arg. Perigee [deg]
                                     M,                            ... % Mean anomoly, M at time given t0 [deg]
                                     toNodeAttributes.t_dep); % Time at which Mo is given [MJD2000]  

toNodeAttributes.kep_trans = curr_departure_orbit;

veinlength = deltaV_Total;
else
    error('The low-thrust functionality is currently not yet available.')
    DepBody = StardustTool.KeplerianElementsAt2(fromNode.attributes.kep_trans, toNodeAttributes.t_dep);
    ArrBody = StardustTool.KeplerianElementsAt2(toNodeAttributes.kep_body, toNodeAttributes.t_arr);
    TOF = toNodeAttributes.tof;
    t_dep = toNodeAttributes.t_dep;
    global mu_S
    mu_S = mu;
    
    [thetaVec,R,Phi,u_mat,DV,isLoopConverged,cont,TPrime1,Ucart,Vcart] = SphericalShapingWithFSolve(DepBody,ArrBody,t_dep,TOF,1,1e-3,2*pi/200,0);
end


