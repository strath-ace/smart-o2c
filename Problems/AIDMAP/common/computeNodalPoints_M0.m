function [Mnode1, Mnode2, error_status, theta1, E1, theta2, E2 ] = computeNodalPoints_M0( keporb1, keporb2, mu)
%% computeNodalPoints_M0: computes the Mean Anomaly of the mutal nodal points between two orbits
% 
%% Inputs:
%  keporb          : Keplerian Elements of the orbits. The keporb has
%                    to be a KeplerianElements Object
%  mu              : Standard gravitational parameter of the main body
% 
%% Outputs:
%  Mnode           : Mean anomaly of the nodal points
%  error_status    : Error indication  
%  theta           : True anomaly of the nodal points
%  E               : Eccentric anomaly of the nodal points
% 
%% Author(s):
%    Juan Manuel Romero Martin (2014), Marilena Di Carlo (2014)
% 
% REFERENCE
% 
% 
% MODIFICATIONS:
% 
% --------------------------------------------------------------------------

% Set error flag to true
error_status = 1;

% Sanity Check: First Obirt is given in Keplerian Object
if ( ( isa(keporb1, 'double') == 0 ) && ( isa(keporb1, 'CelestialBody') == 0 ) )     
    error('The first orbit must be given in interger (Planet ID or KeplerianElements Object (the rest of celestial body).')
end

% Sanity Check: Second Obirt is given in Keplerian Object
if ( ( isa(keporb2, 'double') == 0 ) && ( isa(keporb2, 'CelestialBody') == 0 ) )    
    error('The second orbit must be given in interger (Planet ID or KeplerianElements Object (the rest of celestial body).')
end   


% Get the Keplerian Orbital Elements [a e i Om om theta] [km, rad]
if isa(keporb1, 'CelestialBody')
    keporb1_arr = keporb1.getKeplerianElements().getKepArray_KM_RAD();
elseif isa(keporb1, 'double')
    keporb1_arr  = EphSS_kep(keporb1, 0); 
end

% Convert from keplearian to Cartesian 
out = kep2cart(keporb1_arr, mu);
r1 = out(1:3); 
v1 = out(4:6);    

% Get the Keplerian Orbital Elements [a e i Om om theta] [rad]
if isa(keporb2, 'CelestialBody')
    keporb2_arr = keporb2.getKeplerianElements().getKepArray_KM_RAD();
elseif isa(keporb2, 'double')
    keporb2_arr  = EphSS_kep(keporb2, 0); 
end

% Convert from keplearian to Cartesian 
out = kep2cart(keporb2_arr, mu);
r2 = out(1:3); 
v2 = out(4:6); 

% Direction of the Angular momentum
tmp_vec = cross(r1, v1);
angular_moment_orb1 = tmp_vec/norm3(tmp_vec);

% Direction of the Angular momentum
tmp_vec = cross(r2, v2);
angular_moment_orb2 = tmp_vec/norm3(tmp_vec);

% Compute the Mutual Inclination
imutual = acos(dotprod(angular_moment_orb1, angular_moment_orb2'));

% Sanity check
if imutual == 0
    error('The Mutual inclination is Zero, so the mutual nodal points are not defined.')    
end

% Compute the unit vector pointing to the Ascending Mutual Node
tmp_vec = cross(angular_moment_orb1, angular_moment_orb2);
asc_vec = tmp_vec / norm3(tmp_vec); 

% [a e i Om om theta] [km, rad]
ecc = keporb2_arr(2);
inc = keporb2_arr(3);
OM  = keporb2_arr(4);
om  = keporb2_arr(5);

% Compute the unit position vector of the second orbit corresponding to
% eccentricity anomaly equal to ZERO. (the direction of the periapsis) 
x0 = [cos(OM)*cos(om) - sin(OM)*sin(om)*cos(inc), ... 
      sin(OM)*cos(om) + cos(OM)*sin(om)*cos(inc), ... 
      sin(om)*sin(inc) ];
 
% Compute the Periapsis Argument of the Second Orbit in the Mutual
% reference frame
cos_om_node = dotprod(asc_vec, x0');
sin_om_node = norm3(cross(asc_vec, x0'))/(norm3(asc_vec)*norm3(x0));

temp = atan2(sin_om_node, cos_om_node) * 180/pi;
theta = temp;

% Sanity Check: Quadrant Check
if x0(3) > 0
   theta  = 360 - theta;
end  

% Compute the First Node --------------------------------------------------

% Eccentric Anomaly has to be computed using theta and not sin_om_node and
% cos_om_node (as was previously done)
cos_E = (cosd(theta) + ecc ) / (1 + ecc*cosd(theta));

beta  =  sqrt(1 - ecc*ecc);
sin_E =  beta*sind(theta) / (1 + ecc*cosd(theta)); 


% Compute the Eccentricity Anomaly Value
E0 = atan2(sin_E, cos_E);
E0 = mod(E0, 2*pi);

% DEBUGGING JMRM
E1 = E0 * 180/pi;
theta1 = theta;

% Compute the Mean Anomaly
M0 =  E0 - ecc * sin_E;

% Store the Second Nodal 
Mnode1 = M0 * 180/pi;
Mnode1 = mod(Mnode1, 360);

% Compute the Second Node -------------------------------------------------

% The Mean Anomaly of the 2nd node is not Mnode1 + 180

% True anomaly second node
theta2 = mod(theta1 + 180, 360);

% Eccentric anomaly second node                                           
cos_E2 = (cosd(theta2) + ecc ) / (1 + ecc*cosd(theta2));

beta  =  sqrt(1 - ecc*ecc);
sin_E2 =  beta*sind(theta2) / (1 + ecc*cosd(theta2)); 

% Compute the Eccentricity Anomaly Value
E2 = atan2(sin_E2, cos_E2);

E2 = mod(E2, 2*pi);


% Mean anomaly second node
Mnode2 = E2 - ecc * sin(E2);
Mnode2 = Mnode2 * 180/pi;

E2 = E2 * 180/pi;

% Set error flag to false
error_status = 0;

end

