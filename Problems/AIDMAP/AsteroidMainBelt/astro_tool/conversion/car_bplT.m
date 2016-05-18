function x_bpl = car_bplT(x_car,U_car,vp_car)

% Vector reference frame Transformation.
% Cartesian reference frame to b-plane reference frame.
%
%	x_bpl = car_bplT(x_car,U_car,vp_car)
%
% car reference frame: {x,y,z}
%   inertial reference frame
% The b-plane is the plane perpendicular to the incoming relative velocity
% of the small body at the planet arrival, containing the planet.
% b-plane reference frame: {xi,eta,zeta}
%   eta-axis: as the incoming relative velocity of the small body on its
%             arrival.
%   zeta-axis: in the b-plane, direction opposite to the projection of the
%              heliocentric velocity of the planet on the b-plane.
%   xi-axis: in the b-plane, completes the reference frame.
%
% INPUT:
%        x_car = vector expressed in {x,y,z}
%        U_car = velocity of the small body relative to the planet,
%                expressed in {x,y,z}
%        vp_car = orbital velocity of the planet expressed in {x,y,z}
%
% OUTPUT:
%        x_bpl = vector expressed in {xi,eta,zeta}
%
% FUNCTIONS CALLED: none
%
% - Camilla Colombo - 04/05/2007
% - Revised by Matteo Ceriotti - 21/05/2007
% - Matteo Ceriotti - 11/02/2008: Help improved.
%
% ------------------------- - SpaceART Toolbox - --------------------------

x_car = x_car(:);

nn = U_car/norm(U_car);
ee = cross(vp_car,nn)/norm(cross(vp_car,nn));
cc = cross(ee,nn);

x_bpl = [ee nn cc]'*x_car;

return