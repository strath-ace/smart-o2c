function [kep, E, M, dt, p]=cart2kep(in, mu)

% Convertion from cartesian position and velocity to keplerian elements.
%
%   kep = cart2kep(in, mu)
%
% All the units to be consistent, angles in radians.
%
% INPUT
%   in      State in cartesian coordinates (position, velocity).
%   mu      Planetary constant.
%
% OUTPUT
%   kep     Vector of Keplerian elements:
%               kep = [a, e, i, Om, om, th]
%                   0 <= i <= pi [rad]
%                   0 <= Om <= 2*pi [rad]
%                   0 <= om <= 2*pi [rad]
%                   0 <= th <= 2*pi [rad]
%   E       Eccentric anomaly, hyperbolic anomaly or parabolic anomaly (for
%           definitions see Vallado pag. 49).
%   M       Mean anomaly.
%   dt      Time from the pericentre passage.
%   p       Parameter.
%
% FUNCTIONS CALLED
%   (none)
%
% Original function: ca2kep, (c) Massimiliano Vasile 2002
% Matteo Ceriotti, 08-02-2007
% Revised by: Camilla Colombo, Matteo Ceriotti, 17-04-2007
% Modified by: Daniel Novak, 06-11-2007: numerical roundoff error in the
%                                        calculation of theta (line 75) 
%
% ------------------------- - SpaceART Toolbox - --------------------------

r=in(1:3);
v=in(4:6);

nr=sqrt(r(1)^2+r(2)^2+r(3)^2); % Norm of r
h=[r(2)*v(3)-r(3)*v(2), r(3)*v(1)-r(1)*v(3), r(1)*v(2)-r(2)*v(1)]; % cross(r, v): Angular momentum vector
nh=sqrt(h(1)^2+h(2)^2+h(3)^2); % Norm of h
p=nh^2/mu;
e=1/mu*[v(2)*h(3)-v(3)*h(2), v(3)*h(1)-v(1)*h(3), v(1)*h(2)-v(2)*h(1)]-r/nr; % Eccentricity vector
ne=sqrt(e(1)^2+e(2)^2+e(3)^2); % Eccentricity
if ne==1 % Parabola
    a=Inf;
    warning('cart2kep:parabola', 'Parabola. Semi-major axis is Inf.');
else
    a=p/(1-ne^2); % Semi-major axis
end
i=acos(h(3)/nh); % Inclination

% keyboard

% Line of nodes vector
if i~=0 && i~=pi % Inclined orbit
    n=[-h(2), h(1), 0]/sqrt(h(1)^2+h(2)^2); % n=cross([0 0 1], h); n=n/norm(n);
    %n=n/sqrt(n(1)^2+n(2)^2+n(3)^2); % Normalisation to 1 of nv
else % Zero-inclination orbit: n is not defined
    n=[1, 0, 0]; % Arbitrary choic
end

% Argument of the ascending node
Om=acos(n(1));
if n(2)<0
    Om=2*pi-Om;
end


if ne~=0 % Non circular orbit
    % Argument of the pericentre

    om=acos((n(1)*e(1)+n(2)*e(2)+n(3)*e(3))/ne); % acos(dot(n, e)/ne)

%         if e(3)<0         % 06-11-2007 version 
    if e(3)<-10^(-9)        % 16-07-2014 version (evaluate if it is close to the value)
        om=2*pi-om;
    end
%     
    % Check for equatorial orbit (16-07-2014):
    % Ref. Vallado 3rd ed., pag. 121 (Elliptical equatorial orbits)
    if i == 0 
        om = acos(e(1)/ne);
        if e(2)<0
             om=2*pi-om;
        end
    end
%     


    % True anomaly
    th=acos(min(max((e(1)*r(1)+e(2)*r(2)+e(3)*r(3))/ne/nr, -1), 1)); % acos(dot(e, r)/ne/nr);
    if dot(r, v)<0
        th=2*pi-th;
    end
else % Circular orbit: e is not defined
    e_conv=n; % Arbitrary eccentricity vector
    om=0;
    % True anomaly
    th=acos((e_conv(1)*r(1)+e_conv(2)*r(2)+e_conv(3)*r(3))/nr); % acos(dot(e_conv, r)/nr);
    if dot(r, v)<0
        th=2*pi-th;
    end
end

kep=[a, ne, i, Om, om, th];

% The following formulas have been found in Vallado chapter 2.

if nargout>1 % E (= eccentric or hyperbolic anomaly) required as output
    if ne>0 && ne<1 % Ellipse
        sinE=sin(th)*sqrt(1-ne^2)/(1+ne*cos(th));
        cosE=(ne+cos(th))/(1+ne*cos(th));
        E=atan2(sinE, cosE);
        if E<0
            E=E+2*pi;
        end
    elseif ne==0 % Circumference
        E=th;
    elseif ne>1 % Hyperbola
        sinhH=sin(th)*sqrt(ne^2-1)/(1+ne*cos(th));
        % coshH=(e+cos(th))/(1+e*cos(th)); % Not needed
        H=asinh(sinhH);
        E=H;
    elseif ne==1 % Parabola
        B=tan(th/2); % Parabolic anomaly
        E=B;
    end
    
    if nargout>2 % M (= mean anomaly) required as output
        if ne>0 && ne<1 % Ellipse
            M=E-ne*sinE;
        elseif ne>1 % Hyperbola
            M=ne*sinhH-H;
        elseif ne==1 % Parabola
            M=B+B^3/3;
        elseif ne==0 % Circumference
            M=E;
        end
        
        if nargout>3 % dt (= time from the pericentre passage) required as output
            if ne>=0 && ne<1 % Ellipse or circumference
                n=sqrt(mu/a^3);
            elseif ne>1 % Hyperbola
                n=sqrt(-mu/a^3);
            elseif ne==1 % Parabola
                n=2*sqrt(mu/p^3);
            end
            dt=M/n;
        end
    end
end

return