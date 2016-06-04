function [X_ge] = KeplElem2rvMatrix(Omega, i, w, X_pf)

% KeplElem2rvMatrix: Function that calculates the position vector in a cartesian, inertial
% reference frame centered on the attraction body given the position vector
% in the perifocal plane (plane of the orbit) and the angles that define
% an orbit in 3D space
%
% INPUT
% Kepler parameters Omega,w,i and X_pf (matrix containing the row
% vectors [X;Y;Z]

% OUTPUT
% X_ge: matrix containign the vectors in the cartesian inertial
% frame centered on the attraction body

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation

% Compute the rotation matrices
R1p = [ cos(Omega)  -sin(Omega)  0;            
         sin(Omega)  cos(Omega)  0;             
             0        0     1];        


R2p = [1       0          0;
       0   cos(i)  -sin(i);
        0  sin(i)  cos(i)];


R3p = [ cos(w)  -sin(w)    0 ;
         sin(w)  cos(w)  0;
           0       0     1];
       
% Rotation matrix
T_pf2ge = R1p*R2p*R3p;    

% Obtain the first vector
x_ge = T_pf2ge*X_pf(:,1);
X_ge = [x_ge];

% If there are more than 2 elements, then continue and transform all the
% vectors contained in the matrix and put them inside a new matrix
if min(size(X_pf)) > 2
    n = max(size(X_pf));
    for i=2:n      
        x_ge = T_pf2ge*X_pf(:,i);
        X_ge = [X_ge x_ge];     
    end
end

end

