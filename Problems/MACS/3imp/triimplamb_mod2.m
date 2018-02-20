% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [JV]=triimplamb_mod2(yav)
% [JV,C]=triimplamb(yav)
if nargin<1
    yav=load('yav.dat');
    vuoto=1;
else
    vuoto=0;
end
%yav=[1.67518e-05 6.33367e-01 3.66363e-01 4.07518e-01 6.21401e-01 3.81441e-01];
% Valutazione delle funzioni di costo per un trasferimento complanare
% Circ. LEO -> high ecc. (Molnyia-like) orbit


% Parametri pianeta
mu = 3.986e5;  % = GM [km^3/s^2] parametro gravitazionale del pianeta

% Parametri orbita LEO (partenza - Orbita 1)
hLEO = 350; % [km]
rE = 6371;   % [km]
rLEO = 7000;
wLEO = sqrt(mu/rLEO^3);
TLEO = 2*pi/wLEO/3600; % [h]

% Parametri orbita H.e. (arrivo - Orbita 2)
mEd1=sqrt(mu/42164^3);
eHEO = 0.;  % [km]
THEO = 2*pi/mEd1/3600; % [h]
wHEO = 2*pi/THEO/3600;
aHEO = (mu/wHEO^2)^(1/3);   % [km]
pHEO = aHEO*(1-eHEO^2);       % [km]

the0_C = 0;  % Posizione chaser al tempo t0 su Orbita 1
the0_T = 0;  % Posizione target al tempo t0 su Orbita 2

p = [mu rLEO the0_C wLEO aHEO eHEO pHEO the0_T  THEO];
NwaitMAX = 0.9*THEO/TLEO;           % N max orbite di attesa
DincrW = 0.025*TLEO;      % Time increment for waiting time (fraction of LEO period)
%Twait = 0:DincrW:NwaitMAX*TLEO;

DincrT = 0.02*TLEO;      % Time increment for tranfer time (fraction of h.ecc. period)
%Ttran = DincrT:DincrT:0.9*THEO;

% ll=[        0         DincrT
%     NwaitMAX*TLEO    0.9*THEO];
ll=[ 0.0       (THEO-  DincrT)*0.0+DincrT      0           rLEO                   (THEO-  DincrT)*0.0+DincrT
     TLEO      (THEO-  DincrT)*1.0+DincrT    2*pi          2.5*aHEO               (THEO-  DincrT)*1.0+DincrT];
% Parametri di ricerca
J=[];
% C=[];
for i = 1:size(yav,1)
    ya=yav(i,:);
    y=(ll(2,:)-ll(1,:)).*ya+ll(1,:);
    [DvTOT,eccTransf,rTransfMIN] = DvTransf_3a(y,p);
    if isnan(DvTOT)
        DvTOT = realmax;%1.e14;
        eccTransf = realmax;%10;
        rTransfMIN = realmax;%0;
    end
    J=[J; DvTOT, y(1,1)+y(1,2)+y(1,5)];
    
%     C=[C; eccTransf, rTransfMIN];
end
if vuoto==1
save JV.dat -ascii -double J
% save C.dat -ascii -double C
else
    JV=J(1:2);
end