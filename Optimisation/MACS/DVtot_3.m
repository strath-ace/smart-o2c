function [DVtot]=DVtot_3(the1,Dtheta1,the2,aTransfa,eccTransfa,ome0Transfa,aTransfb,eccTransfb,ome0Transfb)


orb_elem_sat_iniz = [7000e3,  0, 0*pi/180,  the1, 0*pi/180, 0*pi/180];

orb_elem_sat_iniztra1 = [aTransfa,  eccTransfa, 0*pi/180,  the1, 0*pi/180, ome0Transfa];
orb_elem_sat_fintra1 = [aTransfa,  eccTransfa, 0*pi/180,  the1+Dtheta1, 0*pi/180, ome0Transfa];

orb_elem_sat_iniztra2 = [aTransfb,  eccTransfb, 0*pi/180,  the1+Dtheta1, 0*pi/180, ome0Transfb];
orb_elem_sat_fintra2 = [aTransfb,  eccTransfb, 0*pi/180,  the2, 0*pi/180, ome0Transfb];

orb_elem_sat_fin =  [ 42164e3,  0, 0*pi/180, the2, 0*pi/180, 0*pi/180];

planet = 3.986e14; %earth
%planet = 4.2828e13; %  4.2828e13; %mars
mu = planet(1);

rv1 = oe2rv(orb_elem_sat_iniz,planet);

rvit1=oe2rv(orb_elem_sat_iniztra1,planet);
rvft1=oe2rv(orb_elem_sat_fintra1,planet);

rvit2=oe2rv(orb_elem_sat_iniztra2,planet);
rvft2=oe2rv(orb_elem_sat_fintra2,planet);

rv2 = oe2rv(orb_elem_sat_fin,planet);

% [rv1 rvit1 rvft1 rvit2 rvft2 rv2];
DVtot=norm(rvit1(4:6,1)-rv1(4:6,1))+norm(rvit2(4:6,1)-rvft1(4:6,1))+norm(rv2(4:6,1)-rvft2(4:6,1));