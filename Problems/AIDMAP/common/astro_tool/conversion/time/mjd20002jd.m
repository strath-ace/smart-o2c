function jd = mjd20002jd(mjd2000)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
<<<<<<< HEAD
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------
=======
%-----------Copyright (C) 2016 University of Strathclyde-------------
>>>>>>> 5b7361d93c9119cf1d2e9e6c885bed93f924d71b
%
%
%
% MJD20002JD Julian day number from Modified Julian day 2000 number.
% 
% JD = MJD20002JD(MJD2000) returns the Julian day number corresponding to
% the given modified Julian day 2000 number.
% 
% INPUT
% 	mjd2000     Date in MJD 2000. MJD2000 is defined as the number of days
%              since 01-01-2000, 12:00 noon. 
% 
% OUTPUT
%  jd          Date in Julian Day. The JD (Julian day) count is from 0 at
%              12:00 noon, 1 January -4712 (4713 BC), Julian proleptic
%              calendar. The corresponding date in Gregorian calendar is
%              12:00 noon, 24 November -4713.
% 
% See also JD2MJD2000.
% 
% FUNCTIONS CALLED
%  none
% 
% - Nicolas Croisard - 16/02/2008
% - Revised by Matteo Ceriotti - 20/02/2008
% 
% ------------------------- - SpaceART Toolbox - --------------------------


mjd = mjd2000 + 51544.5;
jd  = mjd + 2400000.5;


return
