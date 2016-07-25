%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERAL DESCRIPTION : 10 Objective functions (2+5+3)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

There are three folders : DED Codes, ELD Codes, Hydrothermal Codes. 
Each folder contain an MS Word document containing the problem formulation.

Each folder contains objective functions (function scripts) which always start with "fn_" in the name.
Each objective function is capable of operating on a row vector of values and returns the first output 
as a single value which indicates the cost (inclusive of penalty added). The objective function can be called
with an optional second argument upon which the cost and penalty details are displayed.

Each folder contains a function script with name "MY_FUNCTION" which evaluates a matrix (a population or a collection of
row vectors to provide). It returns a column vector of objective values (inclusive of penalty obviously) and also the 
count of function evaluations done in that particular call. Each value in the population is expected to have a precision 
upto four decimal digits only. Any extra decimal digits need to be rounded off. The optimization algorithm is supposed to
pass its population to "MY_FUNCTION" and retrieve the cost values.

Each folder contains a Solution Script which has the data of lower limits and upper limits of the row vectors to be passed
to "MY_FUNCTION". That corresponds to the search space boundaries of the decision variables.

--------------------------------------------------------------------------------------------------------------------------

1)   Readme.txt

DESCRIPTION : This very file that is currently being looked at.

--------------------------------------------------------------------------------------------------------------------------

2)   DED Codes\MY_FUNCTION.m

DESCRIPTION : Takes input a matrix and sends each row to the appropriate
objective function mentioned within. Outputs are a column matrix of 
objective values and the count of rows evaluated. Optimization algorithm is
expected to call this function.

**************************************************************************

3)   DED Codes\fn_DED_5.m

DESCRIPTION : The objective function for evaluating a row vector of (24*5)
values.

**************************************************************************

4)   DED Codes\fn_DED_10.m

DESCRIPTION : The objective function for evaluating a row vector of (24*10)
values.

**************************************************************************

5)   DED Codes\My_Solutions.m

DESCRIPTION : The lower limits and upper limits for the row vectors of the
population. 

--------------------------------------------------------------------------------------------------------------------------

6)   ELD Codes\MY_FUNCTION.m

DESCRIPTION : Takes input a matrix and sends each row to the appropriate
objective function mentioned within. Outputs are a column matrix of 
objective values and the count of rows evaluated. Optimization algorithm is
expected to call this function.

**************************************************************************

7)   ELD Codes\fn_ELD_6.m

DESCRIPTION : The objective function for evaluating a row vector of 6 values.

**************************************************************************

8)   ELD Codes\fn_ELD_13.m

DESCRIPTION : The objective function for evaluating a row vector of 13 values.

**************************************************************************

9)   ELD Codes\fn_ELD_15.m

DESCRIPTION : The objective function for evaluating a row vector of 15 values.

**************************************************************************

10)   ELD Codes\fn_ELD_40.m

DESCRIPTION : The objective function for evaluating a row vector of 40 values.

**************************************************************************

11)   ELD Codes\fn_ELD_140.m

DESCRIPTION : The objective function for evaluating a row vector of 140 values.

**************************************************************************

12)   ELD Codes\My_Solutions.m

DESCRIPTION : The lower limits and upper limits for the row vectors of the
population. 

--------------------------------------------------------------------------------------------------------------------------

13)   Hydrothermal Code\MY_FUNCTION.m

DESCRIPTION : Takes input a matrix and sends each row to the appropriate
objective function mentioned within. Outputs are a column matrix of 
objective values and the count of rows evaluated. Optimization algorithm is
expected to call this function.

**************************************************************************

14)   Hydrothermal Code\fn_HT_ELD_Case_1.m

DESCRIPTION : The objective function for evaluating a row vector of (24*4)
values.

**************************************************************************

15)   Hydrothermal Code\fn_HT_ELD_Case_2.m

DESCRIPTION : The objective function for evaluating a row vector of (24*4)
values.

**************************************************************************

16)   Hydrothermal Code\fn_HT_ELD_Case_3.m

DESCRIPTION : The objective function for evaluating a row vector of (24*4)
values.

**************************************************************************

17)   Hydrothermal Code\My_Solutions.m

DESCRIPTION : The lower limits and upper limits for the row vectors of the
population.

--------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------