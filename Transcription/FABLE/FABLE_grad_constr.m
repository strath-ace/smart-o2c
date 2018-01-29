function [Ceq,t] = FABLE_grad_constr(eps, alpha, beta, L_prop, start_eq, check_eq, options, constants)

% =========================================================================
% Function for the computation of the gradients of the constraints in the
% multiple shooting and single shooting formulation of FABLE by finite
% difference.
% The input to this functions are modified by a delta value wrt the
% considered solution vector. The new equality constraints is computed
% (Ceq) and is then used in conjunction with the unperturbed Ceq (Ceq_old)
% to compute the gradient, as in:
% gradient = (Ceq - Ceq_old ) / delta
% 
% This function is used in FABLE_transcription_MS.m and
% FABLE_transcription_SS.m.
% 
% Input: eps      -> low-thrust acceleration
%        alpha    -> azimuth angle
%        beta     -> elevation angle
%        L_prop   -> true longitude values for porpagation
%        start_eq -> initial equinoctial elements
%        check_eq -> desired final equinoctial elements (node MS)
%        options  -> structure containing the options for the problem
%        constants-> structure of constants
% 
% Output: Ceq -> equality constraint
%         t   ->
% 
% =========================================================================
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk


% Select propagator based on user-defined model: constant acceleration or
% acceleration changing as 1/r2
if options.fun == 0
    [end_eq, t] = AnEquin_forward_m(L_prop, start_eq, eps, alpha, beta, constants.mu);
elseif options.fun == 1
    [end_eq, t] = AnEquin_forward_m_r2(L_prop, start_eq, eps, alpha, beta, constants.mu);
end

% The final point of the initial coast arc has to coincide with the initial point
% of the first low thrust arc
Ceq = end_eq(1:5,end) - check_eq;












