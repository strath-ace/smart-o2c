function [memories,dd,arch_energy]=add_arch(memories,candidates,dd,arch_energy,lx,mfit,max_arch)

% Addition of elements to archive with optimal disposition
%
% [memories,dd,arc_energy]=add_arch(memories,dd,arch_enerfy,lx,mfit,max_arch)
%
% INPUT
%       memories        :       archive
%       candidates      :       potential new elements of the archive
%       dd              :       pairwise distance matrix
%       arch_energy     :       archive energy
%       lx              :       parameter space dimensions
%       mfit            :       objective function dimensions
%       max_arch        :       max allowed size
%
% OUTPUT
%       memories        :       updated archive
%       dd              :       updated pairwise distance matrix
%       arch_energy     :       updeted archive energy
%
% Written by Lorenzo A. Ricciardi 2015
%
% Assuming memories is self dominated and not dominated by any candidate
% Assuming candidates are self dominated, not dominated by any memory and
% do not dominate any memory



return