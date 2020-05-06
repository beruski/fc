## Copyright (C) 2019 Otavio Beruski
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} gci (@var{f},@var{r},@var{p})
##
##
## @seealso{}
## @end deftypefn

## Author: Otavio Beruski
## Created: 2019-08-20

function [retval] = gci (f,r,p)
% Estimates the Grid Convergence Index for f with grid refinement factors r
% and accuracy order p over two or three grids.

    n = size(f,2); % Number of grids used, influences calculations

    if n == 2
        gci = 0.0; % Grid convergence index
        Fs = 3; % Safety factor in error estimation
        eps = (f(2) - f(1))/f(1); % Difference in numerical solutions
        gci = Fs*abs(eps)/(r^p - 1);
        retval = gci; % Returns only the estimate GCI
    elseif n == 3
        gci12 = 0.0; % GCI between finer and mean grids
        gci23 = 0.0; % GCI between mean and coarser grids
        Fs = 1.25;
        eps12 = (f(2) - f(1))/f(1);
        eps23 = (f(3) - f(2))/f(2);
        gci12 = Fs*abs(eps12)/(r(1)^p - 1);
        gci23 = Fs*abs(eps23)/(r(2)^p - 1);
        R = gci23/(gci12*r(1)^p);
        retval = [gci12 gci23 R]; % Returns both GCIs and the asymptotic range
    end
endfunction
