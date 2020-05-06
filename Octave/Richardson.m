## Copyright (C) 2019 Windows
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} Richardson (@var{a},@var{r},@var{p},@var{f})
##
## See:
## 1 - Annu. Rev. Fluid Mech., 29:123, 1997.
## 2 - AIAA J., 41:595, 2003.
##
## @seealso{}
## @end deftypefn

## Author: Otavio Beruski
## Created: 2019-10-11

function retval = Richardson (a,r,p,f1,f2)
% Calculates the generalized Richardson extrapolation for a variable @var{a},
% with order of convergence @var{p}, based on either a fine (f1 = 1) or a coarse
% grid (f1 = 2), related by the refinement factor @var{r}. The first value in
% @var{a} is that of the fine mesh, with the second being of the coarse mesh.
% The second flag f2 stipulates if the coefficient g and p are to be returned
% along with the extrapolate.

eps_21 = a(2) - a(1); % Difference between fine and coarse meshes.
R = 0.0; % Richardson extrapolate.
g = 0.0; % Series coefficient.

g = eps_21/(r^p - 1);

if f1 == 1
    R = a(1) + eps_21/(1 - r^p);
elseif f1 == 2
    R = a(2) + eps_21*r^p/(1 - r^p);
else
    display('Invalid value for flag. Check input');
    R = NaN;
endif

if f2 == 0
    retval = R;
elseif f2 == 1
    retval = [g;p;R];
endif
endfunction
