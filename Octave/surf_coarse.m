## Copyright (C) 2018 Otavio Beruski
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
## @deftypefn {} {@var{ret} =} surf_coarse (@var{S},@var{n})
##
## @seealso{}
## @end deftypefn

## Author: Otavio Beruski
## Created: 2018-09-13

function [ret] = surf_coarse (S,n)

[ny nx] = size(S);
Nx = nx/n;
Ny = ny/n;

cS = zeros(Ny,Nx);
aux = zeros(ny,Nx);

for i = 1:ny % Coarsening along the columns
    j = 1;
    while j <= Nx
        for k = 0:n-1
            aux(i,j) += S(i,(j-1)*n+1+k);
        endfor
        %aux(i,j) = aux(i,j)/n; % Uncomment if the average is desired
        j += 1;
    endwhile
endfor
for j = 1:Nx % Coarsening along the rows
    i = 1;
    while i <= Ny
        for k = 0:n-1
            cS(i,j) += aux((i-1)*n+1+k,j);
        endfor
        %cS(i,j) = cS(i,j)/n; % Uncomment if the average is desired
        i += 1;
    endwhile
endfor

ret = cS;
endfunction
