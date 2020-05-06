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
## @deftypefn {} {@var{ret} =} surf_smooth (@var{S},@var{n})
##
## @seealso{}
## @end deftypefn

## Author: Otavio Beruski
## Created: 2018-09-13

function [ret] = surf_smooth (S,n)

[ny nx] = size(S);

sS = zeros(ny,nx);

%if n == 1
    sS(1,1) = (S(1,1) + S(2,1) + S(1,2))/3; % Smoothing corners
    sS(ny,1) = (S(ny,1) + S(ny-1,1) + S(ny,2))/3;
    sS(1,nx) = (S(1,nx) + S(2,nx) + S(1,nx-1))/3;
    sS(ny,nx) = (S(ny,nx) + S(ny-1,nx) + S(ny,nx-1))/3;
    for i = 2:(ny-1)
        sS(i,1) = 0.25*(S(i,1) + S(i-1,1) + S(i+1,1) + S(i,2)); % Smoothing 1st and last columns
        sS(i,nx) = 0.25*(S(i,nx) + S(i-1,nx) + S(i+1,nx) + S(i,nx-1));
        for j = 2:(nx-1) % Smoothing the bulk
            sS(i,j) = 0.2*(S(i,j) + S(i-1,j) + S(i+1,j) + S(i,j-1) + S(i,j+1));
        endfor
    endfor
    for j = 2:(nx-1) % Smoothing 1st and last rows
        sS(1,j) = 0.25*(S(1,j) + S(1,j-1) + S(1,j+1) + S(2,j));
        sS(ny,j) = 0.25*(S(ny,j) + S(ny,j-1) + S(ny,j+1) + S(ny-1,j));
    endfor
%endif

ret = sS;
endfunction
