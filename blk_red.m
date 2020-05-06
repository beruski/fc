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
## @deftypefn {} {@var{retval} =} blk_red (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Otavio Beruski
## Created: 2018-09-14

function [ret] = blk_red (A,m,avg)

n = length(A);
r = n/m;

Ap = zeros(n,r);
B = zeros(r,r);

for i = 1:n
    for j = 1:r
        for k = 1:m
            Ap(i,j) = Ap(i,j) + A(i,k+(j-1)*m);
        endfor
    endfor
endfor

for i = 1:r
    for j = 1:r
        for k = 1:m
            B(j,i) = B(j,i) + Ap(k+(j-1)*m,i);
        endfor
    endfor
endfor


if avg == 0
    ret = B;
else
    ret = B/m^2;
endif
endfunction
