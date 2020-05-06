## Copyright (C) 2020 Windows
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
## @deftypefn {} {@var{retval} =} ash (@var{y},@var{t},@var{h},@var{m})
##
## @seealso{histc}
## @end deftypefn

## Author: Otavio Beruski
## Created: 2020-04-27

function [retval] = ash (y,t,h,m)

    L = length(y);
    uy = t(2);
    ly = t(1);
    d = h/m;
    f = zeros(ceil(m*(uy-ly)/h),1);
    nu = zeros(2,1);
    nu0 = 0;

    k = 1;
    x = ly;
    while x < uy
        nu0 = 0;
        for i = 1:m
            lx = x + (i-m)*d;
            ux = x + i*d;
            nu = histc(y,[lx,ux]);
            nu0 = nu0 + nu(1);
        endfor
        f(k) = nu0/m;
        k = k + 1;
        x = x + d;
    endwhile

    retval = f;
endfunction
