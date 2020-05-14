## Copyright (C) 2020 Otavio Beruski
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
## @deftypefn {} {[@var{ret2} @var{ret1}] =} ash (@var{y},@var{t},@var{h},@var{m})
##
## Estimates the probability distribution of a dataset @var{y} using the
## Averaged Shifted Histogram (ASH) method. The data is analyzed between the
## intervals set by @var{t}: [t(1) t(2)), with @var{h} evenly spaced intervals. 
## The number of shifts and histograms to be averaged over is given by @var{m},
## with the weights given simply by the triangular kernel (see ref. [1]).
## The function returns the unnormalized PDF estimate @var{f} and optionally
## the limits of the bins used.
##
## [1] Scott, D.W., Averaged shifted histogram. WIREs Comp. Stat., 2:160-164,
## 2010. Available at:
## https://www.researchgate.net/publication/229760716_Averaged_Shifted_Histogram
## 
## @seealso{histc}
## @end deftypefn

## Author: Otavio Beruski
## Created: 2020-04-27

function [ret1 ret2] = ash (y,t,h,m)

    L = length(y);                      % Length of the data set
    uy = t(2);                          % Upper range limit
    ly = t(1);                          % Lower range limit
    d = h/m;                            % Increment in bin size
    f = zeros(ceil(m*(uy-ly)/h),1);     % Final frequency array
    nu = zeros(2,1);                    % Frequency for each h/m subinterval
    nu0 = 0;                            % Averaged shifted frequency in a bin

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

    ret1 = f;
    ret2 = linspace(ly,uy,length(f))';
endfunction
