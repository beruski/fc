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
## @deftypefn {} {@var{V} =} X_corr (@var{X},@var{H},@var{Y})
##
## Bins the data in @var{X} according to @var{H}, attaching corresponding values
## of @var{Y} alongside it. @var{H} can be either a scalar, in which case the
## number of bins is given by @var{H}, or a vector, which contains the bin
## limits. The binned data, @var{V}, is returned as a cell where each entry is
## arranged as [@var{X} @var{Y}] for a given bin, where the @var{V} data are
## associated to the @var{X} data. The function expects that a given line of
## @var{Y} be associated to values at the same line in @var{X}.
##
## @seealso{YX_corr}
## @end deftypefn

## Author: Otavio Beruski
## Created: 2020-05-09

function V = X_corr (X,H,Y)

    % Variables and initialization
    npts = length(X);                       % Number of points in input array.
    ndim = size(Y,2);                       % Number of variables to correlate.
    if npts > 25000                         % Tentative upper limit on points
        nwish = 25000;                      % in a given bin.
    else
        nwish = npts;
    endif
    if length(H) == 1                       % Getting number of bins and bin
        nbins = H;                          % limits from input data.
        H = linspace(min(X),max(X),nbins)';
    else
        nbins = length(H);
    endif
    data = zeros(nwish,ndim+1);             % Binned data.
    k = 0;                                  % Counter for a given bin.
    for i = 1:(nbins-1)
        V{i} = 0.0;                         % Final data: cell with nbins-1
    endfor                                  % entries.
    clk = time;                             % CPU time for total time estimate.

    % Binning data
    for i = 2:nbins
        data(:,:) = 0.0;
        k = 1;
        for j = 1:npts
            if X(j) < H(i) && X(j) >= H(i-1)
                data(k,1) = X(j);
                data(k,2:(ndim+1)) = Y(j,:);
                k = k + 1;
            endif
        endfor
        V{i-1} = data(1:(k-1),:);
    endfor

    time-clk
endfunction
