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
## @deftypefn {} {[@var{Ybins} @var{YX}] =} YX_corr (@var{V},@var{c},@var{t},@var{m},@var{f_log},@var{f_norm})
##
## Given a cell generated by X_corr() with the structure
##              @var{V}{i} = [X ... Y ...],
## where Y is the @var{c} column of each i-th entry of @var{V}, estimates the
## probability density of Y for each bin of X previously estipulated. The
## output is a surface YX, where each column contains the PDF of Y for the
## respective bin of X in a given range of values
##              [@var{t}(1) @var{t}(2))
## and increment @var{h}, estimated through the Averaged Shifted Histogram
## method. Optionally, the function also returns the bin limits, Ybins, used
## to estimate the PDF via ASH.
##
## The parameters @var{t}, @var{h}, and @var{m} concerns the fucntion ash(),
## which is used to estimated the probability density function. The flags
## @var{f_log} and @var{f_norm} indicate if the binning space is logarithmic
## (f_log = 1) or linear, and if the final frequency surface YX is to be
## normalized (f_norm = 1) or not. As a matter of consistency, it was chosen to
## normalize over the entire surface and following the choice of binning space,
## e.g. if f_log = 1 and f_norm = 1, the normalization is performed over a
## logarithmic binning space.
##
## @seealso{X_corr,ash}
## @end deftypefn

## Author: Otavio Beruski
## Created: 2020-05-09

function [YX Ybins] = YX_corr (V,c,t,h,m,f_log,f_norm)

    nbins = size(V,2);                       % Follows from X_corr data.
    nYbins = length(ash(1,t,h,m));           % Extracting actual array length.
    if f_log == 1                            % Bins in Y:
        Ybins = logspace(t(1),t(2),nYbins)'; % logarithmic
    else                                     % or
        Ybins = linspace(t(1),t(2),nYbins)'; % linear.
    endif
    YX = zeros(nYbins,nbins);                % Final array.
    YXsum = 0.0;                             % Cumulative integral of YX.
    clk = time;                              % CPU time for total time estimate.

    for i = 1:nbins
        % Taking requested data out of the cell.
        if length(V{i}) == 0 % In case there's no entry for a given bin.
            continue
        endif
        Y = zeros(length(V{i}),1); % I really prefer to initialize it first.
        Y(:) = V{i}(:,c);
        % Calculating ASHs
        if f_log == 1 % For the choice of logarithmic binning space.
            YX(:,i) = ash(log10(Y),t,h,m);
        else % For a linear binning space.
            YX(:,i) = ash(Y,t,h,m);
        endif
    endfor
    % Normalization
    if f_norm == 1
        if f_log == 1 % Carrying over the choice of a logarithmic binnign space.
            for i = 1:nbins
                YXsum = YXsum + trapz(log10(Ybins),YX(:,i));
            endfor
        else % For linear binning space.
            for i = 1:nbins
                YXsum = YXsum + trapz(Ybins,YX(:,i));
            endfor
        endif
        YX = YX/YXsum;
    endif

    time - clk
endfunction