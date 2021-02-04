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
## @deftypefn {} {@var{V} =} volCorr (@var{X},@var{bins},@var{F})
##
## Correlates the 3D coordinates matrix @var{X} with the fields @var{F},
## following the sectioning of the volume given by @var{bins}.
##
## The coordinates matrix @var{X} is an Nx3 array containing the coordinates 
## of the discretized volumetric domain where @var{F} is evaluated. @var{bins}
## is a 3x1 cell array, where each cell element is an array containing the
## edges of sections, or bins, for each spatial coordinate. And @var{F} is an
## NxM array where each line contains the M fields corresponding to each 3D
## point (i.e. line) in @var{X}. 

## The output @var{V} is a (lx)x(ly)x(lz) cell array with dimensions given by
## the length of bins given in @var{bins}. Each cell of @var{V} is an
## (nz[i,j,k])xM array contaning the nz points of the i-th bin of coordinate
## var{X}(1), j-th bin of the coordinate @var{X}(2), and k-th bin of coordinate
## @var{X}(3); of the M fields in @var{F} belonging to the [i,j,k] section of
## the discretized volumetric domain.
##
## @seealso{X_corr,YX_corr}
## @end deftypefn

## Author: Otavio Beruski
## Created: 2020-12-11

function V = volCorr (X,bins,F)

N = size(X,1);          % Total number of points to be worked on
M = size(F,2);          % Number of fields defined over X
xBin = bins{1};         % Edges of bins in coordinate X(1), e.g. X
yBin = bins{2};         % Edges of bins in coordinate X(2), e.g. Y
zBin = bins{3};         % Edges of bins in coordinate X(3), e.g. Z
lx = length(xBin)-1;    % Number of bins in coordinate x
ly = length(yBin)-1;    % Number of bins in coordinate y
lz = length(zBin)-1;    % Number of bins in coordinate z

V = cell(lx,ly,lz);     % Points of F associated with each 3D bin of X

nx = 0;                 % Number of points in a given x bin
ny = 0;                 % Number of points in a given y bin
nz = 0;                 % Number of points in a given z bin
xIdx = [];              % Indeces of points belonging to a given x bin
yIdx = [];              % Indeces of points belonging to a given y bin
zIdx = [];              % Indeces of points belonging to a given z bin

p = 1;                  % Counter related to x
q = 1;                  % Counter related to y
r = 1;                  % Counter related to z
xCorr = [];             % Array containing X(2:3) and F in the i-th bin
yCorr = [];             % Array containing X(3) and F in the [i,j]-th bin
zCorr = [];             % Array containing F in the [i,j,k]-th bin

for i = 1:lx    % Loop in coordinate X(1)
    [nx xIdx] = histc(X(:,1),[xBin(i) xBin(i+1)]);  % Counting and labelling
    xCorr = zeros(nx(1),M+2);                       % points in the i-th bin
    p = 1;
    for j = 1:N
        if xIdx(j) == 1     % From histc(), means the point belongs to the bin
            xCorr(p,:) = [X(j,2:3) F(j,:)];
            p = p + 1;
        endif
    endfor
    for j = 1:ly    % Loop in coordinate X(2)
        [ny yIdx] = histc(xCorr(:,1),[yBin(j) yBin(j+1)]); % Counting and
        yCorr = zeros(ny(1),M+1);                          % labelling points
        q = 1;                                             % in [i,j]-th bin
        for k = 1:nx
            if yIdx(k) == 1
                yCorr(q,:) = xCorr(k,2:(M+2));
                q = q + 1;
            endif
        endfor
        for k = 1:lz    % Loop in coordinate X(3)
            [nz zIdx] = histc(yCorr(:,1),[zBin(k) zBin(k+1)]); % Counting an
            zCorr = zeros(nz(1),M);                            % labelling
            r = 1;                                             % points in 
            for m = 1:ny                                       % [i,j,k]-th bin
                if zIdx(m) == 1
                    zCorr(r,:) = yCorr(m,2:(M+1));
                    r = r + 1;
                endif
            endfor
            V{i,j,k} = zCorr;   % Values of F corresponding to the [i,j,k]-th
        endfor                  % bin of the volumetric domain.
    endfor
endfor

endfunction