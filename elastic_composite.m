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
## @deftypefn {} {@var{ret1},@var{ret2} =} elastic_composite (@var{K1},@var{G1},@var{K2},@var{G2},@var{e})
##
## Calculates the bulk modulus K = @var{ret1} and shear modulus G = @var{ret2}
## for a binary composite material with a volume fraction @var{e} of component
## 2. The inputs [@var{K1},@var{G1}] and [@var{K2},@var{G2}] are the bulk and
## shear modulus for components 1 and 2, respectively.
##
## The Self-consistent (SCS) approximation is used to calculate the elastic
## properties iteratively. The equations can be written in a grouped fashion:
##
## (M - M1)/(M2 - M1) = e/(1 + (1-e)*(M2 - M1)/(M1 + F))
##
## where M = K,G and F = F(K1,G1,K2,G2,K,G). For the SCS approximation for a
## binary composite of spheres, F is given by:
##
## for K: F = 4*G/3
##
## for G: F = G*(9*K + 8*G)/(6*(K + 2*G))
##
## Since the equations for K and G are coupled, an iterative procedure is used
## to solve them. In this case, the Hashin-Shtrikman model is used as a first
## guess, in particular the lower bounds established by the model:
##
## if (G2 - G1)*(K2 - K1) >= 0
##
##     for K: F = 4*G1/3
##     for G: F = G1*(9*K1 + 8*K1)/(6*(K1 + 2*G1))
##
## else if (G2 - G1)*(K2 - K1) < 0
##
##     for K: F = 4*G2/3
##     for G: F = G2*(9*K2 + 8*K2)/(6*(K2 + 2*G2))
##
## In the algorithm below, K is first estimated using the appropriate guess,
## and then used in the first estimate G. Two special cases can be pointed out:
##
## If one component is "void": K2 = G2 = 0
##
## If one component is a fluid: K2 != 0, G2 = 0
##
## An appropriate citation for this script would be:
##
## Watt, J. P., Davies, G. F., O'Connel, R. J. The Elastic Properties of
## Composite Materials. Reviews of Geophysics and Space Physics, 14(4), 1976.
##
## See also the reference above for a proper discussion on the models and its
## limitations and special cases.
##
## @seealso{}
## @end deftypefn

## Author: Otavio Beruski
## Created: 2018-08-09

function [ret1 ret2] = elastic_composite (K1,G1,K2,G2,e)

K = 0; # Composite bulk modulus.
G = 0; # Composite shear modulus.
F = 0; # Auxiliary function. See description above.

K0 = 0; # Reference values.
G0 = 0;

THRESH = 1e-6; # Threshold for convergence.

if (G2-G1)*(K2-K1) >= 0
    K = K1 + e*(K2-K1)/(1 + (1-e)*(K2-K1)/(K1+4*G1/3));
    F = G1*(9*K+8*G1)/(6*(K+2*G1));
    G = G1 + e*(G2-G1)/(1 + (1-e)*(G2-G1)/(G1 + F));
else
    K = K1 + e*(K2-K1)/(1 + (1-e)*(K2-K1)/(K1+4*G2/3));
    F = G2*(9*K+8*G2)/(6*(K+2*G2));
    G = G1 + e*(G2-G1)/(1 + (1-e)*(G2-G1)/(G1 + F));
endif

while abs(1-K0/K) >= THRESH && abs(1-G0/G) >= THRESH
    K0 = K;
    G0 = G;
    K = K1 + e*(K2-K1)/(1 + (1-e)*(K2-K1)/(K1+4*G/3));
    F = G*(9*K+8*G)/(6*(K+2*G));
    G = G1 + e*(G2-G1)/(1 + (1-e)*(G2-G1)/(G1 + F));
endwhile

ret1 = K;
ret2 = G;
endfunction
