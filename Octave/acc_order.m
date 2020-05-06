## Copyright (C) 2019 Otavio Beruski
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
## @deftypefn {} {@var{retval} =} acc_order (@var{f1},@var{f2},@var{f3},@var{r12},@var{r23})
##
## Estimates the accuracy order of the numerical solution to a given variable
## f as a function of grid spacing, using the generalized Richardson
## extrapolation as framework. The concept is summarized in the following
## reference:
##
## Annu. Rev. Fluid. Mech., 29:123-160, 1997.
##
## Three grids are required, defined as 1, 2, and 3 from finer to coarser. The
## refinement ratios r_ij are defined as
##
## r_ij = h_j/h_i
##
## where h is a measurement of grid spacing. In the case of unstrucutred grids,
## the following may be used:
##
## r_ij = (N_i/N_j)^(1/D)
##
## where N is the total number of elements and D is the dimensionality. The
## accuracy order is estimated using the generalized Richardson extrapolate as
## a framework, admitting r_ij non-integer and non-constant, through the
## implicit equation:
##
## eps_23/(r_23^p - 1) = r_12^p*[eps_12/(r_12^p - 1)]
##
## where eps_ij is the difference in numerical solution of f between grids
## j and i, and p is the accuracy order. Estimation of p follows via
## substitution iteration using the equations:
##
## p = omega*rho + (1 - omega)*log(beta)/log(r_12)
##
## beta = (r_12^rho - 1)*eps_23/(r_23^rho - 1)*eps_12
##
## where omega is the relaxation factor and rho is the previously calculated
## accuracy order.
##
## @seealso{}
## @end deftypefn

## Author: Otavio Beruski
## Created: 2019-08-20

function [retval] = acc_order (a,r,w)
% Returns the estimated accuracy order for f, using three different grids, from
% coarser to finer: f(3), f(2) and f(1). The refinement factors r need not be
% integers or equal. The relaxation factor w controls the convergence of the
% iterative procedure.
    
    p = 2; % Estimated accuracy order
    
    eps12 = a(2) - a(1); % Difference between fine- and average-grid solutions
    eps23 = a(3) - a(2); % Differecen between average- and coarse-grid solutions
    
    %w = 0.5; % Relaxation factor
    beta = 0.0; % Function of eps and r, part of the iterative solver
    rho = 0; % Previously calculated p in the iterative solver
    
    iter = 0; % Number of iterations.
    
    while abs(p-rho) >= 1E-6
        rho = p;
        beta = (r(1)^rho - 1)*eps23/((r(2)^rho - 1)*eps12);
        p = w*rho + (1-w)*log(beta)/log(r(1));
        iter = iter + 1;
        if iter > 1E4
            display "Number of iterations too high. Try modifying w."
            retval = NaN;
            return
        end
    end
    
    retval = p;
endfunction
