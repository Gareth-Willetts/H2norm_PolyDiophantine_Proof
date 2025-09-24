function [pn_plus_1, zn_minus_1, an_out] = solve_fraction_free_H2_norm_mimo(A, C)
% SOLVE_FRACTION_FREE_H2_NORM_MIMO
%
% @brief:
%   This function implements the logic from Function 6. It implements
%   a fraction-free algorithm to find the terms needed to calculate the 
%   H2 norm of a multi-input, multi-output (MIMO) system. 
%   It returns values only if the system is stable (i.e., all
%   poles are in the open left half-plane).
%
%   The function first computes a common denominator 'a' for all transfer
%   functions, then computes a composite polynomial 'z0'. Finally, it
%   passes these intermediate polynomials to a core routine that determines
%   stability and calculates the final terms.
%
%   **Dependencies:** This function requires the following helper functions
%   to be on the MATLAB path:
%     - get_a_mimo.m
%     - get_p1_p2.m
%     - get_z0_mimo.m
%     - get_fraction_free_terms.m
%
% Input:
%   A - An m1-by-m2 cell array. A{i,j} contains the coefficient vector
%       for the denominator polynomial a_ij(s).
%   C - An m1-by-m2 cell array. C{i,j} contains the coefficient vector
%       for the numerator polynomial c_ij(s).
%
% Output:
%   pn_plus_1 - The term p_{n+1}. Returns [] if the stability test fails.
%   zn_minus_1 - The term z_{n-1}. Returns [] if the stability test fails.
%   an_out - The term a_n. Returns [] if the stability test fails.
%

    % a <- get_a_mimo(m1, m2, A)
    a = get_a_mimo(A);

    % n <- deg(a)
    n = length(a) - 1;

    % (p1, p2) <- get_p1_p2(n, a)
    [p1, p2] = get_p1_p2(n, a);

    % z0 <- get_z0_mimo(m1, m2, A, C, a, n)
    z0 = get_z0_mimo(A, C, a);

    % d1 <- deg(p1), d2 <- deg(p2)
    d1 = length(p1) - 1;
    d2 = length(p2) - 1;

    % an <- coeff(p1, d1)
    an = p1(1); % Leading coefficient of p1

    % s <- 1
    s_flag = 1;

    % if an * coeff(p2, d2) <= 0 then s <- 0 end
    % This is an initial stability check. If p2 is zero or the leading
    % coefficients have different signs (or one is zero), the roots of 'a'
    % cannot all be in the open left half-plane.
    if isempty(p2) || p2(1) == 0 || an * p2(1) <= 0
        s_flag = 0;
    end

    % (s, pn+1, zn-1, an) <- get_fraction_free_terms(...)
    [s_flag, pn_plus_1_val, zn_minus_1_val, an_val] = ...
        get_fraction_free_terms(p1, p2, z0, n, an, d1, d2, s_flag);

    % if s then return pn+1, zn-1, an
    if s_flag
        pn_plus_1 = pn_plus_1_val;
        zn_minus_1 = zn_minus_1_val;
        an_out = an_val;
    else
        % Return nothing (empty arrays in MATLAB) if the stability test failed.
        pn_plus_1 = [];
        zn_minus_1 = [];
        an_out = [];
    end
end