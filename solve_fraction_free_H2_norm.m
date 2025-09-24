function [pn_plus_1, zn_minus_1, an_out] = solve_fraction_free_H2_norm(a, c)
% SOLVE_FRACTION_FREE_H2_NORM 
%
% @brief:
%   This function implements Function 3 to calculate terms
%   used in finding the square of the H2 norm of a transfer function
%   G(s) = c(s)/a(s). 
%   The algorithm is fraction-free. It returns values only if the input 
%   polynomial 'a' has all its roots in the open left half-plane 
%   (a stability condition that is checked within the helper function 
%   'get_fraction_free_terms').
%
%   **Dependencies:** This function requires the following files to be on the
%   MATLAB path:
%     - get_p1_p2.m
%     - get_fraction_free_terms.m
%
% Input:
%   a - A row vector of coefficients for the denominator polynomial a(s).
%   c - A row vector of coefficients for the numerator polynomial c(s).
%   These both need to be input in order of descending powers of s.
%
% Output:
%   pn_plus_1 - The term p_{n+1}. Returns [] if the stability test fails.
%   zn_minus_1 - The term z_{n-1}. Returns [] if the stability test fails.
%   an_out - The term a_n. Returns [] if the stability test fails.
%

    % n <- deg(a)
    n = length(a) - 1;
    [p1, p2] = get_p1_p2(n, a);

    % --- Construct polynomials ce and co from c(s) ---
    % The pseudocode creates ce(s) from even-indexed coefficients of c(s)
    % and co(s) from odd-indexed coefficients.
    % ce(s) = c_0*s^0 + c_2*s^1 + c_4*s^2 + ...
    % co(s) = c_1*s^0 + c_3*s^1 + c_5*s^2 + ...
    c_rev = fliplr(c);
    ce_coeffs_ascending = c_rev(1:2:end);
    co_coeffs_ascending = c_rev(2:2:end);
    
    ce = fliplr(ce_coeffs_ascending);
    co = fliplr(co_coeffs_ascending);
    
    % z0(s) <- (ce(s))^2 - s * (co(s))^2
    ce_squared = conv(ce, ce);
    co_squared = conv(co, co);
    s_times_co_squared = [co_squared, 0];
    
    % Pad with zeros to subtract polynomials of different lengths
    len1 = length(ce_squared);
    len2 = length(s_times_co_squared);
    max_len = max(len1, len2);
    
    ce_squared_padded = [zeros(1, max_len - len1), ce_squared];
    s_co_sq_padded = [zeros(1, max_len - len2), s_times_co_squared];
    
    z0 = ce_squared_padded - s_co_sq_padded;

    % d1 <- deg(p1), d2 <- deg(p2)
    d1 = length(p1) - 1;
    d2 = length(p2) - 1;
    
    % an <- coeff(p1, d1)
    an = p1(1); % Leading coefficient of p1

    % s <- 1
    s_flag = 1;
    
    % if an * coeff(p2, d2) <= 0 then s <- 0 end
    % This is a stability check. If p2 is zero or the leading coefficients
    % have different signs, the algorithm terminates.
    if isempty(p2) || p2(1) == 0 || an * p2(1) <= 0
        s_flag = 0;
    end
    
    % (s, pn+1, zn-1, an) <- get_fraction_free_terms(...)
    [s_flag, pn_plus_1_val, zn_minus_1_val, an_val] = ...
        get_fraction_free_terms(p1, p2, z0, n, an, d1, d2, s_flag);
    
    % if s then return pn+1, zn-1, an end
    if s_flag
        pn_plus_1 = pn_plus_1_val;
        zn_minus_1 = zn_minus_1_val;
        an_out = an_val;
    else
        % Return nothing (empty arrays in MATLAB)
        pn_plus_1 = [];
        zn_minus_1 = [];
        an_out = [];
    end
end