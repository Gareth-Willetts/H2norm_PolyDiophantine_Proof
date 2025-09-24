function [s, pn_plus_1, zn_minus_1, an] = get_fraction_free_terms(p1, p2, z0, n, an, da, db, s)
% GET_FRACTION_FREE_TERMS
%
% @brief:
%   This function implements Function 1 from (ACM paper citation): a
%   fraction-free algorithm for computing p_{n+1} and z_{n-1} from
%   Definition 5.2, following the results that will be established in
%   Theorem 7.1.
%
% Inputs:
%   p1 - A row vector of coefficients for the first polynomial (p_i).
%   p2 - A row vector of coefficients for the second polynomial (p_{i+1}).
%   z0 - A row vector of coefficients for the initial z polynomial (z_{i-1}).
%   Note that these must be entered in descending power.
%
% Outputs:
%   pn_plus_1  - The final polynomial in the sequence of p polynomials (should be a constant).
%   zn_minus_1 - The final polynomial in the sequence of z polynomials (should be a constant).
%   an         - The leading coefficient of the input polynomial p2 (possibly incorrect).

% --- Helper Function to Subtract Polynomials ---
function p_out = poly_subtract(poly1, poly2)
    % Subtracts poly2 from poly1.
    % Padding is done in the main while loop.
    len1 = length(poly1);
    len2 = length(poly2);
    if len1 == len2
        p_out = poly1 - poly2;
    else
        s = 0;
        error('Mismatched polynomial lengths.');
    end
    
    % Remove leading zeros from the result to keep polynomials tidy.
    first_nonzero = find(p_out ~= 0, 1, 'first');
    if isempty(first_nonzero)
        p_out = 0; % The polynomial is the scalar zero.
    else
        p_out = p_out(first_nonzero:end);
    end
end

% ---
% ## Initialization
% ---
i = 1;
pa = p1;  % p_i in Theorem 7.1
pb = p2;  % p_{i+1} in Theorem 7.1
z = z0;   % z_{i-1} in Theorem 7.1
v = 1;    % Initialized to 1

% q determines whether to multiply by s when calculating next p polynomial
q = mod(n, 2);

% ---
% ## Main Algorithm Loop
% ---

while (i <= n-1) && s

    % --- Update z polynomial ---
    % Align pa's degree to z's by multiplying by s^(floor((n-i-1)/2)). 
    % In MATLAB, this is done by appending the appropriate number of zeros 
    % to the coefficient vector.
    pa_aligned = [pa, zeros(1, floor((n-i-1)/2))];

    % This subtraction now correctly reduces the degree of z by one.
    z = poly_subtract(pa(1) * z, z(1) * pa_aligned);

    % --- Update p polynomial ---
    % If q=1, multiply pb by s to align degrees before subtracting.
    % In MATLAB, this is done by appending a zero to the coefficient
    % vector.
    pc = poly_subtract(pb(1) * pa, pa(1) * [pb zeros(1, q)]);
    
    % --- Division step ---
    if i >= 2
        % From the third iteration onwards, divide by the previous v.
        if i > 2
            z = z / v;
            pc = pc / v;
        end
        % Update v for the next division step.
        v = pa(1);
    end

    % --- Update variables for the next iteration ---
    pa = pb;
    pb = pc;
    i  = i + 1;
    q  = ~q; % Toggle q between 0 and 1
    da = db;
    if q
        db = db - 1;
    end
    
    % --- Stability Test ---
    % Checks the signs of the leading coefficients of the p polynomials.
    if (mod(i, 2) == 0 && pb(1) <= 0) || ...
       (mod(i, 2) == 1 && an * pb(1) <= 0)
        s = 0;
    end
end

% ---
% ## Return Results
% ---
if s
    % If stable, return the final constants.
    pn_plus_1  = pb;
    zn_minus_1 = z;
else
    % If the stability test failed, return empty arrays and throw an error.
    pn_plus_1  = [];
    zn_minus_1 = [];
    an         = [];
    error('System stability test failed during algorithm execution.');
end

end