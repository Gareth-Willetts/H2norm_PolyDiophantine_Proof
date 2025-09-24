function a = get_a_mimo(A)
%   GET_A_MIMO
%
% @brief:
%   This function implements Function 4 from (ACM paper citation) and
%   finds the least common multiple of a collection of
%   polynomials a_11(s), a_12(s), ..., a_{m1, m2}(s) stored in a cell array.
%
%   The resulting polynomial a(s) will have all its roots in the open left
%   half-plane if all input polynomials also have this property, as the
%   LCM preserves the set of roots of the inputs.
%
%   **Dependency:** This function requires the **Symbolic Math Toolbox**.
%
% Input:
%   A - An m1-by-m2 cell array, where each cell A{i, j} contains a row
%       vector of coefficients for a single polynomial.
%
% Output:
%   a - A row vector of coefficients for the resulting LCM polynomial.
%
% Example:
%   % Create a 2x2 cell array of polynomials:
%   % p11(s) = s + 1
%   % p12(s) = s + 2
%   % p21(s) = s^2 + 4s + 3  (which is (s+1)(s+3))
%   % p22(s) = 1 (a constant)
%   A = {[1 1], [1 2]; [1 4 3], [1]};
%
%   % The LCM is (s+1)(s+2)(s+3) = s^3 + 6s^2 + 11s + 6.
%   % The expected output is the coefficient vector [1 6 11 6].
%
%   a = get_a_mimo(A);
%

    % Define a symbolic variable 's' for polynomial representation.
    syms s;

    % Get the dimensions (m1, m2) from the input cell array A.
    [m1, m2] = size(A);

    % Initialize the result as a symbolic polynomial with the value 1.
    a_sym = sym(1);

    % Iterate through each polynomial in the input matrix A.
    for i = 1:m1
        for j = 1:m2
            % Convert the current numeric coefficient vector to a symbolic polynomial.
            current_poly_sym = poly2sym(A{i, j}, s);
            
            % Update the result by computing the LCM with the current polynomial.
            a_sym = lcm(a_sym, current_poly_sym);
        end
    end

    % Convert the final symbolic LCM polynomial back to a numeric coefficient vector.
    a = sym2poly(a_sym);
    
end