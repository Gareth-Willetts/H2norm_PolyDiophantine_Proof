function z0 = get_z0_mimo(A, C, a)
%   GET_Z0_MIMO
%
% @brief:
%   This function implements Function 5 from (ACM paper citation) and
%   calculates the polynomial z0(s) based on a multi-input,
%   multi-output (MIMO) system defined by matrices of transfer function
%   polynomials.
%
%   The calculation involves polynomial division, multiplication, and a
%   decomposition of intermediate polynomials into parts with even and
%   odd-powered coefficients.
%
% Input:
%   A - An m1-by-m2 cell array. A{i,j} contains the coefficient vector
%       for the denominator polynomial a_ij(s).
%   C - An m1-by-m2 cell array. C{i,j} contains the coefficient vector
%       for the numerator polynomial c_ij(s).
%   a - A row vector of coefficients for the polynomial a(s), which is a
%       common multiple of all polynomials in A.
%
% Output:
%   z0 - A row vector of coefficients for the resulting polynomial z0(s).
%

    % Get the dimensions (m1, m2) from the input cell array A.
    [m1, m2] = size(A);

    % --- Step 1: Calculate Â and Ĉ matrices ---
    % Initialize cell arrays to hold the new polynomials.
    A_hat = cell(m1, m2);
    C_hat = cell(m1, m2);

    for i = 1:m1
        for j = 1:m2
            % Â_{i,j} = a / a_{i,j}
            % deconv performs polynomial division. We assume the remainder is
            % negligible since 'a' is a common multiple.
            A_hat{i,j} = deconv(a, A{i,j});

            % Ĉ_{i,j} = c_{i,j} * Â_{i,j}
            % conv performs polynomial multiplication.
            C_hat{i,j} = conv(C{i,j}, A_hat{i,j});
        end
    end

    % --- Step 2: Decompose Ĉ into even and odd parts (Ĉe, Ĉo) ---
    C_hat_e = cell(m1, m2);
    C_hat_o = cell(m1, m2);

    for i = 1:m1
        for j = 1:m2
            % Decompose Ĉ_{i,j}(s) into Ĉe_{i,j}(s) and Ĉo_{i,j}(s) where:
            % Ĉe(s) = ĉ_0*s^0 + ĉ_2*s^1 + ĉ_4*s^2 + ...
            % Ĉo(s) = ĉ_1*s^0 + ĉ_3*s^1 + ĉ_5*s^2 + ...
            
            % This is done efficiently by slicing the reversed coefficient vector.
            c_hat_rev = fliplr(C_hat{i,j});
            
            % Extract coefficients for the even-part polynomial.
            C_hat_e{i,j} = fliplr(c_hat_rev(1:2:end));
            
            % Extract coefficients for the odd-part polynomial.
            co_coeffs = fliplr(c_hat_rev(2:2:end));
            
            % Handle cases where there are no odd-powered terms.
            if isempty(co_coeffs)
                C_hat_o{i,j} = 0;
            else
                C_hat_o{i,j} = co_coeffs;
            end
        end
    end

    % --- Step 3: Calculate the final polynomial z0 ---
    % Initialize z0(s) = 0.
    z0 = 0;

    for i = 1:m1
        for j = 1:m2
            % Calculate the term for this element: (Ĉe)² - s * (Ĉo)²
            Ce_sq = conv(C_hat_e{i,j}, C_hat_e{i,j});
            Co_sq = conv(C_hat_o{i,j}, C_hat_o{i,j});
            s_Co_sq = [Co_sq, 0]; % Multiplication by 's'

            % Pad polynomials with leading zeros to make their lengths equal
            % before subtracting.
            len1 = length(Ce_sq);
            len2 = length(s_Co_sq);
            max_len_term = max(len1, len2);
            term_ij = [zeros(1, max_len_term - len1), Ce_sq] - ...
                      [zeros(1, max_len_term - len2), s_Co_sq];

            % Add the resulting term to the running total z0.
            % This also requires padding to handle polynomials of different degrees.
            len_z = length(z0);
            len_t = length(term_ij);
            max_len_sum = max(len_z, len_t);
            z0 = [zeros(1, max_len_sum - len_z), z0] + ...
                 [zeros(1, max_len_sum - len_t), term_ij];
        end
    end
end