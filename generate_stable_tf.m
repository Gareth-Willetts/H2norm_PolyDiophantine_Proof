function [G, num_coeffs, den_coeffs] = generate_stable_tf(num_order, den_order)
% GENERATE_STABLE_TF
%
% @brief:
%   This function generates a random, proper, continuous-time transfer function
%   G(s) that is guaranteed to be stable. Stability is ensured by creating
%   poles (the roots of the denominator) that all have negative real parts.
%
%   The zeros of the transfer function are placed randomly in the complex
%   plane. The function ensures that any complex poles or zeros appear in
%   conjugate pairs to yield real-valued polynomial coefficients.
%
% Inputs:
%   num_order - The desired order of the numerator polynomial.
%   den_order - The desired order of the denominator polynomial. For a
%               proper transfer function, den_order must be >= num_order.
%
% Outputs:
%   G          - The generated transfer function as a 'tf' object.
%   num_coeffs - A row vector of the numerator coefficients.
%   den_coeffs - A row vector of the denominator coefficients.
%
% Example:
%   % Generate a 3rd-order stable system with one zero
%   [G, num, den] = generate_stable_tf(1, 3);
%   
%   % Plot the pole-zero map to verify stability
%   figure;
%   pzmap(G);
%   title('Pole-Zero Map of Random Stable System');
%   grid on;

% --- Input Validation ---
if num_order > den_order
    error('Numerator order cannot be greater than denominator order for a proper transfer function.');
end
if nargin < 2
    error('Both numerator and denominator order must be specified.');
end

% --- Generate Stable Poles (Roots of the Denominator) ---
% To ensure stability, all poles must have negative real parts.
poles = [];
i = 0;
while length(poles) < den_order
    % Randomly decide to add a real pole or a complex conjugate pair
    is_real_pole = rand() > 0.5;
    
    if is_real_pole || (den_order - length(poles) < 2)
        % Add a single real pole.
        % Generate a random number and ensure it's negative.
        p_real = -abs(10 * randn());
        poles = [poles; p_real];
    else
        % Add a complex conjugate pair.
        % Ensure the real part is negative for stability.
        p_complex_real = -abs(10 * randn());
        p_complex_imag = 10 * randn(); % Imaginary part can be anything
        
        poles = [poles; p_complex_real + 1j*p_complex_imag];
        poles = [poles; p_complex_real - 1j*p_complex_imag];
    end
end


% --- Generate Zeros (Roots of the Numerator) ---
% Zeros can be anywhere in the complex plane, so no stability constraints.
zeros_ = [];
i = 0;
while length(zeros_) < num_order
    % Randomly decide to add a real zero or a complex conjugate pair
    is_real_zero = rand() > 0.5;
    
    if is_real_zero || (num_order - length(zeros_) < 2)
        % Add a single real zero.
        z_real = 20 * rand() - 10; % Random value between -10 and 10
        zeros_ = [zeros_; z_real];
    else
        % Add a complex conjugate pair.
        z_complex_real = 20 * rand() - 10;
        z_complex_imag = 20 * rand() - 10;
        
        zeros_ = [zeros_; z_complex_real + 1j*z_complex_imag];
        zeros_ = [zeros_; z_complex_real - 1j*z_complex_imag];
    end
end


% --- Convert Roots to Polynomial Coefficients ---
% The poly() function converts a vector of roots into a polynomial.
% We take the real part to eliminate any tiny imaginary numerical errors.
num_coeffs = real(poly(zeros_));
den_coeffs = real(poly(poles));

% If the numerator is order 0 (a constant), poly() returns 1. 
% We can add a random gain in this case.
if num_order == 0
    num_coeffs = randi([1, 10]);
end


% --- Create the Transfer Function Object ---
G = tf(num_coeffs, den_coeffs);

end


% =========================================================================
% ==                  EXAMPLE USAGE SCRIPT                             ==
% =========================================================================
% 
% To use the function, save it as 'generate_stable_tf.m' in your MATLAB
% path. Then you can run the following script.

%{

% --- Settings ---
NUM_ORDER = 4;  % Desired numerator order
DEN_ORDER = 6;  % Desired denominator order

% --- Generate the System ---
[G, num_poly, den_poly] = generate_stable_tf(NUM_ORDER, DEN_ORDER);

% --- Display Results ---
fprintf('Generated a stable %d-order system with %d zeros.\n\n', DEN_ORDER, NUM_ORDER);
disp('Transfer Function G(s):');
disp(G);

% --- Verify Stability with a Pole-Zero Plot ---
% All poles (marked with 'x') should be in the left-half of the plane.
figure;
pzmap(G);
title(sprintf('Pole-Zero Map (Order %d/%d)', NUM_ORDER, DEN_ORDER));
grid on;
sgrid; % Shows damping and natural frequency lines

%}
