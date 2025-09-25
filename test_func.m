% ------ Basic Testing Code ------ %
% Testing SISO implementation
err = 0;
while err ~= 1
    % Create random degree TF, max degree of 10
    deg = randi(10);
    [~, cn, cd] = generate_stable_tf(deg-1, deg);
    % Calculate H2 norm
    [pn_plus_1, zn_minus_1, an] = solve_fraction_free_H2_norm(cd, cn);
    matlabH2n = norm(tf(cn,cd),2)^2;
    paperH2n = zn_minus_1/(2*an*pn_plus_1);
    fprintf("Difference: %.6f, H2ns: %.6f\n", abs(paperH2n - matlabH2n), paperH2n)
    % Compare to MATLAB implementation, error if not matching
    if ~(paperH2n - matlabH2n <= 0.001*matlabH2n)
        disp("ERROR")
        err = 1;
    end
end

%% ------ Testing MIMO implementation ------ %%
% Example taken from paper
% Define the numerator and denominator coefficients for the 4 transfer functions
c1 = [1 1 1]; a1 = [1 3 3 1];
c2 = [1 1 1 1]; a2 = [1 4 6 4 1];
c3 = [1 1 1 1 1]; a3 = [1 5 10 10 5 1];
c4 = [1 1 1 1 1 1]; a4 = [1 6 15 20 15 6 1];
A = {a1, a2; a3, a4};
C = {c1, c2; c3, c4};
% z0 matches paper
[pn_plus_1, zn_minus_1, an_out] = solve_fraction_free_H2_norm_mimo(A,C);
newMIMOH2n = zn_minus_1/(2*an_out*pn_plus_1);

% Create the SISO transfer function objects
G11 = tf(c1, a1);
G12 = tf(c2, a2);
G21 = tf(c3, a3);
G22 = tf(c4, a4);

% Assemble the 2x2 MIMO transfer function object G(s)
G_MIMO = [G11, G12; 
          G21, G22];

% Compute the H2 Norm of the MIMO system
H2_norm_MIMO = norm(G_MIMO, 2);

% Display the result
fprintf('The MIMO H2 Norm of the 2x2 system is: %.6f\n', H2_norm_MIMO);