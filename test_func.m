% ------ Basic Testing Code ------ %
% Testing SISO implementation
syms x
err = 0;
while err ~= 1
    % Create random degree TF, max degree of 10
    deg = randi(10);
    [~, cn, cd] = generate_stable_tf(deg-1, deg);
    % Calculate H2 norm
    [pn_plus_1, zn_minus_1, an] = solve_fraction_free_H2_norm(cd, cn);
    % Compare to MATLAB implementation, error if not matching
    if ~(zn_minus_1/(2*an*pn_plus_1) - norm(tf(cn,cd),2)^2 <= 1E-1)
        disp("ERROR")
        err = 1;
    end
end

%% ------ Testing MIMO implementation ------ %%
% Example taken from paper
A = {[1 3 3 1], [1 4 6 4 1]; [1 5 10 10 5 1], [1 6 15 20 15 6 1]};
C = {[1 1 1], [1 1 1 1]; [1 1 1 1 1], [1 1 1 1 1 1]};
% z0 matches paper
[pn_plus_1, zn_minus_1, an_out] = solve_fraction_free_H2_norm_mimo(A,C);
% Output correct
zn_minus_1/(2*an_out*pn_plus_1)