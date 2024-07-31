% inputs: state space model of system
% outputs: 
%   par_facs - mode index, frequency and damping of modes in
%   percentage
%   P - Participation factor matrix containing modes of interest (under 200
%   Hz)

function [par_facs,P] = small_signal_assessment(sys)

A = sys.A;
k = min(size(A));

[V, D] = eig(A);
[W, ~] = eig(A.');

% Normalize the right eigenvectors (V)
for i = 1:size(V, 2)
    V(:, i) = V(:, i) / norm(V(:, i));
end

% Normalize the left eigenvectors (W)
for i = 1:size(W, 1)
    W(i, :) = W(i, :) / norm(W(i, :));
end

P = round(abs(V*W'),1);

eigenv = eig(A);
f_a = zeros(k,1);
d_a = zeros(k,1);

for i = 1:length(eigenv)
    if imag(eigenv(i)) == 0
        continue
    else
        f_a(i) = imag(eigenv(i))/(2*pi);
        d_a(i) = -100*real(eigenv(i))/(sqrt(real(eigenv(i))^2 + imag(eigenv(i))^2));
    end
end

ss_freq = 200;
% ss_ind = find(f_a < ss_freq  & f_a ~= 0 & f_a > -ss_freq);
ss_freq_ind = find(f_a < ss_freq  & f_a >= 0);

pf_ind = find(real(eigenv)>-100 & imag(eigenv)/(2*pi) < ss_freq & imag(eigenv)/(2*pi) > 0);

ss_eigs = eigenv(pf_ind);
ss_fa = f_a(pf_ind);
ss_da = d_a(pf_ind);

par_facs = [pf_ind,ss_fa,ss_da];

P = P(:,pf_ind);

end