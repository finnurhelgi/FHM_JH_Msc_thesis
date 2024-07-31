%================================================
%=            P-factor analysis                 =
%================================================ 


% Read found A matrix 

A_scr_1_5 = load("A_scr_1_5_newest.mat").A_scr_1_5
A_scr_2_5 = load("A_scr_2_5_newest.mat").A_scr_2_5

k = min(size(A_scr_2_5));

[V, D] = eig(A_scr_2_5);
[W, ~] = eig(A_scr_2_5.');

% Normalize the right eigenvectors (V)
for i = 1:size(V, 2)
    V(:, i) = V(:, i) / norm(V(:, i));
end

% Normalize the left eigenvectors (W)
for i = 1:size(W, 1)
    W(i, :) = W(i, :) / norm(W(i, :));
end

P = (V*W');

eigenv = eig(A_scr_1_5);
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


% To color a minimum amoing of damping ratio 
zeta = 0.02% Minimum damping ratio

% Eigenvalues of the two systems 
eigenv_1_5 = eig(A_scr_1_5);
eigenv_2_5 = eig(A_scr_2_5);

omega_1_5 = imag(eigenv_1_5)
omega_2_5 = imag(eigenv_2_5)
Hz_1_5 = omega_1_5./(2*pi)
Hz_2_5 = imag(eigenv_2_5)./(2*pi)
damping_1_5 = real(eigenv_1_5)
damping_2_5 = real(eigenv_2_5)

damping_ratio_1_5 = -damping_1_5./(sqrt(damping_1_5.^2.+omega_1_5.^2))
damping_ratio_2_5 = -damping_2_5./(sqrt(damping_2_5.^2.+omega_2_5.^2))


sigma = linspace(-1000, 0, 1000); % range of real parts 
omega = sqrt(((sigma) / zeta).^2 - (sigma).^2)*2*pi; % corresponding imaginary parts
Hz_omega = omega/(2*pi)

k = min(size(A_scr_2_5));

[V, D] = eig(A_scr_2_5);
[W, ~] = eig(A_scr_2_5.');

% Normalize the right eigenvectors (V)
for i = 1:size(V, 2)
    V(:, i) = V(:, i) / norm(V(:, i));
end

% Normalize the left eigenvectors (W)
for i = 1:size(W, 1)
    W(i, :) = W(i, :) / norm(W(i, :));
end

P = V*W'
P = round(abs(V*W'),3);


P_norm = normalizeParticipationFactors(P)

figure;
mark = 10
plot(damping_2_5,Hz_2_5,'X',color='blue',MarkerSize=mark)

hold on
plot(damping_1_5,Hz_1_5,'d',color=[1,0.5,0],MarkerSize=mark)
% Adding FDM
% plot(real(eigenvalue_FDM),imag(eigenvalue_FDM./(2*pi)),'d',MarkerSize=mark)


plot(sigma,Hz_omega,color = [1, 0.5, 0])
plot(sigma,-Hz_omega,color = [1, 0.5, 0])

X = [sigma, fliplr(sigma)];
Y = [Hz_omega, fliplr(-Hz_omega)];

% Fill the area below the line
% Define the area to fill (below the line and within plot limits)
x_fill = [sigma, -1000000*ones(1, length(sigma))]; % Extend x-values to a lower limit (-10)
y_fill = [-Hz_omega, Hz_omega];     % Close the polygon by adding the y-values as zeros
fill(X,Y, 'cyan', 'FaceAlpha', 0.2, 'EdgeColor', 'none');% Fill the area with a color (e.g., light blue)
% y_fill = [omega, -omega];     % Close the polygon by adding the y-values as zeros
% fill(x_fill, y_fill, 'cyan', 'FaceAlpha', 0.3, 'EdgeColor', 'none');% Fill the area with a color (e.g., light blue)
hold off
xlim([-20, 1]);
ylim([-820,820]);
xlabel('Real sec^{-1}', 'Interpreter', 'tex',FontSize=14);
ylabel('Frequency [Hz]', 'Interpreter', 'tex',FontSize = 14);
grid on

zeta_text = sprintf('Zeta above %1.0f%%',zeta*100)
legend('VF SCR 2.5 System','VF SCR 1.5 System','','',zeta_text,fontsize=12)




function P_norm = normalizeParticipationFactors(P)

    % Calculate the magnitudes of the participation factors
    magnitudes = abs(P);

    % Sum the magnitudes for each column
    column_sums = sum(magnitudes, 1);

    % Normalize each participation factor by its column sum
    P_norm = P ./ column_sums;

end