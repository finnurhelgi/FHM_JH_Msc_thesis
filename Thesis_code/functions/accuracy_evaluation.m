%Used in chapter 4.2.4

 % Input: 
    % sys_fom_eigs - white box eigenvalues
    % sys_fit - model eigenvalues
    % range_opt - option to filter eigenvalues based on range input
    % real_range - real  component of target area
    % imag_range - imag. component of targe area
    % damp_tol - real tolerance of "hit box"
    % freq_tol - imag tolerance of "hit box"
    % min_tol - minimum tolerance of "hit box"

    % Output:
    % total_values - number of white box values in target area
    % id_values - number of identified values
    % box_width, box_height, x_box, y_box - "hit box" dimensions and
    % coordinates
    % rel_error - relative error output
    % model_value_nr - total number of model values in area of interest
    
function [rmserr, total_values, id_values, box_width, box_height, x_box, y_box,rel_error,model_value_nr] = accuracy_evaluation(sys_fom_eigs, sys_fit, range_opt, real_range, imag_range, damp_tol, freq_tol,min_tol)

% Extract eigenvalues
p_n = sys_fom_eigs;
c_n = eig(sys_fit);


% Filter eigenvalues based on the range if required
if range_opt == 1
    real_min = real_range(1);
    real_max = real_range(2);
    imag_min = imag_range(1);
    imag_max = imag_range(2);
    is_in_range = @(z) real(z) >= real_min & real(z) <= real_max & imag(z) >= imag_min & imag(z) <= imag_max;
    p_n = p_n(arrayfun(is_in_range, p_n));
    c_n = c_n(arrayfun(is_in_range, c_n));
end

model_value_nr = length(c_n);
d_n = zeros(size(p_n));
dist = zeros(size(p_n));
p_n2 = p_n;

x_box = zeros(5, length(p_n));
y_box = zeros(5, length(p_n));

for i = 1:length(p_n)
    % Calculate the distances in the real and imaginary parts between p_n(i) and all elements in c_n
    real_distances = abs(real(p_n(i)) - real(c_n));
    imag_distances = abs(imag(p_n(i)) - imag(c_n));

    tol_damp = abs((2*damp_tol/100)*real(p_n(i)));
    tol_freq = abs((2*freq_tol/100)*imag(p_n(i)));
   
    box_width = min_tol + tol_damp;
    box_height = min_tol*2*pi + tol_freq;

    % Check if any eigenvalue is within the predefined box
    in_box = (real_distances <= box_width/2) & (imag_distances <= box_height/2);
    
%     Uncomment of you want to see the hit box of every eigenvalue
%     x_center = real(p_n(i));
%     y_center = imag(p_n(i));
%     x_box(:, i) = [x_center - box_width/2, x_center + box_width/2, x_center + box_width/2, x_center - box_width/2, x_center - box_width/2];
%     y_box(:, i) = [y_center - box_height/2, y_center - box_height/2, y_center + box_height/2, y_center + box_height/2, y_center - box_height/2];   
    if any(in_box)
        % Find the index of the closest eigenvalue within the box
        euc_distances = sqrt(real_distances.^2 + imag_distances.^2);
        [min_dist, minIndex] = min(euc_distances);
        % Assign the corresponding value from c_n to d_n
        d_n(i) = c_n(minIndex);
        dist(i) = min_dist;
        
        % Define the coordinates of the box corners
        x_center = real(p_n(i));
        y_center = imag(p_n(i));
        x_box(:, i) = [x_center - box_width/2, x_center + box_width/2, x_center + box_width/2, x_center - box_width/2, x_center - box_width/2];
        y_box(:, i) = [y_center - box_height/2, y_center - box_height/2, y_center + box_height/2, y_center + box_height/2, y_center - box_height/2];
    end
end

index = find(d_n == 0);
p_n2(index) = [];
d_n = d_n(d_n ~= 0);
dist = dist(dist ~= 0);

total_values = length(p_n);
id_values = length(p_n2);

rel_error =sum(dist)/id_values;

end
