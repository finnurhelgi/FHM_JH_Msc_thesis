%Used by the PC-CF method to identify pole clusters
 % Input: 
    % complexNumbers - array of complex numbers
    % tolerancePercent - tolerance percentage for clustering
    % minimum tolerance - minimum tolerance value for clustering
    
    % Output:
    % result - array of clustered complex numbers
function [result,x_box,y_box] = cluster_poles(poles, tolerancePercent, mintol)

    
    x_box = zeros(5, length(poles));
    y_box = zeros(5, length(poles));


%     x_box = [];
%     y_box = [];
    % Initialize result array
    result = [];
    
    % Create a boolean array to mark complex numbers that have been clustered
    clustered = false(size(poles));
    
    % Special handling for numbers close to the zero point within 10^-5 tolerance
%     zeroTolerance = 1e-5;
%     closeToZero = abs(poles) <= zeroTolerance;
%     if any(closeToZero)
%         result = [result; poles(closeToZero)];
%         % Mark all instances close to 0+0i as clustered
%         clustered(closeToZero) = true;
%     end
    
    % Iterate through each complex number
    for i = 1:length(poles)
        if clustered(i)
            continue; % Skip if already clustered
        end
        
        % Get the current complex number
        z = poles(i);
        
        % Define the tolerance region
%         realTolerance = max(abs(real(z)) * tolerancePercent / 100, mintol); box_width = realTolerance*2;
%         imagTolerance = max(abs(imag(z)) * tolerancePercent / 100, mintol); box_height = imagTolerance*2;

        realTolerance = mintol + (abs(real(z)) * tolerancePercent / 100); box_width = realTolerance*2;
        imagTolerance = mintol*2*pi + (abs(imag(z)) * tolerancePercent / 100); box_height = imagTolerance*2;
        
        % Find numbers within the tolerance rectangular area
        inCluster = abs(real(poles) - real(z)) <= realTolerance & ...
                    abs(imag(poles) - imag(z)) <= imagTolerance;

        
        % Check if there are at least two other complex numbers in the cluster
        if sum(inCluster) >= 2
            % Mark the clustered numbers
            clustered(inCluster) = true;

            % Extract the clustered numbers
            clusterNumbers = poles(inCluster);
            
            % Find the median of the cluster
            medianReal = median(real(clusterNumbers));
            medianImag = median(imag(clusterNumbers));
            medianComplex = medianReal + 1i * medianImag; 
            
            % Find the number closest to the median
            [~, closestIndex] = min(abs(clusterNumbers - medianComplex));
            representative = clusterNumbers(closestIndex);
            
            x_center = real(representative);
            y_center = imag(representative);
            x_box(:, i) = [x_center - box_width/2, x_center + box_width/2, x_center + box_width/2, x_center - box_width/2, x_center - box_width/2];
            y_box(:, i) = [y_center - box_height/2, y_center - box_height/2, y_center + box_height/2, y_center + box_height/2, y_center - box_height/2]; 
         
            % Add the representative to the result
            result = [result; representative];
        end
    end

end
