function [detection_table, coeff_table] = AnalyzeDetectionTable(input_table)
    % c(1) = rate of change, c(2) = x-offset, c(3) = multiplier, c(4) = offset
    sigfun = @(c,x) (c(3) .* (1./(1 + exp(-c(1).*(x-c(2)))))) + c(4); 

    x = input_table.TestStimAmp;
%     x = input_table.CondStimAmp;
    y = strcmpi(input_table.Response, 'correct');
    [ux, ~, ic] = unique(x);
    detection_vector = zeros(size(ux));
    for d = 1:length(detection_vector)
        if ux(d) == 0
            detection_vector(d) = mean(~y(ic == d));
        else
            detection_vector(d) = mean(y(ic == d));
        end
    end
    
    % Make d'
    if detection_vector(1) < 1e-3
        z_fa = norminv(1e-3);
    else
        z_fa = norminv(detection_vector(1));
    end
    z_hit = norminv(detection_vector);
    z_hit(isinf(z_hit)) = norminv(1-1e-3);
    dprime = z_hit - z_fa;

    detection_table = table(ux, detection_vector, dprime, 'VariableNames', {'StimAmp',  'pDetect', 'dPrime'});
    
    dt_coeffs = lsqcurvefit(sigfun,... % Function to fit
                            [1.7/mean(detection_table.StimAmp),...% c(1)
                            mean(detection_table.StimAmp),...% c(2)
                            max(detection_table.pDetect),...% c(3)
                            min(detection_table.pDetect)],... % c(4)
                            detection_table.StimAmp, detection_table.pDetect,... % x & y
                            [0, 0, 0, 0],... % Minimum coefficient allowed
                            [1, 100, 1, 1],... % Maximum coefficient allowed
                            optimset('Display','off'));

    dp_coeffs = lsqcurvefit(sigfun,... % Function to fit
                            [1.7/mean(detection_table.dPrime),...% c(1)
                            mean(detection_table.dPrime),...% c(2)
                            max(detection_table.pDetect),...% c(3)
                            min(detection_table.pDetect)],... % c(4)
                            detection_table.StimAmp, detection_table.dPrime,... % x & y
                            [0, 0, 0, 0],... % Minimum coefficient allowed
                            [1, 100, 4, 1],... % Maximum coefficient allowed
                            optimset('Display','off'));
    coeff_table = table(dt_coeffs, dp_coeffs, 'VariableNames', {'pDetect', 'dPrime'});
end