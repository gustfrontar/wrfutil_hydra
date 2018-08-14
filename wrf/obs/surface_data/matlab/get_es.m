function [e_s]=get_es(T)
    % in:         T = temperature in Kelvin,
    % out:      e_s = saturation water vapour pressure over water in Pa 

  e_s = 100. .* 10.0 .^(23.832241 - 5.02808 .* log10(T) - 1.3816e-7 .* ...
                       10.0 .^(11.344 - 0.0303998 .* T) + 8.1328e-3 .* ...
                       10.0 .^(3.49149 - 1302.8844 ./ T) - 2949.076 ./ T);



