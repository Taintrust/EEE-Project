clc
%results_da = runpf('case15da'); % low-voltage distribution system studies.
%Base MVA = 1


results_nbr = runpf('case15nbr'); % Run power flow with updated options
    
%Base MVA = 100
%These systems are similar except for the values of R and X


