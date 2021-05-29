function kappa_f = getkappa_f(functionname,parameter)

    switch functionname
        
        case 'log'
            kappa_f= @weight_log;
            
        case 'bottonup'
            kappa_f = @(nu) weight_bottonup(nu,parameter);
            
        case 'lp'
            kappa_f = @(nu) weight_lp(nu,parameter);
            
        case 'l2'
            kappa_f = @weight_l2;
            
        case 'exp'
            kappa_f = @(nu) weight_exp(nu,parameter);
            
    end
end