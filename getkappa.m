function kappa = getkappa(functionname,parameter)
    
    if(exist('parameter'))
        kappa_f = getkappa_f(functionname,parameter);
    else
        kappa_f = getkappa_f(functionname);
    end
    
    switch functionname
        
        case 'log'
            rho_f= @(nu) log(abs(nu) + eps);
            
        case 'lp'
            rho_f = @(nu) abs(nu).^parameter;
            
        case 'l2'
            rho_f = @(nu) abs(nu).^2;
            
        case 'exp'
            rho_f = @(nu) - parameter* exp(-(nu.*nu)/2/parameter);
            
    end
    
    switch functionname
        
        case 'log'
            alpha_f= @(val) 1.0 + 1.0/val;
            
        case 'lp'
            alpha_f = @(rho) 1.0/(eps + rho);
            
        case 'l2'
            alpha_f = @(rho) 1.0/(eps + rho);
            
        case 'exp'
            alpha_f = @(rho) - parameter* exp(-(nu.*nu)/2/parameter);
            
    end
    kappa={kappa_f rho_f alpha_f};
end