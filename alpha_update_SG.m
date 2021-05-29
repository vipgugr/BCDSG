function [alpha, WTfilterednu] = alpha_update_SG(CTfilterednu,DnutDnu,kappa,SigmaC,epsW)

    ns = size(CTfilterednu,1);
    tamm = size(CTfilterednu,2);
    [m, n] = size(DnutDnu);
    
    WTfilterednu = zeros(ns,tamm);
    
    kappa_f=kappa{1};
    alpha_f=kappa{3};
    
    alpha = zeros(ns,1);
    for s=1:ns
        tmp = reshape(SigmaC(:,s),m,n);
        tmp = tmp .* DnutDnu;
        traza = sum(tmp(:));
        
        u(s,:) = epsW + ( abs ( CTfilterednu(s,:).*CTfilterednu(s,:) + traza/tamm) ).^0.5;
        
        W = kappa_f (u(s,:));
        val=mean(W(:)) + eps;
        alpha(s)=alpha_f(val);
         
        WTfilterednu(s,:) = W(:)';
    end

end
