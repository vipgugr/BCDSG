function [CT,SigmaC]  = conc_update_SG(YT,CT,nfilters,M,SigmaM,filterTfilter,filters,filtersT,beta,alpha,m,n,WTfiltered)

    ns = size(CT,1);
    tamm = size(YT,2);
    
    eps_cs=1.0e-06;
    itmax_cs=1000;
    
  
    for s=1:ns 
        
        expM2 = M(:,s)' * M(:,s)+ 3.0 * SigmaM(s);
        
        auxSigmaC = 0;
        [zminus, ~] = computingEsZs(YT,CT,M);
        for nu=1:nfilters            
            Wfiltered{nu} = reshape(WTfiltered{nu}(s,:)',m,n);
            auxSigmaC = auxSigmaC + beta * expM2 + alpha{nu}(s) * mean(Wfiltered{nu}(:)) * filterTfilter{nu};          
            
            alphas{nu} = alpha{nu}(s);
        end
        SigmaC(:,s) = 1.0 ./auxSigmaC(:);
        
        inv_cov_fix = beta * expM2;
        indep_term = beta *zminus(:,s);
        
        fprintf('\tGradiente conjugado s = %d\n', s);
        
        inv_cov = get_mic_handle(inv_cov_fix,nfilters,filters,filtersT,alphas,Wfiltered,m,n);
        
        [cs, flag_cs,relres_cs,iter_cs ]=pcg(inv_cov, indep_term, eps_cs, itmax_cs,[], [], CT(s,:)' );
        fprintf('\tFlag: %d, RelRes: %d, Iters: %d\t\n',flag_cs,relres_cs,iter_cs);
        
        CT(s,:) = cs(:)';
        %CT(CT < eps) = eps; %%%%% Yo no forzaba no negatividad en cada vuelta. La forzaba al final, en todo caso
    end
end

function h = get_mic_handle(inv_cov_fix,nfilters,filters,filtersT,alphas,Wfiltered,nr,nc)
            h = @mic;
            function SigmaCscs = mic(cs)
                SigmaCscs = multiply_by_invcov(cs, inv_cov_fix,nfilters,filters,filtersT,alphas,Wfiltered,nr,nc);
            end
end

function SigmaCscs = multiply_by_invcov(cs,inv_cov_fix,nfilters,filters,filtersT,alphas,Wfiltered,nr,nc)

        F1 = inv_cov_fix* cs;
        
        %% Prior
        
        cs2df=fft2(reshape(cs,nr,nc));   
        priorTerm = zeros(nr,nc);
        
        for nu = 1:nfilters
%             priorTerm = priorTerm + alphas{nu} * imfilter( Wfiltered{nu}.* imfilter(cs2d,filters{nu}) , filtersT{nu} );
            temp = fft2( Wfiltered{nu} .* ifft2(filters{nu} .* cs2df) );
            priorTerm = priorTerm + alphas{nu} * filtersT{nu} .* temp;
        end
        
        priorTerm = reshape( ifft2(priorTerm) , nr*nc, 1);

        SigmaCscs = F1 + priorTerm;
    
    end

