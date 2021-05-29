function [CT, M, alpha, beta, gamma ] = BCDHElpnf(I, p, filtersetname, RM)
%BCDHE Bayesian TV Color Deconvolution of Hematoxilyn, Eosin slides
% Inputs: 
%   I: RGB observed image
%   p: value of the parameter p, usually 1
%   filtersetname: A choice of filters to obtain the differences
%          available options: 'none', 'fohv', 'fo'
%   RM: Reference color vector matrix of size 3xns
% Outputs:
%   CT: Stain concentration matrix of size nsxn_pixels
%   M: Estimated color vector matrix
%   alpha, beta, gamma: model parameters

    [m,n,nc] = size(I);
    tamm = m*n;
    ns = size(RM,2); %number of stains

    if (nc ~= 3 )
        error('Input image does not have 3 channels');
    end

    y2d = intensities2OD( I ); % dimension (m,n,nc)
    YT=reshape(y2d,m*n,nc)';

    
    % filtered images
    
    filters = getfilters(filtersetname);
    nfilters=prod(size(filters));    
    
    % Initial values
    CT = RM \ YT;
    M = RM;

    for s=1:ns
        SigmaM(s) = 0 ;
    end
    SigmaC = zeros(m*n,ns);    

    % Some stuff
    
    clear I y2d
    
    term=1.e-03;
    nitermin=5;
    nitermax=100;
    epsW=mean(CT(:))*1.0e-6;
    
    for nu=1:nfilters
        filtersT{nu} = flip(flip(filters{nu},1),2);
        filterTfilter{nu} = conv2(filtersT{nu},filters{nu});
        filterTfilter{nu} = psf2otf(filterTfilter{nu}, [m, n]) ; % otf
        filters{nu} = psf2otf(filters{nu}, [m, n]) ; % otf
        filtersT{nu} = psf2otf(filtersT{nu}, [m, n]) ; % otf
    end 
    
    kappa = getkappa('lp',p);
     

    % Primero trabajamos en el dominio filtrado para obtener la Matriz de
    % color
    CT0 = CT;
    iter = 1; convH = term +1.0; convE = term +1.0;
    %Iterations
    while ( (iter <= nitermin) || (((convH > term) || (convE > term)) && (iter <= nitermax)) )
        
        CT2d = reshape(CT,ns,m,n);
        for nu=1:nfilters
            for s=1:ns
%                 temp = imfilter(CT2d(s,:,:),filters{nu});
                temp = fft2 (reshape(CT2d(s,:,:),m,n) );
                temp(:,:) = ifft2( filters{nu} .* temp );
                CTfiltered{nu}(s,:) = reshape(temp,1,tamm);
            end
        end
 %       clear temp CT2d

        % Parameters update
        
        % beta
        beta = beta_update(YT,CT,SigmaC,M,SigmaM);
        fprintf('iter: %3d\t beta: %f\n',iter,beta)
        
        for nu=1:nfilters
 
            % alpha
            [alpha{nu}, WTfiltered{nu}] = alpha_update_SG(CTfiltered{nu},filterTfilter{nu},kappa,SigmaC,epsW);
            
            fprintf('nu %d\n',nu);
            fprintf('H\t alpha: %f\n',alpha{nu}(1))
            fprintf('E\t alpha: %f\n',alpha{nu}(2))
        
        end

        % gamma
        gamma = gamma_update(M,SigmaM,RM);

        fprintf('H\t gamma: %f\n',gamma(1))
        fprintf('E\t gamma: %f\n',gamma(2))

        % Color vector filtered update
        [M, SigmaM] = color_vector_update(YT,CT,SigmaC,M,RM,beta,gamma);

        % Concentration update
        
        [CT,SigmaC]  = conc_update_SG(YT,CT,nfilters,M,SigmaM,filterTfilter,filters,filtersT,beta,alpha,m,n,WTfiltered);
        
        convH = sum((CT(1,:)- CT0(1,:)).*(CT(1,:)- CT0(1,:))) / sum(CT0(1,:).*CT0(1,:));
        convE = sum((CT(2,:)- CT0(2,:)).*(CT(2,:)- CT0(2,:))) / sum(CT0(2,:).*CT0(2,:));
        CT0 = CT;
        
        fprintf('H\t conv: %e\n',convH)
        fprintf('E\t conv: %e\n',convE)
        M
        
        iter = iter +1;

    end       
    
    CT(CT < eps) = eps;

end


