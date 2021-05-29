function [M, SigmaM] = color_vector_updatefiltered(nfilters,YTfiltered,CTfiltered,SigmaC,M,RM,beta,gamma)


nc = size(YTfiltered{1},1);
ns = size(CTfiltered{1},1);
tamm = size(YTfiltered{1},2);

SigmaM = zeros(ns,1);


for s=1:ns
    expC2 = 0.0;
    expCeminus = zeros(nc,1);
    for nu=1:nfilters
        [~, eminus{nu}] = computingEsZs(YTfiltered{nu},CTfiltered{nu},M);
        expC2 = expC2 + beta{nu} * (CTfiltered{nu}(s,:) * CTfiltered{nu}(s,:)' + sum(SigmaC{nu}(:,s)));
        expCeminus = expCeminus + beta{nu} * (CTfiltered{nu}(s,:)*eminus{nu}(:,:,s))';
    end
        
    SigmaM(s) = 1.0 ./( expC2 + gamma(s));
    M(:,s) = SigmaM(s) .* (expCeminus  + gamma(s)*RM(:,s)); %%% CUIDADO AQUI. ESTO HAY QUE COMPROBARLO CON LA TEORIA
    % M(M<0)=0;
    scalefactor = norm(M(:,s));
    M(:,s) = M(:,s) / scalefactor;
    SigmaM(s) = SigmaM(s) / (scalefactor.^2);  %%%% IMPORTANTE: Si se normaliza la media hay que normalizar tambi�n la varianza
end

end
