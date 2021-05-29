% This example shows how to use the method proposed in:
% Fernando Pérez-Bueno, Miguel Vega, Valery Naranjo, Rafael Molina, Aggelos K. Katsaggelos,
%
%
% 
%% Load image and reference vectors
clc,clear all
I = imread('histWB.jpg');
load 'MLandini' RM;
[m,n,nc] = size(I);
subplot(241),imshow(I)
title('Original H&E Image')
%% Deconvolution
ns=2; %number of stains
p=1; %prior parameter
prior= 'lp' %'lp' or 'log'
filtersetname='fo'; %high pass filters ('none', 'fohv', 'fo') 

if strcmp(prior,'lp')
    [CT, M, alpha, beta, gamma] = BCDHElpnf(im2double(I), p, filtersetname, RM(:,1:ns));
elseif strcmp(prior,'log')
    [CT, M, alpha, beta, gamma] = BCDHElognf(im2double(I), filtersetname, RM(:,1:ns));
else
    disp('prior value not valid')
end
disp('completed')

%% Band visualization (OD space)

ns = size(M,2)
concentrations = reshape(CT',m,n,ns);

%figure()
subplot(242),imshow(concentrations(:,:,1))
title('OD H Band')
subplot(246),imshow(concentrations(:,:,2))
title('OD E Band')


%% Band reconstruction (RGB space)
Hrec_OD = reshape((M(:,1)*CT(1,:))',m,n,nc);
Hrec_RGB = OD2intensities(Hrec_OD);

Erec_OD = reshape((M(:,2)*CT(2,:))',m,n,nc);
Erec_RGB = OD2intensities(Erec_OD);

%figure()
subplot(243),imshow(Hrec_RGB)
title('RGB H Band')
subplot(247),imshow(Erec_RGB)
title('RGB E Band')

%% Image Normalization

Iref = imread('Reference.jpg');
%Deconvolution of the reference image
if strcmp(prior,'lp')
    [Cref, Mref, ~,~,~] = BCDHElpnf(im2double(Iref), p, filtersetname, RM(:,1:ns));
elseif strcmp(prior,'log')
    [Cref, Mref, ~,~,~,] = BCDHElognf(im2double(Iref), filtersetname, RM(:,1:ns));
else
    disp('prior value not valid')
end
disp('completed')

%% Range adjustment
CT_Rmax = prctile(CT',99)
Cref_Rmax= prctile(Cref',99)
norm_fac=Cref_Rmax./CT_Rmax
CT_norm=CT.*norm_fac';

%Reconstruction
Yrec_norm=Mref(:,1:ns)*CT_norm;
Y2d_norm=reshape(Yrec_norm',m,n,nc);
Irec_norm=OD2intensities(Y2d_norm);

subplot(244),imshow(Irec_norm)
title('Normalized image')
subplot(248),imshow(Iref)
title('Reference')
