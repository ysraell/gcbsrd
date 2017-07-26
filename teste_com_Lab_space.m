clear all
close all
clc

% The main idea fo this code is calculate eq. (4) for K colors.

%K for K-means++
K = 85;
% m nearest neighbors (eq. (4), m = M)
M = 20;

img = imread('Imgs/818.jpg');

% ground truth, to compare
%gtr = imread('Imgs/818.png');

% distance metric between pixels (or colors)
% simple L-2:
% D= @(a,b) norm(a-b);
% take the luminance in account:
% RGB -> Y   Y=0.2126R+0.7152G+0.0722B
% Y= @(a) [0.2126 0.7152 0.0722]*[a(1) a(2) a(3)]';
% D= @(a,b) norm(a-b)+norm(Y(a)-Y(b));

% Eq. (3):

% first, do a qunatization!
[ll,cc,~] = size(img);
m_rgb = [reshape(img(:,:,1),ll*cc,1)...
         reshape(img(:,:,2),ll*cc,1)...
         reshape(img(:,:,3),ll*cc,1)];
[idx,C] = kmeanspp(double(m_rgb)',K);
m_qua = C(:,idx)';
hist_qua = hist(idx',K)';
hist_relative = hist_qua./sum(hist_qua);

[L_cie,a_cie,b_cie] = RGB2Lab(C');

D= @(l,j) norm([L_cie(l) a_cie(l) b_cie(l)]'-[L_cie(j) a_cie(j) b_cie(j)]');

%S(c)
Sc = zeros(K,1);
for l=1:K
    for j=1:K
        Sc(l) = Sc(l)+hist_relative(j)*D(l,j);
    end
end

SI = reshape(Sc(idx),ll,cc); 
figure;
title('Without Color space Smoothing')
imagesc(SI)
colormap(gray)

% To achive the m nearest color, first we most know the all distances.
Dlab = nan(K,K);
idm = nan(M,K);
T = zeros(K,1);
Lab = [L_cie a_cie b_cie];
%the loops bellow is not the smartest way, but I'm feel so lazy.
%ok I'll try something more intelligent
for i=1:K
    for j=i:K
        Dlab(i,j) = norm(Lab(i,:)-Lab(j,:));
        Dlab(j,i) = Dlab(i,j);
    end
    
    [~,tmp] = sort(Dlab(:,i),'ascend');
    for m=1:M
        T(i)=T(i)+D(i,tmp(m+1));
        idm(m,i)=tmp(m+1);
    end
    
end

Sc2 = zeros(K,1);
for i=1:K   
    ST = 0;
    for m=1:M
        ST=ST+(T(i)-D(i,idm(m,i)))*Sc(idm(m,i));
    end
    Sc2(i) = ST/(T(i)*(m-1));
end


SI = reshape(Sc2(idx),ll,cc); 
figure;
title('With Color space Smoothing')
imagesc(SI)
colormap(gray)






