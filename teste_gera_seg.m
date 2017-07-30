clear all
close all
clc

% The main idea fo this code is calculate eq. (4) for K colors.

%K for K-means++
K = 85;
% m nearest neighbors (eq. (4), m = M)
M = round(K/2);

% With segmentation from http://cs.brown.edu/~pff/segment/
% command_srt = ['./segment/segment ' num2str(sigma_seg) num2str(K_seg) num2str(min_seg) 'Smap.ppm Smap_seg.ppm'];
sigma_seg = 0.6; % choose sigma <5
K_seg = 5000; % choose <1000
min_seg = 200; %choose I dont know!

%show parcial pics?
pimage = 1;
listing = dir('Imgs/*.jpg');
N = size(listing,1);
filename = ['Imgs/' listing(2).name];
% original image
img = imread(filename);
% ground truth, to compare
filename = [filename(1:end-3) 'png'];
gtr = imread(filename);
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
if pimage
    figure;
    imagesc(SI)
    colormap(gray)
    title('Without Color space Smoothing')
end

% To achive the m nearest color, first we most know the all distances.
Dlab = zeros(K,K);
idm = zeros(M,K);
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


Smap = reshape(Sc2(idx),ll,cc); 
if pimage
    figure;
    imagesc(Smap)
    colormap(gray)
    title('With Color space Smoothing')
end

Smap = uint8(normalizar(Smap).*255);
addpath ./segment
%usage: ./segment sigma k min input(ppm) output(ppm)
%  sigma = 0.5, K = 500, min = 50.
imwrite(Smap,'Smap.ppm')
command_srt = ['./segment/segment ' num2str(sigma_seg) ' ' num2str(K_seg) ' ' num2str(min_seg) ' ' 'Smap.ppm Smap_seg.ppm'];
system(command_srt);

Smap_seg = rgb2gray(imread('Smap_seg.ppm'));
% 
if pimage
    figure;
    imagesc(Smap_seg)
    colormap(gray)
    title('After graph segmentation.')
end

save teste_seg_gerado.mat


%%%