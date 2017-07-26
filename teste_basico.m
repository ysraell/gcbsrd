clear all
close all
clc

%K for K-means++
K = 85;

img = imread('Imgs/818.jpg');
gtr = imread('Imgs/818.png');

% distance metric between pixels (or colors)
% simple L-2:
% D= @(a,b) norm(a-b);
% take the luminance in account:
% RGB -> Y   Y=0.2126R+0.7152G+0.0722B
Y= @(a) [0.2126 0.7152 0.0722]*[a(1) a(2) a(3)]';
D= @(a,b) norm(a-b)+norm(Y(a)-Y(b));

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

%S(c)
Sc = zeros(K,1);
for l=1:K
    for j=1:K
        Sc(l) = Sc(l)+hist_relative(j)*D(C(:,l)',C(:,j)');
    end
end

SI = reshape(Sc(idx),ll,cc);


figure;
imagesc(SI)
colormap(gray)






