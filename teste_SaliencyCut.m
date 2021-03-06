clear all
close all
clc

% To use GrabCut, clone from https://github.com/xiumingzhang/grabcut
addpath ~/Documents/CMP165/trab3/grabcut

% The main idea fo this code is calculate eq. (4) for K colors.

%K for K-means++
K = 85;
% m nearest neighbors (eq. (4), m = M)
M = round(K/2);
%sigma_s is the spatial distance weighting (eq. (7))
sigma_s = 10;

% With segmentation from http://cs.brown.edu/~pff/segment/
% command_srt = ['./segment/segment ' num2str(sigma_seg) num2str(K_seg) num2str(min_seg) 'Smap.ppm Smap_seg.ppm'];
sigma_seg = 0.5; % choose sigma <5
K_seg = 5000; % choose <1000
min_seg = 50; %choose I dont know!

% At final, cut the background regions (or try it)
threshold_final_cut = 0.8;

%show parcial pics?
pimage = 0;

% original image
img = imread('Imgs/818.jpg');
% ground truth, to compare
gtr = imread('Imgs/818.png');

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


values_r = unique(Smap_seg);
T_r = size(values_r,1);
hist_seg = hist(double(Smap_seg(:)),T_r);

% To implement Dr (eq. (6) it is important to get a bin frame by region
Segs = uint8(zeros(ll,cc,T_r));
wr = zeros(T_r,1);
for r=1:T_r
    Segs(:,:,r) = uint8(Smap_seg==values_r(r));
    wr(r) = sum(sum(Segs(:,:,r)));
    
    
%     pause
%     figure;
%     title('After graph segmentation.')
%     imagesc(Segs(:,:,r))
%     colormap(gray)
% %     
%     B = bwboundaries(255-Segs(:,:,r).*255);
%     
%     pause
%     figure;
%     title('After graph segmentation.')
%     hold on
%     imagesc(255-Segs(:,:,r).*255)
%     colormap(gray)
%     for k = 1:length(B)
%         boundary = B{k};
%         plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
%     end
%     hold off
    
end

% Calulating eq. (6)

Dr = zeros(T_r,T_r);
Ds =  Dr; %for eq. (7)
IDX= reshape(idx,ll,cc);
ws = zeros(T_r,1);
for rk=1:T_r
    for ri=rk:T_r
        if ri~=rk
        
            % for Dr
            rn =[rk ri];
            N = zeros(2,1);
            fcki = [];
            idxn = [];
            for r=1:2
                Temp = double(Segs(:,:,rn(r))).*IDX;
                temp = unique(Temp(:));
                idxn{r} = temp(temp>0);
                N(r) = size(idxn{r},1);
                temp = hist(Temp(:),N(r))';
                fcki{r} = temp./sum(temp);
            end
            
            for n1=1:N(1)
                for n2=1:N(2)
                    Dr(rk,ri)=Dr(rk,ri)+ fcki{1}(n1)*fcki{2}(n2)*D(idxn{1}(n1),idxn{2}(n2));
                end
            end
            
            Dr(ri,rk) = Dr(rk,ri);
            
            
            %for Ds
            Cr = [];
            for r=1:2
                sl =0;
                tl=0;
                sc = 0;
                tc =0;
                for l=1:ll
                    for c=1:cc
                        sl=sl+l.*double(Segs(l,c,rn(r)));
                        tl = tl+double(Segs(l,c,rn(r)));
                        sc=sc+c.*double(Segs(l,c,rn(r)));
                        tc = tc+double(Segs(l,c,rn(r)));
                    end
                end             
                Cr{r} = [sl/tl sc/tc];
            end
            Ds(rk,ri) = norm(Cr{1}-Cr{2});
            Ds(ri,rk) = Ds(rk,ri);
            
        else

            sl =0;
            tl=0;
            sc = 0;
            tc =0;
            for l=1:ll
                for c=1:cc
                    sl=sl+l.*double(Segs(l,c,rk));
                    tl = tl+double(Segs(l,c,rk));
                    sc=sc+c.*double(Segs(l,c,rk));
                    tc = tc+double(Segs(l,c,rk));
                end
            end             
            Cr = [sl/tl sc/tc];
            dk = norm((Cr-[l/2 c/2])./[l c]);
            ws(rk) = exp(-9*dk^2);            
        end
    end
end


% Calulating eq. (5)
Sr = zeros(T_r,1);
wr = wr./sum(wr);
for rk=1:T_r
    for ri=1:T_r
        if ri~=rk
            Sr(rk) = Sr(rk)+wr(ri)*Dr(rk,ri);
        end
    end
end

Result_seg = zeros(ll,cc);
for r=1:T_r
    Result_seg=Result_seg+ double(Segs(:,:,r)).*floor(Sr(r));
end

if pimage
    figure;
    imagesc(Result_seg)
    colormap(gray)
    title('Eq. (5) S_rk')
end

%eq. (7)

Srk = zeros(T_r,1);
%wr = wr./sum(wr);
for rk=1:T_r
    for ri=1:T_r
        if ri~=rk
            Srk(rk) = Srk(rk)+exp(-Ds(rk,ri)/(sigma_s^2))*wr(ri)*Dr(rk,ri);
%             [Sr(rk) Ds(rk,ri) exp(-Ds(rk,ri)/1000) wr(ri) Dr(rk,ri)]
        end
    end
    Srk(rk) = ws(rk)*Srk(rk);
end

Result = zeros(ll,cc);
for r=1:T_r 
    Result=Result+ double(Segs(:,:,r)).*floor(Srk(r));
end

if pimage
    figure;
    imagesc(Result)
    colormap(gray)
    title('Eq. (7) S_rk')
end

% Just do a threshold
Trh = max(Result(:));
R_trh = uint8(Result>threshold_final_cut*Trh);

% figure;
% title('Simple Threshold (T<max(Regions))')
% imagesc(R_trh)
% colormap(gray)

% I dont konw how use this:
% R_gc = grabcut(Result);
% figure;
% title('After GrabCut')
% imagesc(R_gc)
% colormap(gray)

gtr_bin = gtr/255;
err = abs(double(R_trh)-double(gtr_bin));

figure;
subplot(2,2,1)
imshow(R_trh.*255)
colormap(gray)
title('Our')
subplot(2,2,2)
imshow(gtr)
colormap(gray)
title('Ground Truth')
subplot(2,2,3)
imshow(img)
colormap(gray)
title('Original image')
subplot(2,2,4)
imshow(err.*255)
colormap(gray)
title('Error (Our - GndTr)')
