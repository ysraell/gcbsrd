function [prec,recal]=func_SaliencyCut(filename)


% The main idea fo this code is calculate eq. (4) for K colors.

%K for K-means++
K = 85;
% m nearest neighbors (eq. (4), m = M)
M = round(K/4);

% With segmentation from http://cs.brown.edu/~pff/segment/
% command_srt = ['./segment/segment ' num2str(sigma_seg) num2str(K_seg) num2str(min_seg) 'Smap.ppm Smap_seg.ppm'];
% sigma_seg = 0.6; % choose sigma <5
% K_seg = 5000; % choose <1000
% min_seg = 200; %choose I dont know!

sigma_seg = 0.5; % choose sigma <5
K_seg = 2000; % choose <1000
min_seg = 100; %choose I dont know!

%show parcial pics?
% pimage = 0;

% original image
img = imread(filename);
% ground truth, to compare
filename = [filename(1:end-3) 'png'];
gtr = imread(filename);

% distance metric between pixels (or colors)
% simple L-2:
% D= @(a,b) norm(a-b);
% take the luminance in account:
% RGB -> Y   Y=0.2126R+img = imread0.7152G+0.0722B
% Y= @(a) [0.2126 0.7152 0.0722]*[a(1) a(2) a(3)]';
% D= @(a,b) norm(a-b)+norm(Y(a)-Y(b));

% Eq. (3):

% first, do a qunatization!
[ll,cc,~] = size(img);
m_rgb = [reshape(img(:,:,1),ll*cc,1)...
         reshape(img(:,:,2),ll*cc,1)...
         reshape(img(:,:,3),ll*cc,1)];
[idx,C] = kmeanspp(double(m_rgb)',K);
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

% SI = reshape(Sc(idx),ll,cc); 
% if pimage
%     figure;
%     imagesc(SI)
%     colormap(gray)
%     title('Without Color space Smoothing')
% end

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
% if pimage
%     figure;
%     imagesc(Smap)
%     colormap(gray)
%     title('With Color space Smoothing')
% end

Smap = uint8(normalizar(Smap).*255);
addpath ./segment
%usage: ./segment sigma k min input(ppm) output(ppm)
%  sigma = 0.5, K = 500, min = 50.
file_temp = num2str(randi(999999),'%.6d');
file_ppm = [file_temp  '.ppm'];
imwrite(Smap,file_ppm)
command_srt = ['./segment/segment ' num2str(sigma_seg) ' ' num2str(K_seg) ' ' num2str(min_seg) ' ' file_ppm ' ' file_temp '_seg.ppm'];
[~,~] = system(command_srt);

Smap_seg = rgb2gray(imread([file_temp '_seg.ppm']));
command_srt = ['rm -f ' file_ppm ' ' file_temp '_seg.ppm'];
[~,~] = system(command_srt);
% 
% if pimage
%     figure;
%     imagesc(Smap_seg)
%     colormap(gray)
%     title('After graph segmentation.')
% end

%sigma_s is the spatial distance weighting (eq. (7))
%sigma_s = [0.1:0.1:20];
sigma_s = [0.1 0.5 1:2:20];
% At final, cut the background regions (or try it)
%threshold_final_cut = [0.1:0.005:0.98]';

threshold_final_cut = 0.05:0.01:0.99;

values_r = unique(Smap_seg);
T_r = size(values_r,1);

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

% if pimage
%     figure;
%     imagesc(Result_seg)
%     colormap(gray)
%     title('Eq. (5) S_rk')
% end

gtr_bin = gtr/255;
T_T = max(size(threshold_final_cut));
T_S = max(size(sigma_s));
err_p = zeros(T_T,T_S);

for s=1:T_S
    %eq. (7)

    %% Bloco
    Srk = zeros(T_r,1);
    %wr = wr./sum(wr);
    for rk=1:T_r
        for ri=1:T_r
            if ri~=rk
                Srk(rk) = Srk(rk)+exp(-Ds(rk,ri)/(sigma_s(s)^2))*wr(ri)*Dr(rk,ri);
    %             [Sr(rk) Ds(rk,ri) exp(-Ds(rk,ri)/1000) wr(ri) Dr(rk,ri)]
            end
        end
        Srk(rk) = ws(rk)*Srk(rk);
    end

    Result = zeros(ll,cc);
    for r=1:T_r 
        Result=Result+ double(Segs(:,:,r)).*floor(Srk(r));
    end

%     if pimage
%         figure;
%         imagesc(Result)
%         colormap(gray)
%         title('Eq. (7) S_rk')
%     end
    
    Trh = max(Result(:));
    
    %%
    for t=1:T_T
        R_trh = uint8(Result>threshold_final_cut(t)*Trh);
        err = abs(double(R_trh)-double(gtr_bin));
        err_p(t,s) = 1-sum(err(:))/(ll*cc);
    end
end



%% Presision
[prec,temp] = max(err_p(:));
[t,~]=ind2sub(size(err_p),temp);
recal =1-threshold_final_cut(t);


end

%EOF