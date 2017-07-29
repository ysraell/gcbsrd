clear all
close all
clc

tic;
load teste_seg_gerado_ok

%sigma_s is the spatial distance weighting (eq. (7))
%sigma_s = [0.1:0.1:20];
sigma_s = [0.1:0.1:20];
% At final, cut the background regions (or try it)
%threshold_final_cut = [0.1:0.005:0.98]';

threshold_final_cut = [0.05:0.01:0.99];

%show parcial pics?
pimage = 0;

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

gtr_bin = gtr/255;
T_T = max(size(threshold_final_cut));
T_S = max(size(sigma_s));
err_p = zeros(T_T,T_S);
err_n = zeros(T_T,T_S);
err_s = zeros(T_T,T_S);
normones = norm(ones(ll,cc));

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

    if pimage
        figure;
        imagesc(Result)
        colormap(gray)
        title('Eq. (7) S_rk')
    end
    
    Trh = max(Result(:));
    
    %%
    for t=1:T_T
        R_trh = uint8(Result>threshold_final_cut(t)*Trh);
        err = abs(double(R_trh)-double(gtr_bin));
        err_p(t,s) = 1-sum(err(:))/(ll*cc);
        err_n(t,s) = norm(err)/normones;
        err_s(t,s) = (1-ssim(double(R_trh),double(gtr_bin)))/2;
    end
end



%% norm
[tempv,temp] = min(err_n(:));
[t,s]=ind2sub(size(err_n),temp);

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

    if pimage
        figure;
        imagesc(Result)
        colormap(gray)
        title('Eq. (7) S_rk')
    end
    
    Trh = max(Result(:));
    R_trh_n = uint8(Result>threshold_final_cut(t)*Trh);
    Err_n = abs(double(R_trh_n)-double(gtr_bin));

%% fig norm
figure;
subplot(2,2,1)
imshow(R_trh_n.*255)
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
imshow(Err_n.*255)
colormap(gray)
temp = ['Relative error $\|$Our - GndTr$|/\|1_{l \times c}\| =$ ' num2str(tempv)];
title(temp,'Interpreter','LaTex')

%% ssim
[tempv,temp] = min(err_s(:));
[t,s]=ind2sub(size(err_s),temp);

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

    if pimage
        figure;
        imagesc(Result)
        colormap(gray)
        title('Eq. (7) S_rk')
    end
    
    Trh = max(Result(:));
    R_trh_s = uint8(Result>threshold_final_cut(t)*Trh);
    Err_s = abs(double(R_trh_s)-double(gtr_bin));

    
    %% fig ssim
figure;
subplot(2,2,1)
imshow(R_trh_s.*255)
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
imshow(Err_s.*255)
colormap(gray)
temp = ['Error (Our - GndTr). SSIM = ' num2str(tempv)];
title(temp)

%% Presision
[tempv,temp] = max(err_p(:));
[t,s]=ind2sub(size(err_p),temp);

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

    if pimage
        figure;
        imagesc(Result)
        colormap(gray)
        title('Eq. (7) S_rk')
    end
    
    Trh = max(Result(:));
    R_trh_p = uint8(Result>threshold_final_cut(t)*Trh);
    Err_p = abs(double(R_trh_p)-double(gtr_bin));

    
    %% fig ssim
figure;
subplot(2,2,1)
imshow(R_trh_p.*255)
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
imshow(Err_p.*255)
colormap(gray)
temp = ['Error (Our - GndTr). Precision = ' num2str(tempv)];
title(temp)

%%
recall = 1-threshold_final_cut;

figure;
mesh(sigma_s,recall,err_n)

figure;
mesh(sigma_s,recall,err_s)

figure;
mesh(sigma_s,recall,err_p)

time =toc;


save Limiar_sigma_var.mat