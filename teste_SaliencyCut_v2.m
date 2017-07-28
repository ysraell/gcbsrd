clear all
close all
clc

load teste_seg_gerado

%sigma_s is the spatial distance weighting (eq. (7))
sigma_s = 10;

% At final, cut the background regions (or try it)
threshold_final_cut = 0.6;

%show parcial pics?
pimage = 1;

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
