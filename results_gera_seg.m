clear all
close all
clc

listing = dir('Imgs/*.jpg');
N = size(listing,1);

P = zeros(N,1);
R = zeros(N,1);
tic
parfor i=1:N
    filename = ['Imgs/' listing(i).name];
    [prec,recall] = func_SaliencyCut(filename);
    P(i)=prec;
    R(i)=recall;
    disp([i N prec recall])
end
time = toc;
disp([mean(P) std(P) min(P) max(P)])
disp([mean(R) std(R) min(R) max(R)])


save results_gera_seg_FullDB.mat

disp('poweroff')
pause(60)
system('poweroff')