clear all
close all
clc

load Limiar_sigma_var_ok.mat

figure;
plot(recall,err_p')

break
figure;
mesh(sigma_s,recall,err_n)

figure;
mesh(sigma_s,recall,err_s)

figure;
mesh(sigma_s,recall,err_p)

