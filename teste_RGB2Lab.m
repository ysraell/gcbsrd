clear all
close all
clc

img = imread('Imgs/818.jpg');
gtr = imread('Imgs/818.png');



[L,a,b] = RGB2Lab(img);
