function [SNR]=ComputePredictedSNRPoint(c,delta,gamma)

b=(c+1-delta)/gamma;
b1=b^2;
b2=-4*c/gamma^2;
SNR=10*log10(1/2*(b + sqrt(b1+b2)));