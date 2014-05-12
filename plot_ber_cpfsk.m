function [AV_NF_BER AV_SMF_BER1 AV_SMF_BER2 AV_MF_BER]   = plot_ber_cpfsk()
%plot of BER for unfiltered and filtered case
% FSK signal corrupted by strong narrowband interference noise and little white noise
%AV_NF_BER -- average BER for no filtering case for given SINR
%AV_MF_BER -- average BER for plain median filtering case for given SINR
%AV_SMF_BER1 -- average BER for selective median filtering case for given SINR - pass 1
%AV_SMF_BER2 -- average BER for selective median filtering case for given SINR - pass 2
%ber_ncFSK_awgn - theoretical BER for non-coherent BFSK for given SINR

%SINR = [ 0 2.5 5 7.5 10 12.5 15]; %in dB
eps=1e-10;
fs=100000; % Sampling frequency is 180 KHz
fc= 25000;% CPFSK center frequency 1is 30 KHz
numbit = 1024;
bitperiod = 0.0005;
%h = 1.0 ; % modulation index
h = 2.0 ; % modulation index
%h = 4.0 ; % modulation index
SINR = [ 0 2.5 5  7.5  10 12.5  15]; %in dB
num_sinr = length(SINR);
num_run=10; %Number of simulation runs for one SINR, use 50 or 100 for statistical reliability 
AV_NF_BER = zeros(1,num_sinr);%average BER of non-filtered FSK signal after num_run runs for given SINR
AV_SMF_BER1 = zeros(1,num_sinr);%average BER of filtered FSK signal after num_run runs for given SINR (pass 1)
AV_SMF_BER2 = zeros(1,num_sinr);%average BER of filtered FSK signal after num_run runs for given SINR (pass 2)
AV_MF_BER = zeros(1,num_sinr);
for ii=1:num_sinr
    SINR_dB = SINR(ii);    
    BER_noise_sum = 0;
    BER_filter1_sum = 0;
    BER_filter2_sum = 0;
    BER_plainmedian_sum = 0;
    for jj=1:num_run
%compute average BER for different SINR  
 
%    [Y r Z fZ ber_nf_av ber_mf_av ]  =  fsk_signal(SINR_dB, num_run);
      [BER_clean BER_noise BER_filter1 BER_filter2 BER_plainmedian]= cpfsk_real_mod_demod(numbit,fs, fc, bitperiod, h, SINR_dB);
      
      BER_noise_sum =  BER_noise_sum +  BER_noise;
      BER_filter1_sum =BER_filter1_sum +  BER_filter1;
      BER_filter2_sum = BER_filter2_sum +  BER_filter2;
      BER_plainmedian_sum = BER_plainmedian_sum + BER_plainmedian;
    end
    %[Y r Z fZ ber_nf_av ber_mf_av ber_mf2_av]  =  fsk_signal2(SINR_dB, num_run);    
    AV_NF_BER(ii) = BER_noise_sum/num_run; %average BER of non-filtered FSK signal for given SINR after num_run runs
    if (AV_NF_BER(ii) == 0) 
          AV_NF_BER(ii) = eps;
    end
    AV_SMF_BER1(ii)= BER_filter1_sum/num_run;%average BER of selectively filtered FSK signal for given SINR after num_run runs (pass 1)
    if (AV_SMF_BER1(ii) == 0) 
          AV_SMF_BER1(ii) = eps;
    end     
    AV_SMF_BER2(ii)= BER_filter2_sum/num_run;%average BER of selectively filtered FSK signal for given SINR after num_run runs (pass 2)
    if (AV_SMF_BER2(ii) == 0) 
          AV_SMF_BER2(ii) = eps;
    end      
    AV_MF_BER(ii)= BER_plainmedian_sum/num_run;%average BER of filtered FSK signal for given SINR after num_run runs
    if (AV_MF_BER(ii) == 0) 
          AV_MF_BER(ii) = eps;
    end      

end
%av_bwr_nf= [31.7871 21.5586 11.0830 5.3643  1.0068 0.0566 0];
%%(unfiltered); TYPICAL VALUES
%av_ber_mf = [8.6680 5.2070 1.6504  0.5020 0.0537  0.0039 0]; %(selective median filtered); TYPICAL VALUES

 
SINR_ab= power(10, SINR./10);
av_ber1_smf = AV_SMF_BER1;
av_ber2_smf = AV_SMF_BER2;
av_ber_mf = AV_MF_BER;
av_ber_nf = AV_NF_BER;
ber_ncFSK_awgn = 0.5*exp(-SINR_ab./2);
figure
loglog(SINR_ab,av_ber_nf./100, 'r*:');
hold on;
loglog(SINR_ab,av_ber1_smf./100, 'bx:' );
hold on;
loglog(SINR_ab,av_ber2_smf./100, 'gd:' )
hold on;
loglog(SINR_ab,av_ber_mf./100, 'ks:' );
hold on
loglog(SINR_ab, ber_ncFSK_awgn, 'm+:' );
title('Comparison of BER for Unfiltered and  Three Filtered Version of Corrupted CPFSK Signal');
legend('Red--Unfiltered','Blue--Selective Median Filter (pass1)', 'Green--Selective Median Filter (pass2)','Black-- Median Filter','Magenta - Non-coherent AWGN BFSK', 'Location','SouthWest');
%legend('Red--Unfiltered','Blue--Selectively Filtered', 'Black-- Filtered','Magenta - Non-coherent AWGN FSK','Location','SouthWest' );
xlabel('Signal to Interference and Noise Ratio in Log scale');
ylabel('Bit Error Rate in Log scale');
hold off;