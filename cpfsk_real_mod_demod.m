function [BER_clean BER_noise BER_fil1 BER_fil2 BER_plainmedian ]=cpfsk_real_mod_demod(numbit,fs, fc, bitperiod, h, SINR_dB)
% Make it a function like -- [BER_clean BER_noise BER_fil]= fsk_cpfsk_mod_demod(numbit,fs, bitperiod, hparam, fc, SINR_dB);
%clear all
close all
%numbit = 1024;   %numbit -- number of bits
%fs=100000;       %fs -- sampling frequency
%fc=25000;   %center frequency
%bitperiod=.0005;   %bitperiod = duration of one bit -> 2000 bps 
%bit = randsrc(1,numbit); % bit stream of random 1's and 0's
%h is the modulation index
bit = round(rand(1,numbit)); 
% bit(find(bit==0)) = -1;
bit = 2 * bit - 1;% bit stream of random 1's and  -1's (for 0's)
delt=1/fs;
%dev=.5/bitperiod;    %frequency deviation
hparam= h;
dev=hparam /(2*bitperiod);    %frequency deviation  
% create the FSK signal
L=fs*bitperiod; %number of samples in one bit, 50 now
L1=numbit*L;  % total number of samples per frame
t=[0:delt:(L1-1)*delt];     %total time array
bitup=upsample(bit,50);
bitup=conv(bitup,ones(1,50));   %rectangular modulation
theta=2*pi*dev*cumsum(bitup)*delt;     %phase process, piecewise linear and continuous-phase
% plot(theta)
% title('piecewise linear phase function')
    %now the real bandpass transmitted signal
transmitted=cos(2*pi*fc*t+theta(1:51200));%51200 is equal to L1
% figure(2)
% plot(transmitted)
% title('CP FSK Signal')
% figure
% temp=abs(fft(transmitted)).^2;
% %temp=abs(fft(transmitted));
% %psd=temp(1:25600);  %plot half the spectrum
% psd=temp(1:51200);  %plot half the spectrum
% title('Power spectral density of FSK signal only')
% semilogy(psd)
%  xlabel('Frequency (Hz)')
%  ylabel('Magnitude of FFT Normalized by Length')
%  title(' Magnitude Spectrum of  Signal'); 
  
 % Now add interference and noise
 %SINR_dB=0;
 SINR=10.^(SINR_dB/10); 
% INR=100; %interference to noise ratio, interference dominated
 INR=25; %interference to noise ratio, interference dominated, 
 %according to Prof. Wilson
 SNR=SINR*(1+INR); %signal-to-noise ratio
 SIR=SNR./(INR+eps); %signal-to-interference ratio
 SIRdB = 10 * log10(SIR); % Signal to interference in dB
 
 P=1;
 num_int=30; %number of interfering sinusoids 
 sigma_z=sqrt(L*P./SNR); %thermal noise standard deviation
 sigma_i=sqrt(L*P./SIR); %interference standard deviation
 z=sigma_z*(randn(1,L1));  %generate white additive random noise
 fn=0.5 * fs *[.3*rand(10,1);.3+.4*rand(10,1);.7+.3*rand(10,1)]; %interferer frequencies
 A=rand(1,num_int);
 A=sigma_i * A/norm(A); %interferer amplitudes
% %Create interference if present
in=zeros(1,length(t));
for m=1:num_int
    in=in+A(m)*cos(2*pi*fn(m)*t);
end
%Now, create the corrupted signal 
 r=transmitted +in +z; % FINALLY got the corrupted signal
%  temp=abs(fft(r)).^2;
%  %temp=abs(fft(r)); 
%  psd2=temp(1:25600);
%  figure
%  semilogy(psd2)
%  xlabel('Frequency (Hz)')
%  ylabel('Magnitude of FFT Normalized by Length')
%  title(' Magnitude Spectrum of  Noise Corrupted Signal'); 
% now filter the signal, 21-point median filter is used
 [sig_mod] = narrow_band_int_cpfsk_fil(r, fc, dev,L1,fs,10,25, 0, 40.0); %first pass of selective median filter
 [sig_mod1] = narrow_band_int_cpfsk_fil(sig_mod, fc, dev,L1,fs,10,25, 0, 40.0);%second pass of selective median filter
 [sig_mod2] = narrow_band_int_cpfsk_plainmedian(r, fc, dev,L1,fs,10,25, 0);% plain median filtering
% rF= fft(r);
% %length(rF)
% figure
% plot(fscale,2*abs(rF(1:L1/2 + 1))/L1) 
% xlabel('Frequency (Hz)')
% ylabel('Magnitude of FFT Normalized by Length')
% title('Single-Sided Magnitude Spectrum of  Signal+Noise+Interference');

%%%%%%%%%%%%%%%%%DEMODULATION BELOW HERE

% CPFSK Demodulation using non-coherent quadratic receiver

%decode clean data
 [demod_output] = cpfsk_demod(L,delt,fc,dev,transmitted, numbit);
%  figure 
%  stem(demod_output);
%  title('Demodulated Bits');
 BER_clean = sum(abs(demod_output - bit))/(2*numbit); %factor 2 needed
% Decode noisy data
 [demod_output] = cpfsk_demod(L,delt,fc,dev,r, numbit) ;
%  figure 
%  stem(demod_output);
%  title('Demodulated Bits');
 BER_noise = sum(abs(demod_output - bit))/(2*numbit);%factor 2 needed

%decode filtered data after first pass of selective median filtering
 [demod_output] = cpfsk_demod(L,delt,fc,dev,sig_mod, numbit);
%  figure 
%  stem(demod_output);
%  title('Demodulated Bits');
 BER_fil1 = sum(abs(demod_output - bit))/(2*numbit); %factor 2 needed
 
%decode filtered data after second pass of selective median filtering 
 [demod_output] = cpfsk_demod(L,delt,fc,dev,sig_mod1, numbit);
 BER_fil2 = sum(abs(demod_output - bit))/(2*numbit); %factor 2 needed
 
%decode filtered data after plain median filtering 
 [demod_output] = cpfsk_demod(L,delt,fc,dev,sig_mod2, numbit);
 BER_plainmedian = sum(abs(demod_output - bit))/(2*numbit); %factor 2 neede  
end
% 
function [demod_output] = cpfsk_demod(L,delt,fc,dev,transmitted, numbit)
%Optimal Non-coherent Detection Using Quadratic Receiver of CP-BFSK signal
%   Detailed explanation goes here
% CP-BFSK Demodulation
tb=[0:delt:(L-1)*delt];     %time array for 1 bit
len_bit=size(tb, 2);
bit_pone = ones(1, len_bit); %Numbers of 1's within a bit of identity 1
bit_none = -1 * ones(1, len_bit); %Numbers of -1's within a bit of identity -1
theta1=2*pi*dev*cumsum(bit_pone)*delt;
theta2=2*pi*dev*cumsum(bit_none)*delt;
%Create four basis functions of a bit
cosfc1 = transpose(cos(2*pi*fc*tb+theta1));
cosfc2 = transpose(cos(2*pi*fc*tb+theta2));
sinfc1 = transpose(sin(2*pi*fc*tb+theta1));
sinfc2 = transpose(sin(2*pi*fc*tb+theta2));


%Decode transmitted data
 for x = 1:numbit
     s1=transmitted((x-1)*L+1: x*L) * cosfc1;
     s2=transmitted((x-1)*L+1: x*L) * sinfc1;
     S1 = sqrt((s1)^2 + (s2)^2);
     s3=transmitted((x-1)*L+1: x*L) * cosfc2;
     s4=transmitted((x-1)*L+1: x*L) * sinfc2;
     S2 = sqrt((s3)^2 + (s4)^2);     

     if S1 > S2
        demod_output(x)=1;
     else
         demod_output(x)=-1;
     end
 end

end




