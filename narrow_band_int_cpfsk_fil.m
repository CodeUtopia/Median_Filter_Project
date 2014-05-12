function [sig_mod] = narrow_band_int_cpfsk_fil(sig, fc, dev,Nspf,Nspsec,span,safetybuf, pad_index, peakThrParam)

% *** narrow_band_int_real_fil.m ***
%
% Usage: [sig_mod] = narrow_band_int_fil(sig, freq1,freq2,Nspf,Nspsec,span,safetybuf, pad_index)
%pad_index takes care of silence and preamble waveform -- future extension
% where:
%
%       sig = input time-series
% dev=  freq (separation frequency,
%freq1 and freq2 are two frequencies 
%to be preserved even if peaks are detected
% Nspf -- no. of smaples per frame 
% Nspsec -- no. of samples per second        
%span -- window length for filtering
%   sig_mod = output time-series
%safetybuf -- a safety buffer around freq1 and freq2 to make sure freq1 and
%freq2 are preserved after filtering

%peakThrParam  -- use this parameter to compute peakThr= floor(Nspf_mod/peakThrParam);
% Returns a time-series whose phase spectrum is the same as the input
% time-series, but whose magnitude spectrum is replaced by some filtered
% magnitude spectrum.  


N1 = length ( sig );
fftsig = fft ( sig );
rs = abs ( fftsig ); %magnitude of FFT of signal
t = angle ( fftsig );%phase of FFT of signal

  
 % span = input('\nEnter the length of the moving average filter >> ');
  
  if span < 1 || span ~= round ( span ),
    error('length must be a positive integer ... ')
  end
  Nspf_mod=Nspf + pad_index; %In case you want to change later

  f = Nspsec/2*linspace(0,1,Nspf_mod/2+1);

% Plot single-sided amplitude spectrum.
  figure
  plot(f,2*abs(rs(1:Nspf_mod/2+1))/Nspf_mod) 
  xlabel('Frequency (Hz)')
  ylabel('Magnitude of FFT Normalized by Length')
  title('Single-Sided Magnitude Spectrum of  Signal+Noise+Interference');    
%   % Plot double-sided amplitude spectrum.
%   f1 = Nspsec*linspace(0,1,Nspf_mod);
%   figure
%   plot(f1,2*abs(rs(1:Nspf_mod))/Nspf_mod) 
%   xlabel('Frequency (Hz)')
%   ylabel('Magnitude of FFT Normalized by Length')
%   title('Double-Sided Magnitude Spectrum of  Signal+Noise+Interference'); 


  
  freq1= floor(fc + dev);
  freq2= floor(fc - dev);

  freq1= floor(Nspf_mod*(freq1/Nspsec)); % FSK frequency in terms of FFT index 
  freq2= floor(Nspf_mod*(freq2/Nspsec)); % FSK frequency in terms of FFT index
  peakThr= floor(Nspf_mod/peakThrParam);%otherwise use a range from 500-2000
 % peakThr= floor(Nspf_mod/25.0);%otherwise use a range from 500-2000
  P=findpeaks_freq(rs, 0.01, peakThr, 5,5); %got this program and renamed it
%   peakThr
%   length(P)
  %P contains all peak values and its position
  %So, P conatins the position (frequency) of all interfer
  size1=size(P(:,2));
  % Inteference removal
  %Make sure that FSK signal frequencies are not touched
  for j=1:size1
          for jj=-safetybuf:safetybuf
            if P(j,2) == freq1 - jj  || P(j,2) == Nspf_mod - freq1 - jj %Do not touch freq 1 and its reflection
              P(j,2)=0;
            end
          end

          for jj=-safetybuf:safetybuf
            if P(j,2) == freq2 - jj || P(j,2) == Nspf_mod - freq2 - jj  %Do not touch freq 2 and its reflection
              P(j,2)=0;
            end
          end           
 
  end
  
                    

  
  for i = span+1:Nspf_mod-span
      for j=1:size1
          
          if P(j,2)==i 

           R=0;
           for kk=-span:span
             r=rs(1,i-kk);
             R=[R r];
           end
           rs(1,i)=median(R);
% Not doing median filtering on phase as it is harmful
%            R=0;
%            for kk=-span:span
%              r=t(1,i-kk);
%              R=[R r];
%            end
%           t(1,i)=median(R); %filtering the phase, CAREFUL
          end
      end
  end


% Plot single-sided amplitude spectrum.
  f = Nspsec/2*linspace(0,1,Nspf_mod/2+1);
  figure
  plot(f,2*abs(rs(1:Nspf_mod/2+1))/Nspf_mod) 
  xlabel('Frequency (Hz)')
  ylabel('Magnitude of FFT Normalized by Length')
  title('Single-Sided Magnitude Spectrum of Filtered Signal+Noise+Interference');
% 
%   % Plot double-sided amplitude spectrum.
%   f1 = Nspsec*linspace(0,1,Nspf_mod);
%   figure
%   plot(f1,2*abs(rs(1:Nspf_mod))/Nspf_mod) 
%   xlabel('Frequency (Hz)')
%   ylabel('Magnitude of FFT Normalized by Length')
%   title('Double-Sided Magnitude Spectrum of Filtered Noisy Signal');
    sig_mod = ifft ( rs.* exp ( 1i * t ) );
%   
%   freq1
%   freq2
if isreal ( sig ), sig_mod = real ( sig_mod ); end


