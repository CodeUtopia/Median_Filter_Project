function [sig_mod] = narrow_band_int_cpfsk_plainmedian(sig, fc, dev,Nspf,Nspsec,span,safetybuf, pad_index)

% *** narrow_band_int_fil.m -- plain median filtering***
%
% Usage: [sig_mod] = narrow_band_int_fil(sig, freq,,Nspf,Nspsec,span,pad_index)
%pad_index takes care of silence and preamble waveform -- future extension
% where:
%
%       sig = input time-series
%freq1 and freq2 are two frequencies (fc+dev) and (fc-dev)
%to be preserved even if peaks are detected
% Nspf -- no. of smaples per frame 
% Nspsec -- no. of samples per second        
%span -- window length for filtering
%safetybuf -- a safety buffer around freq1 and freq2 to make sure freq1 and
%freq2 are preserved after filtering

%No peak detection required in this implementation
% ALL data points (except FSK signal frequencies) belonging to magnitude spectrum will be filtered
%
% Returns a time-series whose phase spectrum is the same as the input
% time-series, but whose magnitude spectrum is replaced by some filtered
% magnitude spectrum.  


N1 = length ( sig );
% N2 = fix ( ( N + 1 ) / 2 );
fftsig = fft ( sig );
rs = abs ( fftsig );
t = angle ( fftsig );
  
 % span = input('\nEnter the length of the moving average filter >> ');
  
  if span < 1 || span ~= round ( span ),
    error('length must be a positive integer ... ')
  end
  Nspf_mod=Nspf + pad_index;

%axis([ 1 count -0.001 -0.25]);
  freq1= floor(fc + dev);
  freq2= floor(fc - dev);

% compute frequencies of FSK in terms of FFT indices
  freq1= floor(Nspf_mod*(freq1/Nspsec)); % FSK frequency in terms of FFT index  
  freq2= floor(Nspf_mod*(freq2/Nspsec)); % FSK frequency in terms of FFT index   
  freq1fb = Nspf_mod - freq1;
  freq2fb = Nspf_mod - freq2;
                    

  P=1; %index that should be one if the data point is to be filtered
  Q=1;%index that should be one if the data point is to be filtered
  for i = span+1:Nspf_mod-span
  
          for jj=-safetybuf:safetybuf
            if (i == freq1 - jj || i == freq2 - jj) %Do not touch freq 1 and freq 2
              P=0;
              break;
            else
              P=1;
            end
          end
         for jj=-safetybuf:safetybuf
            if (i == freq1fb - jj || i == freq2fb - jj ) %Do not touch freq 1 and freq 2 from end (fs)
              Q=0;
              break;
            else
              Q=1;
            end
          end         

        if (P==1 && Q ==1)
           R=0;
           for kk=-span:span
             r=rs(1,i-kk);
             R=[R r];
           end
           rs(1,i)=median(R);
           %temp = D_filter(R, 2*span+1);
           %rs(1,i) = temp(1, span+1);
           R=0;
           for kk=-span:span
             r=t(1,i-kk);
             R=[R r];
           end
%           t(1,i)=median(R); %filtering the phase, CAREFUL
        end

  end
  

%   f1 = Nspsec*linspace(0,1,Nspf_mod);
%   figure
%   plot(f1,2*abs(rs(1:Nspf_mod))/Nspf_mod) 
%   xlabel('Frequency (Hz)')
%   ylabel('Magnitude of FFT Normalized by Length')
%   title('Double-Sided Magnitude Spectrum of Plain Median Filtered Noisy Signal'); 
   sig_mod = ifft ( rs.* exp ( 1i * t ) );
  
%   freq1
%   freq2
if isreal ( sig ), sig_mod = real ( sig_mod ); end


