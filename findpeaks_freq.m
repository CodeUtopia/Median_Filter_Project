function P=findpeaks_freq(s, SlopeThreshold,AmpThreshold,smoothwidth,peakgroup)
% Function to locate the positive peaks in a noisy x-y data set
% Detects peaks by looking for downward zero-crossings
% in the first derivative that exceed SlopeThreshold.
% Returns list (P) containing peak number and
% position, height, and width of each peak. SlopeThreshold,
% AmpThreshold, and smoothwidth control sensitivity
% Higher values will neglect smaller features. Peakgroup
% is the number of points around the "top part" of the peak.
% T. C. O'Haver, 1995.  Version 2  Last revised Oct 27, 2006
% Data x and y are generated from input data s (or loaded  file)

%If you are loading the data, load the data in and call it variable 's' as shown

%s = load(filename);
max1 = max(s(1,:));
% max1  %checking the max
size(s) % checking the size
max1=1;
%option for normalization
xs=size(s,2);
for cl=1:size(s,2)
     y(cl) = s(1,cl)/max1;
     if ( y(cl) < 0)
        y(cl) = 0; 
     end
end

x=zeros(1,xs);
x(1)=1;
for cl=2:size(s,2)
     x(cl) = x(cl-1)+1;
end

%At this point, y contains the signal and X contains the index

smoothwidth=round(smoothwidth);
peakgroup=round(peakgroup);
d=fastsmooth(deriv(y),smoothwidth);
%d=-fastsmooth(deriv(d1),smoothwidth);
%d=fastsmooth(deriv(d2),smoothwidth);
n=round(peakgroup/2+1);
P=[0 0 0 0];

vectorlength=length(y);
%vectorlength=length(y)  %check out vector length
peak=1;
AmpTest=AmpThreshold;
%AmpTest=AmpThreshold   %check out threshold 
for j=smoothwidth:length(y)-smoothwidth,
   if sign(d(j)) > sign (d(j+1)), % Detects zero-crossing
     if d(j)-d(j+1) > SlopeThreshold*y(j), % if slope of derivative is larger than SlopeThreshold
        if y(j) > AmpTest,  % if height of peak is larger than AmpThreshold
          for k=1:peakgroup, % Create sub-group of points near peak
              groupindex=j+k-n+1;
              if groupindex<1, groupindex=1;end
              if groupindex>vectorlength, groupindex=vectorlength;end
            xx(k)=x(groupindex);yy(k)=y(groupindex);
          end
          [coef,S,MU]=polyfit(xx,log(abs(yy)),2);  % Fit parabola to log10 of sub-group with centering and scaling
          c1=coef(3);c2=coef(2);c3=coef(1);
          PeakX=-((MU(2).*c2/(2*c3))-MU(1));   % Compute peak position and height of fitted parabola
          PeakY=exp(c1-c3*(c2/(2*c3))^2);
 %         print PeakX, PeakY;
          MeasuredWidth=norm(MU(2).*2.35703/(sqrt(2)*sqrt(-1*c3)));
          
          if peakgroup<7,
             PeakY=max(yy);
             pindex=val2ind(yy,PeakY);
             PeakX=xx(pindex(1));
          end
          
          % Construct matrix P. One row for each peak 
          % detected, containing the peak number, peak 
          % position (x-value) and peak height (y-value).
          P(peak,:) = [round(peak) PeakX PeakY MeasuredWidth];
          peak=peak+1;
        end
      end
   end
end

function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);

function SmoothY=fastsmooth(Y,smoothwidth)
%  fastsmooth(Y,w) smooths vector Y by triangular
% smooth of width = smoothwidth. Works well with signals up to 
% 100,000 points in length and smooth widths up to 1000 points. 
% Faster than tsmooth for smooth widths above 600 points.
% Example: fastsmooth([0 0 0 0 9 0 0 0 0],3) yields [0 0 1 2 3 2 1 0 0]
%  T. C. O'Haver, 2006.
w=round(smoothwidth);
SumPoints=sum(Y(1:w));
s=zeros(size(Y));
halfw=round(w/2);
for k=1:length(Y)-w,
   s(k+halfw-1)=SumPoints;
   SumPoints=SumPoints-Y(k);
   SumPoints=SumPoints+Y(k+w);
end
s=s./w;
SumPoints=sum(s(1:w));
SmoothY=zeros(size(s));
for k=1:length(s)-w,
   SmoothY(k+halfw-1)=SumPoints;
   SumPoints=SumPoints-s(k);
   SumPoints=SumPoints+s(k+w);
end
SmoothY=SmoothY./w;

function d=deriv(a)
% First derivative of vector using 2-point central difference.
% Example: deriv([1 1 1 2 3 4]) yeilds [0 0 .5 1 1 1]
%  T. C. O'Haver, 1988.
n=length(a);
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
  d(j)=(a(j+1)-a(j-1)) ./ 2;
end