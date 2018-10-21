function [h,ww] = freqs(b,a,w)
% FREQS Laplace-transform (s-domain) frequency response.  
%   H = FREQS(B,A,W) returns the complex frequency response vector H 
%   of the filter B/A:
%                        nb-1         nb-2
%            B(s)   b(1)s     +  b(2)s     + ... +  b(nb)
%     H(s) = ---- = -------------------------------------
%                        na-1         na-2
%            A(s)   a(1)s     +  a(2)s     + ... +  a(na)
%
%   given the numerator and denominator coefficients in vectors B and A. 
%   The frequency response is evaluated at the points specified in 
%   vector W (in rad/s).  The magnitude and phase can be graphed by
%   calling FREQS(B,A,W) with no output arguments.
%
%   [H,W] = FREQS(B,A) automatically picks a set of 200 frequencies W on 
%   which the frequency response is computed.  FREQS(B,A,N) picks N 
%   frequencies. 
%
%   See also LOGSPACE, POLYVAL, INVFREQS, and FREQZ.

% 	Author(s): J.N. Little, 6-26-86
%   	   T. Krauss, 3-19-93, default plots and frequency vector
%   Copyright 1988-2004 The MathWorks, Inc.
%   Revision: 1.12.4.2 $  $Date: 2004/12/26 22:15:54 $


error(nargchk(2,3,nargin));
error(nargoutchk(0,2,nargout));

if ~any(size(b)<=1) || ~any(size(a)<=1),
    error('The numerator and denominator must be vectors.');
end

if nargin == 2,
    w = 200;
end
if length(w) == 1,
    wlen = w;
    w_long = freqint(b,a,wlen);
    % need to interpolate long frequency vector:
    xi = linspace(1,length(w_long),wlen).';
    w = 10.^interp1(1:length(w_long),log10(w_long),xi,'linear');
end

s = j*w;
hh = polyval(b,s) ./ polyval(a,s);

if nargout == 0,
    newplot;
    mag = abs(hh);   phase = angle(hh)*180/pi;
    subplot(211),loglog(w,mag),set(gca,'xgrid','on','ygrid','on') 
    set(gca,'xlim',[w(1) w(length(w))])
    xlabel('Frequency (rad/s)')
    ylabel('Magnitude')
    ax = gca;
    subplot(212), semilogx(w,phase),set(gca,'xgrid','on','ygrid','on') 
    set(gca,'xlim',[w(1) w(length(w))])
    xlabel('Frequency (rad/s)')
    ylabel('Phase (degrees)')
    axes(ax)
elseif nargout == 1,
    h = hh;
elseif nargout == 2,
    h = hh;
    ww = w;
end
% end freqs

function w=freqint(a,b,c,d,npts)
%FREQINT Auto-ranging algorithm for Bode frequency response
%   W=FREQINT(A,B,C,D,Npts)
%   W=FREQINT(NUM,DEN,Npts)

%   Andy Grace  7-6-90
%   Was Revision: 1.9,  Date: 1996/07/25 16:43:37

% Generate more points where graph is changing rapidly.
% Calculate points based on eigenvalues and transmission zeros. 

na = size(a, 1);

if (nargin==3) && (na==1),		% Transfer function form.
  npts=c;
  ep=roots(b);
  tz=roots(a);
else				% State space form
  if nargin==3, npts=c; [a,b,c,d] = tf2ss(a,b); end
  ep=eig(a);
  tz=tzero(a,b,c,d);
end

if isempty(ep), ep=-1000; end

% Note: this algorithm does not handle zeros greater than 1e5
ez=[ep(find(imag(ep)>=0));tz(find(abs(tz)<1e5&imag(tz)>=0))];

% Round first and last frequencies to nearest decade
integ = abs(ez)<1e-10; % Cater for systems with pure integrators
highfreq=round(log10(max(3*abs(real(ez)+integ)+1.5*imag(ez)))+0.5);
lowfreq=round(log10(0.1*min(abs(real(ez+integ))+2*imag(ez)))-0.5);

% Define a base range of frequencies
diffzp=length(ep)-length(tz);
w=logspace(lowfreq,highfreq,npts+diffzp+10*(sum(abs(imag(tz))<abs(real(tz)))>0));
ez=ez(imag(ez)>abs(real(ez)));

% Oscillatory poles and zeros
if ~isempty(ez)
  f=w;
  npts2=2+8/ceil(abs((diffzp+eps)/10));
  [dum,ind]=sort(-abs(real(ez))); %#ok
  z=[];
  for i=ind'
    r1=max([0.8*imag(ez(i))-3*abs(real(ez(i))),10^lowfreq]);
    r2=1.2*imag(ez(i))+4*abs(real(ez(i)));
    z=z(z>r2|z<r1);
    f=f(f>r2|f<r1);
    z=[z,logspace(log10(r1),log10(r2),sum(w<=r2&w>=r1)+npts2)];
  end
  w=sort([f,z]);
end
w = w(:);

% end freqint
