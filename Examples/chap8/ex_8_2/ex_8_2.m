% Example 8.2
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber

% make sure we have a clean environment
clear;
randn('seed',0);

% set up the basic problem parameters
noise=0.05;
N=211;
M=211;
dt=.5;
dF=1/dt;
t=linspace(-5,100,N);
F=linspace(0,dF/2,N-1);

%instrument response is a critically-damped pulse
sigi=10;
for i=1:N-1
  if (t(i)<0)
    g(i) = 0;
  else
    g(i) = t(i)*exp(-t(i)/sigi);
  end
end
gmax=max(g);
g=g'/gmax;

%true signal is two pulses of sig standard deviation
sig=2;
mtrue = (exp(-((t(1:N-1)-8).^2/(sig^2*2)))+...
    0.5*exp(-((t(1:N-1)-25).^2/(sig^2*2))))';
%the model should be unit height
mtrue=mtrue/max(mtrue);

% get the true and noisy data
d=conv(g,mtrue);
dn=d+noise*randn(length(d),1);

% plot the true model
figure(1)
clf
plot(t(1:N-1),mtrue,'k')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
bookfonts
axis tight

disp('Displaying the true model (fig. 1)')

%plot the impulse response
figure(2)
clf
plot(t(1:N-1),g,'k')
xlabel('Time (s)')
ylabel('V')
bookfonts
axis tight

disp('Displaying the instrument impulse response (fig. 2)')

% plot the noise free data
figure(3)
clf
plot(t(1:N-1),d(1:N-1),'k')
xlabel('Time (s)')
ylabel('V')
bookfonts
axis tight

disp('Displaying the noise free data (fig. 3)')

% plot the noisy data
figure(4)
clf
plot(t(1:N-1),dn(1:N-1),'k')
xlabel('Time (s)')
ylabel('V')
bookfonts
axis tight

disp('Displaying the noisy data (fig. 4)')

%generate spectra
% pad g to 2*(N-1) points to avoid wrap-around effects
% append a 0 to the data series since they are one shy of this length
dspec=fft([d; 0]);
dnspec=fft([dn; 0]);
gspec=fft([g; zeros(length(g),1)]);

% do the inverse convolution
mperf=real(ifft(dspec./gspec));
mn=real(ifft(dnspec./gspec));

% plot the model recovered from the noise free data
figure(5)
clf
plot(t(1:N-1),mperf(1:N-1),'k')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
bookfonts
axis tight

disp('Displaying the model recovered from the noise free data (fig. 5)')

% plot the model recovered from the noisy data
figure(6)
clf
plot(t(1:N-1),mn(1:N-1),'k')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
bookfonts
axis tight

disp('Displaying the model recovered from the noisy data (fig. 6)')

% Get the nonrepeated parts of the spectra
gs=abs(gspec(1:length(gspec)/2));
ds=abs(dspec(1:length(dspec)/2));
dns=abs(dnspec(1:length(dspec)/2));

% plot the spectra for the impulse response, noise free and noisy data
figure(7)
clf
loglog(F,gs,'k',F,ds,'k-.',F,dns,'k-');
xlabel('f (Hz)')
ylabel('Spectral Amplitude')
bookfonts
ylim([1e-3 1e2]);
xlim([1/N 1]);

disp('Displaying the spectra for data and impulse response (fig. 7)')

% plot the spectra of the noisy data divided by the spectra of the instrument
figure(8)
clf
loglog(F,dns./gs,'k-')
xlabel('f (Hz)')
ylabel('Spectral Amplitude')
bookfonts
axis tight

disp(['Displaying the spectra of the noisy data divided by'...
    ' the instrument (fig. 8)'])


disp('Displaying an animation of zeroth-level level regularization (fig. 9)')
wcount=0;
for alpha=10.^(-2:.02:1)
  wcount=wcount+1;

  %apply zeroth order Tikhonov regularization
  M=(conj(gspec).*dnspec)./(conj(gspec).*gspec+alpha^2*ones(size(gspec)));

  % get the temporal model and store the necessary information
  mw=real(ifft(M));
  aval(wcount)=alpha;
  mnorm(wcount)=norm(M);
  resid(wcount)=norm(gspec.*M-dnspec);
  mstore(:,wcount)=mw;

  %animation of the evolving solution (not saved)
  figure(9)
  clf
  plot(t(1:N-1),mw(1:N-1),'k')
  xlabel('Time (s)')
  ylabel('Acceleration (m/s^2)')
  bookfonts
  drawnow;
  pause(0.1);
end

% plot the L-curve with respect to zeroth-order regularization
figure(10)
clf
plot(resid,mnorm,'k-')
xlabel('Residual Norm ||GM-D||_2');
ylabel('Solution Norm ||M||_2');
bookfonts
% label some of the alphas
for i=[20,85,140]
  H=text(resid(i)+5,mnorm(i)+10,num2str(aval(i),'%3.1f'));
  set(H,'Fontsize',18);
end
hold on
H=plot(resid(85),mnorm(85),'ko');
set(H,'MarkerSize',14);
hold off
axis tight

disp(['Displaying the L-curve for zeroth-order Tikhonov regularization'...
    ' (fig.  10)'])
print -deps2 c8fmalphatradeo_t0.eps

% plot the suite of solutions
figure(11)
clf
hold on
for i=1:10:length(aval)
  plot(t(1:N-1),mstore(1:N-1,i)/10+log10(aval(i)),'k')
  plot(t(1:N-1),mtrue/10+log10(aval(i)),'-.k')
end
% highlight the selected solution
plot(t(1:N-1),mstore(1:N-1,85)/10+log10(aval(85)),'k','LineWidth',2);
xlabel('Time (s)')
ylabel('log_{10} (\alpha)')
bookfonts
axis tight

disp('Displaying the suite of zeroth-order regularized models (fig. 11)')
print -deps2 c8fmalpharange_t0.eps

%plot the best solution from L-curve plot
figure(12)
clf
plot(t(1:N-1),mstore(1:N-1,85),'k')
hold on
H=plot(t(1:N-1),mtrue(1:N-1),'-.k');
set(H,'LineWidth',1);
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
bookfonts
axis tight

disp('Displaying the preferred zeroth order regularized model (fig. 12)')
print -deps2 c8fmalphabest_t0.eps

% second order regularization
wcount=0;
%find frequencies in spectra for Tikhonov weighting
f=fftshift((0:length(gspec)-1)'/length(gspec))-0.5;

disp('Displaying an animation of second-level level regularization (fig. 13)')
for alpha=sqrt(2*pi)*10.^(-1:.02:3)
  wcount=wcount+1;
  
  %apply second order Tikhonov regularization (the factor of 2 pi is just a
  %scaling on alpha, effectively.
  M=(conj(gspec).*dnspec)./(conj(gspec).*gspec+alpha^2*(f.^4));

  % get the temporal model and store the necessary information
  mw=real(ifft(M));
  aval(wcount)=alpha;
  mnorm(wcount)=norm(M);
  resid(wcount)=norm(gspec.*M-dnspec);
  mstore(:,wcount)=mw;
  
  %animation of the evolving solution (not saved)
  figure(13)
  clf
  plot(t(1:N-1),mw(1:N-1),'k')
  xlabel('Time (s)')
  ylabel('Acceleration (m/s^2)')
  bookfonts
  drawnow;
  pause(0.1);
end

% plot the L-curve with respect to second-order regularization
figure(14)
clf
plot(resid,mnorm,'k-')
xlabel('Residual Norm ||GM-D||_2');
ylabel('Model Seminorm ||K^4 M||_2');
bookfonts
% label some of the alphas
for i=[10,125,180]
  H=text(resid(i)*1.03,mnorm(i)*1.1,num2str(aval(i),'%3.1f'));
  set(H,'Fontsize',18);
end
hold on
H=plot(resid(125),mnorm(125),'ko');
hold off
set(H,'MarkerSize',14);
axis tight

disp(['Displaying the L-curve for second-order Tikhonov regularization'...
    ' (fig.  14)'])
print -deps2 c8fmalphatradeo_t2.eps

% plot the suite of solutions
figure(15)
clf
hold on
for i=1:10:length(aval)
  plot(t(1:N-1),mstore(1:N-1,i)/10+log10(aval(i)),'k')
  plot(t(1:N-1),mtrue/10+log10(aval(i)),'-.k')
end
% highlight the selected solution
plot(t(1:N-1),mstore(1:N-1,125)/10+log10(aval(125)),'k','LineWidth',2);
hold off
xlabel('Time (s)')
ylabel('log_{10} (\alpha)')
bookfonts
axis tight

disp('Displaying the suite of second-order regularized models (fig. 15)')
print -deps2 c8fmalpharange_t2.eps

%plot the best solution from L-curve plot
figure(16)
clf
plot(t(1:N-1),mstore(1:N-1,125),'k')
hold on
H=plot(t(1:N-1),mtrue(1:N-1),'-.k');
set(H,'LineWidth',1);
hold off
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
bookfonts
axis tight

disp('Displaying the preferred second order regularized model (fig. 16)')
print -deps2 c8fmalphabest_t2.eps
