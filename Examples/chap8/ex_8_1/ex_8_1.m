% Example 8.1
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

% instrument response is a critically-damped pulse
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

% true signal is two pulses of sig standard deviation
sig=2;
mtrue = (exp(-((t(1:N-1)-8).^2/(sig^2*2)))+...
    0.5*exp(-((t(1:N-1)-25).^2/(sig^2*2))))';
% the model should be unit height
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

% generate spectra
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
print -deps2 c8fspectra1.eps

% plot the spectra of the noisy data divided by the spectra of the instrument
figure(8)
clf
loglog(F,dns./gs,'k-')
xlabel('f (Hz)')
ylabel('Spectral Amplitude')
bookfonts
axis tight

disp(['Displaying the spectra of the noisy data divided by'...
    ' the instrument (fig.  8)'])
print -deps2 c8fspectra2.eps

disp('Displaying an animation of the water level regularization (fig. 9)')
% apply water level regularization to the model
wcount=0;
for wlev=10.^(-1:.1:2)
  wcount=wcount+1;

  % apply water level damping
  for i=1:length(gspec)
    if abs(gspec(i)) < wlev
      gwspec(i,1)=wlev*(gspec(i)./abs(gspec(i)));
    else
      gwspec(i,1)=gspec(i);
    end
  end

  % get the temporal model and store the necessary information
  mw=real(ifft(dnspec./gwspec));
  wval(wcount)=wlev;
  mnorm(wcount)=norm(mw);
  resid(wcount)=norm(real(ifft((dnspec./gwspec).*gspec-dnspec)));
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

% plot the L-curve with respect to water level regularization
figure(10)
clf
plot(resid,mnorm,'k')
xlabel('Residual Norm ||Gm -d||_2');
ylabel('Solution Norm ||m||_2');
bookfonts
% label some of the alphas
for i=2:7:wcount
  H=text(resid(i)+1,mnorm(i)+.3,num2str(wval(i),'%3.1f'));
  set(H,'Fontsize',18);
end

disp('Displaying the water level L-curve (fig. 10)')
print -deps2 c8fmwtradeo.eps

% plot suite of models recovered
figure(11)
clf
hold on
for i=1:31
  plot(t(1:N-1),mstore(1:N-1,i)+10*log10(wval(i)),'k')
  plot(t(1:N-1),mtrue+10*log10(wval(i)),'-.k')
end
% highlight the selected solution
plot(t(1:N-1),mstore(1:N-1,16)+10*log10(wval(16)),'k','LineWidth',2);
xlabel('Time (s)')
ylabel('10 log_{10} (w)')
bookfonts
axis tight

disp('Displaying the suite of regularized models (fig. 11)')
print -deps2 c8fmwrange.eps

% plot the desired regularized model
figure(12)
clf
plot(t(1:N-1),mstore(1:N-1,16),'k')
hold on
plot(t(1:N-1),mtrue(1:N-1),'-.k')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')

disp('Displaying the preferred regularized model (fig. 12)')
bookfonts
axis tight
print -deps2 c8fmwbest.eps
