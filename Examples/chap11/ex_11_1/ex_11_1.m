% Example 11.1 posterior plotting program
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
%
m=(5.0:0.1:15.0)';
for i=1:length(m)
  p1(i)=(1/sqrt(2*pi))*exp(-(10.3-m(i))^2/2);
  p2(i)=(1/sqrt(pi))*exp(-(10.2-m(i))^2);
end

figure(1)
clf
plot(m,p1,'k');
xlabel('m ({\mu}g)');
ylabel('q(m|d=10.3 {\mu}g)');
bookfonts
axis([5 15 0 0.7]);
print -deps c11fpost1.eps

figure(2)
clf
plot(m,p2,'k');
xlabel('m ({\mu}g)');
ylabel('q(m|d_{1}=10.3 \mu g, d_{2}=10.1 {\mu}g)');
bookfonts
axis([5 15 0 0.7]);
print -deps2 c11fpost2.eps

