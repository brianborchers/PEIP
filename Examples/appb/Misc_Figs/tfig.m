%demonstrate t-distribution
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
clear
x=-4:.01:4;
figure(1)
clf
bookfonts

for nu=[3 20]
    plot(x,tpdf(x,nu),'k')
    hold on
end
plot(x,normpdf(x,0,1),'k--');
hold off
ylabel('f_t(x)')
xlabel('x')
text(-2.8,0.01,'\nu=20','Fontsize',20)
text(-2.8,0.06,'\nu=3','Fontsize',20)
text(0.5,0.38,'N(0, 1)','Fontsize',20)

print -deps2 abftpdf.eps
