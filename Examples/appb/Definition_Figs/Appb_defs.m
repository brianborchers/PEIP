%generate distribution figures for Appendix B
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
clear
figure(1)
clf
bookfonts
x=-0.5:0.01:1.5;
f_U=zeros(length(x));
ind=find(x>=0 & x <=1);
f_U(ind)=1;
plot(x,f_U,'k')
xlabel('x')
ylabel('f_U (x)')
ylim([0 1.2])
print -deps2 abfunipdf.eps

figure(2)
clf
bookfonts
x=-4:0.01:4;
f_N=(1/sqrt(2*pi))*exp(-0.5*x.^2);
plot(x,f_N,'k')
xlabel('x')
ylabel('f_N (x)')
ylim([0 0.5])
print -deps2 abfnormpdf.eps

figure(3)
clf
bookfonts
x=-1:0.01:4;
f_exp=exp(-x);
ind=find(x<0)
f_exp(ind)=0;
plot(x,f_exp,'k')
xlabel('x')
ylabel('f_{exp} (x)')
ylim([0 1.2])
print -deps2 abfexppdf.eps

figure(4)
clf
bookfonts
x=-4:0.01:4;
f_dexp=(1/sqrt(2))*exp(-sqrt(2)*abs(x));
plot(x,f_dexp,'k')
xlabel('x')
ylabel('f_{dexp} (x)')
ylim([0 1])
print -deps2 abfdexppdf.eps

figure(5)
clf
bookfonts
x=0:0.01:15;
f_chi2=zeros(length(x),4);
i=1;
for nu=[3,5,7,9]
f_chi2(:,i)=chi2pdf(x,nu)
i=i+1;
end
for i=1:4
    plot(x,f_chi2(:,i),'k');
    hold on
end
hold off
xlabel('x')
ylabel('f_{\chi^2}(x)')
H=text(3,0.2,'\nu=3');
set(H,'FontSize',14);
H=text(4.5,0.15,'\nu=5');
set(H,'FontSize',14);
H=text(6,0.13,'\nu=7');
set(H,'FontSize',14);
H=text(8.5,0.11,'\nu=9');
set(H,'FontSize',14);
print -deps2 abfchi2pdf.eps

figure(6)
clf
bookfonts
x=-4:0.01:4;
f_t=zeros(length(x),3);
i=1;
for nu=[1,3,100]
f_t(:,i)=tpdf(x,nu)
i=i+1;
end
for i=1:3
    plot(x,f_t(:,i),'k');
    hold on
end
hold off
xlabel('x')
ylabel('f_t(x)')
H=text(-0.2,0.33,'\nu=3');
set(H,'FontSize',14);
H=text(-3.9,0.04,'\nu=5');
set(H,'FontSize',14);
H=text(-2.3,0.02,'\nu=100');
set(H,'FontSize',14);
ylim([0 0.5])
print -deps2 abftpdf.eps





