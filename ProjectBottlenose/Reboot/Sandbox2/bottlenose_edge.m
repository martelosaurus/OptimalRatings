% preliminary calls
close all
set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);

% parameters
n = 1.e2;

% grid
q = (1/(n+1)):(1/(n+1)):(1-1/(n+1));
m = tan(pi*(2*q-1)./2);
plot(q,m,'-k','linewidth',2); hold on;
h = .05;
for i=1:100
    plot(q,m+i*h,'-','color',ones(1,3)*(i/100),'linewidth',2); hold on;
    plot(q,m-i*h,'-','color',ones(1,3)*(i/100),'linewidth',2); hold on;
end
plot(q,m,'-w','linewidth',2); hold on;
axis([0,1,-10,10]); grid on;