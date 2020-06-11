% preliminary calls
close all
set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);

% parameters
ep = .1;
al = 2*tan(5*pi/12);
n = 1.e3;

% message function
q = linspace(-ep,1+ep,n);
qm = linspace(0,1,n);
M = @(t) 1.2*(atan(al*(t-.5))+pi/2)./pi-.1;
L = @(t) fzero(@(s) t-M(s),t);
m = M(q); mm = M(qm);
mu = m+ep; mum = mm+ep;
md = m-ep; mdm = mm-ep;

% grid
%A1 = area(qm,min([mum;ones(size(mum))],[],1)); hold on;
A1 = area(qm,mum); hold on;
A1(1).FaceColor = .8*ones(1,3);
A1(1).LineStyle = 'none';
A2 = area(qm,mdm,-ep); hold on;
A2(1).FaceColor = 'w';
A2(1).LineStyle = 'none';

plot(q,zeros(size(q)),'-k'); hold on;
plot(q,ones(size(q)),'-k'); hold on;
plot(zeros(size(q)),linspace(-2*ep,1+2*ep,n),'-k'); hold on;
plot(ones(size(q)),linspace(-2*ep,1+2*ep,n),'-k'); hold on;
plot(q,(1+ep)*ones(size(q)),'-','color',.6*ones(1,3)); hold on;
plot(q,-ep*ones(size(q)),'-','color',.6*ones(1,3)); hold on;
plot(q,(1-ep)*ones(size(q)),'-','color',.6*ones(1,3)); hold on;
plot(q,ep*ones(size(q)),'-','color',.6*ones(1,3)); hold on;
qu = L(.8+ep); ql = L(.8-ep); q0 = q(q<qu); n0 = length(q0);
plot(q0,.8*ones(1,n0),'--k','linewidth',.5); hold on;
plot(qu*ones(1,n0),linspace(-2*ep,.8,n0),'--k'); hold on;
plot(ql*ones(1,n0),linspace(-2*ep,.8,n0),'--k'); hold on;

plot(q,m,'-k','linewidth',2); hold on;
plot(q,mu,'-k'); hold on;
plot(q,md,'-k'); hold on;

axis([-ep,1+ep,-2*ep,1+2*ep]);
xlabel('$q$','interpreter','latex');
ylabel('$\tilde{m}$','interpreter','latex');
set(gca,'ticklabelinterpreter','latex',...
    'xtick',[0,ql,qu,1],'xticklabels',{'$0$','$\underline{q}(\tilde{m}_{0})$','$\overline{q}(\tilde{m}_{0})$','$1$'},...
    'ytick',[-ep,0,ep,.8,1-ep,1,1+ep],'yticklabels',{'$-\bar{\epsilon}$','$0$','$\bar{\epsilon}$','$\tilde{m}_{0}$','$1-\bar{\epsilon}$','$1$','$1+\bar{\epsilon}$'});
axis square; 