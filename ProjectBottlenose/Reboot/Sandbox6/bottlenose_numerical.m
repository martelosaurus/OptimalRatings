% preliminary calls
close all
set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);

% parameters
ep = .1;
n = 8;
N = 1.e3;
tol = 1.e-6;

% grid
q = linspace(0,1,N); 

% solution
[a,e] = bottlenose_driver(n);
disp(strcat(['error = ',num2str(e)]));
m = zeros(N,n-1);
for i = 1:n-1
    m(:,i) = q(:).^i;
end
m = m*a;

% plot
plot(q,zeros(size(q)),'-k'); hold on;
plot(q,ones(size(q)),'-k'); hold on;
plot(zeros(size(q)),linspace(-2*ep,1+2*ep,N),'-k'); hold on;
plot(ones(size(q)),linspace(-2*ep,1+2*ep,N),'-k'); hold on;
plot(q,(1+ep)*ones(size(q)),'-k'); hold on;
plot(q,-ep*ones(size(q)),'-k'); hold on;

plot(q,m,'-b','linewidth',2); hold on;
plot(q,m+ep,'--b'); hold on;
plot(q,m-ep,'--b'); hold on;
plot(q,exp(-10*(q-.5).^2),'-r'); hold on;

axis([-ep,1+ep,-2*ep,1+2*ep]);
axis square; 