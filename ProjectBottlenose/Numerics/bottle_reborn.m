% preliminary calls
close all;
set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);

% parameters
m = 20; % #of messages
ep = 1/(2*m); % noise parameter
nep = 100;
n = nep/ep; % #of knots
N = 50; % #of nodes
tol = 1.e-5; 

% importance function
I = @(r) r.^2;
I0 = @(s,t) 1; 
I1 = @(s,t) (s+t)./2;

% knots and initial guess
m = (1/n):(1/n):1; m = m(:);
iuu = knnsearch(m,min([m+2*ep,ones(n,1)],[],2));
idd = knnsearch(m,max([m-2*ep,zeros(n,1)],[],2));
iu = knnsearch(m,min([m+ep,ones(n,1)],[],2));
id = knnsearch(m,max([m-ep,zeros(n,1)],[],2));
% iuu = wshift('1D',1:n,2*nep);
% idd = wshift('1D',1:n,2*nep);
% iu = wshift('1D',1:n,nep);
% id = wshift('1D',1:n,nep);
L0 = m.^2; L = zeros(n,1); A = zeros(n,1);
    
% main loop
R = sqrt(2)*[1,1;-1,1]./2;
er = 2*tol; c = 0;
while er>tol
    c = c+1;
    L = .5*I1(L0,L0(iuu))./I0(L0,L0(iuu))+.5*I1(L0(idd),L0)./I0(L0(idd),L0);
    L(1) = 0; L(n) = 1;
    A = I1(L0(id),L0(iu))./I0(L0(id),L0(iu));
    % update
    er = norm(L0-L)
    L0 = L;
    if mod(1.e3,c)==0
        plot(L,m,m,m,m,I(m),m,A,'linewidth',2); grid on; axis square; 
        legend('message','straight talk','importance','action','location','northwest');
        drawnow;
    end
end