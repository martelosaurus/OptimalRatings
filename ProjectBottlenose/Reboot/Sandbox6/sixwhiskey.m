close all

% parameters
nm = 20;
ep = 1/nm; % noise parameter
nep = 1.e3;
n = nep*nm+1; % #of knots
N = 50; % #of nodes
tol = 1.e-5;

% grid
m = 0:(ep/nep):1; m = m(:);
iu = [(nep+1):n,n*ones(1,nep)]; iu = iu(:);
id = [ones(1,nep),1:(n-nep)]; id = id(:);
L0 = m.^3; L = zeros(n,1);

% main loop
er = 2*tol; c = 0;
while er>tol
    c = c+1;
    L = 2*((L0(iu).^2+L0(iu).*L0(id)+L0(id).^2)./(L0(iu)+L0(id)))./3;

    % boundaries
    L(1) = 0; L(n) = 1;
    
    % plot
    er = norm(L0-L); 
    
    % update
    L0 = L;
    
    if mod(c,1.e1)==0
        disp(er);
        plot(m,L,'linestyle','-','marker','none','linewidth',2);
        axis([0,1,0,1]);
        set(gca,'ytick',0:2*ep:1,'yticklabel','','xtick',[],'xticklabel','');
        grid on; axis square; 
        legend('message','location','northwest');
        drawnow;
    end
end