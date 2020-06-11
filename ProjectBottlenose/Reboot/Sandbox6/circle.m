close all

% parameters
nm = 20;
ep = 2*pi/nm; % noise parameter
nep = 1.e3;
n = nep*nm; % #of knots
N = 50; % #of nodes
tol = 1.e-5;

% grid
m = 0:(ep/nep):2*pi; m = m(:);
iu = wshift('1D',1:(n+1),+2*nep); iu = iu(:);
id = wshift('1D',1:(n+1),-2*nep); id = id(:);
L0 = m; L = zeros((n+1),1);

% main loop
er = 2*tol; c =0;
while er>tol
    c = c+1;
    Lu = (L0(iu)+2*pi).*(m+2*ep>2*pi)+L0(iu).*(m+2*ep<=2*pi);
    Lm = L0;
    Ld = (L0(id)-2*pi).*(m-2*ep<0)+L0(id).*(m-2*ep>=0);
    Tu = (.5*Lu.^2+sin(Lu)-Lu.*cos(Lu)-.5*Lm.^2-sin(Lm)+Lm.*cos(Lm))./(Lu-cos(Lu)-Lm+cos(Lm));
    Td = (.5*Lm.^2+sin(Lm)-Lm.*cos(Lm)-.5*Ld.^2-sin(Ld)+Ld.*cos(Ld))./(Lm-cos(Lm)-Ld+cos(Ld));
    L = (Tu+Td)./2; %L(n+1) = 2*pi; %L(1) = 0;
    
    % update
    er = norm(L0-L)
    L0 = L;
    if mod(c,1.e2)==0
    plot(L,m,m,1+sin(m),'linestyle','-','marker','none','linewidth',2);
    axis([0,2*pi,0,2*pi]);
    set(gca,'ytick',0:2*ep:2*pi,'yticklabel','','xtick',[],'xticklabel','');
    grid on; axis square; 
    legend('message','straight talk','location','northwest');
    drawnow;
    end
end