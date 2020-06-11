function bottlenose_jump

    % preliminary calls
    close all;
    set(0,'defaultaxesfontsize',12);
    set(0,'defaulttextfontsize',12);

    global X;

    % parameters
    X.ep = .1; % noise parameter
    X.n = 100; % #of knots
    X.N = 100; % #of nodes
    X.tol = 1.e-8; 

    % importance function
    I = @(r) r.*(1-r);
    X.I0 = @(s,t) integral(@(r) I(r),s,t);
    X.I1 = @(s,t) integral(@(r) r.*I(r),s,t);

    % knots
    X.t = linspace(0,1,X.n); % knots
    
    % solve for inverse of message
    m0 = X.t(:); er = 2*X.tol;
    while er>X.tol
        m1 = F(m0);
        er = norm(m1-m0);
        m0 = m1;
    end
    
    % grid
    q = linspace(0,1,X.N);
    m = dsearchn(m1,q(:));
        
    % plot
    figure;
    plot(q,m,'-k','linewidth',2); hold on;
    plot(q,I(q),'--k'); hold on;
    
    A1 = area(q,m+X.ep); hold on;
    A1(1).FaceColor = .75*ones(1,3);
    A1(1).LineStyle = 'none';
    A2 = area(q,m-X.ep,-X.ep); hold on;
    A2(1).FaceColor = 'w';
    A2(1).LineStyle = 'none';
    
    plot(q,zeros(size(q)),'-k'); hold on;
    plot(q,ones(size(q)),'-k'); hold on;
    plot(zeros(size(q)),linspace(-2*X.ep,1+2*X.ep,X.N),'-k'); hold on;
    plot(ones(size(q)),linspace(-2*X.ep,1+2*X.ep,X.N),'-k'); hold on;
    plot(q,(1+X.ep)*ones(size(q)),'-k'); hold on;
    plot(q,-X.ep*ones(size(q)),'-k'); hold on;

    plot(q,m,'-k','linewidth',2); hold on;
    plot(q,I(q),'--k'); hold on;
    
    legend('message','importance','location','northwest');

    axis([-X.ep,1+X.ep,-2*X.ep,1+2*X.ep]);
    xlabel('$q$','interpreter','latex');
    ylabel('$m$','interpreter','latex');
    set(gca,'ticklabelinterpreter','latex','xtick',[0,1],'ytick',[-X.ep,0,1,1+X.ep],'yticklabels',{'$-\bar{\epsilon}$','$0$','$1$','$1+\bar{\epsilon}$'});
    axis square; 
    
end

function m1 = F(m0)
    global X;
    m1 = zeros(X.n,1);
    % Euler-Lagrange
    for i = 2:X.n-1
        qu = min([m0(dsearchn(X.t(:),X.t(i)+2*X.ep)),1]);
        qm = m0(i);
        ql = max([m0(dsearchn(X.t(:),X.t(i)-2*X.ep)),0]);
        m1(i) = .5*(X.I1(qm,qu)./X.I0(qm,qu)+X.I1(ql,qm)./X.I0(ql,qm));
    end
    % endpoints
    m1(1) = 0; m1(X.n) = 1;
end