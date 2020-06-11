function bottlenose

    % preliminary calls
    close all;
    set(0,'defaultaxesfontsize',12);
    set(0,'defaulttextfontsize',12);

    global X;

    % parameters
    X.ep = .1; % noise parameter
    X.n = 20; % #of knots
    X.N = 50; % #of nodes
    X.tol = 1.e-8; 

    % importance function
    I = @(r) 2*r.*(1-r);
    X.I0 = @(s,t) integral(@(r) I(r),s,t);
    X.I1 = @(s,t) integral(@(r) r.*I(r),s,t);

    % knots
    X.t = (cos(pi*((X.n-1):-1:0)./(X.n-1))+1)./2; % chebychev knots
    
    % test
    figure;
    t_plot = linspace(0,1,X.N);
    demo = @(t) 1.5*(atan(2*tan(pi/3)*(t-.5))+pi/2)./pi-.25;
    demo_plot = Mp(t_plot,Mc(X.t,demo(X.t)));
    plot(t_plot,demo_plot,'-r',t_plot,demo(t_plot),'-b',X.t,demo(X.t),'ob');
    
    % solve for inverse of message
    a0 = [0;1;zeros(X.n-2,1)]; er = 2*X.tol;
    while er>X.tol
        er
        a1 = F(a0);
        er = norm(a1-a0);
        a0 = a1;
    end
    
    % grid
    q = linspace(0,1,X.N);
    m = zeros(X.N,1);
    for i = 1:X.N
        m(i) = fzero(@(s) q(i)-Mp(s,a1),q(i));
    end
        
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

function a1 = F(a0)
    global X;
    y = zeros(X.n,1);
    % Euler-Lagrange
    for i = 2:X.n-1
        qu = min([Mp(X.t(i)+2*X.ep,a0),1]);
        qm = Mp(X.t(i),a0);
        ql = max([Mp(X.t(i)-2*X.ep,a0),0]);
        y(i) = .5*(X.I1(qm,qu)./X.I0(qm,qu)+X.I1(ql,qm)./X.I0(ql,qm));
    end
    % endpoints
    y(1) = 0; y(X.n) = 1;
    a1 = Mc(X.t,y);
end

% chebychev polynomial
function y = Mp(q,a)
    global X;
    z = zeros(length(q),X.n);
    z(:,1) = ones(length(q),1); z(:,2) = q(:);
    for j = 3:X.n
        z(:,j)=2*q(:).*z(:,j-1)-z(:,j-2);
    end
    y = z*a;
end

% chebychev coefficients
function a = Mc(q,y)
    global X;
    z = zeros(length(q),X.n);
    z(:,1) = ones(length(q),1); z(:,2) = q(:);
    for j = 3:X.n
        z(:,j)=2*q(:).*z(:,j-1)-z(:,j-2);
    end
    a = z\y(:);
end