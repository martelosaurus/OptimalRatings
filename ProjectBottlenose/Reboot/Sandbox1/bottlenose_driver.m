function bottlenose_driver

    % preliminary calls
    close all;
    set(0,'defaultaxesfontsize',12);
    set(0,'defaulttextfontsize',12);

    global X;

    % parameters
    X.ep = .1; % noise parameter
    X.n = 4; % #of knots
    X.N = 50; % #of nodes
    X.tol = 1.e-10; 

    % importance function
    I = @(r) 2*r.*(1-r);
    X.I0 = @(s,t) integral(@(r) I(r),s,t,'AbsTol',X.tol,'RelTol',X.tol);
    X.I1 = @(s,t) integral(@(r) r.*I(r),s,t,'AbsTol',X.tol,'RelTol',X.tol);

    % knots
    X.Q = (cos(pi*((X.n-1):-1:0)./(X.n-1))+1)./2; % chebychev knots
    
    % solve
    a0 = [0;1;zeros(X.n-2,1)];
    [a,e] = fsolve(@(s) F(s),a0,optimset('MaxFunEvals',1.e4));
    disp(norm(e));
    
    % grid
    q = linspace(0,1,X.N); 
    m = M(q,a);
    
    % plot
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

function y = F(a)
    global X;
    y = zeros(X.n-1,1);
    % Euler-Lagrange
    for i = 2:X.n-1
        q = X.Q(i);
        qu = min([fzero(@(s)M(s,a)-M(q,a)-2*X.ep,q,optimset('TolX',X.tol)),1]);
        qd = max([fzero(@(s)M(s,a)-M(q,a)+2*X.ep,q,optimset('TolX',X.tol)),0]);
        y(i) = X.I1(qu,q)./X.I0(qu,q)+X.I1(q,qd)./X.I0(q,qd)-2*q;
    end
    % endpoints
    y(1) = M(0,a);
    y(X.n) = 1-M(1,a);
end

function y = M(q,a)
    global X;
    z = zeros(length(q),X.n);
    z(:,1) = ones(length(q),1); z(:,2) = q(:);
    for j = 3:X.n
        z(:,j)=2*q(:).*z(:,j-1)-z(:,j-2);
    end
    y = z*a;
end