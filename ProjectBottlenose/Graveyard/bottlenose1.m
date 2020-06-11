function [a,x]=bottlenose1

    close all;
    set(0,'defaultaxesfontname','consolas');
    set(0,'defaulttextfontname','consolas');
    set(0,'defaultaxesfontweight','bold');
    set(0,'defaulttextfontweight','bold');
    set(0,'defaultaxesfontsize',10);
    set(0,'defaulttextfontsize',10);

    % global 
    global p;

    % parameters
    p.nplot = 100; % #of plotting points
    p.L = @(t) .5*t.^2; % loss function
    p.l = @(t) t; % derivative of loss function
    p.g = @(t) ones(size(t)); % state-of-the-world density

    %% rude #1
    p.n = 2; % #of messages/actions
    
    % solve
    p.I = @(t) .1+.8*t; % importance function
    prob.objective = @(s) f(s);
    prob.nonlcon = @(s) g(s);
    prob.x0 = [p.I((1:p.n-1)./p.n),1,(1:p.n-1)./p.n];
    prob.solver = 'fmincon';
    prob.options = optimset('display','off');
    s = fmincon(prob);
    x = s(p.n+1:2*p.n-1)
    
    % rude
    figure;
    xplot = linspace(0,1,p.nplot);
    plot(xplot,xplot,'-k'); hold on; % plot straight-talk
    plot(xplot,p.I(xplot),'--k'); hold on; % plot importance function 
    mplot = interp1([0,x,1],(0:p.n)./p.n,xplot,'linear','extrap');
    plot(xplot,mplot,'-k','linewidth',2); hold on;
%     mplot = interp1([x,1],(1:p.n)./p.n,xplot,'next','extrap');
%     plot(xplot,mplot,'-k','linewidth',2); hold on;
%     for j = 1:p.n-1
%         temp = linspace(0,(1+j)./p.n,p.nplot);
%         plot(x(j)*ones(1,p.nplot),temp,'--k','linewidth',2); hold on; 
%     end
    csvwrite('rude2.csv',[xplot',p.I(xplot)',mplot']);
    legend('straight-talk','importance function','message function',...
        'location','northwest');
    axis([0,1,0,1]); axis square;
    xlabel('q');
    title('rude');
    
    %% polite #1
    p.n = 5; % #of messages/actions
    
    % solve
    p.I = @(t) .9-.8*t; % importance function
    prob.objective = @(s) f(s);
    prob.nonlcon = @(s) g(s);
    prob.x0 = [(1:p.n-1)./p.n,1,(1:p.n-1)./p.n];
    prob.solver = 'fmincon';
    prob.options = optimset('display','off');
    s = fmincon(prob);
    x = s(p.n+1:2*p.n-1)
    
    % polite
    figure;
    xplot = linspace(0,1,p.nplot);
    plot(xplot,xplot,'-k'); hold on; % plot straight-talk
    plot(xplot,p.I(xplot),'--k'); hold on; % plot importance function 
    %mplot = interp1([0,x,1],(0:p.n)./p.n,xplot,'linear','extrap');
    %plot(xplot,mplot,'-k','linewidth',2); hold on;
    mplot = interp1([x,1],(1:p.n)./p.n,xplot,'next','extrap');
    plot(xplot,mplot,'-k','linewidth',2); hold on;
    for j = 1:p.n-1
        temp = linspace(0,(1+j)./p.n,p.nplot);
        plot(x(j)*ones(1,p.nplot),temp,'--k','linewidth',2); hold on; 
    end
    csvwrite('polite2.csv',[xplot',p.I(xplot)',mplot']);
    legend('straight-talk','importance function','message function',...
        'location','northwest');
    axis([0,1,0,1]); axis square;
    xlabel('q');
    title('polite');
    
    %% rude #2
    p.n = 3; % #of messages/actions
    
    % solve
    p.I = @(t) .1+.8*t; % importance function
    prob.objective = @(s) f(s);
    prob.nonlcon = @(s) g(s);
    prob.x0 = [p.I((1:p.n-1)./p.n),1,(1:p.n-1)./p.n];
    prob.solver = 'fmincon';
    prob.options = optimset('display','off');
    s = fmincon(prob);
    x = s(p.n+1:2*p.n-1)
    
    % rude
    figure;
    xplot = linspace(0,1,p.nplot);
    plot(xplot,xplot,'-k'); hold on; % plot straight-talk
    plot(xplot,p.I(xplot),'--k'); hold on; % plot importance function 
    mplot = interp1([0,x,1],(0:p.n)./p.n,xplot,'linear','extrap');
    plot(xplot,mplot,'-k','linewidth',2); hold on;
%     mplot = interp1([x,1],(1:p.n)./p.n,xplot,'next','extrap');
%     plot(xplot,mplot,'-k','linewidth',2); hold on;
%     for j = 1:p.n-1
%         temp = linspace(0,(1+j)./p.n,p.nplot);
%         plot(x(j)*ones(1,p.nplot),temp,'--k','linewidth',2); hold on; 
%     end
    csvwrite('rude3.csv',[xplot',p.I(xplot)',mplot']);
    legend('straight-talk','importance function','message function',...
        'location','northwest');
    axis([0,1,0,1]); axis square;
    xlabel('q');
    title('rude');
    
    %% polite #2
    p.n = 3; % #of messages/actions
    
    % solve
    p.I = @(t) .9-.8*t; % importance function
    prob.objective = @(s) f(s);
    prob.nonlcon = @(s) g(s);
    prob.x0 = [p.I((1:p.n-1)./p.n),1,(1:p.n-1)./p.n];
    prob.solver = 'fmincon';
    prob.options = optimset('display','off');
    s = fmincon(prob);
    x = s(p.n+1:2*p.n-1)
    
    % polite
    figure;
    xplot = linspace(0,1,p.nplot);
    plot(xplot,xplot,'-k'); hold on; % plot straight-talk
    plot(xplot,p.I(xplot),'--k'); hold on; % plot importance function 
    mplot = interp1([0,x,1],(0:p.n)./p.n,xplot,'linear','extrap');
    plot(xplot,mplot,'-k','linewidth',2); hold on;
%     mplot = interp1([x,1],(1:p.n)./p.n,xplot,'next','extrap');
%     plot(xplot,mplot,'-k','linewidth',2); hold on;
%     for j = 1:p.n-1
%         temp = linspace(0,(1+j)./p.n,p.nplot);
%         plot(x(j)*ones(1,p.nplot),temp,'--k','linewidth',2); hold on; 
%     end
    csvwrite('polite3.csv',[xplot',p.I(xplot)',mplot']);
    legend('straight-talk','importance function','message function',...
        'location','northwest');
    axis([0,1,0,1]); axis square;
    xlabel('q');
    title('polite');
    
end

function v = f(s)

    % load the global
    global p;

    % actions and breakpoints
    a = s(1:p.n); 
    x = [0,s(p.n+1:2*p.n-1),1];
    
    % compute value
    v = 0;
    for j = 1:p.n
        v = v+integral(@(t)p.L(a(j)-t).*p.I(t).*p.g(t),x(j),x(j+1));
    end

end

function [c,ceq] = g(s)
    
    % load the global
    global p;

    % actions and breakpoints
    a = s(1:p.n); 
    x = [0,s(p.n+1:2*p.n-1),1];
    ceq = zeros(1,p.n);
   
    c = [];
    for j = 1:p.n
        ceq(j)=integral(@(t) p.l(a(j)-t).*p.I(t).*p.g(t),x(j),x(j+1));
    end

end