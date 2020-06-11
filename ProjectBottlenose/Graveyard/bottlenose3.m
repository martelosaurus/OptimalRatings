function bottlenose3
    
    % preliminary calls
    close all;
    set(0,'defaultaxesfontname','consolas');
    set(0,'defaulttextfontname','consolas');
    set(0,'defaultaxesfontweight','bold');
    set(0,'defaulttextfontweight','bold');
    set(0,'defaultaxesfontsize',10);
    set(0,'defaulttextfontsize',10);
    
    % WFA data
    if 0
        x0 = [.05,.25,.48,.65,.79,.89,.95,.97,.99];
    end 
    % Yelp
    if 1
        x0 = [.07,.15,.33,.68];
    end
    
    
    % parameters
    nplot = 1.e3; % #of plotting points
    n = length(x0)+1; % #of messages/actions
    
    % solve
    c0 = [zeros(n+1,1);.5*x0(1)];
    c_star = fsolve(@(c) f(c,x0,n),c0,optimset('MaxFunEvals',1.e5,'MaxIter',1.e4));
    c_star(1:n+1)
    
    % plot
    figure;
    xplot = linspace(0,1,nplot);
    plot(xplot,xplot,'-k'); hold on; % plot straight-talk
    plot(xplot,polyval(c_star(1:n+1),xplot),'--k'); hold on; % plot importance function 
    mplot = interp1([x0,1],(1:n)./n,xplot,'linear','extrap');
    plot(xplot,mplot,'-k','linewidth',2); hold on;
%     for j = 1:n-1
%         temp = linspace(0,(1+j)./n,nplot);
%         plot(x0(j)*ones(1,nplot),temp,'--k','linewidth',2); hold on; 
%     end
    legend('straight-talk','importance function','message function',...
        'location','southoutside');
    axis([0,1,0,1]); axis square;
    xlabel('q');

end

function y = f(c,x,n)

    a = [c(n);zeros(n-1,1)];
    for j = 1:n-1
        a(j+1)=2*x(j)-a(j);
    end
    x = [0,x,1];
    y = zeros(n+2,1);
    for j = 1:n
        y(j)=integral(@(t) (a(j)-t).*polyval(c(1:n+1),t),x(j),x(j+1));
    end
    y(n+1)=polyval(c(1:n+1),1);
    y(n+2)=polyval(c(1:n+1),0)-1;
    %c(1:n)'
    
end