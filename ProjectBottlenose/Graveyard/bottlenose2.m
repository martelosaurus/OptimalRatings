function [b,sol]=bottlenose2

    close all;
    set(0,'defaultaxesfontname','consolas');
    set(0,'defaulttextfontname','consolas');
    set(0,'defaultaxesfontweight','bold');
    set(0,'defaulttextfontweight','bold');
    set(0,'defaultaxesfontsize',10);
    set(0,'defaulttextfontsize',10);

    n = 101;
    b = linspace(0,1,n);
    sol = zeros(1,n);
    for i = 1:n
        temp=fsolve(@(t)f(t,b(i)),[.25,.75,.5]);
        sol(i)=temp(3);
    end
    
    plot(b,sol,'-k','linewidth',2);
    axis([0,1,.25,.75]); 
    ax = gca;
    ax.XTick = [0,.25,.5,.75,1];
    ax.YTick = [.25,.375,.5,.675,.75];
    ax.XTickLabel = {'0','rude','1/2','polite','1'};
    ax.YTickLabel = {'1/4','big H bucket','1/2','big L bucket','3/4'};
    axis square; grid on;
    xlabel('b');
    ylabel('x');
    
end

function y = f(x,b)

    y(1)=integral(@(t)(x(1)-t).*(b+(1-2*b).*t),0,x(3));
    y(2)=integral(@(t)(x(2)-t).*(b+(1-2*b).*t),x(3),1);
    y(3)=x(1)+x(2)-2*x(3); 

end