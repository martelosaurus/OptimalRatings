function kiwi

    % globals
    global H; 

    % parameters
    Np = 1.e3;
    IN = [1:10,20:10:50,100];
    N = max(IN); 

    % handles
    H.I = @(t) t.^2;
    H.T = @(t1,t2) integral(@(t)t.*I(t),t1,t2)./integral(I,t1,t2);

    % main loop
    xsol0 = .5;
    for n = 2:N
        xsol1 = fsolve(@F,xsol0);
        xsol0 = sort([xsol1,.5]);
    end
    
    % midpoint estimation
    xp = linspace(0,1,Np);
    cs = spline([xsol1,1],(0:N)./N,xp); 
    
    % plot
    for j = 1:length(IN)
        n = IN(j); 
        
        if any(n==IN)
            saveas(gcf,strcat(['sol',num2str(n),'.eps']),'epsc'); 
        end
    end
end

function y = F(x,n)
    global H;
    y = zeros(1,n-1);
    y(1) = 2*x(1)-H.T(x(1),x(2))-H.T(0,x(1));
    for k = 1:n-2, 
        y(k) = 2*x(k+1)-H.T(x(k+1),x(k+2))-H.T(x(k),x(k+1)); 
    end
    y(n-1) = 2*x(n-1)-H.T(x(n-1),1)-H.T(x(n-2),x(n-1));
end