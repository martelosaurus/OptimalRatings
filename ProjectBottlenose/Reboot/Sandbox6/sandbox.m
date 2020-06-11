function z = sandbox(c)
    n = 10;
    N = 10;
    x = [];
    y = [];
    y1 = 1;
    for i = 0:N-1
        x0 = linspace(.1^(i+1),.1^i,n);
        x = [x0,x];
        if mod(i,2)
            y = [c*x0+(y1-c*.1^i),y];
        else
            y = [y1*ones(1,n),y];
        end
        y1 = y(1);
    end
    loglog(x,y,'linewidth',2);
    axis([-.5,1.5,-.5,1.5]);
    grid on;
    drawnow;
    z = y(1);
end