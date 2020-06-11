function BottleDeck

  % preliminary calls
  close all
  set(0,'defaultaxesfontsize',18);
  set(0,'defaulttextfontsize',18);
  
  % global
  global C; 
  C.Nplot = 1.e2; 
  C.cmap = colormap('lines'); 
  
  BottlePlotter('Identity','MessagePlot');
  BottlePlotter('Identity','ActionPlot');
  BottlePlotter('Discrete','MessagePlot');
  BottlePlotter('Discrete','ActionPlot');
  ErrorPlotter('ContinuousError');
  ErrorPlotter('DiscreteError'); 
  ErrorPlotter('BothError'); 
  ObjectiveFunctions(.25);
  
  % loop over Crawford Sobel plots
  CrawfordSobel(.025); 

end

function ObjectiveFunctions(b)

  % load global
  global C; 
  
  cmap = colormap('lines'); 
  figure; 
  xplot = linspace(0,1,C.Nplot);
  plot(xplot,xplot+b,'linewidth',4); hold on;
  plot(xplot,xplot,'linestyle','-','linewidth',4); 
  axis([0,1,0,1]); axis square;
  set(gca,'ticklabelinterpreter','latex');
  set(gca,'xtick',[0,1],'xticklabels',{'$0$','$1$'}); 
  set(gca,'ytick',[0,b,1],'yticklabels',{'$0$','$b$','$1$'}); 
  xlabel('$q$ (state of the world)','interpreter','latex');
  ylabel('$a$ (receiver''s action)','interpreter','latex');
  % annotate
  xlim=get(gca,'XLim');
  ylim=get(gca,'YLim');
  text(0.75*xlim(1)+0.25*xlim(2),0.55*ylim(1)+0.15*ylim(2),...
    'receiver''s preferred action','color',cmap(2,:),'Rotation',45,...
    'FontSize',20,'Interpreter','latex');
  text(0.55*xlim(1)+0.15*xlim(2),0.5*ylim(1)+0.5*ylim(2),...
    'sender''s preferred action','color',cmap(1,:),'Rotation',45,...
    'FontSize',20,'Interpreter','latex');
  % save
  fname = strcat(['ObjectiveFunctions','.eps']);
  saveas(gcf,fname,'epsc'); 
  system(strcat(['epstopdf ',fname])); 
  system(strcat(['del ',fname]));  

end

function CrawfordSobel(b)

  % preliminary calls
  figure;
  global C; 
    
  % breakpoints
  N = ceil(-.5+.5*sqrt(1+2/b)); 
  A = full(spdiags(repmat([1,-2,1],N+1,1),0:2,N-1,N+1));
  A = [1,zeros(1,N);A;zeros(1,N),1]; 
  y = [0;4*b*ones(N-1,1);1]; % right-hand-side
  x = A\y; % breakpoints
  
  for j = 1:N % loop over messages
    xplot = linspace(x(j),x(j+1),C.Nplot); 
    % message
    plot(xplot,ones(1,C.Nplot).*(j/N),'color',C.cmap(1,:),'linewidth',4); hold on;
    % break point
    if N<10
      yplot = linspace(0,1,C.Nplot); 
      plot(x(j)*ones(1,C.Nplot),yplot,'linestyle','--','color',C.cmap(1,:)); hold on; 
    end
  end
  axis([0,1,0,1]); axis square;
  set(gca,'ticklabelinterpreter','latex');
  t = {'$0$','$1/4$','$1/2$','$3/4$','$1$'}; 
  set(gca,'xtick',[0,x(2),x(3),x(4),1],'xticklabels',{'$0$','$x_{1}$','$x_{2}$','$x_{3}$','$1$'}); 
  set(gca,'ytick',0:.25:1,'yticklabels',t); 
  xlabel('$q$ (state of the world)','interpreter','latex');
  ylabel('$\widetilde{m}$ (message)','interpreter','latex');
  fname = 'CrawfordSobel.eps';
  saveas(gcf,fname,'epsc'); 
  system(strcat(['epstopdf ',fname])); 
  system(strcat(['del ',fname]));  
  
end

function BottlePlotter(ModelType,PlotType)

  % global
  global C; 
  figure;

  % parameters
  ep = 1/8;
  m0 = .5;
  n = 1.e3;

  % message function
  q1 = linspace(-ep,1+ep,n);
  q2 = linspace(0,1,n);

  % message region
  switch(ModelType)
    case 'Identity'
      fill([0,0,1,1],[-ep,ep,1+ep,1-ep],.75*ones(1,3)+.25*C.cmap(1,:),'linestyle','none'); hold on;
    case 'Discrete'
      N = round(1/(2*ep))+1; % number of messages 
      % breakpoints
      A = full(spdiags(repmat([1,-2,1],N+1,1),0:2,N-1,N+1));
      A = [1,zeros(1,N);A;zeros(1,N),1]; 
      y = [0;zeros(N-1,1);1]; % right-hand-side
      x = A\y; % breakpoints
      for j = 1:N % loop over messages
        fill([x(j),x(j),x(j+1),x(j+1)],[-ep,ep,ep,-ep]+((j-1)/(N-1)),...
          .75*ones(1,3)+.25*C.cmap(1,:),'linestyle','none'); hold on;
      end
    otherwise
  end

  % custom grid lines
  plot(q1,zeros(size(q1)),'-k'); hold on;
  plot(q1,ones(size(q1)),'-k'); hold on;
  plot(zeros(size(q1)),linspace(-2*ep,1+2*ep,n),'-k'); hold on;
  plot(ones(size(q1)),linspace(-2*ep,1+2*ep,n),'-k'); hold on;
  plot(q1,(1+ep)*ones(size(q1)),'-','color',.6*ones(1,3)); hold on;
  plot(q1,-ep*ones(size(q1)),'-','color',.6*ones(1,3)); hold on;
  plot(q1,(1-ep)*ones(size(q1)),'-','color',.6*ones(1,3)); hold on;
  plot(q1,ep*ones(size(q1)),'-','color',.6*ones(1,3)); hold on;

  switch(ModelType)
    case 'Identity'
      switch(PlotType)
        case 'MessagePlot'
          qu = m0+ep; ql = m0-ep; q0 = q1(q1<qu); n0 = length(q0);
          plot(q0,m0*ones(1,n0),'--k','linewidth',.5); hold on;
          plot(qu*ones(1,n0),linspace(-2*ep,m0,n0),'--k'); hold on;
          plot(ql*ones(1,n0),linspace(-2*ep,m0,n0),'--k'); hold on;
          plot(q2,q2,'color',C.cmap(1,:),'linewidth',4); hold on;
          set(gca,'xtick',[0,ql,qu,1],'xticklabels',{'$0$','$q_{-}(\widetilde{m}_{0})$','$q_{+}(\widetilde{m}_{0})$','$1$'});
        case 'ActionPlot'
          qu = m0+ep; ql = m0-ep; q0 = q1(q1<.5*(qu+ql)); n0 = length(q0);
          plot(q0,m0*ones(1,n0),'--k','linewidth',.5); hold on;
          yl = linspace(0,ep,C.Nplot); ml = 2*yl-ep; 
          yr = linspace(1-ep,1,C.Nplot); mr = 2*(yr-1)+1+ep; 
          ym = linspace(ep,1-ep,C.Nplot); mm = ym;
          y = [yl,ym,yr]; m = [ml,mm,mr];
          plot(.5*(ql+qu)*ones(1,n0),linspace(-2*ep,m0,n0),'--k'); hold on;
          plot(y,m,'color',C.cmap(5,:),'linewidth',4); hold on;
          set(gca,'xtick',[0,(ql+qu)/2,1],'xticklabels',{'$0$','$a^{*}(\widetilde{m}_{0})$','$1$'});
        otherwise
      end
    case 'Discrete'
      switch(PlotType)
        case 'MessagePlot'
          for j = 1:N
            xplot = linspace(x(j),x(j+1),C.Nplot); 
            plot(xplot,ones(1,C.Nplot).*((j-1)/(N-1)),'color',C.cmap(1,:),'linewidth',4); hold on;
          end
          % THIS ISN'T RIGHT
          qu = m0+.5/N; ql = m0-.5/N; q0 = q1(q1<qu); n0 = length(q0);
          m0 = m0+.5*ep;
          plot(q0,m0*ones(1,n0),'--k','linewidth',m0); hold on;
          plot(qu*ones(1,n0),linspace(-2*ep,m0,n0),'--k'); hold on;
          plot(ql*ones(1,n0),linspace(-2*ep,m0,n0),'--k'); hold on;
          set(gca,'xtick',[0,ql,qu,1],'xticklabels',{'$0$','$q_{-}(\widetilde{m}_{0})$','$q_{+}(\widetilde{m}_{0})$','$1$'}); 
        case 'ActionPlot'
          qu = m0+.5/N; ql = m0-.5/N; q0 = q1(q1<.5*(qu+ql)); n0 = length(q0);
          m0 = m0+.5*ep;
          for j = 1:N % loop over messages
            fill([x(j),x(j),x(j+1),x(j+1)],[-ep,ep,ep,-ep]+((j-1)/(N-1)),...
              .75*ones(1,3)+.25*C.cmap(1,:),'linestyle','none'); hold on;
            plot(.5*(x(j)+x(j+1))*ones(1,C.Nplot),linspace(-ep,ep,C.Nplot)+((j-1)/(N-1)),...
              'color',C.cmap(5,:),'linewidth',4); hold on;
          end
          plot(.5*(ql+qu)*ones(1,n0),linspace(-2*ep,m0,n0),'--k'); hold on;
          plot(q0,m0*ones(1,n0),'--k','linewidth',.5); hold on;
          set(gca,'xtick',[0,(ql+qu)/2,1],'xticklabels',{'$0$','$a^{*}(\widetilde{m}_{0})$','$1$'});
        otherwise
      end
    otherwise
  end
  axis([-ep,1+ep,-2*ep,1+2*ep]);
  xlabel('$q$ (state of the world)','interpreter','latex');
  ylabel('$\widetilde{m}$ (received message)','interpreter','latex');
  set(gca,'ticklabelinterpreter','latex',...
      'ytick',[-ep,0,ep,m0,1-ep,1,1+ep],'yticklabels',{'$-\bar{\epsilon}$','$0$','$\bar{\epsilon}$','$\widetilde{m}_{0}$','$1-\bar{\epsilon}$','$1$','$1+\bar{\epsilon}$'});
  axis square;
  
  % save
  fname = strcat([ModelType,PlotType,'.eps']); 
  saveas(gcf,fname,'epsc'); 
  system(strcat(['epstopdf ',fname])); 
  system(strcat(['del ',fname]));  

end

function ErrorPlotter(model)

  % global
  global C;
  figure; 

  % parameters
  ep = .1;
  n = 1.e3;

  % message function
  q = linspace(-ep,1+ep,n);
  c = .75*ones(1,3)+.25*C.cmap(2,:);
  xcplot = [-ep,ep,1-ep,1+ep,1-ep,ep,-ep];
  ycplot = [0,ep,ep,0,-ep,-ep,0];
  xdplot = [-ep,-ep,1+ep,1+ep,-ep];
  ydplot = ep*[-1,1,1,-1,-1]./(1+2*ep);
  switch(model)
    case 'ContinuousError'
      fill(xcplot,ycplot,c); hold on;
      plot(xcplot,ycplot,'color',C.cmap(2,:),'linewidth',4); hold on;
    case 'DiscreteError'
      fill(xdplot,ydplot,c); hold on;
      plot(xdplot,ydplot,'color',C.cmap(2,:),'linewidth',4); hold on;
    otherwise
      fill(xdplot,ydplot,c); hold on;
      plot(xdplot,ydplot,'color',C.cmap(2,:),'linewidth',4); hold on;    
      plot(xcplot,ycplot,'color',C.cmap(2,:),'linewidth',4); hold on;    
  end
  
  % grid
  plot(q,zeros(size(q)),'-k'); hold on;
  plot(ep*ones(size(q)),linspace(-2*ep,2*ep,n),'-k'); hold on;
  plot((1-ep)*ones(size(q)),linspace(-2*ep,2*ep,n),'-k'); hold on;
  plot(q,ep*ones(size(q)),'-','color',.6*ones(1,3)); hold on;
  plot(q,-ep*ones(size(q)),'-','color',.6*ones(1,3)); hold on;
  
  axis([-ep,1+ep,-1.5*ep,1.5*ep]); axis square; 
  xlabel('$\widetilde{m}$ (received message)','interpreter','latex');
  ylabel('$a(\widetilde{m})-q$ (error)','interpreter','latex');
  set(gca,'ticklabelinterpreter','latex',...
      'xtick',[-ep,ep,1-ep,1+ep],'xticklabels',{'$-\bar{\epsilon}$','$0$','$1-\bar{\epsilon}$','$1+\bar{\epsilon}$','$1$','$1+\bar{\epsilon}$'},...
      'ytick',[-ep,0,ep],'yticklabels',{'$-\bar{\epsilon}$','$0$','$\bar{\epsilon}$',});
  axis square;
  
  % save
  %fname = strcat([model,num2str(k),'.eps']);
  fname = strcat([model,'.eps']);
  saveas(gcf,fname,'epsc'); 
  system(strcat(['epstopdf ',fname])); 
  system(strcat(['del ',fname]));  

end