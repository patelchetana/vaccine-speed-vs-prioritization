function figure_5(nasem,random,optimal,const_mit,x_lim,y_lim)
figure;
hold on;
plot(nasem.theta(1,:).^2,'r-','LineWidth',1);
plot(random.theta(1,:).^2,'b-','LineWidth',1);
plot(optimal.theta(1,:).^2,'k-','LineWidth',1);
plot(repelem(const_mit.^2,1,length(optimal.theta(1,:))),'k--','LineWidth',1);
hold off;
xlim(x_lim);
ylim(y_lim);
lgd = legend('NASEM','None','Optimal','Sustained');
lgd.Title.String = 'Prioritization Policy';
xlabel('Days');
ylabel('Community Mitigation Parameter');
end