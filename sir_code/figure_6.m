function figure_6(selection_range,N_mitigation,...
    nasem,random,optimal,baseline,...
    int,mitigation_labels)
x.nasem = (selection_range.nasem-1)*int;
x.random = (selection_range.random-1)*int;
x.optimal = (selection_range.optimal-1)*int;
figure;
t=tiledlayout(2,N_mitigation,'TileSpacing','Compact');

% Plot years of life lost
for i=1:N_mitigation
    ax = nexttile;
    hold on;
    grid on;
    plot(x.nasem,100*(baseline.YLL-nasem.YLL(selection_range.nasem,i))/baseline.YLL,'r-','LineWidth',1);
    plot(x.random,100*(baseline.YLL-random.YLL(selection_range.random,i))/baseline.YLL,'b-','LineWidth',1);
    plot(x.optimal,100*(baseline.YLL-optimal.YLL(selection_range.optimal,i))/baseline.YLL,'k-','LineWidth',1);
    hold off;
    ylim([0 100]);
    title(mitigation_labels(i));
    set(gca,'xticklabel',[]);
    if i==1
        ylabel('Reduction in YLL (%)');
    else
        set(gca,'yticklabel',[]);
    end
end

% Plot infections
for i=1:N_mitigation
    ax = nexttile;
    hold on;
    grid on;
    plot(x.nasem,100*(baseline.cases-nasem.cases(selection_range.nasem,i))/baseline.cases,'r-','LineWidth',1);
    plot(x.random,100*(baseline.cases-random.cases(selection_range.random,i))/baseline.cases,'b-','LineWidth',1);
    plot(x.optimal,100*(baseline.cases-optimal.cases(selection_range.optimal,i))/baseline.cases,'k-','LineWidth',1);
    hold off;
    ylim([0 100]);
    if i==1
        ylabel('Reduction in Total Cases (%)');
    else
        set(gca,'yticklabel',[]);   
    end
end


lgd = legend(ax,'NASEM','None','Optimal');
lgd.Title.String = 'Prioritization Policy';
lgd.Layout.Tile = 'east';
xlabel(t,'Millions Vaccinated per Month');
end

