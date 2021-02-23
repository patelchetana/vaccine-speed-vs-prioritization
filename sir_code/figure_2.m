function figure_2(categories,share,...
    nasem1,random1,optimal1,...
    nasem2,random2,optimal2,...
    x_lim,y_lims,mitigation_labels,category_labels)

nasem1_cases = nasem1.I+nasem1.Ix+nasem1.Iv+nasem1.E+nasem1.Ex+nasem1.Ev;
random1_cases = random1.I+random1.Ix+random1.Iv+random1.E+random1.Ex+random1.Ev;
optimal1_cases = optimal1.I+optimal1.Ix+optimal1.Iv+optimal1.E+optimal1.Ex+optimal1.Ev;

nasem2_cases = nasem2.I+nasem2.Ix+nasem2.Iv+nasem2.E+nasem2.Ex+nasem2.Ev;
random2_cases = random2.I+random2.Ix+random2.Iv+random2.E+random2.Ex+random2.Ev;
optimal2_cases = optimal2.I+optimal2.Ix+optimal2.Iv+optimal2.E+optimal2.Ex+optimal2.Ev;

figure;
N_categories = length(categories);
t = tiledlayout(2,2*N_categories,'TileSpacing','Compact');
for i=1:N_categories 
    range = categories(i,1):categories(i,2); % Range to subset
    agg_mult = share(range)/sum(share(range)); % Multiplier to get weighted average
    
    % Policy 1 Infections
    nexttile(2*i-1);
    hold on;
    grid on;
    plot(100*(nasem1_cases(range,:))'*agg_mult,'r-','LineWidth',1);
    xline(nasem1.TTHI,'r:','LineWidth',1);
    plot(100*(random1_cases(range,:))'*agg_mult,'b-','LineWidth',1);
    xline(random1.TTHI,'b:','LineWidth',1);
    plot(100*optimal1_cases(range,:)'*agg_mult,'k-','LineWidth',1);
    xline(optimal1.TTHI,'k:','LineWidth',1);
    hold off;
    xlim(x_lim);
    ylim(y_lims(1,:));
    if i==1
        ylabel("Exposed + Infected (%)");
    else
        set(gca,'yticklabels',[]);
    end
    set(gca,'xticklabels',[]);
    title(sprintf("%s \n (%s)",mitigation_labels{1},category_labels(i)));
    
    % Policy 2 Infections
    nexttile(2*i);
    hold on;
    grid on;
    plot(100*(nasem2_cases(range,:))'*agg_mult,'r-','LineWidth',1);
    xline(nasem2.TTHI,'r:','LineWidth',1);
    plot(100*(random2_cases(range,:))'*agg_mult,'b-','LineWidth',1);
    xline(random2.TTHI,'b:','LineWidth',1);
    plot(100*optimal2_cases(range,:)'*agg_mult,'k-','LineWidth',1);
    xline(optimal2.TTHI,'k:','LineWidth',1);
    hold off;
    xlim(x_lim);
    ylim(y_lims(1,:));
    set(gca,'yticklabels',[]);
    set(gca,'xticklabels',[]);
    title(sprintf("%s \n (%s)",mitigation_labels{2},category_labels(i)));
    
    % Policy 1 Deaths
    nexttile(2*(N_categories+i)-1);
    hold on;
    grid on;
    plot(10^5*(nasem1.D(range,:)+nasem1.Dx(range,:)+nasem1.Dv(range,:))'*agg_mult,'r-','LineWidth',1);
    plot(10^5*(random1.D(range,:)+random1.Dx(range,:)+random1.Dv(range,:))'*agg_mult,'b-','LineWidth',1);
    plot(10^5*(optimal1.D(range,:)+optimal1.Dx(range,:)+optimal1.Dv(range,:))'*agg_mult,'k-','LineWidth',1);
    hold off;
    xlim(x_lim);
    ylim(y_lims(2*i,:));
    if i==1
        ylabel("Cumulative Mortality per 100K");
    else
%         set(gca,'yticklabels',[]);
    end
    
    % Policy 2 Deaths
    ax=nexttile(2*(N_categories+i));
    hold on;
    grid on;
    plot(10^5*(nasem2.D(range,:)+nasem2.Dx(range,:)+nasem2.Dv(range,:))'*agg_mult,'r-','LineWidth',1);
    plot(10^5*(random2.D(range,:)+random2.Dx(range,:)+random2.Dv(range,:))'*agg_mult,'b-','LineWidth',1);
    plot(10^5*(optimal2.D(range,:)+optimal2.Dx(range,:)+optimal2.Dv(range,:))'*agg_mult,'k-','LineWidth',1);
    hold off;
    xlim(x_lim);
    ylim(y_lims(2*i+1,:));
%     set(gca,'yticklabels',[]);

%     % Policy 1 YLL
%     nexttile(2*(2*N_categories+i)-1);
%     hold on;
%     grid on;
%     plot(sum(nasem1.YLL(range,:)),'r-','LineWidth',1);
%     plot(sum(random1.YLL(range,:)),'b-','LineWidth',1);
%     plot(sum(optimal1.YLL(range,:)),'k-','LineWidth',1);
%     hold off;
%     xlim(x_lim);
% %     ylim(y_lims(end,:));
%     if i==1
%         ylabel("Cumulative Years of Life Lost");
%     else
% %         set(gca,'yticklabels',[]);
%     end
%     
%     % Policy 2 YLL
%     ax=nexttile(2*(2*N_categories+i));
%     hold on;
%     grid on;
%     plot(sum(nasem1.YLL(range,:)),'r-','LineWidth',1);
%     plot(sum(random1.YLL(range,:)),'b-','LineWidth',1);
%     plot(sum(optimal1.YLL(range,:)),'k-','LineWidth',1);
%     hold off;
%     xlim(x_lim);
% %     ylim(y_lims(end,:));
% %     set(gca,'yticklabels',[]);
end
lgd = legend(ax,'NASEM','None','Optimal');
lgd.Title.String = 'Prioritization Policy';
lgd.Layout.Tile = 'east';
xlabel(t,'Days');
end

