function figure_1(categories,nasem_V1,random_V1,opt_V1,...
    nasem_V2,random_V2,opt_V2,...
    contact_matrix_raw,ifr,share,x_lim,category_labels,type_labels)
% Collapse contact matrix and compute IFRs
N_cat = size(categories,1);
contact_matrix_temp = zeros(N_cat,size(contact_matrix_raw,1));
for i=1:N_cat
    range = categories(i,1):categories(i,2);
    contact_matrix_temp(i,:) = contact_matrix_raw(:,range)*share(range)/sum(share(range));
end
contact_matrix = zeros(N_cat,N_cat);
ifr_rates = repelem(0,1,N_cat);
for i=1:N_cat
    range = categories(i,1):categories(i,2);
    contact_matrix(i,:) = contact_matrix_temp(:,range)*share(range)/sum(share(range));
    ifr_rates(i) = ifr(range)'*share(range)/sum(share(range));
end

figure;
t = tiledlayout(N_cat,3,'TileSpacing','Compact');
for i=1:N_cat
    for j=1:3
        ax=nexttile;
        % Contact Matrix
        if j==1
            if i<N_cat
                b = bar(categorical(category_labels),contact_matrix(i,:),'Facecolor','flat');
                set(gca,'XTick',[],'XTickLabel',[]);
            else
                b = bar(categorical(category_labels),contact_matrix(i,:),'Facecolor','flat');
                xlabel('Age Groups');
            end  
            % b.CData = colors;
            ylim([0 4]);
            if i==1
                title('Daily Contacts');
            end
            ylabel({['Age ' category_labels{i}],['IFR: ' num2str(ifr_rates(i),'%.3f')]},'FontSize',10,'FontWeight','bold');
            
        % Vaccination Policies
        else
            if i==1
                title(sprintf(type_labels(j-1)));
            end
            % Plot the graphs
            hold on;
            grid on;
            range = categories(i,1):categories(i,2); % Range to subset
            agg_mult = share(range)/sum(share(range)); % Multiplier to get weighted average
            if j==2
                plot(100*nasem_V1(range,:)'*agg_mult,'r-','LineWidth',1);
                plot(100*random_V1(range,:)'*agg_mult,'b-','LineWidth',1);
                plot(100*opt_V1(range,:)'*agg_mult,'k-','LineWidth',1);
            else
                plot(100*nasem_V2(range,:)'*agg_mult,'r-','LineWidth',1);
                plot(100*random_V2(range,:)'*agg_mult,'b-','LineWidth',1);
                plot(100*opt_V2(range,:)'*agg_mult,'k-','LineWidth',1);
            end
            if i==N_cat
                xlabel('Days');
            else
                set(gca,'XTickLabel',[]);
            end
            if j~=2
                set(gca,'YTickLabel',[]);
            end
            if j==2 && i==ceil(N_cat/2)
                ylabel('Percent Vaccinated (Cumulative)','FontSize',12);
            end
            xlim(x_lim);
            ylim([0 100]);
            hold off;
        end
    end
end
lgd = legend(ax,'NASEM','None','Optimal');
lgd.Title.String = 'Prioritization Policy';
lgd.Layout.Tile = 'east';
end

