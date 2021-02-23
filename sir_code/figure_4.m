function figure_4(selection_range,N_mitigation,nasem,random,int,line_params,legend_labels)
% Compute threshold increase by deaths
addition_D = zeros(length(selection_range),N_mitigation);
for t=1:N_mitigation
    for i=selection_range
        for j=1:size(random.D,1)-1
            if nasem.D(i,t) <= random.D(j,t) && nasem.D(i,t) >= random.D(j+1,t)
                increase = int*(j-i) + int*(random.D(j,t)-nasem.D(i,t))/(random.D(j,t)-random.D(j+1,t));
                assert((random.D(j,t)-nasem.D(i,t))/(random.D(j,t)-random.D(j+1,t)) < 1);
                addition_D(i-min(selection_range)+1,t) = increase;
                break
            end
        end
    end
end
% Compute threshold increase by YLL
addition_YLL = zeros(length(selection_range),N_mitigation);
for t=1:N_mitigation
    for i=selection_range
        for j=1:size(random.D,1)-1
            if nasem.YLL(i,t) == random.YLL(j,t)
                addition_YLL(i-min(selection_range)+1,t) = int*(j-i);
            elseif nasem.YLL(i,t) < random.YLL(j,t) && nasem.YLL(i,t) > random.YLL(j+1,t)
                increase = int*(j-i) + int*(random.YLL(j,t)-nasem.YLL(i,t))/(random.YLL(j,t)-random.YLL(j+1,t));
                assert((random.YLL(j,t)-nasem.YLL(i,t))/(random.YLL(j,t)-random.YLL(j+1,t)) < 1);
                addition_YLL(i-min(selection_range)+1,t) = increase;
                break
            end
        end
    end
end
% % Compute threshold increase by cases
% addition_cases = zeros(length(selection_range),N_mitigation);
% for t=1:N_mitigation
%     for i=selection_range
%         for j=1:size(random.D,1)-1
%             if nasem.cases(i,t) == random.cases(j,t)
%                 addition_cases(i-min(selection_range)+1,t) = int*(j-i);
%             elseif nasem.cases(i,t) < random.cases(j,t) && nasem.cases(i,t) > random.cases(j+1,t)
%                 increase = int*(j-i) + int*(random.cases(j,t)-nasem.cases(i,t))/(random.cases(j,t)-random.cases(j+1,t));
%                 assert((random.cases(j,t)-nasem.cases(i,t))/(random.cases(j,t)-random.cases(j+1,t)) < 1);
%                 addition_cases(i-min(selection_range)+1,t) = increase;
%                 break
%             end
%         end
%     end
% end
top = max([addition_D addition_YLL],[],'all');

figure;
p=tiledlayout(1,2,'TileSpacing','Compact');
x = (selection_range-1)*int;

nexttile;
hold on;
grid on;
for t=1:N_mitigation
    plot(x,addition_D(:,t),line_params(t),'LineWidth',1);
end
hold off;
ylim([0 top]);
title('Holding Deaths Constant');

ax=nexttile;
hold on;
grid on;
for t=1:N_mitigation
    plot(x,addition_YLL(:,t),line_params(t),'LineWidth',1);
end
hold off;
set(gca,'yticklabel',[]);
ylim([0 top]);
title('Holding Years of Life Lost Constant');

% ax=nexttile;
% hold on;
% grid on;
% for t=1:N_mitigation
%     plot(x,addition_cases(:,t),line_params(t),'LineWidth',1);
% end
% hold off;
% title('Holding Total Cases Constant');

lgd = legend(ax,legend_labels);
lgd.Title.String = 'Mitigation Scenario';
lgd.Layout.Tile = 'east';
xlabel(p,'Millions Vaccinated Per Month Under NASEM');
ylabel(p,{'Increase in Millions Vaccinated Per Month','Under No Prioritization'});
end

