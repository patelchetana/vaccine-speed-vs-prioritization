clear;

% Main specification
flags.no_hes       = false;
flags.JJ_eff       = false;
flags.high_mit     = false;
flags.fifteen_mil  = false;

make(flags,"baseline/","baseline/");

% No hesitation
flags.no_hes = true;
make(flags,"no_hes/","no_hes/");
flags.no_hes = false;

% Lower efficacy
flags.JJ_eff = true;
make(flags,"JJ_eff/","JJ_eff/");
flags.JJ_eff = false;

% Higher mitigation
flags.high_mit = true;
make(flags,"rt_13/","rt_13/");
flags.high_mit = false;

% 15M
flags.fifteen_mil = true;
make(flags,"baseline/","15_mil/");
flags.fifteen_mil = false;

function make(flags,in_folder,out_folder)
% Read data
file_paths = readtable('../file_paths.csv','ReadRowNames',true,'Delimiter',',');
mat_in       = strcat(file_paths.path{'outputs'},in_folder);
exhibits_out = strcat(file_paths.path{'exhibits'},out_folder);
baseline_data = readtable([file_paths.path{'inputs'} 'national_baseline.csv']);
if flags.no_hes
    vaccine_data   = readtable([file_paths.path{'inputs'} 'vaccination_flow_strategy_nasem_unreserved_with_100p_uptake_param.csv']);
else
    vaccine_data  = readtable([file_paths.path{'inputs'} 'vaccination_flow_strategy_nasem_unreserved_with_malik_uptake_param.csv']);
end
init_state    = readtable([file_paths.path{'inputs'} 'seir_model_nat_lvl_infection_starting_points_by_group_on_dec14.csv']);
contact_matrix = csvread([file_paths.path{'inputs'} 'age_group_contact_matrix.csv'],1,1);

N_group = height(baseline_data);
share = baseline_data.sh_group;

% Load data for particular policies
if flags.fifteen_mil
    flow = 15;
else
    flow = 30;
end
flow1 = flow*10; % 30 million per month
flow2 = flow*10;
mit1 = "const";
mit2 = "cal";


load(strcat(mat_in,sprintf("nasem_f%d_%s_mit.mat",flow1,mit1)));
load(strcat(mat_in,sprintf("random_f%d_%s_mit.mat",flow1,mit1)));
load(strcat(mat_in,sprintf("optimal_f%d_%s_mit.mat",flow1,mit1)));
nasem1 = obj_nasem;
random1 = obj_random;
opt1 = obj_opt.calc_cases;

load(strcat(mat_in,sprintf("nasem_f%d_%s_mit.mat",flow2,mit2)));
load(strcat(mat_in,sprintf("random_f%d_%s_mit.mat",flow2,mit2)));
load(strcat(mat_in,sprintf("optimal_f%d_%s_mit.mat",flow2,mit2)));
nasem2 = obj_nasem;
random2 = obj_random;
opt2 = obj_opt.calc_cases;

%% Figure 1
% Types of policies used
type_labels = [sprintf("Sustained Mitigation \n %dM Vaccinated per month",flow),...
    sprintf("Calibrated Mitigation \n %dM Vaccinated per month",flow)];
% Age groups used
categories = [7 9;4 6;1 3];
category_labels = {'60+','30-59','0-29'};

figure_1(categories,nasem1.V_admin,random1.V_admin,opt1.V_admin,...
    nasem2.V_admin,random2.V_admin,opt2.V_admin,...
    contact_matrix',baseline_data.c19_ifr_group,share,[0 365],...
    category_labels,type_labels);
set(gcf,'Position',[100 100 1000 500]);
saveas(gcf,strcat(exhibits_out,'figure1.png'));

%% Figure 2
mitigation_labels = ["Sustained Mitigation","Calibrated Mitigation"]; % Types of mitigation used
% Age groups used
categories = [1 9;1 6;7 9];
category_labels = ["All","Less than 60 years old","60 years or older"];
figure_2(categories,share,nasem1,random1,opt1,...
                          nasem2,random2,opt2,...
                          [1 365], [0 2.5;0 300;0 300;0 50;0 50;0 1200;0 1200;0 0.4],...
                          mitigation_labels,category_labels);
set(gcf,'Position',[100 100 1500 500]);
saveas(gcf,strcat(exhibits_out,'figure2.png'));

%%
% Read in objects and compute deaths and cases
int = 2.5;
N_flows = 100/int+1;
nasem.D         = zeros(N_flows,3);
nasem.TTHI       = zeros(N_flows,3);
nasem.YLL       = zeros(N_flows,3);
nasem.cases = zeros(N_flows,3);
random  = nasem;
optimal = nasem;

policies = ["no","cal","const"];
for f=1:N_flows
    for t=1:3
        load(strcat(mat_in,sprintf("nasem_f%d_%s_mit.mat",(f-1)*int*10,policies(t))));
        load(strcat(mat_in,sprintf("random_f%d_%s_mit.mat",(f-1)*int*10,policies(t))));
        
        if isfile(strcat(mat_in,sprintf("optimal_f%d_%s_mit.mat",(f-1)*int*10,policies(t))))
            load(strcat(mat_in,sprintf("optimal_f%d_%s_mit.mat",(f-1)*int*10,policies(t))));
            optimal.TTHI(f,t) = obj_opt.TTHI;
            optimal.YLL(f,t) = sum(obj_opt.YLL(:,end));
            optimal.D(f,t) = (obj_opt.D(:,end)+obj_opt.Dx(:,end)+obj_opt.Dv(:,end))'*share;
            optimal.cases(f,t) = obj_opt.cases(:,end)'*share;
        end
        
        % TTHI
        nasem.TTHI(f,t) = obj_nasem.TTHI;
        random.TTHI(f,t) = obj_random.TTHI;

        % YLL
        nasem.YLL(f,t) = sum(obj_nasem.YLL(:,end));
        random.YLL(f,t) = sum(obj_random.YLL(:,end));

        % Deaths
        nasem.D(f,t) = (obj_nasem.D(:,end)+obj_nasem.Dx(:,end)+obj_nasem.Dv(:,end))'*share;
        random.D(f,t) = (obj_random.D(:,end)+obj_random.Dx(:,end)+obj_random.Dv(:,end))'*share;
        
        % Incidence/Infections
        nasem.cases(f,t) = obj_nasem.cases(:,end)'*share;
        random.cases(f,t) = obj_random.cases(:,end)'*share;
    end
end

% Select range to plot for figures 3 and 4
low = 15;
high = 40;
selection_range = (low/int+1):(high/int+1);

%% Figure 3
fig_3_range.nasem   = selection_range;
fig_3_range.random  = selection_range;
fig_3_range.optimal = (low/int+1):2:(high/int+1);
fig_3_mit_labels = {'No Mitigation', 'Calibrated Mitigation','Sustained Mitigation'};
baseline.D         = random.D(1,1);
baseline.YLL       = random.YLL(1,1);
baseline.cases = random.cases(1,1);
figure_3(fig_3_range,3,nasem,random,optimal,baseline,int,fig_3_mit_labels)
set(gcf,'Position',[100 100 1000 500]);
saveas(gcf,strcat(exhibits_out,'figure3.png'));

%% Figure 4
fig_4_mit_labels = {'Calibrated', 'Sustained'};
line_params = ["k-","k--"];
fig_4_nasem.D = nasem.D(:,2:3);
fig_4_nasem.YLL = nasem.YLL(:,2:3);
fig_4_nasem.cases = nasem.cases(:,2:3);
fig_4_random.D = random.D(:,2:3);
fig_4_random.YLL = random.YLL(:,2:3);
fig_4_random.cases = random.cases(:,2:3);
figure_4(selection_range,2,fig_4_nasem,fig_4_random,int,line_params,fig_4_mit_labels)
set(gcf,'Position',[100 100 1000 500]);
saveas(gcf,strcat(exhibits_out,'figure4.png'));

%% Figure 5
figure_5(nasem2,random2,opt2,random1.theta(1,1),[0 365],[0.3 1.05]);
set(gcf,'Position',[100 100 1000 500]);
saveas(gcf,strcat(exhibits_out,'figure5.png'));

%% Figure 6
figure_6(fig_3_range,3,nasem,random,optimal,baseline,int,fig_3_mit_labels)
set(gcf,'Position',[100 100 1000 500]);
saveas(gcf,strcat(exhibits_out,'figure6.png'));

end
