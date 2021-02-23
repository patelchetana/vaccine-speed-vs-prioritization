clear;

% Main specification
if ~exist('flags')
    flags = struct();
end

if ~isfield(flags,'run_baseline')
    flags.run_baseline = true;
end
if ~isfield(flags,'run_optimal')
    flags.run_optimal = true;
end
if ~isfield(flags,'no_hes')
    flags.no_hes = false;
end
if ~isfield(flags,'JJ_eff')
    flags.JJ_eff = false;
end
if ~isfield(flags,'high_mit')
    flags.high_mit = false;
end

% Set the mitigation parameter
flags.mit_theta = sqrt(1.5/2.6);

% Baseline
folder = "baseline/";

% Optimization Flows
opt_flows = [30 15 20 25 35 40]; % flows to run optimization for
if flags.JJ_eff
    opt_flows = [30];
    folder    = "JJ_eff/";
end
if flags.no_hes
    opt_flows = [30];
    folder    = "no_hes/";
end
if flags.high_mit
    flags.mit_theta = sqrt(1.3/2.6);
    opt_flows = [30];
    folder    = "rt_13/";
end

generate(flags,folder,opt_flows);

function generate(flags,folder,opt_flows)
%%% Generate the data needed for all the figures needed
% Read in parameters for the simulation
file_paths = readtable('../file_paths.csv','ReadRowNames',true,'Delimiter',',');
baseline_data = readtable([file_paths.path{'inputs'} 'national_baseline.csv']);
init_state    = readtable([file_paths.path{'inputs'} 'seir_model_nat_lvl_infection_starting_points_by_group_on_dec14.csv']);
contact_matrix = csvread([file_paths.path{'inputs'} 'age_group_contact_matrix.csv'],1,1);
if flags.no_hes
    options.uptake = [0 0.40497 1 1 1 1 1 1 1]'; 
    vaccine_data   = readtable([file_paths.path{'inputs'} 'vaccination_flow_strategy_nasem_unreserved_with_100p_uptake_param.csv']);
else
    options.uptake = baseline_data.vax_uptake_census;
    vaccine_data  = readtable([file_paths.path{'inputs'} 'vaccination_flow_strategy_nasem_unreserved_with_census_uptake_param.csv']);
end
if flags.JJ_eff
    options.vax_eff = 0.8;
else
    options.vax_eff = baseline_data.average_2vax_efficacy;
end

% Define other parameters
T_basic = 800; % Number of days to run for NASEM and random
T_optimal = 365; % Number of days to run for optimizing
high = 100; % Highest vaccine flow to simulate
int = 2.5; % Spacing between vaccine flows simulated
smoothing_param = 5e-2; % Smoothing parameter for optimal

head = strcat(file_paths.path{'outputs'},folder); % Folder directory has to be created beforehand

if flags.run_baseline
    % Run simulations for random and NASEM
    tic;
    for flow=0:int:high
        for mit_type=1:3
            simul = struct('constV',false,'T',T_basic,'vflow',flow*1e6/30);

            obj = seirclass(baseline_data,vaccine_data,init_state,contact_matrix,simul,options);
            % Constant mitigation
            if mit_type==1
                cal_mit = 0;
                obj.theta = repelem(flags.mit_theta,obj.G,obj.T);
                out_file_tail = sprintf("f%d_const_mit.mat",10*flow);
            % No mitigation
            elseif mit_type==2
                cal_mit = 0;
                out_file_tail = sprintf("f%d_no_mit.mat",10*flow);
            % Calibrated
            else
                cal_mit = 1;
                out_file_tail = sprintf("f%d_cal_mit.mat",10*flow);
            end

            % Run
            obj_nasem  = obj.sir_sim(2,cal_mit);
            obj_random = obj.sir_sim(3,cal_mit);

            % Save the objects
            save(strcat(head,"nasem_",out_file_tail),'obj_nasem');
            save(strcat(head,"random_",out_file_tail),'obj_random');
        end
    end
    toc;
end

% Run optimal simulations
if flags.run_optimal
    parfor i=1:length(opt_flows)
        display(['Running Flow: ' num2str(opt_flows(i))]);
        for mit_type=1:3
            flow = opt_flows(i);
            simul = struct('constV',false,'T',T_optimal,'vflow',flow*1e6/30);
            obj = seirclass(baseline_data,vaccine_data,init_state,contact_matrix,simul,options);

            % Constant mitigation
            if mit_type==1
                cal_mit = 0;
                obj.theta = repelem(flags.mit_theta,obj.G,obj.T);
                out_file_tail = sprintf("f%d_const_mit.mat",10*flow);
            % No mitigation
            elseif mit_type==2
                cal_mit = 0;
                out_file_tail = sprintf("f%d_no_mit.mat",10*flow);
            % Calibrated
            else
                cal_mit = 1;
                out_file_tail = sprintf("f%d_cal_mit.mat",10*flow);
            end

            % Optimize
            obj.Vdot = obj.optimize(obj.Vdot'*obj.share,@(x)(obj.death_rate(x,cal_mit)),smoothing_param);
            obj_opt = obj.sir_sim(2,cal_mit);

            % Save the objects
            save_optimal(strcat(head,"optimal_",out_file_tail),obj_opt);
        end
    end
end
end

% Helper function to save the optimal object in a parfor loop
function save_optimal(fname,obj_opt)
  save(fname, 'obj_opt')
end