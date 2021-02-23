% Documentation for the code
% Pass init_state = [] or do not pass if we want to start with a
% fully susceptible population


clear;
file_paths = readtable('../file_paths.csv','ReadRowNames',true);
baseline_data = readtable([file_paths.path{'inputs'} 'national_baseline.csv']);
vaccine_data  = readtable([file_paths.path{'inputs'} 'vaccination_flow_strategy_nasem_unreserved_with_malik_uptake_param.csv']);
init_state    = readtable([file_paths.path{'inputs'} 'nat_lvl_infection_starting_points_by_group_on_dec14.csv']);
contact_matrix = csvread([file_paths.path{'inputs'} 'age_group_contact_matrix.csv'],1,1);

simul.constV = false;
simul.T = 750;
simul.vflow = 15*1e6/30;

% vsirm = vsirmclass(baseline_data,vaccine_data,init_state,contact_matrix,simul);
% vsirm = vsirm.sir_sim(3,1);
% sum(vsirm.D'*vsirm.share)
seir = seirclass(baseline_data,vaccine_data,init_state,contact_matrix,simul);
seir = seir.sir_sim(3,1);
% seir.D(:,end)'*seir.share



