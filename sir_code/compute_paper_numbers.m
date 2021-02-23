clear;
file_paths = readtable('../file_paths.csv','ReadRowNames',true,'Delimiter',',');
mat_in       = strcat(file_paths.path{'outputs'},"baseline/");
baseline_data = readtable([file_paths.path{'inputs'} 'national_baseline.csv']);
nat_pop = baseline_data.nat_population(1);

policies = ["cal","const"];
old_range = 7:9;
young_range = 1:6;
diff = zeros(2,1);
mort_old_opt = zeros(2,1);
mort_old_other = zeros(2,2);
mort_young = zeros(2,3);
for i=1:2
    load(strcat(mat_in,sprintf("nasem_f300_%s_mit.mat",policies(i))));
    load(strcat(mat_in,sprintf("random_f300_%s_mit.mat",policies(i))));
    load(strcat(mat_in,sprintf("optimal_f300_%s_mit.mat",policies(i))));
    nasem_deaths = obj_nasem.all_D;
    random_deaths = obj_random.all_D;
    opt_deaths = obj_opt.all_D;
    diff(i) = (nasem_deaths(:,end)'*obj_nasem.share - opt_deaths(:,end)'*obj_opt.share)*nat_pop;
    agg_mult = obj_opt.share(old_range)/sum(obj_opt.share(old_range))*1e5;
    mort_old_opt(i) = opt_deaths(old_range,end)'*agg_mult;
    mort_old_other(i,1) = nasem_deaths(old_range,end)'*agg_mult;
    mort_old_other(i,2) = random_deaths(old_range,end)'*agg_mult;
    agg_mult = obj_opt.share(young_range)/sum(obj_opt.share(young_range))*1e5;
    mort_young(i,1) = opt_deaths(young_range,end)'*agg_mult;
    mort_young(i,2) = nasem_deaths(young_range,end)'*agg_mult;
    mort_young(i,3) = random_deaths(young_range,end)'*agg_mult;
    
end
fprintf("Difference in deaths between random and NASEM -- min: %f, max: %f \n",min(diff),max(diff));
fprintf("Mortality for 60+ under optimal (per 100K) -- min: %f, max: %f \n",min(mort_old_opt),max(mort_old_opt));
fprintf("Mortality for 60+ under NASEM and random (per 100K) -- min: %f, max: %f \n",min(mort_old_other,[],'all'),max(mort_old_other,[],'all'));
fprintf("Mortality for 0-59 (per 100K) -- min: %f, max: %f \n",min(mort_young,[],'all'),max(mort_young,[],'all'));


load(strcat(mat_in,"optimal_f400_no_mit.mat"));
opt_deaths = obj_opt.all_D;
opt_deaths = opt_deaths(:,end)'*obj_opt.share;
int = 2.5;
N_flows = 100/int+1;
for f=1:N_flows
    load(strcat(mat_in,sprintf("random_f%d_cal_mit.mat",(f-1)*int*10)));
    random_deaths = obj_random.all_D;
    random_deaths = random_deaths(:,end)'*obj_random.share;
    if random_deaths <= opt_deaths
        fprintf("Random with %2.1f vpm and calibrated mitigation does no worse than optimal with 40 vpm under no mitigation \n",(f-1)*int);
        break
    end
end
for f=1:N_flows
    load(strcat(mat_in,sprintf("random_f%d_const_mit.mat",(f-1)*int*10)));
    random_deaths = obj_random.all_D;
    random_deaths = random_deaths(:,end)'*obj_random.share;
    if random_deaths <= opt_deaths
        fprintf("Random with %2.1f vpm and sustained mitigation does no worse than optimal with 40 vpm under no mitigation \n",(f-1)*int);
        break
    end
end

load(strcat(mat_in,"nasem_f300_cal_mit.mat"));
nasem_deaths = obj_nasem.all_D;
nasem_deaths = nasem_deaths(:,end)'*obj_nasem.share;
int = 2.5;
N_flows = 100/int+1;
for f=1:N_flows
    load(strcat(mat_in,sprintf("random_f%d_const_mit.mat",(f-1)*int*10)));
    random_deaths = obj_random.all_D;
    random_deaths = random_deaths(:,end)'*obj_random.share;
    if random_deaths <= nasem_deaths
        fprintf("Random with %2.1f vpm and sustained mitigation does no worse than NASEM with 30 vpm under calibrated mitigation \n",(f-1)*int);
        break
    end
end

load(strcat(mat_in,"nasem_f300_cal_mit.mat"));
load(strcat(mat_in,"random_f300_cal_mit.mat"));
nasem_deaths = obj_nasem.all_D;
nasem_deaths = nasem_deaths(:,end)'*obj_nasem.share;
random_deaths = obj_random.all_D;
random_deaths = random_deaths(:,end)'*obj_random.share;
fprintf("Abandoning prioritization results in %f more deaths under calibrated mitigation \n",(random_deaths-nasem_deaths)*nat_pop);

load(strcat(mat_in,"nasem_f300_const_mit.mat"));
load(strcat(mat_in,"random_f300_const_mit.mat"));
nasem_deaths = obj_nasem.all_D;
nasem_deaths = nasem_deaths(:,end)'*obj_nasem.share;
random_deaths = obj_random.all_D;
random_deaths = random_deaths(:,end)'*obj_random.share;
fprintf("Abandoning prioritization results in %f more deaths under sustained mitigation \n",(random_deaths-nasem_deaths)*nat_pop);

% Read in objects and compute deaths and cases
int = 2.5;
N_flows = 100/int+1;
random_deaths = zeros(N_flows,2);

policies = ["const","cal"];
for f=1:N_flows
    for t=1:2
        load(strcat(mat_in,sprintf("random_f%d_%s_mit.mat",(f-1)*int*10,policies(t))));
        random_deaths(f,t) = (obj_random.D(:,end)+obj_random.Dx(:,end)+obj_random.Dv(:,end))'*obj_random.share;
    end
end

load(strcat(mat_in,"nasem_f300_const_mit.mat"))
nasem_deaths = (obj_nasem.D(:,end)+obj_nasem.Dx(:,end)+obj_nasem.Dv(:,end))'*obj_nasem.share;
for j=1:size(random_deaths,1)-1
    if nasem_deaths <= random_deaths(j,1) && nasem_deaths >= random_deaths(j+1,1)
        increase = int*(j-13) + int*(random_deaths(j,1)-nasem_deaths)/(random_deaths(j,1)-random_deaths(j+1,1));
        fprintf("Increase in speed under random to keep deaths equivalent to NASEM at 30 vpm with sustained mitigation: %f vpm, %f %% increase \n",increase,100*increase/30);
        break
    end
end

load(strcat(mat_in,"nasem_f300_cal_mit.mat"))
nasem_deaths = (obj_nasem.D(:,end)+obj_nasem.Dx(:,end)+obj_nasem.Dv(:,end))'*obj_nasem.share;
for j=1:size(random_deaths,1)-1
    if nasem_deaths <= random_deaths(j,2) && nasem_deaths >= random_deaths(j+1,2)
        increase = int*(j-13) + int*(random_deaths(j,2)-nasem_deaths)/(random_deaths(j,2)-random_deaths(j+1,2));
        fprintf("Increase in speed under random to keep deaths equivalent to NASEM at 30 vpm with sustained mitigation: %f vpm, %f %% increase \n",increase,100*increase/30);
        break
    end
end
