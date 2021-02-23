% Class Definition: MG-VSIR-M Model
classdef seirclass
% Parameters and objects
    properties
        G                   % Number of groups
        group_names         % Group identifier
        share               % Share in each group -- vector size 1xG
        gamma               % 1/expected incubation time E->I (size 1x1)
        nu                  % 1/expected recovery time I->R (size 1x1)
        beta                % Contact parameter 1x1
        rho                 % Contact matrix GxG
        T                   % Number of simulation days
        V                   % Fully protected population
        V_admin             % Population that received a vaccine
        S                   % Susceptible GxT
        Sx                  % Vaccine Unprotected/Hesitant S
        Sv                  % Vaccine ineffective S
        E                   % Exposed GxT
        Ex                  % Vaccine Unprotected/Hesitant E
        Ev                  % Vaccine ineffective E
        I                   % Infected GxT
        Ix                  % Vaccine Unprotected/Hesitant I
        Iv                  % Vaccine ineffective I
        R                   % Resolved GxT
        Rx                  % Vaccine Unprotected/Hesitant R
        Rv                  % Vaccine ineffective R
        D                   % Dead GxT
        Dx                  % Vaccine Unprotected/Hesitant D
        Dv                  % Vaccine ineffective D
        YLL                 % Years of lives lost
        yll_times_n         % Group-specific years of lives lost x N
        cases               % Total Cases
        uptake              % Fraction of each group that will take vaccine
        ifr                 % IFR -- Gx1
        vax_eff             % Vaccine efficacy -- Gx1
        vax_ifr             % Vaccinated IFR -- Gx1
        S0                  % Initial susceptible -- Gx1
        E0                  % Initial susceptible -- Gx1
        I0                  % Initial infected -- Gx1
        R0                  % Initial recovered -- Gx1
        D0                  % Initial dead -- Gx1
        Rnaught             % Rnaught -- 1x1
        Vdot                % Flow vaccination rate -- G x T
        theta               % Mitigation parameters -- 1 x T
        vaccinate           % {0,1,2}, 0: no vaccines, 1: S+E, 2: S+E+R, 3: Random S+E+R
        mitigate            % Boolean, mitigation to R=1 or not
        constV              % Boolean, constant vaccination rate
        date                % Date the object was created (timestamp)
        TTHI                % Time till herd immunity
    end
    
    methods
        % Initialize class properties
        function [obj] = seirclass(data,vaccine_data,init_state,contact_matrix,simul,options)
            if ~exist('options','var')
                options = struct();
            end
            if ~exist('init_state','var')
                init_state = [];
            end
            if ~exist('contact_matrix','var')
                contact_matrix = [];
            end

            % Basic parameters
            obj.date = date;
                

            % Read the necessary parameters
            obj.G           = height(data); 
            obj.group_names = data.group; 
            obj.share       = data.sh_group;
            obj.ifr         = data.c19_ifr_group;
            obj.uptake      = options.uptake;
            obj.vax_eff     = options.vax_eff;
            obj.yll_times_n = data.yll.*data.n_group;

            % Set vaccine ifr
            if ~isfield(data,'vax_ifr_group')
                obj.vax_ifr = 0;
            else
                obj.vax_ifr = data.vax_ifr_group;
            end
            
            if ~isfield(options,'R')
                obj.Rnaught     = 2.6;
            else
                obj.Rnaught     = options.Rnaught;
            end
            if ~isfield(options,'gamma')
                obj.gamma     = 1/4;
                options.gamma = 1/4;
            else
                obj.gamma = options.gamma;
            end
            if ~isfield(options,'nu')
                obj.nu     = 1/9;
                options.nu = 1/9;
            else
                obj.nu = options.nu;
            end 
            
            % Contact matrix
            if isempty(contact_matrix)
                obj.rho    = repmat(obj.share',obj.G,1);
            else
                obj.rho    = contact_matrix;
            end

            % Set beta
            obj.beta   = data.sus_to_inf * obj.Rnaught/max(abs(eig(diag(data.sus_to_inf)*contact_matrix*1/obj.nu)));

            if isempty(init_state)
                obj.S0     = ones(obj.G,1)-10^-6;
                obj.E0     = zeros(obj.G,1);
                obj.I0     = ones(obj.G,1)*10^-6;
            else
                obj.S0     = init_state.s0;
                obj.E0     = init_state.e0;
                obj.I0     = init_state.i0;
            end
            resolved = 1-(obj.S0+obj.E0+obj.I0);
            obj.R0   = resolved.*(1-obj.ifr/100);
            obj.D0   = resolved.*(obj.ifr/100);

            obj.T      = simul.T;
            
            % Calculate Vdot unless
            if simul.constV
                obj.Vdot  = repmat(simul.v,1,simul.T);
            else
                obj.Vdot  = calc_Vdot(obj.group_names,vaccine_data,simul);
            end
            obj.constV = simul.constV;
            
            % Initialize
            obj.S  = zeros(obj.G,obj.T);
            obj.Sx = zeros(obj.G,obj.T);
            obj.Sv = zeros(obj.G,obj.T);
            
            obj.E  = zeros(obj.G,obj.T);
            obj.Ex = zeros(obj.G,obj.T);
            obj.Ev = zeros(obj.G,obj.T);
            
            obj.I  = zeros(obj.G,obj.T);
            obj.Ix = zeros(obj.G,obj.T);
            obj.Iv = zeros(obj.G,obj.T);
            
            obj.R  = zeros(obj.G,obj.T);
            obj.Rx = zeros(obj.G,obj.T);
            obj.Rv = zeros(obj.G,obj.T);
            
            obj.D  = zeros(obj.G,obj.T);
            obj.Dx = zeros(obj.G,obj.T);
            obj.Dv = zeros(obj.G,obj.T);
            
            obj.V  = zeros(obj.G,obj.T);
            obj.V_admin = zeros(obj.G,obj.T);
            
            obj.theta = ones(obj.G,obj.T-1);
            
        end
        
        % SIR Simulation backbone
        function obj = sir_sim(obj,vaccinate,mitigate)
            % Store init
            obj.S(:,1)  = obj.uptake.*obj.S0;
            obj.Sx(:,1) = (1-obj.uptake).*obj.S0;
            obj.E(:,1)  = obj.uptake.*obj.E0;
            obj.Ex(:,1) = (1-obj.uptake).*obj.E0;
            obj.I(:,1)  = obj.uptake.*obj.I0;
            obj.Ix(:,1) = (1-obj.uptake).*obj.I0;
            obj.R(:,1)  = obj.uptake.*obj.R0;
            obj.Rx(:,1) = (1-obj.uptake).*obj.R0;
            obj.D(:,1)  = obj.uptake.*obj.D0;
            obj.Dx(:,1) = (1-obj.uptake).*obj.D0;
            
            % Loop
            for tt = 2:obj.T
                if mitigate == 1
                     thetasq = solve_theta_R_zero(obj.S(:,tt-1)+obj.Sx(:,tt-1)+obj.Sv(:,tt-1),...
                         obj.I(:,tt-1)+obj.Ix(:,tt-1)+obj.Iv(:,tt-1),...
                         obj.share, obj.nu, obj.beta, obj.rho);
                     obj.theta(:,tt-1) = sqrt(thetasq);
                else
                    thetasq = obj.theta(:,tt-1).^2;
                end
                
                % Change in infections, using group interaction model
                lambda = thetasq.*obj.beta.*(obj.rho*(obj.I(:,tt-1)+obj.Ix(:,tt-1)+obj.Iv(:,tt-1)));
                
                R_dot  = obj.nu*obj.I(:,tt-1).*(1-obj.ifr/100);
                Rx_dot = obj.nu*obj.Ix(:,tt-1).*(1-obj.ifr/100);
                Rv_dot = obj.nu*obj.Iv(:,tt-1).*(1-obj.vax_ifr/100);
                
                D_dot  = obj.nu*obj.I(:,tt-1).*(obj.ifr/100);
                Dx_dot = obj.nu*obj.Ix(:,tt-1).*(obj.ifr/100);
                Dv_dot = obj.nu*obj.Iv(:,tt-1).*(obj.vax_ifr/100);
                
                S_dot  = max(-lambda.*obj.S(:,tt-1) ,-obj.S(:,tt-1));
                Sx_dot = max(-lambda.*obj.Sx(:,tt-1),-obj.Sx(:,tt-1));
                Sv_dot = max(-lambda.*obj.Sv(:,tt-1),-obj.Sv(:,tt-1));
                
                E_dot  = -S_dot - obj.gamma*obj.E(:,tt-1);
                Ex_dot = -Sx_dot - obj.gamma*obj.Ex(:,tt-1);
                Ev_dot = -Sv_dot - obj.gamma*obj.Ev(:,tt-1);
                
                I_dot  = obj.gamma*obj.E(:,tt-1) - (R_dot + D_dot);
                Ix_dot = obj.gamma*obj.Ex(:,tt-1) - (Rx_dot + Dx_dot);
                Iv_dot = obj.gamma*obj.Ev(:,tt-1) - (Rv_dot + Dv_dot);

                % Update
                obj.S(:,tt)  = obj.S(:,tt-1) + S_dot;
                obj.Sx(:,tt) = obj.Sx(:,tt-1) + Sx_dot;
                obj.Sv(:,tt) = obj.Sv(:,tt-1) + Sv_dot;
                
                obj.E(:,tt)  = obj.E(:,tt-1) + E_dot;
                obj.Ex(:,tt) = obj.Ex(:,tt-1) + Ex_dot;
                obj.Ev(:,tt) = obj.Ev(:,tt-1) + Ev_dot;
                
                obj.I(:,tt)  = obj.I(:,tt-1) + I_dot;
                obj.Ix(:,tt) = obj.Ix(:,tt-1) + Ix_dot;
                obj.Iv(:,tt) = obj.Iv(:,tt-1) + Iv_dot;
                
                obj.R(:,tt)  = obj.R(:,tt-1) + R_dot;
                obj.Rx(:,tt) = obj.Rx(:,tt-1) + Rx_dot;
                obj.Rv(:,tt) = obj.Rv(:,tt-1) + Rv_dot;
                
                obj.D(:,tt)  = obj.D(:,tt-1) + D_dot;
                obj.Dx(:,tt) = obj.Dx(:,tt-1) + Dx_dot;
                obj.Dv(:,tt) = obj.Dv(:,tt-1) + Dv_dot;
                
                % Vaccinate only S+E
                % Post vaccination update
                if vaccinate~=0
                    if vaccinate == 1 % Do not vaccinate R (serotesting)
                        vax_pop = obj.S(:,tt) + obj.E(:,tt);
                    else
                        vax_pop = obj.S(:,tt) + obj.E(:,tt) + obj.R(:,tt);
                    end

                    % Random/No Prioritization
                    if vaccinate == 3 
                        obj.Vdot(:,tt) = obj.Vdot(:,tt)'*obj.share*(vax_pop ./ (vax_pop' * obj.share));
                    end

                    V_dot     = min(obj.Vdot(:,tt), vax_pop);
                    V_S_dot   = V_dot.*obj.S(:,tt)./vax_pop;
                    V_E_dot   = V_dot.*obj.E(:,tt)./vax_pop;

                    % Deal with zero vax_pop groups
                    V_S_dot(vax_pop==0) = 0;
                    V_E_dot(vax_pop==0) = 0;
                    
                    if vaccinate~=1
                        V_R_dot   = V_dot.*obj.R(:,tt)./vax_pop;
                        V_R_dot(vax_pop==0) = 0;
                        R_dot     = - V_R_dot;
                        Rx_dot    = V_R_dot;
                    else
                        R_dot = 0;
                        Rx_dot = 0;                        
                    end

                    % V_dot effective
                    V_dot_eff = obj.vax_eff.*V_S_dot;
                    S_dot     = - V_S_dot;
                    Sv_dot    = V_S_dot - V_dot_eff;
                    E_dot     = - V_E_dot;
                    Ex_dot    = V_E_dot;

                    % Post vaccination updates
                    obj.S(:,tt)  = obj.S(:,tt)  + S_dot;
                    obj.Sv(:,tt) = obj.Sv(:,tt) + Sv_dot;
                    obj.E(:,tt)  = obj.E(:,tt)  + E_dot;
                    obj.Ex(:,tt) = obj.Ex(:,tt) + Ex_dot;
                    obj.R(:,tt)  = obj.R(:,tt)  + R_dot;
                    obj.Rx(:,tt) = obj.Rx(:,tt) + Rx_dot;
                
                    obj.V(:,tt) = obj.V(:,tt-1) + V_dot_eff;
                    obj.V_admin(:,tt) = obj.V_admin(:,tt-1) + V_dot;
                end                    
            end
            obj = obj.calc_TTHI;
            obj = obj.calc_YLL;
            obj = obj.calc_cases;
            
            % Output type
            obj.vaccinate = vaccinate;
            obj.mitigate  = mitigate;
        end
        
        function deaths = all_D(obj)
            deaths = obj.D+obj.Dx+obj.Dv;
        end
        
        function cases = all_EI(obj)
            cases = obj.I+obj.Ix+obj.Iv+obj.E+obj.Ex+obj.Ev;
        end
        
        % Time to herd immunity
        function obj = calc_TTHI(obj)
            % Change in infections, using group the interaction model
            lambda = obj.beta.*(obj.rho*(obj.I+obj.Ix+obj.Iv));
            StoE = lambda.*(obj.S+obj.Sx+obj.Sv);
            ItoR = obj.nu.*(obj.I+obj.Ix+obj.Iv);
            avEIdot = obj.share'*(StoE-ItoR);
            if all(avEIdot>=0)
                obj.TTHI = NaN;
            else
                index    = find(avEIdot<0,1);
                obj.TTHI = index + avEIdot(index)/(avEIdot(index-1) - avEIdot(index));
            end
        end
        
        % Total years of lives lost helper function
        function obj = calc_YLL(obj)
            obj.YLL   = (obj.D+obj.Dx+obj.Dv).*obj.yll_times_n;
        end

        % Total incidence helper function
        function obj  = calc_cases(obj)
            obj.cases = 1 - (obj.S+obj.Sx+obj.Sv+obj.V);
        end
        
        % Total incidence objective function
        function cases = incidence_rate(obj, Vdot, mitigation)
            obj.Vdot = reshape(Vdot, size(obj.Vdot));
            obj = obj.sir_sim(2,mitigation);
            cases = 1 - (obj.S(:,end)+obj.Sx(:,end)+obj.Sv(:,end)+obj.V(:,end))'*obj.share;
        end
        
        % Total deaths objective function
        function deaths = death_rate(obj, Vdot, mitigation)
            obj.Vdot = reshape(Vdot, size(obj.Vdot));
            obj = obj.sir_sim(2,mitigation);
            deaths = (obj.D(:,end)+obj.Dx(:,end)+obj.Dv(:,end))'*obj.share;
        end
        
        % Optimize using knitro
        function [V,D] = optimize(obj, vflow, goal, alpha)
            % T x (GxT) dim matrix. Vacc in period t constrained by
            % A*vacc = vflow where vflow is T dim vec
            A = kron(eye(obj.T),obj.share'); % Constrains daily vacc flow
            lb = zeros(obj.T*obj.G, 1);
            vflow_mat = reshape(repmat(vflow,1,obj.G)',obj.G*obj.T,1);
            b = vflow;
            % Initial condition given by random allocation
            [V,D] = knitro_nlp(@(x)(goal(x) + alpha*mean((x - vflow_mat).^2)),...
                vflow_mat,[],[],A,b,lb,[],[],[],[],'knitro.opt');
            V = reshape(V, size(obj.Vdot));
        end
    end
end


% Read in vaccination rates
function Vdot = calc_Vdot(group_names,vaccine_data,simul)
    Vdot = zeros(size(group_names,1),simul.T);
    for gg = 1:length(group_names)
        % Subset to this group
        this_group = vaccine_data(strcmp(vaccine_data.group,group_names{gg}),:);
        this_group = sortrows(this_group,'cumul_n_accepted_nationally');
        for tt = 1:simul.T
            % Isolate the row to look at
            which_row   = all([this_group.cumul_n_accepted_nationally-1e-6>tt*simul.vflow...
                             (this_group.cumul_n_accepted_nationally-this_group.n_accepted_nationally)<=tt*simul.vflow],2);

            if ~all(~which_row)
                % Calculate the additional share vaccinated on each day
                Vdot(gg,tt) = this_group.flow_vaccinated(which_row)*(simul.vflow/this_group.n_accepted_nationally(which_row)); 
            end 
        end
    end
end


% Solve for R=1 mitigation strategy
function thetasq = solve_theta_R_zero(S,I,shares,nu,beta,rho)
    StoE = shares'*(beta.*(rho*I).*S);
    ItoR = shares'*(nu*I);
    thetasq      = min(ItoR/StoE,1);
end
