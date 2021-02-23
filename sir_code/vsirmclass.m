% Class Definition: MG-VSIR-M Model
classdef vsirmclass
% Parameters and objects
    properties
        G                   % Number of groups
        group_names         % Group identifier
        share               % Share in each group -- vector size 1xG
        gamma               % 1/expected resolution time (size G)
        beta                % Contact parameter 1x1
        rho                 % Contact matrix GxG
        T                   % Number of simulation days
        S                   % Susceptible GxT
        I                   % Infected GxT
        R                   % Resolved GxT
        D                   % Dead GxT
        V                   % Vaccinated GxT
        ifr                 % IFR -- Gx1
        S0                  % Initial susceptible -- Gx1
        I0                  % Initial infected -- Gx1
        Rnaught             % Rnaught -- 1x1
        Vdot                % Flow vaccination rate -- G x T
        eff_Vdot            % Flow vaccination rate to susceptible -- G x T
        theta               % Mitigation parameters -- 1 x T
        vaccinate           % {0,1,2}, 0: no vaccines, 1: vaccinate only S, 2: vaccinate both S and R
        mitigate            % Boolean, mitigation to R=1 or not
        constV              % Boolean, constant vaccination rate
        state               % Which state is being used
        date                % Date the object was created (timestamp)
        TTHI                % Time till herd immunity
    end
    
    methods
        % Initialize class properties
        function [obj] = vsirmclass(data,vaccine_data,init_state,contact_matrix,simul,options)
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
            obj.date  = date;
                

            % Read the necessary parameters
            obj.G     = height(data); 
            obj.group_names = data.group; 
            obj.share = data.sh_group;
            obj.ifr   = data.c19_ifr_group;
            
            % Remove NaNs
            obj.ifr(isnan(obj.ifr)) = 0;
            
            if ~isfield(options,'R')
                obj.Rnaught     = 2.3;
            else
                obj.Rnaught     = options.Rnaught;
            end
            if ~isfield(options,'gamma')
                obj.gamma     = ones(size(obj.share))*1/10;
                options.gamma = 1/10;
            else
                obj.gamma = ones(size(obj.share))*options.gamma;
            end 
            
            % Contact matrix
            if isempty(contact_matrix)
                obj.rho    = repmat(obj.share',obj.G,1);
            else
                obj.rho    = contact_matrix;
            end

            % Set beta_g
            if ~isfield(options,'beta') & isempty(contact_matrix)
                obj.beta  = obj.Rnaught*obj.gamma;
            elseif isempty(contact_matrix)
                obj.beta  = options.beta*ones(obj.G,1);
            else
                % Calibrate beta to set the Rnaught value using
                % eigenvalue formaulation
                obj.beta   = obj.Rnaught/max(eig(contact_matrix*diag(1./obj.gamma)))*ones(obj.G,1);
            end

            if isempty(init_state)
                obj.S0     = ones(size(obj.share))-10^-6;
                obj.I0     = ones(size(obj.share))*10^-6;                
            else
                obj.S0     = init_state.s0;
                obj.I0     = init_state.i0;
            end

            obj.T      = simul.T;
            
            % Calculate Vdot unless
            if simul.constV
                obj.Vdot  = repmat(simul.v,1,simul.T);
            else
                obj.Vdot  = calc_Vdot(obj.group_names,vaccine_data,simul);
            end
            obj.constV = simul.constV;
            
            % Initialize
            obj.S        = zeros(obj.G,obj.T);
            obj.I        = zeros(obj.G,obj.T);
            obj.R        = zeros(obj.G,obj.T);
            obj.D        = zeros(obj.G,obj.T);
            obj.V        = zeros(obj.G,obj.T);
            obj.theta    = ones(obj.G,obj.T-1);
            obj.eff_Vdot = zeros(obj.G,obj.T);
            
            % Size checks
            assert(isequal(size(obj.share),size(obj.gamma)));
            assert(isequal(size(obj.share),[obj.G,1]));
            assert(isequal(size(obj.rho),[obj.G,obj.G]));
            assert(isequal(size(obj.I0),size(obj.gamma)));
            assert(isequal(size(obj.S0),size(obj.gamma)));
            assert(isequal(size(obj.Vdot),size(obj.S)));
        end
        
        % SIR Simulation backbone
        function obj = sir_sim(obj,vaccinate,mitigate)
            % Store init
            obj.S(:,1) = obj.S0;
            obj.I(:,1) = obj.I0;
            resolved0  = 1 - obj.S0 - obj.I0;
            obj.R(:,1) = resolved0.*(1-obj.ifr/100);
            obj.D(:,1) = resolved0.*obj.ifr/100;
            % Set vaccine init
            if vaccinate == 1
                obj.eff_Vdot = obj.Vdot;
            end
            
            % Loop
            for tt = 2:obj.T
                if mitigate == 1
                     thetasq = solve_theta_R_zero(obj.S(:,tt-1), ...
                                                  obj.rho,obj.I(:,tt-1), ...
                                                  obj.share, ...
                                                  obj.gamma, obj.beta);
                     obj.theta(:,tt-1) = sqrt(thetasq);
                else
                    thetasq = obj.theta(:,tt-1).^2;
                end
                % Change in infections, using group the interaction model
                Rdot        = obj.gamma.*obj.I(:,tt-1).*(1-obj.ifr/100);
                Ddot        = obj.gamma.*obj.I(:,tt-1).*obj.ifr/100;
                Idot        = thetasq.*obj.beta.*obj.S(:,tt-1).*(obj.rho*obj.I(:,tt-1))-(Rdot+Ddot);

                obj.I(:,tt) = obj.I(:,tt-1)+Idot;
                obj.D(:,tt) = obj.D(:,tt-1)+Ddot;
                
                % Change in susceptible, including
                % vaccination if desired
                % Model vaccinating both the susceptible and recovered
                % by vaccinating only a proportion of the susceptible pop
                if vaccinate == 2
                    healthy  = obj.S(:,tt-1) + obj.R(:,tt-1);
                    sus_frac = obj.S(:,tt-1) ./ healthy;
                    obj.eff_Vdot(:,tt-1) = obj.Vdot(:,tt-1).*sus_frac;
                % Model random distribution of vaccinnes
                % by proportional weighting
                elseif vaccinate == 3
                    healthy  = obj.S(:,tt-1) + obj.R(:,tt-1);
                    sus_frac = obj.S(:,tt-1) ./ healthy;
                    Vfrac    = healthy ./ (healthy' * obj.share);
                    obj.Vdot(:,tt-1)     = (obj.Vdot(:,tt-1)'*obj.share)*Vfrac;
                    obj.eff_Vdot(:,tt-1) = obj.Vdot(:,tt-1).*sus_frac;
                end
                Sdot        = -Idot - Rdot - obj.eff_Vdot(:,tt-1);
                Sdot        = max(Sdot,-obj.S(:,tt-1));
                obj.S(:,tt) = obj.S(:,tt-1) + Sdot;
                obj.R(:,tt) = obj.R(:,tt-1) + Rdot - (obj.Vdot(:,tt-1)-obj.eff_Vdot(:,tt-1));
            end
            
            % Compute vaccination rate in the final period
            if vaccinate == 3
                healthy  = obj.S(:,tt) + obj.R(:,tt);
                Vfrac    = (healthy.*obj.share) ./ (healthy' * obj.share);
                obj.Vdot(:,tt)     = (obj.Vdot(:,tt)'*obj.share)*Vfrac;
            end
            
            % Output type
            obj.vaccinate = vaccinate;
            obj.mitigate  = mitigate;
        end
        
        % Total incidence helper function
        function cases = compute_incidence(obj)
            cases = 1 - (obj.S'*obj.share + cumsum(obj.eff_Vdot'*obj.share));
        end
        
        % Total incidence objective function
        function cases = incidence_rate(obj, Vdot, mitigation)
            obj.Vdot = reshape(Vdot, size(obj.Vdot));
            obj = obj.sir_sim(2,mitigation);
            cases = 1 - (obj.S(:,end)'*obj.share + sum(obj.eff_Vdot'*obj.share));
        end
        
        % Total deaths objective function
        function deaths = death_rate(obj, Vdot, mitigation)
            obj.Vdot = reshape(Vdot, size(obj.Vdot));
            obj = obj.sir_sim(2,mitigation);
            deaths = obj.D(:,end)'*obj.share;
        end
        
        % Optimize using knitro
        function [V,D] = optimize(obj, vflow, goal, alpha)
            V0 = zeros(obj.G*obj.T,1);
            A  = kron(eye(obj.T),obj.share');
            lb = zeros(obj.T*obj.G, 1);
            [V,D] = knitro_nlp(@(x)(goal(x) + alpha*mean(diff(x).^2)),...
                V0,[],[],A,vflow,lb,[],[],[],[],'knitro.opt');
            V = reshape(V, size(obj.Vdot));
        end

        % Time to herd immunity
        function obj = calc_TTHI(obj)
        % Change in infections, using group the interaction model
            Rdot        = obj.gamma.*obj.I.*(1-obj.ifr/100);
            Ddot        = obj.gamma.*obj.I.*obj.ifr/100;
            Idot        = obj.beta.*obj.S.*(obj.rho*obj.I)-(Rdot+Ddot);
            avIdot      = obj.share'*Idot;
            index       = find(avIdot<0,1);
            if all(avIdot>=0)
                obj.TTHI = NaN;
            else
                obj.TTHI = index + avIdot(index)/(avIdot(index) - avIdot(index+1));
            end
        end
    end
end


% Read in vaccination rates
function Vdot = calc_Vdot(group_names,vaccine_data,simul)
    Vdot = zeros(size(group_names,1),simul.T);
    for gg = 1:length(group_names)
        % Subset to this group
        this_group = vaccine_data(strcmp(vaccine_data.group,group_names{gg}),:);
        this_group = sortrows(this_group,'cumul_n_offered_nationally');
        for tt = 1:simul.T
            % Isolate the row to look at
            which_row   = all([this_group.cumul_n_offered_nationally>tt*simul.vflow...
                             (this_group.cumul_n_offered_nationally-this_group.n_offered_nationally)<=tt*simul.vflow],2);

            if ~all(~which_row)
                % Calculate the additional share vaccinated on each day
                Vdot(gg,tt) = this_group.flow_vaccinated(which_row)*(simul.vflow/this_group.n_offered_nationally(which_row)); 
            end 
        end
    end
end


% Solve for R=1 mitigation strategy
function thetasq = solve_theta_R_zero(S,rho,I,shares,gamma,beta)
    recovery     = shares'*(gamma.*I)./beta;
    newinfect    = shares'*(S.*(rho*I));
    thetasq      = min(recovery/newinfect,1);
end
