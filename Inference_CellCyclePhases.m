%% LOAD DATA ON CELL-CYCLE PHASE FREQUENCIES:
freq_b = xlsread('freq_CellCyclePhases.xlsx','Sheet 1');
freq_bAll_Mut = freq_b(1:4,4)';
freq_bAll_WT = freq_b(5:8,4)';
freq_bEdU_Mut = freq_b(9:12,4)';
freq_bEdU_WT = freq_b(13:16,4)';

% averge EdU+ cells per FoV:
n_FoV_Mut = freq_b(9:12,5)';
n_FoV_WT = freq_b(13:16,5)';

%% PROBABILITY DISTRIBUTION FITTING INDIVIDUAL CELL-CYCLE PERIODS (from Piedrafita et al, 2020, Nat. commun.)
lambda = 2.9; % division rate (week-1)
tlag = 1/7; % refractory period (week)
GamShape = 4;

% GamScale (Gamma-dist scale param. that fits the observed average division rate - considering the tlag)
GamScale = (1/lambda - tlag) ./ GamShape;

% Check average cell cycle period:
tcc_avg = tlag + GamScale*GamShape % (week)
1/tcc_avg == lambda

% Distribution of cell cycle periods:
x = [0:0.01:1]';
y = gampdf(x,GamShape,GamScale);
figure(1)
plot((tlag+x)*7,y/sum(y))
ylabel('Probability')
xlabel('Time (days)')
xlim([0 10])

%% TIME SPENT IN DIFFERENT CELL-CYCLE PHASES:
% Assumption: Cells stay in S-, G2- and M- phases a fixed amount of time,
% given by the fraction of basal cells found in those phases times the
% average cell cycle period, while time in G1 phase is variable (will be
% shorter or longer depending on individual cell cycle period).
tS_WT = freq_bAll_WT(2)/100 * tcc_avg; % (week)
tG2_WT = freq_bAll_WT(3)/100 * tcc_avg; % (week)
tM_WT = freq_bAll_WT(4)/100 * tcc_avg; % (week)

% Check time in S+G2+M is always smaller than any individual cell cycle period:
tSG2M_WT = tS_WT + tG2_WT + tM_WT % (week)
tSG2M_WT < tlag

%% SIMULATION OF INDIVIDUAL CELL-CYCLE PERIODS & TRANSITIONS THROUGH CELL-CYCLE PHASES:

% Vector of times when to collect phase info:
rtime = [0:0.05:10]/7; % (week)

% Number of simulated basal cells:
indiv = 100000;
nx = zeros(indiv,length(rtime),4); % will store phase info

for it = 1:indiv
    
    it
    
    % Initial condition:
    time = 0;
    multiplo = 1;
    %   G1 S G2 M
    x = [0 1 0 0]; % cell in S-phase (EdU labelled)
    S_ini = 1;
    
    while time <= rtime(end)

        % Calculate time to new phase:
        if S_ini == 1 % initial S-phase progression (random, up to completion)
            t_next = time + rand * tS_WT; % initially, cells are randomly picked at different stages of completion of the S-phase
            x_next = [0 0 1 0]; % next is G2
            %disp([sprintf('%.2f',time), ' --- initial S-phase --- ', sprintf('%.2f',t_next)])
            S_ini = 0;
        else
            if x(2) == 1 % S-phase progression (fixed)
                t_next = time + tS_WT;
                x_next = [0 0 1 0]; % next is G2
                %disp([sprintf('%.2f',time), ' --- S-phase --- ', sprintf('%.2f',t_next)])
            elseif x(3) == 1 % G2-phase progression (fixed)
                t_next = time + tG2_WT;
                x_next = [0 0 0 1]; % next is M
                %disp([sprintf('%.2f',time), ' --- G2-phase --- ', sprintf('%.2f',t_next)])
            elseif x(4) == 1 % M-phase progression (fixed)
                t_next = time + tM_WT;
                x_next = [1 0 0 0]; % next is G1
                %disp([sprintf('%.2f',time), ' --- M-phase --- ', sprintf('%.2f',t_next)])
            elseif x(1) == 1 % G1-phase progression (random)
                tcc = gamrnd(GamShape, GamScale) + tlag; % random total cell-cycle period
                t_next = time + (tcc - tSG2M_WT);
                x_next = [0 1 0 0]; % next is S
                %disp([sprintf('%.2f',time), ' --- G1-phase --- ', sprintf('%.2f',t_next)])
            end
        end

        if t_next < rtime(1,multiplo)
            time = t_next;
            x = x_next;
        else
            while (multiplo < length(rtime)) && (t_next >= rtime(1,multiplo))
                nx(it,multiplo,:) = x;
                multiplo = multiplo + 1;
            end
            if (multiplo == length(rtime)) && (t_next >= rtime(1,multiplo))
                nx(it,multiplo,:) = x;
            end
            time = t_next;
            x = x_next;
        end

    end

end

%% SUMMARY OF INFERRED CELL-CYCLE PHASE FREQUENCIES AT DIFFERENT TIME POINTS:

% Inferred frequencies:
ifreq_bEdU_WT = zeros(length(rtime),4);
ifreq_bEdU_WT(:,:) = sum(nx,1) ./ indiv .* 100;

% Inferred cumulative frequencies:
icumfreq_bEdU_WT = cumsum(ifreq_bEdU_WT,2);

% Inferred frequencies within batches limited to realistic sample sizes:
nx_sample = [];
No_batches = 100;
for iter = 1:No_batches % No. of batches
    nx_sample(iter,:,:,:) = nx( randperm(indiv, sum(n_FoV_WT)), :,:);
end
ifreq_bEdU_WT_sample = zeros(No_batches, length(rtime), 4);
ifreq_bEdU_WT_sample(:,:,:) = sum(nx_sample,2) ./ sum(n_FoV_WT) .* 100;

% Standard deviation of frequencies between batches:
ifreq_bEdU_WT_sd(:,:) = std(ifreq_bEdU_WT_sample, 0, 1);

%% Plot changes over time:
figure(2)
mycol = [0.72 0.27 1;... %green
     0.13 0.65 0.73;... %cream
    0.47 0.67 0.19;...
    0.98 0.39 0]; %red
hold on
fill([rtime, fliplr(rtime)].*7, [icumfreq_bEdU_WT(:,1); flipud(zeros(length(rtime),1))], mycol(1,:), 'EdgeColor', 'none')
for iter = 2:4
    fill([rtime, fliplr(rtime)].*7, [icumfreq_bEdU_WT(:,iter); flipud(icumfreq_bEdU_WT(:,iter-1))], mycol(iter,:), 'EdgeColor', 'none')
end
line([2 2],[0 100], 'LineStyle', '--', 'Color', 'k')
line([7 7],[0 100], 'LineStyle', '--', 'Color', 'k')
hold off
xlim([0 10]); ylim([0 100])
yticks([0:25:100])
xlabel('Time post-EdU labelling (days)')
ylabel('Frequency from EdU+ cells')
legend({'G1','S','G2','M'}, 'Location', 'northoutside', 'Orientation', 'horizontal'); legend('boxoff')

% Retrieve statistics at given times:
mytime = [2 7]./7; % (weeks)
myifreq_avg = [];
myifreq_sd = [];
for iter = 1:length(mytime)
    loc = find(rtime >= mytime(iter), 1);
    myifreq_avg(iter,:) = ifreq_bEdU_WT(loc,:);
    myifreq_sd(iter,:) = ifreq_bEdU_WT_sd(loc,:);
end

% Plot at given times with errorbars:
figure(3)
%XX = categorical({'2d' , '7d'}); % group names
b = bar(myifreq_avg,'stacked');
[ngroups,nbars] = size(myifreq_avg);
x = nan(nbars, ngroups);
y = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
    y(i,:) = b(i).YEndPoints; % for the data values
end
hold on
errorbar(x',y',myifreq_sd,'k','linestyle','none','HandleVisibility','off');  
legend({'G1', 'S', 'G2', 'M'}, 'Location', 'eastoutside');
box off
xticklabels({'2d','7d'})
ylim([0 110]); yticks([0:25:100])
ylabel('Frequency from EdU+ cells')
