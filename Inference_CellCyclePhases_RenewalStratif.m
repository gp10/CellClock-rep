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

%% OTHER CELL KINETIC PARAMETERS (from Piedrafita et al, 2020, Nat. commun.)
r = 0.10; % symmetric division probability
dens = 0.65; % fraction of basal progenitor cells (others being basal resting cells -G0- transitioning towards differentiation)
gamma = 5.4; % stratification rate (week-1)
m = 1.25; % suprabasal-to-basal cell ratio
mu = 1.5; % shedding rate (week-1)

%% TIME SPENT IN DIFFERENT CELL-CYCLE PHASES:
% Assumption: Cells stay in S-, G2- and M- phases a fixed amount of time,
% given by the fraction of progenitor basal cells found in those phases times the
% average cell cycle period, while time in G1 phase is variable (will be
% shorter or longer depending on individual cell cycle period).

% Fractions in S-, G2- and M- have to be corrected given that not all basal cells are cycling progenitors
freq_bProg_WT = freq_bAll_WT./dens; % fractions referred exclusively to cycling progenitors (excluding resting basal cells that are differentiating)
tS_WT = freq_bProg_WT(2)/100 * tcc_avg; % (week)
tG2_WT = freq_bProg_WT(3)/100 * tcc_avg; % (week)
tM_WT = freq_bProg_WT(4)/100 * tcc_avg; % (week)

% Check time in S+G2+M is always smaller than any individual cell cycle period:
tSG2M_WT = tS_WT + tG2_WT + tM_WT % (week)
tSG2M_WT < tlag

%% SIMULATION OF CELL-CYCLE TRANSITIONS AND CELLULAR RENEWAL:

% Vector of times when to collect phase info:
rtime = [0:0.05:10]/7; % (week)

% Number of simulated basal clones (initiated as EdU-labelled single cells):
indiv = 100000;

% Call Monte-Carlo simulator of clone dynamics:
% (assumption that division waiting time - the period between consecutive
% divisions - and shedding waiting time - the period in suprabasal layer -
% follow realistic, delayed gamma distributions)
[nx_basal,nx_total, nx, X,Tini,Tend] = MCsimulator_dynamics_EdU_SP_total(rtime,dens,lambda,r,gamma,mu,m,indiv,tlag,GamShape,tS_WT,tG2_WT,tM_WT);


%% SUMMARY OF INFERRED CELL-CYCLE PHASE FREQUENCIES AT DIFFERENT TIME POINTS:

% Inferred frequencies:
ifreq_bEdU_WT = zeros(length(rtime),4);
preifreq_bEdU_WT(:,:) = sum(nx(:,:,1:5),1) ./ sum(sum(nx(:,:,1:5),1),3) .* 100;
ifreq_bEdU_WT = [preifreq_bEdU_WT(:,1)+preifreq_bEdU_WT(:,5) preifreq_bEdU_WT(:,2) preifreq_bEdU_WT(:,3) preifreq_bEdU_WT(:,4)];

% Inferred cumulative frequencies:
icumfreq_bEdU_WT = cumsum(ifreq_bEdU_WT,2);

% Inferred frequencies within batches limited to realistic sample sizes:
nx_sample = [];
No_batches = 100;
for iter = 1:No_batches % No. of batches
    rndId = randperm(indiv);
    for bas = 1:length(rtime)
        cont1 = 0;
        cont2 = 1;
        myx = [];
        % add to the list cells matching the given time frame within new randomly selected clones until reaching threshold n_FoV_WT:
        while cont1 < sum(n_FoV_WT)
            xlocs = find( (Tini{1,rndId(cont2)} <= rtime(1,bas)) & (rtime(1,bas) < Tend{1,rndId(cont2)}) );
            myx = [myx; X{1,rndId(cont2)}(xlocs,:)]; % add cells from that clone
            cont1 = sum(sum(myx(:,1:5),1));
            cont2 = cont2 + 1;
        end
        % convert to frequencies and pile up
        nx_sample(iter,bas,:) = sum(myx,1);
    end
end
ifreq_bEdU_WT_sample = zeros(No_batches, length(rtime), 4);
for aja = 1:No_batches
    preifreq_bEdU_WT_sample(:,:) = nx_sample(aja,:,:) ./ sum(nx_sample(aja,:,1:5),3) .* 100;    
    ifreq_bEdU_WT_sample(aja,:,:) = [preifreq_bEdU_WT_sample(:,1)+preifreq_bEdU_WT_sample(:,5) preifreq_bEdU_WT_sample(:,2) preifreq_bEdU_WT_sample(:,3) preifreq_bEdU_WT_sample(:,4)];
end

% Standard deviation of frequencies between batches:
ifreq_bEdU_WT_sd(:,:) = std(ifreq_bEdU_WT_sample, 0, 1);

%% PLOT CHANGES IN CELL-CYCLE PHASE FREQUENCIES OVER TIME:
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

%% SUMMARY OF INFERRED NUMBER OF EdU+ CELLS AT DIFFERENT TIME POINTS:

% Inferred number of basal & total EdU+ cells:
inum_bEdU_WT = zeros(1,length(rtime));
inum_tEdU_WT = zeros(1,length(rtime));
inum_bEdU_WT(1,:) = sum(nx_basal,1);
inum_tEdU_WT(1,:) = sum(nx_total,1);

% Inferred number of basal & total EdU+ cells within batches limited to a realistic number of clones per FoV:
sizefactor = inum_bEdU_WT(1,find(rtime >= 2./7, 1)) ./ indiv; % population of basal EdU+ cells grew from t0 to 2d, defining a growth factor
nClones_FoV_WT = round( sum(n_FoV_WT) / sizefactor ); % initial avg. density of clones estimate (from avg. number of basal EdU+ cells at 2d)
inum_bEdU_WT_sample = zeros(No_batches,length(rtime));
inum_tEdU_WT_sample = zeros(No_batches,length(rtime));
No_batches = 100;
for iter = 1:No_batches
    rndId = randperm(indiv);
    for bas = 1:length(rtime)
        myx = [];
        % add to the list cells matching the given time frame within {nClones_FoV_WT} randomly selected clones:
        for bes = 1:nClones_FoV_WT
            xlocs = find( (Tini{1,rndId(bes)} <= rtime(1,bas)) & (rtime(1,bas) < Tend{1,rndId(bes)}) );
            myx = [myx; X{1,rndId(bes)}(xlocs,:)]; % add cells from that clone
        end
        % convert to frequencies and pile up
        inum_bEdU_WT_sample(iter,bas) = sum(sum(myx(:,1:5),1),2);
        inum_tEdU_WT_sample(iter,bas) = sum(sum(myx,1),2);
    end
end

% Standard deviation of the number of basal & total EdU+ cells between batches:
inum_bEdU_WT_sd = zeros(1,length(rtime));
inum_tEdU_WT_sd = zeros(1,length(rtime));
inum_bEdU_WT_sd(1,:) = std(inum_bEdU_WT_sample, 0, 1);
inum_tEdU_WT_sd(1,:) = std(inum_tEdU_WT_sample, 0, 1);

% Standard deviation of the sb/total EdU+ cell ratio:
iratio_sbtEdU_WT_sd = std( (inum_tEdU_WT_sample-inum_bEdU_WT_sample)./inum_tEdU_WT_sample ,0, 1);

% Include (for comparison) condition without shedding (full accumulation of suprabasal cells):
mu = 0.000001; % no shedding approximation
[nx_basal_ref,nx_total_ref, nx_ref, X_ref,Tini_ref,Tend_ref] = MCMC_split_GAMdiv_GAMshed_EdU_SP_total(rtime,dens,lambda,r,gamma,mu,m,indiv,tlag,GamShape,tS_WT,tG2_WT,tM_WT);

%% PLOT CHANGES IN %EdU+ CELLS OVER TIME, AS A MEASURE OF RENEWAL RATE: 
figure(4)
hold on
plot(rtime.*7, sum(nx_total_ref,1) ./ indiv, 'k-.')
plot(rtime.*7, inum_tEdU_WT ./ indiv, 'b')
plot(rtime.*7, inum_bEdU_WT ./ indiv, 'r')
%errorbar(rtime.*7, inum_tEdU_WT ./ indiv, inum_tEdU_WT_sd ./ nClones_FoV_WT, 'b')
%errorbar(rtime.*7, inum_bEdU_WT ./ indiv, inum_bEdU_WT_sd ./ nClones_FoV_WT, 'r')
line([2 2],[0 4], 'LineStyle', '--', 'Color', 'k')
hold off
legend('total (no shedding)','total','basal')
xlabel('Time (days)')
ylabel('fraction of EdU+ cells (relative to t0)')
ylim([0 4])

%% PLOT CHANGES IN FRACTION OF SB EdU+ CELLS OVER TIME, AS A MEASURE OF STRATIFICATION RATE:
figure(5)
hold on
plot(rtime.*7, (sum(nx_total_ref,1)-sum(nx_basal_ref,1)) ./ sum(nx_total_ref,1), 'k-.')
plot(rtime.*7, (inum_tEdU_WT-inum_bEdU_WT) ./ inum_tEdU_WT, 'b')
%errorbar(rtime.*7, (inum_tEdU_WT-inum_bEdU_WT) ./ inum_tEdU_WT,  iratio_sbtEdU_WT_sd, 'b')
line([2 2],[0 1], 'LineStyle', '--', 'Color', 'k')
hold off
legend('ratio (no shedding)','ratio')
xlabel('Time (days)')
ylabel('sb EdU / total EdU cells')
ylim([0 1])
