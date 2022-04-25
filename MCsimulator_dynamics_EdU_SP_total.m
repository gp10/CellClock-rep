function [nx_basal,nx_total, nx, X,Tini,Tend] = MCsimulator_dynamics_EdU_SP_total(rtime,dens,lambda,r,gamma,mu,m,indiv,tlag,GamShape,tS,tG2,tM)
%% Non-Markovian Monte-Carlo simulator of Single Progenitor (SP) model dynamics
% Considers a gamma distribution for both the division rate and the shedding
% rate (time spent in suprabasal layers) and an uneven initial induction
% condition where all cells are proliferating and quasi-synchronized (at S-phase).

%% Initial definition of parameters:
tic
% General initial parameters:
if (nargin < 8)
    indiv=10000;
    tlag = 0;
    GamShape = 4;
    tS = 0.05;
    tG2 = 0.02;
    tM = 0.005;
end
tSG2M = tS + tG2 + tM;
timelim=rtime(1,end); %round(365/7); % weeks
nx=zeros(indiv,length(rtime),6);

nx_basal=zeros(indiv,length(rtime));
nx_total=zeros(indiv,length(rtime));
nxss=zeros(indiv,6);

% GamScale (Gamma-dist scale param. that fits the observed average division rate - considering the tlag)
GamScale = (1/lambda - tlag) ./ GamShape;
% GamScale (Gamma-dist scale param. that fits the observed average shedding rate - considering the tlag)
tlagShed = tlag*10;
GamShapeShed = GamShape;
GamScaleShed = (1/mu - tlagShed) ./ GamShapeShed;

%% ITERATION FOR DIFFERENT INDIVIDUAL CLONES
for it=1:indiv

    % Initial cell content:
    %    - Prolif -    D  SB
    %     G1 S  G2 M   G0 SB
    x = [ 0  1  0  0   0  0]; % only single-cell clones with just one proliferating cell in S-phase (EdU labelled)

    % Initial variables & cell attributes:
    Id_Cell = 1;
    tini = 0;
    tend = 0;
    % Synchronous initial cond given by EdU staining: all cells in S-phase (defining a random value for the 'time-to-next-division' of initial cells at time 0)

    % ITERATION FOR EACH SINGLE CELL:
    while Id_Cell <= size(x,1)

        if tini(Id_Cell,1) <= timelim
            
            % Calculate time to finish current phase:
            if (x(Id_Cell,2)~=0) && (Id_Cell==1) % initial S-phase progression (random, up to completion)
                tend(Id_Cell,1) = tini(Id_Cell,1) + rand * tS; % initially, cells are randomly picked at different stages of completion of the S-phase
            else
                if x(Id_Cell,2)~=0 % S-phase progression (fixed)
                    tend(Id_Cell,1) = tini(Id_Cell,1) + tS;
                elseif x(Id_Cell,3)~=0 % G2-phase progression (fixed)
                    tend(Id_Cell,1) = tini(Id_Cell,1) + tG2;
                elseif x(Id_Cell,4)~=0 % M-phase progression (fixed)
                    tend(Id_Cell,1) = tini(Id_Cell,1) + tM;
                elseif x(Id_Cell,1)~=0 % G1-phase progression (random)
                    tend(Id_Cell,1) = tini(Id_Cell,1) + gamrnd(GamShape,GamScale) + tlag - tSG2M; % random total cell-cycle period minus fixed portion in S-G2-M
                elseif x(Id_Cell,5)~=0 % Differentiating (D) cell progression (random)
                    tau = -(1./gamma)*log(rand); % exponential random with avg rate gamma
                    tend(Id_Cell,1) = tini(Id_Cell,1) + tau;
                else % Suprabasal (SB) cell progression (random)
                    %tau = -(1./mu)*log(rand); % exponential random with avg rate mu
                    %tend(Id_Cell,1) = tini(Id_Cell,1) + tau;
                    tend(Id_Cell,1) = tini(Id_Cell,1) + gamrnd(GamShapeShed,GamScaleShed) + tlagShed;
                end
            end
            
            if tend(Id_Cell,1) < timelim
                % Update to new phase-state:
                myphase = find(x(Id_Cell,:)~=0,1);
                if (myphase==2)
                    x = [x; 0 0 1 0 0 0]; % S -> G2
                    tini = [tini; tend(Id_Cell,1)];
                elseif (myphase==3)
                    x = [x; 0 0 0 1 0 0]; % G2 -> M
                    tini = [tini; tend(Id_Cell,1)];
                elseif (myphase==4)
                    % Calculation of single-event probabilities:
                    p(1) = lambda*r; % symetric division
                    p(2) = lambda*(1-2*r); % asymetric division
                    p(3) = lambda*r; % symetric division
                    % Calculation of total probability of event:
                    pt=sum(p);
            
                    % Random number generator:
                    r2n=rand*pt;
                    
                    % Event selection:
                    event = find(cumsum(p)>r2n,1);
                    if (event==1)
                        x = [x; 1 0 0 0 0 0; 1 0 0 0 0 0]; % M -> G1 + G1 (both proliferating)
                        tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                    elseif (event==2)
                        x = [x; 1 0 0 0 0 0; 0 0 0 0 1 0]; % M -> G1 + G0 (proliferating + differentiating)
                        tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                    elseif (event==3)
                        x = [x; 0 0 0 0 1 0; 0 0 0 0 1 0]; % M -> G0 + G0 (both differentiating)
                        tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                    end
                
                elseif (myphase==1)
                    x = [x; 0 1 0 0 0 0]; % G1 -> S
                    tini = [tini; tend(Id_Cell,1)];
                elseif (myphase==5)
                    x = [x; 0 0 0 0 0 1]; % G0 -> SB
                    tini = [tini; tend(Id_Cell,1)];
                %elseif (myphase==6) | loss is implicit in the time update
                %x(6)=x(6)-1;
                end
            end

        end
        Id_Cell = Id_Cell + 1;

    end

    % Save individual cell properties:
    X{1,it} = x;
    Tini{1,it} = tini;
    Tend{1,it} = tend;

    % Save the populations of cells at certain time points:
    for bas = 1:size(rtime,2)
        nx(it,bas,1) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,1)==1 ),1);
        nx(it,bas,2) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,2)==1 ),1);
        nx(it,bas,3) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,3)==1 ),1);
        nx(it,bas,4) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,4)==1 ),1);
        nx(it,bas,5) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,5)==1 ),1);
        nx(it,bas,6) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,6)==1 ),1);
    end

    % Final values of the variables:
    nxss(it,:) = nx(it,end,:);

end

%% Sum all types of basal cells to get basal-layer clone sizes:
nx_basal = nx(:,:,1)+nx(:,:,2)+nx(:,:,3)+nx(:,:,4)+nx(:,:,5);
% Sum all types of cells to get total clone sizes:
nx_total = nx_basal + nx(:,:,6);
toc
