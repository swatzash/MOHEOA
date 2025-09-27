clear all
clc
close all
drawing_flag = 1;

TestProblem = 'UF2';
nVar = 10;

fobj = cec09(TestProblem);  % Assuming cec09 returns a multi-objective function

xrange = xboundary(TestProblem, nVar);

% Lower bound and upper bound
lb = xrange(:, 1)';
ub = xrange(:, 2)';

VarSize = [1 nVar];
Max_iter = 100;  % Maximum Number of Iterations
he_num = 100;
Archive_size = 50;   % Repository Size

N = 50;  % Number of search agents

grid_infl = 0.1;  % Grid Inflation Parameter
nGrid = 10;  % Number of Grids per each Dimension
leadsel = 4;  % Leader Selection Pressure Parameter
todel = 2;    % Extra (to be deleted) Repository Member Selection Pressure

humanevalution = CreateEmptyParticle(he_num);
for i = 1:he_num
    humanevalution(i).Velocity = 0;
    humanevalution(i).Position = zeros(1, nVar);
    for j = 1:nVar
        humanevalution(i).Position(1, j) = unifrnd(lb(j), ub(j), 1);
    end
    humanevalution(i).Cost = fobj(humanevalution(i).Position')';
    humanevalution(i).Best.Position = humanevalution(i).Position;
    humanevalution(i).Best.Cost = humanevalution(i).Cost;
end

humanevalution = DetermineDomination(humanevalution);

Archive = GetNonDominatedParticles(humanevalution);

Archive_costs = GetCosts(Archive);
G = CreateHypercubes(Archive_costs, nGrid, grid_infl);

for i = 1:numel(Archive)
    [Archive(i).GridIndex, Archive(i).GridSubIndex] = GetGridIndex(Archive(i), G);
end

A = 0.6;  % Warning value
LN = 0.4; % Percentage of leaders
EN = 0.4; % Percentage of explorers
FN = 0.1; % Percentage of followers

LNNumber = round(N * LN);  % Number of leaders
ENNumber = round(N * EN);  % Number of explorers
FNNumber = round(N * FN);  % Number of followers

jump_factor = abs(lb(1, 1) - ub(1, 1)) / 1000;

if (max(size(ub)) == 1)
    ub = ub .* ones(1, nVar);
    lb = lb .* ones(1, nVar);
end

% Initialization
X0 = initializationLogistic(N, nVar, ub, lb);  % Logistic chaotic mapping initialization
X = X0;


% Start search
for it = 1:Max_iter
    for i = 1:he_num
        clear rep2
        clear rep3
        leader = SelectLeader(Archive, leadsel);

        R = rand(1);

        for j = 1:size(X, 1)
            % Human exploration stage-------------------------------------------------------------------------------------
            if i <= (1 / 4) * Max_iter
               humanevalution(i).Position = leader.Position.* (1 - i / Max_iter) + (mean(humanevalution(j).Position) - leader.Position) * floor(rand() / jump_factor) * jump_factor + 0.2 * (1 - i / Max_iter) * (humanevalution(j).Position - leader.Position) .* Levy(nVar);  % Exploration stage
            else
                % Human development stage-------------------------------------------------------------------------------------
                for j = 1:LNNumber  % Leaders
                    if (R < A)
                       humanevalution(i).Position = 0.2 * cos(pi / 2 * (1 - (i / Max_iter))) * humanevalution(j).Position * exp((-i * randn(1)) / (rand(1) * Max_iter));
                    else
                        humanevalution(i).Position = 0.2 * cos(pi / 2 * (1 - (i / Max_iter))) * humanevalution(j).Position+ randn() * ones(1, nVar);
                    end
                end

                for j = LNNumber + 1:LNNumber + ENNumber  % Explorers
                    humanevalution(i).Position = randn() .* exp((humanevalution(N).Position - humanevalution(j).Position) / j^2);
                end

                for j = LNNumber + ENNumber + 1:LNNumber + ENNumber + FNNumber  % Followers
                   humanevalution(i).Position = humanevalution(j).Position + 0.2 * cos(pi / 2 * (1 - (i / Max_iter))) * rand(1, nVar) .* (humanevalution(1).Position - humanevalution(j).Position);  % Move with a randomly adaptive step size in the direction of the current best leader
                end

                for j = LNNumber + ENNumber + FNNumber + 1:N  % Losers
                   humanevalution(i).Position = leader.Position + (leader.Position - humanevalution(j).Position) * randn(1);
                end
            end
            humanevalution(i).Position=min(max(humanevalution(i).Position,lb),ub);

            humanevalution(i).Cost = fobj(humanevalution(i).Position')';
        end
    end

    humanevalution = DetermineDomination(humanevalution);
    non_dominated_soln = GetNonDominatedParticles(humanevalution);

    Archive = [Archive; non_dominated_soln];

    Archive = DetermineDomination(Archive);
    Archive = GetNonDominatedParticles(Archive);

    for i = 1:numel(Archive)
        [Archive(i).GridIndex, Archive(i).GridSubIndex] = GetGridIndex(Archive(i), G);
    end

    if numel(Archive) > Archive_size
        EXTRA = numel(Archive) - Archive_size;
        Archive = DeleteFromRep(Archive, EXTRA, todel);

        Archive_costs = GetCosts(Archive);
        G = CreateHypercubes(Archive_costs, nGrid, grid_infl);
    end

    disp(['In iteration ' num2str(it) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
    save results

    % Results

    costs = GetCosts(humanevalution);
    Archive_costs = GetCosts(Archive);

    if drawing_flag == 1
        hold off
        plot(costs(1, :), costs(2, :), 'k.');
        hold on
        plot(Archive_costs(1, :).*10, Archive_costs(2, :).*10, 'ro','markerfacecolor','b','markersize',8);
        legend('Successful solutions', 'Non-dominated solutions');
        drawnow
    end
end

% Levy search strategy
function o = Levy(d)
    beta = 1.5;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v).^(1 / beta);
    o = step;
end

% Logistic chaotic mapping initialization
function Positions = initializationLogistic(pop, dim, ub, lb)
    Boundary_no = length(ub);  % number of boundaries

    for i = 1:pop
        for j = 1:dim
            x0 = rand;
            a = 4;
            x = a * x0 * (1 - x0);
            if Boundary_no == 1
                Positions(i, j) = (ub - lb) * x + lb;
                if Positions(i, j) > ub
                    Positions(i, j) = ub;
                end
                if Positions(i, j) < lb
                    Positions(i, j) = lb;
                end
            else
                Positions(i, j) = (ub(j) - lb(j)) * x + lb(j);
                if Positions(i, j) > ub(j)
                    Positions(i, j) = ub(j);
                end
                if Positions(i, j) < lb(j)
                    Positions(i, j) = lb(j);
                end
            end
            x0 = x;
        end
    end
end