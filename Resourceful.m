% Resourceful.m
% A 'game of life' where aging is an emergent property of resource limits.
% (c) Jack Peterson, 4/21/2012

% Setup:
% A single type of organism with a single resource type.
% Things live on an env_size x env_size grid.
% env is the environmental grid.
% haz is the hazards grid. Each time step, there is a probability haz(i,j)
% that an organism living at tile (i,j) will be killed.
% res is the resource grid. There is infinite resource available
% everywhere.
% pop is the population grid. Initial population is a single organism,
% located on a randomly selected tile.
% age is the age grid.
% rep is the replication grid. Initial replication probabilities are
% assigned from a uniform random distribution.
% fin is the final age of organisms once the simulation ends. If an
% organism has died on the site, then its age at death is reported. If
% multiple organisms have lived (sequentially) on the same site, then their
% median age is reported.
env_size = 10;
env = zeros(env_size);
haz = 0.01;
pop = zeros(env_size);
res = zeros(env_size);
[init_pop_row,init_pop_row] = max(rand(env_size,1));
[init_pop_col,init_pop_col] = max(rand(env_size,1));
pop(init_pop_row,init_pop_col) = 1;
max_pop = pop;
rep = rand(env_size);
median_age = zeros(env_size);
median_fin = zeros(env_size);

% Simulation
t_max = 100;
age = zeros(env_size,env_size) - 1;
bio = zeros(env_size,env_size);
fin = zeros(env_size,env_size) - 1;
bio(:,:,1) = double(~~pop);
total_pop = nan(t_max,1);
res = res + t_max + 1;
for tt  = 1:t_max
    % Every time step, each organism can either replicate or not.
    % If an organism replicates, it produces one new organism in a tile
    % either at the same position or adjacent to it.
    max_depth = max(max_pop(:));
    age(:,:,max_depth + 1) = zeros(env_size,env_size) - 1;
    fin(:,:,max_depth + 1) = zeros(env_size,env_size) - 1;
    bio(:,:,max_depth + 1) = zeros(env_size,env_size);
%     this_rep = (rand(env_size,env_size,max_depth+1) < repmat(rep,[1 1 max_depth+1])).*bio;
    this_rep = zeros(env_size,env_size,max_depth);
    for ii = 1:env_size
        for jj = 1:env_size
            for kk = 1:pop(ii,jj)
                this_rep(ii,jj,kk) = (rand(1) < rep(ii,jj))*bio(ii,jj,kk);
            end
        end
    end
    find_rep_row = cell(max_depth,1);
    find_rep_col = cell(max_depth,1);
    for jj = 1:max_depth
        [find_rep_row{jj},find_rep_col{jj}] = find(this_rep(:,:,jj));
        num_rep = length(find_rep_row{jj});
        for ii = 1:num_rep
            % Select an adjacent tile
            % Element 1 is row, element 2 is column
            select_tile = round(rand(2,1)*2) - 1;
            % Check for boundaries
            % If the selection is over the edge, then pick the other direction
            while select_tile(1) + find_rep_row{jj}(ii) > env_size
                select_tile(1) = round(rand(1)*2) - 1;
            end
            while select_tile(1) + find_rep_row{jj}(ii) <= 0
                select_tile(1) = round(rand(1)*2) - 1;
            end
            while select_tile(2) + find_rep_col{jj}(ii) > env_size
                select_tile(2) = round(rand(1)*2) - 1;
            end
            while select_tile(2) + find_rep_col{jj}(ii) <= 0
                select_tile(2) = round(rand(1)*2) - 1;
            end
            % Increment the population on the selected tile
            inc_row = find_rep_row{jj}(ii) + select_tile(1);
            inc_col = find_rep_col{jj}(ii) + select_tile(2);
            pop(inc_row,inc_col) = pop(inc_row,inc_col) + 1;
            max_pop(inc_row,inc_col) = max_pop(inc_row,inc_col) + 1;
            bio(inc_row,inc_col,max_pop(inc_row,inc_col)) = 1;
            age(inc_row,inc_col,max_pop(inc_row,inc_col)) = 0;
        end
    end
    % Each organism can be killed by a hazard
    depth = size(age);
    if numel(depth) == 2
        depth(3) = 1;
    end
    die = (rand(env_size,env_size,depth(3)) < haz).*bio;
    bio = bio - die;
    for ii = 1:env_size
        for jj = 1:env_size
            pop(ii,jj) = sum(bio(ii,jj,:));
            % If an organism died, record its age at death
            for kk = 1:depth(3)
                if die(ii,jj,kk)
                    fin(ii,jj,kk) = age(ii,jj,kk);
                    age(ii,jj,kk) = -1;
                end
            end
        end
    end
    % Increment the age of each surviving organism by 1
    age = (age + 2).*bio - 1;
    for ii = 1:env_size
        for jj = 1:env_size
            age_gt0 = age(ii,jj,:);
            age_gt0(age_gt0 < 0) = [];
            median_age(ii,jj) = median(age_gt0);
            fin_gt1 = fin(ii,jj,:);
            fin_gt1(fin_gt1 < 0) = [];
            median_fin(ii,jj) = median(fin_gt1);
        end
    end
    total_pop(tt) = sum(pop(:));
    if ~mod(tt,10)
        
        subplot(2,2,1)
        surf(pop)
        view(2)
        caxis([0 t_max/10])
        title(strcat('Population (t=',num2str(tt),')'))
        
        subplot(2,2,2)
        surf(median_age)
        view(2)
        caxis([0 t_max/10])
        title('Median age')

        subplot(2,2,3)
        plot(total_pop)
        xlabel('Time')
        ylabel('Population')

        subplot(2,2,4)
        present_age = age;
        present_age(present_age < 0) = [];
        if isempty(present_age)
            continue;
        end
        [y,x] = hist(present_age(:),0:max(present_age(:)));
        plot(x,y)
        xlabel('Age')
        ylabel('Count')
        
        pause(0.01)
    end
end