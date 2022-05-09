%% K-nearest neighbor / Cross-fold validation step

% Import data as a matrix with features adn class labels from a .csv file
CellFeat = xlsread('CellFeatureData.xlsx');
% Values: [Area, Perimeter, Avg Intensity, Max Intensity, Nucleus Area,
% True Class Value]

% Create fold validation groups
Group = {};
TrueClass = {};
ind = [1, 51, 101, 201, 301, 351, 401];
count = [5, 5, 10, 10, 5, 5, 10];
for i = 1:10
    Group{i} = [];
    TrueClass{i} = [];
    for j = 1:length(ind)
        Group{i} = [Group{i}; CellFeat(ind(j):ind(j)+(count(j)-1), 1:5)];
        TrueClass{i} = [TrueClass{i}; CellFeat(ind(j):ind(j)+(count(j)-1), 6)];
        ind(j) = ind(j) + count(j);
    end
end

% 10-Fold Cross Validation: 1 data group = testing data, the other 9 will
% be the training data that you're comparing to. Conduct K-nn for all
% combinations and observe results.

% Set up the 10-fold cross validation test and training groups
TestData = {};
TrainData = {};
TrainClass = {};
for fold = 1:10
    TestData{fold} = [];
    TrainData{fold} = [];
    TrainClass{fold} = [];
    for i = 1:10
        if i == fold
            TestData{fold} = Group{i};           
        else
            TrainData{fold} = [TrainData{fold}; Group{i}];
            TrainClass{fold} = [TrainClass{fold}; TrueClass{i}];
        end
    end
end
%% K-Nearest Neighbor Euclidean Distance Detection
k = 2; 
acc = [];
for i = 1:10
    acc(i) = 0; % accuracy for each k value
    TestD = TestData{i};
    TrainD = TrainData{i};
    TrainC = TrainClass{i};
    TrueC = TrueClass{i};
    PredictC = zeros(1,50);
    for n = 1:length(TestD)
        dist= sqrt(sum((TrainD - TestD(n,:)).^2')); 
            % distance between Test and Train
        [dist2, ind] = sort(dist, 'ascend'); % find minimum distances
        TrainMaj = TrainC(ind); % resorting Training class values
        class_pos = TrainMaj(1:k); % find k minimum values
        PredictC(n) = mode(class_pos); % estimate mode(Test class)

        if (PredictC(n) == TrueC(n))
            acc(i) = acc(i) + 1;
        end
    end
    acc(i) = acc(i)/length(TestD);
end

avg_acc = mean(acc);




