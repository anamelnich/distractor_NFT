fisher_data = classEpochs(256:308,:,:);
data_mean = squeeze(mean(fisher_data,1));

%%
data_mean = data_mean';

%%
class1_data = data_mean(trainLabels ==1, :);
class0_data = data_mean(trainLabels ==0, :);

%%
mean_class1 = mean(class1_data,1);
mean_class0 = mean(class0_data,1);

var_class1 = var(class1_data,0,1);
var_class0 = var(class0_data,0,1);

%%
fisher_score = (mean_class1-mean_class0).^2./(var_class1+var_class0)