function metrics = calcClassMetrics(behdata)
    % Overall classification metrics
    TP = sum((behdata.trial_type == 1) & (behdata.class == 1));
    FP = sum((behdata.trial_type == 0) & (behdata.class == 1));
    TN = sum((behdata.trial_type == 0) & (behdata.class == 0));
    FN = sum((behdata.trial_type == 1) & (behdata.class == 0));

    % Compute overall metrics (protect against division by zero using eps)
    TPR = TP / max((TP + FN), eps);  % Sensitivity/Recall
    TNR = TN / max((TN + FP), eps);  % Specificity
    accuracy = (TP + TN) / length(behdata.trial_type);
    precision = TP / max((TP + FP), eps);
    F1 = 2 * TP / max((2*TP + FP + FN), eps);

    % Store overall metrics
    metrics.overall.TP = TP;
    metrics.overall.FP = FP;
    metrics.overall.TN = TN;
    metrics.overall.FN = FN;
    metrics.overall.TPR = TPR;
    metrics.overall.TNR = TNR;
    metrics.overall.accuracy = accuracy;
    metrics.overall.precision = precision;
    metrics.overall.F1 = F1;

    % Run-by-run classification metrics stored as vectors
    uniqueRuns = unique(behdata.runNumber);
    nRuns = length(uniqueRuns);

    % Pre-allocate vectors for each metric
    TP_runs = zeros(nRuns,1);
    FP_runs = zeros(nRuns,1);
    TN_runs = zeros(nRuns,1);
    FN_runs = zeros(nRuns,1);
    TPR_runs = zeros(nRuns,1);
    TNR_runs = zeros(nRuns,1);
    accuracy_runs = zeros(nRuns,1);
    precision_runs = zeros(nRuns,1);
    F1_runs = zeros(nRuns,1);
    runNumbers = uniqueRuns;  % This will be your x-axis when plotting

    for r = 1:nRuns
        runIdx = (behdata.runNumber == uniqueRuns(r));
        
        % Calculate counts for the current run
        TP_run = sum((behdata.trial_type(runIdx) == 1) & (behdata.class(runIdx) == 1));
        FP_run = sum((behdata.trial_type(runIdx) == 0) & (behdata.class(runIdx) == 1));
        TN_run = sum((behdata.trial_type(runIdx) == 0) & (behdata.class(runIdx) == 0));
        FN_run = sum((behdata.trial_type(runIdx) == 1) & (behdata.class(runIdx) == 0));
        
        % Compute metrics for the current run (with protection against division by zero)
        TPR_run = TP_run / max((TP_run + FN_run), eps);
        TNR_run = TN_run / max((TN_run + FP_run), eps);
        accuracy_run = (TP_run + TN_run) / sum(runIdx);
        precision_run = TP_run / max((TP_run + FP_run), eps);
        F1_run = 2 * TP_run / max((2*TP_run + FP_run + FN_run), eps);
        
        % Store the results in the pre-allocated vectors
        TP_runs(r) = TP_run;
        FP_runs(r) = FP_run;
        TN_runs(r) = TN_run;
        FN_runs(r) = FN_run;
        TPR_runs(r) = TPR_run;
        TNR_runs(r) = TNR_run;
        accuracy_runs(r) = accuracy_run;
        precision_runs(r) = precision_run;
        F1_runs(r) = F1_run;
    end

    % Store run-by-run metrics as vectors in the structure for easy plotting later.
    metrics.by_run.runNumbers = runNumbers;
    metrics.by_run.TP = TP_runs;
    metrics.by_run.FP = FP_runs;
    metrics.by_run.TN = TN_runs;
    metrics.by_run.FN = FN_runs;
    metrics.by_run.TPR = TPR_runs;
    metrics.by_run.TNR = TNR_runs;
    metrics.by_run.accuracy = accuracy_runs;
    metrics.by_run.precision = precision_runs;
    metrics.by_run.F1 = F1_runs;
end
