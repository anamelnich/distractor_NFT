function ndf_main()

global stream ndf ID ids idm

% warning('off', 'all');
% Include any required toolboxes
ndf_include(); %adds paths to CNBI toolkit and eegc3
addpath(genpath('../decoder'));
addpath(genpath('../functions'));

skip_iterations = true;
% Prepare and enter main loop
try
    load('./decoder.mat');
    disp('Decoder Updated at');
    disp(decoder.datetime);

    ndf_initialization(); %sets up ndf configuration, should automatically setup ndf with 64 ch based on incoming data
    decoder = initializeParams(decoder);
    cleanupObj = onCleanup(@() ndf_down(decoder));

    tid_attach(ID);
    disp('[ndf] Receiving NDF frames...');

    %%
    update_flag = false;
    %% Main Loop %%
    while(true)
        tic
        [ndf.frame, ndf.size] = ndf_read(ndf.sink, ndf.conf, ndf.frame); %read data, outputs frame = data and count = data size
        time_frame = 1000*toc;

        % Acquisition is down, exit
        if(ndf.size == 0)
            disp('[ndf] Broken pipe');
            break;
        end

        if ((time_frame < 20) && (skip_iterations)) %skip iteration if reading data takes less than 20 msec 
            disp(['Skipping iteration. Time was: ' num2str(time_frame) ' ms']);
        else
            skip_iterations = false;

            eeg_input = ndf.frame.eeg; %(samples x eeg_channels)
            eog_input = ndf.frame.exg; %(samples x exg_channels)
            trigger_input = ndf.frame.tri; %(samples x tri_channels)

            ndf_store_signals([eeg_input, eog_input], trigger_input, decoder); %stores EEG and trigger data in stream, includes bandpass filter based on spatial filter in decoder, also has EOG filter (commented out)

            if (~any(isnan(stream.eeg(:))))
                first_index = find(ismember(stream.trigger, [102 104 100 110]), 1, 'first'); %returns sample (out of 512) where one of these triggers is found
                %disp(first_index)
                if (first_index >= round(0.2*decoder.fsamp)) & (first_index <= round(0.3*decoder.fsamp))
                    label_value = stream.trigger(first_index);
                    fprintf('Label value at first_index (%d): %d\n', first_index, label_value);
    
                    ex_posterior = singleClassification(decoder, stream.eeg((first_index - round(0.2*decoder.fsamp)):end, decoder.eegChannels), label_value, 1,decoder.leftElectrodeIndices,decoder.rightElectrodeIndices); %data from first index to end of buffer
                    %stream.eeg is 768 x 68, first_index = 154, label_value
                    %110
                    % eeg will rougly be 728x64
                    disp(['Time Frame: ' num2str(time_frame, '%.2f') ' Posteriors: ' num2str(ex_posterior, ' %.2f')]);
                    stream.trigger(first_index) = 0;
                    sendTiD(1 + (ex_posterior > decoder.decision_threshold)); % sends 1 if below threshold, 2 if above
                end
            end

            if (receiveTiD() == 20)
                break;
            end
        end
    end

catch exception
    ndf_printexception(exception);
end
end