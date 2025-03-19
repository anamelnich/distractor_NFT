function decoder = initializeParams(decoder)

global stream ndf

stream.fsamp = ndf.conf.sf;
stream.frame_size = ndf.conf.samples;
stream.cycle_freq = round(stream.fsamp / stream.frame_size);
stream.cycle_time = (stream.frame_size / stream.fsamp);
stream.num_channels = length([decoder.eegChannels decoder.eogChannels]);

stream.flags.initFilter = true;

% max_sample = 1.0*stream.fsamp;
max_sample = 1.5*stream.fsamp;

dummy_type = 2;
dummy_label = [102];

signalLength = ceil(max_sample/stream.frame_size)*stream.frame_size;
singleClassificationNew(decoder, rand(signalLength, length(decoder.eegChannels)), dummy_label, dummy_type, decoder.leftElectrodeIndices,decoder.rightElectrodeIndices);
stream.eeg = nan(signalLength, stream.num_channels);
stream.trigger = nan(signalLength, 1);
