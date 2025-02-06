function decoder = initializeDecoder(params, decoder)

decoder.eegChannels = params.eegChannels; 
decoder.eogChannels = params.eogChannels;
decoder.spectralFilter = params.spectralFilter;
decoder.channelInterp = params.channelInterp;

decoder.triggerChannels = 1;
decoder.eogFilter = params.eogFilter;

decoder.asynchronous.waitSample = params.asynchronous.waitSample;

