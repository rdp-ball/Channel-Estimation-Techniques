%%
k = 4;    % number of input bits
M = 2^k;  % number of possible input symbols
n = 7;    % number of channel uses
EbNo = 3; % Eb/No in dB

% Convert Eb/No to channel Eb/No values using the code rate
R = k/n;
EbNoChannel = EbNo + 10*log10(R);

wirelessAutoencoder = [
  featureInputLayer(M,"Name","One-hot input","Normalization","none")
  
  fullyConnectedLayer(M,"Name","fc_1")
  reluLayer("Name","relu_1")
  
  fullyConnectedLayer(n,"Name","fc_2")
  
  helperAEWNormalizationLayer("Method", "Energy", "Name", "wnorm")
  
  helperAEWAWGNLayer("Name","channel", ...
    "NoiseMethod","EbNo", ...
    "EbNo",EbNoChannel, ...
    "BitsPerSymbol",2, ... % channel use per channel symbol
    "SignalPower",1)
  
  fullyConnectedLayer(M,"Name","fc_3")
  reluLayer("Name","relu_2")
  
  fullyConnectedLayer(M,"Name","fc_4")
  softmaxLayer("Name","softmax")
  
  classificationLayer("Name","classoutput")]
%%
n = 2;                      % number of channel uses
k = 2;                      % number of input bits
EbNo = 3;                   % dB
normalization = "Energy";   % Normalization "Energy" | "Average power"

[txNet(1),rxNet(1),infoTemp,wirelessAutoEncoder(1)] = ...
  helperAEWTrainWirelessAutoencoder(n,k,normalization,EbNo);
infoTemp.n = n;
infoTemp.k = k;
infoTemp.EbNo = EbNo;
infoTemp.Normalization = normalization;
info = infoTemp;
%%
figure
helperAEWPlotTrainingPerformance(info(1))
%%
figure
tiledlayout(2,2)
nexttile([2 1])
plot(wirelessAutoEncoder(1))
title('Autoencoder')
nexttile
plot(txNet(1))
title('Encoder/Tx')
nexttile
plot(rxNet(1))
title('Decoder/Rx')
%%
simParams.EbNoVec = 0:0.5:8;
simParams.MinNumErrors = 10;
simParams.MaxNumFrames = 300;
simParams.NumSymbolsPerFrame = 10000;
simParams.SignalPower = 1;

%%
R = k/n;
EbNoChannelVec = simParams.EbNoVec + 10*log10(R);
M = 2^k;
txConst = comm.ConstellationDiagram(ShowReferenceConstellation=false, ...
  ShowLegend=true, ChannelNames={'Tx Constellation'});
rxConst = comm.ConstellationDiagram(ShowReferenceConstellation=false, ...
  ShowLegend=true, ChannelNames={'Rx Constellation'});
BLER = zeros(size(EbNoChannelVec));
%parfor trainingEbNoIdx = 1:length(EbNoChannelVec)

% Initialize arrays to store channel parameters
hEstimates = zeros(size(EbNoChannelVec)); % Assuming EbNoChannelVec is defined

for trainingEbNoIdx = 1:length(EbNoChannelVec)
  EbNo = EbNoChannelVec(trainingEbNoIdx);
  chan = comm.AWGNChannel("BitsPerSymbol",2, ...
    "EbNo", EbNo, "SamplesPerSymbol", 1, "SignalPower", 1);

  numBlockErrors = 0;
  frameCnt = 0;
  while (numBlockErrors < simParams.MinNumErrors) ...
      && (frameCnt < simParams.MaxNumFrames)

    d = randi([0 M-1],simParams.NumSymbolsPerFrame,1);    % Random information bits
    x = helperAEWEncode(d,txNet(1));                      % Encoder
    txConst(x)
    y = chan(x);                                          % Channel
    
    % Channel estimation (assuming function estimateChannel exists)
    hEst = estimateChannel(y, x);

    % Store the estimated channel gain
    hEstimates(trainingEbNoIdx) = hEst;
    rxConst(y)

    % Channel estimation
    hEst = estimateChannel(y, x);
    
    % Use the estimated channel to adjust the received signal
    yEst = y / hEst;

    dHat = helperAEWDecode(yEst,rxNet(1));                   % Decoder

    numBlockErrors = numBlockErrors + sum(d ~= dHat);
    frameCnt = frameCnt + 1;
  end
  BLER(trainingEbNoIdx) = numBlockErrors / (frameCnt*simParams.NumSymbolsPerFrame);
end

fprintf('Estimated Channel Gains:\n');
disp(hEstimates);

figure;
plot(EbNoChannelVec, hEstimates, '-o');
xlabel('Eb/No (dB)');
ylabel('Estimated Channel Gain');
title('Estimated Channel Gain vs. Eb/No');
grid on;

function hEst = estimateChannel(receivedSignal, transmittedSignal)
    % Simple channel estimation: least squares estimate of channel gain
    hEst = (transmittedSignal' * receivedSignal) / (transmittedSignal' * transmittedSignal);
end
