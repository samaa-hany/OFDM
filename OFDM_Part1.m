%%  Samaa Hany Seif Elyazal            %%%%%%%%%%%%%%%%%
%%  Wireless Communication, Intake 42  %%%%%%%%%%%%%%%%%
%%  LTE, LAB 1                         %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializtion
clc
clear all

% System Paramters
EbNo_range=0:20;
N_OFDMsymbols=1e4;
ModulationOrder=64;
N_Subcarriers=512;
N_DataSubcarriers=300;    % 300 Subcarriers corresponding to 5MHz BW
SamplingRate=7.68e6;
CP_length=4.7e-6;         %normal CP=4.7us, extended=16.67us
N_ofCPbits = int64(SamplingRate*CP_length);


%Derived Paramters
N_Bits=N_DataSubcarriers*log2(ModulationOrder);
N_Bits=N_Bits*N_OFDMsymbols;


%% Transmitter Side
% Bit Stream Generation
T_Bits = randi([0 1],N_Bits,1);
  
% Symbol Mapper
Symbols=qammod(T_Bits,ModulationOrder,'InputType','Bit');
   
% Series to Parallel
SymbolsParallel=reshape(Symbols,N_DataSubcarriers,N_OFDMsymbols);
   
% Guard add
GuardSide=(N_Subcarriers-N_DataSubcarriers)/2;
InputIFFT=[zeros(GuardSide,N_OFDMsymbols);SymbolsParallel;zeros(GuardSide,N_OFDMsymbols)];

% IFFT
OutputIFFT=ifft(InputIFFT);
    
% CP Insertion
[NoOfRows NoOfCols] = size(OutputIFFT); 
OFDMsymbols=[OutputIFFT((NoOfRows-N_ofCPbits+1):NoOfRows,:); OutputIFFT];
    
% Parallel to Series  
txSig=reshape(OFDMsymbols,548*1e4,1);
[r c] = size(txSig);

%% AWGN
Noise = randn(r,1)+1i*randn(r,1);

BER = [];
for EbNo=EbNo_range 
    Eb = ((ModulationOrder-1)*2^2)/(6*log2(ModulationOrder)*N_Subcarriers);
    No = Eb/(10^(EbNo/10));
    
    
    %% Receiver Side
    Rceived_Sig=txSig+(sqrt(No/2)*Noise);
    
    R_Sig_Parallel = reshape(Rceived_Sig,N_Subcarriers+N_ofCPbits,N_OFDMsymbols);
    
    % CP Removement
    CP_removement=R_Sig_Parallel(N_ofCPbits+1:NoOfRows+N_ofCPbits,:);
    
    % FFT
    FFT_Output=fft(CP_removement);
    
    % Guard removement
    zeros_removement=FFT_Output(GuardSide+1:NoOfRows-GuardSide,:);
    R_Symbols = reshape(zeros_removement,N_DataSubcarriers*N_OFDMsymbols,1);
    
    % Symbol Demapper
    R_Bits = qamdemod(R_Symbols,ModulationOrder,'OutputType','Bit');
    
    % Bit Error Rate
    [ber berRatio] = biterr(R_Bits,T_Bits);
    
    BER = [BER berRatio];
    
end

semilogy(EbNo_range,BER)
hold on
title('Bit Error Rate of OFDM')
xlabel('SNR')
ylabel('BER')