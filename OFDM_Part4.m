%%  Samaa Hany Seif Elyazal            
%%  Wireless Communication, Intake 42  
%%  LTE, LAB 4                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializtion
clc
clear all

%% System parameters
% Actual parameter
EbNo_range=0:20;
N_OFDMsymbols=1000;
ModulationOrder=4;
N_Subcarriers=512;
N_DataSubcarriers=300;    % 300 Subcarriers corresponding to 5MHz BW
SamplingRate=7.68e6;
CP_length=4.7e-6;         %normal CP=4.7us, extended=16.67us
MaxDelaySpread=3*1e-6;
No_ofPaths=round(MaxDelaySpread*SamplingRate);
N_ofCPbits = int64(SamplingRate*CP_length);

%Derived Paramters
N_Bits=N_DataSubcarriers*log2(ModulationOrder);
%N_Bits=N_Bits*N_OFDMsymbols;

BER = [];
ber1=0;
%% OFDM Tx
for EbNo=EbNo_range %EbNo_range    
    for itteration = 0 : 1e3
        
        % Bit Stream Generation
        Bits = randi([0 1],N_Bits,1);
        
        % Symbol Mapper
        Symbols=qammod(Bits,ModulationOrder,'InputType','Bit');
        Tx_Pilot=Symbols(6:6:end,:);
               
        % Guard Subcarriers
        GuardSide=(N_Subcarriers-N_DataSubcarriers)/2;
        InputIFFT=[zeros(GuardSide,1);Symbols;zeros(GuardSide,1)];
        
        % IFFT
        OutputIFFT=ifft(InputIFFT);
        
        % CP Insertion
        [NoOfRows NoOfCols] = size(OutputIFFT);
        OFDMsymbols=[OutputIFFT((NoOfRows-N_ofCPbits+1):NoOfRows,:); OutputIFFT];
        
        % P/S and transmition
        txSig=reshape(OFDMsymbols,548,1);
        
        % AWGN
        realAWGN = randn(N_Subcarriers+N_ofCPbits,1);
        imagAWGN = i*randn(N_Subcarriers+N_ofCPbits,1);
        AWGN = realAWGN+imagAWGN;
        
        channel=(1/sqrt(2*No_ofPaths))*(randn(No_ofPaths,1)+1i*randn(No_ofPaths,1));
        tx=conv(txSig,channel);
        transmitted= tx(1:NoOfRows+N_ofCPbits);
        Eb = ((ModulationOrder-1)*2^2)/(6*log2(ModulationOrder)*N_Subcarriers);        
        No = Eb/(10^(EbNo/10));
        
        %% OFDM Rx
        % Receiving with AWGN
        rxSig=transmitted + sqrt(No/2).*AWGN;
        
        % CP removement
        CP_removement=rxSig(N_ofCPbits+1:NoOfRows+N_ofCPbits);
        H=fft(channel,512);
        FFT_Output=(fft(CP_removement)./H);
        
        % Zero remonement
        zeros_removement=FFT_Output(GuardSide+1:NoOfRows-GuardSide,:);
        
        % Channel Estimation
        Rx_Pilot=zeros_removement(6:6:end,:);
        Ch_estiation=Rx_Pilot./Tx_Pilot;
        Interpolation=interp1(6:6:300,Ch_estiation,1:300,'linear','extrap').'; 
        zeros_removement=zeros_removement./Interpolation;
        rxSymbols = reshape(zeros_removement,N_DataSubcarriers,1);
        
        rxBits = qamdemod(rxSymbols,ModulationOrder,'OutputType','Bit');
                
        % Bit Error Rate
        [ber berRatio] = biterr(Bits,rxBits);
        ber1= ber1+berRatio;
    end   
    BER = [BER ber1/itteration];
    ber1=0;    
end

semilogy(EbNo_range,BER)
hold on