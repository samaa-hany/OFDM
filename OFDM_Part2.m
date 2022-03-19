%%  Samaa Hany Seif Elyazal            %%%%%%%%%%%%%%%%%
%%  Wireless Communication, Intake 42  %%%%%%%%%%%%%%%%%
%%  LTE, LAB 1                         %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializtion
clc
clear all

% System Paramters
EbNo_range=0:20;
N_OFDMsymbols=1;
ModulationOrder=64;
N_Subcarriers=512;
N_DataSubcarriers=300;    % 300 Subcarriers corresponding to 5MHz BW
SamplingRate=7.68e6;
CP_length=4.7e-6;         %normal CP=4.7us, extended=16.67us
maxDelay=3e-6;
iteration=1e3;
ber1=0;
N_CP_Bits =round(SamplingRate*CP_length);


% Derived Paramters
N_Bits=N_DataSubcarriers*log2(ModulationOrder);
N_Bits=N_Bits*N_OFDMsymbols;
no_Path=round(maxDelay*SamplingRate);

BER = [];
for EbNo=EbNo_range 
    for x=1:iteration
    
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
        [Rows ,Cols] = size(OutputIFFT); 
        OFDMsymbols=[OutputIFFT((Rows-N_CP_Bits+1):Rows,:); OutputIFFT];
        h=sqrt(1/(2*no_Path))*(randn(no_Path,1)+1i*randn(no_Path,1));

        % Parallel to Series  
        Transmitted_Sig=reshape(OFDMsymbols,(Rows+N_CP_Bits),1);

        %% AWGN
        Eb = ((ModulationOrder-1)*2^2)/(6*log2(ModulationOrder)*N_Subcarriers);
        No = Eb/(10^(EbNo/10));


        %% Receiver Side
        Received_Sig=(conv(h,Transmitted_Sig));
        Received_Signal=Received_Sig(1:N_Subcarriers+N_CP_Bits,:);
        [R ,C] = size(Received_Signal);
        Noise = randn(R,1)+1i*randn(R,1);
        Received_Signal=(Received_Signal+sqrt(No/2)*Noise);

        % CP Removement
        CP_removement=Received_Signal(N_CP_Bits+1:end-N_CP_Bits,:);
    
        % FFT
        OutputFFT=fft(CP_removement,N_Subcarriers);

        % Equalization
        H=fft(h,N_Subcarriers);
        %H=H(:,N_CP_Bits+1:end-N_CP_Bits);
        R_Signal=OutputFFT./H;

        % Guard removement
        R_Symbols=R_Signal(GuardSide+1:Rows-GuardSide,:);

        % Symbol Demapper
        R_Bits1 = qamdemod(R_Symbols,ModulationOrder,'OutputType','Bit');
        [r ,c]=size(R_Bits1);
        R_Bits = reshape(R_Bits1,c*r,1);
        [ber , berRatio] = biterr(R_Bits,T_Bits);
        ber1=ber1+berRatio;

    end
    % Bit Error Rate
    BER = [BER ber1/iteration];
    ber1=0;
end

semilogy(EbNo_range,BER)
title('Bit Error Rate of OFDM')
xlabel('SNR')
ylabel('BER')