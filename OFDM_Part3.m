%%  Samaa Hany Seif Elyazal            %%%%%%%%%%%%%%%%%
%%  Wireless Communication, Intake 42  %%%%%%%%%%%%%%%%%
%%  LTE, LAB 1                         %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializtion
clc
clear all

EbNo_range=100;
N_OFDMsymbols=138;
ModulationOrder=4;
N_subcarriers=512;
N_Datasubcarriers=300;
SamplingRate=7.68e6;
CP_length=4.7e-6;
u=[25;29;34]; %Root sequence

%Generate PSS
n1=0:30;
d11=reshape(exp((-1i.*pi.*n1.*(n1+1).*u(1))/63),[],1);
d21=reshape(exp((-1i.*pi.*n1.*(n1+1).*u(2))/63),[],1);
d31=reshape(exp((-1i.*pi.*n1.*(n1+1).*u(3))/63),[],1);
n2=31:61;
d12=reshape(exp((-1i.*pi.*(n2+1).*(n2+2).*u(1))/63),[],1);
d22=reshape(exp((-1i.*pi.*(n2+1).*(n2+2).*u(2))/63),[],1);
d32=reshape(exp((-1i.*pi.*(n2+1).*(n2+2).*u(3))/63),[],1);
d1=[d11;d12];
d2=[d21;d22];
d3=[d31;d32];
Zero_size=(N_subcarriers-length(d1))/2;
D1=[zeros(Zero_size,1);d1;zeros(Zero_size,1)];
D2=[zeros(Zero_size,1);d2;zeros(Zero_size,1)];
D3=[zeros(Zero_size,1);d3;zeros(Zero_size,1)];
N_Bits=N_Datasubcarriers*log2(ModulationOrder)*N_OFDMsymbols;
Ebno=100;
Eb=((ModulationOrder-1)*4)/(6*log2(ModulationOrder)*N_subcarriers);
No=Eb/db2pow(Ebno);
%%Generate Bits
Bits=randi([0 1],N_Bits,1);
%%Generate Symbols
Symbols=qammod(Bits,ModulationOrder,'InputType','bit','UnitAveragePower',true);
%%Serial TO parralel
Parralel_Symbol=reshape(Symbols,N_Datasubcarriers,N_OFDMsymbols);
%%Guard
GuardSide=zeros((N_subcarriers-N_Datasubcarriers)/2,N_OFDMsymbols);
SIG=[GuardSide;Parralel_Symbol;GuardSide];
SIG=[SIG(:,1:6),D2,SIG(:,7:75),D2,SIG(:,77:end)];
%%IFFT
sig=ifft(SIG);
%%CyclicPrefix
CP_samples=round(SamplingRate*CP_length);
sig_with_cp=[sig(N_subcarriers-CP_samples+1:N_subcarriers,:);sig];
%%serialization
sig_with_cp=reshape(sig_with_cp,1,size(sig_with_cp,1)*size(sig_with_cp,2));
%%AWGN
AWGN=(sqrt(No/2))*(randn(size(sig_with_cp))+1i*randn(size(sig_with_cp)));
Tx_signal=sig_with_cp+AWGN;
%% RECIEVER
%%MatchedFilter
Y1=fliplr(conj(ifft(D1)));
Y2=fliplr(conj(ifft(D2)));
Y3=fliplr(conj(ifft(D3)));
y1=conv(Y1,Tx_signal);
y2=conv(Y2,Tx_signal);
y3=conv(Y3,Tx_signal);
subplot(3,1,1)
plot(abs(y1));
subplot(3,1,2)
plot(abs(y2));
subplot(3,1,3)
plot(abs(y3));