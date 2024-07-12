clear all;close all;
Path0='.\';
Interlace=1; %This value is 1 if image slice is arranged in interlaced fashion.
SaveData=1;
ExpNum=[45];
OmitFirst=0;Seg=1;BaseCorr=1;FreqVol=5480.58;NavPts=0;
Omit1=6;Omit2=6;%Omit 10 points at beginning and 10 points at the end of each k-line
OmitFID=0; %This value is normally set to 0, but we still keep it. 
Omit=0; %Omit k lines at the beginning and at the end of the 2D k-space

FOV=45;
directory=strcat(Path0,mat2str(ExpNum),'\pdata\1\');
reco=fopen(strcat(directory,'reco'));
tline = fgetl(reco);
tlines = cell(0,1);
while ischar(tline)
    tlines{end+1,1} = tline;
    tline = fgetl(reco);
end
fclose(reco);

Slice=char(tlines(12));
SliceChar=size(Slice);
SlicePar=str2num(Slice(29:SliceChar(2)));%Line 12 indicates slice number
NSlice=SlicePar;

method=fopen(strcat(Path0,mat2str(ExpNum),'\method'));
tline = fgetl(method);
tlines = cell(0,1);
while ischar(tline)
    tlines{end+1,1} = tline;
    tline = fgetl(method);
end
fclose(method);

fsLine=char(tlines(10));
fsLineSize=size(fsLine);
fsLine2=fsLine(15:20);
fs=str2num(fsLine2);

Mtx=char(tlines(75)); %Line 75 indicates matrix size
MtxSize=str2num(Mtx);
NRead0=MtxSize(1);
NPhase=MtxSize(2);

NRead=NRead0;

Rep=char(tlines(27)); %Line 27 indicates slice number
RepChar=size(Rep);
RepPar=str2num(Rep(21:RepChar(2)));
Rep=RepPar;

Seg=char(tlines(14)); %Line 14 indicates segment number
SegChar=size(Seg);
SegPar=str2num(Seg(14:SegChar(2)));
Seg=SegPar;

TR_String=char(tlines(15)); %Line 15 indicates TR
TRChar=size(TR_String);
TRPar=str2double(TR_String(23:TRChar(2)));
TR=TRPar/1000;

NCoil=2; %Double sample
SegSize=round(NPhase/Seg); 
FID_Start=1;FID_End=NRead0;K_Start=FID_Start;NReadAct=FID_End-FID_Start+1; %Actual number of acquisition points
fdemMatSqM=zeros(NReadAct,NCoil,NPhase/Seg,NSlice,Seg,Rep);
Image5D=zeros(NPhase,NReadAct-OmitFID,NCoil,NSlice,Rep);
SNR4D=zeros(NPhase,NReadAct-OmitFID,NSlice,Rep);
SNR1_4D=SNR4D;SNR2_4D=SNR4D;RmKSpike=1;

fdemMat=zeros(NReadAct*NPhase,NCoil);

fid1=fopen(strcat(Path0,mat2str(ExpNum(1)),'\fid'));
data=fread(fid1,'int32');
fclose(fid1);
data00=complex(data(1:2:end),data(2:2:end)); %(data(1:2:end) means index 1, 3, 5, ?, separate real and imaginary parts)
data0=reshape(data00,(NavPts+NRead0*NPhase*NCoil/Seg)*NSlice*Seg,Rep);
SizeData0=size(data0);
fdemWhole=zeros(SizeData0(1),Rep);

for RR=1:Rep
data1=data0(:,RR);
data1S=reshape(data1,(NavPts+NRead0*NPhase*NCoil/Seg),NSlice,Seg);
data2S=data1S(NavPts+1:NavPts+NRead0*NPhase*NCoil/Seg,:,:);
Coils=reshape(data2S,NRead0,NCoil,NPhase/Seg,NSlice,Seg);
phase00=angle(data1);
phase00=unwrap(phase00);
data1Size=size(data1);
data1SegSize=data1Size(1)/Seg;
fdem0=zeros(data1Size(1),1);
for SS=1:Seg
    fdemk=zeros(data1SegSize-1-OmitFirst,1);
    for LL=data1SegSize*(SS-1)+1+OmitFirst:data1SegSize*SS-1
        fdemk(LL-OmitFirst-data1SegSize*(SS-1))=(phase00(LL+FID_Start)-phase00(LL+FID_Start-1))/2/pi*fs;        
    end
    fdemk(LL-OmitFirst-data1SegSize*(SS-1)+1)=fdemk(LL-OmitFirst-data1SegSize*(SS-1));
    fdemWhole((data1SegSize-OmitFirst)*(SS-1)+1:SS*(data1SegSize-OmitFirst),RR)=fdemk;   
    if BaseCorr==1
        %fdemks=smooth(fdemk,8);
        fdemks=smooth(fdemk,8);
        fdem0(data1SegSize*(SS-1)+1:data1SegSize*SS)=fdemk-fdemks;
        fdemS=fdemks;        
    else
        fdem0(data1SegSize*(SS-1)+1:data1SegSize*SS)=fdemk;
    end
end
NavPts=0;NRead0=NRead;
fdem0S=reshape(fdem0,(NavPts+NRead0*NPhase*NCoil/Seg),NSlice,Seg);
fdem0S2=fdem0S(NavPts+1:NavPts+NRead0*NPhase*NCoil/Seg,:,:);
fdem01=reshape(fdem0S2,NRead0,NCoil,NPhase/Seg,NSlice,Seg);

data1S=reshape(data1,(NavPts+NRead0*NPhase*NCoil/Seg),NSlice,Seg);
data2S=data1S(NavPts+1:NavPts+NRead0*NPhase*NCoil/Seg,:,:);
Coils=reshape(data2S,NRead0,NCoil,NPhase/Seg,NSlice,Seg);
FID_Start=1;FID_End=NRead0;
    Coil1M=zeros(NReadAct,NPhase,NCoil,NSlice);
    fdemM=zeros(NReadAct,NPhase,NCoil,NSlice);

SS=1;

%-------------NSlice rearrange matrix---
for M=1:NSlice
for L=1:NCoil

Coil0=squeeze(Coils(FID_Start:FID_End,L,:,M,SS));
fdem02=squeeze(fdem01(FID_Start:FID_End,L,:,M,SS));
for k0=1:NPhase/Seg
    k=k0+(SS-1)*SegSize;
    Coil1M(:,ceil(k/SegSize)+(mod(k-1,SegSize))*Seg,L,M)=Coil0(:,k0); 
    fdemM(:,ceil(k/SegSize)+(mod(k-1,SegSize))*Seg,L,M)=fdem02(:,k0); 
    fdemMatSqM(Omit1+1:NRead-Omit2,L,k0,M,SS,RR)=fdem02(Omit1+1:NRead-Omit2,k0);
end

end
end %For NSlice
%-------------NSlice rearrange matrix

%%-------------NSlice
for M=1:NSlice %Modify this sentence to display the particular slice.
for L=1:NCoil %Segment number
Image2D=zeros(NReadAct-OmitFID,NPhase);
Image2D2=Image2D;
Image2D3=Image2D;
Image2D4=Image2D;

phase=zeros(NReadAct,NPhase);
phase2=phase;
phase0=zeros(NReadAct,NPhase);
phase02=zeros(NReadAct,NPhase);

Coil1=squeeze(Coil1M(:,:,L,M));
fdem=squeeze(fdemM(:,:,L,M));

for k=1:NPhase
phase(1:NReadAct,k)=angle(Coil1(1:NReadAct,k));
phase(1:NReadAct,k)=unwrap(phase(1:NReadAct,k));

phase0(1:NReadAct,k)=angle(Coil1(1:NReadAct,k));
phase0(1:NReadAct,k)=unwrap(phase0(1:NReadAct,k));
%------------------------------
    
    phase(K_Start:NReadAct,k)=detrend(phase(K_Start:NReadAct,k)); %h.	If there is linear phase variation, it means the osc freq is not exactly on center. That is why we need detrend
    phase(1:K_Start-1,k)=zeros(K_Start-1,1);
    
    phase02(:,k)=resample(phase0(1:NReadAct-2,k),NReadAct,NReadAct-2);
end    

%------Calculate initial 2D image without phase correction---------
for k=1:NPhase
fdem1=fdem(1:NReadAct-OmitFID,k);
if RmKSpike==1
fdem1(1:Omit1)=zeros(Omit1,1);
fdem1(NReadAct-OmitFID-Omit2+1:NReadAct-OmitFID)=zeros(Omit2,1);
end
sf=fft(fdem1);
sf=circshift(sf,round((NReadAct-OmitFID)/2));
Image2D(:,k)=sf(1:NReadAct-OmitFID);
end
Image2DF=fft((Image2D'));

%----Remove the spike along phase dimension-----------
Ave=mean(mean(abs(Image2DF)));
for LL=1:NPhase
    for MM=1:NReadAct-OmitFID
        if abs(Image2DF(LL,MM))>Ave*100
            Image2DF(LL,MM)=0;
        end
    end
end

Image2D_IF=ifft(Image2DF);
Image2D_IF=Image2D_IF';
%----------------------

for k=1+Omit:NPhase-Omit

    sf=Image2D_IF(:,k);
    sf=circshift(sf,-round((NRead-OmitFID)/2));
    fdem1=ifft(sf);
    fdem2=fdem1;
    for LL=1:NReadAct-OmitFID
        fdem1(LL)=fdem1(LL)*exp(-j*phase0(LL,k));             
    end
    sf0=fft(fdem(:,k));
    sf1=fft(fdem1);
    sf1=circshift(sf1,round((NReadAct-OmitFID)/2));%-round(IndArray(k)));

    Image2D3(:,k)=(sf1);
end
Image2DF_RM3=fft((Image2D3'));
Image2DF_RM3=circshift(Image2DF_RM3,[round(NPhase/2) 0]);
Image5D(:,:,L,M,RR)=abs(Image2DF_RM3);%+abs(Image2DF3);
fdemMat0=(fdem02(1:NReadAct,:));
fdemMat(NRead0*NPhase/Seg*(SS-1)+1:NRead0*NPhase/Seg*SS,L)=squeeze(fdemMat0(:));
end %For L=1:NCoil (2) for double sample
SNR1=abs(Image5D(:,:,1,M,RR));
SNR2=abs(Image5D(:,NReadAct-OmitFID:-1:1,2,M,RR));
SNR2=circshift(SNR2,[0 1 0]);
Diff=(SNR1-SNR2);
Noise_Amp=std2(Diff(10:30,10:30));
SNR1_4D(:,:,M,RR)=SNR1/Noise_Amp;
SNR2_4D(:,:,M,RR)=SNR2/Noise_Amp;
SNR=(SNR1+SNR2)/2;
SNR4D(:,:,M,RR)=SNR;
%%-------------NSlice
end %for M=1:NSlice
 end %RR repetition

 fdemWhole2=squeeze(fdemWhole(:));
 fdemWhole2_LP=lowpass(fdemWhole2,200,fs);

 %--Display the first repetition---
%Interlace=1;
 PlotImage=1;PlotRep=14;
if PlotImage==1
%for M=1:NSlice
for M=14:14
    if Interlace==0
        M2=M;
    else
        if M/2~=round(M/2)
            M2=ceil(M/2);
        else
            M2=ceil(NSlice/2)+M/2;
        end
    end
fig=figure;imagesc(abs(SNR4D(:,:,M2,PlotRep)));colormap(gray);axis square;title(strcat('SNR ',mat2str(M)));pbaspect([1 NPhase/NRead 1]);
end
end

if SaveData==1
save(strcat(Path0,'\EPIEEG',mat2str(ExpNum)),'fdemWhole2_LP','NRead','NPhase','NCoil','NSlice','Rep','TR','fs', 'Interlace','ExpNum','Path0');
save(strcat(Path0,'\EPIEEGS',mat2str(ExpNum)),'SNR4D');
end