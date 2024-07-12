clear all;close all;
Path0='.\';
WithMedFilt1=0;
WithMedFilt2=0;
WithMedFilt3=1;
ExpNum=[45];
load(strcat(Path0,'EPIEEG',mat2str(ExpNum(1))));
FreqVol=5480.58;

SecondRemove=1; %Do subtraction for the second time? 
SaveSeparateData=1; %Separately save data for MRI and EEG. 
UseHighPass=1; %After the second removal of baseline, use high pass for the overall data. To reduce computaion time, set this value to 0.
HighPassEach=0; %Normally, no need to change this to 1.

%---In the following, we extract EEG---
SaveFig=0;
PlotPart=0;
if WithMedFilt1==1
GradTraceS=medfilt1(fdemWhole2_LP,NRead*NCoil*2);
else
GradTraceS=(fdemWhole2_LP);
end
    
figure;plot(GradTraceS);title('GradTraceS medfilt1 fdemWhole2 LP');
%GradTraceS2=smooth(GradTraceS,129*2*43);
%GradTraceD=GradTraceS-GradTraceS2;
%figure;plot(GradTraceD);title('GradTrace-TradTraceS');

%GradTraceDS=smooth(GradTraceD,129*43*2);
%figure;plot(GradTraceDS);
SlicePerExcite=1; %The number of slices per excitation pulse
%SlicePerExcite=4; %The number of slices per excitation pulse 06/19/23
GroupNum=NSlice*Rep/SlicePerExcite;
Grad=reshape(GradTraceS,NRead*NCoil*NPhase*SlicePerExcite,GroupNum);
GradSize=size(Grad);
Grad2=zeros(GradSize(1),GradSize(2));
Grad20=Grad2;

%---Remove baseline based on average the curved baseline---
NAve=5; %The number of repetitions to be averaged
Span=NSlice*NAve/2;
Ratio=zeros(GroupNum,1);
GradSumArr=Grad2;
for K=1:GroupNum

    %---First, calculate the average baseline within the Span---
    K1=K-Span/2;
    K2=K+Span/2;
    if K1<=0
        K1=1;
        K2=1+Span;
    end
    if K2>GroupNum
        K1=GroupNum-Span;
        K2=GroupNum;
    end
    GradSum=0;
    for MM=K1:K2
        GradSum=GradSum+Grad(:,MM);
    end
    GradSum=GradSum/max(K2-K1+1);
    %---Finish calculating the average baseline---
    %Ratio(K,1)=max(abs(Grad(:,K)))/max(abs(GradSum));
    %Ratio(K,1)=max(abs(Grad(1:5000,K)))/max(abs(GradSum(1:5000)));
    %Ratio(K,1)=mean(Grad(:,K)./GradSum(:,1));
    GradSumArr(:,K)=GradSum;
%    if abs(Ratio(K,1))<0.995
%        Grad2(:,K)=Grad(:,K)-GradSum*Ratio(K,1);
%    else      
        Grad2(:,K)=Grad(:,K)-GradSum;
%    end
    if HighPassEach==1      
        Grad20(:,K)=Grad2(:,K);
        Grad2(:,K)=highpass(Grad20(:,K),NSlice/TR,fs);
    end
end    
Grad2Size=size(Grad2);
TimeInd2=[1:Grad2Size(1)*Grad2Size(2)]/(Grad2Size(1)*Grad2Size(2))*TR*Rep;
Grad2Sq=squeeze(Grad2(:));
if PlotPart==1
    figure;plot(TimeInd2(1:5000000),Grad2Sq(1:5000000));title('Grad2Sq baseline removal');
else
    figure;plot(TimeInd2,Grad2Sq);title('Grad2 baseline removal');
end
if WithMedFilt2==1
    %Grad2SqM=medfilt1(Grad2Sq,NRead*NCoil);
    Grad2SqM=medfilt1(Grad2Sq,NRead);   
else
    Grad2SqM=Grad2Sq;
end    
if PlotPart==1
    figure;plot(TimeInd2(1:5000000),Grad2SqM(1:5000000));title('Grad2 baseline removal after medfilt1');
else
    figure;plot(TimeInd2,Grad2SqM);title('Grad2 baseline removal after medfilt1');
end

%---Separately remove the negative and positive spikes in second baseline removal---
if SecondRemove==0
    Grad4=reshape(Grad2SqM,NRead*NCoil*NPhase*SlicePerExcite,GroupNum);
else

Grad3=reshape(Grad2SqM,NRead*NCoil*NPhase*SlicePerExcite,GroupNum);
%Grad4=zeros(GradSize(1),GradSize(2));
Grad4=Grad3;
PosPeaks=0;
PosMean=0;
NegPeaks=0;
NegMean=0;
GroupNum1=GroupNum;
%---Calculate the mean of positive peaks
PosThres=250;
NegThres=250;
for K=1:GroupNum1
%for K=1:6
    [Peak,Ind]=max(abs(Grad3(100:round(GradSize(1)/15),K)));
    Grad3Ave=mean(Grad3(GradSize(1)/2:GradSize(1),K));
    Peak=Grad3(100+Ind-1,K); %Confirm this is a long-tail peak, rather than short burst
    if abs(Peak)>3.5*abs(Grad3Ave) 
        if Peak>PosThres
        PosPeaks=PosPeaks+1;
        PosMean=PosMean+Grad3(:,K);
        elseif Peak<-NegThres
        NegPeaks=NegPeaks+1;
        NegMean=NegMean+Grad3(:,K);       
        end
    end
end    
if PosPeaks>0
PosMean=PosMean/PosPeaks;
PosMeanMax=max(abs(PosMean(100:round(GradSize(1)/15))));
end
if NegPeaks>0
NegMean=NegMean/NegPeaks;
NegMeanMax=-max(abs(NegMean(100:round(GradSize(1)/15))));
end
%---Finish calculating positive peaks---
%---Subtract positive peaks with the mean of positive---

for K=1:GroupNum1
%for K=1:6
    [Peak,Ind]=max(abs(Grad3(100:round(GradSize(1)/15),K)));
    Grad3Ave=mean(Grad3(GradSize(1)/2:GradSize(1),K));
    Peak=Grad3(100+Ind-1,K);    
    if abs(Peak)>3.5*abs(Grad3Ave) 
        if Peak>PosThres
%        Grad4(:,K)=Grad3(:,K)-PosMean;
        Ratio(K)=Peak/PosMeanMax;
        Grad4(:,K)=Grad3(:,K)-PosMean*Ratio(K);
        elseif Peak<-NegThres
        Ratio(K)=Peak/NegMeanMax;
        Grad4(:,K)=Grad3(:,K)-NegMean*Ratio(K);            
        end
    end
end 
%---Finish subtraction---
PosPeaks
figure;plot(PosMean);title('PosMean');
NegPeaks
figure;plot(NegMean);title('NegMean');

if PlotPart==1
    figure;plot(TimeInd2(1:5000000),Grad4(1:5000000));title('Grad4 after second baseline removal of Grad2');
else
    figure;plot(TimeInd2,squeeze(Grad4(:)));title('Grad4 after second baseline removal of Grad2');
end  
Grad40=Grad4;
Grad4Sq=squeeze(Grad40(:));
if WithMedFilt3==1
    Grad4SqM=medfilt1(Grad4Sq,NRead*NCoil);
%    Grad4SqM=medfilt1(Grad2Sq,NRead);       
%Grad4SqM=medfilt1(Grad4Sq);
else
    Grad4SqM=Grad4Sq;
end    

%---High pass filter to flatten baseline.
if UseHighPass==1
Grad4SqM0=Grad4SqM;
Grad4SqM=highpass(Grad4SqM0,NSlice/TR,fs);
end
%---End of high pass---

Grad4=reshape(Grad4SqM,NRead*NCoil*NPhase*SlicePerExcite,GroupNum);
end %For SecondRemove

%---End of baseline removal---

%---Incorporate delays---
%Before=(0.04+0.207+0.138+0.025+1.4+0.12+0.142+0.12+0.148)/1000;
Before=6.44/1000;
%After=(0.199744+0.519233+0.138+0.199278)/1000;
BefPts=round(Before*fs);
dT=TR/NSlice-NRead*NPhase*NCoil/fs;
dNum=round(dT*fs);
NTran=round(TR/NSlice*fs); %The acquisition points for each EPI transient
AfterPts=NTran-NRead*NPhase*NCoil-BefPts;
fdemDelay=zeros(NTran,GroupNum); % The signal trace including delay
%fdemWhole2_LPM=reshape(fdemWhole2_LP,NRead*NPhase*NCoil,GroupNum);
for k=1:GroupNum
    fdemDelay(BefPts+1:BefPts+NRead*NPhase*NCoil,k)=Grad4(1:NRead*NPhase*NCoil,k);
    if k==1
        V3=Grad4(NRead*NPhase*NCoil,k);
        V4=Grad4(1,k+1);
        Inc2=(V4-V3)/(AfterPts+1);
        for LL=1:AfterPts
            fdemDelay(BefPts+NRead*NPhase*NCoil+LL,k)=V3+Inc2*LL;
        end       
    end
    if k>1 && k<GroupNum
%---Interpolate before acquisition time---        
        V1=Grad4(NRead*NPhase*NCoil,k-1);
        V2=Grad4(1,k);
        Inc=(V2-V1)/(BefPts+1);
        for LL=1:BefPts
            fdemDelay(LL,k)=V1+Inc*LL;
        end
%---Interpolate after acquisition time
        V3=Grad4(NRead*NPhase*NCoil,k);
        V4=Grad4(1,k+1);
        Inc2=(V4-V3)/(AfterPts+1);
        for LL=1:AfterPts
            fdemDelay(BefPts+NRead*NPhase*NCoil+LL,k)=V3+Inc2*LL;
        end
    end
    if k==GroupNum
        V1=Grad4(NRead*NPhase*NCoil,k-1);
        V2=Grad4(1,k);
        Inc=(V2-V1)/(BefPts+1);
        for LL=1:BefPts
            fdemDelay(LL,k)=V1+Inc*LL;
        end
    end
end

TimeInd3=[1:NTran*Grad2Size(2)]'/(NTran*Grad2Size(2))*TR*Rep;
if PlotPart==1
    figure;plot(TimeInd3(1:1000000),fdemDelay(1:1000000)/FreqVol);title('fdemDelay incorporate delays into Grad4');
else
    figure;plot(TimeInd3,squeeze(fdemDelay(:))/FreqVol);title('fdemDelay incorporate delays into Grad4');
end    

%---In the following, we add peaks together. First, define time points of stimulation---
PulseNum=20;
RepEachEpoch=15;
PInterval=200/1000;
RepWithEEG=ceil(PulseNum*PInterval/TR/NPhase);
RepNum=Rep;
EpochNum=RepNum/RepEachEpoch;

pks0=zeros(EpochNum,PulseNum); % 8 stands for epoch number, 24 stands for the number of pulses. 
pks0_Size=size(pks0);
Step=TR*RepEachEpoch; %The interval between adjacent epoch.
%Step=TR*NPhase*NCoil/9;
PInterval=1/5;
%Start=Step-1.12;
%Start=PInterval;
Start=0;
for k=1:pks0_Size(1)
    for m=1:pks0_Size(2)
        pks0(k,m)=Start+(k-1)*Step+(m-1)*PInterval;
    end
end
pksT=pks0';
pksT=squeeze(pksT(:));
pksT=pksT';

%Make a program to calculate pks
pks=round(pksT*fs);
pkSize=size(pks);
EleNum1=0;
EleNum2=round(fs*0.03);
Shape=zeros(EleNum1+EleNum2+1,1);
ShapeRaw=Shape;
fdemMatPreSq=squeeze(fdemDelay(:));

PickPul=[2 3 4 5 6 7 9 10 11 12 14 15];
%PickPul=[16 18 19 20 21 23 24 25 26 27 29 30 31 32 34 35]; % For 070823 E60
NumCurves = max(size(PickPul));
ShapeAll=zeros(EleNum1+EleNum2+1,NumCurves);

colormap('parula'); % You can choose any other colormap as well
colorMap = colormap;
colorIndices = round(linspace(1, size(colorMap, 1), NumCurves));
figure;
hold on
%for k=1:pkSize(2)
for L=1:NumCurves
%for L=1:1
    k=PickPul(L);
    ShapeAll(:,L)=fdemMatPreSq(pks(k)-EleNum1:pks(k)+EleNum2,1);
    Shape=Shape+fdemMatPreSq(pks(k)-EleNum1:pks(k)+EleNum2,1);
%    ShapeRaw=ShapeRaw+fdemMatSq(pks(k)-EleNum:pks(k)+EleNum,1);
    %fdemMatPreSq_Detrend=detrend(fdemMatPreSq(pks(k)-EleNum:pks(k)+EleNum,1));
    curveColor = colorMap(colorIndices(L), :);
    fdemMatPreSq_Detrend=(fdemMatPreSq(pks(k)-EleNum1:pks(k)+EleNum2,1));
    %plot([-EleNum:EleNum]/fs,fdemMatSq_Base2_Detrend);
    TimePart=[-EleNum1:EleNum2]/fs;
    plot(TimePart,fdemMatPreSq_Detrend,'Color',curveColor);
    %plot([-EleNum:EleNum]/fs,fdemSq_Peak(pks(k)-EleNum:pks(k)+EleNum,1));title('Sum EEG Shape');
end
Shape=Shape/NumCurves;

[ShapeMin,Xmin]=min(Shape);
[ShapeMax,Xmax]=max(Shape);
ShapeLen=max(size(Shape));
ShapeBase=mean(Shape(1:round(ShapeLen/10)));
ShapeAmp=ShapeMax-ShapeMin;      %Use this for amp
Tmin=(-EleNum1+Xmin-1)/fs*1000
Tmax=(-EleNum1+Xmax-1)/fs*1000

figure;plot(TimePart,Shape/FreqVol);title('Sum EEG Shape');

if SaveSeparateData==1
save(strcat(Path0,'\EEG'),'fdemDelay','TimeInd3');
end
