clear all;close all;
PlotImage=0;
PlotCorr=1;
SmoothImg=0;
Interlace=1;

ExpNum=[45];
Path0='.\';
data0S=0;
ExpNumSize=size(ExpNum);

for MM=1:ExpNumSize(2)
load(strcat(Path0,'EPIEEGS',mat2str(ExpNum(MM))));
data0P=SNR4D;
data0S=data0S+data0P;
end

SNR4DSize=size(SNR4D);
SliceNum=1;
PhaseNum=SNR4DSize(1);
ReadNum=SNR4DSize(2);
NSlice=SNR4DSize(3);
RepNum=SNR4DSize(4);
ReadNumPlot=ReadNum; %The actual number of read points to be plotted 

Size4D=size(SNR4D);

%---Save to NIFTI format---
NIFTI4D=zeros(Size4D(1),Size4D(2),Size4D(3),Size4D(4));
NSlice=Size4D(3);
for K=1:Size4D(4)
for M=1:Size4D(3)
    if Interlace==0
        M2=M;
    else
        if M/2~=round(M/2)
            M2=ceil(M/2);
        else
            M2=ceil(NSlice/2)+M/2;
        end
    end
    NIFTI4D(:,:,M,K)=(SNR4D(:,:,M2,K));
end
end
NIFTI4DP=permute(NIFTI4D,[2 1 3 4]);
NIFTI4DP2=NIFTI4DP;
niftiwrite(NIFTI4DP2,strcat(Path0,'EPI_Whole'));