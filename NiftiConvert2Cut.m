clear all;close all;
Factor=1;
FOVx=26;FOVy=26;SliceThick=0.6;

Path0='.\';
FileName=strcat(Path0,'EPI_Whole','.nii');
info1=niftiinfo(FileName);
Img1=niftiread(FileName);

info1.PixelDimensions(1)=info1.PixelDimensions(1);
info1.PixelDimensions(2)=info1.PixelDimensions(2);
info1.PixelDimensions(3)=info1.PixelDimensions(3);

%Cut the center 1/3 portion, fit onto standard template
Img3=Img1;
info4=info1;
info4.ImageSize(1)=round(info1.ImageSize(1)/3);
info4.ImageSize(2)=info1.ImageSize(2);
info4.ImageSize(3)=info1.ImageSize(3);
Img4=Img3(info4.ImageSize(1)+1:info4.ImageSize(1)*2,info4.ImageSize(2):-1:1,:,:);
%Img4=Img3(info4.ImageSize(1)+1:info4.ImageSize(1)*2,:,:,:);
Img5=permute(Img4,[1 3 2 4]);
info5=info4;
info5.ImageSize(1)=info4.ImageSize(1);
info5.ImageSize(2)=info4.ImageSize(3);
info5.ImageSize(3)=info4.ImageSize(2);
info5.PixelDimensions(1)=FOVx/info4.ImageSize(1);
info5.PixelDimensions(2)=SliceThick;
info5.PixelDimensions(3)=FOVy/info4.ImageSize(2);
info5.TimeUnits='Second';
FileName5=strcat(Path0,'EPI_OneThird','.nii');
niftiwrite(Img5,FileName5,info5);


