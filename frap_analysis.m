function [D] = frap_analysis(image_path)
%frap_analysis Summary of this function goes here
%   Detailed explanation goes here
%%
for k=1:180
    image(:,:,k)=imread(image_path,k);
end
gaussfilt=imgaussfilt(image(:,:,1),5);
T=graythresh(gaussfilt);
BW=~imbinarize(gaussfilt,T);
BW= imclearborder(BW,4);
cc = bwconncomp(BW);
props = regionprops(cc,'Area','pixelidxlist','Centroid','MajorAxisLength','MinorAxisLength');
roi = [props.Area]<2000;
idx = vertcat(props(roi).PixelIdxList);
region_image=BW;
region_image(idx)=false;
roi = [props.Area]>2000;
centers=props(roi).Centroid;
diameters = mean([props(roi).MajorAxisLength props(roi).MinorAxisLength],2);
radii=diameters/2;
imshow(image(:,:,1));
hold on
viscircles(centers,radii);
%How to check background circle is in image
viscircles(centers+4*radii,radii);

A=mean_intensity(centers, radii, size(image(:,:,1)),image);
C=mean_intensity(centers+4*radii, radii, size(image(:,:,1)),image);
%Need A(1) to be larger than A so put minus 1 for now
corrected_bleach=((A.*(C(1)./C))./A(1))-1;
t=[1:180];
customFitType=fittype('a*(1-exp(-b*x))+c','independent','x','dependent','y');
options=fitoptions('Method','NonlinearLeastSquares');
f=fit(t',corrected_bleach',customFitType,options);
disp(f);
figure;plot(f,t,corrected_bleach);

D=(-0.88*((radii/4.4563)^2)*f.b)/(4*log(0.5));
end
%%
function [A] = mean_intensity(centers, radii, dim,image)
mask = circles2mask(centers,radii,dim);
A=zeros(1,180);
for k=1:180
    timepoint=image(:,:,k);
    A(k)=mean(timepoint(mask));
end
end