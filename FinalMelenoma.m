

clear all
close all

ei = 25;
st = 35;

k=ei*st;
img = imread('Pictures/Melanoma.jpg');

h = ones(ei,st)/k;
I = imfilter(img,h,'symmetric');

figure
subplot(2,2,1), imshow(img), title('Orginal');
subplot(2,2,2), imshow(I), title('Filtered Image');

%Converting to BW
Igray = rgb2gray(I);
I1 = imadjust(Igray,stretchlim(Igray),[]);
level = graythresh(I1);
BWj = im2bw(I1,level);
dim = size(BWj);
IN = ones(dim(1),dim(2));
BW = xor(BWj,IN);  %Inverting
subplot(2,2,3), imshow(BW), title('Black and White');

%Finding of an initial point
row = round(dim(1)/2);
col = min(find(BW(row,:)));

%Tracing
boundary = bwtraceboundary(BW, [row, col], 'W');
subplot(2,2,4),imshow(img), title('Traced');
hold on;

%Display traced boundary
plot(boundary(:,2),boundary(:,1),'g','LineWidth',2);
hold off
%
%

nn = size(boundary);
KM=zeros(dim(1),dim(2));
    ii=0;
    %Create new matrix with boundary points.
    %Other distortions outside boundaries
    while ii<nn(1)
        ii=ii+1;
        KM(boundary(ii,1),boundary(ii,2))=1;
    end
    figure
    subplot(2,2,1),plot(boundary(:,2),boundary(:,1),'black','LineWidth',2);
    subplot(2,2,2),imshow(KM)
    
%Fill inner boundaries where lesion is located
KM2 = imfill(KM,'holes');
subplot(2,2,3),imshow(KM2)
KM1 = xor(KM2,IN);
%subplot


%Geometricxal center
IVx = [1:dim(2)];
IVy = [1:dim(1)];
IMx = ones(dim(1),1)*IVx;
IMy = ones(dim(2),1)*IVy;
IMy = imrotate(IMy,-90);
Koordx = IMx.*KM2;
Koordy = IMy.*KM2;
xmean = mean(Koordx,2);
yc = round(sum(xmean.*IMy(:,1))/sum(xmean));
ymean = mean(Koordy);
xc = round(sum(ymean.*IVx)/sum(ymean));
figure, imshow(img), hold on
plot(boundary(:,2),boundary(:,1),'green','LineWidth',2);
hold on
plot(xc,1:dim(1),'red','LineWidth',2);
plot(1:dim(2),yc,'red','LineWidth',2);
hold off

% id = im2double(img);
ID1(:,:,1) = im2double(img(:,:,1));
ID1(:,:,2) = im2double(img(:,:,2));
ID1(:,:,3) = im2double(img(:,:,3));
figure, subplot(2,2,1), imshow(ID1);
subplot(2,2,2), imshow(ID1(:,:,1));
hold on
plot(xc,1:dim(1),'red','LineWidth',2);
plot(1:dim(2),yc,'red','LineWidth',2);
hold off
subplot(2,2,3), imshow(ID1(:,:,2));
subplot(2,2,4), imshow(ID1(:,:,3));





