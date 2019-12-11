%David Meyer Final Project:Skin Cancer
%CS390S
%11/26/2019 to 12/0#/2019

%Can't work on all skin cancer, so work on just one small branch to start
%and then work on the next one...



%Decided to startfigu working on just Melnoma skin cancer because it has a lot
%of different visual ques to show that it needs to be tested by a doctor.
%Marking this, even with a 100%, with this program only shows that it is
%"most likely" not a gaurantee to be malignant(Cancerous).



%% RESEARCH AREA
   %Picture used form: https://upload.wikimedia.org/wikipedia/commons/6/6c/Melanoma.jpg
   
   %According to my incomplete nursing degree, Melanoma is a skin cancer
   %that you can detect using ASBCDE acronym.  
        %A(Asymmetric-one side lop-sided compared to the other.)
        %B(Border irregularity-Growths with irregular, notched or scalloped
            %borders)
        %C (Color Changes-multiple colors or uneven distribution)
        %D(Diameter-Bigger than a pencil eraser(about 1/4 inches or 6 
            %millemeters))
        %E (Evolving-Anything else, such as a mole that changes size, 
            %color, or shape).


%% Programing Talk For Self
    %So, I can look for A,B,C, and D from above list.  Run either  a
    %running total, or booleans.  Probably best to do a sum to determine
    %likelyhood to go to the doctor by saying it's cancerous.  Else, if
    %boolean, then I need something that does all the booleans. 


%% ACTUAL PROGRAMMING
clc
clear all
close all

totalSum = 0;
picOrig = imread('Pictures/Melanoma.jpg'); %Might need original later
bwPic = rgb2gray(picOrig);
skinPic = skinDetect(picOrig);




%% A (Asymmetric)
  %Approach One: Divied in half and then determine if one side is bigger
        %Then rotate to see if it's bigger top and bottom...
  %First, HAVE/NEED to get rid of anything that isn't "Black".  THenhreshold,
  %skin detection... Except change the skin to white.

avgPic = avgFilter(skinPic);
figure,montage({picOrig,skinPic,avgPic})

        %Do the math to determine A
        %image(x,y,z(color));
left = avgPic(:,1:end/2,:);
right = avgPic(:,end/2+1:end,:);
%figure,imshow(left)
%figure,imshow(right)
symmQuestion = left - right;
sumOf = sum(symmQuestion,'all');
%NEED to create a better IF-Statement? Good for now....
if((sumOf >= 10) || (sumOf <= -10)) 
    totalSum = totalSum + 1;
end



%% B (Border Irregularity) Might want to create a rgb2hsv to get better
%borders
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





%% C (Color change) 
%HSV might be best for this? or YCrCb....








%% D (Diameter)Find the center of the spot, then calculate the distance to the
%border?


%% ALL FUNCTIONS
function avgPicOut = avgFilter(usePic)
    avgFilt=[1 1 1;
             1 1 1;
             1 1 1];%Divide matrix by 9 for the average of them all
    %xDone=3;
    %for x=0:xDone%Should I be looping this to get a better smoothing?
        avgFilt = avgFilt * (1/9);
    
        %avgFilt=[2 2 2;2 2 2;2 2 2]/9;
        %avgPicOut = filter2(usePic,avgFilt,'full');
        %mesh(avgPicOut)
        avgPicOut = conv2(usePic,avgFilt);
    %end
end

function med = medFilt(usePic)%There has to be a better way to implement this...
    [rows, col] = size(usePic);%size of grayscale image
    picSalt = usePic;%Unsure why I did this, to lazy to change now. 
    pad=zeros(rows,col);
    for i=2:rows-1
        for j=2:col-1
            %Make 3x3 mask
            filter = [picSalt(i-1,j-1),picSalt(i-1,j),picSalt(i-1,j+1),picSalt(i,j-1),picSalt(i,j),picSalt(i,j+1),picSalt(i+1,j-1),picSalt(i+1,j),picSalt(i+1,j+1)];
            pad(i,j)= median(filter);%function that just return median value
        end
    end
    med = pad;
end

function sobelKernel = sobFilter(usePic)
    sobFiltOne=[-1 0 +1;
                -2 0 +2;
                -1 0 +1];%xDirection is done
    sobFiltTwo=[+1 +2 +1;
                 0 0 0;
                -1 -2 -1];%yDirection is done
    %sobelFiltOneImage = imfilter(usePic,sobFiltOne,'same');
    %sobelFiltTwoImage = imfilter(usePic,sobFiltTwo,'same');
        %Working is above.
    sobelFiltOneImage = conv2(usePic,sobFiltOne);%Another function to loop through picture to apply Kernel
    sobelFiltTwoImage = conv2(usePic,sobFiltTwo);
    %sobelFilt#Image = conv2(usePic,sobfilt#,'full');
        %Why does filter2 flip the image upside down?
            %conv2(pic, mask) = filter2(rot90(mask,2), pic)
            %conv2 is a bit faster than fitler, no reason for me to use
            %filter2 here.  
    
    sobelKernel = sqrt(sobelFiltOneImage.^2 + sobelFiltTwoImage.^2);
      %Edges detected better due to reading the image as a double from the
      %begging. 
end

%this program is to get the normalized histogram of the image

function [his]=myhist(x)
    x=floor(x);
    r=size(x,1);
    c=size(x,2);
    his=zeros(1,256);
    %calculate the histogram
    for i=1:r
        for j=1:c
            his(x(i,j)+1)=his(x(i,j)+1)+1;
        end
    end
    his=his./c./r;
end





%After working
    %Make an all encompassing window to show everything and the definitive
    %thoughts for users to have an easier experience.
    
    %Maybe use PCA to find the characteristics of each type of cancer
    %instead of doing all this?