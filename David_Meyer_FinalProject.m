%David Meyer Final Project:Skin Cancer
%CS390S
%11/26/2019 to 12/0#/2019

%Can't work on all skin cancer, so work on just one small branch to start
%and then work on the next one...



%Decided to start working on just Melnoma skin cancer because it has a lot
%of different visual ques to show that it needs to be tested by a doctor.
%Marking this, even with a 100%, with this program only shows that it is
%"most likely" not a gaurantee to be malignant(Cancerous).



%% RESEARCH AREA
   %Picture used form: https://upload.wikimedia.org/wikipedia/commons/6/6c/Melanoma.jpg
   
   %According to my incomplete nursing degree, Melanoma is a skin cancer
   %that you can detect using ASBCDE acronym.  
        %A(Asymmetric-one side darker
            %than the other is an example)
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
%A (Asymmetric)
  %Approach One: Divied in half and then determine if one side is bigger
        %Then rotate to see if it's bigger top and bottom...
  %First, HAVE/NEED to get rid of anything that isn't "Black".  Threshold,
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



%B (Border Irregularity) Might want to create a rgb2hsv to get better
%borders






%C (Color change) 
%HSV might be best for this? or YCrCb....








%D (Diameter)Find the center of the spot, then calculate the distance to the
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


