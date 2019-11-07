

filename='test.lsm';

%% Read file into matrix using imreadBF (requires loci-tools.jar)

% get lsm/tif info and create 16bit (double) array
imageinfo=imreadBFmeta(filename);
x=imageinfo.width;
y=imageinfo.height;
z=imageinfo.zsize;
c=imageinfo.channels;
t=imageinfo.nframes;
phasechannel=2; 

imgarray=zeros(y,x,z); % annoyingly xy is inverted in matrix manipulation vs image axes

% read image one timeframe at a time to conserve memory
for current_t=1 % but only reading t=1 for this example

    %imreadBF can only handle single channel at a time
        imgarray(:,:,:)= squeeze(imreadBF(filename,1:z,current_t,phasechannel));


    
    
end

%% convert to rescaled grayscale (values between 0 and 1)
 gray = mat2gray(imgarray); % this fixes an issue with the lsm file format, where values are >1 

%% detect gradients of each pixel
%[x, y, z] = size(gray);
lineofbestfit= zeros(x*y,2);

gray = reshape(gray,x*y,z);
meanstacks= zeros(x*y,1);


%% using edited polyfit function with no warning loop - grayscale values won't trigger warnings, improves speed

for index=1:x*y

    lineofbestfit(index,:)=polyfit_nowarning(1:z,squeeze(gray(index,:)),1);
    meanstacks(index,1)=mean(gray(index,:));
end

meanstacks=reshape(meanstacks,x,y,1);

gray = reshape(gray,x,y,z);

%% first round of very crude segmentation to make subsections to improve focusing later

bwcrude=im2bw(meanstacks,graythresh(meanstacks));
bwcrude= imdilate(~bwcrude, strel('disk', 2));
bwcrude=imfill(bwcrude,'holes');
bwcrude= imerode(bwcrude, strel('disk',1));
smallobjectthreshold=20;
bwcrude = bwareaopen(bwcrude, smallobjectthreshold, 4); % removing small objects, threshold will depend on magnification
[bwcrude, ncrudeobjects]=bwlabel(bwcrude,4);

%% using A/y linear regression

%  for index=1:x*y
%      p = [ones(length(gray(index,:)),1) gray(index,:)']\[1:z]';
%      xy(index)=p(2);
%  end
% 

%  for i=1:x
%      for j=1:y
%         p=polyfit(1:k,squeeze(gray(i,j,:)),1);
%         xy(i,j)=p(1);
%      end
%      
%  end


%% convert to rescaled grayscale again
lineofbestfit = reshape(lineofbestfit,x,y,2);
xygray1=mat2gray(lineofbestfit(:,:,1));
xygray2=mat2gray(lineofbestfit(:,:,2));

%imshow(xygray1)
%imshow(xygray2)

%% find focused z slice - threshold high(cells) and low(halo around cells)
% detect focused stacks where these areas have similar values

% NB individual cell focusing could be improved by cutting each slice into
% subsections and evaluating independantly, either using crude mask for
% cell clusters or along y axis, as focus seems to vary along this axis.
% This would allow selection of a single stack for each subsection and
% likely improved segmentation.


level = graythresh(xygray1);
bwhigh = im2bw(xygray1,level); %cells
level = graythresh(xygray2);
bwlow = im2bw(xygray2,level); %halos

%% apply masks to original image z stack and calculate mean value of selected pixels at each slice through the stack
% now modified to use subsections

masked_cells=zeros(ncrudeobjects,z);
masked_halos=zeros(ncrudeobjects,z);
mean_cells=zeros(ncrudeobjects,z);
mean_halos=zeros(ncrudeobjects,z);
min_cell_slice= zeros(ncrudeobjects,1);
min_halo_slice= zeros(ncrudeobjects,1);

%stitched_cells = ones(x,y,1);
%stitched_halos = ones(x,y,1);

min_stitched_cells = ones(x,y,1);
min_stitched_halos = ones(x,y,1);
min_stitched_halos_stack = ones(x,y,1);

for i=1:ncrudeobjects
    [r,c]=find(bwlabel(bwcrude)==i);
    pixels_cells=zeros(length(r),1);
    pixels_halos=zeros(length(r),1);
    for j=1:z
        for k=1:length(r)
            
            %masked_cells(i,j) = sum(sum(gray(r(k),c(k),j).*(1-bwhigh)));
            %masked_halos(i,j) = sum(sum(gray(r(k),c(k),j).*(1-bwlow)));
    
            %mean_cells(i,j) =  mean(mean(gray(r(k),c(k),j).*(1-bwhigh)));
            %mean_halos(i,j) =  mean(mean(gray(r(k),c(k),j).*(1-bwlow)));

            %masked_gray1 = sum(sum(gray(:,:,i).*(1-BWhigh)));
            %masked_gray2 = sum(sum(gray(:,:,i).*(1-BWlow)));
            
            pixels_cells(k)= gray(r(k),c(k),j)*(1-bwhigh(r(k),c(k)));
            pixels_halos(k)= gray(r(k),c(k),j)*(1-bwlow(r(k),c(k)));
        end
        mean_cells(i,j)=sum(pixels_cells)/sum(pixels_cells~=0);
        mean_halos(i,j)=sum(pixels_halos)/sum(pixels_halos~=0);    
    end
    [minvalue_cell min_cell_slice(i)] =min(mean_cells(i,:));
    [minvalue_halo min_halo_slice(i)] =min(mean_halos(i,:));
    for l=1:length(r)
        %stitched_cells(r(l),c(l))=gray(r(l),c(l),min_cell_slice(i));
        %stitched_halos(r(l),c(l))=gray(r(l),c(l),min_halo_slice(i));
        three_cells= min_cell_slice(i)-1:min_cell_slice(i)+1;
        three_halos= min_halo_slice(i)-1:min_halo_slice(i)+1;
        min_stitched_cells(r(l),c(l))=min(gray(r(l),c(l),nonzeros(three_cells(three_cells~=z+1))'));
        min_stitched_halos(r(l),c(l))=min(gray(r(l),c(l),nonzeros(three_halos(three_halos~=z+1))'));
        min_stitched_halos_stack(r(l),c(l))=min(gray(r(l),c(l),1:min_halo_slice(i)));

    end
        
end

%% plotting over time lets us visualise the area of focus - focused stacks are at intersection of lines

%plot(1:z,masked_gray1,1:z,masked_gray2)
%plot(1:z,mean_cells(i,:),1:z,mean_halos(i,:))

%% select z stacks of interest based on focusing and combine
%%i will select 3 stacks around the minimum - depends on how deep through
%%the z stack the cells extend. can be improved as detailed in earlier
%%comment


bwhalo = im2bw(min_stitched_halos,graythresh(min_stitched_halos));

bwcell = im2bw(min_stitched_cells,graythresh(min_stitched_cells));

%% dilate cells and combine with halos to close gaps in outlines

se = strel('disk',1);
dilatedcellmask = imdilate(~bwcell,se);
closedcells=bwhalo&dilatedcellmask;
imshow(closedcells)

%% output to files

Numberframes=1;

for p=1:Numberframes
    outname= ['p_' FileTif(1:end-4) '-p-' num2str(p,'%03d') '.tif'];
    %imwrite(closedcells,outname);
    imwrite(imresize(uint16((2-closedcells)*30000),2),'closed-p-001.tif','Compression','none');
    %outTiff = Tiff(outname,'w');
    %
    %outTiff.write(xygray1);
    %outTiff.close();
    %xygray2=mat2gray(lineofbestfit(:,:,2));
end
rgb=zeros(x,y,3,Numberframes);
rgb(:,:,1,:)=xygray1(:,:,:);
rgb(:,:,2,:)=xygray1(:,:,:);
rgb(:,:,3,:)=xygray1(:,:,:);


% vidObj = VideoWriter('Phasemovie.avi');
% vidObj.FrameRate=10;
% open(vidObj);
% 
% for p=1:Numberframes
%      frames(p) = im2frame(rgb(:,:,:,p));
% end
% vidObj.writeVideo(frames);
% close(vidObj);



% 
% D = bwdist(~bw);
% D = -D;
% D(~bw) = -Inf;
% L = watershed(D);

%% detect gradient again without focused slices
% gray = reshape(gray,x*y,z);
% 
% for index=1:x*y
% 
%     lineofbestfit(index,:)=polyfit_nowarning(1:z,squeeze(gray(index,:)),1);
%          
% end
% 
% gray = reshape(gray,x,y,z);
% 
% 



%% code graveyard

%% libtiff reading tif function - imreadBF is better
% 
% TifLink = Tiff(FileTif, 'r');
% for i=1:NumberImages
%    TifLink.setDirectory(i);
%    imgarray(:,:,i)=TifLink.read();
% end
% TifLink.close();

%%

%masked_xygray2 = gray.*(1-BWlow);

%%remove irrelevant pixels
% masked_xygray1=masked_xygray1(masked_xygray1~=0);
% masked_xygray2=masked_xygray2(masked_xygray2~=0);
%highvalue= zeros(x,y,1);
%zslice(:,:,1)= -lineofbestfit(:,:,1)/lineofbestfit(:,:,2);

% 
% %%invert and saturate
% comp_xy= xygray2;
% %comp_xy= imcomplement(xygray);
% %comp_xy=imadjust(comp_xy);
% 
% %%threshold
% level = graythresh(comp_xy);
% BW = im2bw(comp_xy,level);
% imshow(BW)


% zxy= permute(gray,[3 1 2]);
% for i = 1%size(gray,3)
%     [Gz,Gx] = imgradientxy(zxy(:,:,i));
%     xy(i,:)= Gz(i,:)
% end
