


function [post_vessels, post_micromet] = post_process_func(seg_vessels,seg_nanoparticle,save_dir,sample_name)

% set matlab folder to return at the end
matlab_folder = pwd;
% select folder with the data to read files
data_folder = uigetdir('','Select folder with image files'); 
cd(data_folder);
% input sample name for rename in the terminal
% sample_name = input('Enter the name of the smple: ','s');
sample_name = strcat(cell2mat(inputdlg('Enter the name of the smple:')));
%select save directory
save_dir = uigetdir('','Select save folder');

% select and read vessel file
vessel_file = uigetfile('','Select trained vessel image file');
vessel_info = imfinfo(vessel_file);
seg_vessels = imread(vessel_file,1);
for ii = 2 : size(vessel_info, 1)
    temp_vessels_tiff = imread(vessel_file, ii);
    seg_vessels = cat(3 , seg_vessels, temp_vessels_tiff);
end

% select and read nanoparticle file
%nanoparticle_file = uigetfile('','Select nanoparticle image file');
%nanoparticle_info = imfinfo(nanoparticle_file);
%seg_nanoparticle = imread(nanoparticle_file);
%for ii = 2 : size(vessel_info, 1)
%    temp_nanoparticle_tiff = imread(nanoparticle_file, ii);
%    seg_nanoparticle = cat(3 , seg_nanoparticle, temp_nanoparticle_tiff);
%end


shortfile = sample_name;
display(['Post-processing ' shortfile])


% select and read pre-processed nanoparticle file
pre_nanoparticle_file = uigetfile('','Select the pre-processed nanoparticle image file');
pre_nanoparticle_info = imfinfo(pre_nanoparticle_file);
pre_nanoparticle = imread(pre_nanoparticle_file);
for ii = 2 : size(pre_nanoparticle_info, 1)
    temp_pre_nanoparticle_tiff = imread(pre_nanoparticle_file, ii);
    pre_nanoparticle = cat(3 , pre_nanoparticle, temp_pre_nanoparticle_tiff);
end

    tissues = pre_nanoparticle;

    fill_tissue = fillholes(tissues>40);
    erode_tissue = erosion(fill_tissue,14,'elliptic');
    open_tissue = opening(erode_tissue,10,'elliptic');
    gauss_tissue = gaussf(open_tissue,3);
    block_tissue = label(gauss_tissue>0.5,1,50000000,0); 
    tissue_area = uint8(block_tissue>0);
    tissue_area_bin = tissue_area>0;    

%    num_slices = size(tissue_area,3);

nanoparticle_file = uigetfile('','Select trained nanoparticle image file');
nanoparticle_info = imfinfo(nanoparticle_file);
nanoparticle_ch = imread(nanoparticle_file,1);
for ii = 2 : size(nanoparticle_info, 1)
    temp_nanoparticle_tiff = imread(nanoparticle_file, ii);
    nanoparticle_ch = cat(3 , nanoparticle_ch, temp_nanoparticle_tiff);
end
particle = uint16(nanoparticle_ch);




particle_crop = uint8(particle>=0.2);

Particle_true = pre_nanoparticle.*particle_crop;


%%%%Load in Ilastik segmented image files
%nuclei_seg = seg_nuclei;
vessels_seg = seg_vessels;
%ki67_seg = seg_micromets;

%%%%%Post-processing of Ilastik segmented nuclei channel
%nuclei_seg_bin = nuclei_seg==1;

%nuclei_seg_crop = (nuclei_seg_bin.*tissue_area_bin)>0;


%threshnuc = nuclei_seg_crop>0;
%    threshnuc3 = opening(threshnuc,4);
%    dt_threshnuc = dt(threshnuc3);
    
%    seeds = maxima(dt_threshnuc,2,0);
%    seeds2 = dilation(seeds,4.5)>0;
%    image_out = waterseed(seeds2,max(dt_threshnuc)-dt_threshnuc,1,0,0);
    
%    threshnuc4 = (uint8(threshnuc3))>0;
%    threshnuc4_dil_se = strel('sphere',1);
%    threshnuc4_dil = imdilate(threshnuc4,threshnuc4_dil_se);
%    threshnuc4_dil(image_out) = false; 
    
%nuclabel = label(threshnuc4_dil>0,1,25,0);
%nuc_dilated_label = dip_growregions(nuclabel,[],[],1,2,'low_first');
%final_nuc = nuclabel;
%final_dilated_nuc = nuc_dilated_label;

%nuclei_processed = uint32(final_nuc);
%nuclei_dilated_processed = uint32(final_dilated_nuc);



%%%%%Post-processing of Ilastik segmented blood vessel channel
vessel_seg_bin = vessels_seg >= 0.2;
vessel_seg_crop = vessel_seg_bin.*tissue_area_bin;
se_erode = strel('sphere',2);
vessel_seg_crop_er = imerode(vessel_seg_bin,se_erode);

vessels_processed = vessel_seg_crop_er>0;



%%%%%Post-processing of Ilastik segmented ki67/metastases channel
%ki67_bin = ki67_seg==1;
%ki67_seg_crop = ki67_bin.*tissue_area_bin;
%se_open_ki67 = strel('sphere',3);
%ki67_open = imopen(ki67_seg_crop,se_open_ki67);

%ki67_label = label(ki67_open>0,1,30000,0);
%ki67_processed = uint16(ki67_label);

%%%%%Final segmentation images
%post_nuclei = nuclei_processed;
post_vessels = vessels_processed;
%post_micromet = ki67_processed;
%post_nuclei_dilate = nuclei_dilated_processed;
%%%%%Write post-processed files

%nuclei_processed_name = strcat(shortfile,'_post_processed_nuclei.tif');
%nuclei_dilated_processed_name = strcat(shortfile,'_post_processed_dialted_nuclei.tif');
vessels_processed_name = strcat(shortfile,'_post_processed_vessels.tif');
nanoparticle_processed_name = strcat(shortfile,'_post_processed_nanoparticle.tif');
%ki67_processed_name = strcat(shortfile,'_post_processed_ki67.tif');

num_slices = size(vessels_seg,3);

cd(save_dir)

%imwrite(uint16(nuclei_processed(:,:,1)),nuclei_processed_name);
         
%  for p = 2:num_slices
%          imwrite(uint16(nuclei_processed(:,:,p)),nuclei_processed_name, 'WriteMode','append');
%  end

%imwrite(uint16(nuclei_dilated_processed(:,:,1)),nuclei_dilated_processed_name);
         
%  for p = 2:num_slices
%          imwrite(uint16(nuclei_dilated_processed(:,:,p)),nuclei_dilated_processed_name, 'WriteMode','append');
%  end  

 imwrite(uint8(vessels_processed(:,:,1)),vessels_processed_name);
         
    for p = 2:num_slices
            imwrite(uint8(vessels_processed(:,:,p)),vessels_processed_name, 'WriteMode','append');
    end   


 imwrite(uint8(Particle_true(:,:,1)),nanoparticle_processed_name);
         
    for p = 2:num_slices
            imwrite(uint8(Particle_true(:,:,p)),nanoparticle_processed_name, 'WriteMode','append');
    end    
    
% imwrite(uint8(ki67_processed(:,:,1)),ki67_processed_name);
         
 %   for p = 2:num_slices
 %            imwrite(uint8(ki67_processed(:,:,p)),ki67_processed_name, 'WriteMode','append');
 %   end   

 
cd(matlab_folder)

clc

display(['Pre-processing completed: ' shortfile])
 
end


