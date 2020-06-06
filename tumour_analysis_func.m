function [met_cell_np_int, met_cell_dist, met_dilate_nuclei_label] = tumour_analysis_func(pre_nuclei,post_nuclei,post_vessels,nanoparticle_ch,save_dir,sample_name)

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

shortfile = sample_name;
display(['Processing ' shortfile])
tic

%Trim the tissue boundary using the nuclei channel
pre_nanoparticle_file = uigetfile('','Select pre-processed nanoparticle image file');
pre_nanoparticle_info = imfinfo(pre_nanoparticle_file);
pre_nanoparticle = imread(pre_nanoparticle_file,1);
for ii = 2 : size(pre_nanoparticle_info, 1)
    temp_pre_nanoparticle_tiff = imread(pre_nanoparticle_file, ii);
    pre_nanoparticle = cat(3 , pre_nanoparticle, temp_pre_nanoparticle_tiff);
end
dapi = pre_nanoparticle(:,:,1);

    fill_tissue = fillholes(dapi>40);
    erode_tissue = erosion(fill_tissue,14,'elliptic');
    open_tissue = opening(erode_tissue,10,'elliptic');
    gauss_tissue = gaussf(open_tissue,3);
    block_tissue = label(gauss_tissue>0.5,1,50000000,0); 
    tissue_area = uint8(block_tissue>0);

tissue_area_bin = tissue_area>0;

%nuclei_post_processed = post_nuclei;


vessel_file = uigetfile('','Select post-processed vessel image file');
vessel_info = imfinfo(vessel_file);
post_vessels = imread(vessel_file,1);
for ii = 2 : size(vessel_info, 1)
    temp_vessels_tiff = imread(vessel_file, ii);
    post_vessels = cat(3 , post_vessels, temp_vessels_tiff);
end
vessels_post_processed = post_vessels;

%%% Load particle channel and trim to tissue boundary

nanoparticle_file = uigetfile('','Select original nanoparticle image file');
nanoparticle_info = imfinfo(nanoparticle_file);
nanoparticle_ch = imread(nanoparticle_file,1);
for ii = 2 : size(nanoparticle_info, 1)
    temp_nanoparticle_tiff = imread(nanoparticle_file, ii);
    nanoparticle_ch = cat(3 , nanoparticle_ch, temp_nanoparticle_tiff);
end
particle = uint16(nanoparticle_ch);




particle_crop = particle;%.*uint16(tissue_area_bin);

%%% Analysis of tumour tissue
%nuclabel = label(nuclei_post_processed>0,1,25,0);
%nuc_dilated_label = dip_growregions(nuclabel,[],[],1,2,'low_first');
%nuc_dilated_label = uint32(nuc_dilated_label);

%Ind_tumor_stats = regionprops3(tissue_area_bin,particle_crop,'Centroid','MeanIntensity','Volume','SurfaceArea');

vess_thresh = vessels_post_processed>0;
dt_vessels = round(bwdist(vess_thresh));
vess_lumen_thresh = vessels_post_processed==0;
dt_lumen_vessels = round(bwdist(vess_lumen_thresh));

all_NP_Vessel_int_stats = regionprops3(dt_vessels,particle,'MeanIntensity','VoxelIdxList');
all_cells_int_stats_struct = table2struct(all_NP_Vessel_int_stats);

all_NP_Lumen_int_stats = regionprops3(dt_lumen_vessels,particle,'MeanIntensity','VoxelIdxList');
all_cells_int_stats_struct = table2struct(all_NP_Lumen_int_stats);

all_cells_dist_stats = regionprops3(particle,dt_vessels,'MeanIntensity','VoxelIdxList');
all_cells_dist_stats_struct = table2struct(all_cells_dist_stats);



%max_cell_in_tissue = max(nuc_dilated_label(:));
%min_cell_in_tissue = min(nuc_dilated_label(:))+1;

%new_all_nuc_int =  uint32(nuc_dilated_label);     
    
%    for b = min_cell_in_tissue:max_cell_in_tissue
%            new_all_nuc_int(all_cells_int_stats_struct(b).VoxelIdxList) = all_cells_int_stats_struct(b).MeanIntensity;
%    end

%new_all_nuc_dist =  uint32(nuc_dilated_label);     

%    for b = min_cell_in_tissue:max_cell_in_tissue
%            new_all_nuc_dist(all_cells_dist_stats_struct(b).VoxelIdxList) = all_cells_dist_stats_struct(b).MeanIntensity;
%    end
    
%new_all_nuc_int_crop = new_all_nuc_int;    
%new_all_nuc_dist_crop = new_all_nuc_dist;      
%nuclei_dilated_post_processed_crop = nuc_dilated_label;      

%met_cell_np_int = new_all_nuc_int_crop;
%met_cell_dist = new_all_nuc_dist_crop;
%met_dilate_nuclei_label = nuclei_dilated_post_processed_crop;


%%% Writing tables and images
cd(save_dir)
all_NP_Vessel_int_stats.Properties.VariableNames = {'VoxelIdxList_NPint' 'MeanIntensity_NPint'};
all_NP_Lumen_int_stats.Properties.VariableNames = {'VoxelIdxList_NPint' 'MeanIntensity_NPint'};
all_cells_dist_stats.Properties.VariableNames = {'VoxelIdxList_dist' 'MeanIntensity_dist'};

Tumor_vessel_npint_Int = [all_NP_Vessel_int_stats.MeanIntensity_NPint];% all_cells_dist_stats.MeanIntensity_dist];
Tumor_vessel_table = array2table(Tumor_vessel_npint_Int);
Tumor_vessel_name_nuc = strcat(shortfile,'-Tumor-vessels-all','.csv');
writetable(Tumor_vessel_table,Tumor_vessel_name_nuc);

Tumor_lumen_npint_Int = [all_NP_Lumen_int_stats.MeanIntensity_NPint];% all_cells_dist_stats.MeanIntensity_dist];
Tumor_lumen_table = array2table(Tumor_lumen_npint_Int);
Tumor_lumen_name_nuc = strcat(shortfile,'-Tumor-lumen-all','.csv');
writetable(Tumor_lumen_table,Tumor_lumen_name_nuc);

vessel_regions_name = strcat(shortfile,'_vessel_regions.tif');

num_slices = size(dt_vessels,3);

 imwrite(uint8(dt_vessels(:,:,1)),vessel_regions_name);
         
    for p = 2:num_slices
            imwrite(uint8(dt_vessels(:,:,p)),vessel_regions_name, 'WriteMode','append');
    end  

%all_met_table_name = strcat(shortfile,'-Tumor-area-','.csv');    
%writetable(Ind_tumor_stats, all_met_table_name);

%tumor_mean_np_int_name = strcat(shortfile,'Tumor mean NP Intensity','.tif');
%tumor_mean_dist_name = strcat(shortfile,'Tumor mean distance to vessel','.tif');
%tumor_nuclei_label_name = strcat(shortfile,'Tumor labeled cells','.tif');

%num_slices = size(tissue_area,3);

%imwrite(uint16(new_all_nuc_int_crop(:,:,1)),tumor_mean_np_int_name);
         
%  for p = 2:num_slices
%          imwrite(uint16(new_all_nuc_int_crop(:,:,p)),tumor_mean_np_int_name, 'WriteMode','append');
%  end

%imwrite(uint16(new_all_nuc_dist_crop(:,:,1)),tumor_mean_dist_name);
         
%  for p = 2:num_slices
%          imwrite(uint16(new_all_nuc_dist_crop(:,:,p)),tumor_mean_dist_name, 'WriteMode','append');
%  end

%imwrite(uint16(nuclei_dilated_post_processed_crop(:,:,1)),tumor_nuclei_label_name);
         
%  for p = 2:num_slices
%          imwrite(uint16(nuclei_dilated_post_processed_crop(:,:,p)),tumor_nuclei_label_name, 'WriteMode','append');
%  end
  toc
  
   
cd(matlab_folder)


end