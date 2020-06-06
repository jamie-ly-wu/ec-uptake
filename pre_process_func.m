

function [pre_vessels, pre_nanoparticle] = pre_process_func(vessels,nanoparticle,save_dir,sample_name)

% set matlab folder to return at the end
matlab_folder = pwd;
% select folder with the data to read files
data_folder = uigetdir('','Select folder with image files'); 
cd(data_folder);

% select and read vessel file
vessel_file = uigetfile('','Select vessel image file');
vessel_info = imfinfo(vessel_file);
vessels = imread(vessel_file,1);
for ii = 2 : size(vessel_info, 1)
    temp_vessels_tiff = imread(vessel_file, ii);
    vessels = cat(3 , vessels, temp_vessels_tiff);
end

% select and read nanoparticle file
nanoparticle_file = uigetfile('','Select nanoparticle image file');
nanoparticle_info = imfinfo(nanoparticle_file);
nanoparticle = imread(nanoparticle_file);
for ii = 2 : size(nanoparticle_info, 1)
    temp_nanoparticle_tiff = imread(nanoparticle_file, ii);
    nanoparticle = cat(3 , nanoparticle, temp_nanoparticle_tiff);
end

% input sample name for rename in the terminal
% sample_name = input('Enter the name of the smple: ','s');
sample_name = strcat(cell2mat(inputdlg('Enter the name of the smple:')));
%select save directory
save_dir = uigetdir('','Select save folder');
shortfile = sample_name;
display(['Pre-processing ' shortfile])


%Pre-processing of blood vessels channel
vessels_single = single(vessels);
vessels_local = vessels_single./(0.2*max(vessels_single(:))+(im2mat(gaussf(vessels_single,size(vessels_single,1)/10))));
vessels_loglocal = mat2im(log(vessels_local+0.1));

vessels_processed = uint8(stretch(vessels_loglocal));
pre_vessels = vessels_processed;

%Pre-processing of nanoparticle channel
nanoparticle_single = single(nanoparticle);
nanoparticle_local = nanoparticle_single./(0.2*max(nanoparticle_single(:))+(im2mat(gaussf(nanoparticle_single,size(nanoparticle_single,1)/10))));
nanoparticle_loglocal = mat2im(log(nanoparticle_local+0.1));

nanoparticle_processed = uint8(stretch(nanoparticle_loglocal));
pre_nanoparticle = nanoparticle_processed;

%Pre-processing of ki67 channel
%ki67_single = single(micromets);
%ki67_local = ki67_single./(0.2*max(ki67_single(:))+(im2mat(gaussf(ki67_single,size(ki67_single,1)/10))));
%ki67_loglocal = mat2im(log(ki67_local+0.1));

%ki67_processed = uint8(stretch(ki67_loglocal));
%pre_micromet = ki67_processed;

%Write pre-processed files
cd(save_dir)

nanoparticle_processed_name = strcat(shortfile,'_pre_processed_nanoparticle.tif');
vessels_processed_name = strcat(shortfile,'_pre_processed_vessels.tif');
%ki67_processed_name = strcat(shortfile,'_pre_processed_micromets.tif');

num_slices = size(vessels,3);



imwrite(uint8(vessels_processed(:,:,1)),vessels_processed_name);
        
   for p = 2:num_slices
            imwrite(uint8(vessels_processed(:,:,p)),vessels_processed_name, 'WriteMode','append');
   end
   
   
imwrite(uint8(nanoparticle_processed(:,:,1)),nanoparticle_processed_name);
        
   for p = 2:num_slices
            imwrite(uint8(nanoparticle_processed(:,:,p)),nanoparticle_processed_name, 'WriteMode','append');
   end
   
   
%imwrite(uint8(ki67_processed(:,:,1)),ki67_processed_name);
        
   %for p = 2:num_slices
            %imwrite(uint8(ki67_processed(:,:,p)),ki67_processed_name, 'WriteMode','append');
   %end   

cd(matlab_folder)

clc

display(['Pre-processing completed: ' shortfile])

end
