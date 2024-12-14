## load and process brainreg images with following scripts
## you will need 3 tiff files from brainreg, and the structure csv file for converting label index to strucutre names
## you can find the structure csv here:
## https://github.com/hung-lo/binge_code/blob/main/structures.csv


## first activate your conda env and type "napari" to open the napari window
## then in the napari console (click the most bottom left icon) to open the console window

import tifffile, napari
import pandas as pd
from napari.viewer import current_viewer

mouse_id = 'CHR017962'
# mouse_id = 'DSC024432'
# mouse_id = 'DSC024434'

data_path = f'/Users/hunglo/Documents/inscopix_csv/brainsaw_measures/anterograde_intensity/data/{mouse_id}_OT_eOPN3/'

index_label = 961 ## this is the index for piriform cortex labels

## open files 
## I'm using different channels here since the red/mScarlett channel has more artifacts (eg. shrinking boundary around the injection site).
## using the autofluorescent channel can prevent this

viewer.open(f'{data_path}brainregoutput_green/registered_atlas.tiff', name='registered_atlas')
viewer.open(f'{data_path}brainregoutput_green/registered_hemispheres.tiff', name='registered_hemispheres')
viewer.open(f'{data_path}brainregoutput_red/downsampled.tiff', name='downsampled')

image = viewer.layers['downsampled'].data
label = viewer.layers['registered_atlas'].data
hemi = viewer.layers['registered_hemispheres']

## try several regions together ##This doesn't work yet
# region = label in [36,180,148,187,638,662] 

region = label==index_label

## try to do the calculation without adding to multiple label layers 
## doesn't work yet
# viewer.add_image(image*region*(hemi.data==1),name='left_hemi_crop_img')
# viewer.add_image(image*(region*[hemi.data==2][0].data),name='right_hemi_crop_img')

## select only the 
region = label==index_label
viewer.add_labels(region)

left_hemi = [hemi.data==1][0]
viewer.add_labels(left_hemi,name='left_hemi')
viewer.add_labels(region*left_hemi,name='left_hemi_pc') ## image crop with hemi shape
left_hemi_pc = viewer.layers['left_hemi_pc']
viewer.add_image(image*left_hemi_pc.data,name='left_hemi_pc_img') ## image crop with hemi shape
data_left_pc = viewer.layers['left_hemi_pc_img'].data
right_hemi = [hemi.data==2][0]
viewer.add_labels(right_hemi,name='right_hemi')
viewer.add_labels(region*right_hemi,name='right_hemi_pc') ## image crop with hemi shape
right_hemi_pc = viewer.layers['right_hemi_pc']
viewer.add_image(image*right_hemi_pc.data,name='right_hemi_pc_img') ## image crop with hemi shape
data_right_pc = viewer.layers['right_hemi_pc_img'].data


def label2name(label_list):
    import pandas as pd
    structure_df = pd.read_csv("/Users/hunglo/Documents/inscopix_csv/brainsaw_measures/anterograde_intensity/structures.csv")
    name_list  = [structure_df[structure_df["id"]==i]["name"].values[0] for i in label_list]
    return name_list

# structure_name = 'Piriform area'
structure_name = label2name([index_label])[0]
output_path = '/Users/hunglo/Documents/inscopix_csv/brainsaw_measures/anterograde_intensity/raw_image_cropped'

# for layer in viewer.layers:
for layer in current_viewer().layers:
    print(layer.name)
    if 'img' in layer.name:
        if 'left' in layer.name:
            napari.save_layers(f'{output_path}/{structure_name}_{mouse_id}_left.tif', [layer])
            print(f'saved {structure_name}_left')

        elif 'right' in layer.name:
            napari.save_layers(f'{output_path}/{structure_name}_{mouse_id}_right.tif', [layer])
            print(f'saved {structure_name}_right')

# viewer.close() # use this to close the napari window
print('============done=============')
