# -*- coding: utf-8 -*-
"""
This program map a single file into a cubemap

@author: Brice Dubost

   Copyright 2020 Brice Dubost

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

import transf_mat
import file_loader
import plotting
import os.path

dfolder = "_out"

sfolder = "data_examples"
fname = "20201217_livingroom.tiff"

fpath = os.path.join(sfolder,fname)


#================= CAMERA CUBEMAP DATA =====================
cub_npix = 256
p_cam_data = transf_mat.gen_cubemap_data(cub_npix)



def file_to_cube(fpath):

    fname = os.path.basename(os.path.normpath(fpath))
    print("\n==========================================\n    Loading %s ..." % fname)

    sen_data = file_loader.load_one_file(fpath)
    p_sen_data = sen_data["frames"]
    



    print("\n==========================================\n    Computing transf mat")
    
    #Compute transfer matrix -- do it once if there are several frames
    transf_mat.compute_all_matrices(p_cam_data,p_sen_data,plot = True)
        
    print("\n==========================================\n    Computing image cubemaps")
    #compute the cubemap
    # len(sen_data["framedata"]) = number of images
    # sen_data["framedata"][0].shape = one image shape (V,H)
    # sen_data["frames"] dict containing angle information
    #Conversion matrix information is in p_cam_data
    sen_data["cubemap"] = transf_mat.compute_cubemap(p_cam_data,sen_data["sp"],sen_data["framedata"])

    #plotting
    # If we want the data back ,savetype = ["plot","return"]
    # if we want to save, we need to set a dstfolder and ,savetype = ["plot","fast"] or ,savetype = "fast"
    plotting.plot_cubemap(sen_data["cubemap"],1000,title = fname) 
    
    plotting.plot_cubemap(sen_data["cubemap"],savetype = "fast_tex",dstfolder = dfolder) 
    
    plotting.plot_cubemap(sen_data["cubemap"],2000,title = fname, triangle_disp = True) 


    
    
    
if __name__ == "__main__":
    file_to_cube(fpath)
    
    
    
    
    
    
    
    
    

