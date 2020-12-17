# -*- coding: utf-8 -*-
"""
Display the different computation steps for sensor to cubemap

This file will allow to demonstrate the various computation setp associated with a cubemap transfer calculation

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


import transf_mat, plotting


MLX_SENSORS_XPIX = 32
MLX_SENSORS_YPIX = 24

#================= CUBEMAP DATA =====================
cub_npix = 256
p_cam_data = transf_mat.gen_cubemap_data(cub_npix)


#================= We look for only one sensor position ===============
 
# We set up the angle at which we look at       
p_sen_data = {}
p_sen_data[0] = { "alpha" : 0,  "beta_i" : 0 }

#The we add some extra information about the sensor itself
for sen in p_sen_data:
    p_sen_data[sen]["HFOV"]   = 45.64884636 
    p_sen_data[sen]["VFOV"]   = 33.6037171
    p_sen_data[sen]["gamma_i"] = 7.36132776
    p_sen_data[sen]["distortion"]  = -2.87108226
    p_sen_data[sen]["HNPIX"] = MLX_SENSORS_XPIX
    p_sen_data[sen]["VNPIX"] = MLX_SENSORS_YPIX


#================ Now we study compute all matrices ========================
   
if __name__ == "__main__":
    print("Hi let's look at sensor calculation ")
    
    print("Cubemap info: p_cam_data",p_cam_data,"\n")
    
    print("Sensor data: p_sen_data",p_sen_data,"\n")
    
    print("Now we do the overlapt computation \n")
    transf_mat.compute_all_matrices(p_cam_data,p_sen_data,plot = True, plot_sen = True)
    
    print("Now we have two new keys in the cubemap data: sen_m and sen_nm\n")
    print("\t sen_m is indexed by sensors, and contains, for each pixel of the cubemap the X and Y coordinates of this sensor in this position\n")
    print("\t sen_nm contains the number of sensors at each pixel of this cubemap\n")

    print("Cubemap info new elements from compute_all_matrices: p_cam_data['F'].keys()",p_cam_data['F']['sen_m'],"\n")
    print("Cubemap info new elements from compute_all_matrices: p_cam_data['F'].keys()",p_cam_data['F']['sen_nm'],"\n")
    
    print("Sensor data: p_sen_data",p_sen_data,"\n")
    
    plotting.plot_image(p_cam_data['F']['sen_nm'],title = "Number of overlapping sensors for\nFront view: p_cam_data['F']['sen_nm']",fidx = 1)
    plotting.plot_image(p_cam_data['F']['sen_m'][0]['x'],title = "X coordinates of\nsensor 0 in the front view:\n p_cam_data['F']['sen_m'][0]['x']", fidx = 2)
    plotting.plot_image(p_cam_data['F']['sen_m'][0]['y'],title = "Y coordinates of\nsensor 0 in the front view:\n p_cam_data['F']['sen_m'][0]['y']", fidx = 3)

    #============ Let's look a bit closer on what compute_all_matrices is doing ================        

    #the key line in this function is the call to compute_pixel_mat
    #            (sen_m[sen_i]['x'], sen_m[sen_i]['y']) = compute_pixel_mat(p_cam_data[cub_i],p_sen_data[sen_i], debug = debug)


    #First it computes the transfer matrix between the two cameras
    p_cam_data_front = transf_mat.gen_cubemap_data(cub_npix)['F']

    print("Cubemap data (destination camera)    ", p_cam_data_front)
    print("Sensor data  (source camera)         ",p_sen_data[0],"\n")
    M = transf_mat.compute_transfer_mat(cam_src = p_sen_data[0],cam_dst = p_cam_data_front, debug = True)




    

        
        
        
       

