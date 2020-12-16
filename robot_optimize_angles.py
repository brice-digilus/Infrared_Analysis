# -*- coding: utf-8 -*-
"""
Generate angle list and plot numbers of sensors for each angle on a cubemap
Allow to manually tune overlap for best coverage versus number of measurements

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

import matplotlib
#Adjust your backend here
#matplotlib.use('qt5agg')
import matplotlib.pyplot as plt


import numpy as np

import transf_mat


MLX_SENSORS_XPIX = 32
MLX_SENSORS_YPIX = 24


#================= CUBEMAP DATA =====================
cub_npix = 256
p_cam_data = transf_mat.gen_cubemap_data(cub_npix)


#================= WHO TO OPTIMIZE FOR ===============
#current option "LIDAR", "calibration", "dynamixel1", "dynamixel2"
optimize_for = "dynamixel1"

#===================================== Let's alanlyze lidar ==================
#=============== This part is still wip ======================================
if optimize_for == "LIDAR":

    L_SENSORS_XPIX = 100
    L_SENSORS_YPIX = 100
    
    sen_deg = 3
    
    
    
    p_sen_data = {}
    i_sen = 0
    
    nsen = len(p_sen_data.keys())
    for i_row in range(int(90/sen_deg)):
        beta = i_row * sen_deg 
        corr = np.cos(beta*np.pi/180)
        dalpha = sen_deg/corr
        for i_col in range(int(180/dalpha)+1):
            
            if i_row %2:
                alpha = 180 - dalpha*i_col
            else:
                alpha = dalpha*i_col
            p_sen_data[i_sen + nsen] =  {
                    "alpha" : alpha,
                    "beta_i" : beta,
                }
            i_sen += 1
    
    for i_row in range(int(90/sen_deg)):
        beta = 90+(i_row+1) * sen_deg 
        corr = np.abs(np.cos(beta*np.pi/180))
        dalpha = sen_deg/corr
        
        for i_col in range(int(180/dalpha)+1):
            
            if i_row %2:
                alpha = 180 - dalpha*i_col
            else:
                alpha = dalpha*i_col
            p_sen_data[i_sen + nsen] =  {
                    "alpha" : alpha,
                    "beta_i" : beta,
                }
            i_sen += 1
    
    
    
    for sen in p_sen_data:
        p_sen_data[sen]["HFOV"]  = sen_deg
        p_sen_data[sen]["VFOV"]  = sen_deg
        p_sen_data[sen]["HNPIX"] = L_SENSORS_XPIX
        p_sen_data[sen]["VNPIX"] = L_SENSORS_YPIX
    
    
    p_sen_data_final = {"sen" : p_sen_data,
                         "arranged":p_sen_data.keys()}
                     
              


#=============DYNAMIXEL 360 in alpha (bottom servo, rotation about vertical axis) 180 Beta =======================
        
if optimize_for == "dynamixel1":                  


    p_sen_data = {}
    
    list_xx_bi=[[9,-32,20],[9,-4,0],[9,23,20],[7,51,0],[3,76,0]]
    
    
    for xx,betai,off in list_xx_bi:
        nsen = len(p_sen_data.keys())
        for i_sen in range(xx):
            p_sen_data[i_sen + nsen] =  { "alpha" : (181+i_sen*(360/xx))%360+off,  "beta_i" : betai }
    
    #Potentially we can look more down if we want
        
    
    for sen in p_sen_data:
        p_sen_data[sen]["HFOV"]   = 45.64884636 
        p_sen_data[sen]["VFOV"]   = 33.6037171
        p_sen_data[sen]["gamma_i"] = 7.36132776
        p_sen_data[sen]["distortion"]  = -2.87108226
        p_sen_data[sen]["HNPIX"] = MLX_SENSORS_XPIX
        p_sen_data[sen]["VNPIX"] = MLX_SENSORS_YPIX
    
    
    p_sen_data_final = {"sen" : p_sen_data,
                      "arranged": [0,1,2,3,4,5,6,7,8,17,16,15,14,13,12,11,10,9,18,19,20,21,22,23,24,25,26,33,32,31,30,29,28,27,34,35,36]}


#=============DYNAMIXEL 360 in alpha doubled for overlapping images (bottom servo, rotation about vertical axis) 180 Beta =======================

if optimize_for == "dynamixel2":

    p_sen_data = {}
    
    list_xx_bi=[[9,-32,20],[9,-4,0],[9,23,20],[7,51,0],[3,80,0],
                [9,-32-9.25,32.5],[9,-4-10.25,22.5],[9,23-10.25,27.5],[8,51-10.25,18.5],[5,76-9.25,28.5]]
    
    
    for xx,betai,off in list_xx_bi:
        nsen = len(p_sen_data.keys())
        for i_sen in range(xx):
            p_sen_data[i_sen + nsen] =  { "alpha" : (181+i_sen*(360/xx))%360+off,  "beta_i" : betai }
    
    #Potentially we can look more down if we want
        
    
    for sen in p_sen_data:
        p_sen_data[sen]["HFOV"]   = 45.64884636 
        p_sen_data[sen]["VFOV"]   = 33.6037171
        p_sen_data[sen]["gamma_i"] = 7.36132776
        p_sen_data[sen]["distortion"]  = -2.87108226
        p_sen_data[sen]["HNPIX"] = MLX_SENSORS_XPIX
        p_sen_data[sen]["VNPIX"] = MLX_SENSORS_YPIX
    
    
    p_sen_data_final = {"sen" : p_sen_data,
                      "arranged": [0,1,2,3,4,5,6,7,8,17,16,15,14,13,12,11,10,9,18,19,20,21,22,23,24,25,26,33,32,31,30,29,28,27,34,35,36,\
                                   37,38,39,40,41,42,43,44,45, 54,53,52,51,50,49,48,47,46, 55,56,57,58,59,60,61,62,63, 71,70,69,68,67,66,65,64, 72,73,74,75]}




#============= DYNAMIXEL calibration =======================
        
                      
if optimize_for == "calibration":
    
    p_sen_data = {}
    
    hcenter = -44
    vcenter = -1.5
    
    dh = 40
    dv = 25
    #step = 3 #crude calibration
    step = 1 #slow calibration
    
    i_sen = 0
    
    for ah in np.arange(hcenter-dh/2.,hcenter + dh/2., step):
        for av in np.arange(vcenter-dv/2.,vcenter + dv/2., step):
            p_sen_data[i_sen] =  {
                        "alpha" : ah,
                        "beta_i" : av,
                        }   
            i_sen +=1
        
    
    for sen in p_sen_data:
        p_sen_data[sen]["HFOV"]  = 55
        p_sen_data[sen]["VFOV"]  = 35
        p_sen_data[sen]["HNPIX"] = MLX_SENSORS_XPIX
        p_sen_data[sen]["VNPIX"] = MLX_SENSORS_YPIX
        p_sen_data[sen]["Hoff"] = 45 #Following the adjustment
        p_sen_data[sen]["Voff"] = 0
    
    
    p_sen_data_final = {"sen" : p_sen_data,
                      "arranged": list(range(len(p_sen_data.keys())))}
    


                  
#===================================================================================
import math                  
if __name__ == "__main__":
    print("Hi let's try to optimize sensor orientation ")
    
    p_data = p_sen_data_final
    p_sen_data = p_data["sen"]

    transf_mat.compute_all_matrices(p_cam_data,p_sen_data,plot = True)
    
    print("Here is your angle list, you have %d positions " % len(p_sen_data.keys()))  
    for i_sen in p_sen_data.keys():
        print("Sensor %3d\t Alpha %5.1f\t Beta_i %5.1f"%(i_sen,p_sen_data[i_sen]["alpha"],p_sen_data[i_sen]["beta_i"]))


    if "arranged" in p_data.keys():
        i_sen_l = p_data["arranged"]
        print("\n\n\nHere is your arranged list by hand !!! , you have %d positions " % len(i_sen_l))  
        pa = 0
        pb = 0
        for i_sen in i_sen_l:
            print("Sensor %3d\t Alpha   %5.1f\t Beta_i   %5.1f \t DAlpha   %6.1f\t DBeta_i   %5.1f" 
                  %(i_sen,p_sen_data[i_sen]["alpha"],p_sen_data[i_sen]["beta_i"],
                    math.fmod(p_sen_data[i_sen]["alpha"]-pa,360),math.fmod(p_sen_data[i_sen]["beta_i"]-pb,360)))
            pa = p_sen_data[i_sen]["alpha"]
            pb = p_sen_data[i_sen]["beta_i"]
        print("# ==== Angles for angles.txt =========")
        dstr = ""
        #Expect a pair of Y angle, Xi angle per line. eg "-20 12.34"
        for i_sen in i_sen_l:
            dstr += "%6.2f " %(p_sen_data[i_sen]["alpha"])
            dstr += " %6.2f" %(p_sen_data[i_sen]["beta_i"])
            dstr += "\n"
        print(dstr)

    

        
        
        
       
