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

#TODO:  Calibration information is hardcoded here should be moved in a configuration file

import time
from PIL import Image
import numpy as np



RAW_min_t = -10
RAW_max_t = 50
RAW_nbits = 16



sensors_data = {}

#Fit results, the key correspond to the model present in the TIFF metadata
sensors_data["THDS"] = {}
sensors_data["THDS"]["gamma_i"]     = 7.36132776
sensors_data["THDS"]["HFOV"]        = 45.64884636
sensors_data["THDS"]["VFOV"]        = 33.6037171
sensors_data["THDS"]["distortion"]  = -2.87108226
sensors_data["THDS"]["Hoff"]        = 45 #Following the adjustment of the rotation plate
#Spatial correction tuned by hand
sensors_data["THDS"]["spacial_corr"]  = 0.12



def get_sensor_data_from_tiff(filename):
    
    taginfo = Image.open(filename).tag_v2.get(270,"").replace("\x00","").split("//")
    
    #if the TIFF is not tagged this is a static sensor
    if not len(taginfo):
        return sensors_data["static"]
    
    thermtags = {"MODL":taginfo[0]}
    #We get the first tag, before the underscore to have the model
    #TODO: fix when we will have several to have a serial ID
    sensor_data = sensors_data[thermtags["MODL"].split("_")[0]]

    #now we get the other tags
    for tag in taginfo[1:]:
        thermtags[tag.split(":")[0]] = tag.split(":")[1]
    
    #getting time information
    sensor_data["starttime_str"] = thermtags["TI_SD"] + " " + thermtags["TI_ST"]
    print(sensor_data["starttime_str"])
    sensor_data["starttime_epoch"] = time.mktime(time.strptime(sensor_data["starttime_str"], '%Y%m%d %H%M%S'))
    sensor_data["duration"] = float(thermtags.get("TI_DU",0))
    
    #getting angles data, as well as ambiant temperature
    alpha_l = np.array(thermtags["ANH"].split(";")).astype(float)
    beta_i_l = np.array(thermtags["ANV"].split(";")).astype(float)
    ta_s = np.array(thermtags["TA_S"].split(";")).astype(float)
    Nframes = len(alpha_l)
    sensor_data["frames"] = {}
    for i in range(Nframes):
        sensor_data["frames"][i] = {"alpha" : alpha_l[i]+sensor_data.get("Hoff",0),
                     "beta_i" : beta_i_l[i]+sensor_data.get("Voff",0),
                     "Ta": ta_s[i] }
    
    return sensor_data
        
        


def tiff_get_deg_image(fpath):
    """Take one TIFF file, convert the data to degrees"""

    pic = np.array(Image.open(fpath))
  
    pix_deg = np.array(pic).astype(np.double)
    pix_deg /= (1<<RAW_nbits)-1
    pix_deg *= (RAW_max_t-RAW_min_t)
    pix_deg += (RAW_min_t)

    return pix_deg




def load_one_file(fpath):
    """Load one sensor file in the data struture used by the analysis
    """
    
    print("  Loading sensor data from %s" % fpath)
    
    sensor_data = {}
    sensor_data["sp"] = get_sensor_data_from_tiff(fpath)
    sensor_data["frames"] = sensor_data["sp"].pop("frames")
    Nframes = len(sensor_data["frames"])
    sensor_data["framedata"] = {}
    #getting degrees C data
    pic_deg = tiff_get_deg_image(fpath)
    #getting frame width and height from image information
    sensor_data["sp"]["VNPIX"] = pic_deg.shape[0]
    if pic_deg.shape[1] % Nframes:
        raise Exception("File %s, the width %d is not a multiple of number of frames %d" %(fpath,pic_deg.shape[1], Nframes))
    sensor_data["sp"]["HNPIX"] = int(pic_deg.shape[1]/Nframes)
    #extracting the sub sensor frames
    for i_frame in range(Nframes):
        sensor_data["framedata"][i_frame] = pic_deg[:,(i_frame*sensor_data["sp"]["HNPIX"]):((i_frame+1)*sensor_data["sp"]["HNPIX"])]

    #THE SENSOR IS UPSIDE DOWN IN THE FRAME --> YES
    #   --> Datasheet says, if pin is at the bottom and if we look toward the front
    #       Top left     corner is Row 1  Col 32,
    #       Bottom left  corner is Row 24 Col 32,
    #       Top right    corner is Row 1  Col 1,
    #       Bottom right corner is Row 24 Col 1,
    #I could expect a flip on Y but there is also a flip on X which is sensor specific
    for i_f in sensor_data["framedata"]:
        sensor_data["framedata"][i_f] = np.fliplr(np.flipud(sensor_data["framedata"][i_f]))

    #We copy the global keys in each frame to make the rest of the processing easier
    copy_key_list = sensor_data["sp"].keys()
    for i,frame in sensor_data["frames"].items():
        for key in copy_key_list:
            frame[key] = sensor_data["sp"][key]
            
    
    
    return sensor_data




















