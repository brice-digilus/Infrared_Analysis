# -*- coding: utf-8 -*-
"""
Calibrate the geometric parameters of the camera

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
import matplotlib.pyplot as plt
import numpy as np
import file_loader

# "Low res"
fpath = "data_examples\\Calib_126_img.tiff"

## "high res"
fpath = "data_examples\\Calib_260_img.tiff"


def plot_images(sen_data_all):
    
    max_img = np.max(np.array([sen_data_all['framedata'][i] for i in sen_data_all['framedata']]),axis = 0)
    plt.figure(1)
    plt.imshow(sen_data_all['framedata'][0], cmap = "gray",aspect = "equal")
    plt.figure(2  )
    plt.imshow(max_img, cmap = "gray",aspect = "equal")
    
def plot_crosses(cross_list, fidx = 1,color = "blue"):
    plt.figure(fidx)
    for cross in cross_list:
        plt.scatter(cross[0],cross[1],marker = '+', c = color)

def plot_crosses_pairs(cross_list0,cross_list1, fidx = 1, clear = False, grid = True):
    plt.figure(fidx)
    if clear:
        plt.clf()
    for cross0,cross1 in zip(cross_list0,cross_list1):
        plt.scatter(cross0[0],cross0[1],marker = '+', c = "red")
        plt.scatter(cross1[0],cross1[1],marker = '+', c = "green")
        plt.plot([cross0[0],cross1[0]],[cross0[1],cross1[1]], ":y")
    for x in range(32):
        plt.plot([x+0.5,x+0.5],[0,23], ":k",linewidth=0.5)
    for y in range(24):
        plt.plot([0,31],[y+0.5,y+0.5], ":k",linewidth=0.5)
    plt.axis('equal')
        


def plot_deltas(deltas, fidx = 1):
    plt.figure(fidx)

    for delta in deltas:
        plt.scatter(delta[0],delta[1],marker = '+', c = "red")
    plt.gca().grid(True)
    plt.plot([-0.5,-0.5],[-0.5,0.5], ":k",linewidth=1.5)
    plt.plot([0.5,0.5],[-0.5,0.5], ":k",linewidth=1.5)
    plt.plot([0.5,-0.5],[-0.5,-0.5], ":k",linewidth=1.5)
    plt.plot([0.5,-0.5],[0.5,0.5], ":k",linewidth=1.5)
    plt.axis('equal')


def sensors_find_maxima(sen_data_all):
    """Find the pixel position of the maxima
    It mask the image around the maximum found,
    then use the first moment to get subpixel resolution"""
    deltapx = 3
    sen_data_all['Tmax_coords'] = []
    imshape = sen_data_all['framedata'][0].shape
    gg = np.mgrid[0:imshape[0],0:imshape[1]]
    for ii,image in sen_data_all['framedata'].items():
        initial_guess = np.unravel_index(np.argmax(image, axis=None), image.shape)
        mask = (np.sqrt((gg[0]-initial_guess[0])**2 + (np.abs(gg[1]-initial_guess[1])**2))<deltapx)
        im_masked = image[mask]
        im_masked = im_masked - np.min(im_masked)
        gg_masked = gg[:,mask]
        sen_data_all['Tmax_coords'].append([np.sum(im_masked*gg_masked[0]/np.sum(im_masked)),np.sum(im_masked*gg_masked[1])/np.sum(im_masked)][::-1])

    return sen_data_all



def camera_distortion(x,y,a1):
    r2 = x**2+y**2
    a1 /= 10000
    xd = (1 + a1*r2)*x
    yd = (1 + a1*r2)*y
    return xd,yd



def guess_images_centers_quat(angles_centers_list,par_dh,par_dv,par_dz,par_HFOV,par_VFOV,d_a1):
    HPX = 32
    VPX = 24
    initial_th_guess = []
    #Be very careful on order, way to think about it is rot is qa-1.qb-1.veb.qb.qa so we take vec and undergo qb
    vec = np.array([0,0,1]) #Position of the heat source 
    th = par_dh*np.pi/180
    quatV = np.array([np.cos(th/2),0,1*np.sin(th/2),0]) #vert rotation
    th = par_dv*np.pi/180
    quatH = np.array([np.cos(th/2),1*np.sin(th/2),0,0]) #horiz rotation (up down)
    quatpos = transf_mat.quat_mult(quatV,quatH)
    distance = 0
    deltas = []
    for alpha,beta_i,c_x,c_y in angles_centers_list:
        th = alpha *np.pi/180
        quatalpha = np.array([np.cos(th/2),0,1*np.sin(th/2),0]) #servo v rotation
        th = beta_i *np.pi/180
        quatbeta = np.array([np.cos(th/2),1*np.sin(th/2),0,0]) #servo h rotation
        th = par_dz *np.pi/180
        quatsen = np.array([np.cos(th/2),0,0,1*np.sin(th/2)]) #sensor rotation
        quatservo = transf_mat.quat_mult(transf_mat.quat_mult(quatalpha,quatbeta),quatsen)
        quatrot = transf_mat.quat_mult(quatpos,quatservo)
        vec_rot = transf_mat.quat_rot_vec(quatrot,vec)
        cam_params = {"HFOV":par_HFOV/HPX*2,"VFOV":par_VFOV/VPX*2}#I do not understand this *2
        cam_mat = transf_mat.compute_cam_matrix(cam_params)
        vec_cam = cam_mat.dot(vec_rot)
        px_val = np.array([vec_cam[0]/vec_cam[2],vec_cam[1]/vec_cam[2]])
        
        px_val[0],px_val[1] = camera_distortion(px_val[0],px_val[1],d_a1)

        px_val[0] = -px_val[0] + (HPX/2-0.5) #these signs are still some black magic !
        px_val[1] = px_val[1] + (VPX/2-0.5)
        initial_th_guess.append(px_val)
        distance += np.sqrt((px_val[0]-c_x)**2+(px_val[1]-c_y)**2)
        deltas.append([px_val[0]-c_x,px_val[1]-c_y])

    distance /= len(angles_centers_list)
    print("Average distance",distance)
    return initial_th_guess,distance,deltas

def fitfunc(angles_centers_list,par_dh,par_dv,par_dz,par_HFOV,par_VFOV,d_a1):
    initial_th_guess,distance,d = guess_images_centers_quat(angles_centers_list,par_dh,par_dv,par_dz,par_HFOV,par_VFOV,d_a1)
    return distance


if __name__ == "__main__":


    print("\n==========================================\n    Loading images")
    sen_data_all = file_loader.load_one_file(fpath)


    print("\n==========================================\n    Plotting images")  
    plot_images(sen_data_all)
    print("\n==========================================\n    Finding images maximum positions")
    sensors_find_maxima(sen_data_all)
    print("\n==========================================\n    Initial data guess")
    #store the angles and the maximum positions in a list for the fit, this will be the X vector
    #The Y vector will be the expected distance ie 0
    angles_centers_list = []
    for img_key,frame in sen_data_all["frames"].items():
        angles_centers_list.append([frame['alpha'],frame['beta_i'],\
                                    sen_data_all['Tmax_coords'][img_key][0],sen_data_all['Tmax_coords'][img_key][1]])


    #Compute initial guess of the max positions
    par_dh,par_dv,par_dz,par_HFOV,par_VFOV,d_a1 = (0, 0,0, 55,35.0, 0)               #Blind initial guess
    #Distance before fit: 1.7208454701076064
    #par_dh,par_dv,par_dz,par_HFOV,par_VFOV,d_a1 = (-1.5, 0.71,0, 55,35.0, 0)         #No rotation
    #Distance before fit: 1.4385273309040487
    #par_dh,par_dv,par_dz,par_HFOV,par_VFOV,d_a1 = (-1.5, 0.71 ,7.52, 55,35.0, 0)     #FOV and distortion off for illustration
    #Distance before fit: 0.94678206054299
    #par_dh,par_dv,par_dz,par_HFOV,par_VFOV,d_a1 = (-1.5, 0.71 ,7.52, 46.1,33.70, 0)  #Distortion only off for illustration
    #Distance before fit: 0.34554400892619785
    #Fit result: [-1.56776493  0.61418545  7.44766667 46.10883469 33.76637301 -2.32665921]
    #Average distance 0.1743952202655435
#   historical Fit result: [-0.92263423  1.67898807  7.38385104 45.7282693  33.61861332 -2.77457796]


    initial_th_guess_quat,distance,deltas = guess_images_centers_quat(angles_centers_list,par_dh,par_dv,par_dz,par_HFOV,par_VFOV,d_a1)
    plot_crosses_pairs(initial_th_guess_quat,\
                       sen_data_all['Tmax_coords'],fidx = 2)

    plot_crosses_pairs(initial_th_guess_quat,\
                       sen_data_all['Tmax_coords'],fidx = 4, clear = True)

    plot_deltas(deltas,fidx = 10)
    

    print("Distance before fit:", distance)

    fit = 1
    if fit:
        print("\n==========================================\n    Fitting camera parameters")
        fit_lbounds = [-10,-10,-10,20,10,-10]
        fit_hbounds = [ 10, 10, 10,60,50, 10]

        #prepare fit result vector
        yfit = np.zeros(len(angles_centers_list))
        from scipy.optimize import curve_fit
        
        popt,pcov = curve_fit(fitfunc, angles_centers_list, yfit, p0 = [par_dh,par_dv,par_dz,par_HFOV,par_VFOV,d_a1]\
                              ,ftol = 1e-5,maxfev = 5000, bounds = (fit_lbounds,fit_hbounds))#, method = "trf") 
        print("Fit result:", popt)
        print("Distance before fit:", distance)
    
        par_dh,par_dv,par_dz,par_HFOV,par_VFOV,d_a1 = popt

        fitted_centers_quat,distance,deltas = guess_images_centers_quat(angles_centers_list,par_dh,par_dv,par_dz,par_HFOV,par_VFOV,d_a1)
        plot_crosses_pairs(fitted_centers_quat,\
                           sen_data_all['Tmax_coords'],fidx = 5, clear = True)
    
    
        plot_deltas(deltas,fidx = 11)

    #Go for the fit !
    







