# -*- coding: utf-8 -*-
"""
This file contains rotation matrices, quaternions and stuff to transfer coordinates from a camera to another

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

import numpy as np
import plotting


# ============= QUATERNIONS ================

#https://stackoverflow.com/questions/39000758/how-to-multiply-two-quaternions-by-python-or-numpy
#https://en.wikipedia.org/wiki/Quaternion#Hamilton_product
def quat_mult(quaternion1, quaternion2):
    """multiply quaternion1*quaternion2
    The product of two rotation quaternions[24] will be equivalent to the rotation a2 + b2i + c2j + d2k followed by the rotation a1 + b1i + c1j + d1k.
    """
    a1, b1, c1, d1 = quaternion1
    a2, b2, c2, d2 = quaternion2
    return np.array([-b2 * b1 - c2 * c1 - d2 * d1 + a2 * a1,
                     b2 * a1 - c2 * d1 + d2 * c1 + a2 * b1,
                     b2 * d1 + c2 * a1 - d2 * b1 + a2 * c1,
                     -b2 * c1 + c2 * b1 + d2 * a1 + a2 * d1], dtype=np.float64)


def quat_normz(quat):
    """Normalize a quaternion with the Z value """
    quat = np.array(quat)
    norm = np.sum(quat[1:4]**2)
    if norm<1:
        quat[0] = np.sqrt(1-norm)
    else:
        quat /= np.sqrt(norm)
        quat[0] = 0
    return quat


def quat_norm(quat):
    """Normalize a quaternion"""
    quat = np.array(quat)
    norm = np.sum(quat**2)
    if norm>1e-6:
        quat /= np.sqrt(norm)
    else:
        quat = np.array([1,0,0,0])
        #TODO: raise warning
    return quat

def quat_from_axis_angle(axis,angle):
    """Create a quaternion from a rotation axis and an angle"""
    axis = np.array(axis)
    axis /= np.sqrt(np.sum(axis**2))
    quat = np.zeros(4)
    quat[1:4] = axis * np.sin(angle/2)
    quat[1] = np.cos(angle/2)
    return quat
    
# If we want to convert from a rotation matrix to a quaternion
    #Bar-Itzhack, Itzhack Y. (2000), “New method for extracting the quaternion from a rotation matrix”, AIAA Journal of Guidance, Control and Dynamics 23(6):1085-1087 (Engineering Note), ISSN 0731-5090

def quat_inv(quat):
    #http://graphics.stanford.edu/courses/cs348a-17-winter/Papers/quaternion.pdf
    """Invert a quaternion """
    qinv = np.array(quat)
    qinv[1:4] = -qinv[1:4]
    #if we do not have an unit quat
    if np.abs(np.sum(quat**2)-1)>(1e-3):
        qinv /= np.sum(quat**2)
    return qinv

def quat_rot_vec(quat,vec):
    """Multiply a quaternion by a vector to do a rotation"""
    vec = np.array(vec)
    quat = np.array(quat)
    qvec = np.zeros(4)
    qvec[1:4] = vec
    return quat_mult(quat_inv(quat),quat_mult(qvec, quat))[1:4]
    
# ============== HOMOGRAPHY ================
    
def compute_cam_matrix(cam_data, invert = False):
    """ Compute homography cam matrix """
    HFOV = cam_data["HFOV"]*np.pi/180
    VFOV = cam_data["VFOV"]*np.pi/180
    fx = 1/np.tan(HFOV/2)
    fy = 1/np.tan(VFOV/2)
    if invert:
        fx = 1/fx
        fy = 1/fy
    return(np.array([[fx,0,0],[0,fy,0],[0,0,1]]))

def compute_rot_matrix(axis,theta, invert = False):
    """ Compute rotation matrix """
    
    if invert:
        theta = -theta
    l = axis[0]
    m = axis[1]
    n = axis[2]

    if np.abs(np.sqrt(l*l+m*m+n*n)-1)>1e-7:
        raise Exception("Please normalize your vector !!!" )


    c = np.cos(theta)
    s = np.sin(theta)
    
    A0 = l*l*(1-c) + c
    A1 = m*l*(1-c) - n*s
    A2 = n*l*(1-c) + m*s

    B0 = l*m*(1-c) + n*s
    B1 = m*m*(1-c) + c
    B2 = n*m*(1-c) - l*s

    C0 = l*n*(1-c) - m*s
    C1 = m*n*(1-c) + l*s
    C2 = n*n*(1-c) + c
    
    return np.array([[A0,A1,A2],
                    [B0,B1,B2],
                    [C0,C1,C2]])

def print_33(mat):
    """Display nicely a 3 by 3 matrix with the numbers well aligned"""
    for yy in range(3):
        print("\t", end = "")
        for xx in range(3):
            if mat[yy,xx] >= 0:
                print(" ", end = "")
            print("%.1f " % mat[yy,xx], end = "")
        print("")
        




def compute_rot_mat(cam, invert = False):
    """This function compute a rotation matrix for a given camera with the parameters from the cam dict
        Expected 
            alpha: vertical rotation
            beta XOR beta_i: horizontal rotation
        Optional
            gamma OR gamma_i: gamma_i is the rotation about the camera axis, gamma is the same axis when alpha = 0
    
    """
    #don't know why I need a minus sign on X maybe that is the camera axis versus the world axis, TODO: Think about it
    if "beta" in cam.keys(): #extrinsinc rotation
        beta = cam["beta"]
        axis_beta = [-1,0,0]
    elif "beta_i" in cam.keys(): #intrinsic rotation
        beta = cam["beta_i"]
        ca = np.cos(cam["alpha"]*np.pi/180)
        sa = np.sin(cam["alpha"]*np.pi/180)
        axis_beta = [-ca,0,-sa] #axis of the second rotation AFTER the first rotation
    else: #weird
        print(" You are weird, so I put 0 as beta (X) rotation" )
        beta = 0
        axis_beta = [-1,0,0]
    

    gamma_i = cam.get("gamma_i",0)
    Rk = compute_rot_matrix([0,0,1],gamma_i*np.pi/180, invert = invert)                        #Intrinsic rotation about the sensor/last axis (Z)
    Rk = np.matmul(Rk,compute_rot_matrix([0,1,0],cam["alpha"]*np.pi/180, invert = invert))     #Rotation about Y (vertical axis)
    Rk = np.matmul(Rk,compute_rot_matrix(axis_beta,beta*np.pi/180, invert = invert))           #Rotation about X or Xi (Horizontal axis)
    gamma = cam.get("gamma",0)
    Rk = np.matmul(Rk,compute_rot_matrix([0,0,1],gamma*np.pi/180, invert = invert))            #Extrinsic rotation about the front axis (Z)   
    
    return Rk




def compute_transfer_mat(cam_src,cam_dst, debug = False):
    """Compute the transfer matrix. Do not normalize coords, do not do the homography
    Here src is generally the sensor point of view and dst is the camera/cubemap point of view
    """

    Cm_src = compute_cam_matrix(cam_src)                            #Camera matrix of the source camera
    Rk = compute_rot_mat(cam_src)                                   #Rotation matrix of the source camera
    Rl_i = compute_rot_mat(cam_dst, invert = True)                  #Inverse of the Rotation matrix of the destination camera
    Cm_dst_i = compute_cam_matrix(cam_dst, invert = True)           #Inverse of the Camera matrix of the destination camera

    R = np.matmul(Cm_src,np.matmul(Rk,np.matmul(Rl_i,Cm_dst_i)))    #Now we compute this nice product between all this and we should be able to transfer from one camera to the other
    if debug:
        print("Camera matrix of the source camera  ")
        print_33(Cm_src)
        print("Rotation matrix of the source camera  ")
        print_33(Rk)
        print("Inverse of the Rotation matrix of the destination camera")
        print_33(Rl_i)
        print("Inverse of the Camera matrix of the destination camera")
        print_33(Cm_dst_i)
        print("Final transfer matrix   ")
        print_33(R)
    return R



def compute_inv_distortion_px(distortion,x,y,xpx,ypx):
    #https://stackoverflow.com/questions/6199636/formulas-for-barrel-pincushion-distortion
    #radius computed in pixels as the fit is in pixels
    #Here the distortion is symmetric, very marginal improvement if it is not
    r2 = np.power(x*((xpx-1)/2.),2)+np.power(y*((ypx-1)/2.),2)
    
    k = distortion/10000 #to be consistent with the calibration
    
    fx = 1 + k*r2

    
    #We limit the range of distortion to reasonable bounds as distortion is inverted on unmasked coordinates, r2 can become very high
    #as r2 is very high we can get fx negative and close to one, making the program believe it is valid pixels
    fx[fx<0.5] = 0.5
    x = x/fx
    y = y/fx
    
    return (x,y)





def compute_pixel_mat(cubemap_data,sensor_data, debug = False):
    """Precompute for each cubemap where to get the info from which sensor
    At the end 5 cubemaps * N sensors * 2 coords = 10*N matrices
    
    At each position of the cubemap, we will not the theoretical coordinate of the sensor image, -1 if this place is not covered
    
        Then one will do linear interpolation or any other algorithm between two pix int(val) and int(val) + 1. 
            This is done in another function
    
        NOTE: check if int val + 1 is not out of bounds for the very particular case where int(val) is at the edge
    """
    
    #Get the transfer matrix between the sensor and the cubemap
    M = compute_transfer_mat(cam_src = sensor_data,cam_dst = cubemap_data, debug = debug)
    

    #Initialize the destination arrays containing the sensor X and Y coordinate for this point of the cubemap
    # shape is V,H because images are stored Y,X
    cmap_sen_x = np.ones([cubemap_data['VNPIX'],cubemap_data['HNPIX']]).astype(np.float)*(-1)
    cmap_sen_y = np.ones([cubemap_data['VNPIX'],cubemap_data['HNPIX']]).astype(np.float)*(-1)

    #Vectorized version harder to read but fast
    #Here we apply the transfer matrix onto the cubemap coordinates to get the sensor coordinates
    (hcub,vcub)= np.meshgrid(range(cubemap_data['HNPIX']),range(cubemap_data['VNPIX']))
    #The matrices work with normalized coordinates, because of how the FOV is defined (as pixel independant, full sensor FOV)
    #So we normalize the cubemap coordinates
    hcubnorm = (hcub - cubemap_data["HNPIX"]/2. + 0.5)/(cubemap_data["HNPIX"]/2.)
    vcubnorm = (vcub - cubemap_data["VNPIX"]/2. + 0.5)/(cubemap_data["VNPIX"]/2.)
    #using tensordot, found the right axes = 1 by trying, not thinking !!!
    #https://docs.scipy.org/doc/numpy/reference/generated/numpy.tensordot.html#numpy.tensordot
    #basically what does this do is multiply the matrix M by each coordinate vector
    vout = np.tensordot(M,[hcubnorm,vcubnorm,np.ones(vcubnorm.shape)],axes = 1)
    
    #We get 3D coordinates, we have to apply the homography transformation
    vout[0,:,:] /= vout[2,:,:]
    vout[1,:,:] /= vout[2,:,:]

    #We invert the distortion on the centered normalized coordinates
    distortion = sensor_data.get("distortion",0)
    (vout[0,:,:],vout[1,:,:]) = compute_inv_distortion_px(distortion, vout[0,:,:], vout[1,:,:],sensor_data["HNPIX"],sensor_data["VNPIX"])


    #We check if the pixel is visible by the sensor, are we withing the sensor boundaries (-1 to 1) and are we looking in the same direction as the sensor (Z>0)
    mask = (vout[0,:,:]>=-1) * (vout[0,:,:] < 1) * (vout[1,:,:]>=-1) * (vout[1,:,:] < 1) * (vout[2,:,:] > 0)
    
    #We transform the normalized coordinates into pixel coordinates
    cmap_sen_x[mask] = vout[0,mask]*((sensor_data["HNPIX"]-1)/2.) + (sensor_data["HNPIX"]/2. - 0.5)
    cmap_sen_y[mask] = vout[1,mask]*((sensor_data["VNPIX"]-1)/2.) + (sensor_data["VNPIX"]/2. - 0.5)

    #We put NaNs when the cubemap pixel is not visible by the sensor
    #Here we basically merge both axis conditions, careful on the order
    not_visible = (cmap_sen_x<0) + (cmap_sen_y<0) #boolean OR
    cmap_sen_x[not_visible] = np.nan
    cmap_sen_y[not_visible] = np.nan

    
    return (cmap_sen_x,cmap_sen_y)



def compute_all_matrices(p_cam_data,p_sen_data,plot = False,plot_sen = False, debug = False):
    """Loop over all cubemap orientation and all sensors to get the transfer matrices
    
    In short, compute the whole 5 cubemaps * 2 coordinates * N sensor transfer matrices

    for cub_i in ["F","R","L","B","T"]
        store the coordinates of the sensor for each cubemap pixel
            in p_cam_data[cub_i]["sen_m"][sen_i]["x"] and p_cam_data[cub_i]["sen_m"][sen_i]["y"]
        store the number of sensor visibles at each cubemap pixel in 
            p_cam_data[cub_i]["sen_nm"]
        
    """
    print("Computing the transfer matrices to map the sensor data onto a cubemap")
    for cub_i in p_cam_data.keys():
        print("  Cube map ",cub_i)
        print("      perspective data ",p_cam_data[cub_i])
        sen_m = p_cam_data[cub_i]["sen_m"] = {}
        #store how many sensors are visible for each pixel of the cubemap
        p_cam_data[cub_i]["sen_nm"] = np.zeros([p_cam_data[cub_i]['VNPIX'],p_cam_data[cub_i]['HNPIX']]).astype(np.int)
        nsensor = len(p_sen_data.keys())
        for sen_i in p_sen_data.keys():
            print("\t\tCub ",cub_i,"   Persp for sensor ",sen_i,"  perspective data ",p_sen_data[sen_i])
            sen_m[sen_i] = {}
            (sen_m[sen_i]['x'], sen_m[sen_i]['y']) = compute_pixel_mat(p_cam_data[cub_i],p_sen_data[sen_i], debug = debug)
            p_cam_data[cub_i]["sen_nm"] += (sen_m[sen_i]['x'] >= 0).astype(int)


    #plot the cubemap for the number of sensors to get coverage            
    if plot:
        fidx = 0
        cmap_d = {face:p_cam_data[face]["sen_nm"] for face in p_cam_data.keys()}
        plotting.plot_cubemap(cmap_d,fidx=fidx, title = " N sen Cube maps ", compute_percentile = False, clim = [0,None])

    #plot the cubemap for one sensor            
    if plot_sen:
        fidx = 10
        for sen_i in range(nsensor):
            for axis in ['x','y']:
                cmap_d = {face:p_cam_data[face]["sen_m"][sen_i][axis] for face in p_cam_data.keys()}
                plotting.plot_cubemap(cmap_d,fidx=fidx,  title = " XY Cmap sensor " +str(sen_i) + "  axis " + axis, compute_percentile = False)
                fidx += 1
        


def gen_cubemap_data(cub_npix = 128):
    """Generate the camera for the cubemap views """
    #================= CAMERA CUBEMAP DATA =====================
    #The wide sensor is about 3 degrees per pixel
    #If we put the cubemap at 64 px it will be 1.4 deg/px
    #theta phi in spherical, alpha beta in Y X extrinsic rotations (in this order) Y up X to the right and Z negative relative to front, seems that X is reverted somewhat !
    p_cam_data = {
      "F": {   "alpha" : 0,
                "beta" : 0,},
      "L": {"alpha" : 90,
                "beta" : 0,},
      "R": {"alpha" : 270,
                "beta" : 0,},
      "B": {"alpha" : 180,
                "beta" : 0,},
      "T": {"alpha" : 0,
                "beta" : 90,},
      "Bo": {"alpha" : 0,
                "beta" : -90,},
    }
    
    
    for cam in p_cam_data:
        p_cam_data[cam]["HFOV"] = 90
        p_cam_data[cam]["VFOV"] = 90
        p_cam_data[cam]["HNPIX"] = cub_npix
        p_cam_data[cam]["VNPIX"] = cub_npix
    return p_cam_data


#======================= Cubemap computation ============================
# Here will be 
#    flatfield correction
#    super resolution
#    deconvolution
#    outlier median filtering
#    and all the fun stuff
# For the moment, "boring" mapping with bilinear interpolation + feathering of overlaps
    
#================= Compute cubemaps ============== 
    

def compute_dist_weight_map(HPIX,VPIX):
    """ Compute a map of pixel weight depending on the position in the sensor.
    This is computed once and then looked up.
    This code is used to feather overlapping sensors where the weight will be bigger closer to the center so we favor close to center information when overlapping.

    NOTE: we could also just write a function to be more continuous
    """
    #for more inspiration about feathering one will want to  google "panorama image overlap distance center stitching feathering" 
    #and read section 6.2 of http://www.cs.toronto.edu/~kyros/courses/2530/papers/Lecture-14/Szeliski2006.pdf
    # we use the distance to the nearest invalid pixel, which here is the distance to the edge as a weight


    (xx,yy)= np.meshgrid(range(HPIX),range(VPIX))
    dist = np.min((np.min((xx,HPIX-1-xx),axis = 0),np.min((yy,VPIX-1-yy),axis = 0)),axis = 0)+1
    
    return dist



def bilinear(data,x_s,y_s):
    """Bilinear interpolation"""
    #https://en.wikipedia.org/wiki/Bilinear_interpolation
    x1 = np.floor(x_s).astype(int)
    x2 = x1+1
    y1 = np.floor(y_s).astype(int)
    y2 = y1+1
    f11 = data[y1,x1]
    f12 = data[y2,x1]
    f21 = data[y1,x2]
    f22 = data[y2,x2]
    r = (f11*(x2-x_s)*(y2-y_s)+f21*(x_s-x1)*(y2-y_s)+f12*(x2-x_s)*(y_s-y1)+f22*(x_s-x1)*(y_s-y1))
    return r


def spacial_correction(r,x_s,y_s,s_info):
    """Compute a simple spatial correction in cos()**4
    Basically vigetting
    
    NOTE: This is a very simple vignetting model, for infrared data the sensor lens holder itself is also an emitter... """
    
    if s_info.get("spacial_corr",0):
            
        HFOV = s_info["HFOV"]*np.pi/180
        VFOV = s_info["VFOV"]*np.pi/180
        corr = s_info["spacial_corr"]
        dxa = (x_s - (s_info["HNPIX"]/2-0.5))/s_info["HNPIX"]*HFOV
        dya = (y_s - (s_info["VNPIX"]/2-0.5))/s_info["VNPIX"]*VFOV
        da = np.sqrt(dxa**2+dya**2)
        invfalloff = 1 - corr + corr * np.cos(da)**4
        
        return r * invfalloff
    else:
        return r


def compute_cubemap(p_cam_data,sensor_data,sensor_frames, feathering = True):
    """ Compute the cubemap from a dict of sensor frames and the camera data.
    The camera data should already have been filled with the transformation matrices compute_all_matrices
    This function will map the sensor values to the cubemap parametrized in p_cam_data
    Today we use a bilinear interpolation and we average overlapping sensors
    
    """
    
    if feathering:
        weight_map = compute_dist_weight_map(sensor_frames[0].shape[1],sensor_frames[0].shape[0])
    
    cubemap_frames = {}
    for cub_i in p_cam_data.keys():
        cubemap_frames[cub_i] = c_frame = np.zeros([p_cam_data[cub_i]['VNPIX'],p_cam_data[cub_i]['HNPIX']])
        c_weight = np.zeros([p_cam_data[cub_i]['VNPIX'],p_cam_data[cub_i]['HNPIX']])


        #vectorized version of the bilinear interp
        #For each sensor we get the data and accumulate
        for sen_i,sen_frame in sensor_frames.items():
            x_s0 = p_cam_data[cub_i]["sen_m"][sen_i]['x']
            y_s0 = p_cam_data[cub_i]["sen_m"][sen_i]['y']
            #we take only pixel with valid values,
            #the previous function (compute_all_matrices) have put to nan every point where there is not sensor visible at this
            #cubemap poisition
            mask = np.isfinite(x_s0)
            x_s = x_s0[mask]
            y_s = y_s0[mask]

            #here we can make a distance dependant overlay
            #We multiply the values by distance weight that was precomputed as a LUT (less continuous but simpler)
            #Then we accumulate the weight at each position and divide by the sum after all sensors
            if feathering:
                w_r = bilinear(weight_map,x_s,y_s)
            else: #We just have uniform weights if there is no feathering
                w_r = np.ones(y_s.shape)
            #We accumulate the weights
            c_weight[mask] += w_r
            #and we add the sensor values with their associated weights
            r = bilinear(sen_frame,x_s,y_s)
            #Spatial correction, will be done only if the value exist in the sensor data
            r = spacial_correction(r,x_s,y_s,sensor_data)
            c_frame[mask] += r*w_r

        #Now as we have done all the sensors we
        #normalize by the accumulated weights
        c_frame[c_weight > 0] = c_frame[c_weight > 0]/c_weight[c_weight > 0]
        #and put Nan for pixels without sensors
        c_frame[c_weight == 0] = np.nan

    return cubemap_frames









