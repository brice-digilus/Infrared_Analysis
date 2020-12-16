# -*- coding: utf-8 -*-
"""
This file contains helper functions to do plots of cubemaps and alike

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


import matplotlib.pyplot as plt
import numpy as np
import os

#================= plotting functions =====================

def plot_image(px,fidx=1,title = "", clim =[],save = False, dstfolder = ""):
    plt.figure(fidx)
    if len(clim) == 2:
        plt.imshow(px, cmap = "plasma",aspect = "equal",vmin = clim[0], vmax = clim[1])
    else:
        plt.imshow(px, cmap = "plasma",aspect = "equal")
    plt.title(title)
    plt.colorbar()
    if save:
        plt.savefig(os.path.join(dstfolder,"___FIG_%d.png"%fidx))
        plt.close(fidx)

def plot_2images(px1,px2,fidx=1,title = ""):
    plt.figure(fidx)
    plt.subplot(211)
    plt.imshow(px1, cmap = "bwr",aspect = "equal")
    plt.colorbar()
    plt.xlim([0,px1.shape[1]])
    plt.ylim([px1.shape[0],0])
    plt.title(title)
    plt.subplot(212)
    plt.imshow(px2, cmap = "bwr",aspect = "equal")
    plt.colorbar()
    plt.xlim([0,px1.shape[1]])
    plt.ylim([px1.shape[0],0])
    plt.title(title)
    

def gray_to_color(px,clim =[0,1]):
    #https://stackoverflow.com/questions/15140072/how-to-map-number-to-color-using-matplotlibs-colormap
    import matplotlib.cm
    px_norm = np.array(px) - clim[0]
    px_norm /= (clim[1]-clim[0])
    px_norm[px_norm<0] = 0
    px_norm[px_norm > 1] = 1
    
    colors = matplotlib.cm.plasma(px_norm.flatten())
    colors = colors.reshape([px_norm.shape[0],px_norm.shape[1],4])
    return colors

def save_image(px,fidx=1,clim =[0,1],title = "", dstfolder = "", fname = None):
    """Fast way to save an image with a colormap"""
    from PIL import Image

    colors = gray_to_color(px,clim = clim)
    if fname is None:
        fname = "___FIG_%d.png"%fidx
    dfname = os.path.join(dstfolder,fname)
    
    Image.fromarray((colors*255).astype(np.uint8)).save(dfname)
    if len(title):
        ff = open(dfname + '.txt','w')
        ff.write(title)
        ff.close()



        
def tri_square(sq,tri,nrot = 0):
    """Helper function for triangular wrapping cubemaps ie quasi isocahedral
    
      /\    /\    /\    /\  
     /__\  /__\  /__\  /__\ 
    |    ||    ||    ||    |
    |____||____||____||____|
     \  /  \  /  \  /  \  / 
      \/    \/    \/    \/  
    """
    dsq = sq
    while nrot < 0:
        nrot += 4
    for i in range(nrot):
        dsq = np.rot90(dsq)
    return tri*dsq + (-999)*(1-tri)


def plot_cubemap(cmap_d,fidx = 1,title = "", clim = None,savetype = "plot", dstfolder = "", compute_percentile = True,triangle_disp = False):
    """
        Plot a cubemap.
        savetype define the output. either a string or a list of string, default value "plot"
            "fast"     is just a file (in the folder dstfolder) computed with the quick colormap algorithm
            "saveplot" is just a file (in the folder dstfolder) computed with the matplotlib colormap (slower)
            "plot"     is just a plot with figure idx fidx
            "return"   will return the merged cubemap
        clim is either None or a list of two elements, if one of the two element is None it will be replaced by an automatic value
            if we compute percentiles we will use the 0.5% and 99.5% percentile to set the clim in order to avoid outliers
            
        triangle disp is a fake cubemap where the top and bottom cubes are cut into triangles for human readability
    
    """
    
    #put cubemaps together
    hh,ww = cmap_d["F"].shape
    dcubemap = {}
    if triangle_disp:
        dhh = int(hh/2)
        xx,yy = np.meshgrid(range(hh),range(hh))
        triangle = ((yy[::-1,:]+np.abs(xx-(hh-1)/2))<(hh/2))*1
        tdown = np.rot90(np.rot90(triangle))
    else:
        dhh = hh
    dcubemap = np.ones([hh+dhh,ww*4])*(-999)  
    
    if "Bo" in cmap_d.keys():
        dcubemap = np.ones([hh+2*dhh,ww*4])*(-999)
        dcubemap[(hh+dhh):(hh+2*dhh),ww:(2*ww)]            = cmap_d["Bo"][0:dhh,0:ww]
        
        if triangle_disp:
            for ir in range(4):
                dcubemap[(hh+dhh):(hh+2*dhh),(ww*ir):((ir+1)*ww)]            = tri_square(cmap_d["Bo"],tdown,ir-1)[0:dhh,:]
        

    dcubemap[0:dhh,ww:(2*ww)]              = cmap_d["T"][0:dhh,0:ww]
    dcubemap[dhh:(hh+dhh),(ww):(2*ww)]     = cmap_d["F"]
    dcubemap[dhh:(hh+dhh),0:(ww)]          = cmap_d["L"]
    dcubemap[dhh:(hh+dhh),(2*ww):(3*ww)]   = cmap_d["R"]
    dcubemap[dhh:(hh+dhh),(3*ww):(4*ww)]   = cmap_d["B"]
    if triangle_disp:
        for ir in range(4):
            dcubemap[0:dhh,(ww*ir):((ir+1)*ww)]            = tri_square(cmap_d["T"],triangle,-ir+1)[dhh:hh,:]


    
    dcubemap[dcubemap<-373] = np.nan

    if clim is None:
        clim = [None,None]
    coooold = np.nanmin(dcubemap)
    hoooot  = np.nanmax(dcubemap)
    title = "Cubemap " + title + "\n min %.1f max %.1f " %(coooold,hoooot) 
    if compute_percentile:
        median = np.nanmedian(dcubemap)
        percentile10 = np.nanpercentile(dcubemap,10)
        percentile90 = np.nanpercentile(dcubemap,90)
        percentile2 = np.nanpercentile(dcubemap,2)
        percentile98 = np.nanpercentile(dcubemap,98)
    
        percentile005 = np.nanpercentile(dcubemap,0.5)
        percentile995 = np.nanpercentile(dcubemap,99.5)

        if clim[0] is None:
            clim[0] = int(percentile005)-1-2
        if clim[1] is None:
            clim[1] = int(np.ceil(percentile995))+1+3
        title = title + "median T %.1f\n percentile 10%% T %.1f 80%% T %.1f\n percentile 2%% T %.1f 98%% T %.1f clim %s" %(median,percentile10,percentile90,percentile2,percentile98,str(clim))
    else:
        if clim[0] is None:
            clim[0] = coooold
        if clim[1] is None:
            clim[1] = hoooot

    #display or save cubemaps        
    if type(savetype) == type("fast"):
        savetype = [savetype]
    if "fast" in savetype:
        save_image(dcubemap,fidx,clim = clim, title = title, dstfolder = dstfolder)
    if "saveplot" in savetype:
        plot_image(dcubemap,fidx,clim = clim,title = title,save = True, dstfolder = dstfolder)
    if "plot" in savetype:
        plot_image(dcubemap,fidx,clim = clim,title = title,save = False, dstfolder = dstfolder)
        
    mapping = {"T": "posy","F":"posz","L":"negx","R":"posx","B":"negz","Bo":"negy"}
    if "fast_tex" in savetype:
        for fid in cmap_d.keys():
            save_image(cmap_d[fid],fidx,clim = clim, title = "", dstfolder = dstfolder, fname = mapping[fid]+ ".png")
        
    if "return" in savetype:
        return dcubemap
        
        







