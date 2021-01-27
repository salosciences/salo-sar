#!/usr/bin/env python
# coding: utf-8

# # Forest Stand Height Functions

# In[ ]:

import sys
import numpy as np
import math
import subprocess
import scipy.io as sio
import pdb
import scipy
import xml.etree.ElementTree as ET
import os
import time
import json
from osgeo import gdal, osr
from scipy import interpolate
from scipy import stats
import mpmath as mp
import simplekml
from PIL import Image


def CROP_ISCE():
    xmlfile = "resampOnlyImage.amp.xml"
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    size_array = np.array([])
    for size in root.iter('property'):
        if size.items()[0][1] == 'size':
            size_array = np.append(size_array, int(size.find('value').text))
    width = size_array[0]
    length = size_array[1]

    nanval = 0

    # Read amp files in radar coordinates
    amp_file = np.fromfile("resampOnlyImage.amp", dtype='complex64')
    inty = amp_file.reshape((length,width))

    inty[:176,:] = nanval
    inty[5488:,:] = nanval
    inty[:,:163] = nanval
    inty[:,4846:] = nanval


    # Write output files
    inty.tofile("resampOnlyImage.amp")

# In[ ]:

def crop_ROIPAC():
    directory = sys.argv[1]					# Directory where amp/cor files are
    date1 = sys.argv[2]						# SAR date1
    date2 = sys.argv[3]						# SAR date2

    # Extract ROI_PAC parameters
    amp_rsc_file = date1+"-"+date2+"_2rlks.amp.rsc"
    width = int(read_rsc_data(amp_rsc_file, directory, "WIDTH"))
    length = int(read_rsc_data(amp_rsc_file, directory, "FILE_LENGTH"))
    fullwidth = width*2
    nanval = 0

    # Read cor files in radar coordinates
    cor_file = np.fromfile(os.path.join(directory, date1+"-"+date2+"_2rlks.cor"),np.float32, count=length*fullwidth)
    corr = cor_file.reshape((length, fullwidth))
    mag = corr[:,0:width]
    phs = corr[:,width:fullwidth]

    # Read amp files in radar coordinates
    amp_file = np.fromfile(os.path.join(directory, date1+"-"+date2+"_2rlks.amp"), np.complex64)
    inty = amp_file.reshape((length,width))

    # Creating empty array for cropped square list
    mag[:638,:] = nanval
    mag[3288:,:] = nanval
    mag[:,:84] = nanval
    mag[:,2418:] = nanval

    phs[:638,:] = nanval
    phs[3288:,:] = nanval
    phs[:,:84] = nanval
    phs[:,2418:] = nanval

    inty[:638,:] = nanval
    inty[3288:,:] = nanval
    inty[:,:84] = nanval
    inty[:,2418:] = nanval

    # Creating empty array for square list
    c_out = np.zeros((length,fullwidth))

    # Writing vals
    c_out[:,0:width] = mag
    c_out[:,width:fullwidth] = phs

    # Write output files
    cx = c_out.astype('f4')
    cx.tofile(os.path.join(directory, date1+"-"+date2+"_2rlks_fix.cor"))
    inty.tofile(os.path.join(directory, date1+"-"+date2+"_2rlks_fix.amp"))


# In[ ]:

def arc_sinc(x, c_param):
    # Get rid of extreme values by set all values where x > 1 equal to 1, and x < 0 equal to 0 
    x[(x > 1)] = 1
    x[(x < 0)] = 0

    # Create array of increments between 0 and pi of size pi/100
    XX = np.linspace(0, math.pi, num=100, endpoint=True)

    # Set the first value of XX to eps to avoid division by zero issues -> Paul's suggestion
    XX[0] = np.spacing(1)
    
    # Calculate sinc for XX and save it to YY
    YY = np.sin(XX) / XX

    # Reset the first value of XX to zero and the first value of YY to the corresponding output
    XX[0] = 0
    YY[0] = 1
    
    # Set the last value of YY to 0 to avoid NaN issues
    YY[-1] = 0

    # Flip XX and YY left to right
    XX = XX[::-1]
    YY = YY[::-1]
  

    XX = np.array([ 3.14159265,  3.10985939,  3.07812614,  3.04639288,  3.01465962, 2.98292636,  2.9511931 ,  2.91945984,  2.88772658,  2.85599332, 2.82426006,  2.7925268 ,  2.76079354,  2.72906028,  2.69732703, 2.66559377,  2.63386051,  2.60212725,  2.57039399,  2.53866073, 2.50692747,  2.47519421,  2.44346095,  2.41172769,  2.37999443, 2.34826118,  2.31652792,  2.28479466,  2.2530614 ,  2.22132814, 2.18959488,  2.15786162,  2.12612836,  2.0943951 ,  2.06266184, 2.03092858,  1.99919533,  1.96746207,  1.93572881,  1.90399555, 1.87226229,  1.84052903,  1.80879577,  1.77706251,  1.74532925, 1.71359599,  1.68186273,  1.65012947,  1.61839622,  1.58666296, 1.5549297 ,  1.52319644,  1.49146318,  1.45972992,  1.42799666, 1.3962634 ,  1.36453014,  1.33279688,  1.30106362,  1.26933037, 1.23759711,  1.20586385,  1.17413059,  1.14239733,  1.11066407, 1.07893081,  1.04719755,  1.01546429,  0.98373103,  0.95199777, 0.92026451,  0.88853126,  0.856798  ,  0.82506474,  0.79333148, 0.76159822,  0.72986496,  0.6981317 ,  0.66639844,  0.63466518, 0.60293192,  0.57119866,  0.53946541,  0.50773215,  0.47599889, 0.44426563,  0.41253237,  0.38079911,  0.34906585,  0.31733259, 0.28559933,  0.25386607,  0.22213281,  0.19039955,  0.1586663 , 0.12693304,  0.09519978,  0.06346652,  0.03173326,  0.        ])
    YY = np.array([ 0.        ,  0.01020237,  0.02060472,  0.03120282,  0.04199229, 0.05296859,  0.06412703,  0.07546277,  0.08697083,  0.09864608, 0.11048326,  0.12247694,  0.1346216 ,  0.14691157,  0.15934105, 0.17190411,  0.18459472,  0.19740671,  0.21033383,  0.22336969, 0.23650781,  0.24974161,  0.26306441,  0.27646944,  0.28994984, 0.30349868,  0.31710894,  0.33077352,  0.34448527,  0.35823696, 0.37202131,  0.38583098,  0.39965857,  0.41349667,  0.42733779, 0.44117444,  0.45499906,  0.46880411,  0.48258199,  0.49632512, 0.51002589,  0.52367669,  0.53726993,  0.55079798,  0.56425328, 0.57762824,  0.59091533,  0.60410701,  0.6171958 ,  0.63017424, 0.64303494,  0.65577053,  0.66837371,  0.68083722,  0.69315389, 0.7053166 ,  0.7173183 ,  0.72915204,  0.74081093,  0.75228819, 0.76357711,  0.77467109,  0.78556364,  0.79624836,  0.80671897, 0.81696931,  0.82699334,  0.83678514,  0.84633891,  0.85564901, 0.8647099 ,  0.87351622,  0.88206272,  0.89034433,  0.8983561 , 0.90609326,  0.91355119,  0.92072543,  0.92761169,  0.93420585, 0.94050396,  0.94650224,  0.95219709,  0.9575851 ,  0.96266301, 0.96742778,  0.97187655,  0.97600663,  0.97981554,  0.98330097, 0.98646084,  0.98929323,  0.99179643,  0.99396894,  0.99580945, 0.99731683,  0.99849018,  0.9993288 ,  0.99983218,  1.        ])

    # Run interpolation
    # XX and YY are your original values, x is the query values, and y is the interpolated values that correspond to x
    interp_func = interpolate.interp1d(YY, XX * c_param, kind='slinear') 
    y = interp_func(x)

    # Set all values in y less than 0 equal to 0
    y[(y < 0)] = 0
    # return y
    
    return y


# In[ ]:


# Input parameters are the numbers of scenes, edges, start scene, iterations, the input/output file directory, 
    # averaging numbers in lat and lon for "self" and "pairwise" fitting, bin_size for density calculation in scatter plot fitting, 
    # flag for sparse data cloud filtering.
def auto_mosaicking_new(scenes, edges, start_scene, N, linkarray, directory, Nd_pairwise, Nd_self, bin_size, flag_sparse):


    # Set average S and C parameters (0<s<1, 0<c<20 so s=0.65 and c=13)
    avg_S = 0.65
    avg_C = 13
    
    # Create avg_dp matrix, and fill  with average S and C parameters
    avg_dp = np.zeros(scenes * 2)
    np.put(avg_dp, range(0, scenes * 2, 2), avg_S)
    np.put(avg_dp, range(1, scenes * 2, 2), avg_C)
    
    # Create the dp matrix 
    # the difference of the avg and the initial SC values OR all zeros (avg - avg)
    dp = np.zeros(scenes * 2)

    # Initialize target matrix and fill with K=1, B=0
    target_KB = np.zeros((edges + 1) * 2)
    np.put(target_KB, range(0, (edges + 1) * 2, 2), 1)

    # Run cal_KB()
    Y = cal_KB(dp, edges, start_scene, linkarray, directory, Nd_pairwise, Nd_self, bin_size, flag_sparse)

    # Calculate the residual for cal_KB - target
    res = sum((Y - target_KB)**2)

    # Save dp and the residual as the first iteration output file (using JSON)
    iter_file = open(os.path.join(directory, "output/SC_0_iter.json"), 'w')
    json.dump([dp.tolist(), res], iter_file)
    iter_file.close()


    # For the rest of the iterations run ls_deltaSC() and save to output file (using JSON)
    for i in range(1, N + 1, 1): # this will run from i=1 to i=N
        [dp, res] = ls_deltaSC(dp, edges, scenes, start_scene, linkarray, directory, Nd_pairwise, Nd_self, bin_size, flag_sparse)
        print ("%d iterations completed!\n" % i)
        print (time.strftime("%H:%M:%S"))
        filename = "SC_%d_iter.json" % i
        
        iter_file = open(os.path.join(directory, "output/" + filename), 'w')
        json.dump([dp.tolist(), res], iter_file)
        iter_file.close()

    print ("auto_mosaicking_new finished at " + (time.strftime("%H:%M:%S")))


# In[ ]:


def auto_tree_height_many(scenes, flagfile, directory, numLooks, noiselevel, flag_proc, flag_grad):

    # For each scene name the file, run auto_tree_height_single and save the output to a .json file
    for i in range(scenes):

        # Get the scene data and set the file name and image folder name (f#_o# where # is the frame and orbit numbers, respectively)
        scene_data = flag_scene_file(flagfile, i + 1, directory) # 0 vs 1 indexing
        filename = scene_data[1]
        image_folder = "f" + scene_data[4] + "_o" + scene_data[5] + "/"
        impth = os.path.join(directory, image_folder)

        # Run auto_tree_height_single
        if flag_proc == 0:
            ######## ROI_PAC results
            file_data = auto_tree_height_single_ROIPAC(impth, scene_data[2], scene_data[3], numLooks, noiselevel, flag_grad)
        elif flag_proc == 1:
            ######## ISCE results
            file_data = auto_tree_height_single_ISCE(impth, scene_data[2], scene_data[3], numLooks, noiselevel, flag_grad)
        else:
            print ("Invalid processor provided!!!")

        linkfile = os.path.join(directory, image_folder, filename + '_orig.mat')
        sio.savemat(linkfile,{'corr_vs':file_data[0],'kz':file_data[1],'coords':file_data[2]})

        # Write geodata to a text file (4th - 9th values in file_data) -> this gets stored in the scene specific folder
        geofile = open(os.path.join(directory, image_folder, filename + "_geo.txt"), "w")
        geofile.write("width: %d \n" % file_data[3])
        geofile.write("nlines: %d \n" % file_data[4])
        geofile.write("corner_lat: %f \n" % file_data[5])
        geofile.write("corner_lon: %f \n" % file_data[6])
        geofile.write("post_lat: %f \n" % file_data[7])
        geofile.write("post_lon: %f \n" % file_data[8])
        geofile.close()

    print ("auto_tree_height_many finished at " + (time.strftime("%H:%M:%S")))


# In[ ]:


def auto_tree_height_single_ISCE(directory, date1, date2, numLooks, noiselevel, flag_grad):
    # Extract ISCE parameters
    xmlfile = subprocess.getoutput('find ' + directory + '/int_' + date1 + '_' + date2 + '/ -name *Proc.xml')
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    root_tag = root.tag
    
    range_pixel_res = float(root.findall("./master/instrument/range_pixel_size")[0].text)
    llambda = float(root.findall("./master/instrument/radar_wavelength")[0].text)
    try:
        first_range = float(root.findall("./runTopo/inputs/range_first_sample")[0].text)
    except:
        first_range = float(root.findall("./runTopo/inputs/RANGE_FIRST_SAMPLE")[0].text)
    try:
        num_range_bin = int(root.findall("./runTopo/inputs/width")[0].text)
    except:
        num_range_bin = int(root.findall("./runTopo/inputs/WIDTH")[0].text)
    try:
        num_range_looks = int(root.findall("./runTopo/inputs/number_range_looks")[0].text)
    except:
        num_range_looks = int(root.findall("./runTopo/inputs/NUMBER_RANGE_LOOKS")[0].text)
    center_range = first_range + (num_range_bin/2-1)*range_pixel_res*num_range_looks
    incid_angle = float(root.findall("./master/instrument/incidence_angle")[0].text)
    baseline_top = float(root.findall("./baseline/perp_baseline_top")[0].text)
    baseline_bottom = float(root.findall("./baseline/perp_baseline_bottom")[0].text)
    baseline = (baseline_bottom + baseline_top)/2
    

    xmlfile = os.path.join(directory, "int_" + date1 + "_" + date2 + "/topophase.cor.geo.xml")
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    delta_array = np.array([])
    start_array = np.array([])
    size_array = np.array([], dtype=np.int32)
    for size in root.iter('property'):
        if size.items()[0][1] == 'size':
            size_array = np.append(size_array, int(size.find('value').text))
    for delta_val in root.iter('property'):
        if delta_val.items()[0][1] == 'delta':
            delta_array = np.append(delta_array, float(delta_val.find('value').text))
    for start_val in root.iter('property'):
        if start_val.items()[0][1] == 'startingvalue':
            start_array = np.append(start_array, float(start_val.find('value').text))
    end_array = start_array + size_array * delta_array
    north = max(start_array[1],end_array[1])
    south = min(start_array[1],end_array[1])
    east = max(start_array[0],end_array[0])
    west = min(start_array[0],end_array[0])
    coords = [north, south, west, east]
    geo_width = size_array[0]
    geo_nlines = size_array[1]
    corner_lat = north
    corner_lon = west
    step_lat = delta_array[1]
    step_lon = delta_array[0]

    xmlfile = os.path.join(directory, "int_" + date1 + "_" + date2 + "/resampOnlyImage.amp.geo.xml")
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    delta_array = np.array([])
    start_array = np.array([])
    size_array = np.array([], dtype=np.int32)
    for size in root.iter('property'):
        if size.items()[0][1] == 'size':
            size_array = np.append(size_array, int(size.find('value').text))
    if (size_array[0]<geo_width)|(size_array[1]<geo_nlines):
        for delta_val in root.iter('property'):
            if delta_val.items()[0][1] == 'delta':
                delta_array = np.append(delta_array, float(delta_val.find('value').text))
        for start_val in root.iter('property'):
            if start_val.items()[0][1] == 'startingvalue':
                start_array = np.append(start_array, float(start_val.find('value').text))
        end_array = start_array + size_array * delta_array
        north = max(start_array[1],end_array[1])
        south = min(start_array[1],end_array[1])
        east = max(start_array[0],end_array[0])
        west = min(start_array[0],end_array[0])
        coords = [north, south, west, east]
        geo_width = size_array[0]
        geo_nlines = size_array[1]
        corner_lat = north
        corner_lon = west
        step_lat = delta_array[1]
        step_lon = delta_array[0]


    # Read geolocated amp and cor files

    fid_cor = open(os.path.join(directory, "int_" + date1 + "_" + date2 + "/topophase.cor.geo"), "rb")
    cor_file = np.fromfile(fid_cor, dtype=np.dtype('<f'))
    pdb.set_trace()
    corr = cor_file.reshape(2 * geo_width, -1, order='F')
    corr = corr[:, 0:geo_nlines]
    corr_mag = corr[geo_width:2 * geo_width, :]

    fid_amp = open(os.path.join(directory, "int_" + date1 + "_" + date2 + "/resampOnlyImage.amp.geo"), "rb")
    amp_file = np.fromfile(fid_amp, dtype=np.dtype('<f'))
    inty = amp_file.reshape(2 * geo_width, -1, order='F')
    inty = inty[:,0:geo_nlines]
    inty1 = inty[::2,:]
    inty2 = inty[1::2,:]


    # Operations
    inty1 = np.power(inty1, 2) # Hardcoded based on 2 range looks and 10 azimuth looks
    inty2 = np.power(inty2, 2)

    inty1[inty1 <= 0] = np.NaN
    inty2[inty2 <= 0] = np.NaN
    corr_mag[corr_mag <= 0] = np.NaN

    ####### Noise level for ISCE-processed SAR backscatter power output
    if noiselevel == 0.0:
        if root_tag[0] == 'i':
            ####### ALOS thermal noise level (insarApp)
            N1 = 55.5**2
            N2 = 55.5**2
        elif root_tag[0] == 's':
            ####### ALOS thermal noise level (stripmapApp)
            N1 = (55.5/81)**2
            N2 = (55.5/81)**2
        else:
            raise Exception("invalid *Proc.xml file!!!")
    else:
        N1 = noiselevel
        N2 = noiselevel


    S1 = inty1 - N1
    g_th_1 = np.zeros(S1.shape)
    g_th_1[S1>N1] = np.sqrt(S1[S1>N1] / (S1[S1>N1] + N1))
    g_th_1[np.isnan(S1)] = np.NaN
    g_th_1[S1 <= N1] = np.NaN

    S2 = inty2-N2
    g_th_2 = np.zeros(S2.shape)
    g_th_2[S2>N2] = np.sqrt(S2[S2>N2] / (S2[S2>N2] + N2))
    g_th_2[np.isnan(S2)] = np.NaN
    g_th_2[S2 <= N2] = np.NaN

    g_th = g_th_1 * g_th_2

    corr_mag[corr_mag<0] = 0
    corr_mag[corr_mag>1] = 1
    corr_mag = remove_corr_bias(corr_mag,numLooks)
    corr_mag[corr_mag<0] = 0

    corr_vs = corr_mag / g_th

    # set constants
    pi = math.pi

    # correcting geometric decorrelation related to value compensation of ROI result compared to GAMMA. Caused by baseline/other decorrelation
    gamma_base = 1 - (2 * math.fabs(baseline) * math.cos(incid_angle / 180 * pi) * range_pixel_res / math.sin(incid_angle / 180 * pi) / llambda / center_range)
    gamma_geo = gamma_base
    corr_vs = corr_vs / gamma_geo
    corr_vs[corr_vs>1] = 1

    ##### Simple Radiometric correction of the coherences
    if flag_grad == 1:
        y = np.linspace(1, geo_width, geo_width)
        x = np.linspace(1, geo_nlines, geo_nlines)
        [X, Y] = np.meshgrid(x, y)
        A = np.vstack([X[~np.isnan(corr_vs)], Y[~np.isnan(corr_vs)], np.ones(np.size(corr_vs[~np.isnan(corr_vs)]))]).T
        coeff = np.linalg.lstsq(A, corr_vs[~np.isnan(corr_vs)])[0]
        corr_vs = corr_vs - X*coeff[0] - Y*coeff[1]
        corr_vs[corr_vs>1] = 1
        corr_vs[corr_vs<0] = 0

    kz = -2 * pi * 2 / llambda / center_range / math.sin(incid_angle/180*pi) * baseline
    kz = gdal.GVOT_MIN_TARGET_HEIGHT_FROM_DEM.fabs(kz)

    # Return corr_vs, kz, coords
    return corr_vs, kz, coords, geo_width, geo_nlines, corner_lat, corner_lon, step_lat, step_lon

# In[ ]:

def auto_tree_height_single_ROIPAC(directory, date1, date2, numLooks, noiselevel, flag_grad):

    # Extract ROI_PAC parameters
    amp_rsc_file = "int_"+date1+"_"+date2+"/"+date1+"-"+date2+".amp.rsc"
    range_pixel_res = read_rsc_data(amp_rsc_file, directory, "RANGE_PIXEL_SIZE")
    azimuth_pixel_res = read_rsc_data(amp_rsc_file, directory, "AZIMUTH_PIXEL_SIZE")


    geo_cor_rsc_file = "int_"+date1+"_"+date2+"/"+"geo_"+date1+"-"+date2+"_2rlks.cor.rsc"
    geo_width = int(read_rsc_data(geo_cor_rsc_file, directory, "WIDTH"))
    geo_nlines = int(read_rsc_data(geo_cor_rsc_file, directory, "FILE_LENGTH"))
    corner_lat = read_rsc_data(geo_cor_rsc_file, directory, "Y_FIRST")
    corner_lon = read_rsc_data(geo_cor_rsc_file, directory, "X_FIRST")
    step_lat = read_rsc_data(geo_cor_rsc_file, directory, "Y_STEP")
    step_lon = read_rsc_data(geo_cor_rsc_file, directory, "X_STEP")
    llambda = read_rsc_data(geo_cor_rsc_file, directory, "WAVELENGTH")

    int_rsc_file = "int_"+date1+"_"+date2+"/"+date1+"-"+date2+"-sim_SIM_2rlks.int.rsc"
    range1 = read_rsc_data(int_rsc_file, directory, "RGE_REF1")
    range2 = read_rsc_data(int_rsc_file, directory, "RGE_REF2")
    center_range = (range1 + range2) / 2 * 1000

    amp_4rlks_file = "int_"+date1+"_"+date2+"/"+date1+"-"+date2+"_2rlks.amp.rsc"
    incid_angle = read_rsc_data(amp_4rlks_file, directory, "BEAM")

    baseline_file = "int_"+date1+"_"+date2+"/"+date1+"_"+date2+"_baseline.rsc"
    p_baseline_1 = read_rsc_data(baseline_file, directory, "P_BASELINE_BOTTOM_HDR")
    p_baseline_2 = read_rsc_data(baseline_file, directory, "P_BASELINE_TOP_HDR")
    baseline = (p_baseline_1 + p_baseline_2) / 2

    # Read geolocated amp and cor files
    fid_cor = open(os.path.join(directory, "int_"+date1+"_"+date2+"/geo_"+date1+"-"+date2+"_2rlks.cor"), "rb")
    cor_file = np.fromfile(fid_cor, dtype=np.float32)
    corr_mag = cor_file.reshape(2*geo_width, geo_nlines, order='F')
    corr_mag = corr_mag[geo_width:len(corr_mag),:]

    fid_amp = open(os.path.join(directory, "int_"+date1+"_"+date2+"/geo_"+date1+"-"+date2+"_2rlks.amp"), "rb")
    amp_file = np.fromfile(fid_amp, dtype=np.float32)
    inty = amp_file.reshape(2*geo_width, geo_nlines, order='F')
    inty1 = inty[::2,:]
    inty2 = inty[1::2,:]

    # Set coordinate list
    coords = [corner_lat, corner_lat+(geo_nlines-1)*step_lat, corner_lon, corner_lon+(geo_width-1)*step_lon]
    
    # Operations
    inty1 = np.power(inty1,2) / 20                  # Hardcoded based on 2 range looks and 10 azimuth looks
    inty2 = np.power(inty2,2) / 20

    inty1[inty1 <= 0] = np.NaN
    inty2[inty2 <= 0] = np.NaN
    corr_mag[corr_mag <= 0] = np.NaN  

    ################### Noise level for ROI_PAC-processed SAR backscatter power output
    if noiselevel == 0.0:
        ####### ALOS thermal noise level
        N1 = 0.0192                                     
        N2 = 0.0192
    else:
        N1 = noiselevel
        N2 = noiselevel

    S1 = inty1 - N1
    g_th_1 = np.zeros(S1.shape)
    g_th_1[S1>0] = np.sqrt(S1[S1>0] / (S1[S1>0] + N1))
    g_th_1[np.isnan(S1)] = np.NaN
    g_th_1[S1 <= 0] = np.NaN
    
    S2 = inty2-N2
    g_th_2 = np.zeros(S2.shape)
    g_th_2[S2>0] = np.sqrt(S2[S2>0] / (S2[S2>0] + N2))
    g_th_2[np.isnan(S2)] = np.NaN
    g_th_2[S2 <= 0] = np.NaN
    
    g_th = g_th_1 * g_th_2

    corr_mag[corr_mag<0] = 0
    corr_mag[corr_mag>1] = 1
    corr_mag = remove_corr_bias(corr_mag,numLooks)
    corr_mag[corr_mag<0] = 0
    
    corr_vs = corr_mag / g_th
    
    # set constants
    pi=math.pi
    
    # correcting geometric decorrelation related to value compensation of ROI result compared to GAMMA. Caused by baseline/other decorrelation
    gamma_base = 1 - (2 * math.fabs(baseline) * math.cos(incid_angle / 180 * pi) * range_pixel_res / math.sin(incid_angle / 180 * pi) / llambda / center_range)
    gamma_geo = gamma_base
    corr_vs = corr_vs / gamma_geo                         
    corr_vs[corr_vs>1] = 1

    #################### Simple Radiometric correction of the coherences
    if flag_grad == 1:
        y = np.linspace(1, geo_width, geo_width)
        x = np.linspace(1, geo_nlines, geo_nlines)
        [X, Y] = np.meshgrid(x, y)
        A = np.vstack([X[~np.isnan(corr_vs)], Y[~np.isnan(corr_vs)], np.ones(np.size(corr_vs[~np.isnan(corr_vs)]))]).T
        coeff = np.linalg.lstsq(A, corr_vs[~np.isnan(corr_vs)])[0]
        corr_vs = corr_vs - X*coeff[0] - Y*coeff[1]
        corr_vs[corr_vs>1] = 1
        corr_vs[corr_vs<0] = 0

    kz = -2 * pi * 2 / llambda / center_range / math.sin(incid_angle/180*pi) * baseline
    kz = math.fabs(kz)

    # Return corr_vs, kz, coords
    return corr_vs, kz, coords, geo_width, geo_nlines, corner_lat, corner_lon, step_lat, step_lon
 
# In[ ]:

def read_rsc_data(filename, directory, param):
    # set default output value
    result = 1

    # set filename for file to be searched
    rsc_file = os.path.join(directory, filename)

    # Read parameters from file
    for line in open(rsc_file):
        if line.startswith(param):
            result = float(line.strip().split()[1])
    
    return result

# In[ ]:


def cal_KB(dp, edges, start_scene, link, directory, Nd_pairwise, Nd_self, bin_size, flag_sparse):

    # make output matrix of zeros
    YY = np.zeros((edges + 1) * 2)
    if link.size != 0:
        # for each edge run cal_KB_pairwise_new and put the output into YY
        for i in range(edges):
            print ("cal_KB_pwn started at " + (time.strftime("%H:%M:%S")))
            k_temp, b_temp = cal_KB_pairwise_new(int(link[i, 0]), int(link[i, 1]), dp[int((2*link[i, 0])-2)], dp[int((2*link[i, 0])-1)], dp[int((2*link[i, 1])-2)], dp[int((2*link[i, 1])-1)], directory, Nd_pairwise, bin_size)
            print ("cal_KB_pwn finished at " + (time.strftime("%H:%M:%S")))
            YY[2 * i] = k_temp
            YY[(2 * i) + 1] = b_temp

    # run cal_KB_self_new and put output into YY

    k_temp, b_temp = cal_KB_self_new(dp[int((2 * start_scene) - 2)], dp[int((2 * start_scene) - 1)], directory, Nd_self, bin_size, flag_sparse)
    YY[(2 * (edges + 1)) - 2] = k_temp
    YY[(2 * (edges + 1)) - 1] = b_temp
    
    # return Y
    return YY

# In[ ]:


def cal_KB_pairwise_new(scene1, scene2, deltaS1, deltaC1, deltaS2, deltaC2, directory, Nd_pairwise, bin_size):

    # Set main file name string as scene1_scene2
    file_str = str(scene1) + '_' + str(scene2)
   
    # Load and read data from .mat file
    # Samples and lines are calculated from the shape of the images
    selffile_data = sio.loadmat(os.path.join(directory, "output/", file_str + ".mat"))
    image1 = selffile_data['I1']
    image2 = selffile_data['I2']
    lines = int(image1.shape[0])
    samples = int(image1.shape[1])
    
    # S and C parameters are the average S and C plus the delta value
    S_param1 = 0.65 + deltaS1
    S_param2 = 0.65 + deltaS2
    C_param1 = 13 + deltaC1
    C_param2 = 13 + deltaC2
    
    # Create gamma and run arc_since for image1
    gamma1 = image1.copy()
    gamma1 = gamma1 / S_param1
    image1 = arc_sinc(gamma1, C_param1)
    image1[np.isnan(gamma1)] = np.nan

    # Create gamma and run arc_since for image2     
    gamma2 = image2.copy()
    gamma2 = gamma2 / S_param2
    image2 = arc_sinc(gamma2, C_param2)
    image2[np.isnan(gamma2)] = np.nan
    
    
    # Partition image into subsections for noise suppression (multi-step process)
    # Create M and N which are the number of subsections in each direction; fix() rounds towards zero
    # NX and NY are the subsection dimensions
    NX = Nd_pairwise
    NY = Nd_pairwise
    M = int(np.fix(lines / NY))
    N = int(np.fix(samples / NX))

    # Create JM and JN, which is the remainder after dividing into subsections
    JM = lines % NY
    JN = samples % NX

    # Select the portions of images that are within the subsections
    image1 = image1[0:lines - JM][:, 0:samples - JN]
    image2 = image2[0:lines - JM][:, 0:samples - JN]

    # Split each image into subsections and run mean_wo_nan on each subsection
    
    # Declare new arrays to hold the subsection averages
    image1_means = np.zeros((M, N))
    image2_means = np.zeros((M, N))
    
    # Processing image1
    # Split image into sections with NY number of rows in each
    image1_rows = np.split(image1, M, 0)
    for i in range(M):
        # split each section into subsections with NX number of columns in each
        row_array = np.split(image1_rows[i], N, 1)
        # for each subsection shape take the mean without NaN and save the value in another array
        for j in range(N):
            image1_means[i, j] = mean_wo_nan(row_array[j]) 

    # Processing image2
    # Split image into sections with NY number of rows in each
    image2_rows = np.split(image2, M, 0)
    for i in range(M):
        # split each section into subsections with NX number of columns in each
        row_array = np.split(image2_rows[i], N, 1)
         # for each subsection shape take the mean without NaN and save the value in another array
        for j in range(N):
            image2_means[i, j] = mean_wo_nan(row_array[j])           

    
    # Make an array for each image of where mean > 0 for both images
    IND1 = np.logical_and((image1_means > 0), (image2_means > 0))
    I1m_trunc = image1_means[IND1, ...]
    I2m_trunc = image2_means[IND1, ...]
    
    # Remove the overestimation at low height end (usually subjet to imperfection of the mask 
    # over water bodies, farmlands and human activities) and the saturation points over the forested areas due to logging

    I1m_trunc, I2m_trunc = remove_outlier(I1m_trunc, I2m_trunc, 0.5, 2)
    
    I1m_den = I1m_trunc   
    I2m_den = I2m_trunc

    # Calculate the covariance matrix of the data with outliers removed
    cov_matrix = np.cov(I1m_den, I2m_den)
     
    
    # Calculate the eigenvalues
    dA, vA = np.linalg.eig(cov_matrix)

    # Calculate K and B
    # K is based on whichever value in dA is the largest
    if (dA[0] > dA[1]): # dA[0] is largest
        K = vA[1, 0] / vA[0, 0]
    else: # dA[1] is largest
        K = vA[1, 1] / vA[0, 1]
    B = 2 * np.mean(I1m_den - I2m_den) / np.mean(I1m_den + I2m_den)
    
    return K, B


# In[ ]:


def cal_KB_self_new(deltaS2, deltaC2, directory, Nd_self, bin_size, sparse_lidar_flag):
    
    # Load and read data from .mat file
    # Samples and lines are calculated from the shape of the images
    selffile_data = sio.loadmat(os.path.join(directory, "output/self.mat"))
    image1 = selffile_data['I1']
    image2 = selffile_data['I2']
    lines = int(image1.shape[0])
    samples = int(image1.shape[1])

    # Set the S and C parameters to the average S and C plus the delta values
    S_param2 = 0.65 + deltaS2
    C_param2 = 13 + deltaC2

    # Create gamma and run arc_since for image2   
    gamma2 = image2.copy()
    gamma2 = gamma2 / S_param2
    image2 = arc_sinc(gamma2, C_param2)
    image2[np.isnan(gamma2)] = np.nan

    # Partition image into subsections for noise suppression (multi-step process)
    # Create M and N which are the number of subsections in each direction
    # NX and NY are the subsection dimensions
    NX = Nd_self
    NY = Nd_self
    M = int(np.fix(lines / NY))
    N = int(np.fix(samples / NX))

    # Create JM and JN, which is the remainder after dividing into subsections
    JM = lines % NY
    JN = samples % NX

    # Select the portions of images that are within the subsections
    image1 = image1[0:lines - JM][:, 0:samples - JN]
    image2 = image2[0:lines - JM][:, 0:samples - JN]

    # Split each image into subsections and run mean_wo_nan on each subsection
    
    # Declare new arrays to hold the subsection averages
    image1_means = np.zeros((M, N))
    image2_means = np.zeros((M, N))
    
    # Processing image1
    # Split image into sections with NY number of rows in each
    image1_rows = np.split(image1, M, 0)
    for i in range(M):
        # split each section into subsections with NX number of columns in each
        row_array = np.split(image1_rows[i], N, 1)
        # for each subsection shape take the mean without NaN and save the value in another array
        for j in range(N):
            image1_means[i, j] = mean_wo_nan(row_array[j])

    # Processing image2
    # Split image into sections with NY number of rows in each
    image2_rows = np.split(image2, M, 0)
    for i in range(M):
        # split each section into subsections with NX number of columns in each
        row_array = np.split(image2_rows[i], N, 1)
         # for each subsection shape take the mean without NaN and save the value in another array
        for j in range(N):
            image2_means[i, j] = mean_wo_nan(row_array[j]) 
    
    # Make an array for each image of where mean > 0 for both images
    IND1 = np.logical_and((image1_means > 0), (image2_means > 0))
    I1m_trunc = image1_means[IND1, ...]
    I2m_trunc = image2_means[IND1, ...]
    
    # Remove the overestimation at low height end (usually subjet to imperfection of the mask 
    # over water bodies, farmlands and human activities) and the saturation points over the forested areas due to logging
    IND2 = np.logical_or((I1m_trunc < 5), (I2m_trunc > (math.pi * C_param2 - 1)))
    IND2 = np.logical_not(IND2)


    # Call remove_outlier on these cells when there are only a few of lidar samples that are sparsely distributed
    if sparse_lidar_flag == 1:
        I1m_trunc = I1m_trunc[IND2, ...]
        I2m_trunc = I2m_trunc[IND2, ...]
        # Extract density values from the 2D scatter plot
        I1m_den, I2m_den = extract_scatterplot_density(I1m_trunc, I2m_trunc, bin_size)
    else:
        I1m_trunc, I2m_trunc = remove_outlier(I1m_trunc, I2m_trunc, 0.5, 2)
        I1m_den = I1m_trunc
        I2m_den = I2m_trunc



    # Calculate the covariance matrix of the data with outliers removed
    cov_matrix = np.cov(I1m_den, I2m_den)
    
    # Calculate the eigenvalues
    dA, vA = np.linalg.eig(cov_matrix)

    # Calculate K and B
    # K is based on whichever value in dA is the largest
    if (dA[0] > dA[1]): # dA[0] is largest
        K = vA[1, 0] / vA[0, 0]
    else: # dA[1] is largest
        K = vA[1, 1] / vA[0, 1]
    B = 2 * np.mean(I1m_den - I2m_den) / np.mean(I1m_den + I2m_den)
    
    return K, B


# In[ ]:


def cal_error_metric(dp, edges, start_scene, link, directory, N_pairwise, N_self):

    # make output matrix of zeros
    YY = np.zeros((edges + 1) * 2)
    if link.size != 0:

        # for each edge run cal_error_metric_pairwise and put the output into YY
        for i in range(edges):
            R_temp, RMSE_temp = cal_error_metric_pairwise(int(link[i, 0]), int(link[i, 1]), dp[int((2*link[i, 0])-2)], dp[int((2*link[i, 0])-1)], dp[int((2*link[i, 1])-2)], dp[int((2*link[i, 1])-1)], directory, N_pairwise)
            YY[2 * i] = R_temp
            YY[(2 * i) + 1] = RMSE_temp
    
    # run cal_error_metric_self and put output into YY
    R_temp, RMSE_temp = cal_error_metric_self(dp[int((2 * start_scene) - 2)], dp[int((2 * start_scene) - 1)], directory, N_self)
    YY[(2 * (edges + 1)) - 2] = R_temp
    YY[(2 * (edges + 1)) - 1] = RMSE_temp
    
    # return Y
    return YY


# In[ ]:


def cal_error_metric_pairwise(scene1, scene2, deltaS1, deltaC1, deltaS2, deltaC2, directory, N_pairwise):

    # Set main file name string as scene1_scene2
    file_str = str(scene1) + '_' + str(scene2)

    # Load and read data from .mat file
    # Samples and lines are calculated from the shape of the images
    selffile_data = sio.loadmat(os.path.join(directory, "output/" + file_str + ".mat"))
    image1 = selffile_data['I1']
    image2 = selffile_data['I2']
    lines = int(image1.shape[0])
    samples = int(image1.shape[1])

    # S and C parameters are the average S and C plus the delta value
    S_param1 = 0.65 + deltaS1
    S_param2 = 0.65 + deltaS2
    C_param1 = 13 + deltaC1
    C_param2 = 13 + deltaC2
    
    # Create gamma and run arc_since for image1
    gamma1 = image1.copy()
    gamma1 = gamma1 / S_param1
    image1 = arc_sinc(gamma1, C_param1)
    image1[np.isnan(gamma1)] = np.nan

    # Create gamma and run arc_since for image2     
    gamma2 = image2.copy()
    gamma2 = gamma2 / S_param2
    image2 = arc_sinc(gamma2, C_param2)
    image2[np.isnan(gamma2)] = np.nan

    # Partition image into subsections for noise suppression (multi-step process)
    # Create M and N which are the number of subsections in each direction; fix() rounds towards zero
    # NX and NY are the subsection dimensions
    NX = N_pairwise
    NY = N_pairwise
    M = int(np.fix(lines / NY))
    N = int(np.fix(samples / NX))

    # Create JM and JN, which is the remainder after dividing into subsections
    JM = lines % NY
    JN = samples % NX

    # Select the portions of images that are within the subsections
    image1 = image1[0:lines - JM][:, 0:samples - JN]
    image2 = image2[0:lines - JM][:, 0:samples - JN]

    # Split each image into subsections and run mean_wo_nan on each subsection
    
    # Declare new arrays to hold the subsection averages
    image1_means = np.zeros((M, N))
    image2_means = np.zeros((M, N))
    
    # Processing image1
    # Split image into sections with NY number of rows in each
    image1_rows = np.split(image1, M, 0)
    for i in range(M):
        # split each section into subsections with NX number of columns in each
        row_array = np.split(image1_rows[i], N, 1)
        # for each subsection shape take the mean without NaN and save the value in another array
        for j in range(N):
            image1_means[i, j] = mean_wo_nan(row_array[j]) 

    # Processing image2
    # Split image into sections with NY number of rows in each
    image2_rows = np.split(image2, M, 0)
    for i in range(M):
        # split each section into subsections with NX number of columns in each
        row_array = np.split(image2_rows[i], N, 1)
         # for each subsection shape take the mean without NaN and save the value in another array
        for j in range(N):
            image2_means[i, j] = mean_wo_nan(row_array[j])           


    # Make an array for each image of where mean > 0 for both images
    IND1 = np.logical_and((image1_means > 0), (image2_means > 0))
    I1m_trunc = image1_means[IND1, ...]
    I2m_trunc = image2_means[IND1, ...]
    
    R = np.corrcoef(I1m_trunc,I2m_trunc)
    R = R[0,1]
    RMSE = math.sqrt(sum((I1m_trunc-I2m_trunc)**2)/I1m_trunc.size)
    
    #Export the pair of heights for future scatter plot
    filename = file_str + "_I1andI2.json"
    R_RSME_file = open(os.path.join(directory, "output/" + filename), 'w')
    json.dump([I1m_trunc.tolist(), I2m_trunc.tolist()], R_RSME_file)
    R_RSME_file.close()
    
    return R, RMSE


# In[ ]:


def cal_error_metric_self(deltaS2, deltaC2, directory, N_self):

    # Load and read data from .mat file
    # Samples and lines are calculated from the shape of the images
    selffile_data = sio.loadmat(os.path.join(directory, "output/self.mat"))
    image1 = selffile_data['I1']
    image2 = selffile_data['I2']
    lines = int(image1.shape[0])
    samples = int(image1.shape[1])

    # Set the S and C parameters to the average S and C plus the delta values
    S_param2 = 0.65 + deltaS2
    C_param2 = 13 + deltaC2

    # Create gamma and run arc_since for image2   
    gamma2 = image2.copy()
    gamma2 = gamma2 / S_param2
    image2 = a.arc_sinc(gamma2, C_param2)
    image2[np.isnan(gamma2)] = np.nan

    # Partition image into subsections for noise suppression (multi-step process)
    # Create M and N which are the number of subsections in each direction
    # NX and NY are the subsection dimensions
    NX = N_self
    NY = N_self
    M = int(np.fix(lines / NY))
    N = int(np.fix(samples / NX))

    # Create JM and JN, which is the remainder after dividing into subsections
    JM = lines % NY
    JN = samples % NX

    # Select the portions of images that are within the subsections
    image1 = image1[0:lines - JM][:, 0:samples - JN]
    image2 = image2[0:lines - JM][:, 0:samples - JN]

    # Split each image into subsections and run mean_wo_nan on each subsection
    
    # Declare new arrays to hold the subsection averages
    image1_means = np.zeros((M, N))
    image2_means = np.zeros((M, N))
    
    # Processing image1
    # Split image into sections with NY number of rows in each
    image1_rows = np.split(image1, M, 0)
    for i in range(M):
        # split each section into subsections with NX number of columns in each
        row_array = np.split(image1_rows[i], N, 1)
        # for each subsection shape take the mean without NaN and save the value in another array
        for j in range(N):
            image1_means[i, j] = mean_wo_nan(row_array[j])

    # Processing image2
    # Split image into sections with NY number of rows in each
    image2_rows = np.split(image2, M, 0)
    for i in range(M):
        # split each section into subsections with NX number of columns in each
        row_array = np.split(image2_rows[i], N, 1)
         # for each subsection shape take the mean without NaN and save the value in another array
        for j in range(N):
            image2_means[i, j] = mean_wo_nan(row_array[j]) 
    
    # Make an array for each image of where mean > 0 for both images
    IND1 = np.logical_and((image1_means > 0), (image2_means > 0))
    I1m_trunc = image1_means[IND1, ...]
    I2m_trunc = image2_means[IND1, ...]
    
    R = np.corrcoef(I1m_trunc,I2m_trunc)
    R = R[0,1]
    RMSE = math.sqrt(sum((I1m_trunc-I2m_trunc)**2)/I1m_trunc.size)
    
    #Export the pair of heights for future scatter plot
    filename = "self_I1andI2.json"
    R_RSME_file = open(os.path.join(directory, "output/" + filename), 'w')
    json.dump([I1m_trunc.tolist(), I2m_trunc.tolist()], R_RSME_file)
    R_RSME_file.close()
    
    return R, RMSE


# In[ ]:

def create_mosaic(directory,mosaicfile,listoffiles):

        print ("Create mosaic started at " + (time.strftime("%H:%M:%S")))

        subprocess.getoutput('gdalbuildvrt -separate -srcnodata 255 -overwrite '+directory+'mosaic.vrt '+listoffiles)
        subprocess.getoutput('gdal_translate -of GTiff -a_nodata 255 '+directory+'mosaic.vrt '+directory+'mosaic.tif')

        # Load mosaic.tif and associated parameters - .tif
        driver = gdal.GetDriverByName('GTiff')
        driver.Register()
        img = gdal.Open(os.path.join(directory,'mosaic.tif'))
        ref_data = np.array(img.ReadAsArray())
        refgeotrans = img.GetGeoTransform()
        corner_lon = refgeotrans[0]
        post_lon = refgeotrans[1]
        corner_lat = refgeotrans[3]
        post_lat = refgeotrans[5]
        geo_width = img.RasterXSize
        geo_lines = img.RasterYSize

        ######################## average all of the overlapping pixels at the same area
        ref_data = np.single(ref_data)
        ref_data[ref_data==255] = np.NaN
        avg = np.nanmean(ref_data,axis=0)
        avg[np.isnan(avg)] = 255

        ################## Create the final GeoTiff
        driver = gdal.GetDriverByName('GTiff')

        outRaster = driver.Create(directory+mosaicfile, geo_width, geo_lines)
        outRaster.SetGeoTransform([corner_lon, post_lon, 0, corner_lat, 0, post_lat])
        outband = outRaster.GetRasterBand(1)

        outband.WriteArray(avg)
        outRasterSRS = osr.SpatialReference()
        outRasterSRS.ImportFromEPSG(4326)
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()

        print ("Create mosaic finished at " + (time.strftime("%H:%M:%S")))


        print ("Final mosaic generation done !!!")



# In[ ]:


def extract_scatterplot_density(x, y, bin_size=100, threshold=0.5):

	values, xedges, yedges = np.histogram2d(x, y, bin_size)
	xbin_center = xedges[0:-1] + (xedges[1]-xedges[0])/2
	ybin_center = yedges[0:-1] + (yedges[1]-yedges[0])/2
	max_den = np.max(values)
	threshold_den = max_den * threshold
	[BCX, BCY] = np.meshgrid(xbin_center, ybin_center)
	values = values.transpose()
	IND_den = (values >= threshold_den)
	Hm_den = BCX[IND_den]
	Pm_den = BCY[IND_den]

	return Hm_den, Pm_den


# In[ ]:


def flag_scene_file(flagfilename, flag, directory):
    
    # Open the file
    flagfile = open(os.path.join(directory, flagfilename))
    
    # Set default value for scene_file
    data_array = ["", "", "", "", "", ""]
    
    # For each line in the file compare the line flag with the input flag
    for line in flagfile:
        
        # Set the line values
        line = line.strip().split()
        lineflag = line[0]    
        
        # Compare line and input flags
        if int(lineflag) == flag:
            data_array = list(line)
            
    # Close file
    flagfile.close()
    
    # Print error message is input flag not found
    if(data_array[0] == ""):
        print ("ERROR: Invalid flag number for the given text file")
        
    # Return scene_file
    return data_array


# In[ ]:


def forest_stand_height(scenes, edges, start_scene, iterations, linkfilename, flagfile, ref_file, maskfile, file_directory, filetypes=['gif', 'json', 'kml', 'mat', 'tif'], Nd_pairwise=20, Nd_self=20, N_pairwise=20, N_self=20, bin_size=100, flag_sparse=0, flag_diff=0, flag_error=0, numLooks=20, noiselevel=0.0, flag_proc=0, flag_grad=0):

    print ("Forest Stand Height started at " + (time.strftime("%H:%M:%S")))

    # Set error warnings to ignore "invalid value" warnings caused by NaN values
    np.seterr(invalid='ignore')

    if flag_sparse == 1:
        Nd_self = 1

    subprocess.getoutput('mkdir '+file_directory+'output')

    # Extract the correlation map, kz, and corner coordinates for each scene
    auto_tree_height_many(scenes, flagfile, file_directory, numLooks, noiselevel, flag_proc, flag_grad)

    if linkfilename == '-':
        # Run intermediate_self() (Central scene and LiDAR overlap)
        intermediate_self(start_scene, flagfile, ref_file, maskfile, file_directory)
        edge_array = np.array([])
        print ("Intermediate_self finished at" + (time.strftime("%H:%M:%S")))
    else:
        # Read in the list of edges
        edge_array = read_linkfile(edges, linkfilename, file_directory)
        # Calculate the overlap areas between the different scenes and the LiDAR (or other groundtruth)
        intermediate(edges, start_scene, edge_array, maskfile, flagfile, ref_file, file_directory)

    # Mosaic the interferograms
    auto_mosaicking_new(scenes, edges, start_scene, iterations, edge_array, file_directory, Nd_pairwise, Nd_self, bin_size, flag_sparse)

    # Store the delta S and C values for each scene
    write_deltaSC(scenes, iterations, flagfile, file_directory)

    # Create the tree height map
    write_mapfile_new(scenes, flagfile, maskfile, file_directory, filetypes)

    if flag_diff == 1:
        # Create the diff_height map
        write_diff_height_map(start_scene, ref_file, flagfile, maskfile, file_directory, filetypes)


    # Run cal_error_metric() when error metrics/scatter plots are needed
    if flag_error == 1:
        # Load the dp vector from the final iteration
        filename = "SC_%d_iter.json" % iterations
        ##        filename = "SC_%d_iter.json" % 6
        iter_file = open(os.path.join(file_directory, "output/" + filename))
        file_data = json.load(iter_file)
        iter_file.close()
        dp = np.array(file_data[0])
        # Run cal_error_metric() and create a json file containing all of the "pairwise" and "self" R & RMSE error measures
        Y = cal_error_metric(dp, edges, start_scene, edge_array, file_directory, N_pairwise, N_self)
        output_file = open(os.path.join(file_directory, "output/error_metric.json"), 'w')
        json.dump([Y.tolist()], output_file)
        output_file.close()
        print ("cal_error_metric file written at " + (time.strftime("%H:%M:%S"))) 


# In[ ]:

def intermediate(edges, start_scene, linkarray, maskfile, flagfile, ref_file, directory):
    
    # For each edge run intermediate_pairwise
    for i in range(edges):
        intermediate_pairwise(linkarray[i, 0], linkarray[i, 1], flagfile, maskfile, directory)
        print (("%d edge file(s) created at " % (i + 1)) + (time.strftime("%H:%M:%S")))
        
    # Run intermediate_self() (Central scene and LiDAR overlap)
    intermediate_self(start_scene, flagfile, ref_file, maskfile, directory)
    
    print ("intermediate() complete - overlap areas calculated at " + (time.strftime("%H:%M:%S")))

# In[ ]:


def intermediate_pairwise(flag1, flag2, flagfile, maskfile, directory):

    # Get flag-scene file data
    scene1_data = flag_scene_file(flagfile, flag1, directory)
    scene2_data = flag_scene_file(flagfile, flag2, directory)

    # Set file names based on flags
    filename1 = scene1_data[1]
    filename2 = scene2_data[1]
    
    # Set the image folder names
    image1_folder = "f" + scene1_data[4] + "_o" + scene1_data[5] + "/"
    image2_folder = "f" + scene2_data[4] + "_o" + scene2_data[5] + "/"

    file1 = sio.loadmat(os.path.join(directory, image1_folder, filename1 + "_orig.mat"))
    corr1 = file1['corr_vs']
    kz1 = file1['kz'][0][0]
    coords1 = file1['coords'][0]

    file2 = sio.loadmat(os.path.join(directory, image2_folder, filename2 + "_orig.mat"))
    corr2 = file2['corr_vs']
    kz2 = file2['kz'][0][0]
    coords2 = file2['coords'][0]
    
    # Set D constant --- D = 1 arc second
    D = 2.7777778 * (10**-4)

    # Remove non-forest from both images
    if maskfile != '-':
        corr1 = remove_nonforest(corr1, coords1, maskfile, directory)
        corr2 = remove_nonforest(corr2, coords2, maskfile, directory)
    
    # Set the image boundaries
    north1 = coords1[0]
    south1 = coords1[1]
    west1 = coords1[2]
    east1 = coords1[3]
    north2 = coords2[0]
    south2 = coords2[1]
    west2 = coords2[2]
    east2 = coords2[3]

    # Determine boundaries of the overlap area
    overlap_north = min(north1, north2)
    overlap_south = max(south1, south2)
    overlap_east = min(east1, east2)
    overlap_west = max(west1, west2)

    # Calculate overlap boundaries in coordinates of each image (ex image1[1000-1200] vs image2[0-200])
    xw1 = int(round(((overlap_west - west1) / D) + 1))
    xe1 = int(round(((overlap_east - west1) / D) + 1))
    xn1 = int(round((-(overlap_north - north1) / D) + 1))
    xs1 = int(round((-(overlap_south - north1) / D) + 1))
    xw2 = int(round(((overlap_west - west2) / D) + 1))
    xe2 = int(round(((overlap_east - west2) / D) + 1))
    xn2 = int(round((-(overlap_north - north2) / D) + 1))
    xs2 = int(round((-(overlap_south - north2) / D) + 1))

    # Set overlap sections from each image
    I1 = corr1[xw1-1:xe1][:, xn1-1:xs1]
    I2 = corr2[xw2-1:xe2][:, xn2-1:xs2]
    
    # Set average S and C parameters based on the average S and C (0<s<1, 0<c<20 so s=0.65 and c=13)
    S_param1 = 0.65
    C_param1 = 13
    S_param2 = 0.65
    C_param2 = 13

    # Create grid for image1
    [Dy1, Dx1] = I1.shape
    x1 = np.linspace(0, 1, Dx1)
    y1 = np.linspace(0, 1, Dy1)
    [X1, Y1] = np.meshgrid(x1, y1)
    
    # Create grid for image2
    [Dy2, Dx2] = I2.shape    
    x2 = np.linspace(0, 1, Dx2)
    y2 = np.linspace(0, 1, Dy2)
    [X2, Y2] = np.meshgrid(x2, y2)

    # Set NaN values to -100 to avoid interpolation errors
    I1[np.isnan(I1)] = -100
    I2[np.isnan(I2)] = -100
    
    # Co-register the two images
    I2 = interpolate.griddata((X2.flatten(), Y2.flatten()), I2.flatten(), (X1, Y1), method='nearest')

    # Reset NaN values
    IND1 = (I1 == -100)
    IND2 = (I2 == -100)
    IND = np.logical_or(IND1, IND2)
    I1[IND] = np.NaN
    I2[IND] = np.NaN

    
    # Save link file using MAT
    linkfilename = "%s_%s.mat" % (int(flag1), int(flag2))
    linkfile = os.path.join(directory, "output/" + linkfilename)
    sio.savemat(linkfile,{'I1':I1,'I2':I2})


# In[ ]:


def intermediate_self(start_scene, flagfile, ref_file, maskfile, directory):
   
    # Set scene data, file name, and image folder name
    scene2_data = flag_scene_file(flagfile, start_scene, directory)
    filename2 = scene2_data[1]
    image_folder = "f" + scene2_data[4] + "_o" + scene2_data[5] + "/"

    # Set D constant --- D = 1 arc second, this parameter is based on the use of ALOS data
    #D = 8.3333333 * (10**-4)
    D = 2.77777778 * (10**-4)  #for tracy test case
    
    # Load central image file and associated parameters

    file2 = sio.loadmat(os.path.join(directory, image_folder, filename2 + "_orig.mat"))
    corr2 = file2['corr_vs']
    kz2 = file2['kz'][0][0]
    coords2 = file2['coords'][0]

    # Remove non-forest from the image
    if maskfile != '-':
        corr2 = remove_nonforest(corr2, coords2, maskfile, directory)
    
    # Load LiDAR files and associated parameters - .tif
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()
    img = gdal.Open(os.path.join(directory, ref_file))
    ref_data = np.array(img.ReadAsArray())
    refgeotrans = img.GetGeoTransform()
    corner_lon = refgeotrans[0]
    post_lon = refgeotrans[1]
    corner_lat = refgeotrans[3]
    post_lat = refgeotrans[5]
    width = img.RasterXSize
    lines = img.RasterYSize
    
    # Set LiDAR parameters into correct format
    corr1 = ref_data.transpose()
    corr1[corr1 < 0] = np.NaN   # set margin areas to NaN
    coords1 = np.array([corner_lat, corner_lat + (lines * post_lat), corner_lon, corner_lon + (width * post_lon)])
    
    # Set the image boundaries
    north1 = coords1[0]
    south1 = coords1[1]
    west1 = coords1[2]
    east1 = coords1[3]
    north2 = coords2[0]
    south2 = coords2[1]
    west2 = coords2[2]
    east2 = coords2[3]

    # Determine boundaries of the overlap area
    overlap_north = min(north1, north2)
    overlap_south = max(south1, south2)
    overlap_east = min(east1, east2)
    overlap_west = max(west1, west2)

    # Calculate overlap boundaries in coordinates of image2 (ex image1[1000-1200] vs image2[0-200])
    xw2 = int(round(((overlap_west - west2) / D) + 1))
    xe2 = int(round(((overlap_east - west2) / D) + 1))
    xn2 = int(round((-(overlap_north - north2) / D) + 1))
    xs2 = int(round((-(overlap_south - north2) / D) + 1))
  
    # Set overlap sections for the LiDAR and SAR images
    I1 = corr1.copy()
    I2 = corr2[xw2-1:xe2][:, xn2-1:xs2]


    # Set average S and C parameters based on the average S and C (0<s<1, 0<c<20 so s=0.65 and c=13)
    S_param2 = 0.65
    C_param2 = 13

    # Create grid for image1
    [Dy1, Dx1] = I1.shape
    x1 = np.linspace(0, 1, Dx1)
    y1 = np.linspace(0, 1, Dy1)
    [X1, Y1] = np.meshgrid(x1, y1)
    
    # Create grid for image2
    [Dy2, Dx2] = I2.shape    
    x2 = np.linspace(0, 1, Dx2)
    y2 = np.linspace(0, 1, Dy2)
    [X2, Y2] = np.meshgrid(x2, y2)

    # Set NaN values to -100 to avoid interpolation errors
    I1[np.isnan(I1)] = -100
    I2[np.isnan(I2)] = -100
    

    # Co-register the two images
    I2 = interpolate.griddata((X2.flatten(), Y2.flatten()), I2.flatten(), (X1, Y1), method='nearest')

    # Reset NaN values
    IND1 = (I1 == -100)
    IND2 = (I2 == -100)
    IND = np.logical_or(IND1, IND2)
    I1[IND] = np.NaN
    I2[IND] = np.NaN
    
    # Save link file using JSON
    linkfilename = "self.mat"
    linkfile = os.path.join(directory, "output/" + linkfilename)
    sio.savemat(linkfile,{'I1':I1,'I2':I2})

# In[ ]:


def ls_deltaSC(dp, edges, scenes, start_scene, linkarray, directory, Nd_pairwise, Nd_self, bin_size, flag_sparse):

    # run cal_KB
    y = cal_KB(dp, edges, start_scene, linkarray, directory, Nd_pairwise, Nd_self, bin_size, flag_sparse)

    # Create a blank array for the Jacobi matrix
    jacobi = np.zeros(4 * scenes * (edges + 1)) #(2 * (edges + 1)) * (2 * scenes)
    jacobi = np.reshape(jacobi, (2 * (edges + 1), scenes * 2))
    


    # Fill in the Jacobi matrix
    for i in range(scenes):
        # fill K section
        temp = dp.copy()
        temp[2 * i] = temp[2 * i] + 0.1 # i and i+1 instead of i-1 and i due to index 0 vs 1
        temp = cal_KB(temp, edges, start_scene, linkarray, directory, Nd_pairwise, Nd_self, bin_size, flag_sparse)
        jacobi[:, 2 * i] = np.reshape(((temp - y) / 0.1), (2 * (edges + 1), )) # reshape temp-y part to 1D (a) instead of 2D (ax1)


        # fill B section       
        temp = dp.copy()
        temp[(2 * i) + 1] = temp[(2 * i) + 1] + 1
        temp = cal_KB(temp, edges, start_scene, linkarray, directory, Nd_pairwise, Nd_self, bin_size, flag_sparse)
        jacobi[:, (2 * i) + 1] = np.reshape(((temp - y) / 1), (2 * (edges + 1), ))
    
    # Create matrix of target K and B values (K = 1, B = 0 in order K-B-K-B-K-B...)
    target = np.zeros((edges + 1) * 2)
    target[::2] = 1

    # Calculate the change in S and C 
    changeSC = np.dot(np.dot(np.linalg.inv(np.dot(jacobi.conj().transpose(), jacobi)), jacobi.conj().transpose()), (target - y))
    
    changeSC = changeSC + dp
    YY = cal_KB(changeSC, edges, start_scene, linkarray, directory, Nd_pairwise, Nd_self, bin_size, flag_sparse)
    res = sum((YY - target)**2)

    # return changeSC and res
    return changeSC, res


# In[ ]:


def mean_wo_nan(A):

    # Copy and flatten A
    B = A.copy().flatten('F')
    
    # Remove NaN values from B
    B = B[~np.isnan(B)]
    
    # Return the mean of B
    return np.mean(B)


# In[ ]:

def read_geo_data(coord_file, directory):
    
    # Set filename for file to be searched
    filename = os.path.join(directory, coord_file)
    
    # Read parameters based on file type GeoTIFF or ROI_PAC text file)
    if(coord_file[-3:] == "tif"):
        # Read GeoTIFF
        driver = gdal.GetDriverByName('GTiff')
        driver.Register()
        image = gdal.Open(filename)
        refgeotrans = image.GetGeoTransform()
        corner_long = refgeotrans[0]
        post_long = refgeotrans[1]
        corner_lat = refgeotrans[3]
        post_lat = refgeotrans[5]
        width = image.RasterXSize
        nlines = image.RasterYSize
        
    else:
        # Read ROI_PAC text file
        for line in open(filename):
            if line.startswith("width"):
                width = np.int(line.strip().split()[1])
            elif line.startswith("nlines"):
                nlines = np.int(line.strip().split()[1])
            elif line.startswith("corner_lat"):
                corner_lat = np.float(line.strip().split()[1])
            elif line.startswith("corner_lon"):
                corner_long = np.float(line.strip().split()[1])
            elif line.startswith("post_lat"):
                post_lat = np.float(line.strip().split()[1])
            elif line.startswith("post_lon"):
                post_long = np.float(line.strip().split()[1])
            
    return width, nlines, corner_lat, corner_long, post_lat, post_long


# In[ ]:


def read_linkfile(edges, filename, directory):
    if edges > 0:
        # Open the file
        linkfile = open(os.path.join(directory, filename))
    
        # Create output array
        linkarray = np.zeros(edges * 2).reshape(edges, 2)
    
        # Set line counter
        counter = 0
    
        # For each line in the file compare the line flag with the input flag
        for line in linkfile:
        
            # Set array values from each line
            line = line.strip().split()
            linkarray[counter][0] = line[0]
            linkarray[counter][1] = line[1]
        
            # Increment counter
            counter += 1
            
        # Close file
        linkfile.close()
        
        # Return linkarray
        return linkarray
    else:
        linkarray = []


# In[ ]:


def remove_corr_bias(C,L):
    
    # Set m and D arrays that correspond to ROI_PAC. In Matlab these values can be created using the hypergeom function.

    k = 1.0
    D = np.double(np.arange(0,1,0.01))
    m= np.array([])
    for i in range(D.shape[0]):
        m= np.append(m,mp.gamma(L)*mp.gamma(1+k/2)/mp.gamma(L+k/2)*mp.hyper([1+k/2,L,L],[L+k/2,1],D[i]**2)*(1-D[i]**2)**L)
    m=np.double(m)
    p=np.arange(0,min(m),0.01)
    p=np.double(p)
    m=np.append(np.append(p, m), np.array([1]))
    D=np.append(np.append(0*p, D), np.array([1]))

    # Run interpolation
    set_interp = interpolate.interp1d(m, D, kind='cubic')
    YC = set_interp(C)

    return YC


# In[ ]:


def remove_nonforest(I, func_coords, maskfile, directory):
    
    # Load mask file as a GeoTIFF
    maskfile = gdal.Open(os.path.join(directory, maskfile))
    mask = np.array(maskfile.ReadAsArray())
    
    # Set any NaN values to 1 (aka not a forest)
    mask[np.isnan(mask)] = 1
    
    # Get mask geo parameters
    width = maskfile.RasterXSize
    nlines = maskfile.RasterYSize
    maskgeotrans = maskfile.GetGeoTransform()
    corner_lon = maskgeotrans[0]
    post_lon = maskgeotrans[1]
    corner_lat = maskgeotrans[3]
    post_lat = maskgeotrans[5]
    
    # Transpose mask so that it matches orientation of the radar data
    mask = mask.transpose()
    widthT = nlines
    nlinesT = width
    
    # Set coordinates based on file parameters
    file_coords = np.array([corner_lat, (corner_lat + (nlinesT - 1.0) * post_lat), corner_lon, (corner_lon + (widthT - 1.0) * post_lon)])
    
    # Calculate overlap boundaries in new coordinate system
    xw = int(round(((func_coords[2] - file_coords[2]) / post_lon) + 1))
    xe = int(round(((func_coords[3] - file_coords[2]) / post_lon) + 1))
    xn = int(round(((func_coords[0] - file_coords[0]) / post_lat) + 1))
    xs = int(round(((func_coords[1] - file_coords[0]) / post_lat) + 1))

    # Trim mask
    mask = np.logical_not(mask[xw-1:xe][:, xn-1:xs])
    
    # Get size of image and mask
    [m, n] = I.shape
    [M, N] = mask.shape
    
    # Make range of values from 0-1 based on M and N (not including 1), and run linspace
    x = np.linspace(0, 1, N, endpoint=False)
    y = np.linspace(0, 1, M, endpoint=False)
    [X, Y] = np.meshgrid(x, y)

    # Make range of values from 0-1 based on m and n (not including 1), and run linspace
    xp = np.linspace(0, 1, n, endpoint=False)
    yp = np.linspace(0, 1, m, endpoint=False) 
    [XP, YP] = np.meshgrid(xp, yp)    
    
    # Run interpolation
    O = interpolate.griddata((X.flatten(), Y.flatten()), mask.flatten(), (XP, YP), method='nearest')
    O = np.double(O)
    O[O == 0] = np.NaN
    O = I * O
    return O


# In[ ]:


def remove_outlier(x, y, win_size=0.5, threshold=5):

    # initialize other variables
    outliers_ind = [] # in Matlab code this variable is IND
    ind_x = np.zeros(x.size)
    ind_y = np.zeros(x.size)
    ind = np.zeros(x.size)

    # For each value in x check more or less neighboring points within the given window than the given threshold
    for i in range(x.size):

        # set base equal to a pair of x, y values        
        current_x = x[i]
        current_y = y[i]

        # for each x and y check if they are within +- the window from the current x(i) and y(i)
        # store a list of where both a and y are within the window in the array ind
        ind_x = (x > current_x - win_size) & (x < current_x + win_size) 
        ind_y = (y > current_y - win_size) & (y < current_y + win_size)
        ind = ind_x & ind_y

        # if for the current i there are fewer nearby points (within window) than the threshold
        if sum(ind) <= threshold:
            # then append it to the outliers array
            outliers_ind.append(i)

    # make new copies of x and y and delete all of the outlying points
    XX = np.delete(x, outliers_ind)
    YY = np.delete(y, outliers_ind)

    # return XX and YY (aka x and y without the outliers)
    return XX, YY


# In[ ]:


def write_deltaSC(scenes, N, flagfile, directory):

    # Load dp data from final iteration .json file
    filename = "SC_%d_iter.json" % N

    selffile = open(os.path.join(directory, "output/" + filename))
    selffile_data = json.load(selffile)
    dp = np.array(selffile_data[0])

    # For each scene name the file, create delta S and C, and save them to the file
    for i in range(scenes):

        # Set file name
        scene_data = flag_scene_file(flagfile, i + 1, directory) # 0 vs 1 indexing
        filename = scene_data[1]
        image_folder = "f" + scene_data[4] + "_o" + scene_data[5] + "/"

        # Calculate delta S and C
        DS = dp[2 * i]
        DC = dp[(2 * i) + 1]

        # Save DS and DC to output .json file
        outfile = open(os.path.join(directory, image_folder, filename + '_tempD.json'), "w")
        json.dump([DS, DC], outfile)
        outfile.close()

    print ("write_deltaSC completed at " + (time.strftime("%H:%M:%S")))



# In[ ]:


def write_diff_height_map(start_scene, ref_file, flagfile, maskfile, directory, output_files):
    
    if isinstance(start_scene,int) == 1:

        # Load and read data from .mat file
        # Samples and lines are calculated from the shape of the images
        file_data = sio.loadmat(os.path.join(directory, "output/self.mat"))
        lidar = file_data['I1']
        corr_vs = file_data['I2']
        
    
        # Load and read data from temp .json files
        scene_data = flag_scene_file(flagfile, start_scene, directory) # 0 vs 1 indexing
        filename = scene_data[1]
        image_folder = "f" + scene_data[4] + "_o" + scene_data[5] + "/"
        file_tempD = open(os.path.join(directory, image_folder, filename + "_tempD.json"))
        B = json.load(file_tempD)
        
        # Set S and C paramters based on the default and data from B
        S_param = 0.65 + B[0]
        C_param = 13 + B[1]
        
        # Run sinc model to calculate the heights
        gamma = corr_vs.copy()
        gamma = gamma / S_param
        height = arc_sinc(gamma, C_param)
        
        # Calculate the diff_height map (diviation of the InSAR inverted height away from the lidar height)
        diff_height = lidar - height
        
        # Transpose height to correctly align it (ie so it isn't rotated in relation to an underlying map)
        diff_height = diff_height.transpose()
        
        # Get rid of NaN so future processing software doesn't error

        diff_height[np.isnan(diff_height)] = 255
    
        for filetype in output_files:
            write_file_type(diff_height, "diff_height", filename, os.path.join(directory, image_folder), filetype, 0, ref_file)
            
    print ("all diff_height output files written at " + (time.strftime("%H:%M:%S")))
    
        


# In[ ]:


def write_file_type(data, outtype, filename, directory, filetype, coords, ref_file=""):


    # Use if/else to determine the desired type of file output
    if(filename[-8:] == "_255_255"):
        outfilename = filename[:-4]
    elif((filename[-4:] == "_fsh") or (filename[-5:] == "_diff") or (filename[-4:] == "_255")):
        outfilename = filename
    else:
        if(outtype == "stand_height"):
            outfilename = filename + "_fsh"
        elif(outtype == "diff_height"):
            outfilename = filename + "_diff"

    # Use if/else to determine the desired type of file output
    # Create .gif output
    if(filetype == "gif"):
        # Check if a 0-255 .tif with the same filename already exists, and if not create it.
        if (os.path.isfile(os.path.join(directory,outfilename + "_255.tif")) == True):
            gif_img = Image.open(os.path.join(directory, outfilename + "_255.tif"))
        else:
            # Set array in a 0-255 range for gif/kml
            # Get dimensions of array and then flatten for use with nonzero()
            (row, col) = data.shape
            data = data.flatten()
            # Get the nonzero indices and min/max
            nz_IND = np.nonzero(data)
            nz_min = data[nz_IND[0]].min()
            nz_max = data[nz_IND[0]].max()
            # Set the scaled values
            data255 = data.copy()
            data255[nz_IND[0]] = (data[nz_IND[0]] - nz_min) * (255 / (nz_max - nz_min)) + 1
            # Reshape the array of scaled values
            data255 = np.reshape(data255, (row, col))
            data = np.reshape(data, (row, col))

            # Write 0-255 .tif
            write_file_type(data255, outtype, outfilename + "_255", directory, "tif", coords, ref_file)
            gif_img = Image.open(os.path.join(directory, outfilename + "_255.tif"))

        # Create the .gif
        gif_img.save(os.path.join(directory, outfilename + "_255.gif"), "GIF", transparency=0)

    # Create .json output
    elif(filetype == "json"):
        jsonfile = open(os.path.join(directory, outfilename + '.json'), 'w')
        json.dump([data.tolist()], jsonfile)
        jsonfile.close()

    # Create .kml output
    elif(filetype == "kml"):
        # Determine the realname based on whether or not a single image is being processed or a pair
        if(filename[3] == "_"):   # pair
            realname = filename[:31]
        else:
            realname = filename[:23]

        # Read geo location information in from a text or geotiff file depending on outtype
        if(outtype == "stand_height"):
            (width, lines, north, west, lat_step, long_step) = read_geo_data(realname + "_geo.txt", directory)
            north = coords[0]
            west = coords[2]
            south = coords[1]
            east = coords[3]
        elif(outtype == "diff_height"):
            (width, lines, north, west, lat_step, long_step) = read_geo_data(ref_file, directory[:-10])
            south = north + (lat_step * lines)
            east = west + (long_step * width)


        # Check if a .gif with the same filename does not already exist then create it.
        if (os.path.isfile(os.path.join(directory, outfilename + "_255.gif")) == False):
            write_file_type(data, outtype, outfilename, directory, "gif", coords, ref_file)

        # Create the .kml
        kml = simplekml.Kml()
        arraykml = kml.newgroundoverlay(name=outfilename)
        arraykml.icon.href = os.path.join(directory, outfilename + "_255.gif")
        arraykml.latlonbox.north = north
        arraykml.latlonbox.south = south
        arraykml.latlonbox.east = east
        arraykml.latlonbox.west = west
        kml.save(os.path.join(directory, outfilename + "_255.kml"))

    # Create .mat output
    elif(filetype == "mat"):
        sio.savemat(os.path.join(directory, outfilename + '.mat'), {'data':data})

    # Create .tif output
    elif(filetype == "tif"):
        # Determine the realname based on whether or not a single image is being processed or a pair
        if(filename[3] == "_"):   # pair
            realname = filename[:31]
        else:
            realname = filename[:23]

        # Read geo location information in from a text or geotiff file depending on outtype
        if(outtype == "stand_height"):
            (cols, rows, corner_lat, corner_long, lat_step, long_step) = read_geo_data(realname + "_geo.txt", directory)
            corner_lat = coords[0]
            corner_long = coords[2]
            lat_step = -2.77777777778 * (10**-4)
            long_step = 2.77777777778 * (10**-4)
        elif(outtype == "diff_height"):
            (cols, rows, corner_lat, corner_long, lat_step, long_step) = read_geo_data(ref_file, directory[:-10])
            selffile_data = sio.loadmat(directory[:-10] + "output/" + "self.mat")
            image1 = selffile_data['I1']
            cols = int(image1.shape[0])
            rows = int(image1.shape[1])
            lat_step = -2.77777778 * (10**-4)
            long_step = 2.77777778 * (10**-4)

        # Create the GeoTiff
        driver = gdal.GetDriverByName('GTiff')
        outRaster = driver.Create(directory + outfilename + ".tif", cols, rows)
        outRaster.SetGeoTransform([corner_long, long_step, 0, corner_lat, 0, lat_step])
        outband = outRaster.GetRasterBand(1)
        outband.WriteArray(data)
        outRasterSRS = osr.SpatialReference()
        outRasterSRS.ImportFromEPSG(4326)
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()

    else:
        # Error message
        print ("Error: The selected file type is invalid. Please try again and choose a different output format.")
        print ("You selected %s" % filetype)
        print ("File types available: .gif, .json, .kml, .mat, .tif -- input without the ., such as kml instead of .kml\n")


# In[ ]:


def write_mapfile_new(scenes, flagfile, maskfile, directory, output_files):
    for i in range(scenes):

        # Set the filename
        scene_data = flag_scene_file(flagfile, i + 1, directory) # 0 vs 1 indexing
        filename = scene_data[1]
        image_folder = "f" + scene_data[4] + "_o" + scene_data[5] + "/"

        # Load first image file and associated parameters

        file1 = sio.loadmat(os.path.join(directory, image_folder, filename + "_orig.mat"))
        corr_vs = file1['corr_vs']
        coords = file1['coords'][0]


        # Load and read data from temp .json files
        file_tempD = open(os.path.join(directory, image_folder, filename + "_tempD.json"))
        B = json.load(file_tempD)

        # Set S and C paramters based on the default and data from B
        S_param = 0.65 + B[0]
        C_param = 13 + B[1]

        # Run interpolation to calculate the heights
        gamma = corr_vs.copy()
        gamma = gamma / S_param
        height = arc_sinc(gamma, C_param)
        height[np.isnan(gamma)] = np.nan

        # Mask out non-forest areas
        if maskfile != '-':
            forest_only_height = remove_nonforest(height, coords, maskfile, directory)
        else:
            forest_only_height = height

        # Transpose height to correctly align it (ie so it isn't rotated in relation to an underlying map)
        forest_only_height = forest_only_height.transpose()

        # Get rid of NaN so future processing software doesn't error
        forest_only_height[np.isnan(forest_only_height)] = 255

        #pdb.set_trace()

        # Write all the desired output file types for the forest height map
        for filetype in output_files:
            write_file_type(forest_only_height, "stand_height", filename, os.path.join(directory , image_folder), filetype, coords)

    print ("all tree height map files written at "+ (time.strftime("%H:%M:%S")))

