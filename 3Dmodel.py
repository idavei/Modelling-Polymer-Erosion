"""

    Copyright [2022] [Ishen Dave]

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/

"""

#3-DIMENSIONAL MODEL###########################################################

'Importing required modules'
import numpy as np
import random
import copy
import math
import matplotlib.pyplot as plt
import scipy.stats as stats
from matplotlib import cm 
from operator import sub

#INITIALSING###################################################################

'Initialising conditions for the model'
#Relative proportion of crystalline to amorphous pixels - these
#are fixed at 35:65 for p(CPP-SA)-20:80 (x-ray diffraction)
crys = 0.35
amorph = 1-crys
#Relative rates of erosion for crystalline and amorphous pixels
crys_rate = 0.012 #slow-eroding ordered region
amorph_rate = 1 #fast-eroding tangled region

#FUNCTIONS#####################################################################

'Function to construct the polymer via a 3-D matrix cube'
def grid(thickness,length,height):
    
    grid_crys = 0 #Count for number of crystalline pixels in the grid
    grid_crystallinity = 0
        
    #Initialising Boundaries for polymer matrix cube   
    water = 10*np.zeros([1,length,height]) #Water boundary
    front = 10*np.ones([1,length,height]) #Front boundary
    left = 10*np.ones([thickness+2,1,height]) #Left boundary
    right = 10*np.ones([thickness+2,1,height]) #Right boundary
    upper = 10*np.ones([thickness+2,length+2,1]) #Upper boundary
    lower = 10*np.ones([thickness+2,length+2,1]) #Lower boundary
    
    # #'White-noise' assignment - Primary spatial allocation of amorphous and
    # #crystalline regions randomly throughout the modelled polymeric structure
    
    # init = np.random.rand(thickness,length,height)
    # for i in range(0,thickness):
    #     for j in range(0,length):
    #         for k in range(0,height):
    #             if init[i,j,k] <= crys:
    #                 init[i,j,k] = 5 #Crystalline pixel denoted by a '5'
    #                 grid_crys += 1
    #             else: 
    #                 init[i,j,k] = 8 #Amorphous pixel denoted by an '8'
    
    # grid_crystallinity = grid_crys/(thickness*length*height)
    # pixels = init
    
    #Grid of pixels for spatial-complexity assignment
    #Initialise a random matrix of the dimension we require
    pixels = np.random.rand(thickness,length,height)
    
    #Adding boundaries to pixel grid    
    pixels = np.row_stack((water,pixels))
    pixels = np.row_stack((pixels,front))
    pixels = np.column_stack((left,pixels))
    pixels = np.column_stack((pixels,right))
    pixels = np.dstack((pixels,upper))
    pixels = np.dstack((lower,pixels))
    
    #Spherulites
    max_radius = 12
    min_radius = 8
    
    #Generating a set of sample radii from a truncated normal distribution
    mean = (max_radius+min_radius)/2
    std = (mean/10)
    lower_bound = (min_radius-mean)/std
    upper_bound = (max_radius-mean)/std
    truncated = stats.truncnorm(lower_bound,upper_bound,mean,std)
    sample_radii = truncated.rvs(10000)
    
    # #Basic spherulite spatial-complexity assignment - Initial representation of
    # #crystalline regions as spherulites embedded amongst a 'sea' of amorphous material
    
    # while grid_crystallinity <= crys: #Controlling crystallinity of grid 
    #     t = random.randint(1,thickness)
    #     l = random.randint(1,length)
    #     h = random.randint(1,height)
    #     pixel_mid = [t+0.5,l+0.5,h+0.5]
        
    #     radius = random.uniform(min_radius,max_radius)
    #     for i in range(t-math.ceil(radius),t+math.ceil(radius)+1): 
    #         for j in range(l-math.ceil(radius),l+math.ceil(radius)+1):
    #             for k in range(h-math.ceil(radius),h+math.ceil(radius)+1):
    #                 if 0 < i <= thickness: #No overlap onto or around front, back boundaries
    #                     gen_mid = [i+0.5,j+0.5,k+0.5]
    #                     difference = list(map(sub,pixel_mid,gen_mid))
    #                     y = copy.copy(j)
    #                     z = copy.copy(k)
    #                     #Euclidean distance
    #                     if np.linalg.norm(difference) <= radius:
    #                           #Periodic boundary condition on left and right boundaries
    #                         if y == 0:
    #                             y = y + length
    #                         if y < 0 or y >= length+1:
    #                             y = y%length
    #                         #Periodic boundary condition on lower and upper boundaries
    #                         if z == 0:
    #                             z = z + height
    #                         if z < 0 or z >= height+1:
    #                             z = z%height
                                
    #                         if pixels[i,y,z] != 5:
    #                             pixels[i,y,z] = 5
    #                             grid_crys += 1
    
    # Concentric banded spherulite spatial-complexity assignment - Final representation of
    # a spherulite's interior structure as alternating bands of crystalline and amorphous 
    # material emanating from a central crystalline point    
    
    #Controlling crystallinity of grid 
    while grid_crystallinity <= crys: 
        #Selecting a random spherulite centroid
        t = random.randint(1,thickness)
        l = random.randint(1,length)
        h = random.randint(1,height)
        pixel_mid = [t+0.5,l+0.5,h+0.5]
        centroid_euclideans = []
        
        #Checking Euclidean distance from potential centroid to all existing centroids
        for c in range(0,len(centroids)):
            centroid_manhattan = list(map(sub,pixel_mid,centroids[c]))
            #Finding shortest distances between centroid dimensions 'length' and 'height'
            #Length considering periodic boundary
            delta_l = abs(centroid_manhattan[1])
            delta_lperiodic = length-abs(centroid_manhattan[1])
            length_diff = min(delta_l,delta_lperiodic)
            #Height considering periodic boundary
            delta_h = abs(centroid_manhattan[2])
            delta_hperiodic = height-abs(centroid_manhattan[2])
            height_diff = min(delta_h,delta_hperiodic)
            
            centroid_difference = [centroid_manhattan[0],length_diff,height_diff]
            centroid_euclidean = np.linalg.norm(centroid_difference)
            centroid_euclideans.append(centroid_euclidean)
        
        #If the Euclidean distance to all existing centroids is larger than a defined
        #number, accept centroid and proceed with growth
        if all([x >= 15 for x in centroid_euclideans]): 
            centroids.append(pixel_mid)
            
            # radius = random.uniform(min_radius,max_radius)
            radius = random.choice(sample_radii)
            radii.append(radius)
            for i in range(t-math.ceil(radius),t+math.ceil(radius)+1): 
                for j in range(l-math.ceil(radius),l+math.ceil(radius)+1):
                    for k in range(h-math.ceil(radius),h+math.ceil(radius)+1):
                        if 0 < i <= thickness: #No overlap onto or around left, right boundaries
                            gen_mid = [i+0.5,j+0.5,k+0.5]
                            difference = list(map(sub,pixel_mid,gen_mid))
                            y = copy.copy(j)
                            z = copy.copy(k)
                            #Even radii
                            if math.ceil(radius)%2 == 0:
                                for m in reversed(range(1,(math.ceil(radius/2)+1))):
                                    #Adding concentric bands of crystalline lamellae between
                                    #amorphous regions
                                    if radius-(2*m) <= np.linalg.norm(difference) < radius+(1-(2*m)):
                                        #Periodic boundary condition on left and right boundaries
                                        if y == 0:
                                            y = y + length
                                        if y < 0 or y >= length+1:
                                            y = y%length
                                        #Periodic boundary condition on lower and upper boundaries
                                        if z == 0:
                                            z = z + height
                                        if z < 0 or z >= height+1:
                                            z = z%height        
                                        if pixels[i,y,z] == 5 or pixels[i,y,z] == 3:
                                            break
                                        else:
                                            pixels[i,y,z] = 5
                                            grid_crys += 1
                                            
                                    elif radius-((2*m)+1) <= np.linalg.norm(difference) < radius-(2*m):
                                        #Periodic boundary condition on left and right boundaries
                                        if y == 0:
                                            y = y + length
                                        if y < 0 or y >= length+1:
                                            y = y%length
                                        #Periodic boundary condition on lower and upper boundaries
                                        if z == 0:
                                            z = z + height
                                        if z < 0 or z >= height+1:
                                            z = z%height    
                                        if not pixels[i,y,z] == 5 or pixels[i,y,z] == 3:
                                            pixels[i,y,z] = 3
                                                                        
                            #Odd radii
                            else:
                                for m in reversed(range(0,math.floor(radius/2)+1)):
                                    #Adding concentric bands of crystalline lamellae between amorphous regions
                                    if radius-((2*m)+1)<=np.linalg.norm(difference)<radius-(2*m):
                                        #Periodic boundary condition on left and right boundaries
                                        if y == 0:
                                            y = y + length
                                        if y < 0 or y >= length+1:
                                            y = y%length
                                        #Periodic boundary condition on lower and upper boundaries
                                        if z == 0:
                                            z = z + height
                                        if z < 0 or z >= height+1:
                                            z = z%height 
                                        if pixels[i,y,z] == 5 or pixels[i,y,z] == 3:
                                            break
                                        else:
                                            pixels[i,y,z] = 5
                                            grid_crys += 1
                                            
                                for n in reversed(range(1,math.floor(radius/2)+1)):
                                    if radius-(2*n)<=np.linalg.norm(difference)<radius+(1-(2*n)):
                                        #Periodic boundary condition on left and right boundaries
                                        if y == 0:
                                            y = y + length
                                        if y < 0 or y >= length+1:
                                            y = y%length
                                        #Periodic boundary condition on lower and upper boundaries
                                        if z == 0:
                                            z = z + height
                                        if z < 0 or z >= height+1:
                                            z = z%height
                                        if not pixels[i,y,z] == 5 or pixels[i,y,z] == 3:
                                            pixels[i,y,z] = 3
                                        
        grid_crystallinity = grid_crys/(thickness*length*height)
        
    for i in range(1,thickness+1): #Vertically moving through rows
        for j in range(1,length+1): #Horizontally moving through columns
            for k in range(1,height+1):
                if pixels[i,j,k] != 5:
                    pixels[i,j,k] = 8 #Amorphous pixel denoted by an '8' 
    
    # print("{} {}".format("Grid crystallinity of", grid_crystallinity))
    # print(grid_crys)
    
    return pixels


'Function that finds vulnerable pixels in the polymer grid'
def is_vulnerable(thickness,length,height):
    
    #Using information from last dissolved pixel to update vulnerable pixels
    step = dissolved[-1]
    vulnerable_pixels.remove(step)
    
    back = [step[0]-1,step[1],step[2]]
    if pixels[back[0],back[1],back[2]]!=0:
        if vulnerable_pixels.count(back) < 1:
            vulnerable_pixels.append(back)
        
    front = [step[0]+1,step[1],step[2]]
    if pixels[front[0],front[1],front[2]]!=0 and pixels[front[0],front[1],front[2]]!=10:
        if vulnerable_pixels.count(front) < 1:
            vulnerable_pixels.append(front)
    
    left = [step[0],step[1]-1,step[2]]
    if pixels[left[0],left[1],left[2]]!=0 and pixels[left[0],left[1],left[2]]!=10:
        if vulnerable_pixels.count(left) < 1:
            vulnerable_pixels.append(left)
    
    right = [step[0],step[1]+1,step[2]]
    if pixels[right[0],right[1],right[2]]!=0 and pixels[right[0],right[1],right[2]]!=10:
        if vulnerable_pixels.count(right) < 1:
            vulnerable_pixels.append(right)
                
    up = [step[0],step[1],step[2]+1]
    if pixels[up[0],up[1],up[2]]!=0 and pixels[up[0],up[1],up[2]]!=10:
        if vulnerable_pixels.count(up) < 1:
            vulnerable_pixels.append(up)
                    
    down = [step[0],step[1],step[2]-1]
    if pixels[down[0],down[1],down[2]]!=0 and pixels[down[0],down[1],down[2]]!=10:
        if vulnerable_pixels.count(down) < 1:
            vulnerable_pixels.append(down)
            
    #Periodic boundary conditions on length
    if pixels[left[0],left[1],left[2]] == 10:
        p = [left[0],left[1]+length,left[2]]
        if pixels[p[0],p[1],p[2]]!=0:
            if vulnerable_pixels.count(p) < 1:
                vulnerable_pixels.append(p)
          
    if pixels[right[0],right[1],right[2]] == 10:
        q = [right[0],right[0]%length,right[2]]
        if pixels[q[0],q[1],q[2]]!=0:
            if vulnerable_pixels.count(q) < 1:
                vulnerable_pixels.append(q)
    
    #Periodic boundary conditions on height
    if pixels[up[0],up[1],up[2]] == 10:
        r = [up[0],up[1],up[2]%height]
        if pixels[r[0],r[1],r[2]]!=0:
            if vulnerable_pixels.count(r) < 1:
                vulnerable_pixels.append(r)
                
    if pixels[down[0],down[1],down[2]] == 10:
        s = [down[0],down[1],down[2]+height]
        if pixels[s[0],s[1],s[2]]!=0:
            if vulnerable_pixels.count(s) < 1:
                vulnerable_pixels.append(s)

'Function determining whether a specific vulnerable pixel...'
'...is amorphous or crystalline'
def vulnerable_classifier():
        
    for c in range(0,len(vulnerable_pixels)):
        item = vulnerable_pixels[c]
        if pixels[item[0],item[1],item[2]] == 5:
            crystalline.append(item)
        elif pixels[item[0],item[1],item[2]] == 8:
            amorphous.append(item)
            
'Function to choose and subsequently dissolve a pixel'
def dissolve():
        
    chance_class = random.uniform(0,1)
    if chance_class < amorph_prob:
        dissolve = random.choice(amorphous)
        pixels[dissolve[0],dissolve[1],dissolve[2]] = 0
    else:
        dissolve = random.choice(crystalline)
        pixels[dissolve[0],dissolve[1],dissolve[2]] = 0
         
    dissolved.append(dissolve)


#SIMULATION####################################################################

'TIME ZERO DATA'
'Initialising lists to be used'
#List of currently vulnerable pixels in the polymer grid
vulnerable_pixels = []
#List of currently vulnerable crystalline pixels
crystalline = []
#List of currently vulnerable amorphous pixels
amorphous = []
#Time series of polymer grid through simulation
time_series = []
#Time count through simulation
sim = []
#Time count of all steps
time_steps = []
#Step count
step_count = []
#List of dissolved pixels
dissolved = []
#List of spherulite centroids
centroids = []
#List of spherulite radii
radii = []

'Generating initial data for the simulation'
time = 0
sim.append(time)
time_steps.append(time)
step = 0
step_count.append(step)
thickness = 35
length = thickness
height = thickness
print("Step", step)
print("Time", time)
print()
print("Generating matrix...")
pixels = grid(thickness,length,height)
stationary = copy.copy(pixels)
time_series.append(stationary)
print()

#Initial condition rendering relevant pixels within first layer as vulnerable to erosion
for j in range(1,length+1):
    for k in range(1, height+1):
        vulnerable_pixels.append([1,j,k])
        
vulnerable_classifier()

'Calculating rates for the first step'
#Dissolution rate for amorphous
amorph_total_rate = amorph_rate*len(amorphous)
#Dissolution rate for crystalline
crys_total_rate = crys_rate*len(crystalline)
#Find total rate of dissolution occurring
total_rate = amorph_total_rate + crys_total_rate 
#Probability of an amorphous pixel dissolving
amorph_prob = amorph_total_rate/total_rate
#Probability of a crystalline pixel dissolving
crys_prob = crys_total_rate/total_rate
#Generate time to next event via Exponential distribution
#mu_timestep = -np.log(random.uniform(0,1))/total_rate
mu_timestep = np.random.exponential(1/total_rate)

'Gillespie Algorithm'
total_time = 7
while time < total_time:
    step += 1
    step_count.append(step)
    time += mu_timestep
    time_steps.append(time)
    print("Step", step)
    print("Time", time)
    print("Simulation terminates after Time", total_time)
    print()
    
    #Choosing pixels to dissolve and then dissolving them
    dissolve()
    
    #Saving grid composition at strategic steps to prevent out of memory error
    for m in range(1,2000):
        if time > m and time_steps[-2] < m:
            if sim.count(time) < 1:
                sim.append(time)
                stationary = copy.copy(pixels)
                time_series.append(stationary)

    #Clearing lists
    # vulnerable_pixels.clear()
    crystalline.clear()
    amorphous.clear()
    
    #Identifying and classifying the updated vulnerable pixels
    is_vulnerable(thickness,length,height)
    vulnerable_classifier()
    
    #Updating the new rates
    amorph_total_rate = amorph_rate*len(amorphous)
    crys_total_rate = crys_rate*len(crystalline)
    
    total_rate = amorph_total_rate + crys_total_rate 
    if total_rate == 0:
        break
    else:
        amorph_prob = amorph_total_rate/total_rate
        crys_prob = crys_total_rate/total_rate
    #mu_timestep = -np.log(random.uniform(0,1))/total_rate
    mu_timestep = np.random.exponential(1/total_rate)

#PLOTTING######################################################################

fig = plt.figure()
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')
erosion = cm.get_cmap('jet', 10) #Adjusts colorbar

# #Boundaries
# coords = np.arange(0,length+2,1) #Array of coordinates for each axis
# X,Y,Z = np.meshgrid(coords,coords,coords)

#No boundaries
coords = np.arange(0,length,1) #Array of coordinates for each axis
X,Y,Z = np.meshgrid(coords,coords,coords)
for a in range(0,3):
    pixels = np.delete(pixels,0,a)
for b in range(0,3):
    pixels = np.delete(pixels,length,b)

grid = ax.scatter(X,Y,Z,c=pixels,depthshade=False,cmap=erosion,vmin=0,vmax=10,alpha=0.7)
# plt.colorbar(grid)
ax.set_xlabel('Length (\u03BCm)',fontsize=15)
ax.set_ylabel('Thickness (\u03BCm)',fontsize=15)
ax.set_zlabel('Height (\u03BCm)',fontsize=15)
ax.view_init(elev=24, azim=-30)
ax.set_title('Polymer matrix after Time = {}'.format(total_time))
plt.show()