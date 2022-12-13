"""

    Copyright [2022] [Ishen Dave]

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/

"""

#2-DIMENSIONAL MODEL###########################################################

'Importing required modules'
import numpy as np
import random
import copy
import math
import matplotlib.pyplot as plt
import scipy.stats as stats
from matplotlib import cm 
from operator import sub
from itertools import chain

#INITIALSING###################################################################

'Initialising conditions for the model'
#Relative proportion of crystalline and amorphous regions - these
#are fixed at 35:65 for p(CPP-SA)-20:80 (x-ray diffraction)
crys = 0.35
amorph = 1-crys
#Relative rates of erosion for crystalline and amorphous pixels
crys_rate = 0.012 #slow-eroding ordered region
amorph_rate = 1 #fast-eroding tangled region

#2-DIMENSIONAL MODEL###########################################################
#FUNCTIONS#####################################################################

'Function to construct the polymer via a 2-D grid matrix'
def grid(length,thickness):
   
    grid_crys = 0 #Count for number of crystalline pixels in the grid
    grid_crystallinity = grid_crys/(length*thickness)
    
    #Initialising Boundaries            
    upper = 10*np.ones([1,length+2]) #Upper boundary  
    lower = 10*np.ones([1,length+2]) #Lower boundary  
    right = 10*np.ones([thickness,1]) #Right boundary  
    water = np.zeros([thickness,1]) #Water(left) boundary
    
    # #'White-noise' assignment - Primary spatial allocation of amorphous and
    # #crystalline regions randomly throughout the modelled polymeric structure

    # init = np.random.rand(length,thickness)
    # for i in range(0,length): #Vertically moving through rows
    #     for j in range(0,thickness): #Horizontally moving through columns
    #         if init[i,j] <= crys:
    #             init[i,j] = 5 #Crystalline pixel denoted by a '5'
    #             grid_crys += 1
    #         else: 
    #             init[i,j] = 8 #Amorphous pixel denoted by an '8'
    
    # grid_crystallinity = grid_crys/(length*thickness)
    # pixels = init
                   
    #Grid of pixels for spatial-complexity assignment
    #Initialise a random matrix of the dimension we require
    pixels = np.random.rand(length,thickness)

    #Adding boundaries to pixel grid
    pixels = np.column_stack((water,pixels))
    pixels = np.column_stack((pixels,right))
    pixels = np.row_stack((upper,pixels))
    pixels = np.row_stack((pixels,lower))
    
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
    #     l = random.randint(1,length)
    #     t = random.randint(1,thickness)
    #     pixel_mid = [l+0.5,t+0.5]
    #     centroids.append(pixel_mid)
        
    #     # radius = random.uniform(min_radius,max_radius)
    #     radius = random.choice(sample_radii)
    #     radii.append(radius)
    #     for i in range(l-math.ceil(radius),l+math.ceil(radius)+1): 
    #         for j in range(t-math.ceil(radius),t+math.ceil(radius)+1):
    #           if 0 < j <= thickness: #No overlap onto or around left, right boundaries
    #                 gen_mid = [i+0.5,j+0.5]
    #                 difference = list(map(sub,pixel_mid,gen_mid))
    #                 x = copy.copy(i)
    #                 #Euclidean distance
    #                 if np.linalg.norm(difference) <= radius:
    #                     #Periodic boundary condition on lower and upper boundaries
    #                     if x == 0:
    #                         x = x + length
    #                     if x < 0 or x >= length+1:
    #                         x = x%length                
    #                     if pixels[x,j] != 5:
    #                         pixels[x,j] = 5
    #                         grid_crys += 1
                            
    #Concentric banded spherulite spatial-complexity assignment - Final representation of
    #a spherulite's interior structure as alternating bands of crystalline and amorphous 
    #material emanating from a central crystalline point

    #Controlling crystallinity of grid 
    while grid_crystallinity <= crys: 
        #Selecting a potential spherulite centroid at random
        l = random.randint(1,length)
        t = random.randint(1,thickness)
        pixel_mid = [l+0.5,t+0.5]       
        centroid_manhattans = []
        centroid_differences = []
        
        #Checking Euclidean distance from potential centroid to all existing centroids
        for c in range(0,len(centroids)):
            centroid_manhattan = list(map(sub,pixel_mid,centroids[c]))
            #Shortest distance across periodic boundary
            if abs(centroid_manhattan[0]) > length/2:
                periodic_manhattan = [length-abs(centroid_manhattan[0]),centroid_manhattan[1]]
                centroid_manhattans.append(periodic_manhattan)
                centroid_difference = np.linalg.norm(periodic_manhattan)
            #Shortest distance across polymer bulk
            else:
                centroid_difference = np.linalg.norm(centroid_manhattan)
                centroid_manhattans.append(centroid_manhattan)
            centroid_differences.append(centroid_difference)

        #If the Euclidean distance to all existing centroids is larger than a defined
        #number, accept centroid and proceed with growth
        if all([x >= 15 for x in centroid_differences]): 
            centroids.append(pixel_mid)

            # radius = random.uniform(min_radius,max_radius)
            radius = random.choice(sample_radii)
            radii.append(radius)
            for i in range(l-math.ceil(radius),l+math.ceil(radius)+1): 
                for j in range(t-math.ceil(radius),t+math.ceil(radius)+1):
                    if 0 < j <= thickness: #No overlap onto or around left, right boundaries
                        gen_mid = [i+0.5,j+0.5]
                        difference = list(map(sub,pixel_mid,gen_mid))
                        x = copy.copy(i)
                        #Even radii
                        if math.ceil(radius)%2 == 0:
                            for m in reversed(range(1,(math.ceil(radius/2)+1))):
                                #Adding concentric bands of crystalline lamellae between
                                #amorphous regions
                                if radius-(2*m) <= np.linalg.norm(difference) < radius+(1-(2*m)):
                                    #Periodic boundary condition on lower and upper boundaries
                                    if x == 0:
                                        x = x + length
                                    if x < 0 or x >= length+1:
                                        x = x%length
                                    if pixels[x,j] == 5 or pixels[x,j] == 3:
                                        break
                                    else:
                                        pixels[x,j] = 5
                                        grid_crys += 1
                                        
                                elif radius-((2*m)+1) <= np.linalg.norm(difference) < radius-(2*m):
                                    if x == 0:
                                        x = x + length
                                    if x < 0 or x >= length+1:
                                        x = x%length
                                    if not pixels[x,j] == 5 or pixels[x,j] == 3:
                                        pixels[x,j] = 3
                            
                            if np.linalg.norm(difference) == radius:
                                if x == 0:
                                    x = x + length
                                if x < 0 or x >= length+1:
                                    x = x%length
                                if pixels[x,j] == 5 or pixels[x,j] == 3:
                                    break
                                else:
                                    pixels[x,j] = 5
                                    grid_crys += 1
                                                                    
                        #Odd radii
                        else:
                            if isinstance(radius, int) == True:
                                if np.linalg.norm(difference) == 0:       
                                    if x == 0:
                                        x = x + length
                                    if x < 0 or x >= length+1:
                                        x = x%length
                                    if pixels[x,j] == 5 or pixels[x,j] == 3:
                                        break
                                    else:
                                        pixels[x,j] = 5
                                        grid_crys += 1
                                            
                            for m in reversed(range(0,math.floor(radius/2)+1)):
                                #Adding concentric bands of crystalline lamellae between amorphous regions
                                if radius-(1+(2*m)) < np.linalg.norm(difference) <= radius-(2*m):
                                    #Periodic boundary condition on lower and upper boundaries
                                    if x == 0:
                                        x = x + length
                                    if x < 0 or x >= length+1:
                                        x = x%length             
                                    if pixels[x,j] == 5 or pixels[x,j] == 3:
                                        break
                                    else:
                                        pixels[x,j] = 5
                                        grid_crys += 1
                                        
                            for n in reversed(range(1,math.floor(radius/2)+1)):
                                if radius-(2*n) < np.linalg.norm(difference) <= radius+(1-(2*n)):
                                    if x == 0:
                                        x = x + length
                                    if x < 0 or x >= length+1:
                                        x = x%length
                                    if not pixels[x,j] == 5 or pixels[x,j] == 3:
                                        pixels[x,j] = 3
                                    
        grid_crystallinity = grid_crys/(length*thickness)
                         
    for i in range(1,length+1):
        for j in range(1,thickness+1): 
            if pixels[i,j] != 5:
                    pixels[i,j] = 8 #Amorphous pixel denoted by an '8'  
    
    # print("{} {}".format("Grid crystallinity of", grid_crystallinity))
    # print(grid_crys)
          
    return pixels


'Function that finds vulnerable pixels in the polymer grid, and also classifies...'
'...how many acidic water units are adjacent in the self-inhibition model'
def is_vulnerable(length,thickness):
    
    # #Naive searcher
    # for i in range(1,length+1): 
    #     for j in range(1,thickness+1):
    #         left = pixels[i,j-1]
    #         right = pixels[i,j+1]
    #         up = pixels[i-1,j]
    #         down = pixels[i+1,j]
    #         if left*right*up*down == 0 and pixels[i,j]!=0:
    #             vulnerable_pixels.append([i,j])
    
    #Using information from last dissolved pixel to update vulnerable pixels
    step = dissolved[-1]
    vulnerable_pixels.remove(step)
    if step in one_water:
        one_water.remove(step)
    if step in two_water:
        two_water.remove(step)
    if step in three_water:
        three_water.remove(step)
    if step in four_water:
        four_water.remove(step)
    if step in pure_water:
        pure_water.remove(step)
    if step in pure_one_acidic:
        pure_one_acidic.remove(step)
    if step in pure_two_acidic:
        pure_two_acidic.remove(step)
    if step in pure_three_acidic:
        pure_three_acidic.remove(step)
        
    left = [step[0],step[1]-1]
    if pixels[left[0],left[1]]!=0:   
        if vulnerable_pixels.count(left) < 1:
            vulnerable_pixels.append(left)
            one_water.append(left)
        elif vulnerable_pixels.count(left) == 1:
            if step[1] == 2:
                if pure_two_acidic.count(left) == 1:    
                    pure_two_acidic.remove(left)
                    pure_three_acidic.append(left)
                elif pure_one_acidic.count(left) == 1:
                    pure_one_acidic.remove(left)
                    pure_two_acidic.append(left)
                elif pure_one_acidic.count(left) < 1:
                    pure_one_acidic.append(left)
                    pure_water.remove(left)
            elif one_water.count(left) == 1:
                one_water.remove(left)
                two_water.append(left)
            elif two_water.count(left) == 1:
                two_water.remove(left)
                three_water.append(left)
            elif three_water.count(left) == 1:
                three_water.remove(left)
                four_water.append(left)
                
    right = [step[0],step[1]+1]
    if pixels[right[0],right[1]]!=0 and pixels[right[0],right[1]]!=10:
        if vulnerable_pixels.count(right) < 1:
            vulnerable_pixels.append(right)
            one_water.append(right)
        elif vulnerable_pixels.count(right) == 1:
            if one_water.count(right) == 1:
                one_water.remove(right)
                two_water.append(right)
            elif two_water.count(right) == 1:
                two_water.remove(right)
                three_water.append(right)
            elif three_water.count(right) == 1:
                three_water.remove(right)
                four_water.append(right)
                  
    down = [step[0]-1,step[1]]
    if pixels[down[0],down[1]]!=0 and pixels[down[0],down[1]]!=10:
        if vulnerable_pixels.count(down) < 1: 
            vulnerable_pixels.append(down)
            one_water.append(down)
        elif vulnerable_pixels.count(down) == 1:
            if step[1] == 1:
                if pure_two_acidic.count(down) == 1:    
                    pure_two_acidic.remove(down)
                    pure_three_acidic.append(down)
                elif pure_one_acidic.count(down) == 1:
                    pure_one_acidic.remove(down)
                    pure_two_acidic.append(down)
                elif pure_one_acidic.count(down) < 1:
                    pure_one_acidic.append(down)
                    pure_water.remove(down)
            elif one_water.count(down) == 1:
                one_water.remove(down)
                two_water.append(down)
            elif two_water.count(down) == 1:
                two_water.remove(down)
                three_water.append(down)
            elif three_water.count(down) == 1:
                three_water.remove(down)
                four_water.append(down)
        
    up = [step[0]+1,step[1]]
    if pixels[up[0],up[1]]!=0 and pixels[up[0],up[1]]!=10:
        if vulnerable_pixels.count(up) < 1:
            vulnerable_pixels.append(up)
            one_water.append(up)
        elif vulnerable_pixels.count(up) == 1:
            if step[1] == 1:
                if pure_two_acidic.count(up) == 1:    
                    pure_two_acidic.remove(up)
                    pure_three_acidic.append(up)
                elif pure_one_acidic.count(up) == 1:
                    pure_one_acidic.remove(up)
                    pure_two_acidic.append(up)
                elif pure_one_acidic.count(up) < 1:
                    pure_one_acidic.append(up)
                    pure_water.remove(up)
            elif one_water.count(up) == 1:
                one_water.remove(up)
                two_water.append(up)
            elif two_water.count(up) == 1:
                two_water.remove(up)
                three_water.append(up)
            elif three_water.count(up) == 1:
                three_water.remove(up)
                four_water.append(up)
                  
    #Periodic boundary condition across upper boundary   
    if pixels[up[0],up[1]] == 10:
        p = [up[0]%length,up[1]]
        if pixels[p[0],p[1]]!=0:
            if vulnerable_pixels.count(p) < 1:
                vulnerable_pixels.append(p)
                one_water.append(p)
            elif vulnerable_pixels.count(p) == 1:
                if step[1] == 1:
                    if pure_two_acidic.count(p) == 1:    
                        pure_two_acidic.remove(p)
                        pure_three_acidic.append(p)
                    elif pure_one_acidic.count(p) == 1:
                        pure_one_acidic.remove(p)
                        pure_two_acidic.append(p)
                    elif pure_one_acidic.count(p) < 1:
                        pure_one_acidic.append(p)
                        pure_water.remove(p)                
                elif one_water.count(p) == 1:
                    one_water.remove(p)
                    two_water.append(p)
                elif two_water.count(p) == 1:
                    two_water.remove(p)
                    three_water.append(p)
                elif three_water.count(p) == 1:
                    three_water.remove(p)
                    four_water.append(p)
    
    #Periodic boundary condition across lower boundary
    if pixels[down[0],down[1]] == 10:
        q = [down[0]+length,down[1]]
        if pixels[q[0],q[1]]!=0: 
            if vulnerable_pixels.count(q) < 1:
                vulnerable_pixels.append(q)
                one_water.append(q)
            elif vulnerable_pixels.count(q) == 1:
                if step[1] == 1:
                    if pure_two_acidic.count(q) == 1:    
                        pure_two_acidic.remove(q)
                        pure_three_acidic.append(q)
                    elif pure_one_acidic.count(q) == 1:
                        pure_one_acidic.remove(q)
                        pure_two_acidic.append(q)
                    elif pure_one_acidic.count(q) < 1:
                        pure_one_acidic.append(q)
                        pure_water.remove(q)               
                elif one_water.count(q) == 1:
                    one_water.remove(q)
                    two_water.append(q)
                elif two_water.count(q) == 1:
                    two_water.remove(q)
                    three_water.append(q)
                elif three_water.count(q) == 1:
                    three_water.remove(q)
                    four_water.append(q)


'Function determining whether a specific vulnerable pixel...'
'...is amorphous or crystalline'
def vulnerable_classifier():
        
    for c in range(0,len(vulnerable_pixels)):
        item = vulnerable_pixels[c]
        if pixels[item[0],item[1]] == 5:
            crystalline.append(item)
        elif pixels[item[0],item[1]] == 8:
            amorphous.append(item)
    
    #Self-inhibition model        
    for c in range(0,len(one_water)):
        item = one_water[c]
        if pixels[item[0],item[1]] == 5:
            one_water_crystalline.append(item)
        elif pixels[item[0],item[1]] == 8:
            one_water_amorphous.append(item)
            
    for c in range(0,len(two_water)):
        item = two_water[c]
        if pixels[item[0],item[1]] == 5:
            two_water_crystalline.append(item)
        elif pixels[item[0],item[1]] == 8:
            two_water_amorphous.append(item)
        
    for c in range(0,len(three_water)):
        item = three_water[c]
        if pixels[item[0],item[1]] == 5:
            three_water_crystalline.append(item)
        elif pixels[item[0],item[1]] == 8:
            three_water_amorphous.append(item)
            
    for c in range(0,len(four_water)):
        item = four_water[c]
        if pixels[item[0],item[1]] == 5:
            four_water_crystalline.append(item)
        elif pixels[item[0],item[1]] == 8:
            four_water_amorphous.append(item)
    
    for c in range(0,len(pure_water)):
        item = pure_water[c]
        if pixels[item[0],item[1]] == 5:
            pure_water_crystalline.append(item)
        elif pixels[item[0],item[1]] == 8:
            pure_water_amorphous.append(item)
    
    for c in range(0,len(pure_one_acidic)):
        item = pure_one_acidic[c]
        if pixels[item[0],item[1]] == 5:
            pure_one_acidic_crystalline.append(item)
        elif pixels[item[0],item[1]] == 8:
            pure_one_acidic_amorphous.append(item)
            
    for c in range(0,len(pure_two_acidic)):
        item = pure_two_acidic[c]
        if pixels[item[0],item[1]] == 5:
            pure_two_acidic_crystalline.append(item)
        elif pixels[item[0],item[1]] == 8:
            pure_two_acidic_amorphous.append(item)
    
    for c in range(0,len(pure_three_acidic)):
        item = pure_three_acidic[c]
        if pixels[item[0],item[1]] == 5:
            pure_three_acidic_crystalline.append(item)
        elif pixels[item[0],item[1]] == 8:
            pure_three_acidic_amorphous.append(item)
            

'Function to choose and subsequently dissolve a pixel'
def dissolve():
    
    chance_class = random.uniform(0,1)
    if chance_class < amorph_prob: #Dissolve an amorphous pixel
        # #Self-inhibition model
        # amorphous_list = []
        # amorphous_choice = []
        
        # if pure_water_amorphous or pure_one_acidic_amorphous or pure_two_acidic_amorphous or pure_three_acidic_amorphous: 
        #     while len(amorphous_list) == 0:
        #         amorphous_choice = random.choices([pure_water_amorphous,pure_one_acidic_amorphous,
        #                     pure_two_acidic_amorphous,pure_three_acidic_amorphous,one_water_amorphous,
        #                     two_water_amorphous,three_water_amorphous,four_water_amorphous],
        #                     [p0,p1,p2,p3,a1,a2,a3,a4], k = 1)
        #         amorphous_list = list(chain.from_iterable(amorphous_choice))
                
        #     dissolve = random.choice(amorphous_list)
        
        # else:
        #     while len(amorphous_list) == 0:
        #         amorphous_choice = random.choices([one_water_amorphous,
        #                     two_water_amorphous,three_water_amorphous,four_water_amorphous],
        #                     [a1,a2,a3,a4], k = 1)
        #         amorphous_list = list(chain.from_iterable(amorphous_choice))
                
        #     dissolve = random.choice(amorphous_list)
        
        #No self-inhibition model
        dissolve = random.choice(amorphous)
        
        pixels[dissolve[0],dissolve[1]] = 0
        
    else: #Dissolve a crystalline pixel
        # #Self-inhibition model
        # crystalline_list = []
        # crystalline_choice = []
        
        # if pure_water_crystalline or pure_one_acidic_crystalline or pure_two_acidic_crystalline or pure_three_acidic_crystalline:
        #     while len(crystalline_list) == 0:
        #         crystalline_choice = random.choices([pure_water_crystalline,pure_one_acidic_crystalline,
        #                     pure_two_acidic_crystalline,pure_three_acidic_crystalline,one_water_crystalline,
        #                     two_water_crystalline,three_water_crystalline,four_water_crystalline],
        #                     [p0,p1,p2,p3,a1,a2,a3,a4], k = 1)
        #         crystalline_list = list(chain.from_iterable(crystalline_choice))
              
        #     dissolve = random.choice(crystalline_list)
        
        # else:
        #     while len(crystalline_list) == 0:
        #         crystalline_choice = random.choices([one_water_crystalline,
        #                     two_water_crystalline,three_water_crystalline,four_water_crystalline],
        #                     [a1,a2,a3,a4], k = 1)
        #         crystalline_list = list(chain.from_iterable(crystalline_choice))
              
        #     dissolve = random.choice(crystalline_list)
        
        #No self-inhibition model
        dissolve = random.choice(crystalline)
        
        pixels[dissolve[0],dissolve[1]] = 0
        
    dissolved.append(dissolve)
        
    return pixels


'Function to clear lists'
def clear_lists():

    crystalline.clear()
    amorphous.clear()
    one_water_amorphous.clear()
    one_water_crystalline.clear()
    two_water_amorphous.clear()
    two_water_crystalline.clear()
    three_water_amorphous.clear()
    three_water_crystalline.clear()
    four_water_amorphous.clear()
    four_water_crystalline.clear()
    pure_water_amorphous.clear()
    pure_water_crystalline.clear()
    pure_one_acidic_amorphous.clear()
    pure_one_acidic_crystalline.clear()
    pure_two_acidic_amorphous.clear()
    pure_two_acidic_crystalline.clear()
    pure_three_acidic_amorphous.clear()
    pure_three_acidic_crystalline.clear()
    

'Function to compute total rate of dissolution for all...'
'...classes of amorphous pixels in the self-inhibition model'
def amorph_total_rates():
    
    rate_count = 0
    amorph_total = 0
    
    if pure_water_amorphous:
        rate_count = p0*len(pure_water_amorphous) + rate_count
    if pure_one_acidic_amorphous:
        rate_count = p1*len(pure_one_acidic_amorphous) + rate_count
    if pure_two_acidic_amorphous:
        rate_count = p2*len(pure_two_acidic_amorphous) + rate_count
    if pure_three_acidic_amorphous:
        rate_count = p3*len(pure_three_acidic_amorphous) + rate_count
    if one_water_amorphous:
        rate_count = a1*len(one_water_amorphous) + rate_count
    if two_water_amorphous:
        rate_count = a2*len(two_water_amorphous) + rate_count
    if three_water_amorphous:
        rate_count = a3*len(three_water_amorphous) + rate_count
    if four_water_amorphous:
        rate_count = a4*len(four_water_amorphous) + rate_count
    
    if rate_count == 0:
        pass
    else:
        amorph_total = amorph_rate*rate_count
    
    return amorph_total


'Function to compute the total rate of dissolution for all...'
'...classes of crystalline pixels in the self-inhibition model'   
def crys_total_rates():
    
    rate_count = 0
    crys_total = 0
    
    if pure_water_crystalline:
        rate_count = p0*len(pure_water_crystalline) + rate_count
    if pure_one_acidic_crystalline:
        rate_count = p1*len(pure_one_acidic_crystalline) + rate_count
    if pure_two_acidic_crystalline:
        rate_count = p2*len(pure_two_acidic_crystalline) + rate_count
    if pure_three_acidic_crystalline:
        rate_count = p3*len(pure_three_acidic_crystalline) + rate_count
    if one_water_crystalline:
        rate_count = a1*len(one_water_crystalline) + rate_count
    if two_water_crystalline:
        rate_count = a2*len(two_water_crystalline) + rate_count
    if three_water_crystalline:
        rate_count = a3*len(three_water_crystalline) + rate_count
    if four_water_crystalline:
        rate_count = a4*len(four_water_crystalline) + rate_count
        
    if rate_count == 0:
        pass
    else:
        crys_total = crys_rate*rate_count
        
    return crys_total 
     

#SIMULATION####################################################################

'TIME ZERO DATA'
'Initialising lists to be used'
#List of currently vulnerable pixels in the polymer grid
vulnerable_pixels = []
#List of currently vulnerable crystalline pixels
crystalline = []
#List of currently vulnerable amorphous pixels
amorphous = []
#Lists for currently vulnerable pixels with 1 adjacent acidic water unit 
one_water = []
one_water_amorphous = []
one_water_crystalline = []
#Lists for currently vulnerable pixels with 2 adjacent acidic water units 
two_water = []
two_water_amorphous = []
two_water_crystalline = []
#Lists for currently vulnerable pixels with 3 adjacent acidic water units 
three_water = []
three_water_amorphous = []
three_water_crystalline = []
#Lists for currently vulnerable pixels with 4 adjacent acidic water units 
four_water = []
four_water_amorphous = []
four_water_crystalline = []
#Special cases for first column of pixels
#Lists for currently vulnerable pixels right adjacent only to a pure pH 7 water unit
pure_water = []
pure_water_amorphous = []
pure_water_crystalline = []
#Lists for currently vulnerable pixels right adjacent to a pure water unit, and 1 acidic unit
pure_one_acidic = []
pure_one_acidic_amorphous = []
pure_one_acidic_crystalline = []
#Lists for currently vulnerable pixels right adjacent to a pure water unit, and 2 acidic units
pure_two_acidic = []
pure_two_acidic_amorphous = []
pure_two_acidic_crystalline = []
#Lists for currently vulnerable pixels right adjacent to a pure water unit, and 3 acidic units
pure_three_acidic = []
pure_three_acidic_amorphous = []
pure_three_acidic_crystalline = []
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
length = 500
thickness = length
print("Step", step)
print("Time", time)
print()
print("Generating grid...")
pixels = grid(length,thickness)
stationary = copy.copy(pixels)
time_series.append(stationary)
print()

#Initial condition rendering relevant pixels on first column as vulnerable to erosion
for i in range(1,length+1):
    vulnerable_pixels.append([i,1])
    pure_water.append([i,1])
    
vulnerable_classifier()

'Calculating rates for the first step'
#No self-inhibition model
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

# #Self-inhibition model
# #Dissolution rate for amorphous
# amorph_first_rate = amorph_rate*len(pure_water_amorphous)
# #Dissolution rate for crystalline
# crys_first_rate = crys_rate*len(pure_water_crystalline)
# #Find total rate of dissolution occurring
# total_first_rate = amorph_first_rate + crys_first_rate 
# #Probability of an amorphous pixel dissolving
# amorph_prob = amorph_first_rate/total_first_rate
# #Probability of a crystalline pixel dissolving
# crys_prob = crys_first_rate/total_first_rate
# #Generate time to next event via Exponential distribution
# #mu_timestep = -np.log(random.uniform(0,1))/total_rate
# mu_timestep = np.random.exponential(1/total_first_rate)

'Relative hydrolysis rates for different cases of vulnerable pixel exposure'
#One pH 7 water pixel only
p0 = 0.250
#One pH 7 water pixel and one acidic pH 4.8 pixel
p1 = 0.310
#One pH 7 water pixel and two acidic pH 4.8 pixels
p2 = 0.279
#One pH 7 water pixel and three acidic pH 4.8 pixels
p3 = 0.260
#One acidic pH 4.8 pixel
a1 = 0.060
#Two acidic pH 4.8 pixels
a2 = 0.029
#Three acidic pH 4.8 pixels
a3 = 0.010
#Four acidic pH 4.8 pixels
a4 = 0.003

'Gillespie Algorithm'
total_time = 50
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
    clear_lists()
    
    #Identifying and classifying the updated vulnerable pixels
    is_vulnerable(length,thickness)
    vulnerable_classifier()
    
    #Updating the new rates
    #No self-inhibition model
    amorph_total_rate = amorph_rate*len(amorphous)
    crys_total_rate = crys_rate*len(crystalline)
    
    # #Self-inhibition model
    # amorph_total_rate = amorph_total_rates()
    # crys_total_rate = crys_total_rates()
    
    total_rate = amorph_total_rate + crys_total_rate 
    if total_rate == 0:
        break
    else:
        amorph_prob = amorph_total_rate/total_rate
        crys_prob = crys_total_rate/total_rate
    #mu_timestep = -np.log(random.uniform(0,1))/total_rate
    mu_timestep = np.random.exponential(1/total_rate)

#PLOTTING######################################################################

erosion = cm.get_cmap('jet', 10) #Adjusts colorbar
plt.figure(dpi=125)
plt.pcolormesh(pixels,cmap = erosion)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel('Thickness (\u03BCm)',fontsize=13)
plt.ylabel('Length (\u03BCm)',fontsize=13)
plt.title('Polymer matrix after Time = {}'.format(total_time))
plt.show()