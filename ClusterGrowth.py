"""
Created on Sun Oct 27 17:50:45 2013
"""

# Cluster Growth Models
"""
The following cluster growth models on a 2D square lattice are below. 
1. Eden Cluster
2. Epidemic Cluster
3. Diffusion-Limited Aggregation Cluster (DLA)
4. Percolation Cluster
"""

import random
import numpy as np


# 1
# Eden Cluster Model
#*******************************************************************
#-------------------------------------------------------------------
# N = number of sites to occupy
# seed = start location of occupied site (integers only - ex: [0,0])
#------------------------------------------------------------------- 
def Eden(N, seed):

    # define initial perimeter sites
    p1 = [seed[0]-1,seed[1]]    
    p2 = [seed[0]+1,seed[1]]     
    p3 = [seed[0],seed[1]-1] 
    p4 = [seed[0],seed[1]+1]
    
    # store occupied sites and perimeter sites
    Osites = [seed]
    Psites = [p1,p2,p3,p4]
    
    # start while loop counter
    n = 1    
    
    # loop until N sites occupied
    while n < N:
        # randomly choose from perimeter sites to occupy
        rand = len(Psites)
        indx = int(np.random.randint(0,rand,1)) 
        # get and store new occupied site
        site = Psites.pop(indx)
        Osites.append(site)
        
        # get new perimeter sites from new occupied site
        p1 = [site[0]-1,site[1]]    
        p2 = [site[0]+1,site[1]]     
        p3 = [site[0],site[1]-1] 
        p4 = [site[0],site[1]+1]
        perims = [p1,p2,p3,p4]

        # check that new perimeter sites not already used
        for i in perims:
            if i not in Psites and i not in Osites:
                Psites.append(i)
                
        # increment counter
        n = n+1
        
    # store and return x,y locations of occupied sites
    Ox = [i[0] for i in Osites]
    Oy = [i[1] for i in Osites]    
    
    # returns occupied sites of the cluster
    return Ox,Oy
#*******************************************************************
    
    
# Epidemic Cluster Model
#*******************************************************************
#-------------------------------------------------------------------
# N = number of sites to occupy
# seed = start location of occupied site (integers only - ex: [0,0])
# p = probability to occupy neighbor site (kill site with 1-p) 
#------------------------------------------------------------------- 
def Epidemic(N, seed, p):

    # define initial perimeter sites
    p1 = [seed[0]-1,seed[1]]    
    p2 = [seed[0]+1,seed[1]]     
    p3 = [seed[0],seed[1]-1] 
    p4 = [seed[0],seed[1]+1]
    
    # store occupied sites and perimeter sites
    Osites = [seed]
    Psites = [p1,p2,p3,p4]
    
    # start while loop counter
    n = 1    
    
    # loop until N sites occupied
    while n < N:
        # randomly choose from perimeter sites to occupy or remove
        rand = len(Psites)
        indx = int(np.random.randint(0,rand,1)) 
        
        
        # get random number to determine if site will be occupied or killed
        # based on input probability p
        randGen = np.random.uniform(0,1,1)[0]
        
        # get and store new occupied site based on p
        if randGen >= 0 and randGen <= p:
            site = Psites.pop(indx)
            Osites.append(site)
            
            # get new perimeter sites from new occupied site
            p1 = [site[0]-1,site[1]]    
            p2 = [site[0]+1,site[1]]     
            p3 = [site[0],site[1]-1] 
            p4 = [site[0],site[1]+1]
            perims = [p1,p2,p3,p4]
    
            # check that new perimeter sites not already used
            for i in perims:
                if i not in Psites and i not in Osites:
                    Psites.append(i)
                    
            # increment counter
            n = n+1
        
        # remove site based on p
        else:
            Psites.remove(Psites[indx])
        
    # store and return x,y locations of occupied sites
    Ox = [i[0] for i in Osites]
    Oy = [i[1] for i in Osites]    
    
    # returns occupied sites of the cluster
    return Ox,Oy
#*******************************************************************
    
    
# 3
# Diffusion-Limited Aggregation Cluster Model (DLA)
#*******************************************************************
#-------------------------------------------------------------------
# N = number of sites to occupy
# seed = start location of occupied site (integers only - ex: [0,0])
# R = radius of circle from outtermost perimeter site to start walks
#------------------------------------------------------------------- 
def DLA(N, seed, R):

    # define initial perimeter sites
    p1 = [seed[0]-1,seed[1]]    
    p2 = [seed[0]+1,seed[1]]     
    p3 = [seed[0],seed[1]-1] 
    p4 = [seed[0],seed[1]+1]
    
    # store occupied sites and perimeter sites
    Osites = [seed]
    Psites = [p1,p2,p3,p4]
    perim0 = [p1,p2,p3,p4]
    
    # start while loop counter
    n = 1
    
    # loop until N sites occupied
    while n < N:
        
        # get max distance of perimeter sites for updating radius
        Pdist = [(i[0]**2 + i[1]**2)**0.5 for i in Psites]
        Pmax = int(max(np.round(Pdist)))
        
        # update radius and generate random angle along circle
        # this becomes start point for random walk from outter radius into cluster
        r = Pmax + R
        theta = np.random.uniform(0,2*np.pi,1)
        
        # convert randomly generated point along circle to cartesian coords
        # this is start point for random walk from outter radius into cluster
        x0 = int(r*np.cos(theta))
        y0 = int(r*np.sin(theta))
        xyRand = [x0,y0]
        
        # set distance for loop below
        dist = r
        
        
        # step through random walk as long as particle is inside circle
        while dist <= r:
            
            # perform random walk
            # generate random uniform number between [0,1)
            rand = np.random.uniform(0,1,1)
            
            # depending on random number go right, left, up, or down by 1 unit
            if 0.00 <= rand < 0.25: xyRand[0] = xyRand[0]+1
            if 0.25 <= rand < 0.50: xyRand[0] = xyRand[0]-1
            if 0.50 <= rand < 0.75: xyRand[1] = xyRand[1]+1
            if 0.75 <= rand < 1.00: xyRand[1] = xyRand[1]-1
            
            # set distance to check if walk goes outside radius
            dist = int(round((xyRand[0]**2 + xyRand[1]**2)**0.5))
            
            # if random walk hits a perimeter site it becomes an occupied site
            if xyRand in Psites:
                Osites.append(xyRand)
                Psites.remove(xyRand)
                
                # get new perimeter sites from new occupied site
                p1 = [xyRand[0]-1,xyRand[1]]    
                p2 = [xyRand[0]+1,xyRand[1]]     
                p3 = [xyRand[0],xyRand[1]-1] 
                p4 = [xyRand[0],xyRand[1]+1]
                perims = [p1,p2,p3,p4]
                
                # check that new perimeter sites not already used
                for i in perims:
                    if i not in Psites and i not in Osites:
                        Psites.append(i)
                        perim0.append(i)
                break
            

        # increment counter only if site was occupied
        n = len(Osites)
        
    # store and return x,y locations of occupied sites
    Ox = [i[0] for i in Osites]
    Oy = [i[1] for i in Osites] 
    
    # store and return x,y locations of perimeter sites
    Px = [i[0] for i in perim0]
    Py = [i[1] for i in perim0]
    
    # returns occupied sites of the cluster and unoccupied lattice sites
    return Ox,Oy,Px,Py
#*******************************************************************
    

# 4
# Percolation Cluster Model
#*******************************************************************
#-------------------------------------------------------------------
# N = square lattice size (N x N)
# p = percolation probability (number of sites to occupy, ex: p=0.4)
#------------------------------------------------------------------- 
def Percolation(N, p):
    
    # define number of lattice sites and number of sites to occupy
    Nsites = N*N
    Nocc = int(Nsites*p)
    
    # loop to write all coordinates of lattice to array    
    lattice = [[i,j] for i in range(N) for j in range(N)]
    lattice = np.array(lattice)
    
    # get random sites to occupy without duplicates
    Oindx = random.sample(range(Nsites),Nocc)   
    # store occupied and un-occupied sites
    Osites = lattice[Oindx]
    Psites = np.delete(lattice,Oindx,0)
        
    # store and return x,y locations of occupied sites
    Ox = [i[0] for i in Osites]
    Oy = [i[1] for i in Osites] 
    
    # store and return x,y locations of un-occupied sites
    Px = [i[0] for i in Psites]
    Py = [i[1] for i in Psites]
    
    return Ox,Oy,Px,Py
#*******************************************************************