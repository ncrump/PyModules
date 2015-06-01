"""
Created on Fri Nov 1 20:15:48 2013
"""

#  Cluster Find
"""
The following routines identify clusters on a 2D square lattice
and check for a spanning cluster.  

1. Cluster Find (slow version for plotting small lattices ~ 50x50)
2. Cluster Find (faster version for finding spanning clusters)  
"""

import numpy as np
import matplotlib.pyplot as plt


# 1
# Cluster Find (slow version for plotting)
#*******************************************************************
#-------------------------------------------------------------------
# Ox = array containing x positions of occupied sites
# Oy = array containing y positions of occupied sites
# N = length of square lattice (ex: N=40 for 40x40 grid)
# iters = number of correction scans for labeling (ex: iters = 2)
# txt = turn cluster numbers on/off (entered as: 'on' or 'off')
# vis = turn plotting on/off (only use for small clusters of N<50)
#------------------------------------------------------------------- 
def ClusterFindSlow(Ox, Oy, N, iters, txt='off', vis='off'):
    
    Nsites = len(Ox)    # number of occupied sites
    lattice = range(N)  # length of square lattice
    
    # generate [x,y] position arrays of occupied and unoccupied sites
    Osites = [[Ox[i],Oy[i]] for i in range(Nsites)]
    Psites = [[i,j] for i in lattice for j in lattice]
    
    # define array for storing labels, initialized to integer zeros
    labels = np.zeros(Nsites) 
    labels = labels.astype(int)
    label = 1        
    
    # enter main loop to find clusters and assign cluster labels
    # ----------------------------------------------------------
    for i in range(Nsites):
        # get occupied site
        site = Osites[i]
        
        # if site has no label, assign new integer label
        if labels[i] == 0: 
            labels[i] = label
            label = label + 1  # increment label
            iLabl = label - 1  # unincremented label for checks below
        
        # get neighboring sites left, right, up, down
        lt = [site[0]-1, site[1]]
        rt = [site[0]+1, site[1]]
        up = [site[0], site[1]+1]
        dn = [site[0], site[1]-1]
        pos = [lt,rt,up,dn]
        
        # store site index and label to temp arrays for checks below       
        posIndx = [i]
        posLabl = [iLabl]
        
        # loop through neighbor sites to see if cluster sites
        for j in pos:
            if j in Osites:
                # get neighbor site index and label
                indx = Osites.index(j)
                jLabl = labels[indx]
                
                # if cluster site already labeled, store to array
                if jLabl != 0:
                    posIndx.append(indx)
                    posLabl.append(jLabl)
                
                # if cluster site not labeled, give it new label
                if jLabl == 0:
                    labels[indx] = labels[i]
        
        # get lowest label of neighbor cluster sites
        # relabel neighbor cluster sites to lowest label
        if len(posLabl) > 0:
            minLabl = min(posLabl)
            labels[(posIndx)] = minLabl
    # ----------------------------------------------------------
            
                    
    # enter correction loop to refine cluster labels
    # ----------------------------------------------------------
    # loop through lattice as many times as user input
    for chk in range(iters+1):
        for site in Osites:
            # get neighbor sites again
            lt = [site[0]-1, site[1]]
            rt = [site[0]+1, site[1]]
            up = [site[0], site[1]+1]
            dn = [site[0], site[1]-1]
            pos = [lt,rt,up,dn]
            
            # get lowest label of neighbor cluster sites
            chkIndx = [Osites.index(z) for z in pos if z in Osites]
            chkLabl = [labels[z] for z in chkIndx]
            
            # relabel neighbor cluster sites to lowest label 
            # this refines cluster identification
            if len(chkLabl) > 0:
                minLabl = min(chkLabl)
                labels[(chkIndx)] = minLabl            
    # ----------------------------------------------------------
                

    # scale labels to lowest sequential numbering of clusters
    Lsort = list(set(labels))
    mx = len(Lsort) + 1      
    for i in range(1,mx):
        indx1 = Lsort[i-1]
        indx2 = np.where(labels == indx1)[0]
        labels[indx2] = i
    
    # store x,y locations of unoccupied sites
    Px = [i[0] for i in Psites]
    Py = [i[1] for i in Psites]
        
    # enter color and plotting routine
    # ----------------------------------------------------------
    if vis == 'on':
        # define RGB color array for plotting clusters
        nColors = max(labels) + 1           
        R = np.random.uniform(0.1,1,nColors)  # random R value
        G = np.random.uniform(0.1,1,nColors)  # random G value
        B = np.random.uniform(0.1,1,nColors)  # random B value 
        colors = [[R[i],G[i],B[i]] for i in range(nColors)]
        
        # set marker size from derived visual scaling
        msize = int(128325*(N**-2.179))
        
        # loop to plot cluster sites with custom colors
        for i in range(Nsites):
            siteColor = colors[labels[i]]
            plt.scatter(Ox[i],Oy[i],s=msize,marker='s',c=siteColor,edgecolor='k')
            # plot cluster numbers if user input = 'on'
            if txt == 'on':
                plt.text(Ox[i],Oy[i],str(labels[i]),horizontalalignment='center',verticalalignment='center')
        
        # plot remaining unoccupied sites of lattice
        plt.scatter(Px,Py,s=msize,marker='s',c='none',edgecolor='k')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.xlim(-1.0,N)
        plt.ylim(-1.0,N)
    # ----------------------------------------------------------
        
    # print total clusters found
    clusters = len(set(labels))
    print '\n',clusters,' clusters found'
#*******************************************************************
    
    
# 2
# Cluster Find (faster version for finding spanning cluster)
#*******************************************************************
#-------------------------------------------------------------------
# Ox = array containing x positions of occupied sites
# Oy = array containing y positions of occupied sites
# N = length of square lattice (ex: N=40 for 40x40 grid)
# iters = number of correction scans for labeling (ex: iters = 2)
# msg = set to display message or not (turns off print statements)
#------------------------------------------------------------------- 
def ClusterFindFast(Ox, Oy, N, iters, msg='off'):
    
    # get number of occupied sites
    Nsites = len(Ox)    
    
    # initialize matrix of labels
    labels = np.zeros((N,N))
    labels = labels.astype(int)
    
    # turn [x,y] coords of occupied sites into matrix form
    for i in range(Nsites):
        labels[Oy[i]][Ox[i]] = 1 
    
    # reorder matrix to match x,y representation
    labels = labels[::-1]
    indx = np.where(labels == 1)
    label = 1

    # arrays to store equivalent labels    
    equiv0 = np.array([])
    equiv1 = np.array([])
    
    # find and assign initial cluster labels
    # ----------------------------------------------------------
    for i in range(Nsites):
        # get occupied site
        row = indx[0][i]
        col = indx[1][i]
        
        # get label of site to left
        if col != 0:
            left = labels[row][col-1]
        else: left = 0
        
        # get label of site to right
        if row != 0: 
            up = labels[row-1][col]
        else: up = 0
                
        # check and assign initial labels
        # keep track of equivalent labels
        if left != 0 and left < up:
            labels[row][col] = left
            labels[row-1][col] = left
            equiv0 = np.append(equiv0,up)
            equiv1 = np.append(equiv1,left)
            
        elif up != 0 and up < left:
            labels[row][col] = up
            labels[row][col-1] = up
            equiv0 = np.append(equiv0,left)
            equiv1 = np.append(equiv1,up)
            
        elif left != 0 and left == up:
            labels[row][col] = up
            
        elif left != 0:
            labels[row][col] = left
            
        elif up != 0:
            labels[row][col] = up
            
        elif left == 0 and up == 0:
            labels[row][col] = label
            label = label + 1 
    # ----------------------------------------------------------
        
    # find equivalent cluster labels
    # ----------------------------------------------------------
    equivs = []
    eqvLbl = []
    
    while len(equiv0) >= 1:
        tmp = [equiv0[0],equiv1[0]]
        ndx = [0]
        
        for i in range(iters):
            Neqv = len(equiv1)
            
            for j in range(Neqv):
                if equiv0[j] in tmp or equiv1[j] in tmp:
                    tmp.append(equiv0[j])
                    tmp.append(equiv1[j])
                    ndx.append(j)
                    
            equiv0 = np.delete(equiv0,ndx)
            equiv1 = np.delete(equiv1,ndx)
            tmp = list(set(tmp))
            ndx = []
       
        Lmn = tmp.index(min(tmp))
        eqvLbl.append(tmp.pop(Lmn))
        equivs.append(tmp)
    # ---------------------------------------------------------- 
        
    # scan back through lattice to update final cluster labels
    # ----------------------------------------------------------      
    pos = 0    
    for i in equivs:
        for j in i:
            Lndx = np.where(labels == j)
            labels[Lndx[0],Lndx[1]] = eqvLbl[pos]
        pos = pos + 1      
    # ----------------------------------------------------------      
        
    # check for spanning cluster
    # ----------------------------------------------------------
    top = set(labels[0][::])
    bot = set(labels[N-1][::])
    lft = set(labels[::][0])
    rgt = set(labels[::][N-1])
    
    span = list(set.intersection(top,bot,lft,rgt))
    
    if span[0] == 0:
        span.remove(0)
    
    if len(span) > 0:
        spanStr = 'spanning cluster found'
        spanNum = len(np.where(labels == span[0])[0])
        spanLab = span[0]
    else: 
        spanStr = 'no spanning cluster found'
        spanNum = 0
        spanLab = 0
        
    clusters = len(np.unique(labels)) - 1
    # ----------------------------------------------------------
    
    if msg == 'on':
        # output and return clusters found
        print '\n',clusters,' clusters found'
        print spanStr

    # returns total clusters, sites in spanning cluster, labels, spanning label
    return clusters,spanNum,labels,spanLab
#*******************************************************************