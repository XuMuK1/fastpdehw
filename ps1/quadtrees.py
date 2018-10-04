import numpy as np


#Tree Structure
#class definition

class TreeNode:
    
    def __init__(self,boxsize=0,mass=1,coord=np.zeros([2]),level=0):
    #boxsize is a diameter of the box in 2norm
    #mass and coord -- obvious
    #level is tree depth, level=0 corresponds to root

        self.mass=mass
        self.coord=coord
        self.children=[]
        self.level=level
        self.boxsize=boxsize
        
        
    def addChild(self,child):
        self.children.append(child)
        
    def recomputeValues(self):
        if(len(self.children)>0):
            [self.children[k].recomputeValues for k in np.arange(len(self.children))]#faster
            self.mass = 0
            self.coord = -1
            for child in self.children:
                    try:
                        if self.coord.shape:#if it is not -1
                            self.mass = self.mass + child.mass #really????
                            self.coord = self.coord+ child.coord*child.mass
                    except:
                        self.mass = child.mass
                        self.coord = child.coord * child.mass
                        
            
            self.coord = self.coord / self.mass
        else:
            return self.coord, self.mass
        



#CONSTRUCTION

# Tree Construction

def AssignParticlesToClusters(locs,seps):
    #locs: d X N locations
    #seps: N separators
    clIds = np.zeros([locs.shape[1]])
    for locId in np.arange(locs.shape[1]):
        ind=(locs[:,locId]-seps >=0).astype('int32')
        ind=str(ind.tolist())[1:-1].replace(' ','')
        ind=ind.replace(',','')
        clIds[locId]=int( ind, base=2)
        
    return clIds
    
def ConstructTree(xs,ms,root,level):
    #xs is array d X N, N is number of particles, d is dimension
    
    
    if(len(xs.shape)==1 or xs.shape[1]==1):
        #1 particle, leave it
        root.mass = ms
        root.coord = xs
        root.level = level
        
        return root
    
    #Otherwise go recursively
    
    xMin = np.amin(xs,axis=1)
    xMax = np.amax(xs,axis=1)
    separators = (xMin + xMax)/2 # QuadTree
    clusters = AssignParticlesToClusters(xs,separators)
    
    clIds = np.unique(clusters)
    
    for clId in clIds:
        clKeys = (clusters==clId)
        child = TreeNode(boxsize=np.linalg.norm(xMax-separators,2),level=level+1)
        child=ConstructTree(xs[:,clKeys],ms[clKeys],child,level+1)
        child.recomputeValues()
        root.addChild(child)
    
    return root

def DeepTraverse(root): #depth traversal
    
    if(len(root.children)==0):
        return
    
    for c in root.children:
        print(c.level,c.mass,c.coord)
        DeepTraverse(c)

## working EXAMPLE!!

#xs = np.array([[0,1,0,1,1.2,1.25],[0,0.75,1,1,1,0.9]])
#ms= np.array([1,2,3,4,5,6])

#treeRoot = TreeNode()
#treeRoot=ConstructTree(xs,ms,treeRoot,level=0) 
