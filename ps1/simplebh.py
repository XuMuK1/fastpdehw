import numpy as np
import quadtrees as qtrees

#BARNES HUT
#idea is simple: given particle system with interations,
#sum up the interactions on one particular particle

def SumTraverse(root,loc,m,interFun,farParameter=2):
    #root is QTree root
    #loc is location of particle
    #m is its mass
    #interFun is interactionFunction(m,m0,x,x0) between two particles
    #farParameter=2 means to approximate if dist/boxsize > 2
    #find all interactions
    
    if(len(root.children)>0):
        inter=0
        for child in root.children:
            
            if( np.linalg.norm(child.coord-loc,2)/child.boxsize >farParameter ):
                #if it can be approximated,then
                #print('Approximate!!')
                inter=inter+interFun(m,child.mass,loc,child.coord)
            else:
                #if cannot be approximated
                inter=inter+SumTraverse(child,loc,m,interFun,farParameter)
                
        return inter
    else:
        if(np.linalg.norm(loc-root.coord,2)>1e-10):#if this is not the same particle
            return interFun(m,root.mass,loc,root.coord)
        else:
            return 0
        
        
def BHSummation(xs,ms,interFun,farParameter=2):
    # farParameter = 2  if distance is 2 times larger than box size, then approximate

    treeRoot = qtrees.TreeNode()
    treeRoot = qtrees.ConstructTree(xs,ms,treeRoot,level=0)
    
    inters = np.zeros([xs.shape[1]])
    for i in np.arange(xs.shape[1]):
        inters[i] = SumTraverse(treeRoot,xs[:,i],ms[i],interFun,farParameter)
    return inters 
