from cmath import inf
import scipy
from unicodedata import name
import numpy as np
from ete3 import Tree


def is_symmetric(A, tol=1e-8):
    # return scipy.sparse.linalg.norm(A-A.T, scipy.Inf) < tol;
    
    return scipy.linalg.norm(A-A.T, scipy.Inf) < tol;

def triangularInequality(A):
    n= np.shape(A)[0]
    flag= True
    for i in range(n):
        if(not flag):
             break
        for j in range(n):
            if(not flag):
                 break
            for k in range(n):
                if(A[i][j] + A[j][k] < A[k][i]):
                    flag= False
                    break
    return flag            

def BasicNJ(D,labels):
    #checking whether it is a valid distance matrix
    # if(not is_symmetric(A) or not triangularInequality(A)):
    if(not is_symmetric(D)):
        return -1
    D0 = np.copy(D) #Creating a copy of the original mat so as to not lose it
    t= Tree()
    labels1= labels.copy() 
    for x in labels1:
        temp=t.add_child(name=x)
    print(t) #initial graph

    its=0 #internal node counter
    while(np.shape(D0)[0] > 2):
        its+=1
        r= np.shape(D0)[0]
        
        # we don't really need to maintain a matrix Q(i,j)
        
        imin = -inf #stores the i index of current min element and jmin follows similarly
        jmin = -inf
        dmin = inf # current min value
        distSum= np.sum(D0, axis=1) #stores summmation(d(i,k))  {for each i, k ranges from 0 to r-1}

        for i in range(r):
            for j in range(i+1,r):
                temp= (r-2)*D0[i][j] - distSum[i]- distSum[j]
                if(temp < dmin):
                    dmin=temp
                    imin=i
                    jmin=j # will always be jmin > imin
        # suppose f,g are the selected branches and their parent that we created is u
        # finding len of fu and gu
        dfu= (dmin+ 2*distSum[imin])/(2*(r-2))
        dgu= (dmin+ 2*distSum[jmin])/(2*(r-2))
        #note to self: use vectorisation to add and remove the f and g column and directly add the u column
        u= np.add(D0[imin],D0[jmin])/2- ((dfu+dgu)/2) #u th row
        #adding the uth row and column
        D0= np.r_[D0,[u]]
        # print(D0)

        u=np.append(u,0)
        # print(u)
        D0= np.c_[D0,u]
        F= t&labels1[imin]
        F_detached= F.detach()
        G= t&labels1[jmin]
        G_detached= G.detach()
        newname="("+labels1[imin]+","+labels1[jmin]+")IntNode"+str(its)
        labels1.append(newname)  #creating label for the new node
        F_detached.dist = dfu
        G_detached.dist = dgu
        tmep= t.add_child(name=newname)
        
        tmep.add_child(F_detached)
        tmep.add_child(G_detached)
        labels1.remove(labels1[imin])
        labels1.remove(labels1[jmin-1])     #removing one node after one has been removed hence the -1
        
        #deleting the f g rows and columns from the distance matrix
        D0= np.delete(D0,imin,0)
        D0= np.delete(D0,jmin-1,0)
        D0= np.delete(D0,imin,1)
        D0= np.delete(D0,jmin-1,1)
        # print("Updated Distance Matrix: ")
        # print(D0)
        # print("Updated Labels: ")
        # print(labels1)
        # print(t.get_ascii(attributes=["name","dist"],show_internal=True))
    return t
