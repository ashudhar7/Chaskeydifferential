#!/usr/bin/env python
# coding: utf-8

# In[9]:


def highwaylist():
    alph,bet,gam,pro = [],[],[],[]
    with open('alph.txt') as ip1, open('bet.txt') as ip2, open('gam.txt') as ip3, open('pro.txt') as ip4:
        for i1,i2,i3,i4 in zip(ip1,ip2,ip3,ip4):
            hw.append(NMTS(int(i1.strip()),int(i2.strip()),int(i3.strip()),float(i4.strip()))) 


# In[73]:

def mask_1(a,b,wshift):
    s1=0
    var=0
    while(s1==0):
        a1=(a<<(wshift))&mask
        b1=(b<<(wshift))&mask
        s1=a1+b1
        if(s1==0):
            var=var+1
        wshift=wshift-1
        mask1=(mask<<var)&mask
    return mask1


# In[74]:

import math
def weight1(alpha1,beta1,gamma1): 
    x1= (alpha1 << 1)& mask
    y1= (beta1  << 1)& mask
    z1= (gamma1  << 1)& mask
    not_x=(x1^mask)
    eq=((not_x)^y1)&(not_x^z1)
    s=eq&(alpha1^beta1^gamma1^(y1))
    if(s==0):        
        h=bin((~eq)& mask)
        wt=(h[1:].count("1"))
        return wt,math.pow(2,-wt)
    else:
        return -200,-200


# In[75]:

#weight1(0,2147483648,2147483648)


# In[76]:

#right circular shifts
def ROR(x,r):
    x = ((x << (Chaskey_TYPE - r)) + (x >> r)) & mask
    return x

#left circular shifts
def ROL(x,r):
    x = ((x >> (Chaskey_TYPE - r)) + (x << r)) & mask
    return x


# In[77]:

def PDDT(n,p_thres,k,wt,da,db,dc,ddt,m):
    if(n==k):
        #print(da,db,dc,wt)
        ddt.append(NMTS(da,db,dc,wt))
        return 
    
    for z in range(0,2):
        new_da= (da & ((2 ** m) - 1))
        new_db= (db & ((2 ** m) - 1))
        new_dc= (dc | (z<<k) & mask)
        pw,wt= weight1(new_da,new_db,new_dc)
        if(wt >= p_thres or ((z==1) and (k+1)!=n)):                    
            PDDT(n,p_thres,k+1,wt,da,db,new_dc,ddt,m+1)  
            
    return ddt


# In[78]:

def find_diff_path(st3,st2,st1,st0,Chaskey_ROUNDS,n):
    lp=1   
    pthres=0.1
    temp_wt=0
    tempdec_list=[]
    #print(hex(st1),hex(st0),'findpath')
    for i in range (0,Chaskey_ROUNDS): 
        #condition to check random output decision with valid differential
        inp3=st3
        inp2=st2
        inp1=st1
        inp0=st0
        ddt=[]
        ddt=PDDT(Chaskey_TYPE,pthres,0,1,st3,st2,0,ddt,1)
        if(len(ddt)!=0):
            index=(random.randint(0,len(ddt)-1))
            op1=ddt[index].dz
            w,w1=weight1(st3,st2,op1)
            wt1=w
        else:
            op1=st3^st2
            w,w1=weight1(st3,st2,op1)
            wt1=w
        #print(st3,st2,op1,wt1)
        out1=op1
        
        ddt=[]
        ddt=PDDT(Chaskey_TYPE,pthres,0,1,st1,st0,0,ddt,1)
        if(len(ddt)!=0):
            index=(random.randint(0,len(ddt)-1))
            op2=ddt[index].dz
            w,w2=weight1(st1,st0,op2)
            wt2=w
        else:
            op2=st1^st0
            w,w2=weight1(st1,st0,op2)
            wt2=w
        #print(st2,st1,op2,wt2)
        
        out2=op2
        
        st3=ROL(st3,alpha)
        st3=st3^op1
        st0=ROL(st0,gamma)
        st0=st0^op2
        op1=ROL(op1,theta)
        
        ddt=[]
        ddt=PDDT(Chaskey_TYPE,pthres,0,1,st3,op2,0,ddt,1)
        if(len(ddt)!=0):
            index=(random.randint(0,len(ddt)-1))
            op3=ddt[index].dz
            w,w3=weight1(st3,op2,op3)
            wt3=w
        else:
            op3=st3^op2
            w,w3=weight1(st3,op2,op3)
            wt3=w
            
        out3=op3    
        #print(st1,st0,op3,wt3)
        st3=ROL(st3,beta)
        st3=st3^op3
        op3=ROL(op3,theta)
        
        ddt=[]
        ddt=PDDT(Chaskey_TYPE,pthres,0,1,op1,st0,0,ddt,1)
        if(len(ddt)!=0):
            index=(random.randint(0,len(ddt)-1))
            op4=ddt[index].dz
            w,w4=weight1(op1,st0,op4)
            wt4=w
        else:
            op4=op1^st0
            w,w4=weight1(op1,st0,op4)
            wt4=w
            
        out4=op4
        #print(st1,st0,op3,wt3)
        st0=ROL(st0,delta)
        st0=st0^op4
        st1=op3
        st2=op4
        #print(st1,st0,op3,wt1,wt2,wt3,wt4)
        wt=wt1+wt2+wt3+wt4
        #print(wt)
        tempdec_list.append(NMTS1(inp3,inp2,inp1,inp0,out1,out2,out3,out4,wt1,wt2,wt3,wt4,wt)) 
        temp_wt= temp_wt+tempdec_list[i].wt
        #update state with new valid output differential
    n=n+1
    #print(wt1,wt2,wt3,w4)
    return tempdec_list, temp_wt, n


# In[79]:

def find_best_path(st3,st2,st1,st0, Chaskey_ROUNDS, wt_above, best_wt):
    #print(hex(st1),hex(st0),'best_path')
    #using n as index value for list
    n=0
    for r in range(Chaskey_ROUNDS,0,-1):
        tempdec_list, temp_wt, n = find_diff_path(st3,st2,st1,st0,r, n)
        if((temp_wt+wt_above) < best_wt):
            best_wt= temp_wt+wt_above
            for i,j in zip(range(n-1,22), range(len(tempdec_list))):
                dec_list[i]=(NMTS1(tempdec_list[j].da,tempdec_list[j].db,tempdec_list[j].dc,tempdec_list[j].dd,tempdec_list[j].op1,tempdec_list[j].op2,tempdec_list[j].op3, tempdec_list[j].op4, tempdec_list[j].wt1,tempdec_list[j].wt2,tempdec_list[j].wt3, tempdec_list[j].wt4, tempdec_list[j].wt))
        #print(n)
        if(n<Chaskey_ROUNDS):
            st3=dec_list[n].da
            st2=dec_list[n].db
            st1=dec_list[n].dc
            st0=dec_list[n].dd
        wt_above= wt_above+dec_list[n-1].wt
    return best_wt


# In[ ]:

#alpha beta are for left and right circular shift  
import random
import math 
alpha, beta, gamma, delta, theta = 5,7,8,13,16
Chaskey_ROUNDS=8
Chaskey_TYPE=32
mask = 2 ** 32 - 1
wshift=31
dec_list=[0]*22
wt_above=0
best_wt=9999
s=99999
class NMTS1(object):
    """__init__() functions as the class constructor"""
    def __init__(self, da=None,db=None, dc=None, dd=None, op1=None, op2=None, op3=None, op4=None, wt1=None, wt2=None, wt3=None, wt4=None, wt=None):
        self.da = da
        self.db = db
        self.dc = dc
        self.dd = dd
        self.op1 = op1
        self.op2 = op2
        self.op3 = op3
        self.op4 = op4
        self.wt1 = wt1
        self.wt2 = wt2
        self.wt3 = wt3
        self.wt4 = wt4
        self.wt = wt

class NMTS(object):
    """__init__() functions as the class constructor"""
    def __init__(self, dx=None, dy=None, dz=None, wt=None):
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.wt = wt
#highwaylist()
hw=[]
#nested loop to run number of time
highwaylist()
#ch=0 
while(best_wt>128):
    ch=(random.randint(0,len(hw)))
    ch1=(random.randint(0,len(hw)))
    #best_wt=find_best_path(0x1AC8DA46, 0x73C0D20A, 0x9282B2A3, 0x02947AA1,Chaskey_ROUNDS,wt_above,best_wt)
    best_wt=find_best_path(hw[ch].dx,hw[ch].dy,hw[(ch+ch1)%len(hw)].dx,hw[(ch+ch1)%len(hw)].dy,Chaskey_ROUNDS,wt_above,best_wt)
    #ch=ch+1
    if(best_wt<s):
        s=best_wt
        print(s,ch) 

#print(hex(st1),hex(st0),0)
print("Dec list is: st3,st2,st1,st0,op1,op2,op3,op4,wt1,wt2,wt3,wt4,wt")   
for i in range(0,Chaskey_ROUNDS): 
    print("Starting input of %i round and weight:" %i,hex(dec_list[i].da), hex(dec_list[i].db),hex(dec_list[i].dc),hex(dec_list[i].dd),hex(dec_list[i].op1),hex(dec_list[i].op2),hex(dec_list[i].op3), hex(dec_list[i].op4),(dec_list[i].wt1),(dec_list[i].wt2),(dec_list[i].wt3),(dec_list[i].wt4),(dec_list[i].wt))
    #print("Starting input of %i round and weight:" %i,hex(dec_list[i].da), hex(dec_list[i].db),hex(dec_list[i].dc),hex(dec_list[i].dd),(dec_list[i].wt))


print("Best weight is:",best_wt) 




# In[ ]:




