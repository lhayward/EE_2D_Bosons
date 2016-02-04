import numpy as np
import time
import matplotlib.pyplot as plt

########################################  decimalStr  #########################################
# Converts a number to a string and replaces all commas with hyphens
# (used when naming files, where we don't want periods in the file name)
###############################################################################################
def decimalStr(num):
  res = str(num)
  length = len(res)
  index = res.find('.')
  if index >= 0:
    res = res[0:index] + '-' + res[(index+1):length]
  return res
# end decimalStr(num) function

##################################  getPhiPhiCorrelator_2d  ###################################
# Calculates the phi-phi correlator for a 1D free bosonic system
# Input paramters: 
#    L: length of lattice (integer)
#    bc: boundary conditions ('PBC' or 'APBC')
#    mass: boson mass (float)
#    r, rprime: (x,y) coordinates on lattice of the two phi variables (floats)
###############################################################################################
def getPhiPhiCorrelator_2d(L, bc, mass, r, rprime):
  corr = 0
  
  (x,y) = r
  (xprime, yprime ) = rprime
  
  def omega_2d(kx, ky):
    return np.sqrt( 4*np.sin(kx/2.0)**2 + 4*np.sin(ky/2.0)**2 + mass**2 )
  
  if bc == 'PBC':
    kx = (2*np.array(range(0,L)))*np.pi/L
    ky = (2*np.array(range(0,L)))*np.pi/L
  else: #APBC
    kx = ( 2*np.array(range(0,L)) + 1)*np.pi/L
    ky = ( 2*np.array(range(0,L)) + 1)*np.pi/L
  
  corr = sum( [ sum( np.cos(kx*(x-xprime))*np.cos(kyy*(y-yprime))/omega_2d(kx,kyy) ) for kyy in ky] )
  #corr = sum( [np.cos(kxx*(x-xprime))*np.cos(kyy*(y-yprime))/omega_2d(kxx,kyy) for kxx in kx for kyy in ky] )

  return corr/(2*L**2)

###################################  getPiPiCorrelator_1d  ####################################
# Calculates the pi-pi correlator for a 1D free bosonic system
# Input paramters: 
#    L: length of lattice (integer)
#    bc: boundary conditions ('PBC' or 'APBC')
#    mass: boson mass (float)
#    r, rprime: (x,y) coordinates on lattice of the two phi variables (floats)
###############################################################################################
def getPiPiCorrelator_2d(L, bc, mass, r, rprime):
  corr = 0
  
  (x,y) = r
  (xprime, yprime ) = rprime
  
  def omega_2d(kx, ky):
    return np.sqrt( 4*np.sin(kx/2.0)**2 + 4*np.sin(ky/2.0)**2 + mass**2 )
  
  if bc == 'PBC':
    kx = (2*np.array(range(0,L)))*np.pi/L
    ky = (2*np.array(range(0,L)))*np.pi/L
  else: #APBC
    kx = ( 2*np.array(range(0,L)) + 1)*np.pi/L
    ky = ( 2*np.array(range(0,L)) + 1)*np.pi/L
  
  corr = sum( [ sum( np.cos(kx*(x-xprime))*np.cos(kyy*(y-yprime))*omega_2d(kx,kyy) ) for kyy in ky] )
  #corr = sum( [np.cos(kxx*(x-xprime))*np.cos(kyy*(y-yprime))/omega_2d(kxx,kyy) for kxx in kx for kyy in ky] )

  return corr/(2*L**2)

#########################################  getEE_1d  ##########################################
# Calculates the von Neumann / Renyi entanglement entropy from the phi-phi and pi-pi
# correlators in 1D
# Input paramters: 
#    alpha: Renyi index (float)
#    ell: length of region A (integer)
#    X: matrix of phi-phi correlators
#    P: matrix of pi-pi correlators
###############################################################################################
def getEE_1d(alpha, ell, X, P):
  sitesA = range(ell) #the indices of the sites in region A
  XA = X[sitesA][:,sitesA]
  PA = P[sitesA][:,sitesA]
  
  CA_sq = np.matrix(XA)*np.matrix(PA)
  Ev = np.sqrt(np.linalg.eigvals(CA_sq)) #spectrum of eigenvalues of CA_sq
  
  S_alpha = 0
  Ev_new = np.array([e.real for e in Ev if (e.real - 1./2.)>0])
  #Ev_new = np.array([e for e in Ev if (e - 1./2.)>0])
  if alpha == 1:
    S_alpha = np.sum( (Ev_new+1./2)*np.log(Ev_new+1./2.) - (Ev_new-1./2.)*np.log(Ev_new-1./2) )
    #S_alpha = np.sum( (Ev+1./2)*np.log(Ev+1./2.) - (Ev-1./2.)*np.log(Ev-1./2) )
  else:
    S_alpha = 1.0/(alpha-1.0)*np.sum( np.log( (Ev_new+1./2)**alpha - (Ev_new-1./2.)**alpha ) )
    #S_alpha = 1.0/(alpha-1.0)*np.sum( np.log( (Ev+1./2)**alpha - (Ev-1./2.)**alpha ) )
  return S_alpha

###############################################################################################
###########################################  main  ############################################
###############################################################################################

###### Read input from file: ######
#L     = 6
bc    = 'PBC'
mass  = 0.1
alpha = 2
#L, bc_x, bc_y, mass, alpha = readParams("input.txt")
###################################

t1 = time.clock() #for timing
  
filename = "EE_alpha" + str(alpha) + "_2D_" + bc + "_mass" + decimalStr(mass) + ".txt"
fout = open(filename, 'w')

LList = range(10,11)
for L in LList:
  N = L*L
  print 'L = ' + str(L)

  #Calculate X and P matrices:
#   X = np.matrix(np.zeros((N,N)))
#   P = np.matrix(np.zeros((N,N)))
#   for i in range(N):
#     for j in range(i,N):
#       ix = i%L
#       iy = (i-ix)/L
#       jx = j%L
#       jy = (j-jx)/L
#     
#       X[i,j] = getPhiPhiCorrelator_2d(L, bc, mass, (ix,iy), (jx,jy))
#       X[j,i] = X[i,j]
#       P[i,j] = getPiPiCorrelator_2d(L, bc, mass, (ix,iy), (jx,jy))
#       P[j,i] = P[i,j]
#       #print "X[" + str(i) + "," + str(j) + "] = " + str(X[i,j])

  LA_x = L/2
  LA_y = L
  X_from0 = np.zeros((LA_x,LA_y))
  P_from0 = np.zeros((LA_x,LA_y))
  for x in range(LA_x):
    for y in range(LA_y):
      X_from0[x,y] = getPhiPhiCorrelator_2d(L, bc, mass, (0,0), (x,y))
      P_from0[x,y] = getPiPiCorrelator_2d(L, bc, mass, (0,0), (x,y))

  sitesA = np.zeros(LA_x*LA_y).tolist() #the indices of the sites in region A
  count = 0
  for x in range(0,LA_x):
    for y in range(0,LA_y):
      sitesA[count] = y*L + x
      count = count + 1

#   XA = X[sitesA][:,sitesA]
#   PA = P[sitesA][:,sitesA]
  
  NA = len(sitesA)
  XA = np.zeros((NA,NA))
  PA = np.zeros((NA,NA))
  for iA, sitei in enumerate(sitesA):
    for jA, sitej in enumerate(sitesA): 
      xi = sitei%L
      yi = (sitei-xi)/L
      
      xj = sitej%L
      yj = (sitej-xj)/L
      
      XA[iA,jA] = X_from0[abs(xj-xi),abs(yj-yi)]
      PA[iA,jA] = P_from0[abs(xj-xi),abs(yj-yi)]
  
  
  CA_sq = np.matrix(XA)*np.matrix(PA)
  Ev = np.sqrt(np.linalg.eigvals(CA_sq)) #spectrum of eigenvalues of CA_sq

  S_alpha = 0
  Ev_new = np.array([e.real for e in Ev if (e.real - 1./2.)>0])
  #Ev_new = np.array([e for e in Ev if (e - 1./2.)>0])
  if alpha == 1:
    S_alpha = np.sum( (Ev_new+1./2)*np.log(Ev_new+1./2.) - (Ev_new-1./2.)*np.log(Ev_new-1./2) )
    #S_alpha = np.sum( (Ev+1./2)*np.log(Ev+1./2.) - (Ev-1./2.)*np.log(Ev-1./2) )
  else:
    S_alpha = 1.0/(alpha-1.0)*np.sum( np.log( (Ev_new+1./2)**alpha - (Ev_new-1./2.)**alpha ) )
    #S_alpha = 1.0/(alpha-1.0)*np.sum( np.log( (Ev+1./2)**alpha - (Ev-1./2.)**alpha ) )
  print "  " + str(S_alpha)
  fout.write("%d %.15f"%(L,S_alpha)+'\n')  #Save result to file

#regA_list = np.array(range(1,L))
#EE_list = np.zeros(len(regA_list))
# for i,ell in enumerate(regA_list):
#   #print "ell = " + str(ell)
#   EE_list[i] = getEE_1d(alpha, ell, X, P)
#   fout.write("%d %.15f"%(ell,EE_list[i])+'\n')  #Save result to file
#   #print EE_list[i]
#   #print

# x = np.log( L/np.pi * np.sin(np.pi*regA_list/L) )
# plt.plot(x,EE_list,'o-')
# 
# x2 = np.linspace(0,3.5,10)
# y2 = x2/3.0 + 0.6
# plt.plot(x2,y2)
# plt.show()  

fout.close()

t2 = time.clock()
print "Total elapsed time: " + str(t2-t1) + " sec."
