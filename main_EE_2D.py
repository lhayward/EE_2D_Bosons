import argparse
import numpy as np
import os.path  #to check if file exists
#import sys  #for sys.stdout.flush()
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

#########################################  omega_2d  ##########################################
# Calculates omega for given kx, ky and mass 
###############################################################################################
def omega_2d(kx, ky, mass):
  return np.sqrt( 4*np.sin(kx/2.0)**2 + 4*np.sin(ky/2.0)**2 + mass**2 )

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
  
  if bc == 'PBC':
    kx = (2*np.array(range(0,L)))*np.pi/L
    ky = (2*np.array(range(0,L)))*np.pi/L
  else: #APBC
    kx = ( 2*np.array(range(0,L)) + 1)*np.pi/L
    ky = ( 2*np.array(range(0,L)) + 1)*np.pi/L
  
  corr = sum( [ sum( np.cos(kx*(x-xprime))*np.cos(kyy*(y-yprime))/omega_2d(kx,kyy,mass) ) for kyy in ky] )
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
  
  if bc == 'PBC':
    kx = (2*np.array(range(0,L)))*np.pi/L
    ky = (2*np.array(range(0,L)))*np.pi/L
  else: #APBC
    kx = ( 2*np.array(range(0,L)) + 1)*np.pi/L
    ky = ( 2*np.array(range(0,L)) + 1)*np.pi/L
  
  corr = sum( [ sum( np.cos(kx*(x-xprime))*np.cos(kyy*(y-yprime))*omega_2d(kx,kyy,mass) ) for kyy in ky] )
  #corr = sum( [np.cos(kxx*(x-xprime))*np.cos(kyy*(y-yprime))/omega_2d(kxx,kyy) for kxx in kx for kyy in ky] )

  return corr/(2*L**2)


#########################################  readArray  #########################################
# Takes a string representation of an array and removes the '[' and ']' characters
###############################################################################################
def readArray(line):
  start = max( 0, line.find('[') )
  end = min( len(line), line.find(']') )
  return line[start+1:end] 
####### end readArray(line) function #######

########################################  readParams  #########################################
# Reads in values from input file. Input should have the following form:
#
# L    = ___ (int)
# bc_x = ___ (string: 'PBC' or 'APBC')
# bc_y = ___ (string: 'PBC' or 'APBC')
# mass = ___ (float)
# alpha     = [ ___  ___  ___ ...] (list of floats in square brackets with items separated by spaces)
###############################################################################################
def readParams(filename):
  L    = 1
  bc_x = 'PBC'
  bc_y = 'PBC'
  massterm = 1
  alpha    = [0]
  if os.path.isfile(filename):
    fin = open(filename,'r')
    
    line=fin.readline()
    L = int(line[ max(0,line.find('='))+1:])
    
    line=fin.readline()
    bc_x = line[ max(0,line.find('='))+1:].strip()
    
    line=fin.readline()
    bc_y = line[ max(0,line.find('='))+1:].strip()
    
    line=fin.readline()
    mass = float(line[ max(0,line.find('='))+1:])
    
    line=fin.readline()
    alpha = [float(a) for a in readArray(line).split()]
    
    fin.close()
  return L, bc_x, bc_y, mass, alpha

###############################################################################################
###########################################  main  ############################################
###############################################################################################

parser=argparse.ArgumentParser(description="Code to calculate EE for 2D free bosons")
parser.add_argument('-f', '--file')
args=parser.parse_args()

###### Read input from file: ######
inFile = "input"
if args.file != None:
  inFile = inFile + "_" + args.file
inFile = inFile + ".txt"
print "Input file: %s" %inFile

L     = 10
bc    = 'APBC'
mass  = 0
alpha = 1
L_test, bc_x, bc_y, mass_test, alpha_test = readParams("input.txt")
###################################

print "L_test     = %d" %L_test
print "BC along x = %s" %bc_x
print "BC along y = %s" %bc_y
print "mass_test  = %f" %mass_test
print "alpha_test = " + str(alpha_test)

t1 = time.clock() #for timing

  
filename = "EE_2D_" + bc + "_L" + str(L) + "_mass" + decimalStr(mass) + ".txt"
fout = open(filename, 'w')

#LList = range(6,13)
print 'L = ' + str(L)
for LA_x in range(1,L):
  print '  LA = ' + str(LA_x)

#  N = L*L
#  #Calculate X and P matrices:
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

  #LA_x = L/2
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
  fout.write("%d %.15f"%(LA_x,S_alpha)+'\n')  #Save result to file

fout.close()

t2 = time.clock()
print "Total elapsed time: " + str(t2-t1) + " sec."
