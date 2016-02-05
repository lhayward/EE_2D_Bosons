import argparse
import numpy as np
import os.path  #to check if file exists
import sys  #for sys.stdout.flush()
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

#####################################  getCorrelators_2d  #####################################
# Calculates the phi-phi and pi-pi correlator for a 2D free bosonic system
# Input paramters: 
#    L: length of lattice (integer)
#    bc_x: boundary conditions along x ('PBC' or 'APBC')
#    bc_y: boundary conditions along y ('PBC' or 'APBC')
#    mass: boson mass (float)
#    r, rprime: (x,y) coordinates on lattice of the two phi variables (floats)
###############################################################################################
def getCorrelators_2d(L, bc_x, bc_y, mass, r, rprime):
  bc_err = False
  (x,y) = r
  (xprime, yprime ) = rprime
  
  #BC along x:
  if bc_x == 'PBC':
    kx = (2*np.array(range(0,L)))*np.pi/L
  elif bc_x == 'APBC':
    kx = ( 2*np.array(range(0,L)) + 1)*np.pi/L
  else:
    kx = np.zeros(L)
    print "*** Boundary condition %s along x is not supported ***" %bc_x
    bc_err = True 
  
  #BC along y:
  if bc_y == 'PBC':
    ky = (2*np.array(range(0,L)))*np.pi/L
  elif bc_y == 'APBC':
    ky = ( 2*np.array(range(0,L)) + 1)*np.pi/L 
  else:
    kx = np.zeros(L)
    print "*** Boundary condition %s along y is not supported ***" %bc_y
    bc_err = True 
  
  phiphi = 0
  pipi   = 0
  if not bc_err:
    for kyy in ky:
      omega = omega_2d(kx,kyy,mass)
      phiphi = phiphi + sum( np.cos(kx*(x-xprime))*np.cos(kyy*(y-yprime))/omega )
      pipi   = pipi + sum( np.cos(kx*(x-xprime))*np.cos(kyy*(y-yprime))*omega )

  return phiphi/(2*L**2), pipi/(2*L**2)

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
    alpha = np.array([float(a) for a in readArray(line).split()])
    
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

L, bc_x, bc_y, mass, alpha = readParams(inFile)
###################################

print "L          = %d" %L
print "BC along x = %s" %bc_x
print "BC along y = %s" %bc_y
print "mass       = %f" %mass
print "alpha      = " + str(alpha)

t1 = time.clock() #for timing

filename = "EE_2D_%sx_%sy_L%d_mass%s.txt" %(bc_x,bc_y,L,decimalStr(mass))
fout = open(filename, 'w')

LA_y     = L
LA_x_max = (L+1)/2

#Calculate all needing correlators:
X_from0 = np.zeros((LA_x_max,LA_y))
P_from0 = np.zeros((LA_x_max,LA_y))
for x in range(LA_x_max):
  for y in range(LA_y):
    X_from0[x,y], P_from0[x,y] = getCorrelators_2d(L, bc_x, bc_y, mass, (0,0), (x,y))

#Loop over all cylinders:
for LA_x in range(1,LA_x_max+1):
  print "\nLA = %d" %LA_x
  sys.stdout.flush()

  #Make a list of the sites in this region A:
  sitesA = np.zeros(LA_x*LA_y).tolist() #the indices of the sites in region A
  count = 0
  for x in range(0,LA_x):
    for y in range(0,LA_y):
      sitesA[count] = y*L + x
      count = count + 1
  
  #Calculate XA and PA:
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
  
  #Calculate the matrix CA and its eigenvalues:
  CA_sq = np.matrix(XA)*np.matrix(PA)
  Ev = np.sqrt(np.linalg.eigvals(CA_sq)) #spectrum of eigenvalues of CA_sq
  Ev_new = np.array([e.real for e in Ev if (e.real - 1./2.)>0])
  
  #Calculate the EE for each alpha:
  S_alpha = np.zeros(len(alpha))
  for i, n in enumerate(alpha):
    if n == 1:
      S_alpha[i] = np.sum( (Ev_new+1./2)*np.log(Ev_new+1./2.) - (Ev_new-1./2.)*np.log(Ev_new-1./2) )
    else:
      S_alpha[i] = 1.0/(n-1.0)*np.sum( np.log( (Ev_new+1./2)**n - (Ev_new-1./2.)**n ) )
  #end alpha loop
  print "  " + str(S_alpha)
  
  #Save results to file:
  fout.write("%d" %LA_x)
  for Sn in S_alpha:
    fout.write(" %.15f" %Sn)
  fout.write("\n")
  fout.flush()
#End loop over cylinders

fout.close()

t2 = time.clock()
print "Total elapsed time: " + str(t2-t1) + " sec."
