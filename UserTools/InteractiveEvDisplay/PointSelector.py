### The purpose of this script is 
###

import Store
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.image as img
#from PIL import Image

### constants ####################################
ELLIPSE_CTR_X = 448.
ELLIPSE_CTR_YT = 254.
ELLIPSE_CTR_YB = 620.
ELLIPSE_MAJOR = 88.
ELLIPSE_MINOR = 88.
BULK_TOP_EDGE = 342.
BULK_BOT_EDGE = 532.
BULK_LEFT_EDGE = 172.
BULK_RIGHT_EDGE = 724.

#root constants
size_top_drawing = 0.1
max_y = 1.532699    #m
min_y = -1.631441   #m
rho = 1.137         #m; avg


### function defs ################################
def tellme(s):
  print(s)
  plt.title(s,fontsize=14)
  plt.draw()

def isInsideTankBulk(screen_x, screen_y):
  if ((screen_x > BULK_LEFT_EDGE) and (screen_x < BULK_RIGHT_EDGE) and (screen_y > BULK_TOP_EDGE) and (screen_y < BULK_BOT_EDGE)):
    return True
  else:
    return False

def isInsideTankCaps(screen_x, screen_y):
  c = (ELLIPSE_MAJOR*ELLIPSE_MINOR)**2 + 1.
  cap = ""
  #TODO:multiply both sides by rx^2*ry^2 and compare with
  # this value instead of 1.
  if (screen_y < BULK_TOP_EDGE):   #top; I know, it's reversed in python..
    c = ((screen_x-ELLIPSE_CTR_X)*ELLIPSE_MINOR)**2 + ((screen_y-ELLIPSE_CTR_YT)*ELLIPSE_MAJOR)**2
    cap = "top"
  elif (screen_y > BULK_BOT_EDGE):   #bottom
    c = ((screen_x-ELLIPSE_CTR_X)*ELLIPSE_MINOR)**2 + ((screen_y-ELLIPSE_CTR_YB)*ELLIPSE_MAJOR)**2
    cap = "bottom"
  if (c <= (ELLIPSE_MAJOR*ELLIPSE_MINOR)**2):
    return(True,cap)
  else:
    cap = ""
    return(False,cap)     #return empty string

def convert2Dto3D(screen_x, screen_y, pmt_pos, PYIMG_MAX_X, PYIMG_MAX_Y):
  root_x = screen_x/PYIMG_MAX_X
  root_y = (PYIMG_MAX_Y-screen_y)/PYIMG_MAX_Y
  tank_x = 0.
  tank_y = 0.
  tank_z = 0.

  if (pmt_pos == "top"):
    tank_x = (0.5-root_x)*tank_radius/size_top_drawing
    tank_y = 1.38805
    tank_z = (0.5+(max_y/tank_radius+1)*size_top_drawing-root_y)*tank_radius/size_top_drawing
  elif (pmt_pos == "bottom"):
    tank_x = (0.5-root_x)*tank_radius/size_top_drawing
    tank_y = -1.77609
    tank_z = (-0.5+(abs(min_y)/tank_radius+1)*size_top_drawing+root_y)*tank_radius/size_top_drawing
  else:   #barrel
    phi = (root_x-0.5)/size_top_drawing
    if (screen_x < 308.):   #x>0,z<0
      tank_x = -rho*np.sin(phi)
      tank_z = rho*np.cos(phi)
    elif ((screen_x > 308.) and (screen_x < 448.)):   #x,z>0
      tank_x = -rho*np.sin(phi)
      tank_z = rho*np.cos(phi)
    elif ((screen_x > 448.) and (screen_x < 584.)):   #x<0,z>0
      tank_x = -rho*np.sin(phi)
      tank_z = rho*np.cos(phi)
    elif (screen_x > 584.):   #x,z<0
      tank_x = -rho*np.sin(phi)
      tank_z = rho*np.cos(phi)
    else:
      print("!! Out of bounds!")

    tank_y = (root_y-0.5)*tank_radius/size_top_drawing
  return((tank_x,tank_y,tank_z))

##################################################
def Initialise(pyinit):
  return 1

def Finalise():
  return 1

def Execute():
  plt.ion()

  # Grab event display from root file or Store (is it possible?)
  # for now, just open image file
  img_file = Store.GetStoreVariable('Config','ImageFile')
  print("Using event display found in: " + img_file)
  evd = img.imread(img_file)
  print(evd.shape)  #debug
  PYIMG_MAX_X = evd.shape[1]
  PYIMG_MAX_Y = evd.shape[0]
  plt.imshow(evd)
  plt.show()

  # Grab tank radius from Store
  #tank_radius = Store.GetStoreVariable("ANNIEGeometry", "TankRadius")
  tank_radius = 1.37504  #m

  while True:
    pts = []
    line, = plt.plot(pts, pts,'ro')   #initialize Line2D object
    while (len(pts) < 1):
      tellme("Select candidate exit point")
      pts = np.asarray(plt.ginput(1,timeout=-1))
      #plt.plot(pts[0][0],pts[0][1],'ro')
      line.set_data(pts[0][0],pts[0][1])
      if (len(pts) < 1):
        tellme("Too few points, starting over")
        time.sleep(1)

    #ph = plt.fill(pts[:,0], pts[:,1], 'r', lw=2)

    # if out of bounds, tell user to try again
    screen_x = pts[0][0]
    screen_y = pts[0][1]
    tankCap = isInsideTankCaps(screen_x, screen_y)
    tankBulk = isInsideTankBulk(screen_x, screen_y)
    if (tankBulk != True and tankCap[0] != True):
      print("Oopsies! Out of bounds. Try again.")
      line.remove()
      continue

    tellme("Happy? Press y for yes, click on plot for no")

    if (plt.waitforbuttonpress()):
      break

    #for p in ph:
    #  p.remove()
    line.remove()

  # Convert screen coordinates to physical detector coordinates
  root_x = screen_x/PYIMG_MAX_X
  root_y = (PYIMG_MAX_Y-screen_y)/PYIMG_MAX_Y
  tank_x = 0.
  tank_y = 0.
  tank_z = 0.

  if (tankCap[1] == "top"):
    tank_x = (0.5-root_x)*tank_radius/size_top_drawing
    tank_y = 1.38805
    tank_z = (0.5+(max_y/tank_radius+1)*size_top_drawing-root_y)*tank_radius/size_top_drawing
  elif (tankCap[1] == "bottom"):
    tank_x = (0.5-root_x)*tank_radius/size_top_drawing
    tank_y = -1.77609
    tank_z = (-0.5+(abs(min_y)/tank_radius+1)*size_top_drawing+root_y)*tank_radius/size_top_drawing
  else:   #barrel
    #rho = Store.GetStoreVariable('ANNIEEvent', 'evDisplayRho')
    phi = (root_x-0.5)/size_top_drawing
    if (screen_x < 308.):   #x>0,z<0
      tank_x = -rho*np.sin(phi)
      tank_z = rho*np.cos(phi)
    elif ((screen_x > 308.) and (screen_x < 448.)):   #x,z>0
      tank_x = -rho*np.sin(phi)
      tank_z = rho*np.cos(phi)
    elif ((screen_x > 448.) and (screen_x < 584.)):   #x<0,z>0
      tank_x = -rho*np.sin(phi)
      tank_z = rho*np.cos(phi)
    elif (screen_x > 584.):   #x,z<0
      tank_x = -rho*np.sin(phi)
      tank_z = rho*np.cos(phi)
    else:
      print("!! Out of bounds!")

    tank_y = (root_y-0.5)*tank_radius/size_top_drawing

  print("Physical coordinates: "+str(tank_x)+","+str(tank_y)+","+str(tank_z)+" [m]")

  # Store these coordinates in Store to be used in vertex reco
  #Store.SetStoreVariable('Store','ExitPointX',tank_x)
  #Store.SetStoreVariable('Store','ExitPointY',tank_y)
  #Store.SetStoreVariable('Store','ExitPointZ',tank_z)

  # new syntax: Store.GetStoreVariable('StoreName','varName')
  '''
  a=Store.GetStoreVariable('Config','a') 
  b=Store.GetStoreVariable('Config','b')
  c=Store.GetStoreVariable('Config','c')
  print(a)
  print(b)
  print(c)

  Store.SetStoreVariable('CStore','a',69)
  num = Store.GetStoreVariable('CStore','a')
  print(num)
  '''
  print("Done!")
  return 1
