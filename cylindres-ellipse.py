# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#!/usr/bin/python
# Copyright (C) 2013 FlowKit Ltd, Lausanne, Switzerland
# E-mail contact: contact@flowkit.com
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License, either
# version 3 of the License, or (at your option) any later version.

#
# 2D flow
#

from numpy import *; 
from numpy.linalg import *; 
from tkinter import *
import matplotlib.pyplot as plt; 
from matplotlib import cm
from matplotlib import animation

###### Flow definition #########################################################

maxIter = 29000 # Total number of time iterations initialement 200000.
Re      = 100.0  # Reynolds number.
nx = 780; ny = 360; ly=ny-1.0; q = 9 # Lattice dimensions and populations.
uLB     = 0.04                       # Velocity in lattice units.
theta=0.
nom_anim = "DefaultETU2020"
Cas=0
d=0.
#pert=1.e-4
#pert=0.
Ray = 30
beta = 4
cs = 1

##### Fenêtre de saisie ####

fenetre = Tk()
fenetre.geometry("800x250")
var_cas = IntVar()
var_HP = IntVar()
var_Dcyl=IntVar()
cas_cyl = Radiobutton(fenetre, text="Cylindre elipsoïdal", variable=var_cas, value=0)
cas_helm = Radiobutton(fenetre, text="Helmholtz", variable=var_cas, value=1)
cas_onzecyl=Radiobutton(fenetre, text="5 Cyl. circ.", variable=var_cas, value=2)
cas_vingtdeuxcyl=Radiobutton(fenetre, text="22 Cyl. circ.+Elipse", variable=var_cas, value=3)
cas_cyl.grid(row=1, column=1)
cas_helm.grid(row=1, column=2)
cas_onzecyl.grid(row=1,column=3)
cas_vingtdeuxcyl.grid(row=1,column=4)
cas_HP = Checkbutton(fenetre, text="HP", variable = var_HP)
cas_HP.grid(row=1, column=5)
# champ_label = Label(fenetre, text="Entrer le nom du fichier mp4")
# champ_label.grid(row=2, column=1)
# ligne_texteNF = Entry(fenetre, textvariable = var_texteNF, width=30)
# ligne_texteNF.grid(row=3, column=1)
Re_label = Label(fenetre, text="Entrer le Reynolds")
Re_label.grid(row=4, column=1)
var_texteRe = StringVar(fenetre,"220.")
ligne_texteRe = Entry(fenetre, textvariable = var_texteRe, width=10)
ligne_texteRe.grid(row=4, column=2)
Pert_label = Label(fenetre, text="Entrer l'amplitude de la perturbation du HP: ")
Pert_label.grid(row=5, column=1)
var_textePert = StringVar(fenetre,"0.04")
ligne_textePert = Entry(fenetre, textvariable = var_textePert, width=10)
ligne_textePert.grid(row=5, column=2)
Freq_label = Label(fenetre, text="Entrer la demi période de pulsation HP (multiple de 100): ")
Freq_label.grid(row=6, column=1)
var_texteFreq = StringVar(fenetre,"500")
ligne_texteFreq = Entry(fenetre, textvariable = var_texteFreq, width=10)
ligne_texteFreq.grid(row=6, column=2)
Theta_label = Label(fenetre, text="Entrer l'inclinaison du profil Theta ou la longueur du tube en pixels")
Theta_label.grid(row=7, column=1)
var_texteTheta = StringVar()
ligne_texteTheta = Entry(fenetre, textvariable = var_texteTheta, width=10)
ligne_texteTheta.grid(row=7, column=2)
Ray_label = Label(fenetre, text="Entrer le rayon du tube en pixels (Helmholtz)")
Ray_label.grid(row=8, column=1)
var_texteRay = StringVar(fenetre,"30")
ligne_texteRay = Entry(fenetre, textvariable = var_texteRay, width=10)
ligne_texteRay.grid(row=8, column=2)
def Validrecup():
    global Cas, nom_anim, Thetastr, Restr, theta, Ray, Re, ltpx, AmpHP, HP, nom_fichstr
    Cas=var_cas.get()
    CasHP=var_HP.get()
    if (CasHP==1):
        strHP="HP"
    else:
        strHP=""
        
    if(Cas==1):
        nom_anim = 'HelmD'+strHP
    else:
        nom_anim = 'CylT'+strHP
        if(Cas==2):
            nom_anim='Cyl11Eli'+strHP
        if(Cas==3):
            nom_anim='Cyl22Eli'+strHP
    Restr = ligne_texteRe.get()
    Pertstr = ligne_textePert.get()
    Thetastr = ligne_texteTheta.get()
    Raystr = ligne_texteRay.get()
    Freqstr = ligne_texteFreq.get()
    if Thetastr=="": Thetastr = "0."
    Re = float(Restr)
    AmpHP = float(Pertstr)
    theta = float(Thetastr)
    Ray = int(Raystr)
    if(Cas==1):ltpx = int(Thetastr)
    HP = float(Freqstr)
    nom_fichstr = nom_anim+"_"+Restr+"_"+Raystr+"_"+Thetastr+"_"+Pertstr+"_"+Freqstr+".mp4"
    var_label = Label(fenetre, text = nom_fichstr)
    var_label.grid(row=9, column=1)
bouton_valider = Button(fenetre, text="Valider", command=Validrecup)
bouton_valider.grid(row=10, column=1)
bouton_quitter = Button(fenetre, text="Lancer le calcul", command=fenetre.quit)
bouton_quitter.grid(row=10, column=2)
fenetre.mainloop()
fenetre.destroy()

###### Lattice Constants #######################################################
c = array([(x,y) for x in [0,-1,1] for y in [0,-1,1]]) # Lattice velocities.
t = 1./36. * ones(q)                                   # Lattice weights.
t[asarray([norm(ci)<1.1 for ci in c])] = 1./9.; t[0] = 4./9.
noslip = [c.tolist().index((-c[i]).tolist()) for i in range(q)] 
i1 = arange(q)[asarray([ci[0]<0  for ci in c])] # Unknown on right wall.
i2 = arange(q)[asarray([ci[0]==0 for ci in c])] # Vertical middle.
i3 = arange(q)[asarray([ci[0]>0  for ci in c])] # Unknown on left wall.

###### Function Definitions ####################################################
sumpop = lambda fin: sum(fin,axis=0) # Helper function for density computation.
def equilibrium(rho,u):              # density computation function.
    cu   = 3.0 * dot(c,u.transpose(1,0,2))
    usqr = 3./2.*(u[0]**2+u[1]**2)
    feq = zeros((q,nx,ny))
    for i in range(q): feq[i,:,:] = rho*t[i]*(1.+cu[i]+0.5*cu[i]**2-usqr)
    return feq

###### Setup: obstacles and velocity inlet with perturbation ########
if (Cas == 0):
    ## Cylindre base eliptique##
    cx = nx/5; cy=ny*1/3; r = ny/40;          # Coordinates of the cylinder.
    obstacle = fromfunction(lambda x,y: 0.05*(x*cos(theta)+y*sin(theta)-cx)**2+(y*cos(theta)+ x*sin(theta)-cy)**2 <r**2 , (nx,ny))
if (Cas == 0):
    ## Cylindre base eliptique##
    cx = nx/5; cy=ny*1/3; r = ny/10;          # Coordinates of the cylinder.
    obstacle = fromfunction(lambda x,y: 0.05*(x*cos(theta)+y*sin(theta)-cx)**2+(y*cos(theta)+ x*sin(theta)-cy)**2 <r**2 , (nx,ny))
if ((Cas == 2) or (Cas == 3)):
    ## Cylindres de bases circulaires##
    cx = nx/10; cy1 = ny/6; cy2= ny/3; cy3 = ny/2; cy4 = 2*ny/3; cy5 = 5*ny/6; r = ny/25; cx2 = cx + 8*r; 
    deccy = ny/24; relips = ny/10
    cxelips=nx/2; cyelips=ny/2
    #if (theta < -0.1):
        #cyelips=ny/3
    # Coordinates of the cylinders.
    if(Cas == 2):
        obstacle = fromfunction(lambda x,y: (((x-cx)**2+(y-cy1)**2 <r**2)\
                                         | ((x-cx)**2+(y-cy2)**2 <r**2) | ((x-cx)**2+(y-cy3)**2 <r**2)\
                                         | ((x-cx)**2+(y-cy4)**2 <r**2) | ((x-cx)**2+(y-cy5)**2 <r**2)\
                                         | (y < 5) | (y > 355)), (nx,ny))
    if(Cas == 3):
        obstacle = fromfunction(lambda x,y: (((x-cx)**2+(y-cy1)**2 <r**2)\
                                         | ((x-cx)**2+(y-cy2)**2 <r**2) | ((x-cx)**2+(y-cy3)**2 <r**2)\
                                         | ((x-cx)**2+(y-cy4)**2 <r**2) | ((x-cx)**2+(y-cy5)**2 <r**2)\
                                         | ((x-cx2)**2+(y-cy1-deccy)**2 <r**2)\
                                         | ((x-cx2)**2+(y-cy2-deccy)**2 <r**2) | ((x-cx2)**2+(y-cy3-deccy)**2 <r**2)\
                                         | ((x-cx2)**2+(y-cy4-deccy)**2 <r**2) | ((x-cx2)**2+(y-cy5-deccy)**2 <r**2)\
                                         | (0.05*(x*cos(theta)+y*sin(theta)-cxelips)**2
                                            +(y*cos(theta)+ x*sin(theta)-cyelips)**2 <relips**2)\
                                         | (y < 5) | (y > 355)), (nx,ny))
                                                                                      
if (Cas == 1):
    ## Helmholtz ##    Tx1 = 200; Tx2 = 220; Tx3 = Tx2 + ltpx; Ry1 = 160; Ry2 = 220; Ry3 = 156; Ry4 = 216; Ty2 = 4; Ty3 = 356 
    # Obstacle
    Tx1 = 200; Tx2 = 220; Tx3 = Tx2 + ltpx; Ry1 = 180 - Ray; Ry2 = 180 + Ray; Ry3 = Ry1 - 4; Ry4 = Ry2 - 4; Ty2 = 4; Ty3 = 356 
    # Obstacle
    #r = Tx3-Tx2 # Longueur du tube en pixels ####
    r = Ry2-Ry1 # Internal tube diameter in pixels
    #obstacle = fromfunction(lambda x,y: ((((x <= Tx2)&( x >= Tx1))&((y <= Ry1)|(y >= Ry2)))|((((y <= Ry1)&(y >= Ry3))
    #|((y <= Ry2)&(y >= Ry4)))&((x >= Tx1)&(x < Tx3))))|((((x<=Tx1)&(y<=Ty2))|((x<=Tx1)&(y>=Ty3)))) , (nx,ny))
    obstacle = fromfunction(lambda x,y: ((((x < Tx2)&( x > Tx1))&((y < Ry1)|(y > Ry2)))
                                         |((((y < Ry1)&(y > Ry3))|((y < Ry2)&(y > Ry4)))&\
    ((x > Tx1)&(x < Tx3))))|((((x<Tx2)&(y<Ty2))|((x<Tx2)&(y>Ty3)))) , (nx,ny))
    omega0helmholtz = cs*sqrt(Ray*2/(ltpx*(Ty3-Ty2)*Tx1))
    print ('omega 0 résonnance de helmhotz unités LBM', omega0helmholtz)
nulb = uLB*r/Re; omega = 1.0 / (3.*nulb+0.5); # Relaxation parameter.
vel = fromfunction(lambda d,x,y: (1-d)*uLB,(2,nx,ny)) #*(1.0+pert*sin(beta*y/ly*2*pi))
feq = equilibrium(1.0,vel); fin = feq.copy()

###### Main time loop ##########################################################
fig = plt.figure()
ims = []
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Patrick VIGLIANO'), bitrate=1800)

####################################
nom_fichstr=''
for time in range(maxIter):
    fin[i1,-1,:] = fin[i1,-2,:] # Right wall: outflow condition.
    rho = sumpop(fin)           # Calculate macroscopic density and velocity.
    u = dot(c.transpose(), fin.transpose((1,0,2)))/rho
    u[:,0,:] =vel[:,0,:] # Left wall: compute density from known populations.
    #u[:,0,:] =vel[:,0,:]*(1. + AmpHP * sin(time * pi / HP)) # Left wall: vitesse oscillante.
    #u[:,0,:] =vel[:,0,:]*(1. + AmpHP * sin(time * omega0helmholtz * 2 * pi)) # Left wall: vitesse oscillante pulsation résonnance?
    #u[:,0,:] =vel[:,0,:]*(1. + AmpHP * sin(time * omega0helmholtz))
     # Left wall: vitesse oscillante fréquence résonnance?
    rho[0,:] = 1./(1.-u[0,0,:]) * (sumpop(fin[i2,0,:])+2.*sumpop(fin[i1,0,:]))
    feq = equilibrium(rho,u) # Left wall: Zou/He boundary condition.
    fin[i3,0,:] = fin[i1,0,:] + feq[i3,0,:] - fin[i1,0,:]
    fout = fin - omega * (fin - feq)  # Collision step.
    for i in range(q): fout[i,obstacle] = fin[noslip[i],obstacle]
    for i in range(q): # Streaming step.
        fin[i,:,:] = roll(roll(fout[i,:,:],c[i,0],axis=0),c[i,1],axis=1)
    if ((time >= 1) & (time%100 == 1)): # Visualization
    #if (time%100 == 1): # Visualization
       # im=plt.imshow(((u[0]**2+u[1]**2)+rho*cs**2).transpose(),cmap='rainbow', animated=True) #pression
        im = plt.imshow(sqrt(u[0]**2+u[1]**2).transpose(),cmap=cm.jet, animated=True) #vitesse
        # if (time/100 == 1): plt.savefig("test.png")
        ims.append([im])
        #ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=1000)
        #ani.save(nom_fichstr, writer=writer)
        print(time,' sur ', maxIter)
plt.show()

