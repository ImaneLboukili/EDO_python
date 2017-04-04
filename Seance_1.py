#! usr/bin/env python
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
import functools



def xp(t):
  return (-t)


def xpm(g, tau, W, x, t):
  return (-x+np.tanh(g*np.dot(W,x)))*1./tau
  
  #Fait la meme chose que :
  """
  y = -x
  l = len(x)
  for i in range(l):
    sig = 0
    for j in range(l):
      sig += W[i][j] * x[j]

    y[i] += np.tanh(g * sig)
  """
  #return y/tau

def part_xpm(g, tau, W):
  return functools.partial(xpm, g, tau, W)



def Euler(func, dt, x0, t0=0, tmax=1):
  it = t0 + dt
  i = 0
  
  dots = [(t0, x0)]

  while it <= tmax :
    dots.append((it*1.0, func(dots[i][1],it)*dt + dots[i][1]))
    it += dt
    i += 1

  return dots


def simul(N, g, tau, x0, k):
  W = np.zeros((N,N))
  W.fill(k)

  #Recuperation de la fonction utilisable par la methode d'Euler
  f = part_xpm(g, tau, W)

  #Approximation de x(t) par Euler
  d = Euler(f, 0.001, x0, t0=0, tmax=10)

  
  #Recuperation des donnees pour plot
  x=[]
  y=[[] for _ in range(len(d[0][1]))]
  
  for i in range(len(d)):
    x.append(d[i][0])
    for j in range(len(d[i][1])):
      y[j].append(d[i][1][j])
  
  #plot
  for i in range(len(y)):
    plt.plot(x,y[i])
  plt.title("N="+str(N)+", tau="+str(tau)+", k ="+str(k))
  plt.show()


def simul_W(N, g, tau, x0, W):
  #Recuperation de la fonction utilisable par la methode d'Euler
  f = part_xpm(g, tau, W)

  #Approximation de x(t) par Euler
  d = Euler(f, 0.001, x0, t0=0, tmax=10)

  
  #Recuperation des donnees pour plot
  x=[]
  y=[[] for _ in range(len(d[0][1]))]
  
  for i in range(len(d)):
    x.append(d[i][0])
    for j in range(len(d[i][1])):
      y[j].append(d[i][1][j])
  
  #plot
  for i in range(len(y)):
    plt.plot(x,y[i])
  plt.title("N="+str(N)+", tau="+str(tau)+", g ="+str(g))
  plt.show()



if __name__ == '__main__':
  
  """
  #Test of Euler Method with 
  for x0i in (0,3,5):
    d = Euler(xp, 0.1, x0i, t0=0, tmax=10)
    x=[]
    y=[]
    
    for i in range(len(d)):
      x.append(d[i][0])
      y.append(d[i][1])

    #print d 
    #print x
    #print y
      
  #  plt.plot(x,y)
  
  #plt.show()

  """
  print "\n\nCOUCOU\n\nIl faut quitter les graphiques à la main et les suivants apparaissent à la suite.\n Les valeurs des paramètres sont indiquéees dans le titre.\n"


  N = 10
  g = 1
  tau = 10
  x0 = np.linspace(-0.5,0.5,N)

  #pour differentes valeurs de k
  for k in (-1,-0.5,0,0.5,1):
    simul(N,g,tau,x0,k)

  k=0.5

  #pour differentes valeurs de tau
  for tau in (2,5,10,20,-2,-10,-20):
    simul(N,g,tau,x0,k)

  tau = 10

  #pour W avec des valeurs selon une gaussienne.
  W = np.random.randn(N,N)

  simul_W(N,g,tau,x0,W)

  for g in range(1,10):
    simul_W(N,g,tau,x0,W)




