# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 11:17:19 2016

@author: Eric Schmidt
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy.integrate import quad


Ne = 2153 # Number of particles
N = 1000 # 
T = 1/10 # 
t = np.linspace(0.0, N*T, N)
tf = np.linspace(0.0,1.0/(2.0*T),N/2)
w = ((6/5)*np.pi)*10**6 # (rad/s) Oscillating frequency 
#w = 2*np.pi
    
phi_x = -np.pi/8 # (rad) Phase-offset of x-position mean oscillations
phi_sigma = np.pi/2 # (rad) Phase-offset of beam distribution width oscillations
x_mid = 0 # (m) Middle of the acceptance function
x_0 = 0 # (m) Initial beam distribution mean
sigma_0 = 14.5*10**-3 # (m) Initial beam distribution width
D_0 = 80*10**-3 # (m) Maximum physical width of beam
A_x = 2*10**-3 # (m) Initial beam distribution mean oscillation amplitude
A_sigma = 5.5*10**-3 # (m) Initial beam distribution width oscillation amplitude

k_0 = 0.1 # () Linear acceptance function constant
k_1 = 0.1 # () Quadratic acceptance function constant

# Other constants used:
# AD : Acceptance value
# a : Normalization variable to ensure constant number of particles at each step
# ul : Upper limit when integrating over x-position
# ll : Lower limit when integrating over x-position

def main():
    
    # Pick one of the following
    
#    xlinear(noplot=0)
#    sigmalinear(noplot=0)
    combinedLinear()
#    
#    xquad(noplot=0)
#    sigmaquad(noplot=0)
#    combinedQuad()

def xlinear(noplot):

    # Initialize variables    
    
    AD = np.zeros(len(t))

    x_B = x_0 + A_x*np.cos(w*t + phi_x)
    
    D = None
    D = D_0
    sigma = None
    sigma = sigma_0
    
    a = np.zeros(len(t))
    
    ul = x_B + D/2
    ll = x_B - D/2
    
    i = 0
    
    while i < len(t):
        
        # Finds the new normalization variable
        temp = quad(integranda, ll[[i]], ul[[i]], args=(x_B[i],sigma))
        a[i] = N/(temp[0])
        
        # Integrates over x to find the acceptance
        temp1 = quad(integrandADlinear, ll[[i]], ul[[i]], args=(a[i],x_B[i],sigma))
        AD[i] = temp1[0]
        
#        Nt = quad(integrandN, ll[[i]], ul[[i]], args=(a[i],x_B[i],sigma))
#        print(Nt[0])
        
        i = i + 1
    
    title = "Postion vs. Time"
    stitle = "Gaussian Distribution: Linear Position"
    ylabel = "$x_m$ (m)"
    save_title = "GF_Linear_Position"
    
    if noplot==0:
        plotSingle(AD,x_B,ylabel,title,save_title,stitle)
    
    return AD
    
def sigmalinear(noplot):
    
    AD = np.zeros(len(t))    
    a = np.zeros(len(t))

    sigma = sigma_0 + A_sigma*np.cos(w*t + phi_sigma)
    
    x_B = None
    x_B = x_0
    
    ul = x_B + D_0/2
    ll = x_B - D_0/2
    
    i = 0
    
    while i < len(t):
        
        temp1 = quad(integranda, ll, ul, args=(x_B,sigma[i]))
        a[i] = N/(temp1[0])
        
        temp2 = quad(integrandADlinear, ll, ul, args=(a[i],x_B,sigma[i]))
        AD[i] = temp2[0]
        
#        Nt = quad(integrandN, ll, ul, args=(a[i],x_B,sigma[i]))
#        print(Nt[0])
        
        i = i + 1
    
    title = "Sigma vs. Time"
    ylabel = "$\sigma$ (m)"
    save_title = "GF_Linear_Sigma"
    stitle = "Gaussian Distribution: Linear Sigma"
    
    if noplot==0:
        plotSingle(AD,sigma,ylabel,title,save_title,stitle)
    
    return AD
    
def combinedLinear():
    
    AD = np.zeros(len(t))        
    a = np.zeros(len(t))

    x_B = x_0 + A_x*np.cos(w*t + phi_x)
    sigma = sigma_0 + A_sigma*np.cos(w*t + phi_sigma)
    
    ul = x_B + D_0/2
    ll = x_B - D_0/2
    
    i = 0
    
#    nn = 0 # For verification section just below
    
    while i < len(t):
        
        temp1 = quad(integranda, ll[[i]], ul[[i]], args=(x_B[i],sigma[i]))
        a[i] = N/temp1[0]
        
        # This section verifies the distribution is valid (Not updated for Gaussian)
#        tn = nn/100
#        print(tn)
        
#        if tn.is_integer():
#            x = np.linspace(ll[[i]],ul[[i]],1000)
#            plt.figure(nn*10)
#            plt.plot(x,a[i]*(-(x-x_B[i])**2 + (sigma[i]/2)**2))
#            plt.xlim(ll[[i]],ul[[i]])
#        
#        nn = nn + 1
        ################# End Verification #################
        
        temp2 = quad(integrandADlinear, ll[[i]], ul[[i]], args=(a[i],x_B[i],sigma[i]))
        AD[i] = temp2[0]
        
#        Nt = quad(integrandN, ll[[i]], ul[[i]], args=(a[i],x_B[i],sigma[i]))
#        print(Nt[0])
        
        i = i + 1
        
    datax = xlinear(1)
    datasig = sigmalinear(1)
    
    save_title = "GF_Linear_Combined"
    stitle = "Gaussian Distribution: Linear Combined"
    plotCombined(AD, save_title, datax, datasig, stitle)
    
def xquad(noplot):
    
    AD = np.zeros(len(t))    

    x_B = x_0 + A_x*np.cos(w*t + phi_x)
    
    sigma = None
    sigma = sigma_0
    
    ul = x_B + D_0/2
    ll = x_B - D_0/2
    
    a = np.zeros(len(t))
    
    i = 0
    
    while i < len(t):
        
        temp1 = quad(integranda, ll[[i]], ul[[i]], args=(x_B[i],sigma))
        a[i] = Ne/temp1[0]
        
        temp2 = quad(integrandADquad, ll[[i]], ul[[i]], args=(a[i],x_B[i],sigma))
        AD[i] = temp2[0]
        
#        Nt = quad(integrandN, ll[[i]], ul[[i]], args=(a[i],x_B[i],sigma))
#        print(Nt[0])
        
        i = i + 1
    
    title = "Postion vs. Time"
    ylabel = "$x_m$ (m)"
    save_title = "GF_Quad_Position"
    stitle = "Gaussian Distribution: Quadratic Position"
    
    if noplot==0:
        plotSingle(AD,x_B,ylabel,title,save_title,stitle)
    
    return AD
    
def sigmaquad(noplot):
    
    AD = np.zeros(len(t))    

    sigma = sigma_0 + A_sigma*np.cos(w*t + phi_sigma)
    
    x_B = None
    x_B = x_0
    
    ul = x_B + D_0/2
    ll = x_B - D_0/2
    
    a = np.zeros(len(t))
    
    i = 0  
    
    while i < len(t):
        
        temp1 = quad(integranda, ll, ul, args=(x_B,sigma[i]), epsabs=1.49e-15, epsrel=1.49e-15)
        a[i] = Ne/temp1[0]
        
        temp2 = quad(integrandADquad, ll, ul, args=(a[i],x_B,sigma[i]), epsabs=1.49e-15, epsrel=1.49e-15)
        AD[i] = temp2[0]
        
#        Nt = quad(integrandN, ll, ul, args=(a[i],x_B,sigma[i]))
#        print(Nt[0])
        
        i = i + 1
    
    title = "Sigma vs. Time"
    ylabel = "$\sigma$ (m)"
    save_title = "GF_Quad_Sigma"
    stitle = "Gaussian Distribution: Quadratic Sigma"
    
    if noplot==0:
        plotSingle(AD,sigma,ylabel,title,save_title,stitle)
    
    return AD
    
def combinedQuad():
    
    AD = np.zeros(len(t))        
    a = np.zeros(len(t))

    x_B = x_0 + A_x*np.cos(w*t + phi_x)
    sigma = sigma_0 + A_sigma*np.cos(w*t + phi_sigma)
    
    ll = x_B - D_0/2
    ul = x_B + D_0/2
    
    i = 0
    
    while i < len(t):
        
        temp = quad(integranda, ll[[i]], ul[[i]], args=(x_B[i],sigma[i]), epsabs=1.49e-15, epsrel=1.49e-15)
        a[i] = Ne/temp[0]
        
        temp = quad(integrandADquad, ll[[i]], ul[[i]], args=(a[i],x_B[i],sigma[i]), epsabs=1.49e-15, epsrel=1.49e-15)
        AD[i] = temp[0]
        
#        Nt = quad(integrandN, ll[[i]], ul[[i]], args=(a[i],x_B[i],D[i]))
#        print(Nt[0])
        
        i = i + 1
        
    datax = xquad(noplot = 1)
    datasig = sigmaquad(noplot = 1)
        
    save_title = "GF_Quad_Combined"
    stitle = "Gaussian Distribution: Quadratic Combined"
    plotCombined(AD,save_title,datax,datasig,stitle)
    
# Finds the normalizing constant for the distribution function to ensure the
# total number of particles remains the same
def integranda(x, x_B, sigma):
    return (np.sqrt(2*np.pi*sigma**2)**(-1))*np.exp(-((x-x_B)**2)/(2*sigma**2))
    
# Finds the total number of particles (used for verification only)
def integrandN(x, a, x_B, sigma):
    return a*(np.sqrt(2*np.pi*sigma**2)**(-1))*np.exp(-((x-x_B)**2)/(2*sigma**2))
    
# Used to find the acceptance for the linear term
def integrandADlinear(x, a, x_B, sigma):
    return a*(np.sqrt(2*np.pi*sigma**2)**(-1))*np.exp(-((x-x_B)**2)/(2*sigma**2))*(k_0*x)
    
# Used to find the acceptancs for the quadratic term
def integrandADquad(x, a, x_B, sigma):
    return (a*(np.sqrt(2*np.pi*sigma**2)**(-1))*np.exp(-((x-x_B)**2)/(2*sigma**2)))*(k_1*x**2)
    
def plotSingle(AD,data,ylabel,title,save_title,stitle):

    n = 1

    fig = plt.figure(n)
    st = fig.suptitle("%s"%stitle, fontsize="x-large")
    n = n + 1
    
    plt.subplot(2,1,1)
    plt.plot(t,data)
    plt.ylabel("%s"%(ylabel))
    plt.title("%s"%(title))
    
    plt.subplot(2,1,2)
    plt.plot(t,AD)
    plt.ylabel("$N_D$")
    plt.xlabel("Time")
    plt.title("Particles Detected vs. Time")
    plt.tight_layout()
    
    st.set_y(1.0)
    fig.subplots_adjust(top=0.85)
    
    plt.savefig('Plots/%s.png'%(save_title), bbox_inches='tight', dpi=300)

    fig = plt.figure(n)
    st = fig.suptitle("%s"%stitle, fontsize="x-large")
    n = n + 1
    
    Af = fft(AD)
    plt.plot(tf[1:], 2/N * np.abs(Af[:N/2])[1:])
    plt.xlim(0,5)
    plt.xlabel("Frequency")
    plt.title("DFT of Particles Detected")
    
    st.set_y(1.0)
    fig.subplots_adjust(top=0.85)
    
    plt.savefig('Plots/%s-DFT.png'%(save_title), bbox_inches='tight', dpi=300)
    
def plotCombined(AD,save_title,datax,datasig,stitle):

    n = 1

    fig = plt.figure(n)
    st = fig.suptitle("%s"%stitle, fontsize="x-large")
    n = n + 1
    
    plt.subplot(3,1,1)
    plt.plot(t,datax)
    plt.ylabel("$N_D$")
    plt.title("Particles Detected from changing $x_m$")
    
    plt.subplot(3,1,2)
    plt.plot(t,datasig)
    plt.ylabel("$N_D$")
    plt.title("Particles Detected from changing $\sigma$")
    
    plt.subplot(3,1,3)
    plt.plot(t,AD)
    plt.ylabel("$N_D$")
    plt.xlabel("Time")
    plt.title("Particles Detected vs. Time")
    plt.tight_layout()
    
    st.set_y(1.0)
    fig.subplots_adjust(top=0.85)
    
    plt.savefig('Plots/%s.png'%(save_title), bbox_inches='tight', dpi=300)

    fig = plt.figure(n)
    st = fig.suptitle("%s"%stitle, fontsize="x-large")
    n = n + 1
    
    Af = fft(AD)
    plt.plot(tf[1:], 2/N * np.abs(Af[:N/2])[1:])
    plt.xlim(0,5)
#    plt.ylim(0,2)
    plt.xlabel("Frequency")
    plt.title("DFT of Particles Detected")
    
    st.set_y(1.0)
    fig.subplots_adjust(top=0.85)
    
    plt.savefig('Plots/%s-DFT.png'%(save_title), bbox_inches='tight', dpi=300)
    
main()
