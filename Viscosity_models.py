import numpy as np
import scipy
from lmfit import *

def VFT(T, log_eta_inf, m, Tg):
    "Implementando a equação do modelo VFT"
    
    B = (Tg * (12 - log_eta_inf)**2)/m
    T0 = Tg * (1 - (12 - log_eta_inf)/m)
    log_viscosidade = log_eta_inf + B / (T - T0)
    
    return log_viscosidade

def AM(T, log_eta_inf, m, Tg):
    "Implementando a equação do modelo Avramov-Milchev"
    
    a = m/(12 - log_eta_inf)
    B = Tg * (np.log(10)*(12 - log_eta_inf)) ** ((12 - log_eta_inf)/m)
    log_viscosidade = log_eta_inf + (B/T) ** a
    
    return log_viscosidade

def MYEGA(T, log_eta_inf, m, Tg):
    "Implementendo equação do modelo MYEGA"
    2
    log_viscosidade = log_eta_inf + (12 - log_eta_inf) * (Tg/T) * np.exp((((m/(12 - log_eta_inf)) - 1) * ((Tg/T)-1)))
    
    return log_viscosidade

def CW(T, log_eta_inf, m, Tg):
    "Implementendo equação do modelo Cornelissen-Waterman"
    
    a = m/(12 - log_eta_inf)
    B = (12 - log_eta_inf) * (Tg **(m/(12 - log_eta_inf)))
    log_viscosidade = log_eta_inf + B/(T**a)
    
    return log_viscosidade

def F2(T, log_eta_inf, m, Tg):
    "Equação do modelo de Fulcher: expressão 2"
    
    B = ((Tg**2)/3) * (m + 12 - log_eta_inf)
    C = (24 - 2 * log_eta_inf - m)/(3 * Tg)
    
    log_viscosidade = log_eta_inf + B/(T**2) + C * T
    
    return log_viscosidade

def F3(T, log_eta_inf, m, Tg):
    "Equação do modelo de Fulcher: expressão 3"
    
    C = (24 - 2 * log_eta_inf - m)/(2*np.log10(Tg) + (1/np.log(10)))
    B = (12 - log_eta_inf - C * np.log10(Tg) ) * Tg**2
    
    log_viscosidade = log_eta_inf + B/(T**2) + C * np.log10(T)
    
    return log_viscosidade

def F4(T, log_eta_inf, m, Tg):
    "Equação do modelo de Fulcher: expressão 4"
    
    T0 = Tg * (1 - (24 - 2 * log_eta_inf )/m)
    
    B = 4 * Tg * ((12 - log_eta_inf)**3)/(m**2)
    
    log_viscosidade = log_eta_inf + B / ((T - T0)**2)
    
    return log_viscosidade

def CLU(T, log_eta_inf, m, Tg):
    "Implementando equação de Cukierman-Lane-Uhlmann"
    
    a = Tg/(m + (1/(2*np.log(10))))
    T0 = Tg - (12 - log_eta_inf - (np.log10(Tg)/2)) * a
    B = ((12 - log_eta_inf - (np.log10(Tg)/2))**2) * a
    log_viscosidade = log_eta_inf + (np.log10(T)/2) + B/(T - T0)
    
    return log_viscosidade

'''def BS(T, log_eta_inf, m, Tg):
    "Implementando equação de Blender e Shlesinger"
    
    a = (3*(12 - log_eta_inf))/2*m
    T0 = Tg * (1 - a)
    B = Tg * (12 - log_eta_inf) * (a ** (3/2))
    log_viscosidade = log_eta_inf + B/((T - T0)**(3/2))
    
    return log_viscosidade'''

def Bassler(T, log_eta_inf, m, Tg):
    "Implementando equação de Bässler"
    
    B = Tg * np.sqrt(12 - log_eta_inf)
    log_viscosidade = log_eta_inf + (B/T)**2
    
    return log_viscosidade

def FF(T, log_eta_inf, m, Tg):
    "Modelo de Fan & Fecht"
    
    B = (Tg/2)*(m + 12 - log_eta_inf)
    C = (12 - log_eta_inf - m)/(2 * Tg)
    log_viscosidade = log_eta_inf + B/T  + C * T
    
    return log_viscosidade

def Demetriou_et_al(T, log_eta_inf, m, Tg):
    "Modelo de Demetriou, Harmon, Tao, Duan, Samwer e Johnson"
    
    a = ((m/(12 - log_eta_inf)) - 1)
    Tw = Tg / a
    B = Tg * (12 - log_eta_inf) * np.exp(a)
    log_viscosidade = log_eta_inf + (B/T) * np.exp(-(T/Tw))
    
    return log_viscosidade
