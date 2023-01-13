from kivy.app import App
from kivy.uix.widget import Widget
from kivy.properties import ObjectProperty
from kivy.lang import Builder
from kivy.uix.screenmanager import ScreenManager
from kivy.uix.screenmanager import Screen
import numpy as np
import csv 

def a(Pro, T, P, fluid):
    data=list()
    
    if Pro == 'Density':
        if fluid == 'water':
            f = open('Water_Density.csv')
        elif fluid == 'air':
            f = open('Air_Density.csv')
    elif Pro == 'Viscosity':
        if fluid == 'water':
            f = open('Water_Viscosity.csv')
        elif fluid == 'air':
            f = open('Air_Viscosity.csv')
    elif Pro == 'beta':
        if fluid == 'water':
            f = open('Water_Expansion.csv')
        elif fluid == 'air':
            f = open('Air_Expansion.csv')
    elif Pro == 'Prandtl':
        if fluid == 'water':
            f = open('Water_Prandtl.csv')
        elif fluid == 'air':
            f = open('Air_Prandtl.csv')
    elif Pro == 'cp':
        if fluid == 'water':
            f = open('Water_Cp.csv')
        elif fluid == 'air':
            f = open('Air_Cp.csv')
    elif Pro == 'conductivity':
        if fluid == 'water':
            f = open('Water_Conductivity.csv')
        elif fluid == 'air':
            f = open('Air_Conductivity.csv')        
    elif Pro == 'enthalpy':
        if fluid == 'water':
            f = open('Water_Enthalpy.csv')
        elif fluid == 'air':
            f = open('Air_Enthalpy.csv')
    
    rea = csv.reader(f)
    pro_=[]
    
    for row in rea:
        k=(T-270)/5
        data.append(row[int(k)])
        data.append(row[int(k+1)])
        a=row[int(k)]
        b=row[int(k+1)]
        A = float(a)
        B = float(b)
        pro = A+((B-A)*((k-int(k)/(int(k+1)-int(k)))))
        pro_.append(pro)
    p = int(P/10000)
    c=pro_[p]
    d=pro_[p+1]
    C = float(c)
    D = float(d)
    pro =C+((D-C)*((p-int(p)/(int(p+1)-int(p)))))
    f.close

    return pro


class MainWindow(Screen):
    pass

class Properties(Screen):
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
            
    def result(self, instance):
        T = self.ids.T_input.text
        P = self.ids.P_input.text
        T = float(T) + 273
        P = float(P) * 1000
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        
        try:
            k = round(a('conductivity', T, P, fluid), 3)
            mu = round(a('Viscosity', T, P, fluid), 5)
            rho = round(a('Density', T, P, fluid), 3)
            cp = round(a('cp', T, P, fluid), 3)
            beta = round(a('beta', T, P, fluid), 5)
            enthalpy = round(a('enthalpy', T, P, fluid), 3)/1000
            self.ids.k.text = f'Conductivity [W/mK] \n = {str(k)}'
            self.ids.mu.text = f'Viscosity [m^2/s] \n = {str(mu)}'
            self.ids.rho.text = f'Density [kg/m^3] \n = {str(rho)}'
            self.ids.cp.text = f'Specific Heat [kJ/kgK] \n = {str(cp)}'
            self.ids.beta.text = f'Coefficient of Expansion [1/K] \n = {str(beta)}'
            self.ids.enthalpy.text = f'Enthalpy [kJ/kg] \n = {str(enthalpy)}'
            
        except:
            self.ids.k.text = "Error"
            self.ids.mu.text = "Error"
            self.ids.rho.text = "Error"
            self.ids.cp.text = "Error"
            self.ids.beta.text = "Error"
            self.ids.enthalpy.text = "Error"
            

class Demensionless(Screen):
    pass

class Reynolds(Screen):
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
            
    def result_clicked(self, value):
            self.ids.result.text = f'{value}'
                
    def result(self,instance):
        T = self.ids.T_input.text
        P = self.ids.P_input.text
        V = self.ids.V_input.text
        L = self.ids.L_input.text
        T = float(T) + 273
        P = float(P) * 1000
        V = float(V)
        L = float(L)
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        mu = a('Viscosity', T, P, fluid)
        rho = a('Density', T, P, fluid)
        try:
            answer = rho*V*L/mu
            if 1000 <= answer <= 1000000:
                answer = round(answer/1000, 2)
                self.ids.result.text = f'Reynolds Number =\n{str(answer)} X 10\u00b3'
            elif 1000000 <= answer <= 1000000000:
                answer = round(answer/1000000, 2)
                self.ids.result.text = f'Reynolds Number =\n{str(answer)} X 10\u2076'
            
        except:
            self.ids.result.text = "Error"
            
class Prandtl(Screen):
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
            
    def result(self,instance):
        T = self.ids.T_input.text
        P = self.ids.P_input.text
        T = float(T) + 273
        P = float(P) * 1000
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        
        try:
            Pr = a('Prandtl', T, P, fluid)
            self.ids.result.text = f'Prandtl Number =\n{str(Pr)}'
            
        except:
            self.ids.result.text = "Error"

class Grashof(Screen):
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
                
    def result(self,instance):
        Ti = self.ids.Ti_input.text
        Ts = self.ids.Ts_input.text
        P = self.ids.P_input.text
        L = self.ids.L_input.text
        Ti = float(Ti) + 273
        Ts = float(Ts) + 273
        Tm = (Ti + Ts)/2
        P = float(P) * 1000
        L = float(L)
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        g = 9.81
        beta = a('beta', Tm, P, fluid)
        nu =  a('Viscosity', Tm, P, fluid)/a('Density', Tm, P, fluid)
        try:
            answer = g * beta * (Ts - Ti) * np.power(L,3) / np.power(nu, 2)
            if answer <= 1000:
                self.ids.result.text = f'Grashof Number =\n{str(answer)}'
            elif 1000 <= answer <= 1000000:
                answer = round(answer/1000, 2)
                self.ids.result.text = f'Grashof Number =\n{str(answer)} X 10\u00b3'
            elif 1000000 <= answer <= 1000000000:
                answer = round(answer/1000000, 2)
                self.ids.result.text = f'Grashof Number =\n{str(answer)} X 10\u2076'
            else:
                answer = round(answer/1000000000, 2)
                self.ids.result.text = f'Grashof Number =\n{str(answer)} X 10\u2079'
            
        except:
            self.ids.result.text = "Error"

class Rayleigh(Screen):
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
                
    def result(self,instance):
        Ti = self.ids.Ti_input.text
        Ts = self.ids.Ts_input.text
        P = self.ids.P_input.text
        L = self.ids.L_input.text
        Ti = float(Ti) + 273
        Ts = float(Ts) + 273
        Tm = (Ti + Ts)/2
        P = float(P) * 1000
        L = float(L)
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        g = 9.81
        beta = a('beta', Tm, P, fluid)
        nu =  a('Viscosity', Tm, P, fluid)/a('Density', Tm, P, fluid)
        Gr = g * beta * (Ts - Ti) * np.power(L,3) / np.power(nu, 2)
        Pr = a('Prandtl', Tm, P, fluid)
       
        try:
            Ra = Gr*Pr
            
            if Ra < 1000:
                self.ids.result.text = f'Rayleigh Number =\n{str(Ra)}'
            elif 1000 <= Ra <1000000:
                Ra = round(Ra/1000, 2)
                self.ids.result.text = f'Rayleigh Number =\n{str(Ra)} X 10\u00b3'
            elif 1000000 <= Ra <1000000000:
                Ra = round(Ra/1000000, 2)
                self.ids.result.text = f'Rayleigh Number =\n{str(Ra)} X 10\u2076'
            else:
                Ra = round(Ra/1000000000, 2)
                self.ids.result.text = f'Rayleigh Number =\n{str(Ra)} X 10\u2079'
            
        except:
            self.ids.result.text = "Error"
            
class Convection(Screen):
    pass

class Forced_Convection(Screen):
    pass

class External_Flow(Screen):
    pass

class External_Plate(Screen):
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
                
    def result(self,instance):
        Ti = self.ids.Ti_input.text
        Ts = self.ids.Ts_input.text
        P = self.ids.P_input.text
        L = self.ids.L_input.text
        V = self.ids.V_input.text
        Ti = float(Ti) + 273
        Ts = float(Ts) + 273
        Tm = (Ti + Ts)/2
        P = float(P) * 1000
        L = float(L)
        V = float(V)
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        Pr = a('Prandtl', Tm, P, fluid)
        mu = a('Viscosity', Tm, P, fluid)
        rho = a('Density', Tm, P, fluid)
        Re = rho*V*L/mu
        Re_c=5e5
        k = a('conductivity', Tm, P, fluid)
        if Re < Re_c:
            Nu_x = 0.3387*Re**(1/2)*Pr**(1/3)/((1+(0.0468/Pr)**(2/3))**(1/4))
            # Churchill and Ozoe, valid Pe = Re*Pr > 100
            Nu = 2*Nu_x
        else:
            A = 0.037*Re_c**(4/5) - 0.664*Re_c**(1/2)
            Nu = (0.037*Re**(4/5) - A)*Pr**(1/3) 
        try:
            if Nu < 1000:
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/L, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error"
            
class External_Cylinder(Screen):
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
                
    def result(self,instance):
        Ti = self.ids.Ti_input.text
        Ts = self.ids.Ts_input.text
        P = self.ids.P_input.text
        D = self.ids.D_input.text
        V = self.ids.V_input.text
        Ti = float(Ti) + 273
        Ts = float(Ts) + 273
        Tm = (Ti + Ts)/2
        P = float(P) * 1000
        D = float(D)
        V = float(V)
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        Pr = a('Prandtl', Tm, P, fluid)
        mu = a('Viscosity', Tm, P, fluid)
        rho = a('Density', Tm, P, fluid)
        Re = rho*V*D/mu
        k = a('conductivity', Tm, P, fluid)
        d1 = 0.62*Re**(1/2)*Pr**(1/3)
        d2 = (1 + (0.4/Pr)**(2/3))**(1/4)
        d3 = (1 + (Re/282000)**(5/8))**(4/5)
        Nu = 0.3 + d1/d2*d3
       
        try:
            if Nu < 1000:
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/D, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error"

class External_Sphere(Screen):
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
                
    def result(self,instance):
        Ti = self.ids.Ti_input.text
        Ts = self.ids.Ts_input.text
        P = self.ids.P_input.text
        D = self.ids.D_input.text
        V = self.ids.V_input.text
        Ti = float(Ti) + 273
        Ts = float(Ts) + 273
        Tm = (Ti + Ts)/2
        P = float(P) * 1000
        D = float(D)
        V = float(V)
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        Pr = a('Prandtl', Tm, P, fluid)
        mu = a('Viscosity', Ti, P, fluid)
        mus = a('Viscosity', Ts, P, fluid)
        rho = a('Density', Tm, P, fluid)
        Re = rho*V*D/mu
        k = a('conductivity', Tm, P, fluid)
        
        mu_ratio=mu/mus
        Nu = 2 + (0.4*Re**(1/2) + 0.06*Re**(2/3))*Pr**(0.4)*mu_ratio**(1/4)
       
        try:
            if Nu < 1000:
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/D, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error"
            
class Infinite_Cylinder(Screen):
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
                
    def result(self,instance):
        Ti = self.ids.Ti_input.text
        Ts = self.ids.Ts_input.text
        P = self.ids.P_input.text
        D = self.ids.D_input.text
        V = self.ids.V_input.text
        Ti = float(Ti) + 273
        Ts = float(Ts) + 273
        Tm = (Ti + Ts)/2
        P = float(P) * 1000
        D = float(D)
        V = float(V)
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        Pr = a('Prandtl', Tm, P, fluid)
        mu = a('Viscosity', Tm, P, fluid)
        rho = a('Density', Tm, P, fluid)
        Re = rho*V*D/mu
        k = a('conductivity', Tm, P, fluid)
        
        if 0.4 < Re:
            C = 0.989
            m = 0.330
        elif 4.0 < Re:
            C = 0.911
            m = 0.385
        elif 40 < Re:
            C = 0.683
            m = 0.466
        elif 4000 < Re:
            C = 0.193
            m = 0.618
        elif 40000 < Re < 400000:
            C = 0.027
            m = 0.805
            
        Nu = C*Re**m*Pr**(1/3)
       
        try:
            if Nu < 1000:
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/D, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error"

class Internal_Flow(Screen):
    pass

class Entrance_Area(Screen):
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
                
    def result(self,instance):
        T = self.ids.T_input.text
        x = self.ids.x_input.text
        P = self.ids.P_input.text
        D = self.ids.D_input.text
        V = self.ids.V_input.text
        T = float(T) + 273
        x = float(x)
        P = float(P) * 1000
        D = float(D)
        V = float(V)
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        Pr = a('Prandtl', T, P, fluid)
        mu = a('Viscosity', T, P, fluid)
        rho = a('Density', T, P, fluid)
        Re = rho*V*D/mu
        k = a('conductivity', T, P, fluid)
        Gz = (D/x)*Re*Pr
            
        Nu = 3.66 + 0.0668*Gz/(1 + 0.04*Gz**(2/3))
        try:
            if Nu < 1000:
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/D, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error"
            
class Combined_Entrance_Area(Screen):
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
                
    def result(self,instance):
        T = self.ids.T_input.text
        x = self.ids.x_input.text
        P = self.ids.P_input.text
        D = self.ids.D_input.text
        V = self.ids.V_input.text
        T = float(T) + 273
        x = float(x)
        P = float(P) * 1000
        D = float(D)
        V = float(V)
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        Pr = a('Prandtl', T, P, fluid)
        mu = a('Viscosity', T, P, fluid)
        rho = a('Density', T, P, fluid)
        Re = rho*V*D/mu
        k = a('conductivity', T, P, fluid)
        Gz = (D/x)*Re*Pr
        a1 = 3.66/np.tanh(2.264*Gz^(-1/3) + 1.7*Gz^(-2/3))
        a2 = 0.0499*Gz*np.tanh(1/Gz)
        a3 = np.tanh(2.432*Pr^(1/6)*Gz^(-1/6))
        Nu = (a1 + a2)/a3
        try:
            if Nu < 1000:
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/D, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error"

class Circular_Tube_Turbulent_Flow(Screen):
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
                
    def result(self,instance):
        Ti = self.ids.Ti_input.text
        Ts = self.ids.Ts_input.text
        P = self.ids.P_input.text
        D = self.ids.D_input.text
        V = self.ids.V_input.text
        Ti = float(Ti) + 273
        Ts = float(Ts) + 273
        Tm = (Ti + Ts)/2
        P = float(P) * 1000
        D = float(D)
        V = float(V)
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        Pr = a('Prandtl', Tm, P, fluid)
        mu = a('Viscosity', Tm, P, fluid)
        mus = a('Viscosity', Ts, P, fluid)
        rho = a('Density', Tm, P, fluid)
        Re = rho*V*D/mu
        k = a('conductivity', Tm, P, fluid)
        n = 0.4 if Ts > Tm else 0.3
        mu_over_mus=mu/mus
        
        if 0.6 <= Pr <= 160:
            Nu = 0.023*Re**(4/5)*Pr**n
        elif Pr <= 16700:
            Nu = 0.027*Re**(4/5)*Pr**(1/3)*mu_over_mus**(0.14)
        
        try:
            if Nu < 1000:
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/D, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error"

class Free_Convection(Screen):
    pass

class Vertical_Plate(Screen):
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
                
    def result(self,instance):
        Ti = self.ids.Ti_input.text
        Ts = self.ids.Ts_input.text
        P = self.ids.P_input.text
        L = self.ids.L_input.text
        Ti = float(Ti) + 273
        Ts = float(Ts) + 273
        Tm = (Ti + Ts)/2
        P = float(P) * 1000
        L = float(L)
        g = 9.81
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        Pr = a('Prandtl', Tm, P, fluid)
        mu = a('Viscosity', Tm, P, fluid)
        rho = a('Density', Tm, P, fluid)
        beta = a('beta', Tm, P, fluid)
        nu = mu/rho
        k = a('conductivity', Tm, P, fluid)
        Gr = g*beta*np.abs(Ti-Ts)*L**3/nu**2
        Ra = Gr*Pr
        
        if Ra < 1e9 :
           Nu = 0.68 + 0.670*Ra**(1/4)/(1+(0.492/Pr)**(9/16))**(4/9)
        else :
           Nu = (0.825 + 0.387*Ra**(1/6)/(1 + (0.492/Pr)**(9/16))**(8/27))**2
        
        try:
            if Nu < 1000:
                self.ids.result.text=f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text=f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/L, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error"
            

class Horizontal_Plate(Screen):
    pass

class Horizontal_Plate_Top(Screen):
    
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
                
    def result(self,instance):
        Ti = self.ids.Ti_input.text
        Ts = self.ids.Ts_input.text
        P = self.ids.P_input.text
        L = self.ids.L_input.text
        Ti = float(Ti) + 273
        Ts = float(Ts) + 273
        Tm = (Ti + Ts)/2
        P = float(P) * 1000
        L = float(L)
        g = 9.81
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        Pr = a('Prandtl', Tm, P, fluid)
        mu = a('Viscosity', Tm, P, fluid)
        rho = a('Density', Tm, P, fluid)
        beta = a('beta', Tm, P, fluid)
        nu = mu/rho
        k = a('conductivity', Tm, P, fluid)
        Gr = g*beta*np.abs(Ti-Ts)*L**3/nu**2
        Ra = Gr*Pr
        
        if 1e4 < Ra < 1e7 and Pr >= 0.7:
            Nu = 0.54*Ra**(1/4)
        elif 1e7 <= Ra < 1e11:
            Nu = 0.15*Ra**(1/3)
        else:
            Nu = '계산 불가'
        try:
            if Nu < 1000:
                self.ids.result.text=f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text=f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/L, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error"            

                                   
class Horizontal_Plate_Bottom(Screen):
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
                
    def result(self,instance):
        Ti = self.ids.Ti_input.text
        Ts = self.ids.Ts_input.text
        P = self.ids.P_input.text
        L = self.ids.L_input.text
        Ti = float(Ti) + 273
        Ts = float(Ts) + 273
        Tm = (Ti + Ts)/2
        P = float(P) * 1000
        L = float(L)
        g = 9.81
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        Pr = a('Prandtl', Tm, P, fluid)
        mu = a('Viscosity', Tm, P, fluid)
        rho = a('Density', Tm, P, fluid)
        beta = a('beta', Tm, P, fluid)
        nu = mu/rho
        k = a('conductivity', Tm, P, fluid)
        Gr = g*beta*np.abs(Ti-Ts)*L**3/nu**2
        Ra = Gr*Pr
        Nu = 0.52*Ra**(1/5)
        try:
            if Nu < 1000:
                self.ids.result.text=f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text=f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/L, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error"

class Horizontal_Cylinder(Screen):
    
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
                
    def result(self,instance):
        Ti = self.ids.Ti_input.text
        Ts = self.ids.Ts_input.text
        P = self.ids.P_input.text
        L = self.ids.L_input.text
        Ti = float(Ti) + 273
        Ts = float(Ts) + 273
        Tm = (Ti + Ts)/2
        P = float(P) * 1000
        L = float(L)
        g = 9.81
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        Pr = a('Prandtl', Tm, P, fluid)
        mu = a('Viscosity', Tm, P, fluid)
        rho = a('Density', Tm, P, fluid)
        beta = a('beta', Tm, P, fluid)
        nu = mu/rho
        k = a('conductivity', Tm, P, fluid)
        Gr = g*beta*np.abs(Ti-Ts)*L**3/nu**2
        Ra = Gr*Pr
        Nu = (0.6 + (0.387*Ra**(1/6))/(1+(0.559/Pr)**(9/16))**(8/27))**2
        try:
            if Nu < 1000:
                self.ids.result.text=f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text=f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/L, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error" 

class Sphere(Screen):
    
    def spinner_clicked(self, value):
            self.ids.click_label.text = f'Fluid: {value}'
            self.ids.spinner_id.text = f'{value}'
                
    def result(self,instance):
        Ti = self.ids.Ti_input.text
        Ts = self.ids.Ts_input.text
        P = self.ids.P_input.text
        L = self.ids.L_input.text
        Ti = float(Ti) + 273
        Ts = float(Ts) + 273
        Tm = (Ti + Ts)/2
        P = float(P) * 1000
        L = float(L)
        g = 9.81
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        Pr = a('Prandtl', Tm, P, fluid)
        mu = a('Viscosity', Tm, P, fluid)
        rho = a('Density', Tm, P, fluid)
        beta = a('beta', Tm, P, fluid)
        nu = mu/rho
        k = a('conductivity', Tm, P, fluid)
        Gr = g*beta*np.abs(Ti-Ts)*L**3/nu**2
        Ra = Gr*Pr
        Nu = 2 + (0.589*Ra**(1/4))/(1+(0.469/Pr)**(9/16))**(4/9)
        try:
            if Nu < 1000:
                self.ids.result.text=f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text=f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/L, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error" 
        
        
class fc_horizontal_enclosure(Screen):
    
    def spinner_clicked(self, value):
        self.ids.click_label.text = f'Fluid: {value}'
        self.ids.spinner_id.text = f'{value}'
    
    def result_clicked(self, value):
        self.ids.result.text = f'{value}'
            
    def result(self,instance):
        
        T1 = self.ids.T1_input.text
        T2 = self.ids.T2_input.text
        P = self.ids.P_input.text
        L = self.ids.L_input.text
        T1 = float(T1) + 273
        T2 = float(T2) + 273
        P = float(P) * 1000
        L = float(L)
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        g = 9.81
        T = (T1 + T2)/2
        k = a('conductivity', T, P, fluid)
        mu = a('Viscosity', T, P, fluid)
        rho = a('Density', T, P, fluid)
        beta = a('beta', T, P, fluid)
        nu = mu/rho
        Pr = a('Prandtl', T, P, fluid)
        Gr = g*beta*np.abs(T1-T2)*L**3/nu**2
        Ra = Gr*Pr
        Nu = 0.069*Ra**(1/3)*Pr**(0.074)
        
        try:
            if Nu < 1000:
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/L, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error"   
    

class fc_vertical_cavity(Screen):
    
    def spinner_clicked(self, value):
        self.ids.click_label.text = f'Fluid: {value}'
        self.ids.spinner_id.text = f'{value}'
    
    def result_clicked(self, value):
        self.ids.result.text = f'{value}'
            
    def result(self,instance):
        T1 = self.ids.T1_input.text
        T2 = self.ids.T2_input.text
        P = self.ids.P_input.text
        L = self.ids.L_input.text
        T1 = float(T1) + 273
        T2 = float(T2) + 273
        P = float(P) * 1000
        L = float(L)
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        g = 9.81
        T = (T1 + T2)/2
        k = a('conductivity', T, P, fluid)
        mu = a('Viscosity', T, P, fluid)
        rho = a('Density', T, P, fluid)
        beta = a('beta', T, P, fluid)
        nu = mu/rho
        Pr = a('Prandtl', T, P, fluid)
        Gr = g*beta*np.abs(T1-T2)*L**3/nu**2
        Ra = Gr*Pr
        H = self.ids.H_input.text
        H = float(H)
        
        if 2 <= H/L <= 10 and Pr <= 10**10 and 10**3 <= Ra <= 10**10:
            Nu = 0.22 * ((Pr * Ra / (0.2 + Pr))**0.28) * ((H / L)**(-1/4))
        elif 1 <= H/L <= 2 and 10**(-3) <= Pr <= 10**5 and 10**3 <= Ra * Pr / (0.2 + Pr):
            Nu = 0.18 * ((Pr * Ra / (0.2 + Pr))**0.29)
        elif 10 <= H/L <= 40 and 1 <= Pr <= 2 * 10**4 and 10**4 <= Ra <= 10**7:
            Nu = 0.42*(Ra**(1/4))*(Pr**(0.012))*((H/L)**(-0.3))
        elif 1 <= H/L <= 40 and 1 <= Pr <= 20 and 10**6 <= Ra <= 10**9:
            Nu = 0.046 * (Ra**(1/3))
            
            
        try:
            if Nu < 1000:
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/L, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error"    



class fc_concentric_cylinder (Screen):
    
    def spinner_clicked(self, value):
        self.ids.click_label.text = f'Fluid: {value}'
        self.ids.spinner_id.text = f'{value}'
    
    def result_clicked(self, value):
            self.ids.result.text = f'{value}'
            
    def result(self,instance):
        Ti = self.ids.Ti_input.text
        To = self.ids.To_input.text
        ri = self.ids.ri_input.text
        ro = self.ids.ro_input.text
        P = self.ids.P_input.text
        
        Ti = float(Ti) + 273
        To = float(To) + 273
        ri = float(ri)
        ro = float(ro)
        P = float(P) * 1000
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        Lc = 2*(np.log(ro/ri))**(4/3)/(ri**(-3/5) + ro**(-3/5))**(5/3)
        g = 9.81
        T = (Ti + To)/2
        Pr = a('Prandtl', Ti, P, fluid)
        k = a('conductivity', Ti, P, fluid)
        mu = a('Viscosity', T, P, fluid)
        rho = a('Density', T, P, fluid)
        beta = a('beta', T, P, fluid)
        nu = mu/rho
        Gr = g*beta*np.abs(Ti-To)*Lc**3/nu**2
        T = (Ti + To)/2
        Ra = Gr*Pr
        k_eff = k*0.386*(Pr/(0.861 + Pr))**(1/4)*Ra**(1/4)
        
        if k_eff < k:
            k_eff = k
        k_eff = round(k_eff, 2)
        q = round(4*np.pi*k_eff*(Ti-To)/(1/ri-1/ro), 2)
        
        self.ids.result.text = f'q =\n{str(q)}'
        self.ids.result2.text = f'k_eff =\n{str(k_eff)}'
        self.ids.result3.text = f'Rayleigh =\n{str(Ra)}'


class fc_concentric_sphere (Screen):
    
    def spinner_clicked(self, value):
        self.ids.click_label.text = f'Fluid: {value}'
        self.ids.spinner_id.text = f'{value}'
    
    def result_clicked(self, value):
            self.ids.result.text = f'{value}'
            
    def result(self,instance):
        Ti = self.ids.Ti_input.text
        To = self.ids.To_input.text
        ri = self.ids.ri_input.text
        ro = self.ids.ro_input.text
        P = self.ids.P_input.text
        Ti = float(Ti) + 273
        To = float(To) + 273
        ri = float(ri)
        ro = float(ro)
        P = float(P) * 1000
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        Lc = (1/ri-1/ro)**(4/3)/(2**(1/3)*(ri**(-7/5)+ro**(-7/5))**(5/3))
        g = 9.81
        T = (Ti + To)/2
        Pr = a('Prandtl', Ti, P, fluid)
        k = a('conductivity', Ti, P, fluid)
        mu = a('Viscosity', T, P, fluid)
        rho = a('Density', T, P, fluid)
        beta = a('beta', T, P, fluid)  
        nu = mu/rho
        Gr = g*beta*np.abs(Ti-To)*Lc**3/nu**2
        T = (Ti + To)/2
        Ra = Gr*Pr
        k_eff = k*0.74*(Pr/(0.861 + Pr))**(1/4)*Ra**(1/4)
        
        if k_eff < k:
            k_eff = k
        k_eff = round(k_eff, 2)
        q = round(4*np.pi*k_eff*(Ti-To)/(1/ri-1/ro), 2)
        
        self.ids.result.text = f'q =\n{str(q)}'
        self.ids.result2.text = f'k_eff =\n{str(k_eff)}'
        self.ids.result3.text = f'Rayleigh =\n{str(Ra)}'



class fc_tilted_cavity(Screen):
    def spinner_clicked(self, value):
        self.ids.click_label.text = f'Fluid: {value}'
        self.ids.spinner_id.text = f'{value}'
    
    def result_clicked(self, value):
            self.ids.result.text = f'{value}'
            
    def result(self,instance):
        tau = self.ids.tau_input.text
        T1 = self.ids.T1_input.text
        T2 = self.ids.T2_input.text
        P = self.ids.P_input.text
        L = self.ids.L_input.text
        tau = float(tau)
        T1 = float(T1) + 273
        T2 = float(T2) + 273
        P = float(P) * 1000
        L = float(L)
        fluid = self.ids.spinner_id.text
        fluid = str(fluid)
        
        g = 9.81
        T = (T1 + T2)/2
        Pr = a('Prandtl', T, P, fluid)
        k = a('conductivity', T1, P, fluid)
        mu = a('Viscosity', T, P, fluid)
        rho = a('Density', T, P, fluid)
        beta = a('beta', T, P, fluid)   
        nu = mu/rho
        Gr = g*beta*np.abs(T1-T2)*L**3/nu**2
        Ra = Gr*Pr
        a1 = 1 - 1708/(np.cos(tau)*Ra)
        if a1 < 0:
            a1 = 0
        a2 = 1 - 1708*(np.sin(1.8*tau))**(1.6)/(Ra*np.cos(tau))
        a3 = (np.cos(tau)*Ra/5830)**(1/3) - 1
        if a3 < 0:
            a3 = 0
        Nu = 1 + 1.44*a1*a2 + a3
        
        try:
            if Nu < 1000:
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)}'
            elif 1000 <= Nu <1000000:
                Nu = round(Nu/1000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u00b3'
            elif 1000000 <= Nu <1000000000:
                Nu = round(Nu/1000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2076'
            else:
                Nu = round(Nu/1000000000, 2)
                self.ids.result.text = f'Nusselt Number =\n{str(Nu)} X 10\u2079'
            h = round(Nu*k/L, 2)
            self.ids.result2.text = f'h [W/m2K] =\n{str(h)}'
        except:
            self.ids.result.text = "Error"
            self.ids.result2.text = "Error"

class WindowManager(ScreenManager):
    pass


kv = Builder.load_file('main.kv')



class Calc(App):
    def build(self):
        return kv
    
if __name__ == '__main__':
    Calc().run()