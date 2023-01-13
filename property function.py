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
    

aa = a('Prandtl', 300, 12000, 'water')