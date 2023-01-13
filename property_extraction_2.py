from CoolProp.CoolProp import PropsSI
import pandas as pd

p = []
y = [0]
c = [0]
i = 0
df = []

a = pd.Series(['P'], name = 'P')
p_2 = ['P']

while i < 200:
    i += 1
    P = i * 10 * 1000
    
    p.append(P)
    k = 0
    r = []
    t = []
    c.append(P)
    while k < 100:
        
        k += 1
        T = k * 5 + 270
        t.append(T)
    df = pd.DataFrame(data = t, index = t, columns = y)
d = pd.Series(p, name = 'P')
b = pd.concat([a, d], ignore_index = True)
d = b.to_frame(name = 'P')
d.index = c

i = 0
while i < 200:
    i += 1
    P = i * 10 * 1000
    
    p.append(P)
    k = 0
    r = []
    t = []
    while k < 100:
        
        k += 1
        T = k * 5 + 270
        t.append(T)
        
        rho = PropsSI('D', 'T', T, 'P', P, 'air')
        r.append(rho)
       
    dff = pd.DataFrame(data = r, index = t)
    dff.rename(columns = {0: P}, inplace=True)
    df.insert(i, P, dff)
    
df = df.transpose()
dfff = d.join(df)

dfff.to_csv("C:/Users/GP62/Downloads/4-1/캡스톤/air_density.csv", header = False, index = False)