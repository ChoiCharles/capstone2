import pandas as pd

df = pd.read_csv('C:/Users/GP62/Downloads/4-1/캡스톤/Water_density.csv')

p = []
i = 0
while i < 200:
    i += 1
    P = i * 10 * 1000
    p.append(P)

df.index = p