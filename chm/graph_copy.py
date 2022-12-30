import matplotlib.pyplot as plt
import os.path
from sys import exit 
from sys import argv
from decimal import *

with open("file.csv") as res:
    all_result = [row.strip() for row in res]

name = "Метод Рунге_Кутта p = 3"

result = []
for i, row in enumerate(all_result):
    if i > 1:
        result.append(row.split(","))

result.pop()

x = []
In = []
In2 = []
I = []
S = []
err = []

getcontext().prec = 18
for i, row in enumerate(result):
    print(row[0])
    x.append(Decimal(row[2]))
    In.append(Decimal(row[3]))
    In2.append(Decimal(row[4]))
    I.append(Decimal(row[7]))
    S.append(Decimal(row[5]))
    err.append(Decimal(row[8]))


plt.close()
fig = plt.figure(num = name, figsize=(12, 7))
(x_In_In2_I, x_S_err) = fig.subplots(2, 1)
fig.suptitle(name, fontsize=16, fontweight='bold')

x_In_In2_I.plot(x, In, c = 'blue', label = 'Численное решение')
x_In_In2_I.scatter(x, In, c = 'blue')
x_In_In2_I.set_xlabel('x')
#x_In_In2_I.plot(x, In2, c = 'purple', label = '')
#x_In_In2_I.scatter(x, In2, c = 'purple')
x_In_In2_I.plot(x, I, c = 'purple', label = 'Точное решение')
x_In_In2_I.scatter(x, I, c = 'purple')

x_S_err.plot(x, S, c = 'orange', label = 'S')
x_S_err.scatter(x, S, c = 'orange')
x_S_err.set_xlabel('x')
x_S_err.plot(x, err, c = 'yellow', label = 'err')
x_S_err.scatter(x, err, c = 'yellow')

x_In_In2_I.legend(fontsize = 10)
x_S_err.legend(fontsize = 10)

plt.show()
