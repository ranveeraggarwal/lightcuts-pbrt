import random

#f = open('tlights', 'w')
x1 = 1
x2 = 100
y1 = 1
y2 = 100
z1 = 1
z2 = 100

for x in xrange(1,200):
    print('LightSource "point" "rgb I" [2000 2000 2000] "point from" [' + str(random.uniform(x1, x2)) + ' ' + str(random.uniform(y1, y2)) + ' ' + str(random.uniform(z1, z2)) + ']')
