import random

#f = open('tlights', 'w')
x1 = 1
x2 = 100
y1 = 1
y2 = 200
z1 = 1
z2 = 300

for x in xrange(1,1000):
    print('LightSource "point" "rgb I" [5 5 5] "point from" [' + str(random.uniform(x1, x2)) + ' ' + str(random.uniform(y1, y2)) + ' ' + str(random.uniform(z1, z2)) + ']')
