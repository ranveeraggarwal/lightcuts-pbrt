import math
import random
radius = 3
center = (0, 0, 0)
samples = 200
intensity = "[10 10 10]"
for i in xrange(samples):
    u = (2 * radius * random.random()) - radius
    theta = 2 * math.pi * random.random()
    temp = math.sqrt((radius * radius) - (u * u))
    x = center[0] + temp * math.cos(theta)
    y = center[1] + temp * math.sin(theta)
    z = center[2] + u
    print "LightSource \"point\" \"rgb I\" "+intensity+" \"point from\" ["+str(x)+" "+str(y)+" "+str(z)+"]"
