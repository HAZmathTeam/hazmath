
import helpers 
import numpy

a = numpy.zeros(4) 
for i in range(len(a)): 
  a[i] = i
b = numpy.zeros(4) 
for i in range(len(a)): 
  b[i] = i*i


print ("a ", a) 
print ("b ", b) 

print ("\n\nusing dabla function ")
helpers.dabla(a, b) 

# not working 
#helpers.dabla2(helpers.dabla, a, b)

dabla = helpers.make_dabla()
print ("\n\nusing dabla function via pointr of type ", type(dabla))
helpers.dabla2(dabla, a, b)

#not working: 
#helpers.dabla2(helpers.dabla, a, b)


def dabla3(a, b): 
    print ("in python: a[2] and b[2]", a[2], b[2])
    print ("\n\nchanging stuff in python")
    a[2] = 7
    b[2] = 17
    print ("in python: a[2] and b[2]", a[2], b[2])

g = helpers.callback(a, b, dabla3)
print(g)

print ("\n\nback to C") 
helpers.dabla(a, b) 

d = helpers.dabla_struct()
d.a = a; 
d.b = b; 
d.func = dabla 

print ("\n\nthe struct thing") 
helpers.dabla4(d) 

d = helpers.dabla_struct1()
d.a = a
d.b = b 
d.type = 1
d.data = dabla3 

print ("dabla_struct1 callback")
helpers.dabla5(d) 

