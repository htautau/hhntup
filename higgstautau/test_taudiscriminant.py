from taudiscriminant import MethodBDT, models, Error
a = models.CurrentList()
b = MethodBDT()
for name, value in a.items():
    b.addVariable(name, value, value.type, 10000000, 10000000, Error.PERCENT)
b.build('taudiscriminant/current.jet.bdt.bin')

print b.response()
a.NUMTRACK = 3
print b.response()
print b.response_errors()
print b.get_high_error()
print b.get_low_error()
