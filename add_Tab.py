ls=open('3Dcalc.py','r').readlines()

f=open('3Dcalc_id.py','w')
f.write('def radarFields_3d():\n')
for l in ls:
    f.write(r"    "+l)
f.close()
