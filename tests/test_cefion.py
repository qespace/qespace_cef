import qespace_cef

if __name__=='__main__':
    cfp = qespace_cef.CEFpars(dict(B21=1.203, B40=-0.001, B5m3=0.244), 'meV')
    print(cfp)
    cefion = qespace_cef.CEFion(ion=qespace_cef.Ion('Ce'),cfp=cfp, Hfield=[0,0,0])
    print(cefion.__str__(evecs=True))
