import numpy
#**** A function to convert from Mean Anomaly, M, to True Anomaly, f,
#**** given eccentricity
def convmf(m_in,ecc,lim=1e-15,fact=0.382):
    m_in=numpy.array([m_in])
    ecc=numpy.array([ecc])
    tpi = 2*numpy.pi
    m   = m_in % tpi
    k   = numpy.where(m < 0)
    m[k] = m[k] + tpi
    k   = numpy.where(m > numpy.pi)
    m[k] = tpi - m[k]
    en  = m*(1 + fact*ecc*(numpy.pi - m))
    f   = 1
    while (numpy.max(f) > lim) :
        se = numpy.sin(en)
        es = ecc*se
        f  = en - es - m
        ce = numpy.cos(en)
        ec = ecc*ce
        fp = 1 - ec
        dx = -f/fp
        dx = -f/(fp + 0.5*dx*es)
        dx = -f/(fp + 0.5*dx*es + dx*dx*ec/6)
        en = en + dx
    cen  = numpy.cos(en)
    cth  = (cen - ecc) / (1 - ecc*cen)
    th   = numpy.arccos(cth)
    th[k] = tpi - th[k]
    return th

