from numpy.random import choice

def addSeqError(base, size=None):
    'Returns a random error for a given (DNA) base'

    if not base in ['a', 't', 'c', 'g',
                    'A', 'T', 'C', 'G']:
        raise ValueError('base (%s) should be a, t, c, g, A, T, C or G' % base)

    base2error = {'a':['t', 'c', 'g'],
                  't':['a', 'c', 'g'],
                  'c':['a', 't', 'g'],
                  'g':['a', 't', 'c'],
                  'A':['T', 'C', 'G'],
                  'T':['A', 'C', 'G'],
                  'C':['A', 'T', 'G'],
                  'G':['A', 'T', 'C']}

    error = choice(base2error[base], size=size)

    return(error)
