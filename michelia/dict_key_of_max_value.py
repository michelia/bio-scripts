
def dict_key_of_max_value(adict):
    '''
    Function: Return the key of a dict when the key's value is max
    >>>adict = {'a': 1, 'b': 2, 'c': 3}
    >>>dict_key_of_max_value(adict)
    >>>'a'
    '''
    # reverse adict and the max key's value
    return dict(zip(adict.itervalues(), adict.iterkeys()))[max(adict.itervalues())]

# Function is equal to the follow

# def dict_key_of_max_value(adict):
#     k = ''
#     v = 0
#     for key, value in adict.iteritems():
#         if value > v:
#             v = value
#             k = key
#     return k
