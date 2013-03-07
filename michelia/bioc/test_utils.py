
def testacs(instance, python_name):
    acs = getattr(instance.obj, python_name)
    if callable(acs):
        res = acs()
    else:
        res = acs

class RS4Test_Type(type):
    def __new__(mcs, name, bases, cls_dict):
        
        targetclass = cls_dict['__targetclass__']
        accessors = targetclass.__accessors__
        for rname, where, \
                python_name, as_property, \
                docstring in accessors:

            if python_name is None:
                python_name = rname
                
            cls_dict['test_' + python_name] = \
                             lambda self: testacs(self,
                                                  python_name)
            
        return type.__new__(mcs, name, bases, cls_dict)
