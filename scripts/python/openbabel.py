import sys
if sys.platform.find("linux") != -1:
    import dl
    sys.setdlopenflags(sys.getdlopenflags() | dl.RTLD_GLOBAL)

# This file was created automatically by SWIG 1.3.29.
# Don't modify this file, modify the SWIG interface instead.

import _openbabel
import new
new_instancemethod = new.instancemethod
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'PySwigObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


def _swig_setattr_nondynamic_method(set):
    def set_attr(self,name,value):
        if (name == "thisown"): return self.this.own(value)
        if hasattr(self,name) or (name == "this"):
            set(self,name,value)
        else:
            raise AttributeError("You cannot add attributes to %s" % self)
    return set_attr


class PySwigIterator(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    def __init__(self): raise AttributeError, "No constructor defined"
    __repr__ = _swig_repr
    __swig_destroy__ = _openbabel.delete_PySwigIterator
    def __iter__(self): return self
PySwigIterator.value = new_instancemethod(_openbabel.PySwigIterator_value,None,PySwigIterator)
PySwigIterator.incr = new_instancemethod(_openbabel.PySwigIterator_incr,None,PySwigIterator)
PySwigIterator.decr = new_instancemethod(_openbabel.PySwigIterator_decr,None,PySwigIterator)
PySwigIterator.distance = new_instancemethod(_openbabel.PySwigIterator_distance,None,PySwigIterator)
PySwigIterator.equal = new_instancemethod(_openbabel.PySwigIterator_equal,None,PySwigIterator)
PySwigIterator.copy = new_instancemethod(_openbabel.PySwigIterator_copy,None,PySwigIterator)
PySwigIterator.next = new_instancemethod(_openbabel.PySwigIterator_next,None,PySwigIterator)
PySwigIterator.previous = new_instancemethod(_openbabel.PySwigIterator_previous,None,PySwigIterator)
PySwigIterator.advance = new_instancemethod(_openbabel.PySwigIterator_advance,None,PySwigIterator)
PySwigIterator.__eq__ = new_instancemethod(_openbabel.PySwigIterator___eq__,None,PySwigIterator)
PySwigIterator.__ne__ = new_instancemethod(_openbabel.PySwigIterator___ne__,None,PySwigIterator)
PySwigIterator.__iadd__ = new_instancemethod(_openbabel.PySwigIterator___iadd__,None,PySwigIterator)
PySwigIterator.__isub__ = new_instancemethod(_openbabel.PySwigIterator___isub__,None,PySwigIterator)
PySwigIterator.__add__ = new_instancemethod(_openbabel.PySwigIterator___add__,None,PySwigIterator)
PySwigIterator.__sub__ = new_instancemethod(_openbabel.PySwigIterator___sub__,None,PySwigIterator)
PySwigIterator_swigregister = _openbabel.PySwigIterator_swigregister
PySwigIterator_swigregister(PySwigIterator)

class vectorInt(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __iter__(self): return self.iterator()
    def __init__(self, *args): 
        _openbabel.vectorInt_swiginit(self,_openbabel.new_vectorInt(*args))
    __swig_destroy__ = _openbabel.delete_vectorInt
vectorInt.iterator = new_instancemethod(_openbabel.vectorInt_iterator,None,vectorInt)
vectorInt.__nonzero__ = new_instancemethod(_openbabel.vectorInt___nonzero__,None,vectorInt)
vectorInt.__len__ = new_instancemethod(_openbabel.vectorInt___len__,None,vectorInt)
vectorInt.pop = new_instancemethod(_openbabel.vectorInt_pop,None,vectorInt)
vectorInt.__getslice__ = new_instancemethod(_openbabel.vectorInt___getslice__,None,vectorInt)
vectorInt.__setslice__ = new_instancemethod(_openbabel.vectorInt___setslice__,None,vectorInt)
vectorInt.__delslice__ = new_instancemethod(_openbabel.vectorInt___delslice__,None,vectorInt)
vectorInt.__delitem__ = new_instancemethod(_openbabel.vectorInt___delitem__,None,vectorInt)
vectorInt.__getitem__ = new_instancemethod(_openbabel.vectorInt___getitem__,None,vectorInt)
vectorInt.__setitem__ = new_instancemethod(_openbabel.vectorInt___setitem__,None,vectorInt)
vectorInt.append = new_instancemethod(_openbabel.vectorInt_append,None,vectorInt)
vectorInt.empty = new_instancemethod(_openbabel.vectorInt_empty,None,vectorInt)
vectorInt.size = new_instancemethod(_openbabel.vectorInt_size,None,vectorInt)
vectorInt.clear = new_instancemethod(_openbabel.vectorInt_clear,None,vectorInt)
vectorInt.swap = new_instancemethod(_openbabel.vectorInt_swap,None,vectorInt)
vectorInt.get_allocator = new_instancemethod(_openbabel.vectorInt_get_allocator,None,vectorInt)
vectorInt.begin = new_instancemethod(_openbabel.vectorInt_begin,None,vectorInt)
vectorInt.end = new_instancemethod(_openbabel.vectorInt_end,None,vectorInt)
vectorInt.rbegin = new_instancemethod(_openbabel.vectorInt_rbegin,None,vectorInt)
vectorInt.rend = new_instancemethod(_openbabel.vectorInt_rend,None,vectorInt)
vectorInt.pop_back = new_instancemethod(_openbabel.vectorInt_pop_back,None,vectorInt)
vectorInt.erase = new_instancemethod(_openbabel.vectorInt_erase,None,vectorInt)
vectorInt.push_back = new_instancemethod(_openbabel.vectorInt_push_back,None,vectorInt)
vectorInt.front = new_instancemethod(_openbabel.vectorInt_front,None,vectorInt)
vectorInt.back = new_instancemethod(_openbabel.vectorInt_back,None,vectorInt)
vectorInt.assign = new_instancemethod(_openbabel.vectorInt_assign,None,vectorInt)
vectorInt.resize = new_instancemethod(_openbabel.vectorInt_resize,None,vectorInt)
vectorInt.insert = new_instancemethod(_openbabel.vectorInt_insert,None,vectorInt)
vectorInt.reserve = new_instancemethod(_openbabel.vectorInt_reserve,None,vectorInt)
vectorInt.capacity = new_instancemethod(_openbabel.vectorInt_capacity,None,vectorInt)
vectorInt_swigregister = _openbabel.vectorInt_swigregister
vectorInt_swigregister(vectorInt)

class vvInt(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __iter__(self): return self.iterator()
    def __init__(self, *args): 
        _openbabel.vvInt_swiginit(self,_openbabel.new_vvInt(*args))
    __swig_destroy__ = _openbabel.delete_vvInt
vvInt.iterator = new_instancemethod(_openbabel.vvInt_iterator,None,vvInt)
vvInt.__nonzero__ = new_instancemethod(_openbabel.vvInt___nonzero__,None,vvInt)
vvInt.__len__ = new_instancemethod(_openbabel.vvInt___len__,None,vvInt)
vvInt.pop = new_instancemethod(_openbabel.vvInt_pop,None,vvInt)
vvInt.__getslice__ = new_instancemethod(_openbabel.vvInt___getslice__,None,vvInt)
vvInt.__setslice__ = new_instancemethod(_openbabel.vvInt___setslice__,None,vvInt)
vvInt.__delslice__ = new_instancemethod(_openbabel.vvInt___delslice__,None,vvInt)
vvInt.__delitem__ = new_instancemethod(_openbabel.vvInt___delitem__,None,vvInt)
vvInt.__getitem__ = new_instancemethod(_openbabel.vvInt___getitem__,None,vvInt)
vvInt.__setitem__ = new_instancemethod(_openbabel.vvInt___setitem__,None,vvInt)
vvInt.append = new_instancemethod(_openbabel.vvInt_append,None,vvInt)
vvInt.empty = new_instancemethod(_openbabel.vvInt_empty,None,vvInt)
vvInt.size = new_instancemethod(_openbabel.vvInt_size,None,vvInt)
vvInt.clear = new_instancemethod(_openbabel.vvInt_clear,None,vvInt)
vvInt.swap = new_instancemethod(_openbabel.vvInt_swap,None,vvInt)
vvInt.get_allocator = new_instancemethod(_openbabel.vvInt_get_allocator,None,vvInt)
vvInt.begin = new_instancemethod(_openbabel.vvInt_begin,None,vvInt)
vvInt.end = new_instancemethod(_openbabel.vvInt_end,None,vvInt)
vvInt.rbegin = new_instancemethod(_openbabel.vvInt_rbegin,None,vvInt)
vvInt.rend = new_instancemethod(_openbabel.vvInt_rend,None,vvInt)
vvInt.pop_back = new_instancemethod(_openbabel.vvInt_pop_back,None,vvInt)
vvInt.erase = new_instancemethod(_openbabel.vvInt_erase,None,vvInt)
vvInt.push_back = new_instancemethod(_openbabel.vvInt_push_back,None,vvInt)
vvInt.front = new_instancemethod(_openbabel.vvInt_front,None,vvInt)
vvInt.back = new_instancemethod(_openbabel.vvInt_back,None,vvInt)
vvInt.assign = new_instancemethod(_openbabel.vvInt_assign,None,vvInt)
vvInt.resize = new_instancemethod(_openbabel.vvInt_resize,None,vvInt)
vvInt.insert = new_instancemethod(_openbabel.vvInt_insert,None,vvInt)
vvInt.reserve = new_instancemethod(_openbabel.vvInt_reserve,None,vvInt)
vvInt.capacity = new_instancemethod(_openbabel.vvInt_capacity,None,vvInt)
vvInt_swigregister = _openbabel.vvInt_swigregister
vvInt_swigregister(vvInt)

class vectorDouble(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __iter__(self): return self.iterator()
    def __init__(self, *args): 
        _openbabel.vectorDouble_swiginit(self,_openbabel.new_vectorDouble(*args))
    __swig_destroy__ = _openbabel.delete_vectorDouble
vectorDouble.iterator = new_instancemethod(_openbabel.vectorDouble_iterator,None,vectorDouble)
vectorDouble.__nonzero__ = new_instancemethod(_openbabel.vectorDouble___nonzero__,None,vectorDouble)
vectorDouble.__len__ = new_instancemethod(_openbabel.vectorDouble___len__,None,vectorDouble)
vectorDouble.pop = new_instancemethod(_openbabel.vectorDouble_pop,None,vectorDouble)
vectorDouble.__getslice__ = new_instancemethod(_openbabel.vectorDouble___getslice__,None,vectorDouble)
vectorDouble.__setslice__ = new_instancemethod(_openbabel.vectorDouble___setslice__,None,vectorDouble)
vectorDouble.__delslice__ = new_instancemethod(_openbabel.vectorDouble___delslice__,None,vectorDouble)
vectorDouble.__delitem__ = new_instancemethod(_openbabel.vectorDouble___delitem__,None,vectorDouble)
vectorDouble.__getitem__ = new_instancemethod(_openbabel.vectorDouble___getitem__,None,vectorDouble)
vectorDouble.__setitem__ = new_instancemethod(_openbabel.vectorDouble___setitem__,None,vectorDouble)
vectorDouble.append = new_instancemethod(_openbabel.vectorDouble_append,None,vectorDouble)
vectorDouble.empty = new_instancemethod(_openbabel.vectorDouble_empty,None,vectorDouble)
vectorDouble.size = new_instancemethod(_openbabel.vectorDouble_size,None,vectorDouble)
vectorDouble.clear = new_instancemethod(_openbabel.vectorDouble_clear,None,vectorDouble)
vectorDouble.swap = new_instancemethod(_openbabel.vectorDouble_swap,None,vectorDouble)
vectorDouble.get_allocator = new_instancemethod(_openbabel.vectorDouble_get_allocator,None,vectorDouble)
vectorDouble.begin = new_instancemethod(_openbabel.vectorDouble_begin,None,vectorDouble)
vectorDouble.end = new_instancemethod(_openbabel.vectorDouble_end,None,vectorDouble)
vectorDouble.rbegin = new_instancemethod(_openbabel.vectorDouble_rbegin,None,vectorDouble)
vectorDouble.rend = new_instancemethod(_openbabel.vectorDouble_rend,None,vectorDouble)
vectorDouble.pop_back = new_instancemethod(_openbabel.vectorDouble_pop_back,None,vectorDouble)
vectorDouble.erase = new_instancemethod(_openbabel.vectorDouble_erase,None,vectorDouble)
vectorDouble.push_back = new_instancemethod(_openbabel.vectorDouble_push_back,None,vectorDouble)
vectorDouble.front = new_instancemethod(_openbabel.vectorDouble_front,None,vectorDouble)
vectorDouble.back = new_instancemethod(_openbabel.vectorDouble_back,None,vectorDouble)
vectorDouble.assign = new_instancemethod(_openbabel.vectorDouble_assign,None,vectorDouble)
vectorDouble.resize = new_instancemethod(_openbabel.vectorDouble_resize,None,vectorDouble)
vectorDouble.insert = new_instancemethod(_openbabel.vectorDouble_insert,None,vectorDouble)
vectorDouble.reserve = new_instancemethod(_openbabel.vectorDouble_reserve,None,vectorDouble)
vectorDouble.capacity = new_instancemethod(_openbabel.vectorDouble_capacity,None,vectorDouble)
vectorDouble_swigregister = _openbabel.vectorDouble_swigregister
vectorDouble_swigregister(vectorDouble)

class vVector3(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __iter__(self): return self.iterator()
    def __init__(self, *args): 
        _openbabel.vVector3_swiginit(self,_openbabel.new_vVector3(*args))
    __swig_destroy__ = _openbabel.delete_vVector3
vVector3.iterator = new_instancemethod(_openbabel.vVector3_iterator,None,vVector3)
vVector3.__nonzero__ = new_instancemethod(_openbabel.vVector3___nonzero__,None,vVector3)
vVector3.__len__ = new_instancemethod(_openbabel.vVector3___len__,None,vVector3)
vVector3.pop = new_instancemethod(_openbabel.vVector3_pop,None,vVector3)
vVector3.__getslice__ = new_instancemethod(_openbabel.vVector3___getslice__,None,vVector3)
vVector3.__setslice__ = new_instancemethod(_openbabel.vVector3___setslice__,None,vVector3)
vVector3.__delslice__ = new_instancemethod(_openbabel.vVector3___delslice__,None,vVector3)
vVector3.__delitem__ = new_instancemethod(_openbabel.vVector3___delitem__,None,vVector3)
vVector3.__getitem__ = new_instancemethod(_openbabel.vVector3___getitem__,None,vVector3)
vVector3.__setitem__ = new_instancemethod(_openbabel.vVector3___setitem__,None,vVector3)
vVector3.append = new_instancemethod(_openbabel.vVector3_append,None,vVector3)
vVector3.empty = new_instancemethod(_openbabel.vVector3_empty,None,vVector3)
vVector3.size = new_instancemethod(_openbabel.vVector3_size,None,vVector3)
vVector3.clear = new_instancemethod(_openbabel.vVector3_clear,None,vVector3)
vVector3.swap = new_instancemethod(_openbabel.vVector3_swap,None,vVector3)
vVector3.get_allocator = new_instancemethod(_openbabel.vVector3_get_allocator,None,vVector3)
vVector3.begin = new_instancemethod(_openbabel.vVector3_begin,None,vVector3)
vVector3.end = new_instancemethod(_openbabel.vVector3_end,None,vVector3)
vVector3.rbegin = new_instancemethod(_openbabel.vVector3_rbegin,None,vVector3)
vVector3.rend = new_instancemethod(_openbabel.vVector3_rend,None,vVector3)
vVector3.pop_back = new_instancemethod(_openbabel.vVector3_pop_back,None,vVector3)
vVector3.erase = new_instancemethod(_openbabel.vVector3_erase,None,vVector3)
vVector3.push_back = new_instancemethod(_openbabel.vVector3_push_back,None,vVector3)
vVector3.front = new_instancemethod(_openbabel.vVector3_front,None,vVector3)
vVector3.back = new_instancemethod(_openbabel.vVector3_back,None,vVector3)
vVector3.assign = new_instancemethod(_openbabel.vVector3_assign,None,vVector3)
vVector3.resize = new_instancemethod(_openbabel.vVector3_resize,None,vVector3)
vVector3.insert = new_instancemethod(_openbabel.vVector3_insert,None,vVector3)
vVector3.reserve = new_instancemethod(_openbabel.vVector3_reserve,None,vVector3)
vVector3.capacity = new_instancemethod(_openbabel.vVector3_capacity,None,vVector3)
vVector3_swigregister = _openbabel.vVector3_swigregister
vVector3_swigregister(vVector3)

class vectorMol(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __iter__(self): return self.iterator()
    def __init__(self, *args): 
        _openbabel.vectorMol_swiginit(self,_openbabel.new_vectorMol(*args))
    __swig_destroy__ = _openbabel.delete_vectorMol
vectorMol.iterator = new_instancemethod(_openbabel.vectorMol_iterator,None,vectorMol)
vectorMol.__nonzero__ = new_instancemethod(_openbabel.vectorMol___nonzero__,None,vectorMol)
vectorMol.__len__ = new_instancemethod(_openbabel.vectorMol___len__,None,vectorMol)
vectorMol.pop = new_instancemethod(_openbabel.vectorMol_pop,None,vectorMol)
vectorMol.__getslice__ = new_instancemethod(_openbabel.vectorMol___getslice__,None,vectorMol)
vectorMol.__setslice__ = new_instancemethod(_openbabel.vectorMol___setslice__,None,vectorMol)
vectorMol.__delslice__ = new_instancemethod(_openbabel.vectorMol___delslice__,None,vectorMol)
vectorMol.__delitem__ = new_instancemethod(_openbabel.vectorMol___delitem__,None,vectorMol)
vectorMol.__getitem__ = new_instancemethod(_openbabel.vectorMol___getitem__,None,vectorMol)
vectorMol.__setitem__ = new_instancemethod(_openbabel.vectorMol___setitem__,None,vectorMol)
vectorMol.append = new_instancemethod(_openbabel.vectorMol_append,None,vectorMol)
vectorMol.empty = new_instancemethod(_openbabel.vectorMol_empty,None,vectorMol)
vectorMol.size = new_instancemethod(_openbabel.vectorMol_size,None,vectorMol)
vectorMol.clear = new_instancemethod(_openbabel.vectorMol_clear,None,vectorMol)
vectorMol.swap = new_instancemethod(_openbabel.vectorMol_swap,None,vectorMol)
vectorMol.get_allocator = new_instancemethod(_openbabel.vectorMol_get_allocator,None,vectorMol)
vectorMol.begin = new_instancemethod(_openbabel.vectorMol_begin,None,vectorMol)
vectorMol.end = new_instancemethod(_openbabel.vectorMol_end,None,vectorMol)
vectorMol.rbegin = new_instancemethod(_openbabel.vectorMol_rbegin,None,vectorMol)
vectorMol.rend = new_instancemethod(_openbabel.vectorMol_rend,None,vectorMol)
vectorMol.pop_back = new_instancemethod(_openbabel.vectorMol_pop_back,None,vectorMol)
vectorMol.erase = new_instancemethod(_openbabel.vectorMol_erase,None,vectorMol)
vectorMol.push_back = new_instancemethod(_openbabel.vectorMol_push_back,None,vectorMol)
vectorMol.front = new_instancemethod(_openbabel.vectorMol_front,None,vectorMol)
vectorMol.back = new_instancemethod(_openbabel.vectorMol_back,None,vectorMol)
vectorMol.assign = new_instancemethod(_openbabel.vectorMol_assign,None,vectorMol)
vectorMol.resize = new_instancemethod(_openbabel.vectorMol_resize,None,vectorMol)
vectorMol.insert = new_instancemethod(_openbabel.vectorMol_insert,None,vectorMol)
vectorMol.reserve = new_instancemethod(_openbabel.vectorMol_reserve,None,vectorMol)
vectorMol.capacity = new_instancemethod(_openbabel.vectorMol_capacity,None,vectorMol)
vectorMol_swigregister = _openbabel.vectorMol_swigregister
vectorMol_swigregister(vectorMol)

class vectorBond(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __iter__(self): return self.iterator()
    def __init__(self, *args): 
        _openbabel.vectorBond_swiginit(self,_openbabel.new_vectorBond(*args))
    __swig_destroy__ = _openbabel.delete_vectorBond
vectorBond.iterator = new_instancemethod(_openbabel.vectorBond_iterator,None,vectorBond)
vectorBond.__nonzero__ = new_instancemethod(_openbabel.vectorBond___nonzero__,None,vectorBond)
vectorBond.__len__ = new_instancemethod(_openbabel.vectorBond___len__,None,vectorBond)
vectorBond.pop = new_instancemethod(_openbabel.vectorBond_pop,None,vectorBond)
vectorBond.__getslice__ = new_instancemethod(_openbabel.vectorBond___getslice__,None,vectorBond)
vectorBond.__setslice__ = new_instancemethod(_openbabel.vectorBond___setslice__,None,vectorBond)
vectorBond.__delslice__ = new_instancemethod(_openbabel.vectorBond___delslice__,None,vectorBond)
vectorBond.__delitem__ = new_instancemethod(_openbabel.vectorBond___delitem__,None,vectorBond)
vectorBond.__getitem__ = new_instancemethod(_openbabel.vectorBond___getitem__,None,vectorBond)
vectorBond.__setitem__ = new_instancemethod(_openbabel.vectorBond___setitem__,None,vectorBond)
vectorBond.append = new_instancemethod(_openbabel.vectorBond_append,None,vectorBond)
vectorBond.empty = new_instancemethod(_openbabel.vectorBond_empty,None,vectorBond)
vectorBond.size = new_instancemethod(_openbabel.vectorBond_size,None,vectorBond)
vectorBond.clear = new_instancemethod(_openbabel.vectorBond_clear,None,vectorBond)
vectorBond.swap = new_instancemethod(_openbabel.vectorBond_swap,None,vectorBond)
vectorBond.get_allocator = new_instancemethod(_openbabel.vectorBond_get_allocator,None,vectorBond)
vectorBond.begin = new_instancemethod(_openbabel.vectorBond_begin,None,vectorBond)
vectorBond.end = new_instancemethod(_openbabel.vectorBond_end,None,vectorBond)
vectorBond.rbegin = new_instancemethod(_openbabel.vectorBond_rbegin,None,vectorBond)
vectorBond.rend = new_instancemethod(_openbabel.vectorBond_rend,None,vectorBond)
vectorBond.pop_back = new_instancemethod(_openbabel.vectorBond_pop_back,None,vectorBond)
vectorBond.erase = new_instancemethod(_openbabel.vectorBond_erase,None,vectorBond)
vectorBond.push_back = new_instancemethod(_openbabel.vectorBond_push_back,None,vectorBond)
vectorBond.front = new_instancemethod(_openbabel.vectorBond_front,None,vectorBond)
vectorBond.back = new_instancemethod(_openbabel.vectorBond_back,None,vectorBond)
vectorBond.assign = new_instancemethod(_openbabel.vectorBond_assign,None,vectorBond)
vectorBond.resize = new_instancemethod(_openbabel.vectorBond_resize,None,vectorBond)
vectorBond.insert = new_instancemethod(_openbabel.vectorBond_insert,None,vectorBond)
vectorBond.reserve = new_instancemethod(_openbabel.vectorBond_reserve,None,vectorBond)
vectorBond.capacity = new_instancemethod(_openbabel.vectorBond_capacity,None,vectorBond)
vectorBond_swigregister = _openbabel.vectorBond_swigregister
vectorBond_swigregister(vectorBond)

class vectorResidue(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __iter__(self): return self.iterator()
    def __init__(self, *args): 
        _openbabel.vectorResidue_swiginit(self,_openbabel.new_vectorResidue(*args))
    __swig_destroy__ = _openbabel.delete_vectorResidue
vectorResidue.iterator = new_instancemethod(_openbabel.vectorResidue_iterator,None,vectorResidue)
vectorResidue.__nonzero__ = new_instancemethod(_openbabel.vectorResidue___nonzero__,None,vectorResidue)
vectorResidue.__len__ = new_instancemethod(_openbabel.vectorResidue___len__,None,vectorResidue)
vectorResidue.pop = new_instancemethod(_openbabel.vectorResidue_pop,None,vectorResidue)
vectorResidue.__getslice__ = new_instancemethod(_openbabel.vectorResidue___getslice__,None,vectorResidue)
vectorResidue.__setslice__ = new_instancemethod(_openbabel.vectorResidue___setslice__,None,vectorResidue)
vectorResidue.__delslice__ = new_instancemethod(_openbabel.vectorResidue___delslice__,None,vectorResidue)
vectorResidue.__delitem__ = new_instancemethod(_openbabel.vectorResidue___delitem__,None,vectorResidue)
vectorResidue.__getitem__ = new_instancemethod(_openbabel.vectorResidue___getitem__,None,vectorResidue)
vectorResidue.__setitem__ = new_instancemethod(_openbabel.vectorResidue___setitem__,None,vectorResidue)
vectorResidue.append = new_instancemethod(_openbabel.vectorResidue_append,None,vectorResidue)
vectorResidue.empty = new_instancemethod(_openbabel.vectorResidue_empty,None,vectorResidue)
vectorResidue.size = new_instancemethod(_openbabel.vectorResidue_size,None,vectorResidue)
vectorResidue.clear = new_instancemethod(_openbabel.vectorResidue_clear,None,vectorResidue)
vectorResidue.swap = new_instancemethod(_openbabel.vectorResidue_swap,None,vectorResidue)
vectorResidue.get_allocator = new_instancemethod(_openbabel.vectorResidue_get_allocator,None,vectorResidue)
vectorResidue.begin = new_instancemethod(_openbabel.vectorResidue_begin,None,vectorResidue)
vectorResidue.end = new_instancemethod(_openbabel.vectorResidue_end,None,vectorResidue)
vectorResidue.rbegin = new_instancemethod(_openbabel.vectorResidue_rbegin,None,vectorResidue)
vectorResidue.rend = new_instancemethod(_openbabel.vectorResidue_rend,None,vectorResidue)
vectorResidue.pop_back = new_instancemethod(_openbabel.vectorResidue_pop_back,None,vectorResidue)
vectorResidue.erase = new_instancemethod(_openbabel.vectorResidue_erase,None,vectorResidue)
vectorResidue.push_back = new_instancemethod(_openbabel.vectorResidue_push_back,None,vectorResidue)
vectorResidue.front = new_instancemethod(_openbabel.vectorResidue_front,None,vectorResidue)
vectorResidue.back = new_instancemethod(_openbabel.vectorResidue_back,None,vectorResidue)
vectorResidue.assign = new_instancemethod(_openbabel.vectorResidue_assign,None,vectorResidue)
vectorResidue.resize = new_instancemethod(_openbabel.vectorResidue_resize,None,vectorResidue)
vectorResidue.insert = new_instancemethod(_openbabel.vectorResidue_insert,None,vectorResidue)
vectorResidue.reserve = new_instancemethod(_openbabel.vectorResidue_reserve,None,vectorResidue)
vectorResidue.capacity = new_instancemethod(_openbabel.vectorResidue_capacity,None,vectorResidue)
vectorResidue_swigregister = _openbabel.vectorResidue_swigregister
vectorResidue_swigregister(vectorResidue)

class vectorRing(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __iter__(self): return self.iterator()
    def __init__(self, *args): 
        _openbabel.vectorRing_swiginit(self,_openbabel.new_vectorRing(*args))
    __swig_destroy__ = _openbabel.delete_vectorRing
vectorRing.iterator = new_instancemethod(_openbabel.vectorRing_iterator,None,vectorRing)
vectorRing.__nonzero__ = new_instancemethod(_openbabel.vectorRing___nonzero__,None,vectorRing)
vectorRing.__len__ = new_instancemethod(_openbabel.vectorRing___len__,None,vectorRing)
vectorRing.pop = new_instancemethod(_openbabel.vectorRing_pop,None,vectorRing)
vectorRing.__getslice__ = new_instancemethod(_openbabel.vectorRing___getslice__,None,vectorRing)
vectorRing.__setslice__ = new_instancemethod(_openbabel.vectorRing___setslice__,None,vectorRing)
vectorRing.__delslice__ = new_instancemethod(_openbabel.vectorRing___delslice__,None,vectorRing)
vectorRing.__delitem__ = new_instancemethod(_openbabel.vectorRing___delitem__,None,vectorRing)
vectorRing.__getitem__ = new_instancemethod(_openbabel.vectorRing___getitem__,None,vectorRing)
vectorRing.__setitem__ = new_instancemethod(_openbabel.vectorRing___setitem__,None,vectorRing)
vectorRing.append = new_instancemethod(_openbabel.vectorRing_append,None,vectorRing)
vectorRing.empty = new_instancemethod(_openbabel.vectorRing_empty,None,vectorRing)
vectorRing.size = new_instancemethod(_openbabel.vectorRing_size,None,vectorRing)
vectorRing.clear = new_instancemethod(_openbabel.vectorRing_clear,None,vectorRing)
vectorRing.swap = new_instancemethod(_openbabel.vectorRing_swap,None,vectorRing)
vectorRing.get_allocator = new_instancemethod(_openbabel.vectorRing_get_allocator,None,vectorRing)
vectorRing.begin = new_instancemethod(_openbabel.vectorRing_begin,None,vectorRing)
vectorRing.end = new_instancemethod(_openbabel.vectorRing_end,None,vectorRing)
vectorRing.rbegin = new_instancemethod(_openbabel.vectorRing_rbegin,None,vectorRing)
vectorRing.rend = new_instancemethod(_openbabel.vectorRing_rend,None,vectorRing)
vectorRing.pop_back = new_instancemethod(_openbabel.vectorRing_pop_back,None,vectorRing)
vectorRing.erase = new_instancemethod(_openbabel.vectorRing_erase,None,vectorRing)
vectorRing.push_back = new_instancemethod(_openbabel.vectorRing_push_back,None,vectorRing)
vectorRing.front = new_instancemethod(_openbabel.vectorRing_front,None,vectorRing)
vectorRing.back = new_instancemethod(_openbabel.vectorRing_back,None,vectorRing)
vectorRing.assign = new_instancemethod(_openbabel.vectorRing_assign,None,vectorRing)
vectorRing.resize = new_instancemethod(_openbabel.vectorRing_resize,None,vectorRing)
vectorRing.insert = new_instancemethod(_openbabel.vectorRing_insert,None,vectorRing)
vectorRing.reserve = new_instancemethod(_openbabel.vectorRing_reserve,None,vectorRing)
vectorRing.capacity = new_instancemethod(_openbabel.vectorRing_capacity,None,vectorRing)
vectorRing_swigregister = _openbabel.vectorRing_swigregister
vectorRing_swigregister(vectorRing)

class OBGlobalDataBase(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBGlobalDataBase_swiginit(self,_openbabel.new_OBGlobalDataBase(*args))
    __swig_destroy__ = _openbabel.delete_OBGlobalDataBase
OBGlobalDataBase.Init = new_instancemethod(_openbabel.OBGlobalDataBase_Init,None,OBGlobalDataBase)
OBGlobalDataBase.GetSize = new_instancemethod(_openbabel.OBGlobalDataBase_GetSize,None,OBGlobalDataBase)
OBGlobalDataBase.SetReadDirectory = new_instancemethod(_openbabel.OBGlobalDataBase_SetReadDirectory,None,OBGlobalDataBase)
OBGlobalDataBase.SetEnvironmentVariable = new_instancemethod(_openbabel.OBGlobalDataBase_SetEnvironmentVariable,None,OBGlobalDataBase)
OBGlobalDataBase.ParseLine = new_instancemethod(_openbabel.OBGlobalDataBase_ParseLine,None,OBGlobalDataBase)
OBGlobalDataBase_swigregister = _openbabel.OBGlobalDataBase_swigregister
OBGlobalDataBase_swigregister(OBGlobalDataBase)

class OBElement(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBElement_swiginit(self,_openbabel.new_OBElement(*args))
    __swig_destroy__ = _openbabel.delete_OBElement
OBElement.GetAtomicNum = new_instancemethod(_openbabel.OBElement_GetAtomicNum,None,OBElement)
OBElement.GetSymbol = new_instancemethod(_openbabel.OBElement_GetSymbol,None,OBElement)
OBElement.GetCovalentRad = new_instancemethod(_openbabel.OBElement_GetCovalentRad,None,OBElement)
OBElement.GetVdwRad = new_instancemethod(_openbabel.OBElement_GetVdwRad,None,OBElement)
OBElement.GetMass = new_instancemethod(_openbabel.OBElement_GetMass,None,OBElement)
OBElement.GetMaxBonds = new_instancemethod(_openbabel.OBElement_GetMaxBonds,None,OBElement)
OBElement.GetElectroNeg = new_instancemethod(_openbabel.OBElement_GetElectroNeg,None,OBElement)
OBElement.GetIonization = new_instancemethod(_openbabel.OBElement_GetIonization,None,OBElement)
OBElement.GetElectronAffinity = new_instancemethod(_openbabel.OBElement_GetElectronAffinity,None,OBElement)
OBElement.GetName = new_instancemethod(_openbabel.OBElement_GetName,None,OBElement)
OBElement.GetRed = new_instancemethod(_openbabel.OBElement_GetRed,None,OBElement)
OBElement.GetGreen = new_instancemethod(_openbabel.OBElement_GetGreen,None,OBElement)
OBElement.GetBlue = new_instancemethod(_openbabel.OBElement_GetBlue,None,OBElement)
OBElement_swigregister = _openbabel.OBElement_swigregister
OBElement_swigregister(OBElement)

class OBElementTable(OBGlobalDataBase):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBElementTable_swiginit(self,_openbabel.new_OBElementTable(*args))
    __swig_destroy__ = _openbabel.delete_OBElementTable
OBElementTable.GetNumberOfElements = new_instancemethod(_openbabel.OBElementTable_GetNumberOfElements,None,OBElementTable)
OBElementTable.GetAtomicNum = new_instancemethod(_openbabel.OBElementTable_GetAtomicNum,None,OBElementTable)
OBElementTable.GetSymbol = new_instancemethod(_openbabel.OBElementTable_GetSymbol,None,OBElementTable)
OBElementTable.GetVdwRad = new_instancemethod(_openbabel.OBElementTable_GetVdwRad,None,OBElementTable)
OBElementTable.GetCovalentRad = new_instancemethod(_openbabel.OBElementTable_GetCovalentRad,None,OBElementTable)
OBElementTable.GetMass = new_instancemethod(_openbabel.OBElementTable_GetMass,None,OBElementTable)
OBElementTable.CorrectedBondRad = new_instancemethod(_openbabel.OBElementTable_CorrectedBondRad,None,OBElementTable)
OBElementTable.CorrectedVdwRad = new_instancemethod(_openbabel.OBElementTable_CorrectedVdwRad,None,OBElementTable)
OBElementTable.GetMaxBonds = new_instancemethod(_openbabel.OBElementTable_GetMaxBonds,None,OBElementTable)
OBElementTable.GetElectroNeg = new_instancemethod(_openbabel.OBElementTable_GetElectroNeg,None,OBElementTable)
OBElementTable.GetIonization = new_instancemethod(_openbabel.OBElementTable_GetIonization,None,OBElementTable)
OBElementTable.GetElectronAffinity = new_instancemethod(_openbabel.OBElementTable_GetElectronAffinity,None,OBElementTable)
OBElementTable.GetRGB = new_instancemethod(_openbabel.OBElementTable_GetRGB,None,OBElementTable)
OBElementTable.GetName = new_instancemethod(_openbabel.OBElementTable_GetName,None,OBElementTable)
OBElementTable_swigregister = _openbabel.OBElementTable_swigregister
OBElementTable_swigregister(OBElementTable)

class OBIsotopeTable(OBGlobalDataBase):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBIsotopeTable_swiginit(self,_openbabel.new_OBIsotopeTable(*args))
    __swig_destroy__ = _openbabel.delete_OBIsotopeTable
OBIsotopeTable.GetExactMass = new_instancemethod(_openbabel.OBIsotopeTable_GetExactMass,None,OBIsotopeTable)
OBIsotopeTable_swigregister = _openbabel.OBIsotopeTable_swigregister
OBIsotopeTable_swigregister(OBIsotopeTable)

class OBTypeTable(OBGlobalDataBase):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBTypeTable_swiginit(self,_openbabel.new_OBTypeTable(*args))
    __swig_destroy__ = _openbabel.delete_OBTypeTable
OBTypeTable.SetFromType = new_instancemethod(_openbabel.OBTypeTable_SetFromType,None,OBTypeTable)
OBTypeTable.SetToType = new_instancemethod(_openbabel.OBTypeTable_SetToType,None,OBTypeTable)
OBTypeTable.Translate = new_instancemethod(_openbabel.OBTypeTable_Translate,None,OBTypeTable)
OBTypeTable.GetFromType = new_instancemethod(_openbabel.OBTypeTable_GetFromType,None,OBTypeTable)
OBTypeTable.GetToType = new_instancemethod(_openbabel.OBTypeTable_GetToType,None,OBTypeTable)
OBTypeTable_swigregister = _openbabel.OBTypeTable_swigregister
OBTypeTable_swigregister(OBTypeTable)

class OBResidueData(OBGlobalDataBase):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBResidueData_swiginit(self,_openbabel.new_OBResidueData(*args))
    __swig_destroy__ = _openbabel.delete_OBResidueData
OBResidueData.SetResName = new_instancemethod(_openbabel.OBResidueData_SetResName,None,OBResidueData)
OBResidueData.LookupBO = new_instancemethod(_openbabel.OBResidueData_LookupBO,None,OBResidueData)
OBResidueData.LookupType = new_instancemethod(_openbabel.OBResidueData_LookupType,None,OBResidueData)
OBResidueData.AssignBonds = new_instancemethod(_openbabel.OBResidueData_AssignBonds,None,OBResidueData)
OBResidueData_swigregister = _openbabel.OBResidueData_swigregister
OBResidueData_swigregister(OBResidueData)

FILE_SEP_CHAR = _openbabel.FILE_SEP_CHAR
class OBStopwatch(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBStopwatch_swiginit(self,_openbabel.new_OBStopwatch(*args))
    __swig_destroy__ = _openbabel.delete_OBStopwatch
OBStopwatch.Start = new_instancemethod(_openbabel.OBStopwatch_Start,None,OBStopwatch)
OBStopwatch.Lap = new_instancemethod(_openbabel.OBStopwatch_Lap,None,OBStopwatch)
OBStopwatch.Elapsed = new_instancemethod(_openbabel.OBStopwatch_Elapsed,None,OBStopwatch)
OBStopwatch_swigregister = _openbabel.OBStopwatch_swigregister
OBStopwatch_swigregister(OBStopwatch)

class OBSqrtTbl(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBSqrtTbl_swiginit(self,_openbabel.new_OBSqrtTbl(*args))
    __swig_destroy__ = _openbabel.delete_OBSqrtTbl
OBSqrtTbl.Sqrt = new_instancemethod(_openbabel.OBSqrtTbl_Sqrt,None,OBSqrtTbl)
OBSqrtTbl.Init = new_instancemethod(_openbabel.OBSqrtTbl_Init,None,OBSqrtTbl)
OBSqrtTbl_swigregister = _openbabel.OBSqrtTbl_swigregister
OBSqrtTbl_swigregister(OBSqrtTbl)

class DoubleType(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    hi = property(_openbabel.DoubleType_hi_get, _openbabel.DoubleType_hi_set)
    lo = property(_openbabel.DoubleType_lo_get, _openbabel.DoubleType_lo_set)
    def __init__(self, *args): 
        _openbabel.DoubleType_swiginit(self,_openbabel.new_DoubleType(*args))
    __swig_destroy__ = _openbabel.delete_DoubleType
DoubleType_swigregister = _openbabel.DoubleType_swigregister
DoubleType_swigregister(DoubleType)

DoubleMultiply = _openbabel.DoubleMultiply
DoubleAdd = _openbabel.DoubleAdd
DoubleModulus = _openbabel.DoubleModulus
class OBRandom(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBRandom_swiginit(self,_openbabel.new_OBRandom(*args))
    __swig_destroy__ = _openbabel.delete_OBRandom
OBRandom.Seed = new_instancemethod(_openbabel.OBRandom_Seed,None,OBRandom)
OBRandom.TimeSeed = new_instancemethod(_openbabel.OBRandom_TimeSeed,None,OBRandom)
OBRandom.NextInt = new_instancemethod(_openbabel.OBRandom_NextInt,None,OBRandom)
OBRandom.NextFloat = new_instancemethod(_openbabel.OBRandom_NextFloat,None,OBRandom)
OBRandom_swigregister = _openbabel.OBRandom_swigregister
OBRandom_swigregister(OBRandom)

calc_rms = _openbabel.calc_rms
CleanAtomType = _openbabel.CleanAtomType
OBCompareInt = _openbabel.OBCompareInt
OBCompareUnsigned = _openbabel.OBCompareUnsigned
PI = _openbabel.PI
RAD_TO_DEG = _openbabel.RAD_TO_DEG
DEG_TO_RAD = _openbabel.DEG_TO_RAD
class vector3(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.vector3_swiginit(self,_openbabel.new_vector3(*args))
    __swig_destroy__ = _openbabel.delete_vector3
vector3.Set = new_instancemethod(_openbabel.vector3_Set,None,vector3)
vector3.SetX = new_instancemethod(_openbabel.vector3_SetX,None,vector3)
vector3.SetY = new_instancemethod(_openbabel.vector3_SetY,None,vector3)
vector3.SetZ = new_instancemethod(_openbabel.vector3_SetZ,None,vector3)
vector3.Get = new_instancemethod(_openbabel.vector3_Get,None,vector3)
vector3.__iadd__ = new_instancemethod(_openbabel.vector3___iadd__,None,vector3)
vector3.__isub__ = new_instancemethod(_openbabel.vector3___isub__,None,vector3)
vector3.__idiv__ = new_instancemethod(_openbabel.vector3___idiv__,None,vector3)
vector3.__imul__ = new_instancemethod(_openbabel.vector3___imul__,None,vector3)
vector3.randomUnitVector = new_instancemethod(_openbabel.vector3_randomUnitVector,None,vector3)
vector3.normalize = new_instancemethod(_openbabel.vector3_normalize,None,vector3)
vector3.length = new_instancemethod(_openbabel.vector3_length,None,vector3)
vector3.length_2 = new_instancemethod(_openbabel.vector3_length_2,None,vector3)
vector3.x = new_instancemethod(_openbabel.vector3_x,None,vector3)
vector3.y = new_instancemethod(_openbabel.vector3_y,None,vector3)
vector3.z = new_instancemethod(_openbabel.vector3_z,None,vector3)
vector3.distSq = new_instancemethod(_openbabel.vector3_distSq,None,vector3)
vector3.createOrthoVector = new_instancemethod(_openbabel.vector3_createOrthoVector,None,vector3)
vector3_swigregister = _openbabel.vector3_swigregister
vector3_swigregister(vector3)
ToUpper = _openbabel.ToUpper
ToLower = _openbabel.ToLower
IsNear = _openbabel.IsNear
IsNearZero = _openbabel.IsNearZero

Point2Plane = _openbabel.Point2Plane
class OBFormat(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    def __init__(self): raise AttributeError, "No constructor defined"
    __repr__ = _swig_repr
    __swig_destroy__ = _openbabel.delete_OBFormat
OBFormat.ReadMolecule = new_instancemethod(_openbabel.OBFormat_ReadMolecule,None,OBFormat)
OBFormat.ReadChemObject = new_instancemethod(_openbabel.OBFormat_ReadChemObject,None,OBFormat)
OBFormat.WriteMolecule = new_instancemethod(_openbabel.OBFormat_WriteMolecule,None,OBFormat)
OBFormat.WriteChemObject = new_instancemethod(_openbabel.OBFormat_WriteChemObject,None,OBFormat)
OBFormat.Description = new_instancemethod(_openbabel.OBFormat_Description,None,OBFormat)
OBFormat.TargetClassDescription = new_instancemethod(_openbabel.OBFormat_TargetClassDescription,None,OBFormat)
OBFormat.GetType = new_instancemethod(_openbabel.OBFormat_GetType,None,OBFormat)
OBFormat.SpecificationURL = new_instancemethod(_openbabel.OBFormat_SpecificationURL,None,OBFormat)
OBFormat.GetMIMEType = new_instancemethod(_openbabel.OBFormat_GetMIMEType,None,OBFormat)
OBFormat.Flags = new_instancemethod(_openbabel.OBFormat_Flags,None,OBFormat)
OBFormat.SkipObjects = new_instancemethod(_openbabel.OBFormat_SkipObjects,None,OBFormat)
OBFormat.MakeNewInstance = new_instancemethod(_openbabel.OBFormat_MakeNewInstance,None,OBFormat)
OBFormat_swigregister = _openbabel.OBFormat_swigregister
OBFormat_swigregister(OBFormat)
cvar = _openbabel.cvar
VZero = cvar.VZero
VX = cvar.VX
VY = cvar.VY
VZ = cvar.VZ

class CharPtrLess(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.CharPtrLess_swiginit(self,_openbabel.new_CharPtrLess(*args))
    __swig_destroy__ = _openbabel.delete_CharPtrLess
CharPtrLess.__call__ = new_instancemethod(_openbabel.CharPtrLess___call__,None,CharPtrLess)
CharPtrLess_swigregister = _openbabel.CharPtrLess_swigregister
CharPtrLess_swigregister(CharPtrLess)

class OBConversion(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBConversion_swiginit(self,_openbabel.new_OBConversion(*args))
    __swig_destroy__ = _openbabel.delete_OBConversion
    RegisterFormat = staticmethod(_openbabel.OBConversion_RegisterFormat)
    FindFormat = staticmethod(_openbabel.OBConversion_FindFormat)
    FormatFromExt = staticmethod(_openbabel.OBConversion_FormatFromExt)
    FormatFromMIME = staticmethod(_openbabel.OBConversion_FormatFromMIME)
    GetNextFormat = staticmethod(_openbabel.OBConversion_GetNextFormat)
    Description = staticmethod(_openbabel.OBConversion_Description)
    INOPTIONS = _openbabel.OBConversion_INOPTIONS
    OUTOPTIONS = _openbabel.OBConversion_OUTOPTIONS
    GENOPTIONS = _openbabel.OBConversion_GENOPTIONS
    RegisterOptionParam = staticmethod(_openbabel.OBConversion_RegisterOptionParam)
    GetOptionParams = staticmethod(_openbabel.OBConversion_GetOptionParams)
    GetDefaultFormat = staticmethod(_openbabel.OBConversion_GetDefaultFormat)
    BatchFileName = staticmethod(_openbabel.OBConversion_BatchFileName)
    IncrementedFileName = staticmethod(_openbabel.OBConversion_IncrementedFileName)
OBConversion.GetInStream = new_instancemethod(_openbabel.OBConversion_GetInStream,None,OBConversion)
OBConversion.GetOutStream = new_instancemethod(_openbabel.OBConversion_GetOutStream,None,OBConversion)
OBConversion.SetInStream = new_instancemethod(_openbabel.OBConversion_SetInStream,None,OBConversion)
OBConversion.SetOutStream = new_instancemethod(_openbabel.OBConversion_SetOutStream,None,OBConversion)
OBConversion.SetInAndOutFormats = new_instancemethod(_openbabel.OBConversion_SetInAndOutFormats,None,OBConversion)
OBConversion.SetInFormat = new_instancemethod(_openbabel.OBConversion_SetInFormat,None,OBConversion)
OBConversion.SetOutFormat = new_instancemethod(_openbabel.OBConversion_SetOutFormat,None,OBConversion)
OBConversion.GetInFormat = new_instancemethod(_openbabel.OBConversion_GetInFormat,None,OBConversion)
OBConversion.GetOutFormat = new_instancemethod(_openbabel.OBConversion_GetOutFormat,None,OBConversion)
OBConversion.GetInFilename = new_instancemethod(_openbabel.OBConversion_GetInFilename,None,OBConversion)
OBConversion.GetInPos = new_instancemethod(_openbabel.OBConversion_GetInPos,None,OBConversion)
OBConversion.GetInLen = new_instancemethod(_openbabel.OBConversion_GetInLen,None,OBConversion)
OBConversion.GetTitle = new_instancemethod(_openbabel.OBConversion_GetTitle,None,OBConversion)
OBConversion.GetAuxConv = new_instancemethod(_openbabel.OBConversion_GetAuxConv,None,OBConversion)
OBConversion.SetAuxConv = new_instancemethod(_openbabel.OBConversion_SetAuxConv,None,OBConversion)
OBConversion.IsOption = new_instancemethod(_openbabel.OBConversion_IsOption,None,OBConversion)
OBConversion.GetOptions = new_instancemethod(_openbabel.OBConversion_GetOptions,None,OBConversion)
OBConversion.AddOption = new_instancemethod(_openbabel.OBConversion_AddOption,None,OBConversion)
OBConversion.RemoveOption = new_instancemethod(_openbabel.OBConversion_RemoveOption,None,OBConversion)
OBConversion.SetOptions = new_instancemethod(_openbabel.OBConversion_SetOptions,None,OBConversion)
OBConversion.Convert = new_instancemethod(_openbabel.OBConversion_Convert,None,OBConversion)
OBConversion.FullConvert = new_instancemethod(_openbabel.OBConversion_FullConvert,None,OBConversion)
OBConversion.AddChemObject = new_instancemethod(_openbabel.OBConversion_AddChemObject,None,OBConversion)
OBConversion.GetChemObject = new_instancemethod(_openbabel.OBConversion_GetChemObject,None,OBConversion)
OBConversion.IsLast = new_instancemethod(_openbabel.OBConversion_IsLast,None,OBConversion)
OBConversion.IsFirstInput = new_instancemethod(_openbabel.OBConversion_IsFirstInput,None,OBConversion)
OBConversion.GetOutputIndex = new_instancemethod(_openbabel.OBConversion_GetOutputIndex,None,OBConversion)
OBConversion.SetOutputIndex = new_instancemethod(_openbabel.OBConversion_SetOutputIndex,None,OBConversion)
OBConversion.SetMoreFilesToCome = new_instancemethod(_openbabel.OBConversion_SetMoreFilesToCome,None,OBConversion)
OBConversion.SetOneObjectOnly = new_instancemethod(_openbabel.OBConversion_SetOneObjectOnly,None,OBConversion)
OBConversion.Write = new_instancemethod(_openbabel.OBConversion_Write,None,OBConversion)
OBConversion.WriteString = new_instancemethod(_openbabel.OBConversion_WriteString,None,OBConversion)
OBConversion.WriteFile = new_instancemethod(_openbabel.OBConversion_WriteFile,None,OBConversion)
OBConversion.Read = new_instancemethod(_openbabel.OBConversion_Read,None,OBConversion)
OBConversion.ReadString = new_instancemethod(_openbabel.OBConversion_ReadString,None,OBConversion)
OBConversion.ReadFile = new_instancemethod(_openbabel.OBConversion_ReadFile,None,OBConversion)
OBConversion_swigregister = _openbabel.OBConversion_swigregister
OBConversion_swigregister(OBConversion)
OBConversion_RegisterFormat = _openbabel.OBConversion_RegisterFormat
OBConversion_FindFormat = _openbabel.OBConversion_FindFormat
OBConversion_FormatFromExt = _openbabel.OBConversion_FormatFromExt
OBConversion_FormatFromMIME = _openbabel.OBConversion_FormatFromMIME
OBConversion_GetNextFormat = _openbabel.OBConversion_GetNextFormat
OBConversion_Description = _openbabel.OBConversion_Description
OBConversion_RegisterOptionParam = _openbabel.OBConversion_RegisterOptionParam
OBConversion_GetOptionParams = _openbabel.OBConversion_GetOptionParams
OBConversion_GetDefaultFormat = _openbabel.OBConversion_GetDefaultFormat
OBConversion_BatchFileName = _openbabel.OBConversion_BatchFileName
OBConversion_IncrementedFileName = _openbabel.OBConversion_IncrementedFileName

NOTREADABLE = _openbabel.NOTREADABLE
READONEONLY = _openbabel.READONEONLY
READBINARY = _openbabel.READBINARY
NOTWRITABLE = _openbabel.NOTWRITABLE
WRITEONEONLY = _openbabel.WRITEONEONLY
WRITEBINARY = _openbabel.WRITEBINARY
DEFAULTFORMAT = _openbabel.DEFAULTFORMAT
class OBResidue(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBResidue_swiginit(self,_openbabel.new_OBResidue(*args))
    __swig_destroy__ = _openbabel.delete_OBResidue
OBResidue.AddAtom = new_instancemethod(_openbabel.OBResidue_AddAtom,None,OBResidue)
OBResidue.InsertAtom = new_instancemethod(_openbabel.OBResidue_InsertAtom,None,OBResidue)
OBResidue.RemoveAtom = new_instancemethod(_openbabel.OBResidue_RemoveAtom,None,OBResidue)
OBResidue.Clear = new_instancemethod(_openbabel.OBResidue_Clear,None,OBResidue)
OBResidue.SetName = new_instancemethod(_openbabel.OBResidue_SetName,None,OBResidue)
OBResidue.SetNum = new_instancemethod(_openbabel.OBResidue_SetNum,None,OBResidue)
OBResidue.SetChain = new_instancemethod(_openbabel.OBResidue_SetChain,None,OBResidue)
OBResidue.SetChainNum = new_instancemethod(_openbabel.OBResidue_SetChainNum,None,OBResidue)
OBResidue.SetIdx = new_instancemethod(_openbabel.OBResidue_SetIdx,None,OBResidue)
OBResidue.SetAtomID = new_instancemethod(_openbabel.OBResidue_SetAtomID,None,OBResidue)
OBResidue.SetHetAtom = new_instancemethod(_openbabel.OBResidue_SetHetAtom,None,OBResidue)
OBResidue.SetSerialNum = new_instancemethod(_openbabel.OBResidue_SetSerialNum,None,OBResidue)
OBResidue.GetName = new_instancemethod(_openbabel.OBResidue_GetName,None,OBResidue)
OBResidue.GetNum = new_instancemethod(_openbabel.OBResidue_GetNum,None,OBResidue)
OBResidue.GetNumAtoms = new_instancemethod(_openbabel.OBResidue_GetNumAtoms,None,OBResidue)
OBResidue.GetChain = new_instancemethod(_openbabel.OBResidue_GetChain,None,OBResidue)
OBResidue.GetChainNum = new_instancemethod(_openbabel.OBResidue_GetChainNum,None,OBResidue)
OBResidue.GetIdx = new_instancemethod(_openbabel.OBResidue_GetIdx,None,OBResidue)
OBResidue.GetResKey = new_instancemethod(_openbabel.OBResidue_GetResKey,None,OBResidue)
OBResidue.GetAtoms = new_instancemethod(_openbabel.OBResidue_GetAtoms,None,OBResidue)
OBResidue.GetBonds = new_instancemethod(_openbabel.OBResidue_GetBonds,None,OBResidue)
OBResidue.GetAtomID = new_instancemethod(_openbabel.OBResidue_GetAtomID,None,OBResidue)
OBResidue.GetSerialNum = new_instancemethod(_openbabel.OBResidue_GetSerialNum,None,OBResidue)
OBResidue.GetAminoAcidProperty = new_instancemethod(_openbabel.OBResidue_GetAminoAcidProperty,None,OBResidue)
OBResidue.GetAtomProperty = new_instancemethod(_openbabel.OBResidue_GetAtomProperty,None,OBResidue)
OBResidue.GetResidueProperty = new_instancemethod(_openbabel.OBResidue_GetResidueProperty,None,OBResidue)
OBResidue.IsHetAtom = new_instancemethod(_openbabel.OBResidue_IsHetAtom,None,OBResidue)
OBResidue.IsResidueType = new_instancemethod(_openbabel.OBResidue_IsResidueType,None,OBResidue)
OBResidue.BeginAtom = new_instancemethod(_openbabel.OBResidue_BeginAtom,None,OBResidue)
OBResidue.NextAtom = new_instancemethod(_openbabel.OBResidue_NextAtom,None,OBResidue)
OBResidue.HasData = new_instancemethod(_openbabel.OBResidue_HasData,None,OBResidue)
OBResidue.DeleteData = new_instancemethod(_openbabel.OBResidue_DeleteData,None,OBResidue)
OBResidue.SetData = new_instancemethod(_openbabel.OBResidue_SetData,None,OBResidue)
OBResidue.DataSize = new_instancemethod(_openbabel.OBResidue_DataSize,None,OBResidue)
OBResidue.GetData = new_instancemethod(_openbabel.OBResidue_GetData,None,OBResidue)
OBResidue.BeginData = new_instancemethod(_openbabel.OBResidue_BeginData,None,OBResidue)
OBResidue.EndData = new_instancemethod(_openbabel.OBResidue_EndData,None,OBResidue)
OBResidue_swigregister = _openbabel.OBResidue_swigregister
OBResidue_swigregister(OBResidue)

OB_4RING_ATOM = _openbabel.OB_4RING_ATOM
OB_3RING_ATOM = _openbabel.OB_3RING_ATOM
OB_AROMATIC_ATOM = _openbabel.OB_AROMATIC_ATOM
OB_RING_ATOM = _openbabel.OB_RING_ATOM
OB_CSTEREO_ATOM = _openbabel.OB_CSTEREO_ATOM
OB_ACSTEREO_ATOM = _openbabel.OB_ACSTEREO_ATOM
OB_DONOR_ATOM = _openbabel.OB_DONOR_ATOM
OB_ACCEPTOR_ATOM = _openbabel.OB_ACCEPTOR_ATOM
OB_CHIRAL_ATOM = _openbabel.OB_CHIRAL_ATOM
OB_POS_CHIRAL_ATOM = _openbabel.OB_POS_CHIRAL_ATOM
OB_NEG_CHIRAL_ATOM = _openbabel.OB_NEG_CHIRAL_ATOM
class OBAtom(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBAtom_swiginit(self,_openbabel.new_OBAtom(*args))
    __swig_destroy__ = _openbabel.delete_OBAtom
OBAtom.Clear = new_instancemethod(_openbabel.OBAtom_Clear,None,OBAtom)
OBAtom.SetIdx = new_instancemethod(_openbabel.OBAtom_SetIdx,None,OBAtom)
OBAtom.SetHyb = new_instancemethod(_openbabel.OBAtom_SetHyb,None,OBAtom)
OBAtom.SetAtomicNum = new_instancemethod(_openbabel.OBAtom_SetAtomicNum,None,OBAtom)
OBAtom.SetIsotope = new_instancemethod(_openbabel.OBAtom_SetIsotope,None,OBAtom)
OBAtom.SetImplicitValence = new_instancemethod(_openbabel.OBAtom_SetImplicitValence,None,OBAtom)
OBAtom.IncrementImplicitValence = new_instancemethod(_openbabel.OBAtom_IncrementImplicitValence,None,OBAtom)
OBAtom.DecrementImplicitValence = new_instancemethod(_openbabel.OBAtom_DecrementImplicitValence,None,OBAtom)
OBAtom.SetFormalCharge = new_instancemethod(_openbabel.OBAtom_SetFormalCharge,None,OBAtom)
OBAtom.SetSpinMultiplicity = new_instancemethod(_openbabel.OBAtom_SetSpinMultiplicity,None,OBAtom)
OBAtom.SetType = new_instancemethod(_openbabel.OBAtom_SetType,None,OBAtom)
OBAtom.SetPartialCharge = new_instancemethod(_openbabel.OBAtom_SetPartialCharge,None,OBAtom)
OBAtom.SetCoordPtr = new_instancemethod(_openbabel.OBAtom_SetCoordPtr,None,OBAtom)
OBAtom.SetVector = new_instancemethod(_openbabel.OBAtom_SetVector,None,OBAtom)
OBAtom.SetResidue = new_instancemethod(_openbabel.OBAtom_SetResidue,None,OBAtom)
OBAtom.UnsetAromatic = new_instancemethod(_openbabel.OBAtom_UnsetAromatic,None,OBAtom)
OBAtom.SetClockwiseStereo = new_instancemethod(_openbabel.OBAtom_SetClockwiseStereo,None,OBAtom)
OBAtom.SetAntiClockwiseStereo = new_instancemethod(_openbabel.OBAtom_SetAntiClockwiseStereo,None,OBAtom)
OBAtom.SetPositiveStereo = new_instancemethod(_openbabel.OBAtom_SetPositiveStereo,None,OBAtom)
OBAtom.SetNegativeStereo = new_instancemethod(_openbabel.OBAtom_SetNegativeStereo,None,OBAtom)
OBAtom.UnsetStereo = new_instancemethod(_openbabel.OBAtom_UnsetStereo,None,OBAtom)
OBAtom.SetInRing = new_instancemethod(_openbabel.OBAtom_SetInRing,None,OBAtom)
OBAtom.SetChiral = new_instancemethod(_openbabel.OBAtom_SetChiral,None,OBAtom)
OBAtom.ClearCoordPtr = new_instancemethod(_openbabel.OBAtom_ClearCoordPtr,None,OBAtom)
OBAtom.GetIsotope = new_instancemethod(_openbabel.OBAtom_GetIsotope,None,OBAtom)
OBAtom.GetSpinMultiplicity = new_instancemethod(_openbabel.OBAtom_GetSpinMultiplicity,None,OBAtom)
OBAtom.GetAtomicMass = new_instancemethod(_openbabel.OBAtom_GetAtomicMass,None,OBAtom)
OBAtom.GetExactMass = new_instancemethod(_openbabel.OBAtom_GetExactMass,None,OBAtom)
OBAtom.GetCoordinateIdx = new_instancemethod(_openbabel.OBAtom_GetCoordinateIdx,None,OBAtom)
OBAtom.GetCIdx = new_instancemethod(_openbabel.OBAtom_GetCIdx,None,OBAtom)
OBAtom.GetHeteroValence = new_instancemethod(_openbabel.OBAtom_GetHeteroValence,None,OBAtom)
OBAtom.GetType = new_instancemethod(_openbabel.OBAtom_GetType,None,OBAtom)
OBAtom.GetX = new_instancemethod(_openbabel.OBAtom_GetX,None,OBAtom)
OBAtom.GetY = new_instancemethod(_openbabel.OBAtom_GetY,None,OBAtom)
OBAtom.GetZ = new_instancemethod(_openbabel.OBAtom_GetZ,None,OBAtom)
OBAtom.x = new_instancemethod(_openbabel.OBAtom_x,None,OBAtom)
OBAtom.y = new_instancemethod(_openbabel.OBAtom_y,None,OBAtom)
OBAtom.z = new_instancemethod(_openbabel.OBAtom_z,None,OBAtom)
OBAtom.GetCoordinate = new_instancemethod(_openbabel.OBAtom_GetCoordinate,None,OBAtom)
OBAtom.GetVector = new_instancemethod(_openbabel.OBAtom_GetVector,None,OBAtom)
OBAtom.GetPartialCharge = new_instancemethod(_openbabel.OBAtom_GetPartialCharge,None,OBAtom)
OBAtom.GetResidue = new_instancemethod(_openbabel.OBAtom_GetResidue,None,OBAtom)
OBAtom.GetNewBondVector = new_instancemethod(_openbabel.OBAtom_GetNewBondVector,None,OBAtom)
OBAtom.GetBond = new_instancemethod(_openbabel.OBAtom_GetBond,None,OBAtom)
OBAtom.GetNextAtom = new_instancemethod(_openbabel.OBAtom_GetNextAtom,None,OBAtom)
OBAtom.BeginBonds = new_instancemethod(_openbabel.OBAtom_BeginBonds,None,OBAtom)
OBAtom.EndBonds = new_instancemethod(_openbabel.OBAtom_EndBonds,None,OBAtom)
OBAtom.BeginBond = new_instancemethod(_openbabel.OBAtom_BeginBond,None,OBAtom)
OBAtom.NextBond = new_instancemethod(_openbabel.OBAtom_NextBond,None,OBAtom)
OBAtom.BeginNbrAtom = new_instancemethod(_openbabel.OBAtom_BeginNbrAtom,None,OBAtom)
OBAtom.NextNbrAtom = new_instancemethod(_openbabel.OBAtom_NextNbrAtom,None,OBAtom)
OBAtom.GetDistance = new_instancemethod(_openbabel.OBAtom_GetDistance,None,OBAtom)
OBAtom.GetAngle = new_instancemethod(_openbabel.OBAtom_GetAngle,None,OBAtom)
OBAtom.NewResidue = new_instancemethod(_openbabel.OBAtom_NewResidue,None,OBAtom)
OBAtom.DeleteResidue = new_instancemethod(_openbabel.OBAtom_DeleteResidue,None,OBAtom)
OBAtom.AddBond = new_instancemethod(_openbabel.OBAtom_AddBond,None,OBAtom)
OBAtom.InsertBond = new_instancemethod(_openbabel.OBAtom_InsertBond,None,OBAtom)
OBAtom.DeleteBond = new_instancemethod(_openbabel.OBAtom_DeleteBond,None,OBAtom)
OBAtom.ClearBond = new_instancemethod(_openbabel.OBAtom_ClearBond,None,OBAtom)
OBAtom.CountFreeOxygens = new_instancemethod(_openbabel.OBAtom_CountFreeOxygens,None,OBAtom)
OBAtom.MemberOfRingSize = new_instancemethod(_openbabel.OBAtom_MemberOfRingSize,None,OBAtom)
OBAtom.SmallestBondAngle = new_instancemethod(_openbabel.OBAtom_SmallestBondAngle,None,OBAtom)
OBAtom.AverageBondAngle = new_instancemethod(_openbabel.OBAtom_AverageBondAngle,None,OBAtom)
OBAtom.BOSum = new_instancemethod(_openbabel.OBAtom_BOSum,None,OBAtom)
OBAtom.HtoMethyl = new_instancemethod(_openbabel.OBAtom_HtoMethyl,None,OBAtom)
OBAtom.SetHybAndGeom = new_instancemethod(_openbabel.OBAtom_SetHybAndGeom,None,OBAtom)
OBAtom.HasResidue = new_instancemethod(_openbabel.OBAtom_HasResidue,None,OBAtom)
OBAtom.IsHydrogen = new_instancemethod(_openbabel.OBAtom_IsHydrogen,None,OBAtom)
OBAtom.IsCarbon = new_instancemethod(_openbabel.OBAtom_IsCarbon,None,OBAtom)
OBAtom.IsNitrogen = new_instancemethod(_openbabel.OBAtom_IsNitrogen,None,OBAtom)
OBAtom.IsOxygen = new_instancemethod(_openbabel.OBAtom_IsOxygen,None,OBAtom)
OBAtom.IsSulfur = new_instancemethod(_openbabel.OBAtom_IsSulfur,None,OBAtom)
OBAtom.IsPhosphorus = new_instancemethod(_openbabel.OBAtom_IsPhosphorus,None,OBAtom)
OBAtom.IsHeteroatom = new_instancemethod(_openbabel.OBAtom_IsHeteroatom,None,OBAtom)
OBAtom.IsNotCorH = new_instancemethod(_openbabel.OBAtom_IsNotCorH,None,OBAtom)
OBAtom.IsConnected = new_instancemethod(_openbabel.OBAtom_IsConnected,None,OBAtom)
OBAtom.IsOneThree = new_instancemethod(_openbabel.OBAtom_IsOneThree,None,OBAtom)
OBAtom.IsOneFour = new_instancemethod(_openbabel.OBAtom_IsOneFour,None,OBAtom)
OBAtom.IsCarboxylOxygen = new_instancemethod(_openbabel.OBAtom_IsCarboxylOxygen,None,OBAtom)
OBAtom.IsPhosphateOxygen = new_instancemethod(_openbabel.OBAtom_IsPhosphateOxygen,None,OBAtom)
OBAtom.IsSulfateOxygen = new_instancemethod(_openbabel.OBAtom_IsSulfateOxygen,None,OBAtom)
OBAtom.IsNitroOxygen = new_instancemethod(_openbabel.OBAtom_IsNitroOxygen,None,OBAtom)
OBAtom.IsAmideNitrogen = new_instancemethod(_openbabel.OBAtom_IsAmideNitrogen,None,OBAtom)
OBAtom.IsPolarHydrogen = new_instancemethod(_openbabel.OBAtom_IsPolarHydrogen,None,OBAtom)
OBAtom.IsNonPolarHydrogen = new_instancemethod(_openbabel.OBAtom_IsNonPolarHydrogen,None,OBAtom)
OBAtom.IsAromaticNOxide = new_instancemethod(_openbabel.OBAtom_IsAromaticNOxide,None,OBAtom)
OBAtom.IsChiral = new_instancemethod(_openbabel.OBAtom_IsChiral,None,OBAtom)
OBAtom.IsAxial = new_instancemethod(_openbabel.OBAtom_IsAxial,None,OBAtom)
OBAtom.IsClockwise = new_instancemethod(_openbabel.OBAtom_IsClockwise,None,OBAtom)
OBAtom.IsAntiClockwise = new_instancemethod(_openbabel.OBAtom_IsAntiClockwise,None,OBAtom)
OBAtom.IsPositiveStereo = new_instancemethod(_openbabel.OBAtom_IsPositiveStereo,None,OBAtom)
OBAtom.IsNegativeStereo = new_instancemethod(_openbabel.OBAtom_IsNegativeStereo,None,OBAtom)
OBAtom.HasChiralitySpecified = new_instancemethod(_openbabel.OBAtom_HasChiralitySpecified,None,OBAtom)
OBAtom.HasChiralVolume = new_instancemethod(_openbabel.OBAtom_HasChiralVolume,None,OBAtom)
OBAtom.IsHbondAcceptor = new_instancemethod(_openbabel.OBAtom_IsHbondAcceptor,None,OBAtom)
OBAtom.IsHbondDonor = new_instancemethod(_openbabel.OBAtom_IsHbondDonor,None,OBAtom)
OBAtom.IsHbondDonorH = new_instancemethod(_openbabel.OBAtom_IsHbondDonorH,None,OBAtom)
OBAtom.HasAlphaBetaUnsat = new_instancemethod(_openbabel.OBAtom_HasAlphaBetaUnsat,None,OBAtom)
OBAtom.HasBondOfOrder = new_instancemethod(_openbabel.OBAtom_HasBondOfOrder,None,OBAtom)
OBAtom.CountBondsOfOrder = new_instancemethod(_openbabel.OBAtom_CountBondsOfOrder,None,OBAtom)
OBAtom.HasNonSingleBond = new_instancemethod(_openbabel.OBAtom_HasNonSingleBond,None,OBAtom)
OBAtom.HasSingleBond = new_instancemethod(_openbabel.OBAtom_HasSingleBond,None,OBAtom)
OBAtom.HasDoubleBond = new_instancemethod(_openbabel.OBAtom_HasDoubleBond,None,OBAtom)
OBAtom.HasAromaticBond = new_instancemethod(_openbabel.OBAtom_HasAromaticBond,None,OBAtom)
OBAtom.MatchesSMARTS = new_instancemethod(_openbabel.OBAtom_MatchesSMARTS,None,OBAtom)
OBAtom.HasData = new_instancemethod(_openbabel.OBAtom_HasData,None,OBAtom)
OBAtom.DeleteData = new_instancemethod(_openbabel.OBAtom_DeleteData,None,OBAtom)
OBAtom.SetData = new_instancemethod(_openbabel.OBAtom_SetData,None,OBAtom)
OBAtom.DataSize = new_instancemethod(_openbabel.OBAtom_DataSize,None,OBAtom)
OBAtom.GetData = new_instancemethod(_openbabel.OBAtom_GetData,None,OBAtom)
OBAtom.BeginData = new_instancemethod(_openbabel.OBAtom_BeginData,None,OBAtom)
OBAtom.EndData = new_instancemethod(_openbabel.OBAtom_EndData,None,OBAtom)
OBAtom_swigregister = _openbabel.OBAtom_swigregister
OBAtom_swigregister(OBAtom)

OB_AROMATIC_BOND = _openbabel.OB_AROMATIC_BOND
OB_WEDGE_BOND = _openbabel.OB_WEDGE_BOND
OB_HASH_BOND = _openbabel.OB_HASH_BOND
OB_RING_BOND = _openbabel.OB_RING_BOND
OB_TORUP_BOND = _openbabel.OB_TORUP_BOND
OB_TORDOWN_BOND = _openbabel.OB_TORDOWN_BOND
OB_KSINGLE_BOND = _openbabel.OB_KSINGLE_BOND
OB_KDOUBLE_BOND = _openbabel.OB_KDOUBLE_BOND
OB_KTRIPLE_BOND = _openbabel.OB_KTRIPLE_BOND
OB_CLOSURE_BOND = _openbabel.OB_CLOSURE_BOND
class OBBond(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBBond_swiginit(self,_openbabel.new_OBBond(*args))
    __swig_destroy__ = _openbabel.delete_OBBond
OBBond.SetIdx = new_instancemethod(_openbabel.OBBond_SetIdx,None,OBBond)
OBBond.SetBO = new_instancemethod(_openbabel.OBBond_SetBO,None,OBBond)
OBBond.SetBegin = new_instancemethod(_openbabel.OBBond_SetBegin,None,OBBond)
OBBond.SetEnd = new_instancemethod(_openbabel.OBBond_SetEnd,None,OBBond)
OBBond.SetLength = new_instancemethod(_openbabel.OBBond_SetLength,None,OBBond)
OBBond.Set = new_instancemethod(_openbabel.OBBond_Set,None,OBBond)
OBBond.SetKSingle = new_instancemethod(_openbabel.OBBond_SetKSingle,None,OBBond)
OBBond.SetKDouble = new_instancemethod(_openbabel.OBBond_SetKDouble,None,OBBond)
OBBond.SetKTriple = new_instancemethod(_openbabel.OBBond_SetKTriple,None,OBBond)
OBBond.SetAromatic = new_instancemethod(_openbabel.OBBond_SetAromatic,None,OBBond)
OBBond.SetHash = new_instancemethod(_openbabel.OBBond_SetHash,None,OBBond)
OBBond.SetWedge = new_instancemethod(_openbabel.OBBond_SetWedge,None,OBBond)
OBBond.SetUp = new_instancemethod(_openbabel.OBBond_SetUp,None,OBBond)
OBBond.SetDown = new_instancemethod(_openbabel.OBBond_SetDown,None,OBBond)
OBBond.SetInRing = new_instancemethod(_openbabel.OBBond_SetInRing,None,OBBond)
OBBond.UnsetAromatic = new_instancemethod(_openbabel.OBBond_UnsetAromatic,None,OBBond)
OBBond.UnsetKekule = new_instancemethod(_openbabel.OBBond_UnsetKekule,None,OBBond)
OBBond.GetBondOrder = new_instancemethod(_openbabel.OBBond_GetBondOrder,None,OBBond)
OBBond.GetFlags = new_instancemethod(_openbabel.OBBond_GetFlags,None,OBBond)
OBBond.GetBeginAtomIdx = new_instancemethod(_openbabel.OBBond_GetBeginAtomIdx,None,OBBond)
OBBond.GetEndAtomIdx = new_instancemethod(_openbabel.OBBond_GetEndAtomIdx,None,OBBond)
OBBond.GetBeginAtom = new_instancemethod(_openbabel.OBBond_GetBeginAtom,None,OBBond)
OBBond.GetEndAtom = new_instancemethod(_openbabel.OBBond_GetEndAtom,None,OBBond)
OBBond.GetNbrAtom = new_instancemethod(_openbabel.OBBond_GetNbrAtom,None,OBBond)
OBBond.GetEquibLength = new_instancemethod(_openbabel.OBBond_GetEquibLength,None,OBBond)
OBBond.GetLength = new_instancemethod(_openbabel.OBBond_GetLength,None,OBBond)
OBBond.GetNbrAtomIdx = new_instancemethod(_openbabel.OBBond_GetNbrAtomIdx,None,OBBond)
OBBond.IsRotor = new_instancemethod(_openbabel.OBBond_IsRotor,None,OBBond)
OBBond.IsAmide = new_instancemethod(_openbabel.OBBond_IsAmide,None,OBBond)
OBBond.IsPrimaryAmide = new_instancemethod(_openbabel.OBBond_IsPrimaryAmide,None,OBBond)
OBBond.IsSecondaryAmide = new_instancemethod(_openbabel.OBBond_IsSecondaryAmide,None,OBBond)
OBBond.IsEster = new_instancemethod(_openbabel.OBBond_IsEster,None,OBBond)
OBBond.IsCarbonyl = new_instancemethod(_openbabel.OBBond_IsCarbonyl,None,OBBond)
OBBond.IsSingle = new_instancemethod(_openbabel.OBBond_IsSingle,None,OBBond)
OBBond.IsDouble = new_instancemethod(_openbabel.OBBond_IsDouble,None,OBBond)
OBBond.IsTriple = new_instancemethod(_openbabel.OBBond_IsTriple,None,OBBond)
OBBond.IsKSingle = new_instancemethod(_openbabel.OBBond_IsKSingle,None,OBBond)
OBBond.IsKDouble = new_instancemethod(_openbabel.OBBond_IsKDouble,None,OBBond)
OBBond.IsKTriple = new_instancemethod(_openbabel.OBBond_IsKTriple,None,OBBond)
OBBond.IsUp = new_instancemethod(_openbabel.OBBond_IsUp,None,OBBond)
OBBond.IsDown = new_instancemethod(_openbabel.OBBond_IsDown,None,OBBond)
OBBond.IsWedge = new_instancemethod(_openbabel.OBBond_IsWedge,None,OBBond)
OBBond.IsHash = new_instancemethod(_openbabel.OBBond_IsHash,None,OBBond)
OBBond.IsDoubleBondGeometry = new_instancemethod(_openbabel.OBBond_IsDoubleBondGeometry,None,OBBond)
OBBond.HasData = new_instancemethod(_openbabel.OBBond_HasData,None,OBBond)
OBBond.DeleteData = new_instancemethod(_openbabel.OBBond_DeleteData,None,OBBond)
OBBond.SetData = new_instancemethod(_openbabel.OBBond_SetData,None,OBBond)
OBBond.DataSize = new_instancemethod(_openbabel.OBBond_DataSize,None,OBBond)
OBBond.GetData = new_instancemethod(_openbabel.OBBond_GetData,None,OBBond)
OBBond.BeginData = new_instancemethod(_openbabel.OBBond_BeginData,None,OBBond)
OBBond.EndData = new_instancemethod(_openbabel.OBBond_EndData,None,OBBond)
OBBond_swigregister = _openbabel.OBBond_swigregister
OBBond_swigregister(OBBond)

OB_SSSR_MOL = _openbabel.OB_SSSR_MOL
OB_RINGFLAGS_MOL = _openbabel.OB_RINGFLAGS_MOL
OB_AROMATIC_MOL = _openbabel.OB_AROMATIC_MOL
OB_ATOMTYPES_MOL = _openbabel.OB_ATOMTYPES_MOL
OB_CHIRALITY_MOL = _openbabel.OB_CHIRALITY_MOL
OB_PCHARGE_MOL = _openbabel.OB_PCHARGE_MOL
OB_HYBRID_MOL = _openbabel.OB_HYBRID_MOL
OB_IMPVAL_MOL = _openbabel.OB_IMPVAL_MOL
OB_KEKULE_MOL = _openbabel.OB_KEKULE_MOL
OB_CLOSURE_MOL = _openbabel.OB_CLOSURE_MOL
OB_H_ADDED_MOL = _openbabel.OB_H_ADDED_MOL
OB_PH_CORRECTED_MOL = _openbabel.OB_PH_CORRECTED_MOL
OB_AROM_CORRECTED_MOL = _openbabel.OB_AROM_CORRECTED_MOL
OB_CHAINS_MOL = _openbabel.OB_CHAINS_MOL
OB_TCHARGE_MOL = _openbabel.OB_TCHARGE_MOL
OB_TSPIN_MOL = _openbabel.OB_TSPIN_MOL
OB_CURRENT_CONFORMER = _openbabel.OB_CURRENT_CONFORMER
class OBMol(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBMol_swiginit(self,_openbabel.new_OBMol(*args))
    __swig_destroy__ = _openbabel.delete_OBMol
    ClassDescription = staticmethod(_openbabel.OBMol_ClassDescription)
OBMol.__iadd__ = new_instancemethod(_openbabel.OBMol___iadd__,None,OBMol)
OBMol.ReserveAtoms = new_instancemethod(_openbabel.OBMol_ReserveAtoms,None,OBMol)
OBMol.CreateAtom = new_instancemethod(_openbabel.OBMol_CreateAtom,None,OBMol)
OBMol.CreateBond = new_instancemethod(_openbabel.OBMol_CreateBond,None,OBMol)
OBMol.DestroyAtom = new_instancemethod(_openbabel.OBMol_DestroyAtom,None,OBMol)
OBMol.DestroyBond = new_instancemethod(_openbabel.OBMol_DestroyBond,None,OBMol)
OBMol.AddAtom = new_instancemethod(_openbabel.OBMol_AddAtom,None,OBMol)
OBMol.AddBond = new_instancemethod(_openbabel.OBMol_AddBond,None,OBMol)
OBMol.AddResidue = new_instancemethod(_openbabel.OBMol_AddResidue,None,OBMol)
OBMol.InsertAtom = new_instancemethod(_openbabel.OBMol_InsertAtom,None,OBMol)
OBMol.DeleteAtom = new_instancemethod(_openbabel.OBMol_DeleteAtom,None,OBMol)
OBMol.DeleteBond = new_instancemethod(_openbabel.OBMol_DeleteBond,None,OBMol)
OBMol.DeleteResidue = new_instancemethod(_openbabel.OBMol_DeleteResidue,None,OBMol)
OBMol.NewAtom = new_instancemethod(_openbabel.OBMol_NewAtom,None,OBMol)
OBMol.NewResidue = new_instancemethod(_openbabel.OBMol_NewResidue,None,OBMol)
OBMol.BeginModify = new_instancemethod(_openbabel.OBMol_BeginModify,None,OBMol)
OBMol.EndModify = new_instancemethod(_openbabel.OBMol_EndModify,None,OBMol)
OBMol.GetMod = new_instancemethod(_openbabel.OBMol_GetMod,None,OBMol)
OBMol.IncrementMod = new_instancemethod(_openbabel.OBMol_IncrementMod,None,OBMol)
OBMol.DecrementMod = new_instancemethod(_openbabel.OBMol_DecrementMod,None,OBMol)
OBMol.HasData = new_instancemethod(_openbabel.OBMol_HasData,None,OBMol)
OBMol.DeleteData = new_instancemethod(_openbabel.OBMol_DeleteData,None,OBMol)
OBMol.SetData = new_instancemethod(_openbabel.OBMol_SetData,None,OBMol)
OBMol.DataSize = new_instancemethod(_openbabel.OBMol_DataSize,None,OBMol)
OBMol.GetData = new_instancemethod(_openbabel.OBMol_GetData,None,OBMol)
OBMol.BeginData = new_instancemethod(_openbabel.OBMol_BeginData,None,OBMol)
OBMol.EndData = new_instancemethod(_openbabel.OBMol_EndData,None,OBMol)
OBMol.GetFlags = new_instancemethod(_openbabel.OBMol_GetFlags,None,OBMol)
OBMol.GetTitle = new_instancemethod(_openbabel.OBMol_GetTitle,None,OBMol)
OBMol.NumAtoms = new_instancemethod(_openbabel.OBMol_NumAtoms,None,OBMol)
OBMol.NumBonds = new_instancemethod(_openbabel.OBMol_NumBonds,None,OBMol)
OBMol.NumHvyAtoms = new_instancemethod(_openbabel.OBMol_NumHvyAtoms,None,OBMol)
OBMol.NumResidues = new_instancemethod(_openbabel.OBMol_NumResidues,None,OBMol)
OBMol.NumRotors = new_instancemethod(_openbabel.OBMol_NumRotors,None,OBMol)
OBMol.GetAtom = new_instancemethod(_openbabel.OBMol_GetAtom,None,OBMol)
OBMol.GetFirstAtom = new_instancemethod(_openbabel.OBMol_GetFirstAtom,None,OBMol)
OBMol.GetBond = new_instancemethod(_openbabel.OBMol_GetBond,None,OBMol)
OBMol.GetResidue = new_instancemethod(_openbabel.OBMol_GetResidue,None,OBMol)
OBMol.GetInternalCoord = new_instancemethod(_openbabel.OBMol_GetInternalCoord,None,OBMol)
OBMol.GetTorsion = new_instancemethod(_openbabel.OBMol_GetTorsion,None,OBMol)
OBMol.GetFormula = new_instancemethod(_openbabel.OBMol_GetFormula,None,OBMol)
OBMol.GetEnergy = new_instancemethod(_openbabel.OBMol_GetEnergy,None,OBMol)
OBMol.GetMolWt = new_instancemethod(_openbabel.OBMol_GetMolWt,None,OBMol)
OBMol.GetExactMass = new_instancemethod(_openbabel.OBMol_GetExactMass,None,OBMol)
OBMol.GetTotalCharge = new_instancemethod(_openbabel.OBMol_GetTotalCharge,None,OBMol)
OBMol.GetTotalSpinMultiplicity = new_instancemethod(_openbabel.OBMol_GetTotalSpinMultiplicity,None,OBMol)
OBMol.GetDimension = new_instancemethod(_openbabel.OBMol_GetDimension,None,OBMol)
OBMol.GetCoordinates = new_instancemethod(_openbabel.OBMol_GetCoordinates,None,OBMol)
OBMol.GetSSSR = new_instancemethod(_openbabel.OBMol_GetSSSR,None,OBMol)
OBMol.AutomaticFormalCharge = new_instancemethod(_openbabel.OBMol_AutomaticFormalCharge,None,OBMol)
OBMol.AutomaticPartialCharge = new_instancemethod(_openbabel.OBMol_AutomaticPartialCharge,None,OBMol)
OBMol.SetTitle = new_instancemethod(_openbabel.OBMol_SetTitle,None,OBMol)
OBMol.SetFormula = new_instancemethod(_openbabel.OBMol_SetFormula,None,OBMol)
OBMol.SetEnergy = new_instancemethod(_openbabel.OBMol_SetEnergy,None,OBMol)
OBMol.SetDimension = new_instancemethod(_openbabel.OBMol_SetDimension,None,OBMol)
OBMol.SetTotalCharge = new_instancemethod(_openbabel.OBMol_SetTotalCharge,None,OBMol)
OBMol.SetTotalSpinMultiplicity = new_instancemethod(_openbabel.OBMol_SetTotalSpinMultiplicity,None,OBMol)
OBMol.SetInternalCoord = new_instancemethod(_openbabel.OBMol_SetInternalCoord,None,OBMol)
OBMol.SetAutomaticFormalCharge = new_instancemethod(_openbabel.OBMol_SetAutomaticFormalCharge,None,OBMol)
OBMol.SetAutomaticPartialCharge = new_instancemethod(_openbabel.OBMol_SetAutomaticPartialCharge,None,OBMol)
OBMol.SetAromaticPerceived = new_instancemethod(_openbabel.OBMol_SetAromaticPerceived,None,OBMol)
OBMol.SetSSSRPerceived = new_instancemethod(_openbabel.OBMol_SetSSSRPerceived,None,OBMol)
OBMol.SetRingAtomsAndBondsPerceived = new_instancemethod(_openbabel.OBMol_SetRingAtomsAndBondsPerceived,None,OBMol)
OBMol.SetAtomTypesPerceived = new_instancemethod(_openbabel.OBMol_SetAtomTypesPerceived,None,OBMol)
OBMol.SetChainsPerceived = new_instancemethod(_openbabel.OBMol_SetChainsPerceived,None,OBMol)
OBMol.SetChiralityPerceived = new_instancemethod(_openbabel.OBMol_SetChiralityPerceived,None,OBMol)
OBMol.SetPartialChargesPerceived = new_instancemethod(_openbabel.OBMol_SetPartialChargesPerceived,None,OBMol)
OBMol.SetHybridizationPerceived = new_instancemethod(_openbabel.OBMol_SetHybridizationPerceived,None,OBMol)
OBMol.SetImplicitValencePerceived = new_instancemethod(_openbabel.OBMol_SetImplicitValencePerceived,None,OBMol)
OBMol.SetKekulePerceived = new_instancemethod(_openbabel.OBMol_SetKekulePerceived,None,OBMol)
OBMol.SetClosureBondsPerceived = new_instancemethod(_openbabel.OBMol_SetClosureBondsPerceived,None,OBMol)
OBMol.SetHydrogensAdded = new_instancemethod(_openbabel.OBMol_SetHydrogensAdded,None,OBMol)
OBMol.SetCorrectedForPH = new_instancemethod(_openbabel.OBMol_SetCorrectedForPH,None,OBMol)
OBMol.SetAromaticCorrected = new_instancemethod(_openbabel.OBMol_SetAromaticCorrected,None,OBMol)
OBMol.SetSpinMultiplicityAssigned = new_instancemethod(_openbabel.OBMol_SetSpinMultiplicityAssigned,None,OBMol)
OBMol.SetFlags = new_instancemethod(_openbabel.OBMol_SetFlags,None,OBMol)
OBMol.UnsetAromaticPerceived = new_instancemethod(_openbabel.OBMol_UnsetAromaticPerceived,None,OBMol)
OBMol.UnsetPartialChargesPerceived = new_instancemethod(_openbabel.OBMol_UnsetPartialChargesPerceived,None,OBMol)
OBMol.UnsetImplicitValencePerceived = new_instancemethod(_openbabel.OBMol_UnsetImplicitValencePerceived,None,OBMol)
OBMol.UnsetFlag = new_instancemethod(_openbabel.OBMol_UnsetFlag,None,OBMol)
OBMol.Clear = new_instancemethod(_openbabel.OBMol_Clear,None,OBMol)
OBMol.RenumberAtoms = new_instancemethod(_openbabel.OBMol_RenumberAtoms,None,OBMol)
OBMol.ToInertialFrame = new_instancemethod(_openbabel.OBMol_ToInertialFrame,None,OBMol)
OBMol.Translate = new_instancemethod(_openbabel.OBMol_Translate,None,OBMol)
OBMol.Rotate = new_instancemethod(_openbabel.OBMol_Rotate,None,OBMol)
OBMol.Kekulize = new_instancemethod(_openbabel.OBMol_Kekulize,None,OBMol)
OBMol.PerceiveKekuleBonds = new_instancemethod(_openbabel.OBMol_PerceiveKekuleBonds,None,OBMol)
OBMol.NewPerceiveKekuleBonds = new_instancemethod(_openbabel.OBMol_NewPerceiveKekuleBonds,None,OBMol)
OBMol.DeleteHydrogen = new_instancemethod(_openbabel.OBMol_DeleteHydrogen,None,OBMol)
OBMol.DeleteHydrogens = new_instancemethod(_openbabel.OBMol_DeleteHydrogens,None,OBMol)
OBMol.DeleteNonPolarHydrogens = new_instancemethod(_openbabel.OBMol_DeleteNonPolarHydrogens,None,OBMol)
OBMol.AddHydrogens = new_instancemethod(_openbabel.OBMol_AddHydrogens,None,OBMol)
OBMol.AddPolarHydrogens = new_instancemethod(_openbabel.OBMol_AddPolarHydrogens,None,OBMol)
OBMol.StripSalts = new_instancemethod(_openbabel.OBMol_StripSalts,None,OBMol)
OBMol.ConvertDativeBonds = new_instancemethod(_openbabel.OBMol_ConvertDativeBonds,None,OBMol)
OBMol.CorrectForPH = new_instancemethod(_openbabel.OBMol_CorrectForPH,None,OBMol)
OBMol.AssignSpinMultiplicity = new_instancemethod(_openbabel.OBMol_AssignSpinMultiplicity,None,OBMol)
OBMol.Center = new_instancemethod(_openbabel.OBMol_Center,None,OBMol)
OBMol.SetTorsion = new_instancemethod(_openbabel.OBMol_SetTorsion,None,OBMol)
OBMol.FindSSSR = new_instancemethod(_openbabel.OBMol_FindSSSR,None,OBMol)
OBMol.FindRingAtomsAndBonds = new_instancemethod(_openbabel.OBMol_FindRingAtomsAndBonds,None,OBMol)
OBMol.FindChiralCenters = new_instancemethod(_openbabel.OBMol_FindChiralCenters,None,OBMol)
OBMol.FindChildren = new_instancemethod(_openbabel.OBMol_FindChildren,None,OBMol)
OBMol.FindLargestFragment = new_instancemethod(_openbabel.OBMol_FindLargestFragment,None,OBMol)
OBMol.ContigFragList = new_instancemethod(_openbabel.OBMol_ContigFragList,None,OBMol)
OBMol.Align = new_instancemethod(_openbabel.OBMol_Align,None,OBMol)
OBMol.ConnectTheDots = new_instancemethod(_openbabel.OBMol_ConnectTheDots,None,OBMol)
OBMol.PerceiveBondOrders = new_instancemethod(_openbabel.OBMol_PerceiveBondOrders,None,OBMol)
OBMol.FindTorsions = new_instancemethod(_openbabel.OBMol_FindTorsions,None,OBMol)
OBMol.GetGTDVector = new_instancemethod(_openbabel.OBMol_GetGTDVector,None,OBMol)
OBMol.GetGIVector = new_instancemethod(_openbabel.OBMol_GetGIVector,None,OBMol)
OBMol.GetGIDVector = new_instancemethod(_openbabel.OBMol_GetGIDVector,None,OBMol)
OBMol.Has2D = new_instancemethod(_openbabel.OBMol_Has2D,None,OBMol)
OBMol.Has3D = new_instancemethod(_openbabel.OBMol_Has3D,None,OBMol)
OBMol.HasNonZeroCoords = new_instancemethod(_openbabel.OBMol_HasNonZeroCoords,None,OBMol)
OBMol.HasAromaticPerceived = new_instancemethod(_openbabel.OBMol_HasAromaticPerceived,None,OBMol)
OBMol.HasSSSRPerceived = new_instancemethod(_openbabel.OBMol_HasSSSRPerceived,None,OBMol)
OBMol.HasRingAtomsAndBondsPerceived = new_instancemethod(_openbabel.OBMol_HasRingAtomsAndBondsPerceived,None,OBMol)
OBMol.HasAtomTypesPerceived = new_instancemethod(_openbabel.OBMol_HasAtomTypesPerceived,None,OBMol)
OBMol.HasChiralityPerceived = new_instancemethod(_openbabel.OBMol_HasChiralityPerceived,None,OBMol)
OBMol.HasPartialChargesPerceived = new_instancemethod(_openbabel.OBMol_HasPartialChargesPerceived,None,OBMol)
OBMol.HasHybridizationPerceived = new_instancemethod(_openbabel.OBMol_HasHybridizationPerceived,None,OBMol)
OBMol.HasImplicitValencePerceived = new_instancemethod(_openbabel.OBMol_HasImplicitValencePerceived,None,OBMol)
OBMol.HasKekulePerceived = new_instancemethod(_openbabel.OBMol_HasKekulePerceived,None,OBMol)
OBMol.HasClosureBondsPerceived = new_instancemethod(_openbabel.OBMol_HasClosureBondsPerceived,None,OBMol)
OBMol.HasChainsPerceived = new_instancemethod(_openbabel.OBMol_HasChainsPerceived,None,OBMol)
OBMol.HasHydrogensAdded = new_instancemethod(_openbabel.OBMol_HasHydrogensAdded,None,OBMol)
OBMol.HasAromaticCorrected = new_instancemethod(_openbabel.OBMol_HasAromaticCorrected,None,OBMol)
OBMol.IsCorrectedForPH = new_instancemethod(_openbabel.OBMol_IsCorrectedForPH,None,OBMol)
OBMol.HasSpinMultiplicityAssigned = new_instancemethod(_openbabel.OBMol_HasSpinMultiplicityAssigned,None,OBMol)
OBMol.IsChiral = new_instancemethod(_openbabel.OBMol_IsChiral,None,OBMol)
OBMol.Empty = new_instancemethod(_openbabel.OBMol_Empty,None,OBMol)
OBMol.NumConformers = new_instancemethod(_openbabel.OBMol_NumConformers,None,OBMol)
OBMol.SetConformers = new_instancemethod(_openbabel.OBMol_SetConformers,None,OBMol)
OBMol.AddConformer = new_instancemethod(_openbabel.OBMol_AddConformer,None,OBMol)
OBMol.SetConformer = new_instancemethod(_openbabel.OBMol_SetConformer,None,OBMol)
OBMol.CopyConformer = new_instancemethod(_openbabel.OBMol_CopyConformer,None,OBMol)
OBMol.DeleteConformer = new_instancemethod(_openbabel.OBMol_DeleteConformer,None,OBMol)
OBMol.GetConformer = new_instancemethod(_openbabel.OBMol_GetConformer,None,OBMol)
OBMol.BeginConformer = new_instancemethod(_openbabel.OBMol_BeginConformer,None,OBMol)
OBMol.NextConformer = new_instancemethod(_openbabel.OBMol_NextConformer,None,OBMol)
OBMol.GetConformers = new_instancemethod(_openbabel.OBMol_GetConformers,None,OBMol)
OBMol.BeginAtom = new_instancemethod(_openbabel.OBMol_BeginAtom,None,OBMol)
OBMol.NextAtom = new_instancemethod(_openbabel.OBMol_NextAtom,None,OBMol)
OBMol.BeginBond = new_instancemethod(_openbabel.OBMol_BeginBond,None,OBMol)
OBMol.NextBond = new_instancemethod(_openbabel.OBMol_NextBond,None,OBMol)
OBMol.BeginResidue = new_instancemethod(_openbabel.OBMol_BeginResidue,None,OBMol)
OBMol.NextResidue = new_instancemethod(_openbabel.OBMol_NextResidue,None,OBMol)
OBMol.BeginInternalCoord = new_instancemethod(_openbabel.OBMol_BeginInternalCoord,None,OBMol)
OBMol.NextInternalCoord = new_instancemethod(_openbabel.OBMol_NextInternalCoord,None,OBMol)
OBMol_swigregister = _openbabel.OBMol_swigregister
OBMol_swigregister(OBMol)
OBMol_ClassDescription = _openbabel.OBMol_ClassDescription

class OBInternalCoord(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    _a = property(_openbabel.OBInternalCoord__a_get, _openbabel.OBInternalCoord__a_set)
    _b = property(_openbabel.OBInternalCoord__b_get, _openbabel.OBInternalCoord__b_set)
    _c = property(_openbabel.OBInternalCoord__c_get, _openbabel.OBInternalCoord__c_set)
    _dst = property(_openbabel.OBInternalCoord__dst_get, _openbabel.OBInternalCoord__dst_set)
    _ang = property(_openbabel.OBInternalCoord__ang_get, _openbabel.OBInternalCoord__ang_set)
    _tor = property(_openbabel.OBInternalCoord__tor_get, _openbabel.OBInternalCoord__tor_set)
    def __init__(self, *args): 
        _openbabel.OBInternalCoord_swiginit(self,_openbabel.new_OBInternalCoord(*args))
    __swig_destroy__ = _openbabel.delete_OBInternalCoord
OBInternalCoord_swigregister = _openbabel.OBInternalCoord_swigregister
OBInternalCoord_swigregister(OBInternalCoord)

CartesianToInternal = _openbabel.CartesianToInternal
InternalToCartesian = _openbabel.InternalToCartesian
NewExtension = _openbabel.NewExtension
BUFF_SIZE = _openbabel.BUFF_SIZE
get_rmat = _openbabel.get_rmat
ob_make_rmat = _openbabel.ob_make_rmat
qtrfit = _openbabel.qtrfit
superimpose = _openbabel.superimpose
class OBRTree(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBRTree_swiginit(self,_openbabel.new_OBRTree(*args))
    __swig_destroy__ = _openbabel.delete_OBRTree
OBRTree.GetAtomIdx = new_instancemethod(_openbabel.OBRTree_GetAtomIdx,None,OBRTree)
OBRTree.PathToRoot = new_instancemethod(_openbabel.OBRTree_PathToRoot,None,OBRTree)
OBRTree_swigregister = _openbabel.OBRTree_swigregister
OBRTree_swigregister(OBRTree)
tokenize = _openbabel.tokenize
ThrowError = _openbabel.ThrowError

class OBRing(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    _path = property(_openbabel.OBRing__path_get, _openbabel.OBRing__path_set)
    _pathset = property(_openbabel.OBRing__pathset_get, _openbabel.OBRing__pathset_set)
    def __init__(self, *args): 
        _openbabel.OBRing_swiginit(self,_openbabel.new_OBRing(*args))
    __swig_destroy__ = _openbabel.delete_OBRing
OBRing.findCenterAndNormal = new_instancemethod(_openbabel.OBRing_findCenterAndNormal,None,OBRing)
OBRing.Size = new_instancemethod(_openbabel.OBRing_Size,None,OBRing)
OBRing.PathSize = new_instancemethod(_openbabel.OBRing_PathSize,None,OBRing)
OBRing.IsMember = new_instancemethod(_openbabel.OBRing_IsMember,None,OBRing)
OBRing.IsAromatic = new_instancemethod(_openbabel.OBRing_IsAromatic,None,OBRing)
OBRing.IsInRing = new_instancemethod(_openbabel.OBRing_IsInRing,None,OBRing)
OBRing.SetParent = new_instancemethod(_openbabel.OBRing_SetParent,None,OBRing)
OBRing.GetParent = new_instancemethod(_openbabel.OBRing_GetParent,None,OBRing)
OBRing_swigregister = _openbabel.OBRing_swigregister
OBRing_swigregister(OBRing)

CompareRingSize = _openbabel.CompareRingSize
class OBRingSearch(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBRingSearch_swiginit(self,_openbabel.new_OBRingSearch(*args))
    __swig_destroy__ = _openbabel.delete_OBRingSearch
OBRingSearch.SortRings = new_instancemethod(_openbabel.OBRingSearch_SortRings,None,OBRingSearch)
OBRingSearch.RemoveRedundant = new_instancemethod(_openbabel.OBRingSearch_RemoveRedundant,None,OBRingSearch)
OBRingSearch.AddRingFromClosure = new_instancemethod(_openbabel.OBRingSearch_AddRingFromClosure,None,OBRingSearch)
OBRingSearch.WriteRings = new_instancemethod(_openbabel.OBRingSearch_WriteRings,None,OBRingSearch)
OBRingSearch.SaveUniqueRing = new_instancemethod(_openbabel.OBRingSearch_SaveUniqueRing,None,OBRingSearch)
OBRingSearch.BeginRings = new_instancemethod(_openbabel.OBRingSearch_BeginRings,None,OBRingSearch)
OBRingSearch.EndRings = new_instancemethod(_openbabel.OBRingSearch_EndRings,None,OBRingSearch)
OBRingSearch_swigregister = _openbabel.OBRingSearch_swigregister
OBRingSearch_swigregister(OBRingSearch)

AE_LEAF = _openbabel.AE_LEAF
AE_RECUR = _openbabel.AE_RECUR
AE_NOT = _openbabel.AE_NOT
AE_ANDHI = _openbabel.AE_ANDHI
AE_OR = _openbabel.AE_OR
AE_ANDLO = _openbabel.AE_ANDLO
AL_CONST = _openbabel.AL_CONST
AL_MASS = _openbabel.AL_MASS
AL_AROM = _openbabel.AL_AROM
AL_ELEM = _openbabel.AL_ELEM
AL_HCOUNT = _openbabel.AL_HCOUNT
AL_NEGATIVE = _openbabel.AL_NEGATIVE
AL_POSITIVE = _openbabel.AL_POSITIVE
AL_CONNECT = _openbabel.AL_CONNECT
AL_DEGREE = _openbabel.AL_DEGREE
AL_IMPLICIT = _openbabel.AL_IMPLICIT
AL_RINGS = _openbabel.AL_RINGS
AL_SIZE = _openbabel.AL_SIZE
AL_VALENCE = _openbabel.AL_VALENCE
AL_CHIRAL = _openbabel.AL_CHIRAL
AL_HYB = _openbabel.AL_HYB
AL_CLOCKWISE = _openbabel.AL_CLOCKWISE
AL_ANTICLOCKWISE = _openbabel.AL_ANTICLOCKWISE
class OBSmartsPattern(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    __swig_destroy__ = _openbabel.delete_OBSmartsPattern
    def __init__(self, *args): 
        _openbabel.OBSmartsPattern_swiginit(self,_openbabel.new_OBSmartsPattern(*args))
OBSmartsPattern.NumMatches = new_instancemethod(_openbabel.OBSmartsPattern_NumMatches,None,OBSmartsPattern)
OBSmartsPattern.NumAtoms = new_instancemethod(_openbabel.OBSmartsPattern_NumAtoms,None,OBSmartsPattern)
OBSmartsPattern.NumBonds = new_instancemethod(_openbabel.OBSmartsPattern_NumBonds,None,OBSmartsPattern)
OBSmartsPattern.GetAtomicNum = new_instancemethod(_openbabel.OBSmartsPattern_GetAtomicNum,None,OBSmartsPattern)
OBSmartsPattern.GetBond = new_instancemethod(_openbabel.OBSmartsPattern_GetBond,None,OBSmartsPattern)
OBSmartsPattern.GetCharge = new_instancemethod(_openbabel.OBSmartsPattern_GetCharge,None,OBSmartsPattern)
OBSmartsPattern.GetSMARTS = new_instancemethod(_openbabel.OBSmartsPattern_GetSMARTS,None,OBSmartsPattern)
OBSmartsPattern.GetVectorBinding = new_instancemethod(_openbabel.OBSmartsPattern_GetVectorBinding,None,OBSmartsPattern)
OBSmartsPattern.Empty = new_instancemethod(_openbabel.OBSmartsPattern_Empty,None,OBSmartsPattern)
OBSmartsPattern.IsValid = new_instancemethod(_openbabel.OBSmartsPattern_IsValid,None,OBSmartsPattern)
OBSmartsPattern.Init = new_instancemethod(_openbabel.OBSmartsPattern_Init,None,OBSmartsPattern)
OBSmartsPattern.WriteMapList = new_instancemethod(_openbabel.OBSmartsPattern_WriteMapList,None,OBSmartsPattern)
OBSmartsPattern.Match = new_instancemethod(_openbabel.OBSmartsPattern_Match,None,OBSmartsPattern)
OBSmartsPattern.RestrictedMatch = new_instancemethod(_openbabel.OBSmartsPattern_RestrictedMatch,None,OBSmartsPattern)
OBSmartsPattern.GetMapList = new_instancemethod(_openbabel.OBSmartsPattern_GetMapList,None,OBSmartsPattern)
OBSmartsPattern.GetUMapList = new_instancemethod(_openbabel.OBSmartsPattern_GetUMapList,None,OBSmartsPattern)
OBSmartsPattern.BeginMList = new_instancemethod(_openbabel.OBSmartsPattern_BeginMList,None,OBSmartsPattern)
OBSmartsPattern.EndMList = new_instancemethod(_openbabel.OBSmartsPattern_EndMList,None,OBSmartsPattern)
OBSmartsPattern_swigregister = _openbabel.OBSmartsPattern_swigregister
OBSmartsPattern_swigregister(OBSmartsPattern)

class OBSSMatch(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBSSMatch_swiginit(self,_openbabel.new_OBSSMatch(*args))
    __swig_destroy__ = _openbabel.delete_OBSSMatch
OBSSMatch.Match = new_instancemethod(_openbabel.OBSSMatch_Match,None,OBSSMatch)
OBSSMatch_swigregister = _openbabel.OBSSMatch_swigregister
OBSSMatch_swigregister(OBSSMatch)

SmartsLexReplace = _openbabel.SmartsLexReplace


