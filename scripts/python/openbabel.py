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
    __del__ = lambda self : None;
    def value(*args): return _openbabel.PySwigIterator_value(*args)
    def incr(*args): return _openbabel.PySwigIterator_incr(*args)
    def decr(*args): return _openbabel.PySwigIterator_decr(*args)
    def distance(*args): return _openbabel.PySwigIterator_distance(*args)
    def equal(*args): return _openbabel.PySwigIterator_equal(*args)
    def copy(*args): return _openbabel.PySwigIterator_copy(*args)
    def next(*args): return _openbabel.PySwigIterator_next(*args)
    def previous(*args): return _openbabel.PySwigIterator_previous(*args)
    def advance(*args): return _openbabel.PySwigIterator_advance(*args)
    def __eq__(*args): return _openbabel.PySwigIterator___eq__(*args)
    def __ne__(*args): return _openbabel.PySwigIterator___ne__(*args)
    def __iadd__(*args): return _openbabel.PySwigIterator___iadd__(*args)
    def __isub__(*args): return _openbabel.PySwigIterator___isub__(*args)
    def __add__(*args): return _openbabel.PySwigIterator___add__(*args)
    def __sub__(*args): return _openbabel.PySwigIterator___sub__(*args)
    def __iter__(self): return self
PySwigIterator_swigregister = _openbabel.PySwigIterator_swigregister
PySwigIterator_swigregister(PySwigIterator)

class vectorInt(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def iterator(*args): return _openbabel.vectorInt_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _openbabel.vectorInt___nonzero__(*args)
    def __len__(*args): return _openbabel.vectorInt___len__(*args)
    def pop(*args): return _openbabel.vectorInt_pop(*args)
    def __getslice__(*args): return _openbabel.vectorInt___getslice__(*args)
    def __setslice__(*args): return _openbabel.vectorInt___setslice__(*args)
    def __delslice__(*args): return _openbabel.vectorInt___delslice__(*args)
    def __delitem__(*args): return _openbabel.vectorInt___delitem__(*args)
    def __getitem__(*args): return _openbabel.vectorInt___getitem__(*args)
    def __setitem__(*args): return _openbabel.vectorInt___setitem__(*args)
    def append(*args): return _openbabel.vectorInt_append(*args)
    def empty(*args): return _openbabel.vectorInt_empty(*args)
    def size(*args): return _openbabel.vectorInt_size(*args)
    def clear(*args): return _openbabel.vectorInt_clear(*args)
    def swap(*args): return _openbabel.vectorInt_swap(*args)
    def get_allocator(*args): return _openbabel.vectorInt_get_allocator(*args)
    def begin(*args): return _openbabel.vectorInt_begin(*args)
    def end(*args): return _openbabel.vectorInt_end(*args)
    def rbegin(*args): return _openbabel.vectorInt_rbegin(*args)
    def rend(*args): return _openbabel.vectorInt_rend(*args)
    def pop_back(*args): return _openbabel.vectorInt_pop_back(*args)
    def erase(*args): return _openbabel.vectorInt_erase(*args)
    def __init__(self, *args): 
        _openbabel.vectorInt_swiginit(self,_openbabel.new_vectorInt(*args))
    def push_back(*args): return _openbabel.vectorInt_push_back(*args)
    def front(*args): return _openbabel.vectorInt_front(*args)
    def back(*args): return _openbabel.vectorInt_back(*args)
    def assign(*args): return _openbabel.vectorInt_assign(*args)
    def resize(*args): return _openbabel.vectorInt_resize(*args)
    def insert(*args): return _openbabel.vectorInt_insert(*args)
    def reserve(*args): return _openbabel.vectorInt_reserve(*args)
    def capacity(*args): return _openbabel.vectorInt_capacity(*args)
    __swig_destroy__ = _openbabel.delete_vectorInt
    __del__ = lambda self : None;
vectorInt_swigregister = _openbabel.vectorInt_swigregister
vectorInt_swigregister(vectorInt)

class vvInt(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def iterator(*args): return _openbabel.vvInt_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _openbabel.vvInt___nonzero__(*args)
    def __len__(*args): return _openbabel.vvInt___len__(*args)
    def pop(*args): return _openbabel.vvInt_pop(*args)
    def __getslice__(*args): return _openbabel.vvInt___getslice__(*args)
    def __setslice__(*args): return _openbabel.vvInt___setslice__(*args)
    def __delslice__(*args): return _openbabel.vvInt___delslice__(*args)
    def __delitem__(*args): return _openbabel.vvInt___delitem__(*args)
    def __getitem__(*args): return _openbabel.vvInt___getitem__(*args)
    def __setitem__(*args): return _openbabel.vvInt___setitem__(*args)
    def append(*args): return _openbabel.vvInt_append(*args)
    def empty(*args): return _openbabel.vvInt_empty(*args)
    def size(*args): return _openbabel.vvInt_size(*args)
    def clear(*args): return _openbabel.vvInt_clear(*args)
    def swap(*args): return _openbabel.vvInt_swap(*args)
    def get_allocator(*args): return _openbabel.vvInt_get_allocator(*args)
    def begin(*args): return _openbabel.vvInt_begin(*args)
    def end(*args): return _openbabel.vvInt_end(*args)
    def rbegin(*args): return _openbabel.vvInt_rbegin(*args)
    def rend(*args): return _openbabel.vvInt_rend(*args)
    def pop_back(*args): return _openbabel.vvInt_pop_back(*args)
    def erase(*args): return _openbabel.vvInt_erase(*args)
    def __init__(self, *args): 
        _openbabel.vvInt_swiginit(self,_openbabel.new_vvInt(*args))
    def push_back(*args): return _openbabel.vvInt_push_back(*args)
    def front(*args): return _openbabel.vvInt_front(*args)
    def back(*args): return _openbabel.vvInt_back(*args)
    def assign(*args): return _openbabel.vvInt_assign(*args)
    def resize(*args): return _openbabel.vvInt_resize(*args)
    def insert(*args): return _openbabel.vvInt_insert(*args)
    def reserve(*args): return _openbabel.vvInt_reserve(*args)
    def capacity(*args): return _openbabel.vvInt_capacity(*args)
    __swig_destroy__ = _openbabel.delete_vvInt
    __del__ = lambda self : None;
vvInt_swigregister = _openbabel.vvInt_swigregister
vvInt_swigregister(vvInt)

class vectorDouble(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def iterator(*args): return _openbabel.vectorDouble_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _openbabel.vectorDouble___nonzero__(*args)
    def __len__(*args): return _openbabel.vectorDouble___len__(*args)
    def pop(*args): return _openbabel.vectorDouble_pop(*args)
    def __getslice__(*args): return _openbabel.vectorDouble___getslice__(*args)
    def __setslice__(*args): return _openbabel.vectorDouble___setslice__(*args)
    def __delslice__(*args): return _openbabel.vectorDouble___delslice__(*args)
    def __delitem__(*args): return _openbabel.vectorDouble___delitem__(*args)
    def __getitem__(*args): return _openbabel.vectorDouble___getitem__(*args)
    def __setitem__(*args): return _openbabel.vectorDouble___setitem__(*args)
    def append(*args): return _openbabel.vectorDouble_append(*args)
    def empty(*args): return _openbabel.vectorDouble_empty(*args)
    def size(*args): return _openbabel.vectorDouble_size(*args)
    def clear(*args): return _openbabel.vectorDouble_clear(*args)
    def swap(*args): return _openbabel.vectorDouble_swap(*args)
    def get_allocator(*args): return _openbabel.vectorDouble_get_allocator(*args)
    def begin(*args): return _openbabel.vectorDouble_begin(*args)
    def end(*args): return _openbabel.vectorDouble_end(*args)
    def rbegin(*args): return _openbabel.vectorDouble_rbegin(*args)
    def rend(*args): return _openbabel.vectorDouble_rend(*args)
    def pop_back(*args): return _openbabel.vectorDouble_pop_back(*args)
    def erase(*args): return _openbabel.vectorDouble_erase(*args)
    def __init__(self, *args): 
        _openbabel.vectorDouble_swiginit(self,_openbabel.new_vectorDouble(*args))
    def push_back(*args): return _openbabel.vectorDouble_push_back(*args)
    def front(*args): return _openbabel.vectorDouble_front(*args)
    def back(*args): return _openbabel.vectorDouble_back(*args)
    def assign(*args): return _openbabel.vectorDouble_assign(*args)
    def resize(*args): return _openbabel.vectorDouble_resize(*args)
    def insert(*args): return _openbabel.vectorDouble_insert(*args)
    def reserve(*args): return _openbabel.vectorDouble_reserve(*args)
    def capacity(*args): return _openbabel.vectorDouble_capacity(*args)
    __swig_destroy__ = _openbabel.delete_vectorDouble
    __del__ = lambda self : None;
vectorDouble_swigregister = _openbabel.vectorDouble_swigregister
vectorDouble_swigregister(vectorDouble)

class vVector3(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def iterator(*args): return _openbabel.vVector3_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _openbabel.vVector3___nonzero__(*args)
    def __len__(*args): return _openbabel.vVector3___len__(*args)
    def pop(*args): return _openbabel.vVector3_pop(*args)
    def __getslice__(*args): return _openbabel.vVector3___getslice__(*args)
    def __setslice__(*args): return _openbabel.vVector3___setslice__(*args)
    def __delslice__(*args): return _openbabel.vVector3___delslice__(*args)
    def __delitem__(*args): return _openbabel.vVector3___delitem__(*args)
    def __getitem__(*args): return _openbabel.vVector3___getitem__(*args)
    def __setitem__(*args): return _openbabel.vVector3___setitem__(*args)
    def append(*args): return _openbabel.vVector3_append(*args)
    def empty(*args): return _openbabel.vVector3_empty(*args)
    def size(*args): return _openbabel.vVector3_size(*args)
    def clear(*args): return _openbabel.vVector3_clear(*args)
    def swap(*args): return _openbabel.vVector3_swap(*args)
    def get_allocator(*args): return _openbabel.vVector3_get_allocator(*args)
    def begin(*args): return _openbabel.vVector3_begin(*args)
    def end(*args): return _openbabel.vVector3_end(*args)
    def rbegin(*args): return _openbabel.vVector3_rbegin(*args)
    def rend(*args): return _openbabel.vVector3_rend(*args)
    def pop_back(*args): return _openbabel.vVector3_pop_back(*args)
    def erase(*args): return _openbabel.vVector3_erase(*args)
    def __init__(self, *args): 
        _openbabel.vVector3_swiginit(self,_openbabel.new_vVector3(*args))
    def push_back(*args): return _openbabel.vVector3_push_back(*args)
    def front(*args): return _openbabel.vVector3_front(*args)
    def back(*args): return _openbabel.vVector3_back(*args)
    def assign(*args): return _openbabel.vVector3_assign(*args)
    def resize(*args): return _openbabel.vVector3_resize(*args)
    def insert(*args): return _openbabel.vVector3_insert(*args)
    def reserve(*args): return _openbabel.vVector3_reserve(*args)
    def capacity(*args): return _openbabel.vVector3_capacity(*args)
    __swig_destroy__ = _openbabel.delete_vVector3
    __del__ = lambda self : None;
vVector3_swigregister = _openbabel.vVector3_swigregister
vVector3_swigregister(vVector3)

class vectorMol(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def iterator(*args): return _openbabel.vectorMol_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _openbabel.vectorMol___nonzero__(*args)
    def __len__(*args): return _openbabel.vectorMol___len__(*args)
    def pop(*args): return _openbabel.vectorMol_pop(*args)
    def __getslice__(*args): return _openbabel.vectorMol___getslice__(*args)
    def __setslice__(*args): return _openbabel.vectorMol___setslice__(*args)
    def __delslice__(*args): return _openbabel.vectorMol___delslice__(*args)
    def __delitem__(*args): return _openbabel.vectorMol___delitem__(*args)
    def __getitem__(*args): return _openbabel.vectorMol___getitem__(*args)
    def __setitem__(*args): return _openbabel.vectorMol___setitem__(*args)
    def append(*args): return _openbabel.vectorMol_append(*args)
    def empty(*args): return _openbabel.vectorMol_empty(*args)
    def size(*args): return _openbabel.vectorMol_size(*args)
    def clear(*args): return _openbabel.vectorMol_clear(*args)
    def swap(*args): return _openbabel.vectorMol_swap(*args)
    def get_allocator(*args): return _openbabel.vectorMol_get_allocator(*args)
    def begin(*args): return _openbabel.vectorMol_begin(*args)
    def end(*args): return _openbabel.vectorMol_end(*args)
    def rbegin(*args): return _openbabel.vectorMol_rbegin(*args)
    def rend(*args): return _openbabel.vectorMol_rend(*args)
    def pop_back(*args): return _openbabel.vectorMol_pop_back(*args)
    def erase(*args): return _openbabel.vectorMol_erase(*args)
    def __init__(self, *args): 
        _openbabel.vectorMol_swiginit(self,_openbabel.new_vectorMol(*args))
    def push_back(*args): return _openbabel.vectorMol_push_back(*args)
    def front(*args): return _openbabel.vectorMol_front(*args)
    def back(*args): return _openbabel.vectorMol_back(*args)
    def assign(*args): return _openbabel.vectorMol_assign(*args)
    def resize(*args): return _openbabel.vectorMol_resize(*args)
    def insert(*args): return _openbabel.vectorMol_insert(*args)
    def reserve(*args): return _openbabel.vectorMol_reserve(*args)
    def capacity(*args): return _openbabel.vectorMol_capacity(*args)
    __swig_destroy__ = _openbabel.delete_vectorMol
    __del__ = lambda self : None;
vectorMol_swigregister = _openbabel.vectorMol_swigregister
vectorMol_swigregister(vectorMol)

class vectorBond(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def iterator(*args): return _openbabel.vectorBond_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _openbabel.vectorBond___nonzero__(*args)
    def __len__(*args): return _openbabel.vectorBond___len__(*args)
    def pop(*args): return _openbabel.vectorBond_pop(*args)
    def __getslice__(*args): return _openbabel.vectorBond___getslice__(*args)
    def __setslice__(*args): return _openbabel.vectorBond___setslice__(*args)
    def __delslice__(*args): return _openbabel.vectorBond___delslice__(*args)
    def __delitem__(*args): return _openbabel.vectorBond___delitem__(*args)
    def __getitem__(*args): return _openbabel.vectorBond___getitem__(*args)
    def __setitem__(*args): return _openbabel.vectorBond___setitem__(*args)
    def append(*args): return _openbabel.vectorBond_append(*args)
    def empty(*args): return _openbabel.vectorBond_empty(*args)
    def size(*args): return _openbabel.vectorBond_size(*args)
    def clear(*args): return _openbabel.vectorBond_clear(*args)
    def swap(*args): return _openbabel.vectorBond_swap(*args)
    def get_allocator(*args): return _openbabel.vectorBond_get_allocator(*args)
    def begin(*args): return _openbabel.vectorBond_begin(*args)
    def end(*args): return _openbabel.vectorBond_end(*args)
    def rbegin(*args): return _openbabel.vectorBond_rbegin(*args)
    def rend(*args): return _openbabel.vectorBond_rend(*args)
    def pop_back(*args): return _openbabel.vectorBond_pop_back(*args)
    def erase(*args): return _openbabel.vectorBond_erase(*args)
    def __init__(self, *args): 
        _openbabel.vectorBond_swiginit(self,_openbabel.new_vectorBond(*args))
    def push_back(*args): return _openbabel.vectorBond_push_back(*args)
    def front(*args): return _openbabel.vectorBond_front(*args)
    def back(*args): return _openbabel.vectorBond_back(*args)
    def assign(*args): return _openbabel.vectorBond_assign(*args)
    def resize(*args): return _openbabel.vectorBond_resize(*args)
    def insert(*args): return _openbabel.vectorBond_insert(*args)
    def reserve(*args): return _openbabel.vectorBond_reserve(*args)
    def capacity(*args): return _openbabel.vectorBond_capacity(*args)
    __swig_destroy__ = _openbabel.delete_vectorBond
    __del__ = lambda self : None;
vectorBond_swigregister = _openbabel.vectorBond_swigregister
vectorBond_swigregister(vectorBond)

class vectorResidue(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def iterator(*args): return _openbabel.vectorResidue_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _openbabel.vectorResidue___nonzero__(*args)
    def __len__(*args): return _openbabel.vectorResidue___len__(*args)
    def pop(*args): return _openbabel.vectorResidue_pop(*args)
    def __getslice__(*args): return _openbabel.vectorResidue___getslice__(*args)
    def __setslice__(*args): return _openbabel.vectorResidue___setslice__(*args)
    def __delslice__(*args): return _openbabel.vectorResidue___delslice__(*args)
    def __delitem__(*args): return _openbabel.vectorResidue___delitem__(*args)
    def __getitem__(*args): return _openbabel.vectorResidue___getitem__(*args)
    def __setitem__(*args): return _openbabel.vectorResidue___setitem__(*args)
    def append(*args): return _openbabel.vectorResidue_append(*args)
    def empty(*args): return _openbabel.vectorResidue_empty(*args)
    def size(*args): return _openbabel.vectorResidue_size(*args)
    def clear(*args): return _openbabel.vectorResidue_clear(*args)
    def swap(*args): return _openbabel.vectorResidue_swap(*args)
    def get_allocator(*args): return _openbabel.vectorResidue_get_allocator(*args)
    def begin(*args): return _openbabel.vectorResidue_begin(*args)
    def end(*args): return _openbabel.vectorResidue_end(*args)
    def rbegin(*args): return _openbabel.vectorResidue_rbegin(*args)
    def rend(*args): return _openbabel.vectorResidue_rend(*args)
    def pop_back(*args): return _openbabel.vectorResidue_pop_back(*args)
    def erase(*args): return _openbabel.vectorResidue_erase(*args)
    def __init__(self, *args): 
        _openbabel.vectorResidue_swiginit(self,_openbabel.new_vectorResidue(*args))
    def push_back(*args): return _openbabel.vectorResidue_push_back(*args)
    def front(*args): return _openbabel.vectorResidue_front(*args)
    def back(*args): return _openbabel.vectorResidue_back(*args)
    def assign(*args): return _openbabel.vectorResidue_assign(*args)
    def resize(*args): return _openbabel.vectorResidue_resize(*args)
    def insert(*args): return _openbabel.vectorResidue_insert(*args)
    def reserve(*args): return _openbabel.vectorResidue_reserve(*args)
    def capacity(*args): return _openbabel.vectorResidue_capacity(*args)
    __swig_destroy__ = _openbabel.delete_vectorResidue
    __del__ = lambda self : None;
vectorResidue_swigregister = _openbabel.vectorResidue_swigregister
vectorResidue_swigregister(vectorResidue)

class vectorRing(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def iterator(*args): return _openbabel.vectorRing_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _openbabel.vectorRing___nonzero__(*args)
    def __len__(*args): return _openbabel.vectorRing___len__(*args)
    def pop(*args): return _openbabel.vectorRing_pop(*args)
    def __getslice__(*args): return _openbabel.vectorRing___getslice__(*args)
    def __setslice__(*args): return _openbabel.vectorRing___setslice__(*args)
    def __delslice__(*args): return _openbabel.vectorRing___delslice__(*args)
    def __delitem__(*args): return _openbabel.vectorRing___delitem__(*args)
    def __getitem__(*args): return _openbabel.vectorRing___getitem__(*args)
    def __setitem__(*args): return _openbabel.vectorRing___setitem__(*args)
    def append(*args): return _openbabel.vectorRing_append(*args)
    def empty(*args): return _openbabel.vectorRing_empty(*args)
    def size(*args): return _openbabel.vectorRing_size(*args)
    def clear(*args): return _openbabel.vectorRing_clear(*args)
    def swap(*args): return _openbabel.vectorRing_swap(*args)
    def get_allocator(*args): return _openbabel.vectorRing_get_allocator(*args)
    def begin(*args): return _openbabel.vectorRing_begin(*args)
    def end(*args): return _openbabel.vectorRing_end(*args)
    def rbegin(*args): return _openbabel.vectorRing_rbegin(*args)
    def rend(*args): return _openbabel.vectorRing_rend(*args)
    def pop_back(*args): return _openbabel.vectorRing_pop_back(*args)
    def erase(*args): return _openbabel.vectorRing_erase(*args)
    def __init__(self, *args): 
        _openbabel.vectorRing_swiginit(self,_openbabel.new_vectorRing(*args))
    def push_back(*args): return _openbabel.vectorRing_push_back(*args)
    def front(*args): return _openbabel.vectorRing_front(*args)
    def back(*args): return _openbabel.vectorRing_back(*args)
    def assign(*args): return _openbabel.vectorRing_assign(*args)
    def resize(*args): return _openbabel.vectorRing_resize(*args)
    def insert(*args): return _openbabel.vectorRing_insert(*args)
    def reserve(*args): return _openbabel.vectorRing_reserve(*args)
    def capacity(*args): return _openbabel.vectorRing_capacity(*args)
    __swig_destroy__ = _openbabel.delete_vectorRing
    __del__ = lambda self : None;
vectorRing_swigregister = _openbabel.vectorRing_swigregister
vectorRing_swigregister(vectorRing)

class OBGlobalDataBase(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBGlobalDataBase_swiginit(self,_openbabel.new_OBGlobalDataBase(*args))
    __swig_destroy__ = _openbabel.delete_OBGlobalDataBase
    __del__ = lambda self : None;
    def Init(*args): return _openbabel.OBGlobalDataBase_Init(*args)
    def GetSize(*args): return _openbabel.OBGlobalDataBase_GetSize(*args)
    def SetReadDirectory(*args): return _openbabel.OBGlobalDataBase_SetReadDirectory(*args)
    def SetEnvironmentVariable(*args): return _openbabel.OBGlobalDataBase_SetEnvironmentVariable(*args)
    def ParseLine(*args): return _openbabel.OBGlobalDataBase_ParseLine(*args)
OBGlobalDataBase_swigregister = _openbabel.OBGlobalDataBase_swigregister
OBGlobalDataBase_swigregister(OBGlobalDataBase)

class OBElement(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBElement_swiginit(self,_openbabel.new_OBElement(*args))
    def GetAtomicNum(*args): return _openbabel.OBElement_GetAtomicNum(*args)
    def GetSymbol(*args): return _openbabel.OBElement_GetSymbol(*args)
    def GetCovalentRad(*args): return _openbabel.OBElement_GetCovalentRad(*args)
    def GetVdwRad(*args): return _openbabel.OBElement_GetVdwRad(*args)
    def GetMass(*args): return _openbabel.OBElement_GetMass(*args)
    def GetMaxBonds(*args): return _openbabel.OBElement_GetMaxBonds(*args)
    def GetElectroNeg(*args): return _openbabel.OBElement_GetElectroNeg(*args)
    def GetIonization(*args): return _openbabel.OBElement_GetIonization(*args)
    def GetElectronAffinity(*args): return _openbabel.OBElement_GetElectronAffinity(*args)
    def GetName(*args): return _openbabel.OBElement_GetName(*args)
    def GetRed(*args): return _openbabel.OBElement_GetRed(*args)
    def GetGreen(*args): return _openbabel.OBElement_GetGreen(*args)
    def GetBlue(*args): return _openbabel.OBElement_GetBlue(*args)
    __swig_destroy__ = _openbabel.delete_OBElement
    __del__ = lambda self : None;
OBElement_swigregister = _openbabel.OBElement_swigregister
OBElement_swigregister(OBElement)

class OBElementTable(OBGlobalDataBase):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBElementTable_swiginit(self,_openbabel.new_OBElementTable(*args))
    __swig_destroy__ = _openbabel.delete_OBElementTable
    __del__ = lambda self : None;
    def ParseLine(*args): return _openbabel.OBElementTable_ParseLine(*args)
    def GetNumberOfElements(*args): return _openbabel.OBElementTable_GetNumberOfElements(*args)
    def GetSize(*args): return _openbabel.OBElementTable_GetSize(*args)
    def GetAtomicNum(*args): return _openbabel.OBElementTable_GetAtomicNum(*args)
    def GetSymbol(*args): return _openbabel.OBElementTable_GetSymbol(*args)
    def GetVdwRad(*args): return _openbabel.OBElementTable_GetVdwRad(*args)
    def GetCovalentRad(*args): return _openbabel.OBElementTable_GetCovalentRad(*args)
    def GetMass(*args): return _openbabel.OBElementTable_GetMass(*args)
    def CorrectedBondRad(*args): return _openbabel.OBElementTable_CorrectedBondRad(*args)
    def CorrectedVdwRad(*args): return _openbabel.OBElementTable_CorrectedVdwRad(*args)
    def GetMaxBonds(*args): return _openbabel.OBElementTable_GetMaxBonds(*args)
    def GetElectroNeg(*args): return _openbabel.OBElementTable_GetElectroNeg(*args)
    def GetIonization(*args): return _openbabel.OBElementTable_GetIonization(*args)
    def GetElectronAffinity(*args): return _openbabel.OBElementTable_GetElectronAffinity(*args)
    def GetRGB(*args): return _openbabel.OBElementTable_GetRGB(*args)
    def GetName(*args): return _openbabel.OBElementTable_GetName(*args)
OBElementTable_swigregister = _openbabel.OBElementTable_swigregister
OBElementTable_swigregister(OBElementTable)

class OBIsotopeTable(OBGlobalDataBase):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBIsotopeTable_swiginit(self,_openbabel.new_OBIsotopeTable(*args))
    __swig_destroy__ = _openbabel.delete_OBIsotopeTable
    __del__ = lambda self : None;
    def GetSize(*args): return _openbabel.OBIsotopeTable_GetSize(*args)
    def ParseLine(*args): return _openbabel.OBIsotopeTable_ParseLine(*args)
    def GetExactMass(*args): return _openbabel.OBIsotopeTable_GetExactMass(*args)
OBIsotopeTable_swigregister = _openbabel.OBIsotopeTable_swigregister
OBIsotopeTable_swigregister(OBIsotopeTable)

class OBTypeTable(OBGlobalDataBase):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBTypeTable_swiginit(self,_openbabel.new_OBTypeTable(*args))
    __swig_destroy__ = _openbabel.delete_OBTypeTable
    __del__ = lambda self : None;
    def ParseLine(*args): return _openbabel.OBTypeTable_ParseLine(*args)
    def GetSize(*args): return _openbabel.OBTypeTable_GetSize(*args)
    def SetFromType(*args): return _openbabel.OBTypeTable_SetFromType(*args)
    def SetToType(*args): return _openbabel.OBTypeTable_SetToType(*args)
    def Translate(*args): return _openbabel.OBTypeTable_Translate(*args)
    def GetFromType(*args): return _openbabel.OBTypeTable_GetFromType(*args)
    def GetToType(*args): return _openbabel.OBTypeTable_GetToType(*args)
OBTypeTable_swigregister = _openbabel.OBTypeTable_swigregister
OBTypeTable_swigregister(OBTypeTable)

class OBResidueData(OBGlobalDataBase):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBResidueData_swiginit(self,_openbabel.new_OBResidueData(*args))
    def ParseLine(*args): return _openbabel.OBResidueData_ParseLine(*args)
    def GetSize(*args): return _openbabel.OBResidueData_GetSize(*args)
    def SetResName(*args): return _openbabel.OBResidueData_SetResName(*args)
    def LookupBO(*args): return _openbabel.OBResidueData_LookupBO(*args)
    def LookupType(*args): return _openbabel.OBResidueData_LookupType(*args)
    def AssignBonds(*args): return _openbabel.OBResidueData_AssignBonds(*args)
    __swig_destroy__ = _openbabel.delete_OBResidueData
    __del__ = lambda self : None;
OBResidueData_swigregister = _openbabel.OBResidueData_swigregister
OBResidueData_swigregister(OBResidueData)

FILE_SEP_CHAR = _openbabel.FILE_SEP_CHAR
class OBStopwatch(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def Start(*args): return _openbabel.OBStopwatch_Start(*args)
    def Lap(*args): return _openbabel.OBStopwatch_Lap(*args)
    def Elapsed(*args): return _openbabel.OBStopwatch_Elapsed(*args)
    def __init__(self, *args): 
        _openbabel.OBStopwatch_swiginit(self,_openbabel.new_OBStopwatch(*args))
    __swig_destroy__ = _openbabel.delete_OBStopwatch
    __del__ = lambda self : None;
OBStopwatch_swigregister = _openbabel.OBStopwatch_swigregister
OBStopwatch_swigregister(OBStopwatch)

class OBSqrtTbl(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBSqrtTbl_swiginit(self,_openbabel.new_OBSqrtTbl(*args))
    __swig_destroy__ = _openbabel.delete_OBSqrtTbl
    __del__ = lambda self : None;
    def Sqrt(*args): return _openbabel.OBSqrtTbl_Sqrt(*args)
    def Init(*args): return _openbabel.OBSqrtTbl_Init(*args)
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
    __del__ = lambda self : None;
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
    def Seed(*args): return _openbabel.OBRandom_Seed(*args)
    def TimeSeed(*args): return _openbabel.OBRandom_TimeSeed(*args)
    def NextInt(*args): return _openbabel.OBRandom_NextInt(*args)
    def NextFloat(*args): return _openbabel.OBRandom_NextFloat(*args)
    __swig_destroy__ = _openbabel.delete_OBRandom
    __del__ = lambda self : None;
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
    def Set(*args): return _openbabel.vector3_Set(*args)
    def SetX(*args): return _openbabel.vector3_SetX(*args)
    def SetY(*args): return _openbabel.vector3_SetY(*args)
    def SetZ(*args): return _openbabel.vector3_SetZ(*args)
    def Get(*args): return _openbabel.vector3_Get(*args)
    def __iadd__(*args): return _openbabel.vector3___iadd__(*args)
    def __isub__(*args): return _openbabel.vector3___isub__(*args)
    def __idiv__(*args): return _openbabel.vector3___idiv__(*args)
    def __imul__(*args): return _openbabel.vector3___imul__(*args)
    def randomUnitVector(*args): return _openbabel.vector3_randomUnitVector(*args)
    def normalize(*args): return _openbabel.vector3_normalize(*args)
    def length(*args): return _openbabel.vector3_length(*args)
    def length_2(*args): return _openbabel.vector3_length_2(*args)
    def x(*args): return _openbabel.vector3_x(*args)
    def y(*args): return _openbabel.vector3_y(*args)
    def z(*args): return _openbabel.vector3_z(*args)
    def distSq(*args): return _openbabel.vector3_distSq(*args)
    def createOrthoVector(*args): return _openbabel.vector3_createOrthoVector(*args)
    __swig_destroy__ = _openbabel.delete_vector3
    __del__ = lambda self : None;
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
    def ReadMolecule(*args): return _openbabel.OBFormat_ReadMolecule(*args)
    def ReadChemObject(*args): return _openbabel.OBFormat_ReadChemObject(*args)
    def WriteMolecule(*args): return _openbabel.OBFormat_WriteMolecule(*args)
    def WriteChemObject(*args): return _openbabel.OBFormat_WriteChemObject(*args)
    def Description(*args): return _openbabel.OBFormat_Description(*args)
    def TargetClassDescription(*args): return _openbabel.OBFormat_TargetClassDescription(*args)
    def GetType(*args): return _openbabel.OBFormat_GetType(*args)
    def SpecificationURL(*args): return _openbabel.OBFormat_SpecificationURL(*args)
    def GetMIMEType(*args): return _openbabel.OBFormat_GetMIMEType(*args)
    def Flags(*args): return _openbabel.OBFormat_Flags(*args)
    def SkipObjects(*args): return _openbabel.OBFormat_SkipObjects(*args)
    def MakeNewInstance(*args): return _openbabel.OBFormat_MakeNewInstance(*args)
    __swig_destroy__ = _openbabel.delete_OBFormat
    __del__ = lambda self : None;
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
    def __call__(*args): return _openbabel.CharPtrLess___call__(*args)
    def __init__(self, *args): 
        _openbabel.CharPtrLess_swiginit(self,_openbabel.new_CharPtrLess(*args))
    __swig_destroy__ = _openbabel.delete_CharPtrLess
    __del__ = lambda self : None;
CharPtrLess_swigregister = _openbabel.CharPtrLess_swigregister
CharPtrLess_swigregister(CharPtrLess)

class OBConversion(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBConversion_swiginit(self,_openbabel.new_OBConversion(*args))
    __swig_destroy__ = _openbabel.delete_OBConversion
    __del__ = lambda self : None;
    RegisterFormat = staticmethod(_openbabel.OBConversion_RegisterFormat)
    FindFormat = staticmethod(_openbabel.OBConversion_FindFormat)
    FormatFromExt = staticmethod(_openbabel.OBConversion_FormatFromExt)
    FormatFromMIME = staticmethod(_openbabel.OBConversion_FormatFromMIME)
    GetNextFormat = staticmethod(_openbabel.OBConversion_GetNextFormat)
    Description = staticmethod(_openbabel.OBConversion_Description)
    def GetInStream(*args): return _openbabel.OBConversion_GetInStream(*args)
    def GetOutStream(*args): return _openbabel.OBConversion_GetOutStream(*args)
    def SetInStream(*args): return _openbabel.OBConversion_SetInStream(*args)
    def SetOutStream(*args): return _openbabel.OBConversion_SetOutStream(*args)
    def SetInAndOutFormats(*args): return _openbabel.OBConversion_SetInAndOutFormats(*args)
    def SetInFormat(*args): return _openbabel.OBConversion_SetInFormat(*args)
    def SetOutFormat(*args): return _openbabel.OBConversion_SetOutFormat(*args)
    def GetInFormat(*args): return _openbabel.OBConversion_GetInFormat(*args)
    def GetOutFormat(*args): return _openbabel.OBConversion_GetOutFormat(*args)
    def GetInFilename(*args): return _openbabel.OBConversion_GetInFilename(*args)
    def GetInPos(*args): return _openbabel.OBConversion_GetInPos(*args)
    def GetInLen(*args): return _openbabel.OBConversion_GetInLen(*args)
    def GetTitle(*args): return _openbabel.OBConversion_GetTitle(*args)
    def GetAuxConv(*args): return _openbabel.OBConversion_GetAuxConv(*args)
    def SetAuxConv(*args): return _openbabel.OBConversion_SetAuxConv(*args)
    INOPTIONS = _openbabel.OBConversion_INOPTIONS
    OUTOPTIONS = _openbabel.OBConversion_OUTOPTIONS
    GENOPTIONS = _openbabel.OBConversion_GENOPTIONS
    def IsOption(*args): return _openbabel.OBConversion_IsOption(*args)
    def GetOptions(*args): return _openbabel.OBConversion_GetOptions(*args)
    def AddOption(*args): return _openbabel.OBConversion_AddOption(*args)
    def RemoveOption(*args): return _openbabel.OBConversion_RemoveOption(*args)
    def SetOptions(*args): return _openbabel.OBConversion_SetOptions(*args)
    RegisterOptionParam = staticmethod(_openbabel.OBConversion_RegisterOptionParam)
    GetOptionParams = staticmethod(_openbabel.OBConversion_GetOptionParams)
    def Convert(*args): return _openbabel.OBConversion_Convert(*args)
    def FullConvert(*args): return _openbabel.OBConversion_FullConvert(*args)
    def AddChemObject(*args): return _openbabel.OBConversion_AddChemObject(*args)
    def GetChemObject(*args): return _openbabel.OBConversion_GetChemObject(*args)
    def IsLast(*args): return _openbabel.OBConversion_IsLast(*args)
    def IsFirstInput(*args): return _openbabel.OBConversion_IsFirstInput(*args)
    def GetOutputIndex(*args): return _openbabel.OBConversion_GetOutputIndex(*args)
    def SetOutputIndex(*args): return _openbabel.OBConversion_SetOutputIndex(*args)
    def SetMoreFilesToCome(*args): return _openbabel.OBConversion_SetMoreFilesToCome(*args)
    def SetOneObjectOnly(*args): return _openbabel.OBConversion_SetOneObjectOnly(*args)
    GetDefaultFormat = staticmethod(_openbabel.OBConversion_GetDefaultFormat)
    def Write(*args): return _openbabel.OBConversion_Write(*args)
    def WriteString(*args): return _openbabel.OBConversion_WriteString(*args)
    def WriteFile(*args): return _openbabel.OBConversion_WriteFile(*args)
    def Read(*args): return _openbabel.OBConversion_Read(*args)
    def ReadString(*args): return _openbabel.OBConversion_ReadString(*args)
    def ReadFile(*args): return _openbabel.OBConversion_ReadFile(*args)
    BatchFileName = staticmethod(_openbabel.OBConversion_BatchFileName)
    IncrementedFileName = staticmethod(_openbabel.OBConversion_IncrementedFileName)
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
    __del__ = lambda self : None;
    def AddAtom(*args): return _openbabel.OBResidue_AddAtom(*args)
    def InsertAtom(*args): return _openbabel.OBResidue_InsertAtom(*args)
    def RemoveAtom(*args): return _openbabel.OBResidue_RemoveAtom(*args)
    def Clear(*args): return _openbabel.OBResidue_Clear(*args)
    def SetName(*args): return _openbabel.OBResidue_SetName(*args)
    def SetNum(*args): return _openbabel.OBResidue_SetNum(*args)
    def SetChain(*args): return _openbabel.OBResidue_SetChain(*args)
    def SetChainNum(*args): return _openbabel.OBResidue_SetChainNum(*args)
    def SetIdx(*args): return _openbabel.OBResidue_SetIdx(*args)
    def SetAtomID(*args): return _openbabel.OBResidue_SetAtomID(*args)
    def SetHetAtom(*args): return _openbabel.OBResidue_SetHetAtom(*args)
    def SetSerialNum(*args): return _openbabel.OBResidue_SetSerialNum(*args)
    def GetName(*args): return _openbabel.OBResidue_GetName(*args)
    def GetNum(*args): return _openbabel.OBResidue_GetNum(*args)
    def GetNumAtoms(*args): return _openbabel.OBResidue_GetNumAtoms(*args)
    def GetChain(*args): return _openbabel.OBResidue_GetChain(*args)
    def GetChainNum(*args): return _openbabel.OBResidue_GetChainNum(*args)
    def GetIdx(*args): return _openbabel.OBResidue_GetIdx(*args)
    def GetResKey(*args): return _openbabel.OBResidue_GetResKey(*args)
    def GetAtoms(*args): return _openbabel.OBResidue_GetAtoms(*args)
    def GetBonds(*args): return _openbabel.OBResidue_GetBonds(*args)
    def GetAtomID(*args): return _openbabel.OBResidue_GetAtomID(*args)
    def GetSerialNum(*args): return _openbabel.OBResidue_GetSerialNum(*args)
    def GetAminoAcidProperty(*args): return _openbabel.OBResidue_GetAminoAcidProperty(*args)
    def GetAtomProperty(*args): return _openbabel.OBResidue_GetAtomProperty(*args)
    def GetResidueProperty(*args): return _openbabel.OBResidue_GetResidueProperty(*args)
    def IsHetAtom(*args): return _openbabel.OBResidue_IsHetAtom(*args)
    def IsResidueType(*args): return _openbabel.OBResidue_IsResidueType(*args)
    def BeginAtom(*args): return _openbabel.OBResidue_BeginAtom(*args)
    def NextAtom(*args): return _openbabel.OBResidue_NextAtom(*args)
    def HasData(*args): return _openbabel.OBResidue_HasData(*args)
    def DeleteData(*args): return _openbabel.OBResidue_DeleteData(*args)
    def SetData(*args): return _openbabel.OBResidue_SetData(*args)
    def DataSize(*args): return _openbabel.OBResidue_DataSize(*args)
    def GetData(*args): return _openbabel.OBResidue_GetData(*args)
    def BeginData(*args): return _openbabel.OBResidue_BeginData(*args)
    def EndData(*args): return _openbabel.OBResidue_EndData(*args)
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
    __del__ = lambda self : None;
    def Clear(*args): return _openbabel.OBAtom_Clear(*args)
    def SetIdx(*args): return _openbabel.OBAtom_SetIdx(*args)
    def SetHyb(*args): return _openbabel.OBAtom_SetHyb(*args)
    def SetAtomicNum(*args): return _openbabel.OBAtom_SetAtomicNum(*args)
    def SetIsotope(*args): return _openbabel.OBAtom_SetIsotope(*args)
    def SetImplicitValence(*args): return _openbabel.OBAtom_SetImplicitValence(*args)
    def IncrementImplicitValence(*args): return _openbabel.OBAtom_IncrementImplicitValence(*args)
    def DecrementImplicitValence(*args): return _openbabel.OBAtom_DecrementImplicitValence(*args)
    def SetFormalCharge(*args): return _openbabel.OBAtom_SetFormalCharge(*args)
    def SetSpinMultiplicity(*args): return _openbabel.OBAtom_SetSpinMultiplicity(*args)
    def SetType(*args): return _openbabel.OBAtom_SetType(*args)
    def SetPartialCharge(*args): return _openbabel.OBAtom_SetPartialCharge(*args)
    def SetCoordPtr(*args): return _openbabel.OBAtom_SetCoordPtr(*args)
    def SetVector(*args): return _openbabel.OBAtom_SetVector(*args)
    def SetResidue(*args): return _openbabel.OBAtom_SetResidue(*args)
    def SetAromatic(*args): return _openbabel.OBAtom_SetAromatic(*args)
    def UnsetAromatic(*args): return _openbabel.OBAtom_UnsetAromatic(*args)
    def SetClockwiseStereo(*args): return _openbabel.OBAtom_SetClockwiseStereo(*args)
    def SetAntiClockwiseStereo(*args): return _openbabel.OBAtom_SetAntiClockwiseStereo(*args)
    def SetPositiveStereo(*args): return _openbabel.OBAtom_SetPositiveStereo(*args)
    def SetNegativeStereo(*args): return _openbabel.OBAtom_SetNegativeStereo(*args)
    def UnsetStereo(*args): return _openbabel.OBAtom_UnsetStereo(*args)
    def SetInRing(*args): return _openbabel.OBAtom_SetInRing(*args)
    def SetChiral(*args): return _openbabel.OBAtom_SetChiral(*args)
    def ClearCoordPtr(*args): return _openbabel.OBAtom_ClearCoordPtr(*args)
    def GetFormalCharge(*args): return _openbabel.OBAtom_GetFormalCharge(*args)
    def GetAtomicNum(*args): return _openbabel.OBAtom_GetAtomicNum(*args)
    def GetIsotope(*args): return _openbabel.OBAtom_GetIsotope(*args)
    def GetSpinMultiplicity(*args): return _openbabel.OBAtom_GetSpinMultiplicity(*args)
    def GetAtomicMass(*args): return _openbabel.OBAtom_GetAtomicMass(*args)
    def GetExactMass(*args): return _openbabel.OBAtom_GetExactMass(*args)
    def GetIdx(*args): return _openbabel.OBAtom_GetIdx(*args)
    def GetCoordinateIdx(*args): return _openbabel.OBAtom_GetCoordinateIdx(*args)
    def GetCIdx(*args): return _openbabel.OBAtom_GetCIdx(*args)
    def GetValence(*args): return _openbabel.OBAtom_GetValence(*args)
    def GetHyb(*args): return _openbabel.OBAtom_GetHyb(*args)
    def GetImplicitValence(*args): return _openbabel.OBAtom_GetImplicitValence(*args)
    def GetHvyValence(*args): return _openbabel.OBAtom_GetHvyValence(*args)
    def GetHeteroValence(*args): return _openbabel.OBAtom_GetHeteroValence(*args)
    def GetType(*args): return _openbabel.OBAtom_GetType(*args)
    def GetX(*args): return _openbabel.OBAtom_GetX(*args)
    def GetY(*args): return _openbabel.OBAtom_GetY(*args)
    def GetZ(*args): return _openbabel.OBAtom_GetZ(*args)
    def x(*args): return _openbabel.OBAtom_x(*args)
    def y(*args): return _openbabel.OBAtom_y(*args)
    def z(*args): return _openbabel.OBAtom_z(*args)
    def GetCoordinate(*args): return _openbabel.OBAtom_GetCoordinate(*args)
    def GetVector(*args): return _openbabel.OBAtom_GetVector(*args)
    def GetPartialCharge(*args): return _openbabel.OBAtom_GetPartialCharge(*args)
    def GetResidue(*args): return _openbabel.OBAtom_GetResidue(*args)
    def GetNewBondVector(*args): return _openbabel.OBAtom_GetNewBondVector(*args)
    def GetBond(*args): return _openbabel.OBAtom_GetBond(*args)
    def GetNextAtom(*args): return _openbabel.OBAtom_GetNextAtom(*args)
    def BeginBonds(*args): return _openbabel.OBAtom_BeginBonds(*args)
    def EndBonds(*args): return _openbabel.OBAtom_EndBonds(*args)
    def BeginBond(*args): return _openbabel.OBAtom_BeginBond(*args)
    def NextBond(*args): return _openbabel.OBAtom_NextBond(*args)
    def BeginNbrAtom(*args): return _openbabel.OBAtom_BeginNbrAtom(*args)
    def NextNbrAtom(*args): return _openbabel.OBAtom_NextNbrAtom(*args)
    def GetDistance(*args): return _openbabel.OBAtom_GetDistance(*args)
    def GetAngle(*args): return _openbabel.OBAtom_GetAngle(*args)
    def NewResidue(*args): return _openbabel.OBAtom_NewResidue(*args)
    def DeleteResidue(*args): return _openbabel.OBAtom_DeleteResidue(*args)
    def AddBond(*args): return _openbabel.OBAtom_AddBond(*args)
    def InsertBond(*args): return _openbabel.OBAtom_InsertBond(*args)
    def DeleteBond(*args): return _openbabel.OBAtom_DeleteBond(*args)
    def ClearBond(*args): return _openbabel.OBAtom_ClearBond(*args)
    def CountFreeOxygens(*args): return _openbabel.OBAtom_CountFreeOxygens(*args)
    def ImplicitHydrogenCount(*args): return _openbabel.OBAtom_ImplicitHydrogenCount(*args)
    def ExplicitHydrogenCount(*args): return _openbabel.OBAtom_ExplicitHydrogenCount(*args)
    def MemberOfRingCount(*args): return _openbabel.OBAtom_MemberOfRingCount(*args)
    def MemberOfRingSize(*args): return _openbabel.OBAtom_MemberOfRingSize(*args)
    def SmallestBondAngle(*args): return _openbabel.OBAtom_SmallestBondAngle(*args)
    def AverageBondAngle(*args): return _openbabel.OBAtom_AverageBondAngle(*args)
    def BOSum(*args): return _openbabel.OBAtom_BOSum(*args)
    def KBOSum(*args): return _openbabel.OBAtom_KBOSum(*args)
    def HtoMethyl(*args): return _openbabel.OBAtom_HtoMethyl(*args)
    def SetHybAndGeom(*args): return _openbabel.OBAtom_SetHybAndGeom(*args)
    def HasResidue(*args): return _openbabel.OBAtom_HasResidue(*args)
    def IsHydrogen(*args): return _openbabel.OBAtom_IsHydrogen(*args)
    def IsCarbon(*args): return _openbabel.OBAtom_IsCarbon(*args)
    def IsNitrogen(*args): return _openbabel.OBAtom_IsNitrogen(*args)
    def IsOxygen(*args): return _openbabel.OBAtom_IsOxygen(*args)
    def IsSulfur(*args): return _openbabel.OBAtom_IsSulfur(*args)
    def IsPhosphorus(*args): return _openbabel.OBAtom_IsPhosphorus(*args)
    def IsAromatic(*args): return _openbabel.OBAtom_IsAromatic(*args)
    def IsInRing(*args): return _openbabel.OBAtom_IsInRing(*args)
    def IsInRingSize(*args): return _openbabel.OBAtom_IsInRingSize(*args)
    def IsHeteroatom(*args): return _openbabel.OBAtom_IsHeteroatom(*args)
    def IsNotCorH(*args): return _openbabel.OBAtom_IsNotCorH(*args)
    def IsConnected(*args): return _openbabel.OBAtom_IsConnected(*args)
    def IsOneThree(*args): return _openbabel.OBAtom_IsOneThree(*args)
    def IsOneFour(*args): return _openbabel.OBAtom_IsOneFour(*args)
    def IsCarboxylOxygen(*args): return _openbabel.OBAtom_IsCarboxylOxygen(*args)
    def IsPhosphateOxygen(*args): return _openbabel.OBAtom_IsPhosphateOxygen(*args)
    def IsSulfateOxygen(*args): return _openbabel.OBAtom_IsSulfateOxygen(*args)
    def IsNitroOxygen(*args): return _openbabel.OBAtom_IsNitroOxygen(*args)
    def IsAmideNitrogen(*args): return _openbabel.OBAtom_IsAmideNitrogen(*args)
    def IsPolarHydrogen(*args): return _openbabel.OBAtom_IsPolarHydrogen(*args)
    def IsNonPolarHydrogen(*args): return _openbabel.OBAtom_IsNonPolarHydrogen(*args)
    def IsAromaticNOxide(*args): return _openbabel.OBAtom_IsAromaticNOxide(*args)
    def IsChiral(*args): return _openbabel.OBAtom_IsChiral(*args)
    def IsAxial(*args): return _openbabel.OBAtom_IsAxial(*args)
    def IsClockwise(*args): return _openbabel.OBAtom_IsClockwise(*args)
    def IsAntiClockwise(*args): return _openbabel.OBAtom_IsAntiClockwise(*args)
    def IsPositiveStereo(*args): return _openbabel.OBAtom_IsPositiveStereo(*args)
    def IsNegativeStereo(*args): return _openbabel.OBAtom_IsNegativeStereo(*args)
    def HasChiralitySpecified(*args): return _openbabel.OBAtom_HasChiralitySpecified(*args)
    def HasChiralVolume(*args): return _openbabel.OBAtom_HasChiralVolume(*args)
    def IsHbondAcceptor(*args): return _openbabel.OBAtom_IsHbondAcceptor(*args)
    def IsHbondDonor(*args): return _openbabel.OBAtom_IsHbondDonor(*args)
    def IsHbondDonorH(*args): return _openbabel.OBAtom_IsHbondDonorH(*args)
    def HasAlphaBetaUnsat(*args): return _openbabel.OBAtom_HasAlphaBetaUnsat(*args)
    def HasBondOfOrder(*args): return _openbabel.OBAtom_HasBondOfOrder(*args)
    def CountBondsOfOrder(*args): return _openbabel.OBAtom_CountBondsOfOrder(*args)
    def HasNonSingleBond(*args): return _openbabel.OBAtom_HasNonSingleBond(*args)
    def HasSingleBond(*args): return _openbabel.OBAtom_HasSingleBond(*args)
    def HasDoubleBond(*args): return _openbabel.OBAtom_HasDoubleBond(*args)
    def HasAromaticBond(*args): return _openbabel.OBAtom_HasAromaticBond(*args)
    def MatchesSMARTS(*args): return _openbabel.OBAtom_MatchesSMARTS(*args)
    def HasData(*args): return _openbabel.OBAtom_HasData(*args)
    def DeleteData(*args): return _openbabel.OBAtom_DeleteData(*args)
    def SetData(*args): return _openbabel.OBAtom_SetData(*args)
    def DataSize(*args): return _openbabel.OBAtom_DataSize(*args)
    def GetData(*args): return _openbabel.OBAtom_GetData(*args)
    def BeginData(*args): return _openbabel.OBAtom_BeginData(*args)
    def EndData(*args): return _openbabel.OBAtom_EndData(*args)
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
    __del__ = lambda self : None;
    def SetIdx(*args): return _openbabel.OBBond_SetIdx(*args)
    def SetBO(*args): return _openbabel.OBBond_SetBO(*args)
    def SetBegin(*args): return _openbabel.OBBond_SetBegin(*args)
    def SetEnd(*args): return _openbabel.OBBond_SetEnd(*args)
    def SetLength(*args): return _openbabel.OBBond_SetLength(*args)
    def Set(*args): return _openbabel.OBBond_Set(*args)
    def SetKSingle(*args): return _openbabel.OBBond_SetKSingle(*args)
    def SetKDouble(*args): return _openbabel.OBBond_SetKDouble(*args)
    def SetKTriple(*args): return _openbabel.OBBond_SetKTriple(*args)
    def SetAromatic(*args): return _openbabel.OBBond_SetAromatic(*args)
    def SetHash(*args): return _openbabel.OBBond_SetHash(*args)
    def SetWedge(*args): return _openbabel.OBBond_SetWedge(*args)
    def SetUp(*args): return _openbabel.OBBond_SetUp(*args)
    def SetDown(*args): return _openbabel.OBBond_SetDown(*args)
    def SetInRing(*args): return _openbabel.OBBond_SetInRing(*args)
    def SetClosure(*args): return _openbabel.OBBond_SetClosure(*args)
    def UnsetAromatic(*args): return _openbabel.OBBond_UnsetAromatic(*args)
    def UnsetKekule(*args): return _openbabel.OBBond_UnsetKekule(*args)
    def GetBO(*args): return _openbabel.OBBond_GetBO(*args)
    def GetBondOrder(*args): return _openbabel.OBBond_GetBondOrder(*args)
    def GetFlags(*args): return _openbabel.OBBond_GetFlags(*args)
    def GetBeginAtomIdx(*args): return _openbabel.OBBond_GetBeginAtomIdx(*args)
    def GetEndAtomIdx(*args): return _openbabel.OBBond_GetEndAtomIdx(*args)
    def GetBeginAtom(*args): return _openbabel.OBBond_GetBeginAtom(*args)
    def GetEndAtom(*args): return _openbabel.OBBond_GetEndAtom(*args)
    def GetNbrAtom(*args): return _openbabel.OBBond_GetNbrAtom(*args)
    def GetEquibLength(*args): return _openbabel.OBBond_GetEquibLength(*args)
    def GetLength(*args): return _openbabel.OBBond_GetLength(*args)
    def GetNbrAtomIdx(*args): return _openbabel.OBBond_GetNbrAtomIdx(*args)
    def IsAromatic(*args): return _openbabel.OBBond_IsAromatic(*args)
    def IsInRing(*args): return _openbabel.OBBond_IsInRing(*args)
    def IsRotor(*args): return _openbabel.OBBond_IsRotor(*args)
    def IsAmide(*args): return _openbabel.OBBond_IsAmide(*args)
    def IsPrimaryAmide(*args): return _openbabel.OBBond_IsPrimaryAmide(*args)
    def IsSecondaryAmide(*args): return _openbabel.OBBond_IsSecondaryAmide(*args)
    def IsEster(*args): return _openbabel.OBBond_IsEster(*args)
    def IsCarbonyl(*args): return _openbabel.OBBond_IsCarbonyl(*args)
    def IsSingle(*args): return _openbabel.OBBond_IsSingle(*args)
    def IsDouble(*args): return _openbabel.OBBond_IsDouble(*args)
    def IsTriple(*args): return _openbabel.OBBond_IsTriple(*args)
    def IsKSingle(*args): return _openbabel.OBBond_IsKSingle(*args)
    def IsKDouble(*args): return _openbabel.OBBond_IsKDouble(*args)
    def IsKTriple(*args): return _openbabel.OBBond_IsKTriple(*args)
    def IsClosure(*args): return _openbabel.OBBond_IsClosure(*args)
    def IsUp(*args): return _openbabel.OBBond_IsUp(*args)
    def IsDown(*args): return _openbabel.OBBond_IsDown(*args)
    def IsWedge(*args): return _openbabel.OBBond_IsWedge(*args)
    def IsHash(*args): return _openbabel.OBBond_IsHash(*args)
    def IsDoubleBondGeometry(*args): return _openbabel.OBBond_IsDoubleBondGeometry(*args)
    def HasData(*args): return _openbabel.OBBond_HasData(*args)
    def DeleteData(*args): return _openbabel.OBBond_DeleteData(*args)
    def SetData(*args): return _openbabel.OBBond_SetData(*args)
    def DataSize(*args): return _openbabel.OBBond_DataSize(*args)
    def GetData(*args): return _openbabel.OBBond_GetData(*args)
    def BeginData(*args): return _openbabel.OBBond_BeginData(*args)
    def EndData(*args): return _openbabel.OBBond_EndData(*args)
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
    __del__ = lambda self : None;
    def __iadd__(*args): return _openbabel.OBMol___iadd__(*args)
    def ReserveAtoms(*args): return _openbabel.OBMol_ReserveAtoms(*args)
    def CreateAtom(*args): return _openbabel.OBMol_CreateAtom(*args)
    def CreateBond(*args): return _openbabel.OBMol_CreateBond(*args)
    def DestroyAtom(*args): return _openbabel.OBMol_DestroyAtom(*args)
    def DestroyBond(*args): return _openbabel.OBMol_DestroyBond(*args)
    def AddAtom(*args): return _openbabel.OBMol_AddAtom(*args)
    def AddBond(*args): return _openbabel.OBMol_AddBond(*args)
    def AddResidue(*args): return _openbabel.OBMol_AddResidue(*args)
    def InsertAtom(*args): return _openbabel.OBMol_InsertAtom(*args)
    def DeleteAtom(*args): return _openbabel.OBMol_DeleteAtom(*args)
    def DeleteBond(*args): return _openbabel.OBMol_DeleteBond(*args)
    def DeleteResidue(*args): return _openbabel.OBMol_DeleteResidue(*args)
    def NewAtom(*args): return _openbabel.OBMol_NewAtom(*args)
    def NewResidue(*args): return _openbabel.OBMol_NewResidue(*args)
    def BeginModify(*args): return _openbabel.OBMol_BeginModify(*args)
    def EndModify(*args): return _openbabel.OBMol_EndModify(*args)
    def GetMod(*args): return _openbabel.OBMol_GetMod(*args)
    def IncrementMod(*args): return _openbabel.OBMol_IncrementMod(*args)
    def DecrementMod(*args): return _openbabel.OBMol_DecrementMod(*args)
    def HasData(*args): return _openbabel.OBMol_HasData(*args)
    def DeleteData(*args): return _openbabel.OBMol_DeleteData(*args)
    def SetData(*args): return _openbabel.OBMol_SetData(*args)
    def DataSize(*args): return _openbabel.OBMol_DataSize(*args)
    def GetData(*args): return _openbabel.OBMol_GetData(*args)
    def BeginData(*args): return _openbabel.OBMol_BeginData(*args)
    def EndData(*args): return _openbabel.OBMol_EndData(*args)
    def GetFlags(*args): return _openbabel.OBMol_GetFlags(*args)
    def GetTitle(*args): return _openbabel.OBMol_GetTitle(*args)
    def NumAtoms(*args): return _openbabel.OBMol_NumAtoms(*args)
    def NumBonds(*args): return _openbabel.OBMol_NumBonds(*args)
    def NumHvyAtoms(*args): return _openbabel.OBMol_NumHvyAtoms(*args)
    def NumResidues(*args): return _openbabel.OBMol_NumResidues(*args)
    def NumRotors(*args): return _openbabel.OBMol_NumRotors(*args)
    def GetAtom(*args): return _openbabel.OBMol_GetAtom(*args)
    def GetFirstAtom(*args): return _openbabel.OBMol_GetFirstAtom(*args)
    def GetBond(*args): return _openbabel.OBMol_GetBond(*args)
    def GetResidue(*args): return _openbabel.OBMol_GetResidue(*args)
    def GetInternalCoord(*args): return _openbabel.OBMol_GetInternalCoord(*args)
    def GetTorsion(*args): return _openbabel.OBMol_GetTorsion(*args)
    def GetFormula(*args): return _openbabel.OBMol_GetFormula(*args)
    def GetEnergy(*args): return _openbabel.OBMol_GetEnergy(*args)
    def GetMolWt(*args): return _openbabel.OBMol_GetMolWt(*args)
    def GetExactMass(*args): return _openbabel.OBMol_GetExactMass(*args)
    def GetTotalCharge(*args): return _openbabel.OBMol_GetTotalCharge(*args)
    def GetTotalSpinMultiplicity(*args): return _openbabel.OBMol_GetTotalSpinMultiplicity(*args)
    def GetDimension(*args): return _openbabel.OBMol_GetDimension(*args)
    def GetCoordinates(*args): return _openbabel.OBMol_GetCoordinates(*args)
    def GetSSSR(*args): return _openbabel.OBMol_GetSSSR(*args)
    def AutomaticFormalCharge(*args): return _openbabel.OBMol_AutomaticFormalCharge(*args)
    def AutomaticPartialCharge(*args): return _openbabel.OBMol_AutomaticPartialCharge(*args)
    def SetTitle(*args): return _openbabel.OBMol_SetTitle(*args)
    def SetFormula(*args): return _openbabel.OBMol_SetFormula(*args)
    def SetEnergy(*args): return _openbabel.OBMol_SetEnergy(*args)
    def SetDimension(*args): return _openbabel.OBMol_SetDimension(*args)
    def SetTotalCharge(*args): return _openbabel.OBMol_SetTotalCharge(*args)
    def SetTotalSpinMultiplicity(*args): return _openbabel.OBMol_SetTotalSpinMultiplicity(*args)
    def SetInternalCoord(*args): return _openbabel.OBMol_SetInternalCoord(*args)
    def SetAutomaticFormalCharge(*args): return _openbabel.OBMol_SetAutomaticFormalCharge(*args)
    def SetAutomaticPartialCharge(*args): return _openbabel.OBMol_SetAutomaticPartialCharge(*args)
    def SetAromaticPerceived(*args): return _openbabel.OBMol_SetAromaticPerceived(*args)
    def SetSSSRPerceived(*args): return _openbabel.OBMol_SetSSSRPerceived(*args)
    def SetRingAtomsAndBondsPerceived(*args): return _openbabel.OBMol_SetRingAtomsAndBondsPerceived(*args)
    def SetAtomTypesPerceived(*args): return _openbabel.OBMol_SetAtomTypesPerceived(*args)
    def SetChainsPerceived(*args): return _openbabel.OBMol_SetChainsPerceived(*args)
    def SetChiralityPerceived(*args): return _openbabel.OBMol_SetChiralityPerceived(*args)
    def SetPartialChargesPerceived(*args): return _openbabel.OBMol_SetPartialChargesPerceived(*args)
    def SetHybridizationPerceived(*args): return _openbabel.OBMol_SetHybridizationPerceived(*args)
    def SetImplicitValencePerceived(*args): return _openbabel.OBMol_SetImplicitValencePerceived(*args)
    def SetKekulePerceived(*args): return _openbabel.OBMol_SetKekulePerceived(*args)
    def SetClosureBondsPerceived(*args): return _openbabel.OBMol_SetClosureBondsPerceived(*args)
    def SetHydrogensAdded(*args): return _openbabel.OBMol_SetHydrogensAdded(*args)
    def SetCorrectedForPH(*args): return _openbabel.OBMol_SetCorrectedForPH(*args)
    def SetAromaticCorrected(*args): return _openbabel.OBMol_SetAromaticCorrected(*args)
    def SetSpinMultiplicityAssigned(*args): return _openbabel.OBMol_SetSpinMultiplicityAssigned(*args)
    def SetFlags(*args): return _openbabel.OBMol_SetFlags(*args)
    def UnsetAromaticPerceived(*args): return _openbabel.OBMol_UnsetAromaticPerceived(*args)
    def UnsetPartialChargesPerceived(*args): return _openbabel.OBMol_UnsetPartialChargesPerceived(*args)
    def UnsetImplicitValencePerceived(*args): return _openbabel.OBMol_UnsetImplicitValencePerceived(*args)
    def UnsetFlag(*args): return _openbabel.OBMol_UnsetFlag(*args)
    def DoTransformations(*args): return _openbabel.OBMol_DoTransformations(*args)
    ClassDescription = staticmethod(_openbabel.OBMol_ClassDescription)
    def Clear(*args): return _openbabel.OBMol_Clear(*args)
    def RenumberAtoms(*args): return _openbabel.OBMol_RenumberAtoms(*args)
    def ToInertialFrame(*args): return _openbabel.OBMol_ToInertialFrame(*args)
    def Translate(*args): return _openbabel.OBMol_Translate(*args)
    def Rotate(*args): return _openbabel.OBMol_Rotate(*args)
    def Kekulize(*args): return _openbabel.OBMol_Kekulize(*args)
    def PerceiveKekuleBonds(*args): return _openbabel.OBMol_PerceiveKekuleBonds(*args)
    def NewPerceiveKekuleBonds(*args): return _openbabel.OBMol_NewPerceiveKekuleBonds(*args)
    def DeleteHydrogen(*args): return _openbabel.OBMol_DeleteHydrogen(*args)
    def DeleteHydrogens(*args): return _openbabel.OBMol_DeleteHydrogens(*args)
    def DeleteNonPolarHydrogens(*args): return _openbabel.OBMol_DeleteNonPolarHydrogens(*args)
    def AddHydrogens(*args): return _openbabel.OBMol_AddHydrogens(*args)
    def AddPolarHydrogens(*args): return _openbabel.OBMol_AddPolarHydrogens(*args)
    def StripSalts(*args): return _openbabel.OBMol_StripSalts(*args)
    def ConvertDativeBonds(*args): return _openbabel.OBMol_ConvertDativeBonds(*args)
    def CorrectForPH(*args): return _openbabel.OBMol_CorrectForPH(*args)
    def AssignSpinMultiplicity(*args): return _openbabel.OBMol_AssignSpinMultiplicity(*args)
    def Center(*args): return _openbabel.OBMol_Center(*args)
    def SetTorsion(*args): return _openbabel.OBMol_SetTorsion(*args)
    def FindSSSR(*args): return _openbabel.OBMol_FindSSSR(*args)
    def FindRingAtomsAndBonds(*args): return _openbabel.OBMol_FindRingAtomsAndBonds(*args)
    def FindChiralCenters(*args): return _openbabel.OBMol_FindChiralCenters(*args)
    def FindChildren(*args): return _openbabel.OBMol_FindChildren(*args)
    def FindLargestFragment(*args): return _openbabel.OBMol_FindLargestFragment(*args)
    def ContigFragList(*args): return _openbabel.OBMol_ContigFragList(*args)
    def Align(*args): return _openbabel.OBMol_Align(*args)
    def ConnectTheDots(*args): return _openbabel.OBMol_ConnectTheDots(*args)
    def PerceiveBondOrders(*args): return _openbabel.OBMol_PerceiveBondOrders(*args)
    def FindTorsions(*args): return _openbabel.OBMol_FindTorsions(*args)
    def GetGTDVector(*args): return _openbabel.OBMol_GetGTDVector(*args)
    def GetGIVector(*args): return _openbabel.OBMol_GetGIVector(*args)
    def GetGIDVector(*args): return _openbabel.OBMol_GetGIDVector(*args)
    def Has2D(*args): return _openbabel.OBMol_Has2D(*args)
    def Has3D(*args): return _openbabel.OBMol_Has3D(*args)
    def HasNonZeroCoords(*args): return _openbabel.OBMol_HasNonZeroCoords(*args)
    def HasAromaticPerceived(*args): return _openbabel.OBMol_HasAromaticPerceived(*args)
    def HasSSSRPerceived(*args): return _openbabel.OBMol_HasSSSRPerceived(*args)
    def HasRingAtomsAndBondsPerceived(*args): return _openbabel.OBMol_HasRingAtomsAndBondsPerceived(*args)
    def HasAtomTypesPerceived(*args): return _openbabel.OBMol_HasAtomTypesPerceived(*args)
    def HasChiralityPerceived(*args): return _openbabel.OBMol_HasChiralityPerceived(*args)
    def HasPartialChargesPerceived(*args): return _openbabel.OBMol_HasPartialChargesPerceived(*args)
    def HasHybridizationPerceived(*args): return _openbabel.OBMol_HasHybridizationPerceived(*args)
    def HasImplicitValencePerceived(*args): return _openbabel.OBMol_HasImplicitValencePerceived(*args)
    def HasKekulePerceived(*args): return _openbabel.OBMol_HasKekulePerceived(*args)
    def HasClosureBondsPerceived(*args): return _openbabel.OBMol_HasClosureBondsPerceived(*args)
    def HasChainsPerceived(*args): return _openbabel.OBMol_HasChainsPerceived(*args)
    def HasHydrogensAdded(*args): return _openbabel.OBMol_HasHydrogensAdded(*args)
    def HasAromaticCorrected(*args): return _openbabel.OBMol_HasAromaticCorrected(*args)
    def IsCorrectedForPH(*args): return _openbabel.OBMol_IsCorrectedForPH(*args)
    def HasSpinMultiplicityAssigned(*args): return _openbabel.OBMol_HasSpinMultiplicityAssigned(*args)
    def IsChiral(*args): return _openbabel.OBMol_IsChiral(*args)
    def Empty(*args): return _openbabel.OBMol_Empty(*args)
    def NumConformers(*args): return _openbabel.OBMol_NumConformers(*args)
    def SetConformers(*args): return _openbabel.OBMol_SetConformers(*args)
    def AddConformer(*args): return _openbabel.OBMol_AddConformer(*args)
    def SetConformer(*args): return _openbabel.OBMol_SetConformer(*args)
    def CopyConformer(*args): return _openbabel.OBMol_CopyConformer(*args)
    def DeleteConformer(*args): return _openbabel.OBMol_DeleteConformer(*args)
    def GetConformer(*args): return _openbabel.OBMol_GetConformer(*args)
    def BeginConformer(*args): return _openbabel.OBMol_BeginConformer(*args)
    def NextConformer(*args): return _openbabel.OBMol_NextConformer(*args)
    def GetConformers(*args): return _openbabel.OBMol_GetConformers(*args)
    def BeginAtom(*args): return _openbabel.OBMol_BeginAtom(*args)
    def NextAtom(*args): return _openbabel.OBMol_NextAtom(*args)
    def BeginBond(*args): return _openbabel.OBMol_BeginBond(*args)
    def NextBond(*args): return _openbabel.OBMol_NextBond(*args)
    def BeginResidue(*args): return _openbabel.OBMol_BeginResidue(*args)
    def NextResidue(*args): return _openbabel.OBMol_NextResidue(*args)
    def BeginInternalCoord(*args): return _openbabel.OBMol_BeginInternalCoord(*args)
    def NextInternalCoord(*args): return _openbabel.OBMol_NextInternalCoord(*args)
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
    __del__ = lambda self : None;
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
    __del__ = lambda self : None;
    def GetAtomIdx(*args): return _openbabel.OBRTree_GetAtomIdx(*args)
    def PathToRoot(*args): return _openbabel.OBRTree_PathToRoot(*args)
OBRTree_swigregister = _openbabel.OBRTree_swigregister
OBRTree_swigregister(OBRTree)
tokenize = _openbabel.tokenize
ThrowError = _openbabel.ThrowError

class OBRing(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    _path = property(_openbabel.OBRing__path_get, _openbabel.OBRing__path_set)
    _pathset = property(_openbabel.OBRing__pathset_get, _openbabel.OBRing__pathset_set)
    def findCenterAndNormal(*args): return _openbabel.OBRing_findCenterAndNormal(*args)
    def __init__(self, *args): 
        _openbabel.OBRing_swiginit(self,_openbabel.new_OBRing(*args))
    def Size(*args): return _openbabel.OBRing_Size(*args)
    def PathSize(*args): return _openbabel.OBRing_PathSize(*args)
    def IsMember(*args): return _openbabel.OBRing_IsMember(*args)
    def IsAromatic(*args): return _openbabel.OBRing_IsAromatic(*args)
    def IsInRing(*args): return _openbabel.OBRing_IsInRing(*args)
    def SetParent(*args): return _openbabel.OBRing_SetParent(*args)
    def GetParent(*args): return _openbabel.OBRing_GetParent(*args)
    __swig_destroy__ = _openbabel.delete_OBRing
    __del__ = lambda self : None;
OBRing_swigregister = _openbabel.OBRing_swigregister
OBRing_swigregister(OBRing)

CompareRingSize = _openbabel.CompareRingSize
class OBRingSearch(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBRingSearch_swiginit(self,_openbabel.new_OBRingSearch(*args))
    __swig_destroy__ = _openbabel.delete_OBRingSearch
    __del__ = lambda self : None;
    def SortRings(*args): return _openbabel.OBRingSearch_SortRings(*args)
    def RemoveRedundant(*args): return _openbabel.OBRingSearch_RemoveRedundant(*args)
    def AddRingFromClosure(*args): return _openbabel.OBRingSearch_AddRingFromClosure(*args)
    def WriteRings(*args): return _openbabel.OBRingSearch_WriteRings(*args)
    def SaveUniqueRing(*args): return _openbabel.OBRingSearch_SaveUniqueRing(*args)
    def BeginRings(*args): return _openbabel.OBRingSearch_BeginRings(*args)
    def EndRings(*args): return _openbabel.OBRingSearch_EndRings(*args)
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
    __del__ = lambda self : None;
    def __init__(self, *args): 
        _openbabel.OBSmartsPattern_swiginit(self,_openbabel.new_OBSmartsPattern(*args))
    def NumMatches(*args): return _openbabel.OBSmartsPattern_NumMatches(*args)
    def NumAtoms(*args): return _openbabel.OBSmartsPattern_NumAtoms(*args)
    def NumBonds(*args): return _openbabel.OBSmartsPattern_NumBonds(*args)
    def GetAtomicNum(*args): return _openbabel.OBSmartsPattern_GetAtomicNum(*args)
    def GetBond(*args): return _openbabel.OBSmartsPattern_GetBond(*args)
    def GetCharge(*args): return _openbabel.OBSmartsPattern_GetCharge(*args)
    def GetSMARTS(*args): return _openbabel.OBSmartsPattern_GetSMARTS(*args)
    def GetVectorBinding(*args): return _openbabel.OBSmartsPattern_GetVectorBinding(*args)
    def Empty(*args): return _openbabel.OBSmartsPattern_Empty(*args)
    def IsValid(*args): return _openbabel.OBSmartsPattern_IsValid(*args)
    def Init(*args): return _openbabel.OBSmartsPattern_Init(*args)
    def WriteMapList(*args): return _openbabel.OBSmartsPattern_WriteMapList(*args)
    def Match(*args): return _openbabel.OBSmartsPattern_Match(*args)
    def RestrictedMatch(*args): return _openbabel.OBSmartsPattern_RestrictedMatch(*args)
    def GetMapList(*args): return _openbabel.OBSmartsPattern_GetMapList(*args)
    def GetUMapList(*args): return _openbabel.OBSmartsPattern_GetUMapList(*args)
    def BeginMList(*args): return _openbabel.OBSmartsPattern_BeginMList(*args)
    def EndMList(*args): return _openbabel.OBSmartsPattern_EndMList(*args)
OBSmartsPattern_swigregister = _openbabel.OBSmartsPattern_swigregister
OBSmartsPattern_swigregister(OBSmartsPattern)

class OBSSMatch(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBSSMatch_swiginit(self,_openbabel.new_OBSSMatch(*args))
    __swig_destroy__ = _openbabel.delete_OBSSMatch
    __del__ = lambda self : None;
    def Match(*args): return _openbabel.OBSSMatch_Match(*args)
OBSSMatch_swigregister = _openbabel.OBSSMatch_swigregister
OBSSMatch_swigregister(OBSSMatch)

SmartsLexReplace = _openbabel.SmartsLexReplace


