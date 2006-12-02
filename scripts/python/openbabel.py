import sys
if sys.platform.find("linux") != -1:
    import dl
    sys.setdlopenflags(sys.getdlopenflags() | dl.RTLD_GLOBAL)

# This file was created automatically by SWIG 1.3.30.
# Don't modify this file, modify the SWIG interface instead.

import _openbabel
import new
new_instancemethod = new.instancemethod
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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

class vectorData(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def iterator(*args): return _openbabel.vectorData_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _openbabel.vectorData___nonzero__(*args)
    def __len__(*args): return _openbabel.vectorData___len__(*args)
    def pop(*args): return _openbabel.vectorData_pop(*args)
    def __getslice__(*args): return _openbabel.vectorData___getslice__(*args)
    def __setslice__(*args): return _openbabel.vectorData___setslice__(*args)
    def __delslice__(*args): return _openbabel.vectorData___delslice__(*args)
    def __delitem__(*args): return _openbabel.vectorData___delitem__(*args)
    def __getitem__(*args): return _openbabel.vectorData___getitem__(*args)
    def __setitem__(*args): return _openbabel.vectorData___setitem__(*args)
    def append(*args): return _openbabel.vectorData_append(*args)
    def empty(*args): return _openbabel.vectorData_empty(*args)
    def size(*args): return _openbabel.vectorData_size(*args)
    def clear(*args): return _openbabel.vectorData_clear(*args)
    def swap(*args): return _openbabel.vectorData_swap(*args)
    def get_allocator(*args): return _openbabel.vectorData_get_allocator(*args)
    def begin(*args): return _openbabel.vectorData_begin(*args)
    def end(*args): return _openbabel.vectorData_end(*args)
    def rbegin(*args): return _openbabel.vectorData_rbegin(*args)
    def rend(*args): return _openbabel.vectorData_rend(*args)
    def pop_back(*args): return _openbabel.vectorData_pop_back(*args)
    def erase(*args): return _openbabel.vectorData_erase(*args)
    def __init__(self, *args): 
        _openbabel.vectorData_swiginit(self,_openbabel.new_vectorData(*args))
    def push_back(*args): return _openbabel.vectorData_push_back(*args)
    def front(*args): return _openbabel.vectorData_front(*args)
    def back(*args): return _openbabel.vectorData_back(*args)
    def assign(*args): return _openbabel.vectorData_assign(*args)
    def resize(*args): return _openbabel.vectorData_resize(*args)
    def insert(*args): return _openbabel.vectorData_insert(*args)
    def reserve(*args): return _openbabel.vectorData_reserve(*args)
    def capacity(*args): return _openbabel.vectorData_capacity(*args)
    __swig_destroy__ = _openbabel.delete_vectorData
    __del__ = lambda self : None;
vectorData_swigregister = _openbabel.vectorData_swigregister
vectorData_swigregister(vectorData)

class OBGlobalDataBase(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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

OpenDatafile = _openbabel.OpenDatafile
FILE_SEP_CHAR = _openbabel.FILE_SEP_CHAR
class DoubleType(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    hi = _swig_property(_openbabel.DoubleType_hi_get, _openbabel.DoubleType_hi_set)
    lo = _swig_property(_openbabel.DoubleType_lo_get, _openbabel.DoubleType_lo_set)
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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

class OBStopwatch(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBSqrtTbl_swiginit(self,_openbabel.new_OBSqrtTbl(*args))
    __swig_destroy__ = _openbabel.delete_OBSqrtTbl
    __del__ = lambda self : None;
    def Sqrt(*args): return _openbabel.OBSqrtTbl_Sqrt(*args)
    def Init(*args): return _openbabel.OBSqrtTbl_Init(*args)
OBSqrtTbl_swigregister = _openbabel.OBSqrtTbl_swigregister
OBSqrtTbl_swigregister(OBSqrtTbl)

class vector3(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.vector3_swiginit(self,_openbabel.new_vector3(*args))
    def Set(*args): return _openbabel.vector3_Set(*args)
    def SetX(*args): return _openbabel.vector3_SetX(*args)
    def SetY(*args): return _openbabel.vector3_SetY(*args)
    def SetZ(*args): return _openbabel.vector3_SetZ(*args)
    def Get(*args): return _openbabel.vector3_Get(*args)
    def AsArray(*args): return _openbabel.vector3_AsArray(*args)
    def __iadd__(*args): return _openbabel.vector3___iadd__(*args)
    def __isub__(*args): return _openbabel.vector3___isub__(*args)
    def __idiv__(*args): return _openbabel.vector3___idiv__(*args)
    def __imul__(*args): return _openbabel.vector3___imul__(*args)
    def randomUnitVector(*args): return _openbabel.vector3_randomUnitVector(*args)
    def normalize(*args): return _openbabel.vector3_normalize(*args)
    def CanBeNormalized(*args): return _openbabel.vector3_CanBeNormalized(*args)
    def length_2(*args): return _openbabel.vector3_length_2(*args)
    def length(*args): return _openbabel.vector3_length(*args)
    def x(*args): return _openbabel.vector3_x(*args)
    def y(*args): return _openbabel.vector3_y(*args)
    def z(*args): return _openbabel.vector3_z(*args)
    def __eq__(*args): return _openbabel.vector3___eq__(*args)
    def __ne__(*args): return _openbabel.vector3___ne__(*args)
    def IsApprox(*args): return _openbabel.vector3_IsApprox(*args)
    def distSq(*args): return _openbabel.vector3_distSq(*args)
    def createOrthoVector(*args): return _openbabel.vector3_createOrthoVector(*args)
    __swig_destroy__ = _openbabel.delete_vector3
    __del__ = lambda self : None;
vector3_swigregister = _openbabel.vector3_swigregister
vector3_swigregister(vector3)

__lshift__ = _openbabel.__lshift__
__add__ = _openbabel.__add__
__div__ = _openbabel.__div__
dot = _openbabel.dot
cross = _openbabel.cross
vectorAngle = _openbabel.vectorAngle
CalcTorsionAngle = _openbabel.CalcTorsionAngle
Point2Plane = _openbabel.Point2Plane
Trim = _openbabel.Trim
class OBGenericData(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBGenericData_swiginit(self,_openbabel.new_OBGenericData(*args))
    def Clone(*args): return _openbabel.OBGenericData_Clone(*args)
    __swig_destroy__ = _openbabel.delete_OBGenericData
    __del__ = lambda self : None;
    def SetAttribute(*args): return _openbabel.OBGenericData_SetAttribute(*args)
    def GetAttribute(*args): return _openbabel.OBGenericData_GetAttribute(*args)
    def GetDataType(*args): return _openbabel.OBGenericData_GetDataType(*args)
    def GetValue(*args): return _openbabel.OBGenericData_GetValue(*args)
OBGenericData_swigregister = _openbabel.OBGenericData_swigregister
OBGenericData_swigregister(OBGenericData)
__sub__ = _openbabel.__sub__
__mul__ = _openbabel.__mul__
cvar = _openbabel.cvar
VZero = cvar.VZero
VX = cvar.VX
VY = cvar.VY
VZ = cvar.VZ
UndefinedData = cvar.UndefinedData
PairData = cvar.PairData
EnergyData = cvar.EnergyData
CommentData = cvar.CommentData
ConformerData = cvar.ConformerData
ExternalBondData = cvar.ExternalBondData
RotamerList = cvar.RotamerList
VirtualBondData = cvar.VirtualBondData
RingData = cvar.RingData
TorsionData = cvar.TorsionData
AngleData = cvar.AngleData
SerialNums = cvar.SerialNums
UnitCell = cvar.UnitCell
SpinData = cvar.SpinData
ChargeData = cvar.ChargeData
SymmetryData = cvar.SymmetryData
ChiralData = cvar.ChiralData
OccupationData = cvar.OccupationData
DensityData = cvar.DensityData
ElectronicData = cvar.ElectronicData
VibrationData = cvar.VibrationData
RotationData = cvar.RotationData
NuclearData = cvar.NuclearData
SetData = cvar.SetData
CustomData0 = cvar.CustomData0
CustomData1 = cvar.CustomData1
CustomData2 = cvar.CustomData2
CustomData3 = cvar.CustomData3
CustomData4 = cvar.CustomData4
CustomData5 = cvar.CustomData5
CustomData6 = cvar.CustomData6
CustomData7 = cvar.CustomData7
CustomData8 = cvar.CustomData8
CustomData9 = cvar.CustomData9
CustomData10 = cvar.CustomData10
CustomData11 = cvar.CustomData11
CustomData12 = cvar.CustomData12
CustomData13 = cvar.CustomData13
CustomData14 = cvar.CustomData14
CustomData15 = cvar.CustomData15

class OBCommentData(OBGenericData):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBCommentData_swiginit(self,_openbabel.new_OBCommentData(*args))
    def Clone(*args): return _openbabel.OBCommentData_Clone(*args)
    def SetData(*args): return _openbabel.OBCommentData_SetData(*args)
    def GetData(*args): return _openbabel.OBCommentData_GetData(*args)
    def GetValue(*args): return _openbabel.OBCommentData_GetValue(*args)
    __swig_destroy__ = _openbabel.delete_OBCommentData
    __del__ = lambda self : None;
OBCommentData_swigregister = _openbabel.OBCommentData_swigregister
OBCommentData_swigregister(OBCommentData)

class OBExternalBond(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBExternalBond_swiginit(self,_openbabel.new_OBExternalBond(*args))
    __swig_destroy__ = _openbabel.delete_OBExternalBond
    __del__ = lambda self : None;
    def GetIdx(*args): return _openbabel.OBExternalBond_GetIdx(*args)
    def GetAtom(*args): return _openbabel.OBExternalBond_GetAtom(*args)
    def GetBond(*args): return _openbabel.OBExternalBond_GetBond(*args)
    def SetIdx(*args): return _openbabel.OBExternalBond_SetIdx(*args)
    def SetAtom(*args): return _openbabel.OBExternalBond_SetAtom(*args)
    def SetBond(*args): return _openbabel.OBExternalBond_SetBond(*args)
OBExternalBond_swigregister = _openbabel.OBExternalBond_swigregister
OBExternalBond_swigregister(OBExternalBond)

class OBExternalBondData(OBGenericData):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBExternalBondData_swiginit(self,_openbabel.new_OBExternalBondData(*args))
    def Clone(*args): return _openbabel.OBExternalBondData_Clone(*args)
    def SetData(*args): return _openbabel.OBExternalBondData_SetData(*args)
    def GetData(*args): return _openbabel.OBExternalBondData_GetData(*args)
    __swig_destroy__ = _openbabel.delete_OBExternalBondData
    __del__ = lambda self : None;
OBExternalBondData_swigregister = _openbabel.OBExternalBondData_swigregister
OBExternalBondData_swigregister(OBExternalBondData)

class OBPairData(OBGenericData):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBPairData_swiginit(self,_openbabel.new_OBPairData(*args))
    def Clone(*args): return _openbabel.OBPairData_Clone(*args)
    def SetValue(*args): return _openbabel.OBPairData_SetValue(*args)
    def GetValue(*args): return _openbabel.OBPairData_GetValue(*args)
    __swig_destroy__ = _openbabel.delete_OBPairData
    __del__ = lambda self : None;
OBPairData_swigregister = _openbabel.OBPairData_swigregister
OBPairData_swigregister(OBPairData)

class OBSetData(OBGenericData):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBSetData_swiginit(self,_openbabel.new_OBSetData(*args))
    def Clone(*args): return _openbabel.OBSetData_Clone(*args)
    def AddData(*args): return _openbabel.OBSetData_AddData(*args)
    def SetData(*args): return _openbabel.OBSetData_SetData(*args)
    def GetData(*args): return _openbabel.OBSetData_GetData(*args)
    def GetBegin(*args): return _openbabel.OBSetData_GetBegin(*args)
    def GetEnd(*args): return _openbabel.OBSetData_GetEnd(*args)
    def DeleteData(*args): return _openbabel.OBSetData_DeleteData(*args)
    __swig_destroy__ = _openbabel.delete_OBSetData
    __del__ = lambda self : None;
OBSetData_swigregister = _openbabel.OBSetData_swigregister
OBSetData_swigregister(OBSetData)

class OBVirtualBond(OBGenericData):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def Clone(*args): return _openbabel.OBVirtualBond_Clone(*args)
    def __init__(self, *args): 
        _openbabel.OBVirtualBond_swiginit(self,_openbabel.new_OBVirtualBond(*args))
    def GetBgn(*args): return _openbabel.OBVirtualBond_GetBgn(*args)
    def GetEnd(*args): return _openbabel.OBVirtualBond_GetEnd(*args)
    def GetOrder(*args): return _openbabel.OBVirtualBond_GetOrder(*args)
    def GetStereo(*args): return _openbabel.OBVirtualBond_GetStereo(*args)
    __swig_destroy__ = _openbabel.delete_OBVirtualBond
    __del__ = lambda self : None;
OBVirtualBond_swigregister = _openbabel.OBVirtualBond_swigregister
OBVirtualBond_swigregister(OBVirtualBond)

class OBRingData(OBGenericData):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBRingData_swiginit(self,_openbabel.new_OBRingData(*args))
    def Clone(*args): return _openbabel.OBRingData_Clone(*args)
    __swig_destroy__ = _openbabel.delete_OBRingData
    __del__ = lambda self : None;
    def SetData(*args): return _openbabel.OBRingData_SetData(*args)
    def PushBack(*args): return _openbabel.OBRingData_PushBack(*args)
    def GetData(*args): return _openbabel.OBRingData_GetData(*args)
OBRingData_swigregister = _openbabel.OBRingData_swigregister
OBRingData_swigregister(OBRingData)

class OBUnitCell(OBGenericData):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    Undefined = _openbabel.OBUnitCell_Undefined
    Triclinic = _openbabel.OBUnitCell_Triclinic
    Monoclinic = _openbabel.OBUnitCell_Monoclinic
    Orthorhombic = _openbabel.OBUnitCell_Orthorhombic
    Tetragonal = _openbabel.OBUnitCell_Tetragonal
    Rhombohedral = _openbabel.OBUnitCell_Rhombohedral
    Hexagonal = _openbabel.OBUnitCell_Hexagonal
    Cubic = _openbabel.OBUnitCell_Cubic
    def __init__(self, *args): 
        _openbabel.OBUnitCell_swiginit(self,_openbabel.new_OBUnitCell(*args))
    def Clone(*args): return _openbabel.OBUnitCell_Clone(*args)
    __swig_destroy__ = _openbabel.delete_OBUnitCell
    __del__ = lambda self : None;
    def SetData(*args): return _openbabel.OBUnitCell_SetData(*args)
    def SetOffset(*args): return _openbabel.OBUnitCell_SetOffset(*args)
    def SetSpaceGroup(*args): return _openbabel.OBUnitCell_SetSpaceGroup(*args)
    def SetLatticeType(*args): return _openbabel.OBUnitCell_SetLatticeType(*args)
    def GetA(*args): return _openbabel.OBUnitCell_GetA(*args)
    def GetB(*args): return _openbabel.OBUnitCell_GetB(*args)
    def GetC(*args): return _openbabel.OBUnitCell_GetC(*args)
    def GetAlpha(*args): return _openbabel.OBUnitCell_GetAlpha(*args)
    def GetBeta(*args): return _openbabel.OBUnitCell_GetBeta(*args)
    def GetGamma(*args): return _openbabel.OBUnitCell_GetGamma(*args)
    def GetOffset(*args): return _openbabel.OBUnitCell_GetOffset(*args)
    def GetSpaceGroup(*args): return _openbabel.OBUnitCell_GetSpaceGroup(*args)
    def GetLatticeType(*args): return _openbabel.OBUnitCell_GetLatticeType(*args)
    def GetCellVectors(*args): return _openbabel.OBUnitCell_GetCellVectors(*args)
    def GetCellMatrix(*args): return _openbabel.OBUnitCell_GetCellMatrix(*args)
    def GetOrthoMatrix(*args): return _openbabel.OBUnitCell_GetOrthoMatrix(*args)
    def GetFractionalMatrix(*args): return _openbabel.OBUnitCell_GetFractionalMatrix(*args)
    def GetSpaceGroupNumber(*args): return _openbabel.OBUnitCell_GetSpaceGroupNumber(*args)
    def GetCellVolume(*args): return _openbabel.OBUnitCell_GetCellVolume(*args)
OBUnitCell_swigregister = _openbabel.OBUnitCell_swigregister
OBUnitCell_swigregister(OBUnitCell)

class OBConformerData(OBGenericData):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBConformerData_swiginit(self,_openbabel.new_OBConformerData(*args))
    def Clone(*args): return _openbabel.OBConformerData_Clone(*args)
    __swig_destroy__ = _openbabel.delete_OBConformerData
    __del__ = lambda self : None;
    def SetDimension(*args): return _openbabel.OBConformerData_SetDimension(*args)
    def SetEnergies(*args): return _openbabel.OBConformerData_SetEnergies(*args)
    def SetForces(*args): return _openbabel.OBConformerData_SetForces(*args)
    def SetVelocities(*args): return _openbabel.OBConformerData_SetVelocities(*args)
    def SetDisplacements(*args): return _openbabel.OBConformerData_SetDisplacements(*args)
    def SetData(*args): return _openbabel.OBConformerData_SetData(*args)
    def GetDimension(*args): return _openbabel.OBConformerData_GetDimension(*args)
    def GetEnergies(*args): return _openbabel.OBConformerData_GetEnergies(*args)
    def GetForces(*args): return _openbabel.OBConformerData_GetForces(*args)
    def GetVelocities(*args): return _openbabel.OBConformerData_GetVelocities(*args)
    def GetDisplacements(*args): return _openbabel.OBConformerData_GetDisplacements(*args)
    def GetData(*args): return _openbabel.OBConformerData_GetData(*args)
OBConformerData_swigregister = _openbabel.OBConformerData_swigregister
OBConformerData_swigregister(OBConformerData)

class OBSymmetryData(OBGenericData):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBSymmetryData_swiginit(self,_openbabel.new_OBSymmetryData(*args))
    def Clone(*args): return _openbabel.OBSymmetryData_Clone(*args)
    __swig_destroy__ = _openbabel.delete_OBSymmetryData
    __del__ = lambda self : None;
    def SetData(*args): return _openbabel.OBSymmetryData_SetData(*args)
    def SetPointGroup(*args): return _openbabel.OBSymmetryData_SetPointGroup(*args)
    def SetSpaceGroup(*args): return _openbabel.OBSymmetryData_SetSpaceGroup(*args)
    def GetPointGroup(*args): return _openbabel.OBSymmetryData_GetPointGroup(*args)
    def GetSpaceGroup(*args): return _openbabel.OBSymmetryData_GetSpaceGroup(*args)
OBSymmetryData_swigregister = _openbabel.OBSymmetryData_swigregister
OBSymmetryData_swigregister(OBSymmetryData)

class OBTorsion(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBTorsion_swiginit(self,_openbabel.new_OBTorsion(*args))
    __swig_destroy__ = _openbabel.delete_OBTorsion
    __del__ = lambda self : None;
    def Clear(*args): return _openbabel.OBTorsion_Clear(*args)
    def Empty(*args): return _openbabel.OBTorsion_Empty(*args)
    def AddTorsion(*args): return _openbabel.OBTorsion_AddTorsion(*args)
    def SetAngle(*args): return _openbabel.OBTorsion_SetAngle(*args)
    def SetData(*args): return _openbabel.OBTorsion_SetData(*args)
    def GetAngle(*args): return _openbabel.OBTorsion_GetAngle(*args)
    def GetBondIdx(*args): return _openbabel.OBTorsion_GetBondIdx(*args)
    def GetSize(*args): return _openbabel.OBTorsion_GetSize(*args)
    def GetBC(*args): return _openbabel.OBTorsion_GetBC(*args)
    def GetADs(*args): return _openbabel.OBTorsion_GetADs(*args)
    def IsProtonRotor(*args): return _openbabel.OBTorsion_IsProtonRotor(*args)
OBTorsion_swigregister = _openbabel.OBTorsion_swigregister
OBTorsion_swigregister(OBTorsion)

class OBTorsionData(OBGenericData):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    def __init__(self): raise AttributeError, "No constructor defined"
    __repr__ = _swig_repr
    def Clone(*args): return _openbabel.OBTorsionData_Clone(*args)
    def Clear(*args): return _openbabel.OBTorsionData_Clear(*args)
    def GetData(*args): return _openbabel.OBTorsionData_GetData(*args)
    def GetSize(*args): return _openbabel.OBTorsionData_GetSize(*args)
    def SetData(*args): return _openbabel.OBTorsionData_SetData(*args)
    def FillTorsionArray(*args): return _openbabel.OBTorsionData_FillTorsionArray(*args)
    __swig_destroy__ = _openbabel.delete_OBTorsionData
    __del__ = lambda self : None;
OBTorsionData_swigregister = _openbabel.OBTorsionData_swigregister
OBTorsionData_swigregister(OBTorsionData)

class OBAngle(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBAngle_swiginit(self,_openbabel.new_OBAngle(*args))
    __swig_destroy__ = _openbabel.delete_OBAngle
    __del__ = lambda self : None;
    def __eq__(*args): return _openbabel.OBAngle___eq__(*args)
    def Clear(*args): return _openbabel.OBAngle_Clear(*args)
    def GetAngle(*args): return _openbabel.OBAngle_GetAngle(*args)
    def SetAngle(*args): return _openbabel.OBAngle_SetAngle(*args)
    def SetAtoms(*args): return _openbabel.OBAngle_SetAtoms(*args)
OBAngle_swigregister = _openbabel.OBAngle_swigregister
OBAngle_swigregister(OBAngle)

class OBAngleData(OBGenericData):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    def __init__(self): raise AttributeError, "No constructor defined"
    __repr__ = _swig_repr
    def Clone(*args): return _openbabel.OBAngleData_Clone(*args)
    def Clear(*args): return _openbabel.OBAngleData_Clear(*args)
    def FillAngleArray(*args): return _openbabel.OBAngleData_FillAngleArray(*args)
    def SetData(*args): return _openbabel.OBAngleData_SetData(*args)
    def GetSize(*args): return _openbabel.OBAngleData_GetSize(*args)
    __swig_destroy__ = _openbabel.delete_OBAngleData
    __del__ = lambda self : None;
OBAngleData_swigregister = _openbabel.OBAngleData_swigregister
OBAngleData_swigregister(OBAngleData)

output = _openbabel.output
input = _openbabel.input
calcvolume = _openbabel.calcvolume
class OBChiralData(OBGenericData):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def GetAtom4Refs(*args): return _openbabel.OBChiralData_GetAtom4Refs(*args)
    def GetAtomRef(*args): return _openbabel.OBChiralData_GetAtomRef(*args)
    def __init__(self, *args): 
        _openbabel.OBChiralData_swiginit(self,_openbabel.new_OBChiralData(*args))
    def Clone(*args): return _openbabel.OBChiralData_Clone(*args)
    __swig_destroy__ = _openbabel.delete_OBChiralData
    __del__ = lambda self : None;
    def Clear(*args): return _openbabel.OBChiralData_Clear(*args)
    def SetAtom4Refs(*args): return _openbabel.OBChiralData_SetAtom4Refs(*args)
    def AddAtomRef(*args): return _openbabel.OBChiralData_AddAtomRef(*args)
    def GetSize(*args): return _openbabel.OBChiralData_GetSize(*args)
OBChiralData_swigregister = _openbabel.OBChiralData_swigregister
OBChiralData_swigregister(OBChiralData)

class OBSerialNums(OBGenericData):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBSerialNums_swiginit(self,_openbabel.new_OBSerialNums(*args))
    def Clone(*args): return _openbabel.OBSerialNums_Clone(*args)
    def GetData(*args): return _openbabel.OBSerialNums_GetData(*args)
    def SetData(*args): return _openbabel.OBSerialNums_SetData(*args)
    __swig_destroy__ = _openbabel.delete_OBSerialNums
    __del__ = lambda self : None;
OBSerialNums_swigregister = _openbabel.OBSerialNums_swigregister
OBSerialNums_swigregister(OBSerialNums)

class OBBase(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    __swig_destroy__ = _openbabel.delete_OBBase
    __del__ = lambda self : None;
    def DoTransformations(*args): return _openbabel.OBBase_DoTransformations(*args)
    ClassDescription = staticmethod(_openbabel.OBBase_ClassDescription)
    def HasData(*args): return _openbabel.OBBase_HasData(*args)
    def DeleteData(*args): return _openbabel.OBBase_DeleteData(*args)
    def SetData(*args): return _openbabel.OBBase_SetData(*args)
    def DataSize(*args): return _openbabel.OBBase_DataSize(*args)
    def GetData(*args): return _openbabel.OBBase_GetData(*args)
    def BeginData(*args): return _openbabel.OBBase_BeginData(*args)
    def EndData(*args): return _openbabel.OBBase_EndData(*args)
    def __init__(self, *args): 
        _openbabel.OBBase_swiginit(self,_openbabel.new_OBBase(*args))
OBBase_swigregister = _openbabel.OBBase_swigregister
OBBase_swigregister(OBBase)
OBBase_ClassDescription = _openbabel.OBBase_ClassDescription

obError = _openbabel.obError
obWarning = _openbabel.obWarning
obInfo = _openbabel.obInfo
obAuditMsg = _openbabel.obAuditMsg
obDebug = _openbabel.obDebug
class OBError:
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OBError, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OBError, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBError_swiginit(self,_openbabel.new_OBError(*args))
    def message(*args): return _openbabel.OBError_message(*args)
    def GetMethod(*args): return _openbabel.OBError_GetMethod(*args)
    def GetError(*args): return _openbabel.OBError_GetError(*args)
    def GetExplanation(*args): return _openbabel.OBError_GetExplanation(*args)
    def GetPossibleCause(*args): return _openbabel.OBError_GetPossibleCause(*args)
    def GetSuggestedRemedy(*args): return _openbabel.OBError_GetSuggestedRemedy(*args)
    def GetLevel(*args): return _openbabel.OBError_GetLevel(*args)
    __swig_destroy__ = _openbabel.delete_OBError
    __del__ = lambda self : None;
OBError_swigregister = _openbabel.OBError_swigregister
OBError_swigregister(OBError)

class OBMessageHandler(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBMessageHandler_swiginit(self,_openbabel.new_OBMessageHandler(*args))
    __swig_destroy__ = _openbabel.delete_OBMessageHandler
    __del__ = lambda self : None;
    def ThrowError(*args): return _openbabel.OBMessageHandler_ThrowError(*args)
    def GetMessagesOfLevel(*args): return _openbabel.OBMessageHandler_GetMessagesOfLevel(*args)
    def StartLogging(*args): return _openbabel.OBMessageHandler_StartLogging(*args)
    def StopLogging(*args): return _openbabel.OBMessageHandler_StopLogging(*args)
    def SetMaxLogEntries(*args): return _openbabel.OBMessageHandler_SetMaxLogEntries(*args)
    def GetMaxLogEntries(*args): return _openbabel.OBMessageHandler_GetMaxLogEntries(*args)
    def ClearLog(*args): return _openbabel.OBMessageHandler_ClearLog(*args)
    def SetOutputLevel(*args): return _openbabel.OBMessageHandler_SetOutputLevel(*args)
    def GetOutputLevel(*args): return _openbabel.OBMessageHandler_GetOutputLevel(*args)
    def SetOutputStream(*args): return _openbabel.OBMessageHandler_SetOutputStream(*args)
    def GetOutputStream(*args): return _openbabel.OBMessageHandler_GetOutputStream(*args)
    def StartErrorWrap(*args): return _openbabel.OBMessageHandler_StartErrorWrap(*args)
    def StopErrorWrap(*args): return _openbabel.OBMessageHandler_StopErrorWrap(*args)
    def GetErrorMessageCount(*args): return _openbabel.OBMessageHandler_GetErrorMessageCount(*args)
    def GetWarningMessageCount(*args): return _openbabel.OBMessageHandler_GetWarningMessageCount(*args)
    def GetInfoMessageCount(*args): return _openbabel.OBMessageHandler_GetInfoMessageCount(*args)
    def GetAuditMessageCount(*args): return _openbabel.OBMessageHandler_GetAuditMessageCount(*args)
    def GetDebugMessageCount(*args): return _openbabel.OBMessageHandler_GetDebugMessageCount(*args)
    def GetMessageSummary(*args): return _openbabel.OBMessageHandler_GetMessageSummary(*args)
OBMessageHandler_swigregister = _openbabel.OBMessageHandler_swigregister
OBMessageHandler_swigregister(OBMessageHandler)

class obLogBuf(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    __swig_destroy__ = _openbabel.delete_obLogBuf
    __del__ = lambda self : None;
    def __init__(self, *args): 
        _openbabel.obLogBuf_swiginit(self,_openbabel.new_obLogBuf(*args))
obLogBuf_swigregister = _openbabel.obLogBuf_swigregister
obLogBuf_swigregister(obLogBuf)

class OBFormat(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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

class CharPtrLess(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __call__(*args): return _openbabel.CharPtrLess___call__(*args)
    def __init__(self, *args): 
        _openbabel.CharPtrLess_swiginit(self,_openbabel.new_CharPtrLess(*args))
    __swig_destroy__ = _openbabel.delete_CharPtrLess
    __del__ = lambda self : None;
CharPtrLess_swigregister = _openbabel.CharPtrLess_swigregister
CharPtrLess_swigregister(CharPtrLess)

class OBConversion(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    def CloseOutFile(*args): return _openbabel.OBConversion_CloseOutFile(*args)
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
ZEROATOMSOK = _openbabel.ZEROATOMSOK
NOTWRITABLE = _openbabel.NOTWRITABLE
WRITEONEONLY = _openbabel.WRITEONEONLY
WRITEBINARY = _openbabel.WRITEBINARY
DEFAULTFORMAT = _openbabel.DEFAULTFORMAT
class OBResidue(OBBase):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    def BeginAtoms(*args): return _openbabel.OBResidue_BeginAtoms(*args)
    def EndAtoms(*args): return _openbabel.OBResidue_EndAtoms(*args)
    def BeginAtom(*args): return _openbabel.OBResidue_BeginAtom(*args)
    def NextAtom(*args): return _openbabel.OBResidue_NextAtom(*args)
OBResidue_swigregister = _openbabel.OBResidue_swigregister
OBResidue_swigregister(OBResidue)

MAXSETNO = _openbabel.MAXSETNO
MAXELEM = _openbabel.MAXELEM
MINELEM = _openbabel.MINELEM
MAXRES = _openbabel.MAXRES
MINRES = _openbabel.MINRES
AA_ALA = _openbabel.AA_ALA
AA_GLY = _openbabel.AA_GLY
AA_LEU = _openbabel.AA_LEU
AA_SER = _openbabel.AA_SER
AA_VAL = _openbabel.AA_VAL
AA_THR = _openbabel.AA_THR
AA_LYS = _openbabel.AA_LYS
AA_ASP = _openbabel.AA_ASP
AA_ILE = _openbabel.AA_ILE
AA_ASN = _openbabel.AA_ASN
AA_GLU = _openbabel.AA_GLU
AA_PRO = _openbabel.AA_PRO
AA_ARG = _openbabel.AA_ARG
AA_PHE = _openbabel.AA_PHE
AA_GLN = _openbabel.AA_GLN
AA_TYR = _openbabel.AA_TYR
AA_HIS = _openbabel.AA_HIS
AA_CYS = _openbabel.AA_CYS
AA_MET = _openbabel.AA_MET
AA_TRP = _openbabel.AA_TRP
class OBInternalCoord(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    _a = _swig_property(_openbabel.OBInternalCoord__a_get, _openbabel.OBInternalCoord__a_set)
    _b = _swig_property(_openbabel.OBInternalCoord__b_get, _openbabel.OBInternalCoord__b_set)
    _c = _swig_property(_openbabel.OBInternalCoord__c_get, _openbabel.OBInternalCoord__c_set)
    _dst = _swig_property(_openbabel.OBInternalCoord__dst_get, _openbabel.OBInternalCoord__dst_set)
    _ang = _swig_property(_openbabel.OBInternalCoord__ang_get, _openbabel.OBInternalCoord__ang_set)
    _tor = _swig_property(_openbabel.OBInternalCoord__tor_get, _openbabel.OBInternalCoord__tor_set)
    def __init__(self, *args): 
        _openbabel.OBInternalCoord_swiginit(self,_openbabel.new_OBInternalCoord(*args))
    __swig_destroy__ = _openbabel.delete_OBInternalCoord
    __del__ = lambda self : None;
OBInternalCoord_swigregister = _openbabel.OBInternalCoord_swigregister
OBInternalCoord_swigregister(OBInternalCoord)
ACIDIC = cvar.ACIDIC
ACYCLIC = cvar.ACYCLIC
ALIPHATIC = cvar.ALIPHATIC
AROMATIC = cvar.AROMATIC
BASIC = cvar.BASIC
BURIED = cvar.BURIED
CHARGED = cvar.CHARGED
CYCLIC = cvar.CYCLIC
HYDROPHOBIC = cvar.HYDROPHOBIC
LARGE = cvar.LARGE
MEDIUM = cvar.MEDIUM
NEGATIVE = cvar.NEGATIVE
NEUTRAL = cvar.NEUTRAL
POLAR = cvar.POLAR
POSITIVE = cvar.POSITIVE
SMALL = cvar.SMALL
SURFACE = cvar.SURFACE
ALPHA_CARBON = cvar.ALPHA_CARBON
AMINO_BACKBONE = cvar.AMINO_BACKBONE
BACKBONE = cvar.BACKBONE
CYSTEINE_SULPHUR = cvar.CYSTEINE_SULPHUR
LIGAND = cvar.LIGAND
NUCLEIC_BACKBONE = cvar.NUCLEIC_BACKBONE
SHAPELY_BACKBONE = cvar.SHAPELY_BACKBONE
SHAPELY_SPECIAL = cvar.SHAPELY_SPECIAL
SIDECHAIN = cvar.SIDECHAIN
SUGAR_PHOSPHATE = cvar.SUGAR_PHOSPHATE
ALA = cvar.ALA
GLY = cvar.GLY
LEU = cvar.LEU
SER = cvar.SER
VAL = cvar.VAL
THR = cvar.THR
LYS = cvar.LYS
ASP = cvar.ASP
ILE = cvar.ILE
ASN = cvar.ASN
GLU = cvar.GLU
PRO = cvar.PRO
ARG = cvar.ARG
PHE = cvar.PHE
GLN = cvar.GLN
TYR = cvar.TYR
HIS = cvar.HIS
CYS = cvar.CYS
MET = cvar.MET
TRP = cvar.TRP
ASX = cvar.ASX
GLX = cvar.GLX
PCA = cvar.PCA
HYP = cvar.HYP
A = cvar.A
C = cvar.C
G = cvar.G
T = cvar.T
U = cvar.U
UPLUS = cvar.UPLUS
I = cvar.I
_1MA = cvar._1MA
_5MC = cvar._5MC
OMC = cvar.OMC
_1MG = cvar._1MG
_2MG = cvar._2MG
M2G = cvar.M2G
_7MG = cvar._7MG
OMG = cvar.OMG
YG = cvar.YG
H2U = cvar.H2U
_5MU = cvar._5MU
PSU = cvar.PSU
UNK = cvar.UNK
ACE = cvar.ACE
FOR = cvar.FOR
HOH = cvar.HOH
DOD = cvar.DOD
SO4 = cvar.SO4
PO4 = cvar.PO4
NAD = cvar.NAD
COA = cvar.COA
NAP = cvar.NAP
NDP = cvar.NDP
AMINO = cvar.AMINO
AMINO_NUCLEO = cvar.AMINO_NUCLEO
COENZYME = cvar.COENZYME
ION = cvar.ION
NUCLEO = cvar.NUCLEO
PROTEIN = cvar.PROTEIN
PURINE = cvar.PURINE
PYRIMIDINE = cvar.PYRIMIDINE
SOLVENT = cvar.SOLVENT
WATER = cvar.WATER

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
OB_ATOM_HAS_NO_H = _openbabel.OB_ATOM_HAS_NO_H
class OBAtom(OBBase):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    Visit = _swig_property(_openbabel.OBAtom_Visit_get, _openbabel.OBAtom_Visit_set)
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
    def SetParent(*args): return _openbabel.OBAtom_SetParent(*args)
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
    def x(*args): return _openbabel.OBAtom_x(*args)
    def GetY(*args): return _openbabel.OBAtom_GetY(*args)
    def y(*args): return _openbabel.OBAtom_y(*args)
    def GetZ(*args): return _openbabel.OBAtom_GetZ(*args)
    def z(*args): return _openbabel.OBAtom_z(*args)
    def GetCoordinate(*args): return _openbabel.OBAtom_GetCoordinate(*args)
    def GetVector(*args): return _openbabel.OBAtom_GetVector(*args)
    def GetPartialCharge(*args): return _openbabel.OBAtom_GetPartialCharge(*args)
    def GetResidue(*args): return _openbabel.OBAtom_GetResidue(*args)
    def GetParent(*args): return _openbabel.OBAtom_GetParent(*args)
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
    def CountRingBonds(*args): return _openbabel.OBAtom_CountRingBonds(*args)
    def SmallestBondAngle(*args): return _openbabel.OBAtom_SmallestBondAngle(*args)
    def AverageBondAngle(*args): return _openbabel.OBAtom_AverageBondAngle(*args)
    def BOSum(*args): return _openbabel.OBAtom_BOSum(*args)
    def KBOSum(*args): return _openbabel.OBAtom_KBOSum(*args)
    def HtoMethyl(*args): return _openbabel.OBAtom_HtoMethyl(*args)
    def SetHybAndGeom(*args): return _openbabel.OBAtom_SetHybAndGeom(*args)
    def ForceNoH(*args): return _openbabel.OBAtom_ForceNoH(*args)
    def HasNoHForced(*args): return _openbabel.OBAtom_HasNoHForced(*args)
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
class OBBond(OBBase):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    Visit = _swig_property(_openbabel.OBBond_Visit_get, _openbabel.OBBond_Visit_set)
    def __init__(self, *args): 
        _openbabel.OBBond_swiginit(self,_openbabel.new_OBBond(*args))
    __swig_destroy__ = _openbabel.delete_OBBond
    __del__ = lambda self : None;
    def SetIdx(*args): return _openbabel.OBBond_SetIdx(*args)
    def SetBO(*args): return _openbabel.OBBond_SetBO(*args)
    def SetBegin(*args): return _openbabel.OBBond_SetBegin(*args)
    def SetEnd(*args): return _openbabel.OBBond_SetEnd(*args)
    def SetParent(*args): return _openbabel.OBBond_SetParent(*args)
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
    def UnsetHash(*args): return _openbabel.OBBond_UnsetHash(*args)
    def UnsetWedge(*args): return _openbabel.OBBond_UnsetWedge(*args)
    def UnsetUp(*args): return _openbabel.OBBond_UnsetUp(*args)
    def UnsetDown(*args): return _openbabel.OBBond_UnsetDown(*args)
    def UnsetAromatic(*args): return _openbabel.OBBond_UnsetAromatic(*args)
    def UnsetKekule(*args): return _openbabel.OBBond_UnsetKekule(*args)
    def GetIdx(*args): return _openbabel.OBBond_GetIdx(*args)
    def GetBO(*args): return _openbabel.OBBond_GetBO(*args)
    def GetBondOrder(*args): return _openbabel.OBBond_GetBondOrder(*args)
    def GetFlags(*args): return _openbabel.OBBond_GetFlags(*args)
    def GetBeginAtomIdx(*args): return _openbabel.OBBond_GetBeginAtomIdx(*args)
    def GetEndAtomIdx(*args): return _openbabel.OBBond_GetEndAtomIdx(*args)
    def GetBeginAtom(*args): return _openbabel.OBBond_GetBeginAtom(*args)
    def GetEndAtom(*args): return _openbabel.OBBond_GetEndAtom(*args)
    def GetNbrAtom(*args): return _openbabel.OBBond_GetNbrAtom(*args)
    def GetParent(*args): return _openbabel.OBBond_GetParent(*args)
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
class OBMol(OBBase):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBMol_swiginit(self,_openbabel.new_OBMol(*args))
    __swig_destroy__ = _openbabel.delete_OBMol
    __del__ = lambda self : None;
    def __iadd__(*args): return _openbabel.OBMol___iadd__(*args)
    def ReserveAtoms(*args): return _openbabel.OBMol_ReserveAtoms(*args)
    def CreateAtom(*args): return _openbabel.OBMol_CreateAtom(*args)
    def CreateBond(*args): return _openbabel.OBMol_CreateBond(*args)
    def CreateResidue(*args): return _openbabel.OBMol_CreateResidue(*args)
    def DestroyAtom(*args): return _openbabel.OBMol_DestroyAtom(*args)
    def DestroyBond(*args): return _openbabel.OBMol_DestroyBond(*args)
    def DestroyResidue(*args): return _openbabel.OBMol_DestroyResidue(*args)
    def AddAtom(*args): return _openbabel.OBMol_AddAtom(*args)
    def AddBond(*args): return _openbabel.OBMol_AddBond(*args)
    def AddResidue(*args): return _openbabel.OBMol_AddResidue(*args)
    def InsertAtom(*args): return _openbabel.OBMol_InsertAtom(*args)
    def DeleteAtom(*args): return _openbabel.OBMol_DeleteAtom(*args)
    def DeleteBond(*args): return _openbabel.OBMol_DeleteBond(*args)
    def DeleteResidue(*args): return _openbabel.OBMol_DeleteResidue(*args)
    def NewAtom(*args): return _openbabel.OBMol_NewAtom(*args)
    def NewBond(*args): return _openbabel.OBMol_NewBond(*args)
    def NewResidue(*args): return _openbabel.OBMol_NewResidue(*args)
    def BeginModify(*args): return _openbabel.OBMol_BeginModify(*args)
    def EndModify(*args): return _openbabel.OBMol_EndModify(*args)
    def GetMod(*args): return _openbabel.OBMol_GetMod(*args)
    def IncrementMod(*args): return _openbabel.OBMol_IncrementMod(*args)
    def DecrementMod(*args): return _openbabel.OBMol_DecrementMod(*args)
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
    def GetAngle(*args): return _openbabel.OBMol_GetAngle(*args)
    def GetFormula(*args): return _openbabel.OBMol_GetFormula(*args)
    def GetSpacedFormula(*args): return _openbabel.OBMol_GetSpacedFormula(*args)
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
    def BeginAtoms(*args): return _openbabel.OBMol_BeginAtoms(*args)
    def EndAtoms(*args): return _openbabel.OBMol_EndAtoms(*args)
    def BeginBonds(*args): return _openbabel.OBMol_BeginBonds(*args)
    def EndBonds(*args): return _openbabel.OBMol_EndBonds(*args)
    def BeginResidues(*args): return _openbabel.OBMol_BeginResidues(*args)
    def EndResidues(*args): return _openbabel.OBMol_EndResidues(*args)
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

CartesianToInternal = _openbabel.CartesianToInternal
InternalToCartesian = _openbabel.InternalToCartesian
NewExtension = _openbabel.NewExtension
BUFF_SIZE = _openbabel.BUFF_SIZE
get_rmat = _openbabel.get_rmat
ob_make_rmat = _openbabel.ob_make_rmat
qtrfit = _openbabel.qtrfit
superimpose = _openbabel.superimpose
class OBRTree(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    _path = _swig_property(_openbabel.OBRing__path_get, _openbabel.OBRing__path_set)
    _pathset = _swig_property(_openbabel.OBRing__pathset_get, _openbabel.OBRing__pathset_set)
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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

class OBSmartsPattern(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
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
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBSSMatch_swiginit(self,_openbabel.new_OBSSMatch(*args))
    __swig_destroy__ = _openbabel.delete_OBSSMatch
    __del__ = lambda self : None;
    def Match(*args): return _openbabel.OBSSMatch_Match(*args)
OBSSMatch_swigregister = _openbabel.OBSSMatch_swigregister
OBSSMatch_swigregister(OBSSMatch)

SmartsLexReplace = _openbabel.SmartsLexReplace
class OBMolAtomIter(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBMolAtomIter_swiginit(self,_openbabel.new_OBMolAtomIter(*args))
    def good(*args): return _openbabel.OBMolAtomIter_good(*args)
    def inc(*args): return _openbabel.OBMolAtomIter_inc(*args)
    def deref(*args): return _openbabel.OBMolAtomIter_deref(*args)
    def __ref__(*args): return _openbabel.OBMolAtomIter___ref__(*args)
    __swig_destroy__ = _openbabel.delete_OBMolAtomIter
    __del__ = lambda self : None;
    Visit = _swig_property(_openbabel.OBMolAtomIter_Visit_get, _openbabel.OBMolAtomIter_Visit_set)
    def Clear(*args): return _openbabel.OBMolAtomIter_Clear(*args)
    def SetIdx(*args): return _openbabel.OBMolAtomIter_SetIdx(*args)
    def SetHyb(*args): return _openbabel.OBMolAtomIter_SetHyb(*args)
    def SetAtomicNum(*args): return _openbabel.OBMolAtomIter_SetAtomicNum(*args)
    def SetIsotope(*args): return _openbabel.OBMolAtomIter_SetIsotope(*args)
    def SetImplicitValence(*args): return _openbabel.OBMolAtomIter_SetImplicitValence(*args)
    def IncrementImplicitValence(*args): return _openbabel.OBMolAtomIter_IncrementImplicitValence(*args)
    def DecrementImplicitValence(*args): return _openbabel.OBMolAtomIter_DecrementImplicitValence(*args)
    def SetFormalCharge(*args): return _openbabel.OBMolAtomIter_SetFormalCharge(*args)
    def SetSpinMultiplicity(*args): return _openbabel.OBMolAtomIter_SetSpinMultiplicity(*args)
    def SetType(*args): return _openbabel.OBMolAtomIter_SetType(*args)
    def SetPartialCharge(*args): return _openbabel.OBMolAtomIter_SetPartialCharge(*args)
    def SetVector(*args): return _openbabel.OBMolAtomIter_SetVector(*args)
    def SetCoordPtr(*args): return _openbabel.OBMolAtomIter_SetCoordPtr(*args)
    def SetResidue(*args): return _openbabel.OBMolAtomIter_SetResidue(*args)
    def SetParent(*args): return _openbabel.OBMolAtomIter_SetParent(*args)
    def SetAromatic(*args): return _openbabel.OBMolAtomIter_SetAromatic(*args)
    def UnsetAromatic(*args): return _openbabel.OBMolAtomIter_UnsetAromatic(*args)
    def SetClockwiseStereo(*args): return _openbabel.OBMolAtomIter_SetClockwiseStereo(*args)
    def SetAntiClockwiseStereo(*args): return _openbabel.OBMolAtomIter_SetAntiClockwiseStereo(*args)
    def SetPositiveStereo(*args): return _openbabel.OBMolAtomIter_SetPositiveStereo(*args)
    def SetNegativeStereo(*args): return _openbabel.OBMolAtomIter_SetNegativeStereo(*args)
    def UnsetStereo(*args): return _openbabel.OBMolAtomIter_UnsetStereo(*args)
    def SetInRing(*args): return _openbabel.OBMolAtomIter_SetInRing(*args)
    def SetChiral(*args): return _openbabel.OBMolAtomIter_SetChiral(*args)
    def ClearCoordPtr(*args): return _openbabel.OBMolAtomIter_ClearCoordPtr(*args)
    def GetFormalCharge(*args): return _openbabel.OBMolAtomIter_GetFormalCharge(*args)
    def GetAtomicNum(*args): return _openbabel.OBMolAtomIter_GetAtomicNum(*args)
    def GetIsotope(*args): return _openbabel.OBMolAtomIter_GetIsotope(*args)
    def GetSpinMultiplicity(*args): return _openbabel.OBMolAtomIter_GetSpinMultiplicity(*args)
    def GetAtomicMass(*args): return _openbabel.OBMolAtomIter_GetAtomicMass(*args)
    def GetExactMass(*args): return _openbabel.OBMolAtomIter_GetExactMass(*args)
    def GetIdx(*args): return _openbabel.OBMolAtomIter_GetIdx(*args)
    def GetCoordinateIdx(*args): return _openbabel.OBMolAtomIter_GetCoordinateIdx(*args)
    def GetCIdx(*args): return _openbabel.OBMolAtomIter_GetCIdx(*args)
    def GetValence(*args): return _openbabel.OBMolAtomIter_GetValence(*args)
    def GetHyb(*args): return _openbabel.OBMolAtomIter_GetHyb(*args)
    def GetImplicitValence(*args): return _openbabel.OBMolAtomIter_GetImplicitValence(*args)
    def GetHvyValence(*args): return _openbabel.OBMolAtomIter_GetHvyValence(*args)
    def GetHeteroValence(*args): return _openbabel.OBMolAtomIter_GetHeteroValence(*args)
    def GetType(*args): return _openbabel.OBMolAtomIter_GetType(*args)
    def GetX(*args): return _openbabel.OBMolAtomIter_GetX(*args)
    def x(*args): return _openbabel.OBMolAtomIter_x(*args)
    def GetY(*args): return _openbabel.OBMolAtomIter_GetY(*args)
    def y(*args): return _openbabel.OBMolAtomIter_y(*args)
    def GetZ(*args): return _openbabel.OBMolAtomIter_GetZ(*args)
    def z(*args): return _openbabel.OBMolAtomIter_z(*args)
    def GetCoordinate(*args): return _openbabel.OBMolAtomIter_GetCoordinate(*args)
    def GetVector(*args): return _openbabel.OBMolAtomIter_GetVector(*args)
    def GetPartialCharge(*args): return _openbabel.OBMolAtomIter_GetPartialCharge(*args)
    def GetResidue(*args): return _openbabel.OBMolAtomIter_GetResidue(*args)
    def GetParent(*args): return _openbabel.OBMolAtomIter_GetParent(*args)
    def GetNewBondVector(*args): return _openbabel.OBMolAtomIter_GetNewBondVector(*args)
    def GetBond(*args): return _openbabel.OBMolAtomIter_GetBond(*args)
    def GetNextAtom(*args): return _openbabel.OBMolAtomIter_GetNextAtom(*args)
    def BeginBonds(*args): return _openbabel.OBMolAtomIter_BeginBonds(*args)
    def EndBonds(*args): return _openbabel.OBMolAtomIter_EndBonds(*args)
    def BeginBond(*args): return _openbabel.OBMolAtomIter_BeginBond(*args)
    def NextBond(*args): return _openbabel.OBMolAtomIter_NextBond(*args)
    def BeginNbrAtom(*args): return _openbabel.OBMolAtomIter_BeginNbrAtom(*args)
    def NextNbrAtom(*args): return _openbabel.OBMolAtomIter_NextNbrAtom(*args)
    def GetDistance(*args): return _openbabel.OBMolAtomIter_GetDistance(*args)
    def GetAngle(*args): return _openbabel.OBMolAtomIter_GetAngle(*args)
    def NewResidue(*args): return _openbabel.OBMolAtomIter_NewResidue(*args)
    def DeleteResidue(*args): return _openbabel.OBMolAtomIter_DeleteResidue(*args)
    def AddBond(*args): return _openbabel.OBMolAtomIter_AddBond(*args)
    def InsertBond(*args): return _openbabel.OBMolAtomIter_InsertBond(*args)
    def DeleteBond(*args): return _openbabel.OBMolAtomIter_DeleteBond(*args)
    def ClearBond(*args): return _openbabel.OBMolAtomIter_ClearBond(*args)
    def CountFreeOxygens(*args): return _openbabel.OBMolAtomIter_CountFreeOxygens(*args)
    def ImplicitHydrogenCount(*args): return _openbabel.OBMolAtomIter_ImplicitHydrogenCount(*args)
    def ExplicitHydrogenCount(*args): return _openbabel.OBMolAtomIter_ExplicitHydrogenCount(*args)
    def MemberOfRingCount(*args): return _openbabel.OBMolAtomIter_MemberOfRingCount(*args)
    def MemberOfRingSize(*args): return _openbabel.OBMolAtomIter_MemberOfRingSize(*args)
    def CountRingBonds(*args): return _openbabel.OBMolAtomIter_CountRingBonds(*args)
    def SmallestBondAngle(*args): return _openbabel.OBMolAtomIter_SmallestBondAngle(*args)
    def AverageBondAngle(*args): return _openbabel.OBMolAtomIter_AverageBondAngle(*args)
    def BOSum(*args): return _openbabel.OBMolAtomIter_BOSum(*args)
    def KBOSum(*args): return _openbabel.OBMolAtomIter_KBOSum(*args)
    def HtoMethyl(*args): return _openbabel.OBMolAtomIter_HtoMethyl(*args)
    def SetHybAndGeom(*args): return _openbabel.OBMolAtomIter_SetHybAndGeom(*args)
    def ForceNoH(*args): return _openbabel.OBMolAtomIter_ForceNoH(*args)
    def HasNoHForced(*args): return _openbabel.OBMolAtomIter_HasNoHForced(*args)
    def HasResidue(*args): return _openbabel.OBMolAtomIter_HasResidue(*args)
    def IsHydrogen(*args): return _openbabel.OBMolAtomIter_IsHydrogen(*args)
    def IsCarbon(*args): return _openbabel.OBMolAtomIter_IsCarbon(*args)
    def IsNitrogen(*args): return _openbabel.OBMolAtomIter_IsNitrogen(*args)
    def IsOxygen(*args): return _openbabel.OBMolAtomIter_IsOxygen(*args)
    def IsSulfur(*args): return _openbabel.OBMolAtomIter_IsSulfur(*args)
    def IsPhosphorus(*args): return _openbabel.OBMolAtomIter_IsPhosphorus(*args)
    def IsAromatic(*args): return _openbabel.OBMolAtomIter_IsAromatic(*args)
    def IsInRing(*args): return _openbabel.OBMolAtomIter_IsInRing(*args)
    def IsInRingSize(*args): return _openbabel.OBMolAtomIter_IsInRingSize(*args)
    def IsHeteroatom(*args): return _openbabel.OBMolAtomIter_IsHeteroatom(*args)
    def IsNotCorH(*args): return _openbabel.OBMolAtomIter_IsNotCorH(*args)
    def IsConnected(*args): return _openbabel.OBMolAtomIter_IsConnected(*args)
    def IsOneThree(*args): return _openbabel.OBMolAtomIter_IsOneThree(*args)
    def IsOneFour(*args): return _openbabel.OBMolAtomIter_IsOneFour(*args)
    def IsCarboxylOxygen(*args): return _openbabel.OBMolAtomIter_IsCarboxylOxygen(*args)
    def IsPhosphateOxygen(*args): return _openbabel.OBMolAtomIter_IsPhosphateOxygen(*args)
    def IsSulfateOxygen(*args): return _openbabel.OBMolAtomIter_IsSulfateOxygen(*args)
    def IsNitroOxygen(*args): return _openbabel.OBMolAtomIter_IsNitroOxygen(*args)
    def IsAmideNitrogen(*args): return _openbabel.OBMolAtomIter_IsAmideNitrogen(*args)
    def IsPolarHydrogen(*args): return _openbabel.OBMolAtomIter_IsPolarHydrogen(*args)
    def IsNonPolarHydrogen(*args): return _openbabel.OBMolAtomIter_IsNonPolarHydrogen(*args)
    def IsAromaticNOxide(*args): return _openbabel.OBMolAtomIter_IsAromaticNOxide(*args)
    def IsChiral(*args): return _openbabel.OBMolAtomIter_IsChiral(*args)
    def IsAxial(*args): return _openbabel.OBMolAtomIter_IsAxial(*args)
    def IsClockwise(*args): return _openbabel.OBMolAtomIter_IsClockwise(*args)
    def IsAntiClockwise(*args): return _openbabel.OBMolAtomIter_IsAntiClockwise(*args)
    def IsPositiveStereo(*args): return _openbabel.OBMolAtomIter_IsPositiveStereo(*args)
    def IsNegativeStereo(*args): return _openbabel.OBMolAtomIter_IsNegativeStereo(*args)
    def HasChiralitySpecified(*args): return _openbabel.OBMolAtomIter_HasChiralitySpecified(*args)
    def HasChiralVolume(*args): return _openbabel.OBMolAtomIter_HasChiralVolume(*args)
    def IsHbondAcceptor(*args): return _openbabel.OBMolAtomIter_IsHbondAcceptor(*args)
    def IsHbondDonor(*args): return _openbabel.OBMolAtomIter_IsHbondDonor(*args)
    def IsHbondDonorH(*args): return _openbabel.OBMolAtomIter_IsHbondDonorH(*args)
    def HasAlphaBetaUnsat(*args): return _openbabel.OBMolAtomIter_HasAlphaBetaUnsat(*args)
    def HasBondOfOrder(*args): return _openbabel.OBMolAtomIter_HasBondOfOrder(*args)
    def CountBondsOfOrder(*args): return _openbabel.OBMolAtomIter_CountBondsOfOrder(*args)
    def HasNonSingleBond(*args): return _openbabel.OBMolAtomIter_HasNonSingleBond(*args)
    def HasSingleBond(*args): return _openbabel.OBMolAtomIter_HasSingleBond(*args)
    def HasDoubleBond(*args): return _openbabel.OBMolAtomIter_HasDoubleBond(*args)
    def HasAromaticBond(*args): return _openbabel.OBMolAtomIter_HasAromaticBond(*args)
    def MatchesSMARTS(*args): return _openbabel.OBMolAtomIter_MatchesSMARTS(*args)
    def DoTransformations(*args): return _openbabel.OBMolAtomIter_DoTransformations(*args)
    def ClassDescription(*args): return _openbabel.OBMolAtomIter_ClassDescription(*args)
    def HasData(*args): return _openbabel.OBMolAtomIter_HasData(*args)
    def DeleteData(*args): return _openbabel.OBMolAtomIter_DeleteData(*args)
    def SetData(*args): return _openbabel.OBMolAtomIter_SetData(*args)
    def DataSize(*args): return _openbabel.OBMolAtomIter_DataSize(*args)
    def GetData(*args): return _openbabel.OBMolAtomIter_GetData(*args)
    def BeginData(*args): return _openbabel.OBMolAtomIter_BeginData(*args)
    def EndData(*args): return _openbabel.OBMolAtomIter_EndData(*args)
OBMolAtomIter_swigregister = _openbabel.OBMolAtomIter_swigregister
OBMolAtomIter_swigregister(OBMolAtomIter)

class OBMolAtomDFSIter(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBMolAtomDFSIter_swiginit(self,_openbabel.new_OBMolAtomDFSIter(*args))
    def __deref__(*args): return _openbabel.OBMolAtomDFSIter___deref__(*args)
    def __ref__(*args): return _openbabel.OBMolAtomDFSIter___ref__(*args)
    __swig_destroy__ = _openbabel.delete_OBMolAtomDFSIter
    __del__ = lambda self : None;
    Visit = _swig_property(_openbabel.OBMolAtomDFSIter_Visit_get, _openbabel.OBMolAtomDFSIter_Visit_set)
    def Clear(*args): return _openbabel.OBMolAtomDFSIter_Clear(*args)
    def SetIdx(*args): return _openbabel.OBMolAtomDFSIter_SetIdx(*args)
    def SetHyb(*args): return _openbabel.OBMolAtomDFSIter_SetHyb(*args)
    def SetAtomicNum(*args): return _openbabel.OBMolAtomDFSIter_SetAtomicNum(*args)
    def SetIsotope(*args): return _openbabel.OBMolAtomDFSIter_SetIsotope(*args)
    def SetImplicitValence(*args): return _openbabel.OBMolAtomDFSIter_SetImplicitValence(*args)
    def IncrementImplicitValence(*args): return _openbabel.OBMolAtomDFSIter_IncrementImplicitValence(*args)
    def DecrementImplicitValence(*args): return _openbabel.OBMolAtomDFSIter_DecrementImplicitValence(*args)
    def SetFormalCharge(*args): return _openbabel.OBMolAtomDFSIter_SetFormalCharge(*args)
    def SetSpinMultiplicity(*args): return _openbabel.OBMolAtomDFSIter_SetSpinMultiplicity(*args)
    def SetType(*args): return _openbabel.OBMolAtomDFSIter_SetType(*args)
    def SetPartialCharge(*args): return _openbabel.OBMolAtomDFSIter_SetPartialCharge(*args)
    def SetVector(*args): return _openbabel.OBMolAtomDFSIter_SetVector(*args)
    def SetCoordPtr(*args): return _openbabel.OBMolAtomDFSIter_SetCoordPtr(*args)
    def SetResidue(*args): return _openbabel.OBMolAtomDFSIter_SetResidue(*args)
    def SetParent(*args): return _openbabel.OBMolAtomDFSIter_SetParent(*args)
    def SetAromatic(*args): return _openbabel.OBMolAtomDFSIter_SetAromatic(*args)
    def UnsetAromatic(*args): return _openbabel.OBMolAtomDFSIter_UnsetAromatic(*args)
    def SetClockwiseStereo(*args): return _openbabel.OBMolAtomDFSIter_SetClockwiseStereo(*args)
    def SetAntiClockwiseStereo(*args): return _openbabel.OBMolAtomDFSIter_SetAntiClockwiseStereo(*args)
    def SetPositiveStereo(*args): return _openbabel.OBMolAtomDFSIter_SetPositiveStereo(*args)
    def SetNegativeStereo(*args): return _openbabel.OBMolAtomDFSIter_SetNegativeStereo(*args)
    def UnsetStereo(*args): return _openbabel.OBMolAtomDFSIter_UnsetStereo(*args)
    def SetInRing(*args): return _openbabel.OBMolAtomDFSIter_SetInRing(*args)
    def SetChiral(*args): return _openbabel.OBMolAtomDFSIter_SetChiral(*args)
    def ClearCoordPtr(*args): return _openbabel.OBMolAtomDFSIter_ClearCoordPtr(*args)
    def GetFormalCharge(*args): return _openbabel.OBMolAtomDFSIter_GetFormalCharge(*args)
    def GetAtomicNum(*args): return _openbabel.OBMolAtomDFSIter_GetAtomicNum(*args)
    def GetIsotope(*args): return _openbabel.OBMolAtomDFSIter_GetIsotope(*args)
    def GetSpinMultiplicity(*args): return _openbabel.OBMolAtomDFSIter_GetSpinMultiplicity(*args)
    def GetAtomicMass(*args): return _openbabel.OBMolAtomDFSIter_GetAtomicMass(*args)
    def GetExactMass(*args): return _openbabel.OBMolAtomDFSIter_GetExactMass(*args)
    def GetIdx(*args): return _openbabel.OBMolAtomDFSIter_GetIdx(*args)
    def GetCoordinateIdx(*args): return _openbabel.OBMolAtomDFSIter_GetCoordinateIdx(*args)
    def GetCIdx(*args): return _openbabel.OBMolAtomDFSIter_GetCIdx(*args)
    def GetValence(*args): return _openbabel.OBMolAtomDFSIter_GetValence(*args)
    def GetHyb(*args): return _openbabel.OBMolAtomDFSIter_GetHyb(*args)
    def GetImplicitValence(*args): return _openbabel.OBMolAtomDFSIter_GetImplicitValence(*args)
    def GetHvyValence(*args): return _openbabel.OBMolAtomDFSIter_GetHvyValence(*args)
    def GetHeteroValence(*args): return _openbabel.OBMolAtomDFSIter_GetHeteroValence(*args)
    def GetType(*args): return _openbabel.OBMolAtomDFSIter_GetType(*args)
    def GetX(*args): return _openbabel.OBMolAtomDFSIter_GetX(*args)
    def x(*args): return _openbabel.OBMolAtomDFSIter_x(*args)
    def GetY(*args): return _openbabel.OBMolAtomDFSIter_GetY(*args)
    def y(*args): return _openbabel.OBMolAtomDFSIter_y(*args)
    def GetZ(*args): return _openbabel.OBMolAtomDFSIter_GetZ(*args)
    def z(*args): return _openbabel.OBMolAtomDFSIter_z(*args)
    def GetCoordinate(*args): return _openbabel.OBMolAtomDFSIter_GetCoordinate(*args)
    def GetVector(*args): return _openbabel.OBMolAtomDFSIter_GetVector(*args)
    def GetPartialCharge(*args): return _openbabel.OBMolAtomDFSIter_GetPartialCharge(*args)
    def GetResidue(*args): return _openbabel.OBMolAtomDFSIter_GetResidue(*args)
    def GetParent(*args): return _openbabel.OBMolAtomDFSIter_GetParent(*args)
    def GetNewBondVector(*args): return _openbabel.OBMolAtomDFSIter_GetNewBondVector(*args)
    def GetBond(*args): return _openbabel.OBMolAtomDFSIter_GetBond(*args)
    def GetNextAtom(*args): return _openbabel.OBMolAtomDFSIter_GetNextAtom(*args)
    def BeginBonds(*args): return _openbabel.OBMolAtomDFSIter_BeginBonds(*args)
    def EndBonds(*args): return _openbabel.OBMolAtomDFSIter_EndBonds(*args)
    def BeginBond(*args): return _openbabel.OBMolAtomDFSIter_BeginBond(*args)
    def NextBond(*args): return _openbabel.OBMolAtomDFSIter_NextBond(*args)
    def BeginNbrAtom(*args): return _openbabel.OBMolAtomDFSIter_BeginNbrAtom(*args)
    def NextNbrAtom(*args): return _openbabel.OBMolAtomDFSIter_NextNbrAtom(*args)
    def GetDistance(*args): return _openbabel.OBMolAtomDFSIter_GetDistance(*args)
    def GetAngle(*args): return _openbabel.OBMolAtomDFSIter_GetAngle(*args)
    def NewResidue(*args): return _openbabel.OBMolAtomDFSIter_NewResidue(*args)
    def DeleteResidue(*args): return _openbabel.OBMolAtomDFSIter_DeleteResidue(*args)
    def AddBond(*args): return _openbabel.OBMolAtomDFSIter_AddBond(*args)
    def InsertBond(*args): return _openbabel.OBMolAtomDFSIter_InsertBond(*args)
    def DeleteBond(*args): return _openbabel.OBMolAtomDFSIter_DeleteBond(*args)
    def ClearBond(*args): return _openbabel.OBMolAtomDFSIter_ClearBond(*args)
    def CountFreeOxygens(*args): return _openbabel.OBMolAtomDFSIter_CountFreeOxygens(*args)
    def ImplicitHydrogenCount(*args): return _openbabel.OBMolAtomDFSIter_ImplicitHydrogenCount(*args)
    def ExplicitHydrogenCount(*args): return _openbabel.OBMolAtomDFSIter_ExplicitHydrogenCount(*args)
    def MemberOfRingCount(*args): return _openbabel.OBMolAtomDFSIter_MemberOfRingCount(*args)
    def MemberOfRingSize(*args): return _openbabel.OBMolAtomDFSIter_MemberOfRingSize(*args)
    def CountRingBonds(*args): return _openbabel.OBMolAtomDFSIter_CountRingBonds(*args)
    def SmallestBondAngle(*args): return _openbabel.OBMolAtomDFSIter_SmallestBondAngle(*args)
    def AverageBondAngle(*args): return _openbabel.OBMolAtomDFSIter_AverageBondAngle(*args)
    def BOSum(*args): return _openbabel.OBMolAtomDFSIter_BOSum(*args)
    def KBOSum(*args): return _openbabel.OBMolAtomDFSIter_KBOSum(*args)
    def HtoMethyl(*args): return _openbabel.OBMolAtomDFSIter_HtoMethyl(*args)
    def SetHybAndGeom(*args): return _openbabel.OBMolAtomDFSIter_SetHybAndGeom(*args)
    def ForceNoH(*args): return _openbabel.OBMolAtomDFSIter_ForceNoH(*args)
    def HasNoHForced(*args): return _openbabel.OBMolAtomDFSIter_HasNoHForced(*args)
    def HasResidue(*args): return _openbabel.OBMolAtomDFSIter_HasResidue(*args)
    def IsHydrogen(*args): return _openbabel.OBMolAtomDFSIter_IsHydrogen(*args)
    def IsCarbon(*args): return _openbabel.OBMolAtomDFSIter_IsCarbon(*args)
    def IsNitrogen(*args): return _openbabel.OBMolAtomDFSIter_IsNitrogen(*args)
    def IsOxygen(*args): return _openbabel.OBMolAtomDFSIter_IsOxygen(*args)
    def IsSulfur(*args): return _openbabel.OBMolAtomDFSIter_IsSulfur(*args)
    def IsPhosphorus(*args): return _openbabel.OBMolAtomDFSIter_IsPhosphorus(*args)
    def IsAromatic(*args): return _openbabel.OBMolAtomDFSIter_IsAromatic(*args)
    def IsInRing(*args): return _openbabel.OBMolAtomDFSIter_IsInRing(*args)
    def IsInRingSize(*args): return _openbabel.OBMolAtomDFSIter_IsInRingSize(*args)
    def IsHeteroatom(*args): return _openbabel.OBMolAtomDFSIter_IsHeteroatom(*args)
    def IsNotCorH(*args): return _openbabel.OBMolAtomDFSIter_IsNotCorH(*args)
    def IsConnected(*args): return _openbabel.OBMolAtomDFSIter_IsConnected(*args)
    def IsOneThree(*args): return _openbabel.OBMolAtomDFSIter_IsOneThree(*args)
    def IsOneFour(*args): return _openbabel.OBMolAtomDFSIter_IsOneFour(*args)
    def IsCarboxylOxygen(*args): return _openbabel.OBMolAtomDFSIter_IsCarboxylOxygen(*args)
    def IsPhosphateOxygen(*args): return _openbabel.OBMolAtomDFSIter_IsPhosphateOxygen(*args)
    def IsSulfateOxygen(*args): return _openbabel.OBMolAtomDFSIter_IsSulfateOxygen(*args)
    def IsNitroOxygen(*args): return _openbabel.OBMolAtomDFSIter_IsNitroOxygen(*args)
    def IsAmideNitrogen(*args): return _openbabel.OBMolAtomDFSIter_IsAmideNitrogen(*args)
    def IsPolarHydrogen(*args): return _openbabel.OBMolAtomDFSIter_IsPolarHydrogen(*args)
    def IsNonPolarHydrogen(*args): return _openbabel.OBMolAtomDFSIter_IsNonPolarHydrogen(*args)
    def IsAromaticNOxide(*args): return _openbabel.OBMolAtomDFSIter_IsAromaticNOxide(*args)
    def IsChiral(*args): return _openbabel.OBMolAtomDFSIter_IsChiral(*args)
    def IsAxial(*args): return _openbabel.OBMolAtomDFSIter_IsAxial(*args)
    def IsClockwise(*args): return _openbabel.OBMolAtomDFSIter_IsClockwise(*args)
    def IsAntiClockwise(*args): return _openbabel.OBMolAtomDFSIter_IsAntiClockwise(*args)
    def IsPositiveStereo(*args): return _openbabel.OBMolAtomDFSIter_IsPositiveStereo(*args)
    def IsNegativeStereo(*args): return _openbabel.OBMolAtomDFSIter_IsNegativeStereo(*args)
    def HasChiralitySpecified(*args): return _openbabel.OBMolAtomDFSIter_HasChiralitySpecified(*args)
    def HasChiralVolume(*args): return _openbabel.OBMolAtomDFSIter_HasChiralVolume(*args)
    def IsHbondAcceptor(*args): return _openbabel.OBMolAtomDFSIter_IsHbondAcceptor(*args)
    def IsHbondDonor(*args): return _openbabel.OBMolAtomDFSIter_IsHbondDonor(*args)
    def IsHbondDonorH(*args): return _openbabel.OBMolAtomDFSIter_IsHbondDonorH(*args)
    def HasAlphaBetaUnsat(*args): return _openbabel.OBMolAtomDFSIter_HasAlphaBetaUnsat(*args)
    def HasBondOfOrder(*args): return _openbabel.OBMolAtomDFSIter_HasBondOfOrder(*args)
    def CountBondsOfOrder(*args): return _openbabel.OBMolAtomDFSIter_CountBondsOfOrder(*args)
    def HasNonSingleBond(*args): return _openbabel.OBMolAtomDFSIter_HasNonSingleBond(*args)
    def HasSingleBond(*args): return _openbabel.OBMolAtomDFSIter_HasSingleBond(*args)
    def HasDoubleBond(*args): return _openbabel.OBMolAtomDFSIter_HasDoubleBond(*args)
    def HasAromaticBond(*args): return _openbabel.OBMolAtomDFSIter_HasAromaticBond(*args)
    def MatchesSMARTS(*args): return _openbabel.OBMolAtomDFSIter_MatchesSMARTS(*args)
    def DoTransformations(*args): return _openbabel.OBMolAtomDFSIter_DoTransformations(*args)
    def ClassDescription(*args): return _openbabel.OBMolAtomDFSIter_ClassDescription(*args)
    def HasData(*args): return _openbabel.OBMolAtomDFSIter_HasData(*args)
    def DeleteData(*args): return _openbabel.OBMolAtomDFSIter_DeleteData(*args)
    def SetData(*args): return _openbabel.OBMolAtomDFSIter_SetData(*args)
    def DataSize(*args): return _openbabel.OBMolAtomDFSIter_DataSize(*args)
    def GetData(*args): return _openbabel.OBMolAtomDFSIter_GetData(*args)
    def BeginData(*args): return _openbabel.OBMolAtomDFSIter_BeginData(*args)
    def EndData(*args): return _openbabel.OBMolAtomDFSIter_EndData(*args)
OBMolAtomDFSIter_swigregister = _openbabel.OBMolAtomDFSIter_swigregister
OBMolAtomDFSIter_swigregister(OBMolAtomDFSIter)

class OBMolAtomBFSIter(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBMolAtomBFSIter_swiginit(self,_openbabel.new_OBMolAtomBFSIter(*args))
    def __deref__(*args): return _openbabel.OBMolAtomBFSIter___deref__(*args)
    def __ref__(*args): return _openbabel.OBMolAtomBFSIter___ref__(*args)
    __swig_destroy__ = _openbabel.delete_OBMolAtomBFSIter
    __del__ = lambda self : None;
    Visit = _swig_property(_openbabel.OBMolAtomBFSIter_Visit_get, _openbabel.OBMolAtomBFSIter_Visit_set)
    def Clear(*args): return _openbabel.OBMolAtomBFSIter_Clear(*args)
    def SetIdx(*args): return _openbabel.OBMolAtomBFSIter_SetIdx(*args)
    def SetHyb(*args): return _openbabel.OBMolAtomBFSIter_SetHyb(*args)
    def SetAtomicNum(*args): return _openbabel.OBMolAtomBFSIter_SetAtomicNum(*args)
    def SetIsotope(*args): return _openbabel.OBMolAtomBFSIter_SetIsotope(*args)
    def SetImplicitValence(*args): return _openbabel.OBMolAtomBFSIter_SetImplicitValence(*args)
    def IncrementImplicitValence(*args): return _openbabel.OBMolAtomBFSIter_IncrementImplicitValence(*args)
    def DecrementImplicitValence(*args): return _openbabel.OBMolAtomBFSIter_DecrementImplicitValence(*args)
    def SetFormalCharge(*args): return _openbabel.OBMolAtomBFSIter_SetFormalCharge(*args)
    def SetSpinMultiplicity(*args): return _openbabel.OBMolAtomBFSIter_SetSpinMultiplicity(*args)
    def SetType(*args): return _openbabel.OBMolAtomBFSIter_SetType(*args)
    def SetPartialCharge(*args): return _openbabel.OBMolAtomBFSIter_SetPartialCharge(*args)
    def SetVector(*args): return _openbabel.OBMolAtomBFSIter_SetVector(*args)
    def SetCoordPtr(*args): return _openbabel.OBMolAtomBFSIter_SetCoordPtr(*args)
    def SetResidue(*args): return _openbabel.OBMolAtomBFSIter_SetResidue(*args)
    def SetParent(*args): return _openbabel.OBMolAtomBFSIter_SetParent(*args)
    def SetAromatic(*args): return _openbabel.OBMolAtomBFSIter_SetAromatic(*args)
    def UnsetAromatic(*args): return _openbabel.OBMolAtomBFSIter_UnsetAromatic(*args)
    def SetClockwiseStereo(*args): return _openbabel.OBMolAtomBFSIter_SetClockwiseStereo(*args)
    def SetAntiClockwiseStereo(*args): return _openbabel.OBMolAtomBFSIter_SetAntiClockwiseStereo(*args)
    def SetPositiveStereo(*args): return _openbabel.OBMolAtomBFSIter_SetPositiveStereo(*args)
    def SetNegativeStereo(*args): return _openbabel.OBMolAtomBFSIter_SetNegativeStereo(*args)
    def UnsetStereo(*args): return _openbabel.OBMolAtomBFSIter_UnsetStereo(*args)
    def SetInRing(*args): return _openbabel.OBMolAtomBFSIter_SetInRing(*args)
    def SetChiral(*args): return _openbabel.OBMolAtomBFSIter_SetChiral(*args)
    def ClearCoordPtr(*args): return _openbabel.OBMolAtomBFSIter_ClearCoordPtr(*args)
    def GetFormalCharge(*args): return _openbabel.OBMolAtomBFSIter_GetFormalCharge(*args)
    def GetAtomicNum(*args): return _openbabel.OBMolAtomBFSIter_GetAtomicNum(*args)
    def GetIsotope(*args): return _openbabel.OBMolAtomBFSIter_GetIsotope(*args)
    def GetSpinMultiplicity(*args): return _openbabel.OBMolAtomBFSIter_GetSpinMultiplicity(*args)
    def GetAtomicMass(*args): return _openbabel.OBMolAtomBFSIter_GetAtomicMass(*args)
    def GetExactMass(*args): return _openbabel.OBMolAtomBFSIter_GetExactMass(*args)
    def GetIdx(*args): return _openbabel.OBMolAtomBFSIter_GetIdx(*args)
    def GetCoordinateIdx(*args): return _openbabel.OBMolAtomBFSIter_GetCoordinateIdx(*args)
    def GetCIdx(*args): return _openbabel.OBMolAtomBFSIter_GetCIdx(*args)
    def GetValence(*args): return _openbabel.OBMolAtomBFSIter_GetValence(*args)
    def GetHyb(*args): return _openbabel.OBMolAtomBFSIter_GetHyb(*args)
    def GetImplicitValence(*args): return _openbabel.OBMolAtomBFSIter_GetImplicitValence(*args)
    def GetHvyValence(*args): return _openbabel.OBMolAtomBFSIter_GetHvyValence(*args)
    def GetHeteroValence(*args): return _openbabel.OBMolAtomBFSIter_GetHeteroValence(*args)
    def GetType(*args): return _openbabel.OBMolAtomBFSIter_GetType(*args)
    def GetX(*args): return _openbabel.OBMolAtomBFSIter_GetX(*args)
    def x(*args): return _openbabel.OBMolAtomBFSIter_x(*args)
    def GetY(*args): return _openbabel.OBMolAtomBFSIter_GetY(*args)
    def y(*args): return _openbabel.OBMolAtomBFSIter_y(*args)
    def GetZ(*args): return _openbabel.OBMolAtomBFSIter_GetZ(*args)
    def z(*args): return _openbabel.OBMolAtomBFSIter_z(*args)
    def GetCoordinate(*args): return _openbabel.OBMolAtomBFSIter_GetCoordinate(*args)
    def GetVector(*args): return _openbabel.OBMolAtomBFSIter_GetVector(*args)
    def GetPartialCharge(*args): return _openbabel.OBMolAtomBFSIter_GetPartialCharge(*args)
    def GetResidue(*args): return _openbabel.OBMolAtomBFSIter_GetResidue(*args)
    def GetParent(*args): return _openbabel.OBMolAtomBFSIter_GetParent(*args)
    def GetNewBondVector(*args): return _openbabel.OBMolAtomBFSIter_GetNewBondVector(*args)
    def GetBond(*args): return _openbabel.OBMolAtomBFSIter_GetBond(*args)
    def GetNextAtom(*args): return _openbabel.OBMolAtomBFSIter_GetNextAtom(*args)
    def BeginBonds(*args): return _openbabel.OBMolAtomBFSIter_BeginBonds(*args)
    def EndBonds(*args): return _openbabel.OBMolAtomBFSIter_EndBonds(*args)
    def BeginBond(*args): return _openbabel.OBMolAtomBFSIter_BeginBond(*args)
    def NextBond(*args): return _openbabel.OBMolAtomBFSIter_NextBond(*args)
    def BeginNbrAtom(*args): return _openbabel.OBMolAtomBFSIter_BeginNbrAtom(*args)
    def NextNbrAtom(*args): return _openbabel.OBMolAtomBFSIter_NextNbrAtom(*args)
    def GetDistance(*args): return _openbabel.OBMolAtomBFSIter_GetDistance(*args)
    def GetAngle(*args): return _openbabel.OBMolAtomBFSIter_GetAngle(*args)
    def NewResidue(*args): return _openbabel.OBMolAtomBFSIter_NewResidue(*args)
    def DeleteResidue(*args): return _openbabel.OBMolAtomBFSIter_DeleteResidue(*args)
    def AddBond(*args): return _openbabel.OBMolAtomBFSIter_AddBond(*args)
    def InsertBond(*args): return _openbabel.OBMolAtomBFSIter_InsertBond(*args)
    def DeleteBond(*args): return _openbabel.OBMolAtomBFSIter_DeleteBond(*args)
    def ClearBond(*args): return _openbabel.OBMolAtomBFSIter_ClearBond(*args)
    def CountFreeOxygens(*args): return _openbabel.OBMolAtomBFSIter_CountFreeOxygens(*args)
    def ImplicitHydrogenCount(*args): return _openbabel.OBMolAtomBFSIter_ImplicitHydrogenCount(*args)
    def ExplicitHydrogenCount(*args): return _openbabel.OBMolAtomBFSIter_ExplicitHydrogenCount(*args)
    def MemberOfRingCount(*args): return _openbabel.OBMolAtomBFSIter_MemberOfRingCount(*args)
    def MemberOfRingSize(*args): return _openbabel.OBMolAtomBFSIter_MemberOfRingSize(*args)
    def CountRingBonds(*args): return _openbabel.OBMolAtomBFSIter_CountRingBonds(*args)
    def SmallestBondAngle(*args): return _openbabel.OBMolAtomBFSIter_SmallestBondAngle(*args)
    def AverageBondAngle(*args): return _openbabel.OBMolAtomBFSIter_AverageBondAngle(*args)
    def BOSum(*args): return _openbabel.OBMolAtomBFSIter_BOSum(*args)
    def KBOSum(*args): return _openbabel.OBMolAtomBFSIter_KBOSum(*args)
    def HtoMethyl(*args): return _openbabel.OBMolAtomBFSIter_HtoMethyl(*args)
    def SetHybAndGeom(*args): return _openbabel.OBMolAtomBFSIter_SetHybAndGeom(*args)
    def ForceNoH(*args): return _openbabel.OBMolAtomBFSIter_ForceNoH(*args)
    def HasNoHForced(*args): return _openbabel.OBMolAtomBFSIter_HasNoHForced(*args)
    def HasResidue(*args): return _openbabel.OBMolAtomBFSIter_HasResidue(*args)
    def IsHydrogen(*args): return _openbabel.OBMolAtomBFSIter_IsHydrogen(*args)
    def IsCarbon(*args): return _openbabel.OBMolAtomBFSIter_IsCarbon(*args)
    def IsNitrogen(*args): return _openbabel.OBMolAtomBFSIter_IsNitrogen(*args)
    def IsOxygen(*args): return _openbabel.OBMolAtomBFSIter_IsOxygen(*args)
    def IsSulfur(*args): return _openbabel.OBMolAtomBFSIter_IsSulfur(*args)
    def IsPhosphorus(*args): return _openbabel.OBMolAtomBFSIter_IsPhosphorus(*args)
    def IsAromatic(*args): return _openbabel.OBMolAtomBFSIter_IsAromatic(*args)
    def IsInRing(*args): return _openbabel.OBMolAtomBFSIter_IsInRing(*args)
    def IsInRingSize(*args): return _openbabel.OBMolAtomBFSIter_IsInRingSize(*args)
    def IsHeteroatom(*args): return _openbabel.OBMolAtomBFSIter_IsHeteroatom(*args)
    def IsNotCorH(*args): return _openbabel.OBMolAtomBFSIter_IsNotCorH(*args)
    def IsConnected(*args): return _openbabel.OBMolAtomBFSIter_IsConnected(*args)
    def IsOneThree(*args): return _openbabel.OBMolAtomBFSIter_IsOneThree(*args)
    def IsOneFour(*args): return _openbabel.OBMolAtomBFSIter_IsOneFour(*args)
    def IsCarboxylOxygen(*args): return _openbabel.OBMolAtomBFSIter_IsCarboxylOxygen(*args)
    def IsPhosphateOxygen(*args): return _openbabel.OBMolAtomBFSIter_IsPhosphateOxygen(*args)
    def IsSulfateOxygen(*args): return _openbabel.OBMolAtomBFSIter_IsSulfateOxygen(*args)
    def IsNitroOxygen(*args): return _openbabel.OBMolAtomBFSIter_IsNitroOxygen(*args)
    def IsAmideNitrogen(*args): return _openbabel.OBMolAtomBFSIter_IsAmideNitrogen(*args)
    def IsPolarHydrogen(*args): return _openbabel.OBMolAtomBFSIter_IsPolarHydrogen(*args)
    def IsNonPolarHydrogen(*args): return _openbabel.OBMolAtomBFSIter_IsNonPolarHydrogen(*args)
    def IsAromaticNOxide(*args): return _openbabel.OBMolAtomBFSIter_IsAromaticNOxide(*args)
    def IsChiral(*args): return _openbabel.OBMolAtomBFSIter_IsChiral(*args)
    def IsAxial(*args): return _openbabel.OBMolAtomBFSIter_IsAxial(*args)
    def IsClockwise(*args): return _openbabel.OBMolAtomBFSIter_IsClockwise(*args)
    def IsAntiClockwise(*args): return _openbabel.OBMolAtomBFSIter_IsAntiClockwise(*args)
    def IsPositiveStereo(*args): return _openbabel.OBMolAtomBFSIter_IsPositiveStereo(*args)
    def IsNegativeStereo(*args): return _openbabel.OBMolAtomBFSIter_IsNegativeStereo(*args)
    def HasChiralitySpecified(*args): return _openbabel.OBMolAtomBFSIter_HasChiralitySpecified(*args)
    def HasChiralVolume(*args): return _openbabel.OBMolAtomBFSIter_HasChiralVolume(*args)
    def IsHbondAcceptor(*args): return _openbabel.OBMolAtomBFSIter_IsHbondAcceptor(*args)
    def IsHbondDonor(*args): return _openbabel.OBMolAtomBFSIter_IsHbondDonor(*args)
    def IsHbondDonorH(*args): return _openbabel.OBMolAtomBFSIter_IsHbondDonorH(*args)
    def HasAlphaBetaUnsat(*args): return _openbabel.OBMolAtomBFSIter_HasAlphaBetaUnsat(*args)
    def HasBondOfOrder(*args): return _openbabel.OBMolAtomBFSIter_HasBondOfOrder(*args)
    def CountBondsOfOrder(*args): return _openbabel.OBMolAtomBFSIter_CountBondsOfOrder(*args)
    def HasNonSingleBond(*args): return _openbabel.OBMolAtomBFSIter_HasNonSingleBond(*args)
    def HasSingleBond(*args): return _openbabel.OBMolAtomBFSIter_HasSingleBond(*args)
    def HasDoubleBond(*args): return _openbabel.OBMolAtomBFSIter_HasDoubleBond(*args)
    def HasAromaticBond(*args): return _openbabel.OBMolAtomBFSIter_HasAromaticBond(*args)
    def MatchesSMARTS(*args): return _openbabel.OBMolAtomBFSIter_MatchesSMARTS(*args)
    def DoTransformations(*args): return _openbabel.OBMolAtomBFSIter_DoTransformations(*args)
    def ClassDescription(*args): return _openbabel.OBMolAtomBFSIter_ClassDescription(*args)
    def HasData(*args): return _openbabel.OBMolAtomBFSIter_HasData(*args)
    def DeleteData(*args): return _openbabel.OBMolAtomBFSIter_DeleteData(*args)
    def SetData(*args): return _openbabel.OBMolAtomBFSIter_SetData(*args)
    def DataSize(*args): return _openbabel.OBMolAtomBFSIter_DataSize(*args)
    def GetData(*args): return _openbabel.OBMolAtomBFSIter_GetData(*args)
    def BeginData(*args): return _openbabel.OBMolAtomBFSIter_BeginData(*args)
    def EndData(*args): return _openbabel.OBMolAtomBFSIter_EndData(*args)
OBMolAtomBFSIter_swigregister = _openbabel.OBMolAtomBFSIter_swigregister
OBMolAtomBFSIter_swigregister(OBMolAtomBFSIter)

class OBMolBondIter(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBMolBondIter_swiginit(self,_openbabel.new_OBMolBondIter(*args))
    def good(*args): return _openbabel.OBMolBondIter_good(*args)
    def inc(*args): return _openbabel.OBMolBondIter_inc(*args)
    def deref(*args): return _openbabel.OBMolBondIter_deref(*args)
    def __ref__(*args): return _openbabel.OBMolBondIter___ref__(*args)
    __swig_destroy__ = _openbabel.delete_OBMolBondIter
    __del__ = lambda self : None;
    Visit = _swig_property(_openbabel.OBMolBondIter_Visit_get, _openbabel.OBMolBondIter_Visit_set)
    def SetIdx(*args): return _openbabel.OBMolBondIter_SetIdx(*args)
    def SetBO(*args): return _openbabel.OBMolBondIter_SetBO(*args)
    def SetBegin(*args): return _openbabel.OBMolBondIter_SetBegin(*args)
    def SetEnd(*args): return _openbabel.OBMolBondIter_SetEnd(*args)
    def SetParent(*args): return _openbabel.OBMolBondIter_SetParent(*args)
    def SetLength(*args): return _openbabel.OBMolBondIter_SetLength(*args)
    def Set(*args): return _openbabel.OBMolBondIter_Set(*args)
    def SetKSingle(*args): return _openbabel.OBMolBondIter_SetKSingle(*args)
    def SetKDouble(*args): return _openbabel.OBMolBondIter_SetKDouble(*args)
    def SetKTriple(*args): return _openbabel.OBMolBondIter_SetKTriple(*args)
    def SetAromatic(*args): return _openbabel.OBMolBondIter_SetAromatic(*args)
    def SetHash(*args): return _openbabel.OBMolBondIter_SetHash(*args)
    def SetWedge(*args): return _openbabel.OBMolBondIter_SetWedge(*args)
    def SetUp(*args): return _openbabel.OBMolBondIter_SetUp(*args)
    def SetDown(*args): return _openbabel.OBMolBondIter_SetDown(*args)
    def SetInRing(*args): return _openbabel.OBMolBondIter_SetInRing(*args)
    def SetClosure(*args): return _openbabel.OBMolBondIter_SetClosure(*args)
    def UnsetHash(*args): return _openbabel.OBMolBondIter_UnsetHash(*args)
    def UnsetWedge(*args): return _openbabel.OBMolBondIter_UnsetWedge(*args)
    def UnsetUp(*args): return _openbabel.OBMolBondIter_UnsetUp(*args)
    def UnsetDown(*args): return _openbabel.OBMolBondIter_UnsetDown(*args)
    def UnsetAromatic(*args): return _openbabel.OBMolBondIter_UnsetAromatic(*args)
    def UnsetKekule(*args): return _openbabel.OBMolBondIter_UnsetKekule(*args)
    def GetIdx(*args): return _openbabel.OBMolBondIter_GetIdx(*args)
    def GetBO(*args): return _openbabel.OBMolBondIter_GetBO(*args)
    def GetBondOrder(*args): return _openbabel.OBMolBondIter_GetBondOrder(*args)
    def GetFlags(*args): return _openbabel.OBMolBondIter_GetFlags(*args)
    def GetBeginAtomIdx(*args): return _openbabel.OBMolBondIter_GetBeginAtomIdx(*args)
    def GetEndAtomIdx(*args): return _openbabel.OBMolBondIter_GetEndAtomIdx(*args)
    def GetBeginAtom(*args): return _openbabel.OBMolBondIter_GetBeginAtom(*args)
    def GetEndAtom(*args): return _openbabel.OBMolBondIter_GetEndAtom(*args)
    def GetNbrAtom(*args): return _openbabel.OBMolBondIter_GetNbrAtom(*args)
    def GetParent(*args): return _openbabel.OBMolBondIter_GetParent(*args)
    def GetEquibLength(*args): return _openbabel.OBMolBondIter_GetEquibLength(*args)
    def GetLength(*args): return _openbabel.OBMolBondIter_GetLength(*args)
    def GetNbrAtomIdx(*args): return _openbabel.OBMolBondIter_GetNbrAtomIdx(*args)
    def IsAromatic(*args): return _openbabel.OBMolBondIter_IsAromatic(*args)
    def IsInRing(*args): return _openbabel.OBMolBondIter_IsInRing(*args)
    def IsRotor(*args): return _openbabel.OBMolBondIter_IsRotor(*args)
    def IsAmide(*args): return _openbabel.OBMolBondIter_IsAmide(*args)
    def IsPrimaryAmide(*args): return _openbabel.OBMolBondIter_IsPrimaryAmide(*args)
    def IsSecondaryAmide(*args): return _openbabel.OBMolBondIter_IsSecondaryAmide(*args)
    def IsEster(*args): return _openbabel.OBMolBondIter_IsEster(*args)
    def IsCarbonyl(*args): return _openbabel.OBMolBondIter_IsCarbonyl(*args)
    def IsSingle(*args): return _openbabel.OBMolBondIter_IsSingle(*args)
    def IsDouble(*args): return _openbabel.OBMolBondIter_IsDouble(*args)
    def IsTriple(*args): return _openbabel.OBMolBondIter_IsTriple(*args)
    def IsKSingle(*args): return _openbabel.OBMolBondIter_IsKSingle(*args)
    def IsKDouble(*args): return _openbabel.OBMolBondIter_IsKDouble(*args)
    def IsKTriple(*args): return _openbabel.OBMolBondIter_IsKTriple(*args)
    def IsClosure(*args): return _openbabel.OBMolBondIter_IsClosure(*args)
    def IsUp(*args): return _openbabel.OBMolBondIter_IsUp(*args)
    def IsDown(*args): return _openbabel.OBMolBondIter_IsDown(*args)
    def IsWedge(*args): return _openbabel.OBMolBondIter_IsWedge(*args)
    def IsHash(*args): return _openbabel.OBMolBondIter_IsHash(*args)
    def IsDoubleBondGeometry(*args): return _openbabel.OBMolBondIter_IsDoubleBondGeometry(*args)
    def DoTransformations(*args): return _openbabel.OBMolBondIter_DoTransformations(*args)
    def ClassDescription(*args): return _openbabel.OBMolBondIter_ClassDescription(*args)
    def HasData(*args): return _openbabel.OBMolBondIter_HasData(*args)
    def DeleteData(*args): return _openbabel.OBMolBondIter_DeleteData(*args)
    def SetData(*args): return _openbabel.OBMolBondIter_SetData(*args)
    def DataSize(*args): return _openbabel.OBMolBondIter_DataSize(*args)
    def GetData(*args): return _openbabel.OBMolBondIter_GetData(*args)
    def BeginData(*args): return _openbabel.OBMolBondIter_BeginData(*args)
    def EndData(*args): return _openbabel.OBMolBondIter_EndData(*args)
OBMolBondIter_swigregister = _openbabel.OBMolBondIter_swigregister
OBMolBondIter_swigregister(OBMolBondIter)

class OBAtomAtomIter(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBAtomAtomIter_swiginit(self,_openbabel.new_OBAtomAtomIter(*args))
    def good(*args): return _openbabel.OBAtomAtomIter_good(*args)
    def inc(*args): return _openbabel.OBAtomAtomIter_inc(*args)
    def deref(*args): return _openbabel.OBAtomAtomIter_deref(*args)
    def __ref__(*args): return _openbabel.OBAtomAtomIter___ref__(*args)
    __swig_destroy__ = _openbabel.delete_OBAtomAtomIter
    __del__ = lambda self : None;
    Visit = _swig_property(_openbabel.OBAtomAtomIter_Visit_get, _openbabel.OBAtomAtomIter_Visit_set)
    def Clear(*args): return _openbabel.OBAtomAtomIter_Clear(*args)
    def SetIdx(*args): return _openbabel.OBAtomAtomIter_SetIdx(*args)
    def SetHyb(*args): return _openbabel.OBAtomAtomIter_SetHyb(*args)
    def SetAtomicNum(*args): return _openbabel.OBAtomAtomIter_SetAtomicNum(*args)
    def SetIsotope(*args): return _openbabel.OBAtomAtomIter_SetIsotope(*args)
    def SetImplicitValence(*args): return _openbabel.OBAtomAtomIter_SetImplicitValence(*args)
    def IncrementImplicitValence(*args): return _openbabel.OBAtomAtomIter_IncrementImplicitValence(*args)
    def DecrementImplicitValence(*args): return _openbabel.OBAtomAtomIter_DecrementImplicitValence(*args)
    def SetFormalCharge(*args): return _openbabel.OBAtomAtomIter_SetFormalCharge(*args)
    def SetSpinMultiplicity(*args): return _openbabel.OBAtomAtomIter_SetSpinMultiplicity(*args)
    def SetType(*args): return _openbabel.OBAtomAtomIter_SetType(*args)
    def SetPartialCharge(*args): return _openbabel.OBAtomAtomIter_SetPartialCharge(*args)
    def SetVector(*args): return _openbabel.OBAtomAtomIter_SetVector(*args)
    def SetCoordPtr(*args): return _openbabel.OBAtomAtomIter_SetCoordPtr(*args)
    def SetResidue(*args): return _openbabel.OBAtomAtomIter_SetResidue(*args)
    def SetParent(*args): return _openbabel.OBAtomAtomIter_SetParent(*args)
    def SetAromatic(*args): return _openbabel.OBAtomAtomIter_SetAromatic(*args)
    def UnsetAromatic(*args): return _openbabel.OBAtomAtomIter_UnsetAromatic(*args)
    def SetClockwiseStereo(*args): return _openbabel.OBAtomAtomIter_SetClockwiseStereo(*args)
    def SetAntiClockwiseStereo(*args): return _openbabel.OBAtomAtomIter_SetAntiClockwiseStereo(*args)
    def SetPositiveStereo(*args): return _openbabel.OBAtomAtomIter_SetPositiveStereo(*args)
    def SetNegativeStereo(*args): return _openbabel.OBAtomAtomIter_SetNegativeStereo(*args)
    def UnsetStereo(*args): return _openbabel.OBAtomAtomIter_UnsetStereo(*args)
    def SetInRing(*args): return _openbabel.OBAtomAtomIter_SetInRing(*args)
    def SetChiral(*args): return _openbabel.OBAtomAtomIter_SetChiral(*args)
    def ClearCoordPtr(*args): return _openbabel.OBAtomAtomIter_ClearCoordPtr(*args)
    def GetFormalCharge(*args): return _openbabel.OBAtomAtomIter_GetFormalCharge(*args)
    def GetAtomicNum(*args): return _openbabel.OBAtomAtomIter_GetAtomicNum(*args)
    def GetIsotope(*args): return _openbabel.OBAtomAtomIter_GetIsotope(*args)
    def GetSpinMultiplicity(*args): return _openbabel.OBAtomAtomIter_GetSpinMultiplicity(*args)
    def GetAtomicMass(*args): return _openbabel.OBAtomAtomIter_GetAtomicMass(*args)
    def GetExactMass(*args): return _openbabel.OBAtomAtomIter_GetExactMass(*args)
    def GetIdx(*args): return _openbabel.OBAtomAtomIter_GetIdx(*args)
    def GetCoordinateIdx(*args): return _openbabel.OBAtomAtomIter_GetCoordinateIdx(*args)
    def GetCIdx(*args): return _openbabel.OBAtomAtomIter_GetCIdx(*args)
    def GetValence(*args): return _openbabel.OBAtomAtomIter_GetValence(*args)
    def GetHyb(*args): return _openbabel.OBAtomAtomIter_GetHyb(*args)
    def GetImplicitValence(*args): return _openbabel.OBAtomAtomIter_GetImplicitValence(*args)
    def GetHvyValence(*args): return _openbabel.OBAtomAtomIter_GetHvyValence(*args)
    def GetHeteroValence(*args): return _openbabel.OBAtomAtomIter_GetHeteroValence(*args)
    def GetType(*args): return _openbabel.OBAtomAtomIter_GetType(*args)
    def GetX(*args): return _openbabel.OBAtomAtomIter_GetX(*args)
    def x(*args): return _openbabel.OBAtomAtomIter_x(*args)
    def GetY(*args): return _openbabel.OBAtomAtomIter_GetY(*args)
    def y(*args): return _openbabel.OBAtomAtomIter_y(*args)
    def GetZ(*args): return _openbabel.OBAtomAtomIter_GetZ(*args)
    def z(*args): return _openbabel.OBAtomAtomIter_z(*args)
    def GetCoordinate(*args): return _openbabel.OBAtomAtomIter_GetCoordinate(*args)
    def GetVector(*args): return _openbabel.OBAtomAtomIter_GetVector(*args)
    def GetPartialCharge(*args): return _openbabel.OBAtomAtomIter_GetPartialCharge(*args)
    def GetResidue(*args): return _openbabel.OBAtomAtomIter_GetResidue(*args)
    def GetParent(*args): return _openbabel.OBAtomAtomIter_GetParent(*args)
    def GetNewBondVector(*args): return _openbabel.OBAtomAtomIter_GetNewBondVector(*args)
    def GetBond(*args): return _openbabel.OBAtomAtomIter_GetBond(*args)
    def GetNextAtom(*args): return _openbabel.OBAtomAtomIter_GetNextAtom(*args)
    def BeginBonds(*args): return _openbabel.OBAtomAtomIter_BeginBonds(*args)
    def EndBonds(*args): return _openbabel.OBAtomAtomIter_EndBonds(*args)
    def BeginBond(*args): return _openbabel.OBAtomAtomIter_BeginBond(*args)
    def NextBond(*args): return _openbabel.OBAtomAtomIter_NextBond(*args)
    def BeginNbrAtom(*args): return _openbabel.OBAtomAtomIter_BeginNbrAtom(*args)
    def NextNbrAtom(*args): return _openbabel.OBAtomAtomIter_NextNbrAtom(*args)
    def GetDistance(*args): return _openbabel.OBAtomAtomIter_GetDistance(*args)
    def GetAngle(*args): return _openbabel.OBAtomAtomIter_GetAngle(*args)
    def NewResidue(*args): return _openbabel.OBAtomAtomIter_NewResidue(*args)
    def DeleteResidue(*args): return _openbabel.OBAtomAtomIter_DeleteResidue(*args)
    def AddBond(*args): return _openbabel.OBAtomAtomIter_AddBond(*args)
    def InsertBond(*args): return _openbabel.OBAtomAtomIter_InsertBond(*args)
    def DeleteBond(*args): return _openbabel.OBAtomAtomIter_DeleteBond(*args)
    def ClearBond(*args): return _openbabel.OBAtomAtomIter_ClearBond(*args)
    def CountFreeOxygens(*args): return _openbabel.OBAtomAtomIter_CountFreeOxygens(*args)
    def ImplicitHydrogenCount(*args): return _openbabel.OBAtomAtomIter_ImplicitHydrogenCount(*args)
    def ExplicitHydrogenCount(*args): return _openbabel.OBAtomAtomIter_ExplicitHydrogenCount(*args)
    def MemberOfRingCount(*args): return _openbabel.OBAtomAtomIter_MemberOfRingCount(*args)
    def MemberOfRingSize(*args): return _openbabel.OBAtomAtomIter_MemberOfRingSize(*args)
    def CountRingBonds(*args): return _openbabel.OBAtomAtomIter_CountRingBonds(*args)
    def SmallestBondAngle(*args): return _openbabel.OBAtomAtomIter_SmallestBondAngle(*args)
    def AverageBondAngle(*args): return _openbabel.OBAtomAtomIter_AverageBondAngle(*args)
    def BOSum(*args): return _openbabel.OBAtomAtomIter_BOSum(*args)
    def KBOSum(*args): return _openbabel.OBAtomAtomIter_KBOSum(*args)
    def HtoMethyl(*args): return _openbabel.OBAtomAtomIter_HtoMethyl(*args)
    def SetHybAndGeom(*args): return _openbabel.OBAtomAtomIter_SetHybAndGeom(*args)
    def ForceNoH(*args): return _openbabel.OBAtomAtomIter_ForceNoH(*args)
    def HasNoHForced(*args): return _openbabel.OBAtomAtomIter_HasNoHForced(*args)
    def HasResidue(*args): return _openbabel.OBAtomAtomIter_HasResidue(*args)
    def IsHydrogen(*args): return _openbabel.OBAtomAtomIter_IsHydrogen(*args)
    def IsCarbon(*args): return _openbabel.OBAtomAtomIter_IsCarbon(*args)
    def IsNitrogen(*args): return _openbabel.OBAtomAtomIter_IsNitrogen(*args)
    def IsOxygen(*args): return _openbabel.OBAtomAtomIter_IsOxygen(*args)
    def IsSulfur(*args): return _openbabel.OBAtomAtomIter_IsSulfur(*args)
    def IsPhosphorus(*args): return _openbabel.OBAtomAtomIter_IsPhosphorus(*args)
    def IsAromatic(*args): return _openbabel.OBAtomAtomIter_IsAromatic(*args)
    def IsInRing(*args): return _openbabel.OBAtomAtomIter_IsInRing(*args)
    def IsInRingSize(*args): return _openbabel.OBAtomAtomIter_IsInRingSize(*args)
    def IsHeteroatom(*args): return _openbabel.OBAtomAtomIter_IsHeteroatom(*args)
    def IsNotCorH(*args): return _openbabel.OBAtomAtomIter_IsNotCorH(*args)
    def IsConnected(*args): return _openbabel.OBAtomAtomIter_IsConnected(*args)
    def IsOneThree(*args): return _openbabel.OBAtomAtomIter_IsOneThree(*args)
    def IsOneFour(*args): return _openbabel.OBAtomAtomIter_IsOneFour(*args)
    def IsCarboxylOxygen(*args): return _openbabel.OBAtomAtomIter_IsCarboxylOxygen(*args)
    def IsPhosphateOxygen(*args): return _openbabel.OBAtomAtomIter_IsPhosphateOxygen(*args)
    def IsSulfateOxygen(*args): return _openbabel.OBAtomAtomIter_IsSulfateOxygen(*args)
    def IsNitroOxygen(*args): return _openbabel.OBAtomAtomIter_IsNitroOxygen(*args)
    def IsAmideNitrogen(*args): return _openbabel.OBAtomAtomIter_IsAmideNitrogen(*args)
    def IsPolarHydrogen(*args): return _openbabel.OBAtomAtomIter_IsPolarHydrogen(*args)
    def IsNonPolarHydrogen(*args): return _openbabel.OBAtomAtomIter_IsNonPolarHydrogen(*args)
    def IsAromaticNOxide(*args): return _openbabel.OBAtomAtomIter_IsAromaticNOxide(*args)
    def IsChiral(*args): return _openbabel.OBAtomAtomIter_IsChiral(*args)
    def IsAxial(*args): return _openbabel.OBAtomAtomIter_IsAxial(*args)
    def IsClockwise(*args): return _openbabel.OBAtomAtomIter_IsClockwise(*args)
    def IsAntiClockwise(*args): return _openbabel.OBAtomAtomIter_IsAntiClockwise(*args)
    def IsPositiveStereo(*args): return _openbabel.OBAtomAtomIter_IsPositiveStereo(*args)
    def IsNegativeStereo(*args): return _openbabel.OBAtomAtomIter_IsNegativeStereo(*args)
    def HasChiralitySpecified(*args): return _openbabel.OBAtomAtomIter_HasChiralitySpecified(*args)
    def HasChiralVolume(*args): return _openbabel.OBAtomAtomIter_HasChiralVolume(*args)
    def IsHbondAcceptor(*args): return _openbabel.OBAtomAtomIter_IsHbondAcceptor(*args)
    def IsHbondDonor(*args): return _openbabel.OBAtomAtomIter_IsHbondDonor(*args)
    def IsHbondDonorH(*args): return _openbabel.OBAtomAtomIter_IsHbondDonorH(*args)
    def HasAlphaBetaUnsat(*args): return _openbabel.OBAtomAtomIter_HasAlphaBetaUnsat(*args)
    def HasBondOfOrder(*args): return _openbabel.OBAtomAtomIter_HasBondOfOrder(*args)
    def CountBondsOfOrder(*args): return _openbabel.OBAtomAtomIter_CountBondsOfOrder(*args)
    def HasNonSingleBond(*args): return _openbabel.OBAtomAtomIter_HasNonSingleBond(*args)
    def HasSingleBond(*args): return _openbabel.OBAtomAtomIter_HasSingleBond(*args)
    def HasDoubleBond(*args): return _openbabel.OBAtomAtomIter_HasDoubleBond(*args)
    def HasAromaticBond(*args): return _openbabel.OBAtomAtomIter_HasAromaticBond(*args)
    def MatchesSMARTS(*args): return _openbabel.OBAtomAtomIter_MatchesSMARTS(*args)
    def DoTransformations(*args): return _openbabel.OBAtomAtomIter_DoTransformations(*args)
    def ClassDescription(*args): return _openbabel.OBAtomAtomIter_ClassDescription(*args)
    def HasData(*args): return _openbabel.OBAtomAtomIter_HasData(*args)
    def DeleteData(*args): return _openbabel.OBAtomAtomIter_DeleteData(*args)
    def SetData(*args): return _openbabel.OBAtomAtomIter_SetData(*args)
    def DataSize(*args): return _openbabel.OBAtomAtomIter_DataSize(*args)
    def GetData(*args): return _openbabel.OBAtomAtomIter_GetData(*args)
    def BeginData(*args): return _openbabel.OBAtomAtomIter_BeginData(*args)
    def EndData(*args): return _openbabel.OBAtomAtomIter_EndData(*args)
OBAtomAtomIter_swigregister = _openbabel.OBAtomAtomIter_swigregister
OBAtomAtomIter_swigregister(OBAtomAtomIter)

class OBAtomBondIter(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBAtomBondIter_swiginit(self,_openbabel.new_OBAtomBondIter(*args))
    def good(*args): return _openbabel.OBAtomBondIter_good(*args)
    def inc(*args): return _openbabel.OBAtomBondIter_inc(*args)
    def deref(*args): return _openbabel.OBAtomBondIter_deref(*args)
    def __ref__(*args): return _openbabel.OBAtomBondIter___ref__(*args)
    __swig_destroy__ = _openbabel.delete_OBAtomBondIter
    __del__ = lambda self : None;
    Visit = _swig_property(_openbabel.OBAtomBondIter_Visit_get, _openbabel.OBAtomBondIter_Visit_set)
    def SetIdx(*args): return _openbabel.OBAtomBondIter_SetIdx(*args)
    def SetBO(*args): return _openbabel.OBAtomBondIter_SetBO(*args)
    def SetBegin(*args): return _openbabel.OBAtomBondIter_SetBegin(*args)
    def SetEnd(*args): return _openbabel.OBAtomBondIter_SetEnd(*args)
    def SetParent(*args): return _openbabel.OBAtomBondIter_SetParent(*args)
    def SetLength(*args): return _openbabel.OBAtomBondIter_SetLength(*args)
    def Set(*args): return _openbabel.OBAtomBondIter_Set(*args)
    def SetKSingle(*args): return _openbabel.OBAtomBondIter_SetKSingle(*args)
    def SetKDouble(*args): return _openbabel.OBAtomBondIter_SetKDouble(*args)
    def SetKTriple(*args): return _openbabel.OBAtomBondIter_SetKTriple(*args)
    def SetAromatic(*args): return _openbabel.OBAtomBondIter_SetAromatic(*args)
    def SetHash(*args): return _openbabel.OBAtomBondIter_SetHash(*args)
    def SetWedge(*args): return _openbabel.OBAtomBondIter_SetWedge(*args)
    def SetUp(*args): return _openbabel.OBAtomBondIter_SetUp(*args)
    def SetDown(*args): return _openbabel.OBAtomBondIter_SetDown(*args)
    def SetInRing(*args): return _openbabel.OBAtomBondIter_SetInRing(*args)
    def SetClosure(*args): return _openbabel.OBAtomBondIter_SetClosure(*args)
    def UnsetHash(*args): return _openbabel.OBAtomBondIter_UnsetHash(*args)
    def UnsetWedge(*args): return _openbabel.OBAtomBondIter_UnsetWedge(*args)
    def UnsetUp(*args): return _openbabel.OBAtomBondIter_UnsetUp(*args)
    def UnsetDown(*args): return _openbabel.OBAtomBondIter_UnsetDown(*args)
    def UnsetAromatic(*args): return _openbabel.OBAtomBondIter_UnsetAromatic(*args)
    def UnsetKekule(*args): return _openbabel.OBAtomBondIter_UnsetKekule(*args)
    def GetIdx(*args): return _openbabel.OBAtomBondIter_GetIdx(*args)
    def GetBO(*args): return _openbabel.OBAtomBondIter_GetBO(*args)
    def GetBondOrder(*args): return _openbabel.OBAtomBondIter_GetBondOrder(*args)
    def GetFlags(*args): return _openbabel.OBAtomBondIter_GetFlags(*args)
    def GetBeginAtomIdx(*args): return _openbabel.OBAtomBondIter_GetBeginAtomIdx(*args)
    def GetEndAtomIdx(*args): return _openbabel.OBAtomBondIter_GetEndAtomIdx(*args)
    def GetBeginAtom(*args): return _openbabel.OBAtomBondIter_GetBeginAtom(*args)
    def GetEndAtom(*args): return _openbabel.OBAtomBondIter_GetEndAtom(*args)
    def GetNbrAtom(*args): return _openbabel.OBAtomBondIter_GetNbrAtom(*args)
    def GetParent(*args): return _openbabel.OBAtomBondIter_GetParent(*args)
    def GetEquibLength(*args): return _openbabel.OBAtomBondIter_GetEquibLength(*args)
    def GetLength(*args): return _openbabel.OBAtomBondIter_GetLength(*args)
    def GetNbrAtomIdx(*args): return _openbabel.OBAtomBondIter_GetNbrAtomIdx(*args)
    def IsAromatic(*args): return _openbabel.OBAtomBondIter_IsAromatic(*args)
    def IsInRing(*args): return _openbabel.OBAtomBondIter_IsInRing(*args)
    def IsRotor(*args): return _openbabel.OBAtomBondIter_IsRotor(*args)
    def IsAmide(*args): return _openbabel.OBAtomBondIter_IsAmide(*args)
    def IsPrimaryAmide(*args): return _openbabel.OBAtomBondIter_IsPrimaryAmide(*args)
    def IsSecondaryAmide(*args): return _openbabel.OBAtomBondIter_IsSecondaryAmide(*args)
    def IsEster(*args): return _openbabel.OBAtomBondIter_IsEster(*args)
    def IsCarbonyl(*args): return _openbabel.OBAtomBondIter_IsCarbonyl(*args)
    def IsSingle(*args): return _openbabel.OBAtomBondIter_IsSingle(*args)
    def IsDouble(*args): return _openbabel.OBAtomBondIter_IsDouble(*args)
    def IsTriple(*args): return _openbabel.OBAtomBondIter_IsTriple(*args)
    def IsKSingle(*args): return _openbabel.OBAtomBondIter_IsKSingle(*args)
    def IsKDouble(*args): return _openbabel.OBAtomBondIter_IsKDouble(*args)
    def IsKTriple(*args): return _openbabel.OBAtomBondIter_IsKTriple(*args)
    def IsClosure(*args): return _openbabel.OBAtomBondIter_IsClosure(*args)
    def IsUp(*args): return _openbabel.OBAtomBondIter_IsUp(*args)
    def IsDown(*args): return _openbabel.OBAtomBondIter_IsDown(*args)
    def IsWedge(*args): return _openbabel.OBAtomBondIter_IsWedge(*args)
    def IsHash(*args): return _openbabel.OBAtomBondIter_IsHash(*args)
    def IsDoubleBondGeometry(*args): return _openbabel.OBAtomBondIter_IsDoubleBondGeometry(*args)
    def DoTransformations(*args): return _openbabel.OBAtomBondIter_DoTransformations(*args)
    def ClassDescription(*args): return _openbabel.OBAtomBondIter_ClassDescription(*args)
    def HasData(*args): return _openbabel.OBAtomBondIter_HasData(*args)
    def DeleteData(*args): return _openbabel.OBAtomBondIter_DeleteData(*args)
    def SetData(*args): return _openbabel.OBAtomBondIter_SetData(*args)
    def DataSize(*args): return _openbabel.OBAtomBondIter_DataSize(*args)
    def GetData(*args): return _openbabel.OBAtomBondIter_GetData(*args)
    def BeginData(*args): return _openbabel.OBAtomBondIter_BeginData(*args)
    def EndData(*args): return _openbabel.OBAtomBondIter_EndData(*args)
OBAtomBondIter_swigregister = _openbabel.OBAtomBondIter_swigregister
OBAtomBondIter_swigregister(OBAtomBondIter)

class OBResidueIter(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBResidueIter_swiginit(self,_openbabel.new_OBResidueIter(*args))
    def good(*args): return _openbabel.OBResidueIter_good(*args)
    def inc(*args): return _openbabel.OBResidueIter_inc(*args)
    def __deref__(*args): return _openbabel.OBResidueIter___deref__(*args)
    def __ref__(*args): return _openbabel.OBResidueIter___ref__(*args)
    __swig_destroy__ = _openbabel.delete_OBResidueIter
    __del__ = lambda self : None;
    def AddAtom(*args): return _openbabel.OBResidueIter_AddAtom(*args)
    def InsertAtom(*args): return _openbabel.OBResidueIter_InsertAtom(*args)
    def RemoveAtom(*args): return _openbabel.OBResidueIter_RemoveAtom(*args)
    def Clear(*args): return _openbabel.OBResidueIter_Clear(*args)
    def SetName(*args): return _openbabel.OBResidueIter_SetName(*args)
    def SetNum(*args): return _openbabel.OBResidueIter_SetNum(*args)
    def SetChain(*args): return _openbabel.OBResidueIter_SetChain(*args)
    def SetChainNum(*args): return _openbabel.OBResidueIter_SetChainNum(*args)
    def SetIdx(*args): return _openbabel.OBResidueIter_SetIdx(*args)
    def SetAtomID(*args): return _openbabel.OBResidueIter_SetAtomID(*args)
    def SetHetAtom(*args): return _openbabel.OBResidueIter_SetHetAtom(*args)
    def SetSerialNum(*args): return _openbabel.OBResidueIter_SetSerialNum(*args)
    def GetName(*args): return _openbabel.OBResidueIter_GetName(*args)
    def GetNum(*args): return _openbabel.OBResidueIter_GetNum(*args)
    def GetNumAtoms(*args): return _openbabel.OBResidueIter_GetNumAtoms(*args)
    def GetChain(*args): return _openbabel.OBResidueIter_GetChain(*args)
    def GetChainNum(*args): return _openbabel.OBResidueIter_GetChainNum(*args)
    def GetIdx(*args): return _openbabel.OBResidueIter_GetIdx(*args)
    def GetResKey(*args): return _openbabel.OBResidueIter_GetResKey(*args)
    def GetAtoms(*args): return _openbabel.OBResidueIter_GetAtoms(*args)
    def GetBonds(*args): return _openbabel.OBResidueIter_GetBonds(*args)
    def GetAtomID(*args): return _openbabel.OBResidueIter_GetAtomID(*args)
    def GetSerialNum(*args): return _openbabel.OBResidueIter_GetSerialNum(*args)
    def GetAminoAcidProperty(*args): return _openbabel.OBResidueIter_GetAminoAcidProperty(*args)
    def GetAtomProperty(*args): return _openbabel.OBResidueIter_GetAtomProperty(*args)
    def GetResidueProperty(*args): return _openbabel.OBResidueIter_GetResidueProperty(*args)
    def IsHetAtom(*args): return _openbabel.OBResidueIter_IsHetAtom(*args)
    def IsResidueType(*args): return _openbabel.OBResidueIter_IsResidueType(*args)
    def BeginAtoms(*args): return _openbabel.OBResidueIter_BeginAtoms(*args)
    def EndAtoms(*args): return _openbabel.OBResidueIter_EndAtoms(*args)
    def BeginAtom(*args): return _openbabel.OBResidueIter_BeginAtom(*args)
    def NextAtom(*args): return _openbabel.OBResidueIter_NextAtom(*args)
    def DoTransformations(*args): return _openbabel.OBResidueIter_DoTransformations(*args)
    def ClassDescription(*args): return _openbabel.OBResidueIter_ClassDescription(*args)
    def HasData(*args): return _openbabel.OBResidueIter_HasData(*args)
    def DeleteData(*args): return _openbabel.OBResidueIter_DeleteData(*args)
    def SetData(*args): return _openbabel.OBResidueIter_SetData(*args)
    def DataSize(*args): return _openbabel.OBResidueIter_DataSize(*args)
    def GetData(*args): return _openbabel.OBResidueIter_GetData(*args)
    def BeginData(*args): return _openbabel.OBResidueIter_BeginData(*args)
    def EndData(*args): return _openbabel.OBResidueIter_EndData(*args)
OBResidueIter_swigregister = _openbabel.OBResidueIter_swigregister
OBResidueIter_swigregister(OBResidueIter)

class OBResidueAtomIter(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.OBResidueAtomIter_swiginit(self,_openbabel.new_OBResidueAtomIter(*args))
    def good(*args): return _openbabel.OBResidueAtomIter_good(*args)
    def inc(*args): return _openbabel.OBResidueAtomIter_inc(*args)
    def deref(*args): return _openbabel.OBResidueAtomIter_deref(*args)
    def __ref__(*args): return _openbabel.OBResidueAtomIter___ref__(*args)
    __swig_destroy__ = _openbabel.delete_OBResidueAtomIter
    __del__ = lambda self : None;
    Visit = _swig_property(_openbabel.OBResidueAtomIter_Visit_get, _openbabel.OBResidueAtomIter_Visit_set)
    def Clear(*args): return _openbabel.OBResidueAtomIter_Clear(*args)
    def SetIdx(*args): return _openbabel.OBResidueAtomIter_SetIdx(*args)
    def SetHyb(*args): return _openbabel.OBResidueAtomIter_SetHyb(*args)
    def SetAtomicNum(*args): return _openbabel.OBResidueAtomIter_SetAtomicNum(*args)
    def SetIsotope(*args): return _openbabel.OBResidueAtomIter_SetIsotope(*args)
    def SetImplicitValence(*args): return _openbabel.OBResidueAtomIter_SetImplicitValence(*args)
    def IncrementImplicitValence(*args): return _openbabel.OBResidueAtomIter_IncrementImplicitValence(*args)
    def DecrementImplicitValence(*args): return _openbabel.OBResidueAtomIter_DecrementImplicitValence(*args)
    def SetFormalCharge(*args): return _openbabel.OBResidueAtomIter_SetFormalCharge(*args)
    def SetSpinMultiplicity(*args): return _openbabel.OBResidueAtomIter_SetSpinMultiplicity(*args)
    def SetType(*args): return _openbabel.OBResidueAtomIter_SetType(*args)
    def SetPartialCharge(*args): return _openbabel.OBResidueAtomIter_SetPartialCharge(*args)
    def SetVector(*args): return _openbabel.OBResidueAtomIter_SetVector(*args)
    def SetCoordPtr(*args): return _openbabel.OBResidueAtomIter_SetCoordPtr(*args)
    def SetResidue(*args): return _openbabel.OBResidueAtomIter_SetResidue(*args)
    def SetParent(*args): return _openbabel.OBResidueAtomIter_SetParent(*args)
    def SetAromatic(*args): return _openbabel.OBResidueAtomIter_SetAromatic(*args)
    def UnsetAromatic(*args): return _openbabel.OBResidueAtomIter_UnsetAromatic(*args)
    def SetClockwiseStereo(*args): return _openbabel.OBResidueAtomIter_SetClockwiseStereo(*args)
    def SetAntiClockwiseStereo(*args): return _openbabel.OBResidueAtomIter_SetAntiClockwiseStereo(*args)
    def SetPositiveStereo(*args): return _openbabel.OBResidueAtomIter_SetPositiveStereo(*args)
    def SetNegativeStereo(*args): return _openbabel.OBResidueAtomIter_SetNegativeStereo(*args)
    def UnsetStereo(*args): return _openbabel.OBResidueAtomIter_UnsetStereo(*args)
    def SetInRing(*args): return _openbabel.OBResidueAtomIter_SetInRing(*args)
    def SetChiral(*args): return _openbabel.OBResidueAtomIter_SetChiral(*args)
    def ClearCoordPtr(*args): return _openbabel.OBResidueAtomIter_ClearCoordPtr(*args)
    def GetFormalCharge(*args): return _openbabel.OBResidueAtomIter_GetFormalCharge(*args)
    def GetAtomicNum(*args): return _openbabel.OBResidueAtomIter_GetAtomicNum(*args)
    def GetIsotope(*args): return _openbabel.OBResidueAtomIter_GetIsotope(*args)
    def GetSpinMultiplicity(*args): return _openbabel.OBResidueAtomIter_GetSpinMultiplicity(*args)
    def GetAtomicMass(*args): return _openbabel.OBResidueAtomIter_GetAtomicMass(*args)
    def GetExactMass(*args): return _openbabel.OBResidueAtomIter_GetExactMass(*args)
    def GetIdx(*args): return _openbabel.OBResidueAtomIter_GetIdx(*args)
    def GetCoordinateIdx(*args): return _openbabel.OBResidueAtomIter_GetCoordinateIdx(*args)
    def GetCIdx(*args): return _openbabel.OBResidueAtomIter_GetCIdx(*args)
    def GetValence(*args): return _openbabel.OBResidueAtomIter_GetValence(*args)
    def GetHyb(*args): return _openbabel.OBResidueAtomIter_GetHyb(*args)
    def GetImplicitValence(*args): return _openbabel.OBResidueAtomIter_GetImplicitValence(*args)
    def GetHvyValence(*args): return _openbabel.OBResidueAtomIter_GetHvyValence(*args)
    def GetHeteroValence(*args): return _openbabel.OBResidueAtomIter_GetHeteroValence(*args)
    def GetType(*args): return _openbabel.OBResidueAtomIter_GetType(*args)
    def GetX(*args): return _openbabel.OBResidueAtomIter_GetX(*args)
    def x(*args): return _openbabel.OBResidueAtomIter_x(*args)
    def GetY(*args): return _openbabel.OBResidueAtomIter_GetY(*args)
    def y(*args): return _openbabel.OBResidueAtomIter_y(*args)
    def GetZ(*args): return _openbabel.OBResidueAtomIter_GetZ(*args)
    def z(*args): return _openbabel.OBResidueAtomIter_z(*args)
    def GetCoordinate(*args): return _openbabel.OBResidueAtomIter_GetCoordinate(*args)
    def GetVector(*args): return _openbabel.OBResidueAtomIter_GetVector(*args)
    def GetPartialCharge(*args): return _openbabel.OBResidueAtomIter_GetPartialCharge(*args)
    def GetResidue(*args): return _openbabel.OBResidueAtomIter_GetResidue(*args)
    def GetParent(*args): return _openbabel.OBResidueAtomIter_GetParent(*args)
    def GetNewBondVector(*args): return _openbabel.OBResidueAtomIter_GetNewBondVector(*args)
    def GetBond(*args): return _openbabel.OBResidueAtomIter_GetBond(*args)
    def GetNextAtom(*args): return _openbabel.OBResidueAtomIter_GetNextAtom(*args)
    def BeginBonds(*args): return _openbabel.OBResidueAtomIter_BeginBonds(*args)
    def EndBonds(*args): return _openbabel.OBResidueAtomIter_EndBonds(*args)
    def BeginBond(*args): return _openbabel.OBResidueAtomIter_BeginBond(*args)
    def NextBond(*args): return _openbabel.OBResidueAtomIter_NextBond(*args)
    def BeginNbrAtom(*args): return _openbabel.OBResidueAtomIter_BeginNbrAtom(*args)
    def NextNbrAtom(*args): return _openbabel.OBResidueAtomIter_NextNbrAtom(*args)
    def GetDistance(*args): return _openbabel.OBResidueAtomIter_GetDistance(*args)
    def GetAngle(*args): return _openbabel.OBResidueAtomIter_GetAngle(*args)
    def NewResidue(*args): return _openbabel.OBResidueAtomIter_NewResidue(*args)
    def DeleteResidue(*args): return _openbabel.OBResidueAtomIter_DeleteResidue(*args)
    def AddBond(*args): return _openbabel.OBResidueAtomIter_AddBond(*args)
    def InsertBond(*args): return _openbabel.OBResidueAtomIter_InsertBond(*args)
    def DeleteBond(*args): return _openbabel.OBResidueAtomIter_DeleteBond(*args)
    def ClearBond(*args): return _openbabel.OBResidueAtomIter_ClearBond(*args)
    def CountFreeOxygens(*args): return _openbabel.OBResidueAtomIter_CountFreeOxygens(*args)
    def ImplicitHydrogenCount(*args): return _openbabel.OBResidueAtomIter_ImplicitHydrogenCount(*args)
    def ExplicitHydrogenCount(*args): return _openbabel.OBResidueAtomIter_ExplicitHydrogenCount(*args)
    def MemberOfRingCount(*args): return _openbabel.OBResidueAtomIter_MemberOfRingCount(*args)
    def MemberOfRingSize(*args): return _openbabel.OBResidueAtomIter_MemberOfRingSize(*args)
    def CountRingBonds(*args): return _openbabel.OBResidueAtomIter_CountRingBonds(*args)
    def SmallestBondAngle(*args): return _openbabel.OBResidueAtomIter_SmallestBondAngle(*args)
    def AverageBondAngle(*args): return _openbabel.OBResidueAtomIter_AverageBondAngle(*args)
    def BOSum(*args): return _openbabel.OBResidueAtomIter_BOSum(*args)
    def KBOSum(*args): return _openbabel.OBResidueAtomIter_KBOSum(*args)
    def HtoMethyl(*args): return _openbabel.OBResidueAtomIter_HtoMethyl(*args)
    def SetHybAndGeom(*args): return _openbabel.OBResidueAtomIter_SetHybAndGeom(*args)
    def ForceNoH(*args): return _openbabel.OBResidueAtomIter_ForceNoH(*args)
    def HasNoHForced(*args): return _openbabel.OBResidueAtomIter_HasNoHForced(*args)
    def HasResidue(*args): return _openbabel.OBResidueAtomIter_HasResidue(*args)
    def IsHydrogen(*args): return _openbabel.OBResidueAtomIter_IsHydrogen(*args)
    def IsCarbon(*args): return _openbabel.OBResidueAtomIter_IsCarbon(*args)
    def IsNitrogen(*args): return _openbabel.OBResidueAtomIter_IsNitrogen(*args)
    def IsOxygen(*args): return _openbabel.OBResidueAtomIter_IsOxygen(*args)
    def IsSulfur(*args): return _openbabel.OBResidueAtomIter_IsSulfur(*args)
    def IsPhosphorus(*args): return _openbabel.OBResidueAtomIter_IsPhosphorus(*args)
    def IsAromatic(*args): return _openbabel.OBResidueAtomIter_IsAromatic(*args)
    def IsInRing(*args): return _openbabel.OBResidueAtomIter_IsInRing(*args)
    def IsInRingSize(*args): return _openbabel.OBResidueAtomIter_IsInRingSize(*args)
    def IsHeteroatom(*args): return _openbabel.OBResidueAtomIter_IsHeteroatom(*args)
    def IsNotCorH(*args): return _openbabel.OBResidueAtomIter_IsNotCorH(*args)
    def IsConnected(*args): return _openbabel.OBResidueAtomIter_IsConnected(*args)
    def IsOneThree(*args): return _openbabel.OBResidueAtomIter_IsOneThree(*args)
    def IsOneFour(*args): return _openbabel.OBResidueAtomIter_IsOneFour(*args)
    def IsCarboxylOxygen(*args): return _openbabel.OBResidueAtomIter_IsCarboxylOxygen(*args)
    def IsPhosphateOxygen(*args): return _openbabel.OBResidueAtomIter_IsPhosphateOxygen(*args)
    def IsSulfateOxygen(*args): return _openbabel.OBResidueAtomIter_IsSulfateOxygen(*args)
    def IsNitroOxygen(*args): return _openbabel.OBResidueAtomIter_IsNitroOxygen(*args)
    def IsAmideNitrogen(*args): return _openbabel.OBResidueAtomIter_IsAmideNitrogen(*args)
    def IsPolarHydrogen(*args): return _openbabel.OBResidueAtomIter_IsPolarHydrogen(*args)
    def IsNonPolarHydrogen(*args): return _openbabel.OBResidueAtomIter_IsNonPolarHydrogen(*args)
    def IsAromaticNOxide(*args): return _openbabel.OBResidueAtomIter_IsAromaticNOxide(*args)
    def IsChiral(*args): return _openbabel.OBResidueAtomIter_IsChiral(*args)
    def IsAxial(*args): return _openbabel.OBResidueAtomIter_IsAxial(*args)
    def IsClockwise(*args): return _openbabel.OBResidueAtomIter_IsClockwise(*args)
    def IsAntiClockwise(*args): return _openbabel.OBResidueAtomIter_IsAntiClockwise(*args)
    def IsPositiveStereo(*args): return _openbabel.OBResidueAtomIter_IsPositiveStereo(*args)
    def IsNegativeStereo(*args): return _openbabel.OBResidueAtomIter_IsNegativeStereo(*args)
    def HasChiralitySpecified(*args): return _openbabel.OBResidueAtomIter_HasChiralitySpecified(*args)
    def HasChiralVolume(*args): return _openbabel.OBResidueAtomIter_HasChiralVolume(*args)
    def IsHbondAcceptor(*args): return _openbabel.OBResidueAtomIter_IsHbondAcceptor(*args)
    def IsHbondDonor(*args): return _openbabel.OBResidueAtomIter_IsHbondDonor(*args)
    def IsHbondDonorH(*args): return _openbabel.OBResidueAtomIter_IsHbondDonorH(*args)
    def HasAlphaBetaUnsat(*args): return _openbabel.OBResidueAtomIter_HasAlphaBetaUnsat(*args)
    def HasBondOfOrder(*args): return _openbabel.OBResidueAtomIter_HasBondOfOrder(*args)
    def CountBondsOfOrder(*args): return _openbabel.OBResidueAtomIter_CountBondsOfOrder(*args)
    def HasNonSingleBond(*args): return _openbabel.OBResidueAtomIter_HasNonSingleBond(*args)
    def HasSingleBond(*args): return _openbabel.OBResidueAtomIter_HasSingleBond(*args)
    def HasDoubleBond(*args): return _openbabel.OBResidueAtomIter_HasDoubleBond(*args)
    def HasAromaticBond(*args): return _openbabel.OBResidueAtomIter_HasAromaticBond(*args)
    def MatchesSMARTS(*args): return _openbabel.OBResidueAtomIter_MatchesSMARTS(*args)
    def DoTransformations(*args): return _openbabel.OBResidueAtomIter_DoTransformations(*args)
    def ClassDescription(*args): return _openbabel.OBResidueAtomIter_ClassDescription(*args)
    def HasData(*args): return _openbabel.OBResidueAtomIter_HasData(*args)
    def DeleteData(*args): return _openbabel.OBResidueAtomIter_DeleteData(*args)
    def SetData(*args): return _openbabel.OBResidueAtomIter_SetData(*args)
    def DataSize(*args): return _openbabel.OBResidueAtomIter_DataSize(*args)
    def GetData(*args): return _openbabel.OBResidueAtomIter_GetData(*args)
    def BeginData(*args): return _openbabel.OBResidueAtomIter_BeginData(*args)
    def EndData(*args): return _openbabel.OBResidueAtomIter_EndData(*args)
OBResidueAtomIter_swigregister = _openbabel.OBResidueAtomIter_swigregister
OBResidueAtomIter_swigregister(OBResidueAtomIter)

class doubleArray(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    def __init__(self, *args): 
        _openbabel.doubleArray_swiginit(self,_openbabel.new_doubleArray(*args))
    __swig_destroy__ = _openbabel.delete_doubleArray
    __del__ = lambda self : None;
    def __getitem__(*args): return _openbabel.doubleArray___getitem__(*args)
    def __setitem__(*args): return _openbabel.doubleArray___setitem__(*args)
    def cast(*args): return _openbabel.doubleArray_cast(*args)
    frompointer = staticmethod(_openbabel.doubleArray_frompointer)
doubleArray_swigregister = _openbabel.doubleArray_swigregister
doubleArray_swigregister(doubleArray)
doubleArray_frompointer = _openbabel.doubleArray_frompointer

def double_array(mylist):
    """Create a C array of doubles from a list."""
    c = doubleArray(len(mylist))
    for i,v in enumerate(mylist):
        c[i] = v
    return c



