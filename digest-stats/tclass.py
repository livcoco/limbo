#!/usr/bin/env python

class DigestResult:
    """Create a dummy digest result"""
    name = 'NilDigest'
    digest_size = 1
    def digest(self):      return '7'
    def hexdigest(self):   return '37'
    
class NilDigest:
    """Create a nil (dummy) digester"""
    def func(buf=[]):
        return DigestResult()
    __call__ = func

class NilDigest2:
    """Create a dummy digest result for 'empty-loop' timing"""
    name = 'NilDigest'
    digest_size = 1
    def digest(self):      return '7'
    def hexdigest(self):   return '37'
    def func(self, buf=[]):
        return self
    __call__ = func

class NilDigest3:
    """Create a dummy digest result for 'empty-loop' timing"""
    def __init__(self, buf=[]):  pass
    def digest(self):            return '7'
    def hexdigest(self):         return '37'
    def func(self):              return self
    name = 'NilDigest'
    digest_size = 1
    __call__ = func

p = NilDigest2()
print p
print p.name, p.digest_size, p.digest(), p.hexdigest()

q = NilDigest3('xyzzy')
print q.name, q.digest_size, q.digest(), q.hexdigest()
