from fractions import Fraction as Fr
class Q2:
    """elements a + b*sqrt(2), a,b in Q"""
    __slots__=("a","b")
    def __init__(self,a=0,b=0):
        self.a=Fr(a); self.b=Fr(b)
    @staticmethod
    def _c(x):
        if isinstance(x,Q2): return x
        return Q2(Fr(x),0)
    def __add__(s,o): o=Q2._c(o); return Q2(s.a+o.a,s.b+o.b)
    __radd__=__add__
    def __sub__(s,o): o=Q2._c(o); return Q2(s.a-o.a,s.b-o.b)
    def __rsub__(s,o): o=Q2._c(o); return Q2(o.a-s.a,o.b-s.b)
    def __neg__(s): return Q2(-s.a,-s.b)
    def __mul__(s,o): o=Q2._c(o); return Q2(s.a*o.a+2*s.b*o.b, s.a*o.b+s.b*o.a)
    __rmul__=__mul__
    def inv(s):
        n=s.a*s.a-2*s.b*s.b
        return Q2(s.a/n,-s.b/n)
    def __truediv__(s,o): o=Q2._c(o); return s*o.inv()
    def __rtruediv__(s,o): o=Q2._c(o); return o*s.inv()
    def __pow__(s,k):
        r=Q2(1,0)
        for _ in range(k): r=r*s
        return r
    def __eq__(s,o):
        o=Q2._c(o); return s.a==o.a and s.b==o.b
    def __ne__(s,o): return not s.__eq__(o)
    def __hash__(s): return hash((s.a,s.b))
    def __repr__(s):
        if s.b==0: return f"{s.a}"
        return f"({s.a}+{s.b}*sqrt2)"
SQ2=Q2(0,1)
