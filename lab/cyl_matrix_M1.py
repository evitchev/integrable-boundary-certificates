"""M1's fixed exact-point observable test, OpenAI Codex, 2026-09-18.

Exit2: failure of a mandatory M1 precondition, not a matrix-class exclusion.
Exit1: failed internal check/control; exit64: usage error.
The h1 input-formula audit is conditional physics and not a stage2 candidate test.
Optional --numerical checks the matched local asymptotic using mpmath only.
"""
from pathlib import Path
import argparse,json,sys,time
import sympy as S

class UsageParser(argparse.ArgumentParser):
    def error(self,message):
        self.print_usage(sys.stderr);self.exit(64,f'{self.prog}: error: {message}\n')

def save(out,name,obj):
    (out/name).write_text(json.dumps(obj,indent=2)+'\n')

def exact(out,tamper=False):
    u,v,k=S.symbols('u v k',positive=True);I=S.I;R=S.Rational
    V=(16*u**2+12*u+1)/(16*(1+u)**2);Lam=u**R(3,2)/(1+u)**R(5,2)
    psi=u**R(1,4)*(1+u)**R(-1,4)*(1+2*u/3)
    residual=S.simplify(u*S.diff(u*S.diff(psi,u),u)/psi-V)
    if residual!=0:raise ArithmeticError('zero-energy solution')
    a0=S.limit(psi/u**R(1,4),u,0);ainf=S.limit(psi/u,u,S.oo)
    if a0!=1 or ainf!=R(2,3):raise ArithmeticError('endpoint normalization')
    local_c=S.limit(v**2*V.subs(u,v-1)/(v-1)**2,v,0)
    nu2=S.simplify(16*(local_c+R(1,4)))
    if local_c!=R(5,16) or nu2!=9:raise ArithmeticError('Bessel order')
    phase=S.exp(I*S.pi/4)
    a_minus=phase/3;w0plus=-R(3,2)*a_minus
    N=2*S.sqrt(2/S.pi)*(2 if tamper else 1);Klead=S.Integer(8)
    if S.simplify(N*S.sqrt(S.pi)/(2*S.sqrt(2)))!=1:
        raise ArithmeticError('canonical WKB prefactor must be one; endpoint normalization tamper detected')
    cy=S.simplify(N*Klead/(4*S.exp(-I*S.pi/4))**3)
    cW=S.simplify(cy*w0plus)
    if S.simplify(cW-S.sqrt(2)/(8*S.sqrt(S.pi)))!=0:raise ArithmeticError('connection prefactor')
    # Verify the local Bessel operator by the chain rule before using its
    # small-argument asymptotic. x=4 c kappa v^(-1/4), c^2=-i.
    x=S.symbols('x');f0,f1,f2=S.symbols('f0 f1 f2')
    leading_v2_ypp=-f0/4+x*f1/16+x**2*f2/16
    local_equation=S.expand(leading_v2_ypp-(x*x/16+local_c)*f0)
    expected_bessel=(x*x*f2+x*f1-(x*x+9)*f0)/16
    if S.simplify(local_equation-expected_bessel)!=0:raise ArithmeticError('local Bessel operator balance')
    # Scalar reduction of [[a,b],[c,-a]] with y1=sqrt(b)*psi.
    z=S.symbols('z');aa=S.Function('a')(z);bb=S.Function('b')(z);cc=S.Function('c')(z);f=S.Function('f')(z)
    yy=S.sqrt(bb)*f
    scalar_raw=S.diff(yy,z,2)-S.diff(bb,z)/bb*S.diff(yy,z)-(S.diff(aa,z)+aa**2+bb*cc-aa*S.diff(bb,z)/bb)*yy
    pot=aa**2+bb*cc+S.diff(aa,z)-aa*S.diff(bb,z)/bb-S.diff(bb,z,2)/(2*bb)+3*S.diff(bb,z)**2/(4*bb**2)
    if S.simplify(scalar_raw/S.sqrt(bb)-(S.diff(f,z,2)-pot*f))!=0:raise ArithmeticError('scalar reduction')
    controls={'wrong_zero_energy_coefficient':str(S.factor(u*S.diff(u*S.diff(u**R(1,4)*(1+u)**R(-1,4)*(1+u/3),u),u)/(u**R(1,4)*(1+u)**R(-1,4)*(1+u/3))-V)),
              'drop_local_inverse_square_term_changes_nu_squared':str(16*R(1,4)),
              'multiply_Y_solution_by_kappa_cubed':'changes leading canonical WKB coefficient by kappa^3; forbidden energy-dependent endpoint renormalization, not an admissible repair',
              'unshifted_Q_product_pole_order':-5,'M1_shifted_quantum_Wronskian_pole_order':-4,'scalar_A1_pole_order':0,'companion_A1_pole_order':-1}
    if controls['wrong_zero_energy_coefficient']=='0':raise ArithmeticError('vacuous exact control')
    result={'status':'EXACT','head':'da1e2df','fibre':{'solution':1,'t':'-1/3','also_solution2_t':'1/3','X':'-2','Y':'0'},
      'V':str(V),'Lambda':str(Lam),'zero_energy_solution':str(psi),'zero_energy_residual':str(residual),'zero_energy_normalization':str(a0),'growing_infinity_coefficient':str(ainf),
      'scalar_A1_W_at_kappa0':'-4/3','at_minus1':{'local_inverse_square_coefficient':str(local_c),'outer_exponents':['5/4','-1/4'],'inner_variable':'4 exp(-i*pi/4) kappa v^(-1/4), v=1+u','Bessel_order_squared':str(nu2),'local_Bessel_operator_residual':'0','canonical_Bessel_prefactor':str(N),'limit_x3_K3':str(Klead),'psi0_negative_exponent_coefficient':str(a_minus),'Wz_psi0_fplus':str(w0plus),'limit_kappa3_psiY_prefactor':str(cy),'limit_kappa3_scalar_W0Y':str(cW),'limit_kappa4_QY':str(cW)},
      'nonvanishing_QX_at_kappa_minus1':{'status':'EXACT analytic positivity argument','potential_numerator':[16,12,1],'potential_denominator':'16(1+u)^2','argument':'For u>0, Lambda+V>0. A common decaying solution at both z ends for kappa^2=1 would make integral(|psi_z|^2+(Lambda+V)|psi|^2) dz=0, impossible. Scalar endpoint normalizations depend only on kappa^2, so QX(-1)=-W_A1(1) !=0.'},
      'M1_small_kappa':'D_M1 ~ [sqrt(2)/(4 sqrt(pi))] QX(kappa=-1)/(eta_X-eta_Y) * kappa^(-4); second shifted product is at most kappa^(-1)',
      'obstruction':'canonical small-kappa pole order of mandatory exact-point quantum Wronskian differs from A1 connection determinant; all B_j vanish at this fibre so no joint coefficient solve can change it',
      'admissible_class':'empty with all sealed preconditions','dimension':'not applicable: empty set, not a zero-dimensional family','gauge_quotient':'empty; allowed endpoint-normalized SL2 gauges preserve connection determinants',
      'controls':controls,'exit':2,'stage2':'NOT DONE: no surviving candidate','stage3':'NOT DONE: earlier gate failed','scope':'sealed M1 observable, endpoint frames and additive E+4 shift only; no general matrix or scalar exclusion'}
    save(out,'stage1_exact.json',result);print('EXACT M1 stage1: fixed exact-point pole-order obstruction; exit2.',flush=True)
    return result

def target_audit(out):
    mu,L=S.symbols('mu L');H=mu*(mu**2-1)/3
    RX=mu*(11*mu**2+7)/72;RY=(11*mu**4-20*mu**2+9)/(72*mu)
    B0=L+S.EulerGamma
    res0=[S.residue(H*B0+R,mu,0) for R in (RX,RY)]
    val1=[S.simplify(-S.diff(H,mu).subs(mu,1)/2+R.subs(mu,1)) for R in (RX,RY)]
    if res0!=[0,S.Rational(1,8)] or val1!=[-S.Rational(1,12),-S.Rational(1,3)]:raise ArithmeticError('target local audit')
    n=S.symbols('n',integer=True,positive=True)
    resn=S.factor(-H.subs(mu,n)/2)
    result={'status':'EXACT identity audit of the displayed formula, not a candidate test','physics_conditionality':'h1 target inherits BLZ II section4 conjectures per brief; not independently verified from literature here',
      'reflection_rewrite':'log kappa-log2 - (digamma(1+mu)+digamma(1-mu))/2',
      'bracket_at_zero':'log kappa-log2+EulerGamma; analytic even function near zero','mu0_residues':{'X':str(res0[0]),'Y':str(res0[1])},'mu1_values':{'X':str(val1[0]),'Y':str(val1[1])},
      'positive_integer_residue_for_n_ge2':str(resn),'sample_residues':{str(j):str(resn.subs(n,j)) for j in range(2,9)},
      'finding':'The displayed Y formula has an extra simple pole at mu_Y=0 of residue1/8. The no-pole-at-zero sentence is false for this formula. Integer n>=2 residue statement and no pole at n=1 hold.'}
    import mpmath as mp
    with mp.workdps(65):
        def value(x,sector):
            h=x*(x*x-1)/3;b=mp.log(3)-mp.log(2)-mp.digamma(x)-mp.pi*mp.cot(mp.pi*x)/2-1/(2*x)
            r=x*(11*x*x+7)/72 if sector=='X' else (11*x**4-20*x*x+9)/(72*x)
            return h*b+r
        eps=mp.mpf('1e-24')
        result['numerical_controls']={'label':'NUMERICAL (mpmath), corroboration only','precision':65,'epsilon':str(eps),'mu_times_X_near0':mp.nstr(eps*value(eps,'X'),35),'mu_times_Y_near0':mp.nstr(eps*value(eps,'Y'),35),'integer_residues_Y':{str(j):mp.nstr(eps*value(j+eps,'Y'),35) for j in range(2,9)}}
        if abs(eps*value(eps,'Y')-mp.mpf(1)/8)>mp.mpf('1e-40'):raise ArithmeticError('numerical zero residue audit')
    save(out,'target_formula_audit.json',result);print('EXACT input correction: Y residue at mu0 is 1/8; no-pole-at-zero claim fails as printed.',flush=True)
    return result

def numeric_check(out):
    """Two-boundary mpmath check of kappa^4 QY, independent of zero-energy solve."""
    import mpmath as mp
    mp.mp.dps=55
    def outer0(kappa,rad=mp.mpf('.5'),N=440):
        coeff=[mp.mpc(1)];qs=[mp.mpc(0)]*(N+1)
        for m in range(1,N//2+1):qs[2*m]=mp.mpf(5)*(m-3)*(-1)**m/16
        for j in range((N-3)//2+1):qs[3+2*j]=kappa*kappa*(-1)**j*mp.rf(mp.mpf('2.5'),j)/mp.factorial(j)
        for n in range(1,N+1):coeff.append(4*sum((qs[j]*coeff[n-j] for j in range(1,n+1)),mp.mpc(0))/(n*(n+1)))
        root=mp.j*mp.sqrt(rad);power=mp.mpc(1);val=mp.mpc(0);dz=mp.mpc(0)
        for n,a in enumerate(coeff):val+=a*power;dz+=(mp.mpf('.25')+mp.mpf(n)/2)*a*power;power*=root
        pref=rad**mp.mpf('.25')*mp.exp(mp.j*mp.pi/4)
        return pref*val,pref*dz
    def inner(kappa,Tabs,order=48,step=mp.mpf('.25')):
        v=(4*kappa/Tabs)**4;arg=Tabs*mp.exp(-mp.j*mp.pi/4);norm=2*mp.sqrt(2/mp.pi)
        kval=mp.besselk(3,arg);kp=-(mp.besselk(2,arg)+mp.besselk(4,arg))/2
        y=norm*mp.sqrt(v)*kval;dy=norm/mp.sqrt(v)*(kval/2-arg*kp/4)
        count=0;end=mp.mpf('.5')
        while v<end:
            hh=min(step,end/v-1);A=1-v;b=v/A
            pp=[v/A*b**j for j in range(order)]
            c1=[mp.rf(mp.mpf('.5'),j)/mp.factorial(j)*b**j for j in range(order)]
            c2=[(-1)**j*mp.rf(mp.mpf('2.5'),j)/mp.factorial(j) for j in range(order)]
            c3=[(1 if j==0 else 0)-mp.mpf(3)/(4*A)*b**j+mp.mpf(j+1)/(16*A*A)*b**j for j in range(order)]
            qq=[]
            for j in range(order):
                kinetic=-mp.j*kappa*kappa/mp.sqrt(v*A)*sum((c1[l]*c2[j-l] for l in range(j+1)),mp.mpc(0))
                potential=sum((c3[l]*(-1)**(j-l)*(j-l+1) for l in range(j+1)),mp.mpf(0))
                qq.append(kinetic+potential)
            cs=[y,v*dy]
            for n in range(order-1):
                cs.append(sum((pp[j]*(n-j+1)*cs[n-j+1]+qq[j]*cs[n-j] for j in range(n+1)),mp.mpc(0))/((n+2)*(n+1)))
            y=sum((c*hh**j for j,c in enumerate(cs)),mp.mpc(0));dy=sum((j*cs[j]*hh**(j-1) for j in range(1,len(cs))),mp.mpc(0))/v
            v*=1+hh;count+=1
        return y,dy,count
    result={'status':'NUMERICAL (mpmath); corroboration, not the analytic proof','precision':55,'outer_terms':440,'Taylor_order':48,'relative_step':'1/4','canonical_Bessel_boundary_approximation':'leading local equation only; two Tabs starts test truncation','expected_limit':mp.nstr(mp.sqrt(2)/(8*mp.sqrt(mp.pi)),45),'rows':[]};start=time.monotonic()
    for ks in ['0.04','0.02','0.01','0.005']:
        kap=mp.mpf(ks);y0,dz0=outer0(kap)
        for Ts in ['32','48']:
            y,dy,steps=inner(kap,mp.mpf(Ts));W=y0*(-mp.mpf('.5')*dy)-dz0*y
            scaled=kap**3*W
            row={'kappa':ks,'Tabs':Ts,'steps':steps,'kappa4_QY_real':mp.nstr(mp.re(scaled),40),'kappa4_QY_imag':mp.nstr(mp.im(scaled),40),'absolute_error_to_limit':mp.nstr(abs(scaled-mp.sqrt(2)/(8*mp.sqrt(mp.pi))),15)}
            result['rows'].append(row);save(out,'numerical_connection_check.json',result);print('NUMERICAL',json.dumps(row),flush=True)
    result['elapsed_seconds']=round(time.monotonic()-start,3);save(out,'numerical_connection_check.json',result)
    return result

def main():
    p=UsageParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--numerical',action='store_true');p.add_argument('--tamper-endpoint-prefactor',action='store_true');a=p.parse_args();out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
    stage=exact(out,a.tamper_endpoint_prefactor);audit=target_audit(out)
    if a.numerical:numeric_check(out)
    save(out,'summary.json',{'exit':2,'stage1':'EXACT M1 empty: mandatory exact-point observable has wrong small-kappa pole order','stage2':'NOT DONE: no surviving candidate','stage3':'NOT DONE: earlier gate failed','input_audit':'EXACT displayed Y target has residue1/8 at mu0; physical h1 formula remains conditional','scope':stage['scope']})
    return 2
if __name__=='__main__':sys.exit(main())
