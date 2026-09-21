"""M2 connection observable and origin-response audit, OpenAI Codex, 2026-09-18.

EXACT symbolic identities; optional NUMERICAL Bessel-integral controls (mpmath).
Exit 3: global M2-to-radial connection construction remains unresolved.
This is not a candidate-oper confirmation or an exclusion of M2.
Exit 1: a computation/control fails; exit 64: usage error.
The h1 charge-hierarchy interpretation remains the brief's conditional input.
No certified-charge hold-outs or fitted amplitudes are used.
"""
from pathlib import Path
import argparse
import json
import sys
import sympy as S

R = S.Rational


class UsageParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_usage(sys.stderr)
        self.exit(64, f'{self.prog}: error: {message}\n')


def save(out, name, obj):
    (out / name).write_text(json.dumps(obj, indent=2) + '\n')


def check_zero(expr, label):
    reduced = S.factor(S.cancel(S.expand(expr)))
    if reduced != 0:
        raise ArithmeticError(f'{label}: {reduced}')


def scalar_and_frame(out, tamper):
    z, lam = S.symbols('z lambda', positive=True)
    a, b, c = [S.Function(s)(z) for s in ('a', 'b', 'c')]
    psi = S.Function('psi')(z)
    Q = a*a + b*c + S.diff(a,z) - a*S.diff(b,z)/b - S.diff(b,z,2)/(2*b) + 3*S.diff(b,z)**2/(4*b*b)
    y1 = S.sqrt(b/lam)*psi
    y2 = (S.diff(y1,z)-a*y1)/b
    residual = S.simplify((S.diff(y2,z)-c*y1+a*y2)*S.sqrt(lam*b))
    check_zero(S.simplify(residual - S.diff(psi,z,2)+Q*psi), 'scalarization')
    T = S.Matrix([[S.sqrt(b/lam),0],[(S.diff(b,z)/(2*b)-a)/S.sqrt(lam*b),1/S.sqrt(lam*b)]])
    determinant = S.simplify(T.det())
    frame_factor = (2 if tamper == 'frame' else 1)*lam
    check_zero(frame_factor*determinant-1,'derived matrix-to-scalar Wronskian')
    L,V = S.symbols('Lambda V')
    q0 = Q.subs({a:0,b:lam,c:lam*L+V/lam}).doit()
    check_zero(q0-lam**2*L-V,'A1 exact reduction')
    # Nonlinear admissible subfamily: diagonal f/lambda, same off-diagonals.
    f = S.symbols('f')
    qf = S.simplify(Q.subs({a:f/lam,b:lam,c:lam*L+V/lam}).doit())
    check_zero(qf-(lam**2*L+V+f*f/lam**2),'infinite-dimensional subfamily')
    u = S.symbols('u', positive=True)
    ve = (16*u*u+12*u+1)/(16*(1+u)**2)
    p = u**R(1,4)*(1+u)**R(-1,4)*(1+2*u/3)
    check_zero(S.simplify(u*S.diff(u*S.diff(p,u),u)/p-ve),'exact-point solution')
    norm = S.limit(p/u**R(1,4),u,0)
    growing = S.limit(p/u,u,S.oo)
    check_zero(norm-1,'origin unit coefficient')
    w0 = -2*growing
    check_zero(w0+R(4,3),'A1 finite exact-point Wronskian')
    result = {
        'status':'EXACT', 'scalar_Q':str(Q), 'frame_determinant':str(determinant),
        'observable':'lambda det(v_origin,v_end) = W_z(psi_origin,psi_end); no shifted or multiplied Q functions',
        'Frobenius_pair':'in x=u, chi=sqrt(u) psi, chi_pm ~ x^(1/2 +/- sigma), W_x(chi_plus,chi_minus)=-2 sigma; determinant=-2 sigma times the opposite-line connection coefficient',
        'exact_point_identity':'B_j=0 implies b=lambda, a=0 and c=lambda Lambda+V/lambda; both endpoint solutions and the observable equal A1 with the same normalization',
        'exact_point_example':{'solution':1,'t':'-1/3','X':'-2','Y':'0','also_solution2_t':'1/3','origin_exponents':['1/4','-1/4'],'infinity_exponents':['1','-1'],'scalar_W_at_lambda0':str(w0),'frame_factor':str(frame_factor)},
        'nonempty_subfamily':{'a':'f(t)/lambda','b':'lambda','c':'lambda Lambda+V/lambda','Q':str(qf),'f':'any meromorphic f(t), analytic and zero at applicable exact points; e.g. (t^2-1/9) g(t)','origin_sigma_squared':'-eps/24 + f(t)^2/lambda^2 (on the regular-origin chart eps<0)','inequivalence':'generic different f(t)^2 have different local monodromy exponents, hence cannot be identified by an allowed endpoint-preserving SL2 gauge'},
        'class_dimension':'infinite-dimensional over constants as a space of functions of t; the displayed subfamily alone supplies one arbitrary meromorphic function per component',
        'gauge_quotient_dimension':'infinite-dimensional over constants by the local-monodromy invariant; no finite generic-fibre variety dimension is claimed',
        'dimension_not_done':'number of independent meromorphic coefficient functions after all constraints and the gauge quotient; singular strata at exceptional fibres',
        'generic_fibre_coefficient_bookkeeping':'1890 coefficients before high-order scalar equations and endpoint conditions, at a generic nonresonant t; this is not the dimension of the solution variety',
        'stage1_scope':'algebraic exact-point gate passes. A global anchor matching chart, endpoint confluence, and finite-fibre moduli classification are NOT DONE.'}
    save(out,'stage1_connection.json',result)
    return result


def oscillator_connection(out, tamper):
    x, mu, K = S.symbols('x mu K', positive=True)
    s = S.symbols('s')
    aa = (1+mu)/2+K*K/4
    ff = S.Function('F')(s)
    xx = x*x
    # Verify the confluent hypergeometric reduction, including the energy scale.
    pref = x**(mu+R(1,2))*S.exp(-x*x/2)
    ansatz = pref*ff.subs(s,xx)
    reduced = S.simplify((S.diff(ansatz,x,2)-(x*x+K*K+(mu*mu-R(1,4))/x**2)*ansatz)/(4*pref))
    expected = (s*S.diff(ff,s,2)+(1+mu-s)*S.diff(ff,s)-aa*ff).subs(s,xx)
    check_zero(S.simplify(reduced-expected),'oscillator Kummer reduction')
    # Gamma recurrence in W=-2mu Gamma(mu)/Gamma(a).
    check_zero(S.combsimp(mu*S.gamma(mu)-S.gamma(1+mu)),'origin Wronskian gamma recurrence')
    spectral_scale = 1 if tamper != 'spectral' else 4
    z = S.symbols('z', positive=True)
    bern = S.bernoulli(3,(1+mu)/2)/(2*3)
    # log Gamma denominator: +B3(a)/(6z^2).
    k4 = S.factor(bern*16/spectral_scale**2)
    H = mu*(mu*mu-1)/3
    check_zero(k4-H,'derived first odd oscillator slot')
    odd_checks = []
    for n in range(1,9):
        coeff = S.expand((-1)**n*S.bernoulli(n+1,(1+mu)/2)/(n*(n+1)))
        odd = S.factor((coeff-coeff.subs(mu,-mu))/2)
        if n % 2:
            check_zero(odd, f'oscillator zero odd slot {n}')
        odd_checks.append({'inverse_z_power':n,'beta_odd_coefficient':str(odd)})
    result = {'status':'EXACT',
      'ODE':'psi_xx=[x^2+K^2+(mu^2-1/4)/x^2] psi; E=-K^2',
      'origin_solutions':'x^(1/2 +/- mu) exp(-x^2/2) M((1 +/- mu)/2+K^2/4,1 +/- mu,x^2), unit origin coefficients',
      'recessive_solution':'x^(1/2+mu) exp(-x^2/2) U((1+mu)/2+K^2/4,1+mu,x^2)',
      'infinity_normalization':'exp(-x^2/2) x^(-K^2/2-1/2), leading coefficient one, independent of mu',
      'canonical_W':'-2 Gamma(1+mu)/Gamma((1+mu)/2+K^2/4)',
      'normalization_scope':'The factor -2 is retained in the actual Wronskian. It cancels only in the beta-odd logarithm; no global M2 normalization or K(lambda,t) is inferred.',
      'odd_K0':'-mu log K + mu log 2 + (log Gamma(1+mu)-log Gamma(1-mu))/2',
      'odd_Kminus4':str(k4),'Stirling_parity_checks':odd_checks,
      'h0_singularity_qualification':'The logarithmic block has logarithmic singularities at positive integer mu; Gamma(1-mu) has poles, its reciprocal zeros. These are not the simple meromorphic poles of the h1 coefficient.',
      'global_M2_embedding':'NOT DONE'}
    save(out,'oscillator_connection.json',result)
    return result


def origin_response(out,tamper):
    mu, alpha, v, q, d, g1, ell = S.symbols('mu alpha v q d g1 ell')
    n = S.symbols('n', integer=True, positive=True)
    H = mu*(mu*mu-1)/3
    C = -S.sin(S.pi*mu)/S.pi*2**(2*alpha-1)*S.gamma(alpha+1)**2*S.gamma(alpha+1+mu)*S.gamma(alpha+1-mu)/S.gamma(2*alpha+2)
    # Schlaefli K integral and the changes r=s+t, y=s/(s+t) give this Mellin value.
    J = 2**(2*alpha-1)*S.gamma(alpha+1)**2*S.gamma(alpha+1+mu)*S.gamma(alpha+1-mu)/S.gamma(2*alpha+2)
    anchor_product = mu*(1-mu*mu)*S.pi/S.sin(S.pi*mu)
    anchor = S.cancel(-S.sin(S.pi*mu)/S.pi*anchor_product/3)
    check_zero(anchor-H,'Mellin anchor Gamma product')
    psi_sum = 2*S.digamma(mu)+S.pi*S.cot(S.pi*mu)+1/mu+2/(1-mu*mu)
    log_derivative = S.simplify(S.diff(C,alpha)/C)
    at_one = S.expand_func(log_derivative.subs(alpha,1))
    # Reflection psi(1-mu)=psi(mu)+pi cot(pi mu), with recurrence
    # psi(-mu)=psi(1-mu)+1/mu, accounts for sympy's possible base choices.
    at_one = at_one.subs(S.digamma(1-mu),S.digamma(mu)+S.pi*S.cot(S.pi*mu))
    at_one = at_one.subs(S.digamma(-mu),S.digamma(mu)+S.pi*S.cot(S.pi*mu)+1/mu)
    dalpha = 2*S.log(2)-R(5,3)+psi_sum
    check_zero(at_one-dalpha,'derivative of the source-derived Mellin coefficient')
    source_gain = 2 if tamper == 'source' else 1
    # Free local flow: alpha=1+v h, amplitude=1+g1 h,
    # mu(h)=mu+h(q mu+d/mu). No target R enters here.
    flow = H*(g1+source_gain*v*(dalpha-2*ell))+(q*mu+d/mu)*S.diff(H,mu)
    xflow = S.expand(flow.subs({v:-R(1,2),q:-R(1,8),d:0,g1:0}))
    yflow = S.expand(flow.subs({v:-R(1,2),q:-R(1,8),d:-R(3,8),g1:0}))
    target_bracket = ell-S.log(2)-S.digamma(mu)-S.pi*S.cot(S.pi*mu)/2-1/(2*mu)
    RX = mu*(11*mu*mu+7)/72
    RY = (11*mu**4-20*mu*mu+9)/(72*mu)
    check_zero(xflow-H*target_bracket-RX,'independent origin response X')
    check_zero(yflow-H*target_bracket-RY,'independent origin response Y')
    general_res = S.factor(v*H.subs(mu,n))
    target_res = S.factor(general_res.subs(v,-R(1,2)))
    zero_res = -d/3
    check_zero(target_res+n*(n*n-1)/6,'integer residue law')
    check_zero(zero_res.subs(d,-R(3,8))-R(1,8),'Y zero residue')
    # The dictionary derivative is derived independently before substitution.
    t,X,Y = S.symbols('t X Y')
    al = 2/(t-1)
    N = (t*t-25)/(t*t-1)
    deltas = {'X':-X/2+(N-1)/24,'Y':-Y/2}
    mu2 = {name:S.factor(4*(al+1)*delta+al*al) for name,delta in deltas.items()}
    rates = {name:S.factor(S.diff(value,t).subs(t,3)) for name,value in mu2.items()}
    check_zero(S.diff(al,t).subs(t,3)+R(1,2),'alpha slope')
    check_zero(rates['X']-X,'X mu squared slope')
    check_zero(rates['Y']-(Y-1),'Y mu squared slope')
    # Concrete local families have identical h0 determinants but different h1
    # singularities. They are witnesses for the local class, NOT global M2 members.
    variants = []
    for vv,dd in [(R(-1,2),R(-3,8)),(R(-1,4),R(-3,8)),(R(-1,2),R(0))]:
        variants.append({'alpha_prime':str(vv),'inverse_mu_velocity':str(dd),'positive_integer_residue':str(general_res.subs(v,vv)),'mu0_residue':str(zero_res.subs(d,dd))})
    result = {'status':'EXACT origin-ODE derivation, no BLZ formula used to compute the response',
      'ODE':'psi_xx=[K^2+(mu(h)^2-1/4)/x^2+g(h)x^(2alpha(h))] psi',
      'Green_function':'G_mu(x,x)=x I_mu(Kx) K_mu(Kx); odd G=-sin(pi mu)/pi x K_mu(Kx)^2',
      'Mellin_integral':str(J),'Mellin_convergence':'Re(alpha+1)>abs(Re(mu)); extension elsewhere is meromorphic continuation, not numerical integration through an origin divergence',
      'first_endpoint_coefficient':str(C),
      'spectral_slot':'g C(alpha,mu) K^(-2alpha-2); higher insertions occur at higher multiples of 2alpha+2',
      'anchor_coefficient':str(H),'alpha_derivative_log_coefficient_at1':str(dalpha),
      'general_first_response':str(flow),
      'dictionary_mu_squared':{key:str(value) for key,value in mu2.items()},
      'dictionary_mu_squared_slopes':{key:str(value) for key,value in rates.items()},
      'dictionary_slopes':{'alpha':'-1/2','mu_X':'-mu/8','mu_Y':'-(mu^2+3)/(8mu)','amplitude':'0 (unit power-law coefficient)'},
      'derived_rational_remainders':{'X':str(S.factor(xflow-H*target_bracket)),'Y':str(S.factor(yflow-H*target_bracket))},
      'comparison_residuals':{'X':'0','Y':'0'},
      'general_positive_integer_residue_n_ge2':str(general_res),
      'target_positive_integer_residue_n_ge2':str(target_res),
      'general_mu0_residue':str(zero_res),
      'target_zero_residues':{'X':'0','Y':'1/8'},
      'mu1_pole':'absent: H has a zero, neither rational index velocity is singular at 1',
      'residue_sample_values':{str(j):str(target_res.subs(n,j)) for j in range(2,9)},
      'local_same_h0_counterfamilies':variants,
      'filter_teeth':'In this local power-law class, the specified residues force alpha_prime=-1/2 and inverse_mu_velocity_Y=-3/8. They do not determine g_prime or the coefficient of mu in mu_prime.',
      'conditionality':'This proves the response for the stated radial ODE directly. Its identification with the cylindrical hierarchy is still the brief\'s imported conditional target, not a global oper proof.',
      'global_M2_gate2':'NOT DONE: these local families have not been embedded in the finite three-puncture M2 class. The seal\'s freedom prediction for global M2 is unresolved.'}
    save(out,'origin_response.json',result)
    return result


def numeric(out):
    import mpmath as mp
    with mp.workdps(45):
        def integral(a,m,log_weight=False):
            def integrand(x):
                if not x:
                    return mp.mpf(0)
                val=x**(2*a+1)*mp.besselk(m,x)**2
                return val*(2*mp.log(x) if log_weight else 1)
            return mp.quad(integrand,[0,mp.mpf('.1'),1,4,mp.inf])
        def closed(a,m):
            return 2**(2*a-1)*mp.gamma(a+1)**2*mp.gamma(a+1+m)*mp.gamma(a+1-m)/mp.gamma(2*a+2)
        rows=[]
        for av,mv in [('1','.4'),('1','.7'),('.8','.3'),('1.2','.6')]:
            a,m=mp.mpf(av),mp.mpf(mv)
            val=integral(a,m);ref=closed(a,m)
            err=abs(val-ref)/abs(ref)
            if err>mp.mpf('1e-32'):
                raise ArithmeticError('Bessel quadrature mismatch')
            rows.append({'alpha':av,'mu':mv,'integral':mp.nstr(val,35),'closed':mp.nstr(ref,35),'relative_error':mp.nstr(err,10)})
        a,m=mp.mpf(1),mp.mpf('.4')
        logint=integral(a,m,True)
        ref=mp.diff(lambda aa:closed(aa,m),a)
        err=abs(logint-ref)
        if err>mp.mpf('1e-32'):
            raise ArithmeticError('log source quadrature mismatch')
        rows.append({'alpha':'1','mu':'.4','insertion':'2 log x','integral':mp.nstr(logint,35),'closed_derivative':mp.nstr(ref,35),'absolute_error':mp.nstr(err,10)})
        # Pole controls evaluate the source-derived general coefficient, not the
        # target digamma formula. The difference in alpha is a finite difference
        # of Gamma functions, and mu is kept off resonance.
        def C(aa,mm):
            return -mp.sin(mp.pi*mm)/mp.pi*closed(aa,mm)
        eps=mp.mpf('1e-12')
        residues=[]
        for nv in [2,3,4,7]:
            mm=nv+eps
            vel=-(mm*mm+3)/(8*mm)
            value=-mp.diff(lambda aa:C(aa,mm),1)/2+vel*mp.diff(lambda xx:C(1,xx),mm)
            res=eps*value
            expected=-mp.mpf(nv)*(nv*nv-1)/6
            if abs(res-expected)>mp.mpf('1e-8'):
                raise ArithmeticError('source-derived integer residue')
            residues.append({'n':nv,'epsilon':str(eps),'scaled_response':mp.nstr(res,25),'limit':str(expected)})
        save(out,'numerical_checks.json',{'status':'NUMERICAL, mpmath corroboration only','precision':45,'convergent_Mellin_checks':rows,'meromorphic_continuation_residues':residues})


def main():
    p=UsageParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--numerical',action='store_true')
    p.add_argument('--tamper',choices=('frame','spectral','source'))
    args=p.parse_args()
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    try:
        scalar_and_frame(out,args.tamper)
        oscillator_connection(out,args.tamper)
        origin_response(out,args.tamper)
        if args.numerical:
            numeric(out)
    except (ArithmeticError,ValueError) as exc:
        save(out,'failure.json',{'exit':1,'error':str(exc),'tamper':args.tamper})
        print('FAIL',str(exc))
        return 1
    result={'exit':3,'stage1':'EXACT algebraic observable and exact-point identity; nonempty infinite-dimensional class and quotient as functions of t','stage2':'NOT DONE for global M2; EXACT independent local origin-response derivation reproduces the displayed target for the radial power-law family','stage3':'NOT DONE: global stage2 has not passed; all certified coefficient hold-outs untouched','unresolved_requirement':'A matched anchor chart embedding the derived radial endpoint family in the retained finite M2 spatial/spectral class has not been constructed. Neither emptiness nor a candidate is established.','seal_prediction':'global forced/free question unresolved; freedom proved only in the separate local endpoint family','scope':'exit3 is incomplete global construction, not consistency evidence, not an exclusion, and not a matrix-branch justification','fitting':'none; local flow rates imported from the independently stated dictionary, not from residues'}
    save(out,'summary.json',result)
    print(json.dumps(result,indent=2))
    return 3


if __name__=='__main__':
    sys.exit(main())
