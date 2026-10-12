"""WUW1 post-hoc (labelled): why TAMPER 2 (W_U from the Schur polynomial with exponent a = 21/10) did not break P1/P3.  Prints eigen(S_3(a)) for symbolic
a and checks the field-level / zero-mode properties of the tampered field."""
import os, sys
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'wuw1.py')).read()
pre = src[:src.index('# ================= (1) the weight-3 objects')].replace("open(os.path.join(HERE, 'run_wuw1.log'), 'w')", "open(os.path.join(HERE, 'run_diag_tamper2.log'), 'w')")
exec(pre)
for a in (F(2), F(21, 10), F(3), F(1)):
    Wa = schur(a, U, 3); e = sp.factor(eigen(Wa, internal=False))
    say(f'a = {a}: S_3 field {Wa}; eigen = {e}; field-level Res e^(-U) W = {field_screen(Wa, F(-1), U)}')
say('reading: the cylinder vacuum eigenvalue of S_3(a d^k U/k!) is a^3/6 [P^3 - P/4]: the a-dependence is an overall factor, so the sealed exponent tamper is')
say('VACUOUS for the eigenvalue tests (P1, P3); it is NOT vacuous at field level (screening).')
