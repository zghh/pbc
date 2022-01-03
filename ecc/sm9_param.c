#include <stdarg.h>
#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <stdlib.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_fp.h"
#include "pbc_fieldquadratic.h"
#include "pbc_param.h"
#include "pbc_pairing.h"
#include "pbc_poly.h"
#include "pbc_curve.h"
#include "pbc_memory.h"
#include "pbc_sm9_param.h"
#include "ecc/param.h"

struct sm9_param_s
{
    mpz_t q;    // Curve defined over F_q.
    mpz_t r;    // The order of the curve.
    mpz_t b;    // E: y^2 = x^3 + b
    mpz_t t;    // Curve parameter
    mpz_t beta; // beta is a quadratic nonresidue in Fq, we use F_q^2 = F_q[sqrt(beta)]
    mpz_t alpha0, alpha1;
    // the polynomial x^6 + alpha0 + alpha1 sqrt(beta)
    // is irreducible over F_q^2[x], so
    // we can extend F_q^2 to F_q^12 using the
    // sixth root of -(alpha0 + alpha1 sqrt(beta))
};
typedef struct sm9_param_s sm9_param_t[1];
typedef struct sm9_param_s *sm9_param_ptr;

// TODO: we never use phikonr so don't bother computing it,
// but one day other routines might need it
struct sm9_pairing_data_s
{
    field_t Fq, Fq2, Fq2x, Fq12;
    field_t Eq, Etwist;
    element_t negalpha;
    element_t negalphainv;
    mpz_t tateexp;
    mpz_t t, q;

    // for tate exponentiation speedup:
    // x^{q^k} for various k
    element_t xpowq2, xpowq6, xpowq8;
};
typedef struct sm9_pairing_data_s sm9_pairing_data_t[1];
typedef struct sm9_pairing_data_s *sm9_pairing_data_ptr;

static void sm9_clear(void *data)
{
    sm9_param_ptr fp = data;
    mpz_clear(fp->q);
    mpz_clear(fp->r);
    mpz_clear(fp->b);
    mpz_clear(fp->t);
    mpz_clear(fp->beta);
    mpz_clear(fp->alpha0);
    mpz_clear(fp->alpha1);
    pbc_free(data);
}

static void sm9_out_str(FILE *stream, void *data)
{
    sm9_param_ptr p = data;
    param_out_type(stream, "sm9");
    param_out_mpz(stream, "q", p->q);
    param_out_mpz(stream, "r", p->r);
    param_out_mpz(stream, "b", p->b);
    param_out_mpz(stream, "b", p->t);
    param_out_mpz(stream, "beta", p->beta);
    param_out_mpz(stream, "alpha0", p->alpha0);
    param_out_mpz(stream, "alpha1", p->alpha1);
}

static void sm9_finalpow(element_t out) {}

static void sm9_tangent(element_t out, element_t a, element_t c, sm9_pairing_data_ptr p)
{
    element_ptr ax, ay, cx, cy, cof0x, cof0y, out0, out3, out5;
    element_t lamda, cof0, cof3, cof5, temp;

    ax = element_x(a);
    ay = element_y(a);
    cx = element_x(c);
    cy = element_y(c);

    element_init_same_as(lamda, ax);
    element_square(lamda, ax);
    element_mul_si(lamda, lamda, 3);
    element_init_same_as(temp, ay);
    element_mul_si(temp, ay, 2);
    element_div(lamda, lamda, temp);

    element_init_same_as(cof3, ay);
    element_set(cof3, ay);
    element_mul(temp, lamda, ax);
    element_sub(cof3, ay, temp);

    element_init_same_as(cof5, lamda);
    element_mul_zn(cof5, lamda, cx);

    element_init_same_as(cof0, ax);
    cof0x = element_x(cof0);
    cof0y = element_y(cof0);
    element_neg(cof0x, cy);
    element_set0(cof0y);

    out0 = out->field->item(out, 0);
    out3 = out->field->item(out, 3);
    out5 = out->field->item(out, 5);
    element_set(out0, cof0);
    element_mul(out3, cof3, p->negalphainv);
    element_mul(out5, cof5, p->negalphainv);

    element_clear(lamda);
    element_clear(cof0);
    element_clear(cof3);
    element_clear(cof5);
    element_clear(temp);
}

static void sm9_line(element_t out, element_t a, element_t b, element_t c, sm9_pairing_data_ptr p)
{
    element_ptr ax, ay, bx, by, cx, cy, cof0x, cof0y, out0, out3, out5;
    element_t lamda, cof0, cof3, cof5, temp;

    ax = element_x(a);
    ay = element_y(a);
    bx = element_x(b);
    by = element_y(b);
    cx = element_x(c);
    cy = element_y(c);

    element_init_same_as(lamda, ay);
    element_sub(lamda, ay, by);
    element_init_same_as(temp, ax);
    element_sub(temp, ax, bx);
    element_div(lamda, lamda, temp);

    element_init_same_as(cof3, by);
    element_set(cof3, by);
    element_mul(temp, lamda, bx);
    element_sub(cof3, by, temp);

    element_init_same_as(cof5, lamda);
    element_mul_zn(cof5, lamda, cx);

    element_init_same_as(cof0, ax);
    cof0x = element_x(cof0);
    cof0y = element_y(cof0);
    element_neg(cof0x, cy);
    element_set0(cof0y);

    out0 = out->field->item(out, 0);
    out3 = out->field->item(out, 3);
    out5 = out->field->item(out, 5);
    element_set(out0, cof0);
    element_mul(out3, cof3, p->negalphainv);
    element_mul(out5, cof5, p->negalphainv);

    element_clear(lamda);
    element_clear(cof0);
    element_clear(cof3);
    element_clear(cof5);
    element_clear(temp);
}

static void sm9_frobenius(element_t out, element_t point, sm9_pairing_data_ptr p)
{
    element_t x, r, w, px, py;
    mpz_t temp;

    element_init_same_as(px, element_x(point));
    element_set(px, element_x(point));
    element_init_same_as(py, element_y(point));
    element_set(py, element_y(point));

    element_init(x, p->Fq2);
    element_set0(element_x(x));
    element_set1(element_y(x));
    mpz_init(temp);
    mpz_sub_ui(temp, p->q, 1);
    mpz_divexact_ui(temp, temp, 6);
    element_pow_mpz(x, x, temp);

    element_init_same_as(r, x);
    element_invert(r, x);
    element_init_same_as(w, r);
    element_square(w, r);

    element_neg(element_y(px), element_y(px));
    element_mul(px, px, w);

    element_neg(element_y(py), element_y(py));
    element_mul(py, py, w);
    element_mul(py, py, r);

    element_random(out);
    element_set(element_x(out), px);
    element_set(element_y(out), py);

    element_clear(x);
    element_clear(r);
    element_clear(w);
    element_clear(px);
    element_clear(py);
    mpz_clear(temp);
}

static void sm9_pairing(element_ptr out, element_ptr in1, element_ptr in2,
                        pairing_t pairing)
{
    sm9_pairing_data_ptr p = pairing->data;
    mpz_t a, temp;
    int i;
    element_t t, f, line, q1, q2, negq2;

    element_init(line, p->Fq12);

    mpz_init(a);
    mpz_mul_ui(a, p->t, 6);
    mpz_add_ui(a, a, 2);
    element_init_same_as(t, in2);
    element_set(t, in2);
    element_init(f, p->Fq12);
    element_set1(f);

    for (i = mpz_sizeinbase(a, 2) - 2; i >= 0; i--)
    {
        element_square(f, f);
        sm9_tangent(line, t, in1, p);
        element_mul(f, f, line);
        element_add(t, t, t);
        if (mpz_tstbit(a, i))
        {
            sm9_line(line, t, in2, in1, p);
            element_mul(f, f, line);
            element_add(t, t, in2);
        }
    }

    element_init_same_as(q1, in2);
    sm9_frobenius(q1, in2, p);
    element_init_same_as(q2, q1);
    sm9_frobenius(q2, q1, p);
    sm9_line(line, t, q1, in1, p);
    element_mul(f, f, line);
    element_add(t, t, q1);
    element_init_same_as(negq2, q2);
    element_neg(negq2, q2);
    sm9_line(line, t, negq2, in1, p);
    element_mul(f, f, line);
    element_sub(t, t, q2);
    mpz_init(temp);
    mpz_pow_ui(temp, p->q, 12);
    mpz_sub_ui(temp, temp, 1);
    mpz_divexact(temp, temp, pairing->r);
    element_pow_mpz(out, f, temp);

    mpz_clear(a);
    mpz_clear(temp);
    element_clear(t);
    element_clear(f);
    element_clear(line);
    element_clear(q1);
    element_clear(q2);
    element_clear(negq2);
}

static void sm9_pairing_clear(pairing_t pairing)
{
    field_clear(pairing->GT);
    sm9_pairing_data_ptr p = pairing->data;
    element_clear(p->negalpha);
    element_clear(p->negalphainv);
    mpz_clear(p->tateexp);
    element_clear(p->xpowq2);
    element_clear(p->xpowq6);
    element_clear(p->xpowq8);
    field_clear(p->Etwist);
    field_clear(p->Eq);

    field_clear(p->Fq12);
    field_clear(p->Fq2x);
    field_clear(p->Fq2);
    field_clear(p->Fq);
    mpz_clear(p->t);
    mpz_clear(p->q);
    pbc_free(p);

    mpz_clear(pairing->r);
    field_clear(pairing->Zr);
}

static void sm9_init_pairing(pairing_t pairing, void *data)
{
    sm9_param_ptr param = data;
    sm9_pairing_data_ptr p;
    element_t irred;
    element_t e0, e1, e2;
    p = pairing->data = pbc_malloc(sizeof(sm9_pairing_data_t));
    mpz_init(pairing->r);
    mpz_set(pairing->r, param->r);
    field_init_fp(pairing->Zr, pairing->r);
    field_init_fp(p->Fq, param->q);
    p->Fq->nqr = pbc_malloc(sizeof(element_t));
    element_init(p->Fq->nqr, p->Fq);
    element_set_mpz(p->Fq->nqr, param->beta);
    field_init_quadratic(p->Fq2, p->Fq);
    field_init_poly(p->Fq2x, p->Fq2);
    element_init(irred, p->Fq2x);
    // Call poly_set_coeff1() first so we can use element_item() for the other
    // coefficients.
    poly_set_coeff1(irred, 6);

    element_init(p->negalpha, p->Fq2);
    element_init(p->negalphainv, p->Fq2);
    element_set_mpz(element_x(p->negalpha), param->alpha0);
    element_set_mpz(element_y(p->negalpha), param->alpha1);

    element_set(element_item(irred, 0), p->negalpha);
    field_init_polymod(p->Fq12, irred);
    element_neg(p->negalpha, p->negalpha);
    element_invert(p->negalphainv, p->negalpha);
    element_clear(irred);

    element_init(e0, p->Fq);
    element_init(e1, p->Fq);
    element_init(e2, p->Fq2);

    // Initialize the curve Y^2 = X^3 + b.
    element_set_mpz(e1, param->b);
    field_init_curve_ab(p->Eq, e0, e1, pairing->r, NULL);

    // Initialize the curve Y^2 = X^3 - alpha0 b - alpha1 sqrt(beta) b.
    element_set_mpz(e0, param->alpha0);
    element_neg(e0, e0);
    element_mul(element_x(e2), e0, e1);
    element_set_mpz(e0, param->alpha1);
    element_neg(e0, e0);
    element_mul(element_y(e2), e0, e1);
    element_clear(e0);
    element_init(e0, p->Fq2);
    field_init_curve_ab(p->Etwist, e0, e2, pairing->r, NULL);
    element_clear(e0);
    element_clear(e1);
    element_clear(e2);

    mpz_t ndonr;
    mpz_init(ndonr);
    // ndonr temporarily holds the trace.
    mpz_sub(ndonr, param->q, param->r);
    mpz_add_ui(ndonr, ndonr, 1);
    // TODO: We can use a smaller quotient_cmp, but I have to figure out
    // BN curves again.
    pbc_mpz_curve_order_extn(ndonr, param->q, ndonr, 12);
    mpz_divexact(ndonr, ndonr, param->r);
    mpz_divexact(ndonr, ndonr, param->r);
    field_curve_set_quotient_cmp(p->Etwist, ndonr);
    mpz_clear(ndonr);

    pairing->G1 = p->Eq;
    pairing->G2 = p->Etwist;
    pairing_GT_init(pairing, p->Fq12);
    pairing->finalpow = sm9_finalpow;
    pairing->map = sm9_pairing;
    pairing->clear_func = sm9_pairing_clear;

    mpz_init(p->tateexp);
    /* unoptimized tate exponent
    mpz_pow_ui(p->tateexp, param->q, 12);
    mpz_sub_ui(p->tateexp, p->tateexp, 1);
    mpz_divexact(p->tateexp, p->tateexp, param->r);
    */
    mpz_ptr z = p->tateexp;
    mpz_mul(z, param->q, param->q);
    mpz_sub_ui(z, z, 1);
    mpz_mul(z, z, param->q);
    mpz_mul(z, z, param->q);
    mpz_add_ui(z, z, 1);
    mpz_divexact(z, z, param->r);

    element_init(p->xpowq2, p->Fq2);
    element_init(p->xpowq6, p->Fq2);
    element_init(p->xpowq8, p->Fq2);
    element_t xpowq;
    element_init(xpowq, p->Fq12);

    // there are smarter ways since we know q = 1 mod 6
    // and that x^6 = -alpha
    // but this is fast enough
    element_set1(element_item(xpowq, 1));
    element_pow_mpz(xpowq, xpowq, param->q);
    element_pow_mpz(xpowq, xpowq, param->q);
    element_set(p->xpowq2, element_item(xpowq, 1));

    element_pow_mpz(xpowq, xpowq, param->q);
    element_pow_mpz(xpowq, xpowq, param->q);
    element_pow_mpz(xpowq, xpowq, param->q);
    element_pow_mpz(xpowq, xpowq, param->q);
    element_set(p->xpowq6, element_item(xpowq, 1));

    element_pow_mpz(xpowq, xpowq, param->q);
    element_pow_mpz(xpowq, xpowq, param->q);
    element_set(p->xpowq8, element_item(xpowq, 1));

    element_clear(xpowq);

    mpz_init(p->t);
    mpz_set(p->t, param->t);
    mpz_init(p->q);
    mpz_set(p->q, param->q);
}

static void sm9_init(pbc_param_ptr p)
{
    static pbc_param_interface_t interface = {{
        sm9_clear,
        sm9_init_pairing,
        sm9_out_str,
    }};
    p->api = interface;
    sm9_param_ptr fp = p->data = pbc_malloc(sizeof(*fp));
    mpz_init(fp->q);
    mpz_init(fp->r);
    mpz_init(fp->b);
    mpz_init(fp->t);
    mpz_init(fp->beta);
    mpz_init(fp->alpha0);
    mpz_init(fp->alpha1);
}

// Public interface:

int pbc_param_init_sm9(pbc_param_ptr par, struct symtab_s *tab)
{
    sm9_init(par);
    sm9_param_ptr p = par->data;

    int err = 0;
    err += lookup_mpz(p->q, tab, "q");
    err += lookup_mpz(p->r, tab, "r");
    err += lookup_mpz(p->b, tab, "b");
    err += lookup_mpz(p->t, tab, "t");
    err += lookup_mpz(p->beta, tab, "beta");
    err += lookup_mpz(p->alpha0, tab, "alpha0");
    err += lookup_mpz(p->alpha1, tab, "alpha1");
    return err;
}
