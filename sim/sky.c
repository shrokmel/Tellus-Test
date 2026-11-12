/* sky.c: fast integrator for the contagious-gravity star system.
 *
 * Physics (units AU, years, solar masses):
 *   force between i and j:  G q_i q_j / r^2      (acceleration of i: G q_j (q_i/m_i) / r^2)
 *   pull factor f = q/m
 *   charge flow:  dq_i/dt = -lam * sum_j w(r_ij) mu_ij (f_i - f_j),  w = exp(-(r/d)^2),  mu = m_i m_j/(m_i+m_j)
 *   lam = 0 gives plain Newton.
 * Bodies flagged as test particles (tp = 1) feel all other non-test bodies and exchange charge with them, but exert
 * no gravity and do not interact with each other (for fast scans of many candidate comets).
 *
 * Integrator: Dormand-Prince 5(4), adaptive step, max-norm error control, steps land exactly on sample times.
 *
 * Input (text):  line 1: N G lam d years sample_days rtol
 *                then N lines: name m q x y z vx vy vz tp
 * Output (binary, float64): for each sample: t, then for each body x y z vx vy vz q.
 * Usage: sky input.txt output.bin
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int N; static double G, LAM, D, *M, *ATOLQ; static int *TP;

static void rhs(const double *y, double *dy) {
    const double *X = y, *V = y + 3 * N, *Q = y + 6 * N;
    double *DX = dy, *A = dy + 3 * N, *DQ = dy + 6 * N;
    memcpy(DX, V, sizeof(double) * 3 * N);
    memset(A, 0, sizeof(double) * 3 * N);
    memset(DQ, 0, sizeof(double) * N);
    double inv_d2 = (D > 0) ? 1.0 / (D * D) : 0;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            if (TP[i] && TP[j]) continue;
            double dx = X[3*j] - X[3*i], dy_ = X[3*j+1] - X[3*i+1], dz = X[3*j+2] - X[3*i+2];
            double r2 = dx*dx + dy_*dy_ + dz*dz, inv_r = 1.0 / sqrt(r2), inv_r3 = inv_r * inv_r * inv_r;
            double s = G * Q[i] * Q[j] * inv_r3;          /* force magnitude / r */
            if (!TP[j]) { double ai = s / M[i]; A[3*i] += ai*dx; A[3*i+1] += ai*dy_; A[3*i+2] += ai*dz; }
            if (!TP[i]) { double aj = s / M[j]; A[3*j] -= aj*dx; A[3*j+1] -= aj*dy_; A[3*j+2] -= aj*dz; }
            if (LAM > 0) {
                double w = exp(-r2 * inv_d2);
                if (w > 1e-300) {
                    double mu = M[i] * M[j] / (M[i] + M[j]);
                    double flow = LAM * w * mu * (Q[i] / M[i] - Q[j] / M[j]);
                    DQ[i] -= flow; DQ[j] += flow;
                }
            }
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: sky input.txt output.bin\n"); return 1; }
    FILE *fi = fopen(argv[1], "r"); if (!fi) { perror("input"); return 1; }
    double years, sample_days, rtol;
    if (fscanf(fi, "%d %lf %lf %lf %lf %lf %lf", &N, &G, &LAM, &D, &years, &sample_days, &rtol) != 7) return 2;
    int S = 7 * N;
    M = malloc(sizeof(double) * N); TP = malloc(sizeof(int) * N); ATOLQ = malloc(sizeof(double) * N);
    double *y = malloc(sizeof(double) * S);
    for (int i = 0; i < N; i++) {
        char name[64]; double m, q, x, yy, z, vx, vy, vz; int tp;
        if (fscanf(fi, "%63s %lf %lf %lf %lf %lf %lf %lf %lf %d", name, &m, &q, &x, &yy, &z, &vx, &vy, &vz, &tp) != 10) return 3;
        M[i] = m; TP[i] = tp; ATOLQ[i] = 1e-12 * m;
        y[3*i] = x; y[3*i+1] = yy; y[3*i+2] = z; y[3*N+3*i] = vx; y[3*N+3*i+1] = vy; y[3*N+3*i+2] = vz; y[6*N+i] = q;
    }
    fclose(fi);
    FILE *fo = fopen(argv[2], "wb"); if (!fo) { perror("output"); return 1; }

    static const double c2=1./5, c3=3./10, c4=4./5, c5=8./9;
    static const double a21=1./5, a31=3./40, a32=9./40, a41=44./45, a42=-56./15, a43=32./9,
        a51=19372./6561, a52=-25360./2187, a53=64448./6561, a54=-212./729,
        a61=9017./3168, a62=-355./33, a63=46732./5247, a64=49./176, a65=-5103./18656,
        b1=35./384, b3=500./1113, b4=125./192, b5=-2187./6784, b6=11./84,
        e1=71./57600, e3=-71./16695, e4=71./1920, e5=-17253./339200, e6=22./525, e7=-1./40;
    (void)c2; (void)c3; (void)c4; (void)c5;
    double *k1=malloc(sizeof(double)*S), *k2=malloc(sizeof(double)*S), *k3=malloc(sizeof(double)*S),
           *k4=malloc(sizeof(double)*S), *k5=malloc(sizeof(double)*S), *k6=malloc(sizeof(double)*S),
           *k7=malloc(sizeof(double)*S), *yt=malloc(sizeof(double)*S), *yn=malloc(sizeof(double)*S);
    double dt_s = sample_days / 365.25, t = 0, h = 1e-4;
    long n_samples = (long)floor(years / dt_s + 1e-9) + 1, k_sample = 0, steps = 0, rejects = 0;
    rhs(y, k1);
    fwrite(&t, sizeof(double), 1, fo); fwrite(y, sizeof(double), S, fo); k_sample = 1;
    while (k_sample < n_samples) {
        double t_next = k_sample * dt_s, hh = h; int hit = 0;
        if (t + hh >= t_next) { hh = t_next - t; hit = 1; }
        for (int i=0;i<S;i++) yt[i] = y[i] + hh*(a21*k1[i]); rhs(yt, k2);
        for (int i=0;i<S;i++) yt[i] = y[i] + hh*(a31*k1[i] + a32*k2[i]); rhs(yt, k3);
        for (int i=0;i<S;i++) yt[i] = y[i] + hh*(a41*k1[i] + a42*k2[i] + a43*k3[i]); rhs(yt, k4);
        for (int i=0;i<S;i++) yt[i] = y[i] + hh*(a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]); rhs(yt, k5);
        for (int i=0;i<S;i++) yt[i] = y[i] + hh*(a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]); rhs(yt, k6);
        for (int i=0;i<S;i++) yn[i] = y[i] + hh*(b1*k1[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i]); rhs(yn, k7);
        double err = 0;
        for (int i=0;i<S;i++) {
            double ei = hh*(e1*k1[i] + e3*k3[i] + e4*k4[i] + e5*k5[i] + e6*k6[i] + e7*k7[i]);
            double atol = (i >= 6*N) ? ATOLQ[i-6*N] : 1e-14;
            double sc = atol + rtol * fmax(fabs(y[i]), fabs(yn[i]));
            double r = fabs(ei) / sc; if (r > err) err = r;
        }
        double fac = (err == 0) ? 5.0 : fmin(5.0, fmax(0.2, 0.9 * pow(err, -0.2)));
        if (err <= 1.0) {
            t = hit ? t_next : t + hh; memcpy(y, yn, sizeof(double)*S); memcpy(k1, k7, sizeof(double)*S); steps++;
            if (hit) { fwrite(&t, sizeof(double), 1, fo); fwrite(y, sizeof(double), S, fo); k_sample++; }
            if (!hit || hh * fac < h) h = hh * fac;   /* do not let a short landing step shrink h */
        } else { h = hh * fac; rejects++; }
        if (h < 1e-14) { fprintf(stderr, "step size underflow at t=%g\n", t); return 4; }
    }
    fclose(fo);
    fprintf(stderr, "done: %ld samples, %ld steps, %ld rejected\n", n_samples, steps, rejects);
    return 0;
}
