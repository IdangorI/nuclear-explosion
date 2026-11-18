#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

static int W = 758, H = 630;
static float CX = 379.f, GYP = 557.5f, SC = 200.f, ROLL = 0.f;
static const float TILT = 0.055f;

static inline float clampf(float x, float a, float b) { return x < a ? a : (x > b ? b : x); }
static inline float ss(float a, float b, float x) { float t = clampf((x - a) / (b - a), 0, 1); return t * t * (3 - 2 * t); }
static inline float mixf(float a, float b, float t) { return a + (b - a) * t; }
static inline float smin(float a, float b, float k) { float h = clampf(0.5f + 0.5f * (b - a) / k, 0, 1); return mixf(b, a, h) - k * h * (1 - h); }
static inline float smax(float a, float b, float k) { return -smin(-a, -b, k); }

static inline unsigned hsh(int x, int y, int z) {
    unsigned h = (unsigned)x * 0x8da6b343u ^ (unsigned)y * 0xd8163841u ^ (unsigned)z * 0xcb1ab31fu;
    h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16; return h;
}
static inline float grad(unsigned h, float x, float y, float z) {
    int hh = h & 15; float u = hh < 8 ? x : y; float v = hh < 4 ? y : ((hh == 12 || hh == 14) ? x : z);
    return ((hh & 1) ? -u : u) + ((hh & 2) ? -v : v);
}
static inline float fd(float t) { return t * t * t * (t * (t * 6 - 15) + 10); }
static float pn(float x, float y, float z) {
    int X = (int)floorf(x), Y = (int)floorf(y), Z = (int)floorf(z);
    x -= X; y -= Y; z -= Z;
    float u = fd(x), v = fd(y), w = fd(z);
    float a = mixf(grad(hsh(X, Y, Z), x, y, z), grad(hsh(X + 1, Y, Z), x - 1, y, z), u);
    float b = mixf(grad(hsh(X, Y + 1, Z), x, y - 1, z), grad(hsh(X + 1, Y + 1, Z), x - 1, y - 1, z), u);
    float c = mixf(grad(hsh(X, Y, Z + 1), x, y, z - 1), grad(hsh(X + 1, Y, Z + 1), x - 1, y, z - 1), u);
    float d = mixf(grad(hsh(X, Y + 1, Z + 1), x, y - 1, z - 1), grad(hsh(X + 1, Y + 1, Z + 1), x - 1, y - 1, z - 1), u);
    return 0.7f * mixf(mixf(a, b, v), mixf(c, d, v), w);
}
static float worley(float x, float y, float z) {
    int X = (int)floorf(x), Y = (int)floorf(y), Z = (int)floorf(z);
    float md = 9;
    for (int dz = -1; dz <= 1; dz++) for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
        unsigned h = hsh(X + dx, Y + dy, Z + dz);
        float fx = X + dx + 0.15f + 0.7f * (h & 1023) / 1023.f - x;
        float fy = Y + dy + 0.15f + 0.7f * ((h >> 10) & 1023) / 1023.f - y;
        float fz = Z + dz + 0.15f + 0.7f * ((h >> 20) & 1023) / 1023.f - z;
        float d = fx * fx + fy * fy + fz * fz; if (d < md) md = d;
    }
    return sqrtf(md);
}
static float billow(float x, float y, float z, int oct) {
    float s = 0, a = 0.5f, n = 0;
    for (int i = 0; i < oct; i++) { s += a * fabsf(pn(x + i * 17.1f, y - i * 9.3f, z + i * 5.7f)); n += a; x *= 2.03f; y *= 2.03f; z *= 2.03f; a *= 0.5f; }
    return s / n;
}
static float fbm(float x, float y, float z, int oct) {
    float s = 0, a = 0.5f, n = 0;
    for (int i = 0; i < oct; i++) { s += a * pn(x + i * 31.7f, y + i * 3.1f, z - i * 11.9f); n += a; x *= 2.07f; y *= 2.07f; z *= 2.07f; a *= 0.5f; }
    return s / n;
}

typedef struct {
    float t, yc, ax, ay, m, Rt, rt, tyo, hy0, hr, sc, A, Tc, Ts, sootb, sootw, roll;
    float hs, ytop, rs0, rsb, stemFire, rb, hb, wr, baseFire, base_on;
    float Ik, Ip, amb, fade;
} Par;
static Par P;

static void setpar(float t) {
    P.t = t;
    float rf = 0.56f * (1 - expf(-(t + 0.012f) / 0.10f));
    float k = rf / 0.56f;
    P.yc = 0.26f + 1.18f * (1 - expf(-fmaxf(t - 0.25f, 0) / 1.15f)) + 0.08f * ss(0.0f, 0.4f, t);
    P.ax = (0.56f + 0.30f * ss(0.6f, 3.4f, t)) * k;
    P.ay = (0.56f - 0.17f * ss(0.6f, 2.6f, t)) * k;
    P.m = ss(0.75f, 2.5f, t);
    P.Rt = P.ax * 0.60f;
    P.rt = P.ay * 0.70f * (0.6f + 0.4f * P.m);
    P.tyo = -P.ay * 0.35f;
    P.hy0 = -P.ay * 1.45f;
    P.hr = P.ay * 0.95f * P.m;
    P.sc = fmaxf(P.ax, 0.2f) / 0.86f;
    P.A = 0.28f + 0.72f * ss(0.08f, 1.4f, t);
    P.Tc = 1.8f + 0.75f * expf(-t / 0.8f);
    P.Ts = 1.08f + 1.2f * expf(-t / 0.6f);
    P.sootb = 0.1f * ss(0.6f, 3.6f, t) - 0.28f;
    P.sootw = ss(0.25f, 1.3f, t);
    P.roll = 0.75f * ss(0.75f, 3.6f, t);
    P.hs = 0.04f + 1.7f * powf(ss(0.35f, 2.2f, t), 1.15f);
    P.ytop = fminf(P.hs, P.yc - P.ay * 0.45f);
    P.rs0 = 0.10f + 0.09f * ss(0.8f, 3.2f, t);
    P.rsb = 0.42f * ss(0.25f, 1.6f, t);
    P.stemFire = 1.5f * expf(-fmaxf(t - 0.7f, 0) / 0.9f) * ss(0.35f, 0.8f, t);
    P.base_on = t > 0.05f;
    P.rb = 0.25f + 1.28f * (1 - expf(-fmaxf(t - 0.05f, 0) / 0.5f));
    P.hb = 0.045f + 0.25f * ss(0.06f, 2.8f, t);
    P.wr = 0.12f + 0.20f * ss(0.15f, 2.6f, t);
    P.baseFire = 1.6f * expf(-t / 0.3f);
    P.fade = 1 - ss(3.62f, 3.97f, t);
    P.Ik = 1.0f * P.fade;
    P.Ip = (1.4f + 2.5f * expf(-t / 0.5f)) * P.fade;
    P.amb = 1.0f * P.fade;
}

static float bound(float x, float y, float z) {
    float qy = y - P.yc, rho = sqrtf(x * x + z * z);
    float ex = rho / P.ax, ey = qy / P.ay;
    float b = (sqrtf(ex * ex + ey * ey) - 1) * fminf(P.ax, P.ay) - 0.30f * P.sc;
    float tr = rho - P.Rt, ty = qy - P.tyo;
    b = fminf(b, sqrtf(tr * tr + ty * ty) - P.rt + (1 - P.m) * 0.8f - 0.30f * P.sc);
    if (y < P.ytop + 0.3f) {
        float rs = P.rs0 + P.rsb * expf(-y / 0.16f) + 0.2f;
        b = fminf(b, rho - rs - 0.32f);
    }
    if (P.base_on) {
        float hb = P.hb * 4.0f + 0.12f;
        float d = fmaxf(rho - P.rb * 1.25f - 0.1f, y - hb);
        b = fminf(b, d);
    }
    return b;
}

typedef struct { float d, T, dust, soot; } Smp;

static Smp field(float x, float y, float z, int detail) {
    Smp r = {0, 0, 0, 0};
    float t = P.t;
    float rho = sqrtf(x * x + z * z) + 1e-6f;
    float cth = x / rho, sth = z / rho;
    float dsum = 0, Tsum = 0, dustsum = 0, sootsum = 0;
    int oct = detail == 0 ? 1 : (detail == 2 ? 2 : 4);

    {
        float qy = y - P.yc;
        float ex = rho / P.ax, ey = qy / P.ay;
        float sd = (sqrtf(ex * ex + ey * ey) - 1) * fminf(P.ax, P.ay);
        float tr = rho - P.Rt, ty = qy - P.tyo;
        float sdt = sqrtf(tr * tr + ty * ty) - P.rt + (1 - P.m) * 0.8f;
        sd = smin(sd, sdt, 0.15f);
        float hy = qy - P.hy0;
        sd = smax(sd, -(sqrtf(rho * rho + hy * hy) - P.hr), 0.12f);
        sd = fmaxf(sd, -y);
        if (sd < 0.30f * P.sc) {

            float vx = rho - P.Rt, vy = qy - P.tyo;
            float wgt = expf(-(vx * vx + vy * vy) / (0.42f * 0.42f * P.sc * P.sc));
            float ang = P.roll * wgt;
            float ca = cosf(ang), sa = sinf(ang);
            float rx = vx * ca - vy * sa + P.Rt, ry = vx * sa + vy * ca + P.tyo;
            float nx = rx * cth / P.sc, ny = ry / P.sc, nz = rx * sth / P.sc;

            float wx = pn(nx * 1.6f + 0.35f * t, ny * 1.6f, nz * 1.6f + 1.3f);
            float wy = pn(nx * 1.6f + 7.1f, ny * 1.6f - 0.35f * t, nz * 1.6f);
            float wz = pn(nx * 1.6f, ny * 1.6f + 3.3f, nz * 1.6f + 0.35f * t);
            nx += 0.16f * wx; ny += 0.16f * wy; nz += 0.16f * wz;
            float lump = 1 - worley(nx * 2.7f, ny * 2.7f, nz * 2.7f);
            float lump2 = detail != 1 ? 0.45f : 1 - worley(nx * 7.3f + 3.f, ny * 7.3f, nz * 7.3f - 2.f);
            float bil = billow(nx * 9.f, ny * 9.f, nz * 9.f, oct + (detail == 1));
            float disp = P.A * P.sc * (0.20f * (lump - 0.42f) + 0.115f * (lump2 - 0.42f) + 0.11f * (bil - 0.30f));
            float sdf = sd - disp;
            float d = clampf(-sdf / 0.02f, 0, 1);
            if (d > 0) {
                dsum += d;
                if (detail == 1) {
                    float depth = -sd + 0.6f * disp;
                    float T = P.Ts + (P.Tc - P.Ts) * ss(0.0f, 0.28f, depth);
                    T *= 0.62f + 0.9f * ss(-0.6f, 1.0f, qy / P.ay) * (1 - 0.35f * ss(0.4f, 1.0f, rho / (P.ax + 0.01f)));
                    T *= 0.72f + 0.5f * ss(-0.035f, 0.045f, disp / P.sc);
                    float crev = 1 - ss(-0.03f, 0.03f, disp / P.sc);
                    float so = ss(0.45f, 0.8f, 0.45f + 0.5f * fbm(nx * 3.2f + 5, ny * 3.2f, nz * 3.2f, 2) + 0.45f * crev + P.sootb) * P.sootw;
                    so *= 1 - 0.85f * ss(-0.35f, 0.6f, qy / P.ay);
                    T *= 1 - 0.75f * so;
                    T *= 0.86f + 0.28f * (0.5f + 0.5f * pn(nx * 4.5f + t * 2.2f, ny * 4.5f - t * 1.6f, nz * 4.5f));
                    Tsum += d * T; sootsum += d * so;
                }
            }
        }
    }

    if (P.hs > 0.06f && y < P.ytop + 0.35f) {
        float ytop = P.ytop;
        float ox = 0.09f * pn(1.3f, y * 1.6f - 0.45f * t, 2.1f), oz = 0.09f * pn(5.3f, y * 1.6f - 0.45f * t, 7.7f);
        float sx = x - ox, sz = z - oz;
        float rr = sqrtf(sx * sx + sz * sz);
        float rs = P.rs0 + P.rsb * expf(-y / 0.16f) + 0.20f * expf(-(ytop - y) / 0.2f);
        float sd = smax(rr - rs, y - ytop, 0.12f);
        sd = fmaxf(sd, -y);
        if (sd < 0.3f) {
            float an = 0.9f * t, ca = cosf(an), sa = sinf(an);
            float nx = sx * ca - sz * sa, nz = sx * sa + sz * ca, ny = y - 0.5f * t;
            float wx = pn(nx * 2.2f, ny * 2.2f, nz * 2.2f + 4.1f), wz = pn(nx * 2.2f + 9.f, ny * 2.2f, nz * 2.2f);
            nx += 0.08f * wx; nz += 0.08f * wz;
            float lump = 1 - worley(nx * 6.f + 11.f, ny * 5.f, nz * 6.f);
            float bil = billow(nx * 13.f, ny * 11.f, nz * 13.f, oct);
            float amp = 0.75f * rs + 0.05f;
            float sdf = sd - amp * (0.6f * (lump - 0.42f) + 0.45f * (bil - 0.30f));
            float d = clampf(-sdf / 0.035f, 0, 1) * 0.95f;
            if (d > 0) {
                dsum += d; dustsum += d * 0.7f;
                if (detail == 1) {
                    float T = 0.16f + P.stemFire * (0.35f + 0.9f * expf(-(ytop - y) / 0.35f)) * ss(0.06f, -0.06f, sdf + 0.05f);
                    T += 0.25f * P.stemFire * expf(-y / 0.2f);
                    Tsum += d * T;
                }
            }
        }
    }

    if (P.base_on && y < P.hb * 4.2f + 0.15f) {
        float rb = P.rb * (1 + 0.13f * pn(cth * 1.5f + 3.1f, sth * 1.5f, 0.5f) + 0.06f * pn(cth * 4.f, sth * 4.f, 9.2f));
        float hb = P.hb * (1 + 0.7f * pn(cth * 2.4f + 1.7f, sth * 2.4f, 4.4f + 0.1f * t) + 0.5f * pn(x * 3, z * 3, 8.8f) + 1.6f * fmaxf(pn(x * 1.9f + 2.f, 0.3f * t, z * 1.9f), 0.f));
        hb = fmaxf(hb, 0.02f);
        float wr = P.wr;
        float ex = (rho - (rb - wr)) / wr, ey = y / hb;
        float sdw = (sqrtf(ex * ex + ey * ey) - 1) * fminf(wr, hb);
        float sdl = fmaxf(y - hb * 0.85f, rho - rb + wr * 0.5f);
        float sd = smin(sdw, sdl, 0.08f);
        sd = fmaxf(sd, -y);
        if (sd < 0.3f) {
            float bs = 0.45f + 0.55f * P.rb;
            float nx = x / bs, nz = z / bs, ny = y - 0.12f * t;
            float wx = pn(nx * 2.f, ny * 2.f, nz * 2.f + 1.f), wz = pn(nx * 2.f + 5.f, ny * 2.f, nz * 2.f);
            nx += 0.1f * wx; nz += 0.1f * wz;
            float lump = 1 - worley(nx * 4.6f, ny * 5.5f, nz * 4.6f);
            float bil = billow(nx * 9.f, ny * 10.f, nz * 9.f, oct + (detail == 1));
            float amp = 0.85f * hb + 0.03f;
            float sdf = sd - amp * (0.55f * (lump - 0.42f) + 0.6f * (bil - 0.30f));
            float d = clampf(-sdf / 0.07f, 0, 1) * 0.8f;
            if (d > 0) {
                dsum += d; dustsum += d;
                if (detail == 1) {
                    float T = 0.07f + P.baseFire * expf(-rho / 0.45f) * (0.6f + 0.4f * bil);
                    Tsum += d * T;
                }
            }
        }
    }
    r.d = fminf(dsum, 1.f);
    if (dsum > 0) { r.T = Tsum / dsum; r.dust = dustsum / dsum; r.soot = sootsum / dsum; }
    return r;
}

static const float RT[8] = {0, 0.25f, 0.5f, 0.75f, 1.0f, 1.3f, 1.7f, 2.2f};
static const float RC[8][3] = {{0, 0, 0}, {0.10f, 0.012f, 0.004f}, {0.40f, 0.055f, 0.012f}, {0.95f, 0.24f, 0.04f},
                               {1.55f, 0.62f, 0.12f}, {2.1f, 1.2f, 0.42f}, {2.7f, 2.1f, 1.15f}, {3.3f, 3.1f, 2.7f}};
static void ramp(float T, float *c) {
    if (T <= 0) { c[0] = c[1] = c[2] = 0; return; }
    if (T >= RT[7]) { float k = T / RT[7]; for (int i = 0; i < 3; i++) c[i] = RC[7][i] * k; return; }
    int j = 0; while (T > RT[j + 1]) j++;
    float f = (T - RT[j]) / (RT[j + 1] - RT[j]);
    for (int i = 0; i < 3; i++) c[i] = RC[j][i] + (RC[j + 1][i] - RC[j][i]) * f;
}

#define GX 120
#define GYN 80
#define GZ 120
static const float BX0 = -2.2f, BX1 = 2.2f, BY0 = -0.02f, BY1 = 2.75f;
static float *Dg, *Tg;
static float LX, LY, LZ;
static inline float gget(float *g, float x, float y, float z) {
    float fx = (x - BX0) / (BX1 - BX0) * (GX - 1), fy = (y - BY0) / (BY1 - BY0) * (GYN - 1), fz = (z - BX0) / (BX1 - BX0) * (GZ - 1);
    if (fx < 0 || fy < 0 || fz < 0 || fx >= GX - 1 || fy >= GYN - 1 || fz >= GZ - 1) return g == Tg ? 1.f : 0.f;
    int i = (int)fx, j = (int)fy, k = (int)fz; float u = fx - i, v = fy - j, w = fz - k;
#define G(a, b, c) g[((a) * GYN + (b)) * GZ + (c)]
    float c00 = mixf(G(i, j, k), G(i + 1, j, k), u), c10 = mixf(G(i, j + 1, k), G(i + 1, j + 1, k), u);
    float c01 = mixf(G(i, j, k + 1), G(i + 1, j, k + 1), u), c11 = mixf(G(i, j + 1, k + 1), G(i + 1, j + 1, k + 1), u);
    return mixf(mixf(c00, c10, v), mixf(c01, c11, v), w);
}
static void build_shadow(void) {
    Dg = calloc(GX * GYN * GZ, 4); Tg = calloc(GX * GYN * GZ, 4);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < GX; i++) for (int j = 0; j < GYN; j++) for (int k = 0; k < GZ; k++) {
        float x = BX0 + (BX1 - BX0) * i / (GX - 1), y = BY0 + (BY1 - BY0) * j / (GYN - 1), z = BX0 + (BX1 - BX0) * k / (GZ - 1);
        Dg[(i * GYN + j) * GZ + k] = bound(x, y, z) > 0.03f ? 0 : field(x, y, z, 0).d;
    }
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < GX; i++) for (int j = 0; j < GYN; j++) for (int k = 0; k < GZ; k++) {
        float x = BX0 + (BX1 - BX0) * i / (GX - 1), y = BY0 + (BY1 - BY0) * j / (GYN - 1), z = BX0 + (BX1 - BX0) * k / (GZ - 1);
        float tau = 0, st = 0.035f;
        for (int s = 1; s < 120; s++) {
            float px = x + LX * st * s, py = y + LY * st * s, pz = z + LZ * st * s;
            if (py > BY1 || px < BX0 || px > BX1 || pz < BX0 || pz > BX1) break;
            tau += gget(Dg, px, py, pz) * st;
        }
        Tg[(i * GYN + j) * GZ + k] = expf(-tau * 9.f);
    }
}

int main(int argc, char **argv) {
    float t = atof(argv[1]);
    if (argc > 3) { W = atoi(argv[3]); H = atoi(argv[4]); CX = atof(argv[5]); GYP = atof(argv[6]); SC = atof(argv[7]); ROLL = atof(argv[8]); }
    setpar(t);
    float ln = sqrtf(0.5f * 0.5f + 0.8f * 0.8f + 0.7f * 0.7f); LX = 0.5f / ln; LY = 0.8f / ln; LZ = 0.7f / ln;
    build_shadow();
    float *out = calloc(W * H * 4, 4);
    float sa = sinf(TILT), ca = cosf(TILT);
    float fx = 0, fy = -sa, fz = -ca;
    const float SIG = 85.f, DT = 0.0085f;
    float lp[3] = {0, P.yc, 0};
#pragma omp parallel for schedule(dynamic, 2)
    for (int py = 0; py < H; py++) for (int px = 0; px < W; px++) {
        float sx0 = (px + 0.5f - CX) / SC, sy0 = (GYP - (py + 0.5f)) / SC;
        float sx = sx0 * cosf(ROLL) - sy0 * sinf(ROLL), sy = sx0 * sinf(ROLL) + sy0 * cosf(ROLL);
        float ox = sx, oy = sy * ca + 4 * sa, oz = -sy * sa + 4 * ca;

        float t0 = 0, t1 = 1e9f;
        { float a = (-2.2f - oz) / fz, b = (2.2f - oz) / fz; if (a > b) { float c = a; a = b; b = c; } t0 = fmaxf(t0, a); t1 = fminf(t1, b); }
        if (ox < -2.2f || ox > 2.2f) continue;
        float tg = oy / sa;
        float ytopT = (BY1 - oy) / fy; if (ytopT > t0) t0 = ytopT;
        t1 = fminf(t1, tg);
        if (t0 >= t1) continue;
        unsigned hj = hsh(px, py, (int)(t * 1000));
        float tt = t0 + DT * ((hj & 4095) / 4096.f);
        float Tr = 1, C[3] = {0, 0, 0};
        while (tt < t1 && Tr > 0.004f) {
            float x = ox + fx * tt, y = oy + fy * tt, z = oz + fz * tt;
            float b = bound(x, y, z);
            if (b > 0.02f) { tt += fmaxf(DT, b * 0.8f); continue; }
            Smp s = field(x, y, z, 1);
            if (s.d > 0.001f) {
                float sig = s.d * SIG;
                float a = 1 - expf(-sig * DT);

                float d2 = field(x + LX * 0.035f, y + LY * 0.035f, z + LZ * 0.035f, 2).d;
                float Ls = gget(Tg, x, y, z) * expf(-1.8f * d2);

                float vx = lp[0] - x, vy = lp[1] - y, vz = lp[2] - z;
                float dl = sqrtf(vx * vx + vy * vy + vz * vz) + 1e-4f;
                float Lp = 0;
                if (s.dust > 0.05f) {
                    float d3 = field(x + vx / dl * 0.045f, y + vy / dl * 0.045f, z + vz / dl * 0.045f, 2).d;
                    Lp = P.Ip / (1 + 1.3f * dl * dl) * expf(-1.6f * d3) * s.dust;
                }
                float alb[3] = {mixf(0.75f, 0.62f, s.dust), mixf(0.5f, 0.40f, s.dust), mixf(0.36f, 0.30f, s.dust)};
                float sootk = 1 - 0.7f * s.soot;
                float e[3]; ramp(s.T, e);
                float key[3] = {1.0f, 0.62f, 0.36f}, pl[3] = {1.0f, 0.42f, 0.16f}, am[3] = {0.09f, 0.025f, 0.02f};
                for (int c = 0; c < 3; c++) {
                    float col = alb[c] * sootk * (key[c] * Ls * P.Ik * 0.5f + pl[c] * Lp + am[c] * P.amb) + e[c] * P.fade * 1.3f;
                    C[c] += Tr * a * col;
                }
                Tr *= 1 - a;
            }
            tt += DT;
        }

        if (Tr > 0.004f && tg < 1e8f && tt >= tg - DT) {
            float x = ox + fx * tg, z = oz + fz * tg;
            float dl2 = x * x + z * z + P.yc * P.yc;
            float g = 0.10f * P.Ip / (1 + 1.6f * dl2);
            float gc[3] = {1.0f, 0.42f, 0.16f};
            for (int c = 0; c < 3; c++) C[c] += Tr * g * gc[c];
        }
        float *o = out + (py * W + px) * 4;
        o[0] = C[0]; o[1] = C[1]; o[2] = C[2]; o[3] = 1 - Tr;
    }
    FILE *f = fopen(argv[2], "wb"); fwrite(out, 4, W * H * 4, f); fclose(f);

    float scy = GYP - SC * (P.yc * ca);
    printf("%f %f %f\n", CX, scy, P.ax * SC);
    return 0;
}
