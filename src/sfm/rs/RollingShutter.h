// Rolling-shutter estimation from neighbour matches alone, for photo sequences
// seen roughly straight down. docs/notes/rolling-shutter.md has the model, the
// measurements behind every constant here, and what it does not cover.
#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "sfm/geometry/Fundamental.h"
#include "sfm/geometry/Homography.h"
#include "sfm/geometry/LinAlg.h"
#include "sfm/optim/Ransac.h"

namespace sfm::rs {

// Per image: shift (2), height (1), attitude (3), then the camera's rotation
// during read-out as a cubic in the image row, three axes per order.
constexpr int kPer = 15;
// Shared by a window: camera mounting rotation (3), radial distortion k1.
constexpr int kGlobal = 4;

struct Options {
    double focal = 0;          // px; 0 scans kFocalGrid
    int max_gap = 5;           // sequence neighbours matched against each other
    int per_pair = 200;
    int window = 60, overlap = 12;
    int min_run = 6;           // shorter unbroken runs are left uncorrected
    int min_link = 80;         // loose inliers below which the sequence is cut
    double filter_px = 10;     // loose homography gate; a tight one drops the distortion itself
    double sig_px = 15, huber_px = 4;
    double h0 = 500, sig_xy = 300, sig_hstep = 1.5, sig_habs = 20;
    double sig_att_deg = 3, sig_yaw_deg = 8;
    double sig_rs = 0.006, sig_rs_high = 0.03;  // rad over the frame: constant rate, then the rest
    int threads = 0;
};

constexpr double kFocalGrid[] = {2500, 3500, 4500, 5500, 6500, 8000, 10323, 13000, 17000};

struct ImageModel {
    bool solved = false;
    int width = 0, height = 0;
    double f = 0, k1 = 0;
    double w[9] = {};          // [order][axis], camera axes, radians over the full read-out
    int window = -1;
};

struct PairMatches {
    int i = 0, j = 0;          // sequence positions, i < j
    std::vector<Vec2> p, q;    // loose inliers, full-size pixels
};

struct WindowReport {
    int first = 0, last = 0;   // sequence positions, inclusive
    double neighbours_px = 0, cost = 0, k1 = 0;
    double mount[3] = {};
};

struct Result {
    double focal = 0;
    std::vector<ImageModel> images;      // by sequence position
    std::vector<WindowReport> windows;
    std::vector<int> breaks;             // sequence positions a run ends at
    std::vector<std::pair<double, double>> scan;  // (focal, relative fit), empty when given
    int uncorrected = 0;
};

// ---------------------------------------------------------------------------
// The correction itself
// ---------------------------------------------------------------------------

inline Vec3 readoutRotation(const double* w, double s) {
    const double a = s, b = s * s - 1.0 / 12, c = s * s * s - 0.15 * s;
    return {w[0] * a + w[3] * b + w[6] * c, w[1] * a + w[4] * b + w[7] * c,
            w[2] * a + w[5] * b + w[8] * c};
}

// Second order in the angle; a frame's read-out rotation is a few degrees at most.
inline Vec3 rotateSmall(const Vec3& r, const Vec3& d) {
    const Vec3 c1 = r.cross(d);
    return d + c1 + r.cross(c1) * 0.5;
}

// A pixel of the photo as taken -> where an instant exposure at the mid-row
// pose would have put it. Lens distortion is kept.
inline Vec2 correctPoint(const ImageModel& m, const Vec2& p) {
    if (!m.solved) return p;
    const double cx = m.width * 0.5, cy = m.height * 0.5;
    const double xd = (p.x - cx) / m.f, yd = (p.y - cy) / m.f;
    const double g = 1 + m.k1 * (xd * xd + yd * yd);
    const Vec3 d = rotateSmall(readoutRotation(m.w, (p.y - cy) / m.height), {xd * g, yd * g, 1});
    const double xu = d.x / d.z, yu = d.y / d.z;
    double x = xu, y = yu;
    for (int it = 0; it < 12; it++) {
        const double gg = 1 + m.k1 * (x * x + y * y);
        x = xu / gg;
        y = yu / gg;
    }
    return {m.f * x + cx, m.f * y + cy};
}

// The inverse, for resampling: a pixel of the corrected photo -> where to read
// it in the photo as taken. The row decides the rotation, hence the iteration.
inline Vec2 sourcePoint(const ImageModel& m, const Vec2& p) {
    if (!m.solved) return p;
    const double cx = m.width * 0.5, cy = m.height * 0.5;
    const double xd = (p.x - cx) / m.f, yd = (p.y - cy) / m.f;
    const double g = 1 + m.k1 * (xd * xd + yd * yd);
    const Vec3 ray{xd * g, yd * g, 1};
    double u = p.x, v = p.y;
    for (int it = 0; it < 5; it++) {
        const Vec3 d = rotateSmall(readoutRotation(m.w, (v - cy) / m.height) * -1.0, ray);
        const double xu = d.x / d.z, yu = d.y / d.z;
        double x = xu, y = yu;
        for (int k = 0; k < 12; k++) {
            const double gg = 1 + m.k1 * (x * x + y * y);
            x = xu / gg;
            y = yu / gg;
        }
        u = m.f * x + cx;
        v = m.f * y + cy;
    }
    return {u, v};
}

// ---------------------------------------------------------------------------
// Estimation
// ---------------------------------------------------------------------------

namespace detail {

inline Mat3 expm(const Vec3& r) {
    const double th = r.norm();
    const Mat3 K = crossMatrix(r);
    const Mat3 KK = mul(K, K);
    const double a = th < 1e-8 ? 1.0 - th * th / 6 : std::sin(th) / th;
    const double b = th < 1e-8 ? 0.5 - th * th / 24 : (1 - std::cos(th)) / (th * th);
    Mat3 R = mat3Identity();
    for (int k = 0; k < 9; k++) R[k] += a * K[k] + b * KK[k];
    return R;
}

struct Obs {
    int i, j;                  // window-local
    Vec2 p, q;
    bool neighbours;
};

struct Problem {
    int n = 0, width = 0, height = 0;
    double f = 0;
    std::vector<Obs> obs;
    std::vector<Vec2> pos0;    // chained start, metres at h0 / f per pixel
    std::vector<Mat3> base;    // camera-to-world of a level camera at the chained heading
};

// World: x east, y north, z up, the first photo's image-up called north.
inline Vec2 toGround(const Problem& pb, int i, const double* xi, const double* g, const Vec2& uv,
                     const Mat3& Rwc) {
    const double cx = pb.width * 0.5, cy = pb.height * 0.5;
    double x = (uv.x - cx) / pb.f, y = (uv.y - cy) / pb.f;
    const double gg = 1 + g[3] * (x * x + y * y);
    x *= gg;
    y *= gg;
    const Vec3 d = rotateSmall(readoutRotation(xi + 6, (uv.y - cy) / pb.height), {x, y, 1});
    const Vec3 dw = mul(Rwc, d);
    const double k = xi[2] / -dw.z;
    return {pb.pos0[i].x + xi[0] + dw.x * k, pb.pos0[i].y + xi[1] + dw.y * k};
}

inline Mat3 cameraToWorld(const Problem& pb, int i, const double* xi, const double* g) {
    return mul(mul(pb.base[i], expm({g[0], g[1], g[2]})), expm({xi[3], xi[4], xi[5]}));
}

// Ground mismatch of one match, in pixels over sig_px.
inline void residual(const Problem& pb, const Obs& o, const double* xi, const double* xj,
                     const double* g, const Vec2& gi, const Vec2& gj, double inv_sig, double r[2]) {
    (void)o;
    (void)g;
    const double s = pb.f / (0.5 * (xi[2] + xj[2])) * inv_sig;
    r[0] = (gi.x - gj.x) * s;
    r[1] = (gi.y - gj.y) * s;
}

// In-place Cholesky solve of a dense symmetric system; false when not positive definite.
inline bool choleskySolve(std::vector<double>& A, std::vector<double>& b, int n) {
    for (int j = 0; j < n; j++) {
        double d = A[(size_t)j * n + j];
        for (int k = 0; k < j; k++) d -= A[(size_t)j * n + k] * A[(size_t)j * n + k];
        if (!(d > 0)) return false;
        d = std::sqrt(d);
        A[(size_t)j * n + j] = d;
        for (int i = j + 1; i < n; i++) {
            double s = A[(size_t)i * n + j];
            const double* ri = &A[(size_t)i * n];
            const double* rj = &A[(size_t)j * n];
            for (int k = 0; k < j; k++) s -= ri[k] * rj[k];
            A[(size_t)i * n + j] = s / d;
        }
    }
    for (int i = 0; i < n; i++) {
        double s = b[i];
        for (int k = 0; k < i; k++) s -= A[(size_t)i * n + k] * b[k];
        b[i] = s / A[(size_t)i * n + i];
    }
    for (int i = n - 1; i >= 0; i--) {
        double s = b[i];
        for (int k = i + 1; k < n; k++) s -= A[(size_t)k * n + i] * b[k];
        b[i] = s / A[(size_t)i * n + i];
    }
    return true;
}

struct Solved {
    std::vector<double> x;
    double cost = 1e300, neighbours_px = 0;
};

class Solver {
public:
    Solver(const Problem& pb, const Options& opt) : pb_(pb), opt_(opt) {
        N_ = pb.n * kPer + kGlobal;
        inv_sig_ = 1.0 / opt.sig_px;
        huber_ = opt.huber_px / opt.sig_px;
        const double d2r = 3.14159265358979323846 / 180;
        sig_att_ = opt.sig_att_deg * d2r;
        sig_yaw_ = opt.sig_yaw_deg * d2r;
    }

    Solved run() {
        std::vector<double> x(N_, 0.0);
        for (int i = 0; i < pb_.n; i++) x[(size_t)i * kPer + 2] = opt_.h0;
        double c = cost(x), lam = 1e-3;
        std::vector<double> A, g, At, dx, xt(N_);
        for (int it = 0; it < 60; it++) {
            normalEquations(x, A, g);
            bool ok = false;
            double c2 = c;
            for (int tr = 0; tr < 12 && !ok; tr++) {
                At = A;
                dx.assign(N_, 0.0);
                for (int k = 0; k < N_; k++) {
                    At[(size_t)k * N_ + k] += lam * (A[(size_t)k * N_ + k] + 1e-12);
                    dx[k] = -g[k];
                }
                if (choleskySolve(At, dx, N_)) {
                    for (int k = 0; k < N_; k++) xt[k] = x[k] + dx[k];
                    c2 = cost(xt);
                    ok = c2 < c;
                }
                if (!ok) lam *= 5;
            }
            if (!ok) break;
            const double gain = (c - c2) / c;
            x = xt;
            c = c2;
            lam = std::max(lam / 5, 1e-9);
            if (gain < 1e-6) break;
        }
        Solved s;
        s.x = x;
        s.cost = c;
        s.neighbours_px = neighboursMedian(x);
        return s;
    }

private:
    const Problem& pb_;
    const Options& opt_;
    int N_ = 0;
    double inv_sig_ = 0, huber_ = 0, sig_att_ = 0, sig_yaw_ = 0;

    double priorSigma(int k) const {   // k within an image's kPer
        if (k < 2) return opt_.sig_xy;
        if (k == 2) return opt_.sig_habs;
        // Attitude priors are on the correction in body axes: body roll is
        // -camera y, pitch is camera x, yaw is camera z.
        if (k == 3 || k == 4) return sig_att_;
        if (k == 5) return sig_yaw_;
        return k < 9 ? opt_.sig_rs : opt_.sig_rs_high;
    }

    void grounds(const std::vector<double>& x, std::vector<Mat3>& R) const {
        R.resize(pb_.n);
        const double* g = &x[(size_t)pb_.n * kPer];
        for (int i = 0; i < pb_.n; i++) R[i] = cameraToWorld(pb_, i, &x[(size_t)i * kPer], g);
    }

    double priorCost(const std::vector<double>& x) const {
        double c = 0;
        for (int i = 0; i < pb_.n; i++) {
            const double* xi = &x[(size_t)i * kPer];
            for (int k = 0; k < kPer; k++) {
                const double v = (k == 2 ? xi[2] - opt_.h0 : xi[k]) / priorSigma(k);
                c += 0.5 * v * v;
            }
            if (i + 1 < pb_.n) {
                const double v = (x[(size_t)(i + 1) * kPer + 2] - xi[2]) / opt_.sig_hstep;
                c += 0.5 * v * v;
            }
        }
        return c;
    }

    double cost(const std::vector<double>& x) const {
        std::vector<Mat3> R;
        grounds(x, R);
        const double* g = &x[(size_t)pb_.n * kPer];
        double c = priorCost(x);
        for (const Obs& o : pb_.obs) {
            const double* xi = &x[(size_t)o.i * kPer];
            const double* xj = &x[(size_t)o.j * kPer];
            double r[2];
            residual(pb_, o, xi, xj, g, toGround(pb_, o.i, xi, g, o.p, R[o.i]),
                     toGround(pb_, o.j, xj, g, o.q, R[o.j]), inv_sig_, r);
            const double m = std::sqrt(r[0] * r[0] + r[1] * r[1]);
            c += m > huber_ ? huber_ * m - 0.5 * huber_ * huber_ : 0.5 * m * m;
        }
        return c;
    }

    double neighboursMedian(const std::vector<double>& x) const {
        std::vector<Mat3> R;
        grounds(x, R);
        const double* g = &x[(size_t)pb_.n * kPer];
        std::vector<double> e;
        for (const Obs& o : pb_.obs) {
            if (!o.neighbours) continue;
            const double* xi = &x[(size_t)o.i * kPer];
            const double* xj = &x[(size_t)o.j * kPer];
            double r[2];
            residual(pb_, o, xi, xj, g, toGround(pb_, o.i, xi, g, o.p, R[o.i]),
                     toGround(pb_, o.j, xj, g, o.q, R[o.j]), 1.0, r);
            e.push_back(std::sqrt(r[0] * r[0] + r[1] * r[1]));
        }
        if (e.empty()) return 0;
        std::nth_element(e.begin(), e.begin() + e.size() / 2, e.end());
        return e[e.size() / 2];
    }

    // Forward differences over the 34 unknowns one match touches; metres get a
    // millimetre step, angles 1e-7 rad, k1 1e-6.
    void normalEquations(const std::vector<double>& x, std::vector<double>& A,
                         std::vector<double>& gr) const {
        A.assign((size_t)N_ * N_, 0.0);
        gr.assign(N_, 0.0);
        std::vector<Mat3> R;
        grounds(x, R);
        const int go = pb_.n * kPer;
        const double* g = &x[go];
        double step[kPer + kGlobal];
        for (int k = 0; k < kPer; k++) step[k] = k < 3 ? 1e-3 : 1e-7;
        step[kPer] = step[kPer + 1] = step[kPer + 2] = 1e-7;
        step[kPer + 3] = 1e-6;

        const int M = 2 * kPer + kGlobal;
        double J[2][2 * kPer + kGlobal];
        int cols[2 * kPer + kGlobal];
        double xi[kPer], xj[kPer], gg[kGlobal];
        for (const Obs& o : pb_.obs) {
            for (int k = 0; k < kPer; k++) {
                xi[k] = x[(size_t)o.i * kPer + k];
                xj[k] = x[(size_t)o.j * kPer + k];
                cols[k] = o.i * kPer + k;
                cols[kPer + k] = o.j * kPer + k;
            }
            for (int k = 0; k < kGlobal; k++) {
                gg[k] = g[k];
                cols[2 * kPer + k] = go + k;
            }
            const Vec2 gi = toGround(pb_, o.i, xi, gg, o.p, R[o.i]);
            const Vec2 gj = toGround(pb_, o.j, xj, gg, o.q, R[o.j]);
            double r0[2], r1[2];
            residual(pb_, o, xi, xj, gg, gi, gj, inv_sig_, r0);
            for (int k = 0; k < kPer; k++) {
                const double keep = xi[k];
                xi[k] += step[k];
                const Mat3 Rk = k >= 3 && k < 6 ? cameraToWorld(pb_, o.i, xi, gg) : R[o.i];
                residual(pb_, o, xi, xj, gg, toGround(pb_, o.i, xi, gg, o.p, Rk), gj, inv_sig_, r1);
                xi[k] = keep;
                J[0][k] = (r1[0] - r0[0]) / step[k];
                J[1][k] = (r1[1] - r0[1]) / step[k];
            }
            for (int k = 0; k < kPer; k++) {
                const double keep = xj[k];
                xj[k] += step[k];
                const Mat3 Rk = k >= 3 && k < 6 ? cameraToWorld(pb_, o.j, xj, gg) : R[o.j];
                residual(pb_, o, xi, xj, gg, gi, toGround(pb_, o.j, xj, gg, o.q, Rk), inv_sig_, r1);
                xj[k] = keep;
                J[0][kPer + k] = (r1[0] - r0[0]) / step[k];
                J[1][kPer + k] = (r1[1] - r0[1]) / step[k];
            }
            for (int k = 0; k < kGlobal; k++) {
                const double keep = gg[k];
                gg[k] += step[kPer + k];
                const Mat3 Ri = k < 3 ? cameraToWorld(pb_, o.i, xi, gg) : R[o.i];
                const Mat3 Rj = k < 3 ? cameraToWorld(pb_, o.j, xj, gg) : R[o.j];
                residual(pb_, o, xi, xj, gg, toGround(pb_, o.i, xi, gg, o.p, Ri),
                         toGround(pb_, o.j, xj, gg, o.q, Rj), inv_sig_, r1);
                gg[k] = keep;
                J[0][2 * kPer + k] = (r1[0] - r0[0]) / step[kPer + k];
                J[1][2 * kPer + k] = (r1[1] - r0[1]) / step[kPer + k];
            }
            const double m = std::sqrt(r0[0] * r0[0] + r0[1] * r0[1]);
            const double wgt = m > huber_ ? huber_ / m : 1.0;
            for (int a = 0; a < M; a++) {
                const double ja0 = J[0][a] * wgt, ja1 = J[1][a] * wgt;
                gr[cols[a]] += ja0 * r0[0] + ja1 * r0[1];
                double* row = &A[(size_t)cols[a] * N_];
                for (int b = 0; b < M; b++) row[cols[b]] += ja0 * J[0][b] + ja1 * J[1][b];
            }
        }
        for (int i = 0; i < pb_.n; i++) {
            for (int k = 0; k < kPer; k++) {
                const int c = i * kPer + k;
                const double s = 1.0 / priorSigma(k);
                A[(size_t)c * N_ + c] += s * s;
                gr[c] += s * s * (k == 2 ? x[c] - opt_.h0 : x[c]);
            }
            if (i + 1 < pb_.n) {
                const int a = i * kPer + 2, b = (i + 1) * kPer + 2;
                const double s2 = 1.0 / (opt_.sig_hstep * opt_.sig_hstep), d = x[b] - x[a];
                A[(size_t)a * N_ + a] += s2;
                A[(size_t)b * N_ + b] += s2;
                A[(size_t)a * N_ + b] -= s2;
                A[(size_t)b * N_ + a] -= s2;
                gr[a] -= s2 * d;
                gr[b] += s2 * d;
            }
        }
    }
};

// q = R p + t between two photos, trimmed three times at 2.5 medians.
inline bool rigidFit(const std::vector<Vec2>& p, const std::vector<Vec2>& q, double R[4],
                     double t[2]) {
    const size_t n = p.size();
    if (n < 8) return false;
    std::vector<char> keep(n, 1);
    std::vector<double> e(n), es;
    for (int round = 0; round < 3; round++) {
        double pc[2] = {0, 0}, qc[2] = {0, 0}, cnt = 0;
        for (size_t k = 0; k < n; k++)
            if (keep[k]) {
                pc[0] += p[k].x; pc[1] += p[k].y; qc[0] += q[k].x; qc[1] += q[k].y; cnt++;
            }
        if (cnt < 8) return false;
        for (int a = 0; a < 2; a++) { pc[a] /= cnt; qc[a] /= cnt; }
        double sxx = 0, sxy = 0;   // sum of dot and cross of centred pairs
        for (size_t k = 0; k < n; k++)
            if (keep[k]) {
                const double ax = p[k].x - pc[0], ay = p[k].y - pc[1];
                const double bx = q[k].x - qc[0], by = q[k].y - qc[1];
                sxx += ax * bx + ay * by;
                sxy += ax * by - ay * bx;
            }
        const double th = std::atan2(sxy, sxx), c = std::cos(th), s = std::sin(th);
        R[0] = c; R[1] = -s; R[2] = s; R[3] = c;
        t[0] = qc[0] - (c * pc[0] - s * pc[1]);
        t[1] = qc[1] - (s * pc[0] + c * pc[1]);
        es.clear();
        for (size_t k = 0; k < n; k++) {
            const double dx = R[0] * p[k].x + R[1] * p[k].y + t[0] - q[k].x;
            const double dy = R[2] * p[k].x + R[3] * p[k].y + t[1] - q[k].y;
            e[k] = std::sqrt(dx * dx + dy * dy);
            es.push_back(e[k]);
        }
        std::nth_element(es.begin(), es.begin() + es.size() / 2, es.end());
        const double med = es[es.size() / 2];
        for (size_t k = 0; k < n; k++) keep[k] = e[k] < 2.5 * med;
    }
    return true;
}

// Up to per_pair matches, spread over a 12 x 9 grid of the first photo.
inline std::vector<int> spread(const std::vector<Vec2>& p, int width, int height, int per_pair,
                               unsigned seed) {
    std::vector<int> idx(p.size());
    for (size_t k = 0; k < p.size(); k++) idx[k] = (int)k;
    if ((int)p.size() <= per_pair) return idx;
    std::mt19937 rng(seed);
    std::shuffle(idx.begin(), idx.end(), rng);
    std::vector<int> seen(12 * 9, 0), rank(p.size());
    for (size_t k = 0; k < idx.size(); k++) {
        const Vec2& v = p[idx[k]];
        const int cx = std::min(11, std::max(0, (int)(v.x / (width / 12.0))));
        const int cy = std::min(8, std::max(0, (int)(v.y / (height / 9.0))));
        rank[k] = seen[cx * 9 + cy]++;
    }
    std::vector<int> order(idx.size());
    for (size_t k = 0; k < order.size(); k++) order[k] = (int)k;
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return rank[a] < rank[b]; });
    std::vector<int> out(per_pair);
    for (int k = 0; k < per_pair; k++) out[k] = idx[order[k]];
    return out;
}

struct Window {
    int start = 0, end = 0, own0 = 0, own1 = 0;   // sequence positions, end exclusive
};

inline Problem buildProblem(const Window& w, const std::vector<const PairMatches*>& pairs,
                            int width, int height, double f, const Options& opt) {
    Problem pb;
    pb.n = w.end - w.start;
    pb.width = width;
    pb.height = height;
    pb.f = f;
    std::vector<const PairMatches*> chain(pb.n, nullptr);
    for (const PairMatches* pm : pairs) {
        if (pm->i < w.start || pm->j >= w.end) continue;
        if (pm->j == pm->i + 1) chain[pm->i - w.start] = pm;
        for (int k : spread(pm->p, width, height, opt.per_pair, (unsigned)(pm->i * 131 + pm->j)))
            pb.obs.push_back({pm->i - w.start, pm->j - w.start, pm->p[k], pm->q[k],
                              pm->j == pm->i + 1});
    }
    // Chain: M maps a photo's pixels into the first photo's pixel frame.
    std::vector<std::array<double, 6>> M(pb.n);
    M[0] = {1, 0, 0, 0, 1, 0};
    for (int k = 0; k + 1 < pb.n; k++) {
        double R[4] = {1, 0, 0, 1}, t[2] = {0, 0};
        if (chain[k]) rigidFit(chain[k]->p, chain[k]->q, R, t);
        // inverse of q = R p + t, composed after M[k]
        const double Ri[4] = {R[0], R[2], R[1], R[3]};
        const double ti[2] = {-(Ri[0] * t[0] + Ri[1] * t[1]), -(Ri[2] * t[0] + Ri[3] * t[1])};
        const auto& m = M[k];
        M[k + 1] = {m[0] * Ri[0] + m[1] * Ri[2], m[0] * Ri[1] + m[1] * Ri[3],
                    m[0] * ti[0] + m[1] * ti[1] + m[2],
                    m[3] * Ri[0] + m[4] * Ri[2], m[3] * Ri[1] + m[4] * Ri[3],
                    m[3] * ti[0] + m[4] * ti[1] + m[5]};
    }
    const double gsd = opt.h0 / f, cx = width * 0.5, cy = height * 0.5;
    const double x0 = M[0][0] * cx + M[0][1] * cy + M[0][2], y0 = M[0][3] * cx + M[0][4] * cy + M[0][5];
    pb.pos0.resize(pb.n);
    pb.base.resize(pb.n);
    for (int k = 0; k < pb.n; k++) {
        const auto& m = M[k];
        pb.pos0[k] = {(m[0] * cx + m[1] * cy + m[2] - x0) * gsd,
                      -(m[3] * cx + m[4] * cy + m[5] - y0) * gsd};
        // Bearing of image-up; the level camera has x east-rotated, y the opposite of up, z down.
        const double yaw = std::atan2(-m[1], m[4]);
        const double c = std::cos(yaw), s = std::sin(yaw);
        pb.base[k] = {c, -s, 0, -s, -c, 0, 0, 0, -1};
    }
    return pb;
}

template <class F>
inline void parallelFor(int n, int threads, const F& fn) {
    const unsigned hc = std::thread::hardware_concurrency();
    int nt = threads > 0 ? threads : (hc > 0 ? (int)hc : 1);
    nt = std::max(1, std::min(nt, n));
    std::atomic<int> next{0};
    std::vector<std::thread> pool;
    for (int t = 0; t < nt; t++)
        pool.emplace_back([&] {
            for (int k = next++; k < n; k = next++) fn(k);
        });
    for (std::thread& t : pool) t.join();
}

}  // namespace detail

// Loose homography gate on one pair's putative matches; empty when the pair
// does not hold together.
inline PairMatches filterPair(int i, int j, const std::vector<Vec2>& p, const std::vector<Vec2>& q,
                              const Options& opt) {
    PairMatches out;
    out.i = i;
    out.j = j;
    if (p.size() < 40) return out;
    RansacOptions ro;
    ro.max_error = opt.filter_px;
    ro.confidence = 0.9999;
    ro.max_num_trials = 5000;
    ro.seed = (unsigned)(i * 7919 + j);
    auto fit = [&](const std::vector<int>& idx) { return estimateHomography(p, q, idx); };
    auto res = [&](const Mat3& H, int k) {
        const Vec2 t = applyH(H, p[k]);
        const double dx = t.x - q[k].x, dy = t.y - q[k].y;
        return dx * dx + dy * dy;
    };
    const RansacReport<Mat3> rep = loransac<Mat3>((int)p.size(), 4, fit, fit, res, ro);
    if (!rep.success || rep.num_inliers < 40) return out;
    for (size_t k = 0; k < p.size(); k++)
        if (rep.inlier_mask[k]) {
            out.p.push_back(p[k]);
            out.q.push_back(q[k]);
        }
    return out;
}

// `pairs` between sequence positions at most max_gap apart; images of one size.
inline Result estimate(int num_images, int width, int height,
                       const std::vector<PairMatches>& pairs, const Options& opt) {
    using namespace detail;
    Result res;
    res.images.assign(num_images, ImageModel());
    for (ImageModel& m : res.images) {
        m.width = width;
        m.height = height;
    }
    std::vector<char> linked(std::max(0, num_images - 1), 0);
    std::vector<const PairMatches*> ptr;
    for (const PairMatches& pm : pairs) {
        if (pm.p.empty() || pm.j - pm.i > opt.max_gap) continue;
        ptr.push_back(&pm);
        if (pm.j == pm.i + 1 && (int)pm.p.size() >= opt.min_link) linked[pm.i] = 1;
    }
    std::vector<Window> windows;
    for (int r0 = 0; r0 < num_images;) {
        int r1 = r0 + 1;
        while (r1 < num_images && linked[r1 - 1]) r1++;
        if (r1 < num_images) res.breaks.push_back(r1 - 1);
        if (r1 - r0 < opt.min_run) {
            res.uncorrected += r1 - r0;
        } else {
            const int step = opt.window - opt.overlap;
            for (int w0 = r0, k = 0;; w0 += step, k++) {
                int w1 = std::min(r1, w0 + opt.window);
                if (r1 - w1 < opt.min_run) w1 = r1;
                windows.push_back({w0, w1, k == 0 ? r0 : w0 + opt.overlap / 2,
                                   w1 == r1 ? r1 : w1 - opt.overlap / 2});
                if (w1 == r1) break;
            }
        }
        r0 = r1;
    }
    if (windows.empty()) return res;

    double f = opt.focal;
    if (!(f > 0)) {
        // The three windows with the most matches decide the lens.
        std::vector<std::pair<size_t, int>> load;
        for (size_t k = 0; k < windows.size(); k++) {
            size_t c = 0;
            for (const PairMatches* pm : ptr)
                if (pm->i >= windows[k].start && pm->j < windows[k].end) c += pm->p.size();
            load.push_back({c, (int)k});
        }
        std::sort(load.rbegin(), load.rend());
        const int nw = (int)std::min<size_t>(3, load.size());
        const int ng = (int)(sizeof(kFocalGrid) / sizeof(kFocalGrid[0]));
        std::vector<double> cost((size_t)nw * ng);
        parallelFor(nw * ng, opt.threads, [&](int t) {
            const Problem pb =
                buildProblem(windows[load[t / ng].second], ptr, width, height, kFocalGrid[t % ng], opt);
            cost[t] = Solver(pb, opt).run().cost;
        });
        int best = 0;
        double best_v = 1e300;
        for (int gI = 0; gI < ng; gI++) {
            double v = 0;
            for (int wI = 0; wI < nw; wI++) {
                const double mn = *std::min_element(cost.begin() + (size_t)wI * ng,
                                                    cost.begin() + (size_t)(wI + 1) * ng);
                v += cost[(size_t)wI * ng + gI] / mn / nw;
            }
            res.scan.push_back({kFocalGrid[gI], v});
            if (v < best_v) {
                best_v = v;
                best = gI;
            }
        }
        f = kFocalGrid[best];
    }
    res.focal = f;

    std::vector<Solved> solved(windows.size());
    parallelFor((int)windows.size(), opt.threads, [&](int k) {
        const Problem pb = buildProblem(windows[k], ptr, width, height, f, opt);
        solved[k] = Solver(pb, opt).run();
    });
    for (size_t k = 0; k < windows.size(); k++) {
        const Window& w = windows[k];
        const std::vector<double>& x = solved[k].x;
        const double* g = &x[(size_t)(w.end - w.start) * kPer];
        WindowReport wr;
        wr.first = w.start;
        wr.last = w.end - 1;
        wr.neighbours_px = solved[k].neighbours_px;
        wr.cost = solved[k].cost;
        wr.k1 = g[3];
        for (int a = 0; a < 3; a++) wr.mount[a] = g[a];
        res.windows.push_back(wr);
        for (int s = w.own0; s < w.own1; s++) {
            ImageModel& m = res.images[s];
            m.solved = true;
            m.f = f;
            m.k1 = g[3];
            m.window = (int)k;
            for (int a = 0; a < 9; a++) m.w[a] = x[(size_t)(s - w.start) * kPer + 6 + a];
        }
    }
    return res;
}

}  // namespace sfm::rs
