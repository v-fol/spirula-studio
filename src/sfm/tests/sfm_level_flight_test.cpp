// Level-flight priors (map/LevelFlightPriors.h): relative height factors
// between cameras that are near each other over the ground, for a nadir survey
// that holds its altitude (host only).
//
// Prints PASS/FAIL and returns 0/1. See docs/testing.md.
#include <cmath>
#include <cstdio>
#include <vector>

#include "sfm/ba/Priors.h"
#include "sfm/ba/Problem.h"
#include "sfm/core/Pose.h"
#include "sfm/map/LevelFlightPriors.h"
#include "sfm/tests/TestMain.h"

using namespace sfm;

static int fails = 0;
static void check(bool ok, const char* what) {
    printf("  %-64s %s\n", what, ok ? "ok" : "FAILED");
    if (!ok) fails++;
}

// A camera at C looking straight down -Z of a frame turned by Rw, its image x
// along the flight.
static Pose nadir(const Vec3& C, const Mat3& Rw) {
    // camera axes in the survey frame: x east, y south, z down
    const Mat3 Rc = {1, 0, 0, 0, -1, 0, 0, 0, -1};
    const Mat3 R = mul(Rc, transpose(Rw));
    const Vec3 t = mul(R, C);
    return {R, {-t.x, -t.y, -t.z}};
}

// Three flight lines along X, 60 m apart, a photo every 40 m, 100 m up; line
// `lifted` raised by `dz`. World = survey frame turned by Rw, scaled by s.
static std::vector<PosedImage> survey(const Mat3& Rw, double s, int lifted = -1, double dz = 0) {
    std::vector<PosedImage> imgs;
    uint32_t id = 10;
    for (int line = 0; line < 3; line++)
        for (int k = 0; k < 12; k++) {
            const Vec3 local{40.0 * k, 60.0 * line, 100.0 + (line == lifted ? dz : 0.0)};
            PosedImage p;
            p.image = id++;
            p.pose = nadir(mul(Rw, local) * s, Rw);
            imgs.push_back(p);
        }
    return imgs;
}

// sum A_k c_k - b of a factor, in sigmas along its weighted axis.
static double residual(const PriorCentre& f, const std::vector<PosedImage>& imgs) {
    Vec3 sum{0, 0, 0};
    for (int k = 0; k < f.n; k++)
        for (const PosedImage& p : imgs)
            if (p.image == f.img[k]) sum = sum + mul(f.A[k], cameraCenter(p.pose));
    return (sum.z - f.b.z) / f.sigma.z;
}

static int line_of(uint32_t image) { return (int)((image - 10) / 12); }

int cmdLevelFlightSelftest(int, char**) {
    LevelFlightOptions opt;
    opt.slope = 0.02;
    opt.neighbours = 6;
    const Mat3 turn = angleAxisToRotation(Vec3{0.3, -0.2, 0.9});

    // ---- T1: a level survey: factors, all satisfied --------------------------
    {
        const std::vector<PosedImage> imgs = survey(turn, 0.37);
        PosePriors pr;
        LevelFlightStats st;
        const int n = levelFlightFactors(imgs, opt, pr, &st);
        printf("  T1: %d factors over %d cameras\n", n, st.images);
        check(st.nadir && n > 0 && n == (int)pr.centres.size(), "T1: factors made for a nadir survey");
        bool shape = true, cross = false;
        double worst = 0;
        for (const PriorCentre& f : pr.centres) {
            shape = shape && f.n == 2 && f.sigma.x == 0 && f.sigma.y == 0 && f.sigma.z > 0 &&
                    !f.absolute();
            cross = cross || line_of(f.img[0]) != line_of(f.img[1]);
            worst = std::max(worst, std::fabs(residual(f, imgs)));
        }
        check(shape, "T1: two-camera height factors, level axes free");
        check(cross, "T1: some factors tie neighbouring flight lines");
        check(worst < 1e-9, "T1: a level survey satisfies every factor");
        check(!pr.hasAbsoluteCentres(), "T1: no absolute position stated");
    }

    // ---- T2: a lifted line shows only across the lines, at dz / sigma ---------
    {
        const double s = 0.37, dz = 1.5;
        const std::vector<PosedImage> imgs = survey(turn, s, 1, dz);
        PosePriors pr;
        levelFlightFactors(imgs, opt, pr);
        bool within_zero = true, across_ok = true;
        int across = 0;
        for (const PriorCentre& f : pr.centres) {
            const double r = residual(f, imgs);
            const bool crossing = (line_of(f.img[0]) == 1) != (line_of(f.img[1]) == 1);
            if (!crossing) within_zero = within_zero && std::fabs(r) < 1e-6;
            else {
                across++;
                // d over the ground is s * 60 (or more, for a diagonal neighbour)
                const double expect = dz * s / f.sigma.z;
                across_ok = across_ok && std::fabs(std::fabs(r) - expect) < 1e-6 * expect;
            }
        }
        printf("  T2: %d factors across the lifted line\n", across);
        check(across > 0 && across_ok, "T2: lifted line costs dz / sigma across");
        check(within_zero, "T2: and nothing along the lines");
    }

    // ---- T3: oblique cameras state nothing -----------------------------------
    {
        std::vector<PosedImage> imgs;
        for (int k = 0; k < 12; k++) {
            // an orbit looking at its centre from the side
            const double a = 2 * M_PI * k / 12;
            const Vec3 C{50 * std::cos(a), 50 * std::sin(a), 10};
            const Vec3 f = (Vec3{0, 0, 0} - C).normalized();
            const Vec3 r = f.cross(Vec3{0, 0, 1}).normalized();
            const Vec3 u = f.cross(r);
            const Mat3 R = {r.x, r.y, r.z, u.x, u.y, u.z, f.x, f.y, f.z};
            const Vec3 t = mul(R, C);
            PosedImage p;
            p.image = (uint32_t)k;
            p.pose = {R, {-t.x, -t.y, -t.z}};
            imgs.push_back(p);
        }
        PosePriors pr;
        LevelFlightStats st;
        check(levelFlightFactors(imgs, opt, pr, &st) == 0 && !st.nadir,
              "T3: an orbit gets no factors");
    }

    // ---- T4: the source wraps another one ------------------------------------
    {
        struct Inner : PriorSource {
            bool has(uint32_t) const override { return true; }
            bool relativeRotation(uint32_t, uint32_t, Mat3&, double&) const override { return false; }
            std::vector<uint32_t> neighbours(uint32_t) const override { return {7}; }
            PosePriors factors(const std::vector<PosedImage>&) override {
                PosePriors p;
                PriorCentre c;
                c.n = 1;
                c.img[0] = 10;
                p.centres.push_back(c);
                return p;
            }
        } inner;
        const std::vector<PosedImage> imgs = survey(turn, 1.0);
        LevelFlightPriors alone(nullptr, opt), wrapped(&inner, opt);
        const size_t own = alone.factors(imgs).centres.size();
        const PosePriors both = wrapped.factors(imgs);
        check(!alone.has(10) && wrapped.has(10) && wrapped.neighbours(3) == std::vector<uint32_t>{7},
              "T4: questions go to the wrapped source");
        check(both.centres.size() == own + 1 && both.centres[0].absolute(),
              "T4: its factors first, the level ones added");
        check(wrapped.lastFactors().factors == (int)own, "T4: stats count the level ones");
    }

    // ---- T5: in the BA's prior assembly: zero level, a cost lifted ------------
    {
        auto cost = [&](const std::vector<PosedImage>& imgs) {
            PosePriors pr;
            LevelFlightStats st;
            // Image ids are BA indices here.
            std::vector<PosedImage> ba = imgs;
            for (size_t i = 0; i < ba.size(); i++) ba[i].image = (uint32_t)i;
            levelFlightFactors(ba, opt, pr, &st);
            BAProblem P;
            P.num_images = (uint32_t)ba.size();
            P.identityFrames();
            P.priors = &pr;
            std::vector<double> poses;
            for (const PosedImage& p : ba) {
                const Vec3 aa = rotationToAngleAxis(p.pose.R);
                poses.insert(poses.end(), {aa.x, aa.y, aa.z, p.pose.t.x, p.pose.t.y, p.pose.t.z});
            }
            PriorAssembler A;
            A.init(P);
            return A.cost(P, poses.data(), nullptr);
        };
        const double level = cost(survey(turn, 0.37)), lifted = cost(survey(turn, 0.37, 2, 3.0));
        printf("  T5: prior cost level %.3e, one line 3 m up %.3e\n", level, lifted);
        check(level < 1e-12 && lifted > 1.0, "T5: the solver's prior cost sees the lifted line");
    }

    printf("%s\n", fails ? "FAIL" : "PASS");
    return fails ? 1 : 0;
}

int main(int argc, char** argv) { return sfmTestMain(argc, argv, cmdLevelFlightSelftest); }
