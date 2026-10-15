// Level flight as a prior source (`sfm auto --level-flight SLOPE`): an aircraft
// that holds its altitude puts cameras that are near each other over the ground
// at nearly the same height. Straight-down photos over flat ground tie
// neighbouring flight lines together in height only weakly, so without a
// position source a line and the ground under it can rise or sink together at
// almost no cost in reprojection, and the ground folds where the lines meet.
// These factors make that cost something. docs/notes/level-flight.md.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <mutex>
#include <utility>
#include <vector>

#include "sfm/core/PriorSource.h"

namespace sfm {

struct LevelFlightOptions {
    // Sigma of the height difference between two cameras, per unit of their
    // distance over the ground: 0.02 lets two cameras 50 m apart differ by 1 m
    // at one sigma. Huber past that, so a real climb costs little more.
    double slope = 0.02;
    // Nearest cameras over the ground each camera is tied to, from other passes.
    int neighbours = 6;
    // Shooting order, sequence[image id] (empty: unknown). The images fix heights
    // along a line best, so a camera is tied only to cameras more than `min_gap`
    // photos away: the neighbouring lines, which nothing else holds in height.
    std::vector<int> sequence;
    int min_gap = 10;
    // The cameras' mean viewing direction must be at least this long (1 = all
    // look the same way): below it the capture is not one looking down, and
    // there is no up to hold heights along.
    double min_nadir = 0.8;
    int min_images = 5;
};

struct LevelFlightStats {
    int images = 0;
    int factors = 0;
    bool nadir = false;
};

// The world up of cameras that look down: against their mean viewing direction
// (camera +z, world -> camera rotations). False when they do not look one way.
inline bool levelFlightUp(const std::vector<PosedImage>& imgs, double min_nadir, Vec3& up) {
    Vec3 sum{0, 0, 0};
    for (const PosedImage& p : imgs) sum = sum + mul(transpose(p.pose.R), Vec3{0, 0, 1});
    if (imgs.empty() || sum.norm() < min_nadir * (double)imgs.size()) return false;
    up = sum.normalized() * -1.0;
    return true;
}

// One relative height factor per camera and each of its nearest neighbours over
// the ground, appended to `out`. Returns the number added.
inline int levelFlightFactors(const std::vector<PosedImage>& imgs, const LevelFlightOptions& opt,
                              PosePriors& out, LevelFlightStats* st = nullptr) {
    LevelFlightStats s;
    s.images = (int)imgs.size();
    Vec3 up;
    if ((int)imgs.size() < opt.min_images || !(opt.slope > 0) || opt.neighbours < 1 ||
        !levelFlightUp(imgs, opt.min_nadir, up)) {
        if (st) *st = s;
        return 0;
    }
    s.nadir = true;
    // A level frame: e1, e2 over the ground, up. Only the third row is weighted.
    const Vec3 any = std::fabs(up.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    const Vec3 e1 = up.cross(any).normalized();
    const Vec3 e2 = up.cross(e1);
    const Mat3 M{e1.x, e1.y, e1.z, e2.x, e2.y, e2.z, up.x, up.y, up.z};
    Mat3 negM;
    for (int k = 0; k < 9; k++) negM[k] = -M[k];

    const size_t n = imgs.size();
    std::vector<double> gx(n), gy(n);
    for (size_t i = 0; i < n; i++) {
        const Vec3 c = cameraCenter(imgs[i].pose);
        gx[i] = c.dot(e1);
        gy[i] = c.dot(e2);
    }
    auto seqOf = [&](size_t i) {
        const uint32_t id = imgs[i].image;
        return id < opt.sequence.size() ? opt.sequence[id] : -1;
    };
    // Direction of travel: the shorter step to the posed cameras before and after
    // in the sequence, which at a line's end is along its own line, not the jump
    // to the next; (0, 0) when it has neither.
    std::vector<std::pair<double, double>> dir(n, {0.0, 0.0});
    {
        std::vector<std::pair<int, size_t>> bySeq;
        for (size_t i = 0; i < n; i++)
            if (seqOf(i) >= 0) bySeq.push_back({seqOf(i), i});
        std::sort(bySeq.begin(), bySeq.end());
        for (size_t k = 0; k < bySeq.size(); k++) {
            const size_t i = bySeq[k].second;
            double best = 0, bx = 0, by = 0;
            for (int side = -1; side <= 1; side += 2) {
                const long m = (long)k + side;
                if (m < 0 || m >= (long)bySeq.size()) continue;
                const size_t j = bySeq[(size_t)m].second;
                const double dx = (gx[j] - gx[i]) * side, dy = (gy[j] - gy[i]) * side;
                const double d = std::hypot(dx, dy);
                if (d > 0 && (best == 0 || d < best)) {
                    best = d;
                    bx = dx / d;
                    by = dy / d;
                }
            }
            dir[i] = {bx, by};
        }
    }
    // Another pass: far enough apart in the sequence, and to the side of this
    // camera's track rather than ahead of or behind it on its own line.
    auto otherPass = [&](size_t i, size_t j) {
        const int a = seqOf(i), b = seqOf(j);
        if (a < 0 || b < 0) return true;
        if (std::abs(a - b) <= opt.min_gap) return false;
        const double dx = gx[j] - gx[i], dy = gy[j] - gy[i], d = std::hypot(dx, dy);
        const double side = std::fabs(dir[i].first * dy - dir[i].second * dx);
        return dir[i] == std::make_pair(0.0, 0.0) || side >= 0.5 * d;
    };
    // Nearest neighbours over the ground, both orderings of a pair kept once.
    std::vector<std::pair<size_t, size_t>> pairs;
    std::vector<double> dist;
    std::vector<double> nn(n, 0.0);
    std::vector<std::pair<double, size_t>> cand;
    for (size_t i = 0; i < n; i++) {
        cand.clear();
        for (size_t j = 0; j < n; j++)
            if (j != i && otherPass(i, j))
                cand.push_back({std::hypot(gx[j] - gx[i], gy[j] - gy[i]), j});
        const size_t k = std::min(cand.size(), (size_t)opt.neighbours);
        std::partial_sort(cand.begin(), cand.begin() + k, cand.end());
        nn[i] = k ? cand[0].first : 0.0;
        for (size_t a = 0; a < k; a++) {
            const size_t j = cand[a].second;
            pairs.push_back({std::min(i, j), std::max(i, j)});
            dist.push_back(cand[a].first);
        }
    }
    std::vector<size_t> order(pairs.size());
    for (size_t k = 0; k < order.size(); k++) order[k] = k;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return pairs[a] < pairs[b]; });
    // Cameras on top of each other would get an unbounded weight: half the
    // typical spacing is the least distance a sigma is taken from.
    std::vector<double> nz;
    for (double d : nn)
        if (d > 0) nz.push_back(d);
    double floor_d = 0;
    if (!nz.empty()) {
        std::nth_element(nz.begin(), nz.begin() + nz.size() / 2, nz.end());
        floor_d = 0.5 * nz[nz.size() / 2];
    }
    int added = 0;
    for (size_t k = 0; k < order.size(); k++) {
        const size_t o = order[k];
        if (k > 0 && pairs[o] == pairs[order[k - 1]]) continue;
        const double d = std::max(dist[o], floor_d);
        if (!(d > 0)) continue;
        PriorCentre c;
        c.n = 2;
        c.img[0] = imgs[pairs[o].first].image;
        c.img[1] = imgs[pairs[o].second].image;
        c.A[0] = M;
        c.A[1] = negM;
        c.b = {0, 0, 0};
        c.sigma = {0, 0, opt.slope * d};
        out.centres.push_back(c);
        added++;
    }
    s.factors = added;
    if (st) *st = s;
    return added;
}

// Another source's factors, if there is one, plus the level-flight ones; every
// other question goes to that source unchanged.
class LevelFlightPriors : public PriorSource {
public:
    LevelFlightPriors(PriorSource* inner, LevelFlightOptions opt) : inner_(inner), opt_(opt) {}

    const LevelFlightOptions& options() const { return opt_; }
    LevelFlightStats lastFactors() const {
        std::lock_guard<std::mutex> lk(mu_);
        return stats_;
    }

    bool has(uint32_t img) const override { return inner_ && inner_->has(img); }
    bool relativeRotation(uint32_t i, uint32_t j, Mat3& R_ji, double& sigma) const override {
        return inner_ && inner_->relativeRotation(i, j, R_ji, sigma);
    }
    std::vector<uint32_t> neighbours(uint32_t img) const override {
        return inner_ ? inner_->neighbours(img) : std::vector<uint32_t>{};
    }
    bool position(uint32_t img, Vec3& p) const override {
        return inner_ && inner_->position(img, p);
    }
    bool positionError(uint32_t img, const Pose& pose, const GpsFrame& f,
                       double& metres) const override {
        return inner_ && inner_->positionError(img, pose, f, metres);
    }
    bool declaredUp(uint32_t img, Vec3& u) const override {
        return inner_ && inner_->declaredUp(img, u);
    }
    void disableLevel() override {
        if (inner_) inner_->disableLevel();
    }

    PosePriors factors(const std::vector<PosedImage>& imgs) override {
        PosePriors out = inner_ ? inner_->factors(imgs) : PosePriors{};
        LevelFlightStats st;
        levelFlightFactors(imgs, opt_, out, &st);
        std::lock_guard<std::mutex> lk(mu_);
        stats_ = st;
        return out;
    }

private:
    PriorSource* inner_;
    LevelFlightOptions opt_;
    LevelFlightStats stats_;
    mutable std::mutex mu_;
};

}  // namespace sfm
