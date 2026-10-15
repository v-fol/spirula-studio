// Applying a rolling-shutter estimate (sfm/rs/RollingShutter.h) to a feature
// directory, and to the photos themselves. docs/notes/rolling-shutter.md.
#pragma once

#include <atomic>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <mutex>
#include <string>
#include <vector>

#include "external/stb_image.h"
#include "external/stb_image_write.h"
#include "sfm/core/Features.h"
#include "sfm/core/Matches.h"
#include "sfm/rs/RollingShutter.h"

namespace sfm::rs {

namespace fs = std::filesystem;

// Written beside the corrected feature files; its presence is what says a
// feature directory already holds corrected keypoints.
constexpr const char* kMarker = "rolling_shutter.txt";

inline bool alreadyApplied(const std::string& featdir) {
    std::error_code ec;
    return fs::exists(fs::path(featdir) / kMarker, ec);
}

// "photo9" before "photo10": file names of a capture are rarely zero-padded.
inline bool naturalLess(const std::string& a, const std::string& b) {
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (std::isdigit((unsigned char)a[i]) && std::isdigit((unsigned char)b[j])) {
            size_t i1 = i, j1 = j;
            while (i1 < a.size() && a[i1] == '0') i1++;
            while (j1 < b.size() && b[j1] == '0') j1++;
            size_t i2 = i1, j2 = j1;
            while (i2 < a.size() && std::isdigit((unsigned char)a[i2])) i2++;
            while (j2 < b.size() && std::isdigit((unsigned char)b[j2])) j2++;
            if (i2 - i1 != j2 - j1) return i2 - i1 < j2 - j1;
            const int c = a.compare(i1, i2 - i1, b, j1, j2 - j1);
            if (c != 0) return c < 0;
            i = i2;
            j = j2;
        } else {
            if (a[i] != b[j]) return a[i] < b[j];
            i++;
            j++;
        }
    }
    return a.size() - i < b.size() - j;
}

inline void jpegSink(void* ctx, void* data, int size) {
    auto* v = (std::vector<uint8_t>*)ctx;
    v->insert(v->end(), (uint8_t*)data, (uint8_t*)data + size);
}

// The APP1 segments (EXIF, XMP) of a JPEG file, verbatim.
inline std::vector<uint8_t> app1Segments(const std::vector<uint8_t>& f) {
    std::vector<uint8_t> out;
    if (f.size() < 4 || f[0] != 0xFF || f[1] != 0xD8) return out;
    size_t p = 2;
    while (p + 4 <= f.size() && f[p] == 0xFF) {
        const uint8_t m = f[p + 1];
        if (m == 0xDA || m == 0xD9) break;
        const size_t len = ((size_t)f[p + 2] << 8) | f[p + 3];
        if (p + 2 + len > f.size()) break;
        if (m == 0xE1) out.insert(out.end(), f.begin() + p, f.begin() + p + 2 + len);
        p += 2 + len;
    }
    return out;
}

inline std::vector<uint8_t> readFile(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return std::vector<uint8_t>((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

// Catmull-Rom; black outside the source.
inline void resample(const uint8_t* src, int w, int h, const ImageModel& m, uint8_t* dst) {
    const int L = 16, gw = w / L + 2, gh = h / L + 2;
    std::vector<float> gx((size_t)gw * gh), gy((size_t)gw * gh);
    for (int j = 0; j < gh; j++)
        for (int i = 0; i < gw; i++) {
            const Vec2 s = sourcePoint(m, {(double)i * L, (double)j * L});
            gx[(size_t)j * gw + i] = (float)s.x;
            gy[(size_t)j * gw + i] = (float)s.y;
        }
    auto cubic = [](float t, float wt[4]) {
        const float t2 = t * t, t3 = t2 * t;
        wt[0] = -0.5f * t3 + t2 - 0.5f * t;
        wt[1] = 1.5f * t3 - 2.5f * t2 + 1;
        wt[2] = -1.5f * t3 + 2 * t2 + 0.5f * t;
        wt[3] = 0.5f * t3 - 0.5f * t2;
    };
    for (int y = 0; y < h; y++) {
        const int j = y / L;
        const float fy = (float)(y - j * L) / L;
        for (int x = 0; x < w; x++) {
            const int i = x / L;
            const float fx = (float)(x - i * L) / L;
            const size_t a = (size_t)j * gw + i, b = a + gw;
            const float sx = (gx[a] * (1 - fx) + gx[a + 1] * fx) * (1 - fy) +
                             (gx[b] * (1 - fx) + gx[b + 1] * fx) * fy;
            const float sy = (gy[a] * (1 - fx) + gy[a + 1] * fx) * (1 - fy) +
                             (gy[b] * (1 - fx) + gy[b + 1] * fx) * fy;
            uint8_t* o = dst + ((size_t)y * w + x) * 3;
            if (!(sx >= 1 && sy >= 1 && sx < w - 2 && sy < h - 2)) {
                o[0] = o[1] = o[2] = 0;
                continue;
            }
            const int ix = (int)sx, iy = (int)sy;
            float wx[4], wy[4];
            cubic(sx - ix, wx);
            cubic(sy - iy, wy);
            float acc[3] = {0, 0, 0};
            for (int r = 0; r < 4; r++) {
                const uint8_t* row = src + ((size_t)(iy - 1 + r) * w + (ix - 1)) * 3;
                for (int c = 0; c < 4; c++) {
                    const float k = wy[r] * wx[c];
                    acc[0] += k * row[c * 3];
                    acc[1] += k * row[c * 3 + 1];
                    acc[2] += k * row[c * 3 + 2];
                }
            }
            for (int c = 0; c < 3; c++)
                o[c] = (uint8_t)std::min(255.f, std::max(0.f, acc[c] + 0.5f));
        }
    }
}

inline bool writeImage(const fs::path& src, const fs::path& dst, const ImageModel& m,
                       int quality) {
    if (!m.solved) {
        std::error_code ec;
        fs::copy_file(src, dst, fs::copy_options::overwrite_existing, ec);
        return !ec;
    }
    int w = 0, h = 0, c = 0;
    uint8_t* px = stbi_load(src.string().c_str(), &w, &h, &c, 3);
    if (!px) return false;
    bool ok = w == m.width && h == m.height;
    if (ok) {
        std::vector<uint8_t> out((size_t)w * h * 3), jpg;
        resample(px, w, h, m, out.data());
        ok = stbi_write_jpg_to_func(jpegSink, &jpg, w, h, 3, out.data(), quality) != 0;
        if (ok) {
            const std::vector<uint8_t> app1 = app1Segments(readFile(src));
            // After stb's JFIF APP0 when it wrote one, else straight after SOI.
            size_t at = 2;
            if (jpg.size() > 6 && jpg[2] == 0xFF && jpg[3] == 0xE0)
                at = 4 + (((size_t)jpg[4] << 8) | jpg[5]);
            std::ofstream f(dst, std::ios::binary);
            f.write((const char*)jpg.data(), (std::streamsize)at);
            f.write((const char*)app1.data(), (std::streamsize)app1.size());
            f.write((const char*)jpg.data() + at, (std::streamsize)(jpg.size() - at));
            ok = (bool)f;
        }
    }
    stbi_image_free(px);
    return ok;
}


// `db` holds unverified matches over `featdir`. Feature files are written to
// `outdir`, which may be `featdir` itself; `write_images` empty writes no photos.
inline bool applyToFeatureDir(const MatchesDatabase& db, const std::string& featdir,
                              const std::string& outdir, const Options& opt,
                              const std::string& imagedir, const std::string& write_images,
                              int quality, FILE* log, std::string& err) {
    const auto t0 = std::chrono::steady_clock::now();
    auto secs = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
        const int n = (int)db.images.size();
    std::vector<int> order(n), seq(n);
    for (int i = 0; i < n; i++) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return naturalLess(db.images[a].name, db.images[b].name); });
    for (int s = 0; s < n; s++) seq[order[s]] = s;

    std::vector<FeatureSet> feats(n);
    std::mutex err_mtx;
    std::string first_error;
    detail::parallelFor(n, opt.threads, [&](int i) {
        try {
            feats[i] = readFeatures((fs::path(featdir) / (db.images[i].name + ".bin")).string(), false);
        } catch (const std::exception& e) {
            std::lock_guard<std::mutex> lk(err_mtx);
            if (first_error.empty()) first_error = e.what();
        }
    });
    if (!first_error.empty()) { err = first_error; return false; }
    const int width = feats[0].width, height = feats[0].height;
    for (const FeatureSet& f : feats)
        if (f.width != width || f.height != height)
            { err = "the photos are not all one size; a sequence from one camera is expected"; return false; }

    std::vector<int> cand;
    for (int k = 0; k < (int)db.pairs.size(); k++) {
        const int d = std::abs(seq[db.pairs[k].image1] - seq[db.pairs[k].image2]);
        if (d >= 1 && d <= opt.max_gap) cand.push_back(k);
    }
    std::vector<PairMatches> pairs(cand.size());
    detail::parallelFor((int)cand.size(), opt.threads, [&](int c) {
        const TwoViewMatches& tv = db.pairs[cand[c]];
        uint32_t a = tv.image1, b = tv.image2;
        const bool swap = seq[a] > seq[b];
        if (swap) std::swap(a, b);
        std::vector<Vec2> p, q;
        p.reserve(tv.matches.size());
        q.reserve(tv.matches.size());
        for (const FeatureMatch& m : tv.matches) {
            const Keypoint& ka = feats[a].keypoints[swap ? m.idx2 : m.idx1];
            const Keypoint& kb = feats[b].keypoints[swap ? m.idx1 : m.idx2];
            p.push_back({ka.x, ka.y});
            q.push_back({kb.x, kb.y});
        }
        pairs[c] = filterPair(seq[a], seq[b], p, q, opt);
    });
    size_t kept = 0;
    for (const PairMatches& pm : pairs) kept += !pm.p.empty();
    if (log) std::fprintf(log, "[rs]      %d images, %zu neighbour pairs of %zu matched pairs hold together (%.1f s)\n",
                n, kept, db.pairs.size(), secs());

    const Result res = estimate(n, width, height, pairs, opt);
    if (!res.scan.empty()) {
        if (log) std::fprintf(log, "[rs]      Focal scan (px: fit, 1.00 is best):");
        for (const auto& s : res.scan) if (log) std::fprintf(log, " %.0f: %.2f", s.first, s.second);
        if (log) std::fprintf(log, "\n");
    }
    std::vector<double> fit;
    for (const WindowReport& w : res.windows) fit.push_back(w.neighbours_px);
    std::sort(fit.begin(), fit.end());
    if (log) std::fprintf(log, "[rs]      Focal %.0f px; %zu windows, neighbours agree to %.2f px (median), %.2f (worst); "
                "%zu breaks in the sequence, %d images left uncorrected (%.1f s)\n",
                res.focal, res.windows.size(), fit.empty() ? 0.0 : fit[fit.size() / 2],
                fit.empty() ? 0.0 : fit.back(), res.breaks.size(), res.uncorrected, secs());

    std::error_code ec;
    fs::create_directories(outdir, ec);
    detail::parallelFor(n, opt.threads, [&](int i) {
        try {
            const fs::path rel = db.images[i].name + ".bin";
            FeatureSet fsx = readFeatures((fs::path(featdir) / rel).string(), true);
            const ImageModel& m = res.images[seq[i]];
            for (Keypoint& k : fsx.keypoints) {
                const Vec2 c = correctPoint(m, {k.x, k.y});
                k.x = (float)c.x;
                k.y = (float)c.y;
            }
            const fs::path dst = fs::path(outdir) / rel;
            std::error_code e2;
            fs::create_directories(dst.parent_path(), e2);
            writeFeatures(dst.string(), fsx);
        } catch (const std::exception& e) {
            std::lock_guard<std::mutex> lk(err_mtx);
            if (first_error.empty()) first_error = e.what();
        }
    });
    if (!first_error.empty()) { err = first_error; return false; }
    {
        std::ofstream f(fs::path(outdir) / kMarker);
        f << "# name solved focal_px k1 then 9 read-out rotations (rad over the frame): "
             "order 1 x y z, order 2 x y z, order 3 x y z\n";
        f.precision(9);
        for (int s = 0; s < n; s++) {
            const ImageModel& m = res.images[s];
            f << db.images[order[s]].name << ' ' << (m.solved ? 1 : 0) << ' ' << m.f << ' ' << m.k1;
            for (double v : m.w) f << ' ' << v;
            f << '\n';
        }
    }
    if (log) std::fprintf(log, "[rs]      Wrote %s (%.1f s)\n", outdir.c_str(), secs());

    if (!write_images.empty()) {
        std::map<std::string, fs::path> by_name;
        for (auto it = fs::recursive_directory_iterator(imagedir, ec);
             !ec && it != fs::recursive_directory_iterator(); it.increment(ec)) {
            if (!it->is_regular_file()) continue;
            fs::path rel = it->path().lexically_relative(imagedir);
            by_name[rel.replace_extension().generic_string()] = it->path();
        }
        std::atomic<int> written{0}, missing{0};
        detail::parallelFor(n, opt.threads, [&](int i) {
            auto it = by_name.find(fs::path(db.images[i].name).generic_string());
            if (it == by_name.end()) { missing++; return; }
            const fs::path dst = fs::path(write_images) / it->second.lexically_relative(imagedir);
            std::error_code e2;
            fs::create_directories(dst.parent_path(), e2);
            if (writeImage(it->second, dst, res.images[seq[i]], quality)) written++;
            else missing++;
        });
        if (log) std::fprintf(log, "[rs]      Wrote %d corrected photos to %s, %d could not be written (%.1f s)\n",
                    written.load(), write_images.c_str(), missing.load(), secs());
    }
    return true;
}

}  // namespace sfm::rs
