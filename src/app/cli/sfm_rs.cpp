// `spirula-sfm rs` -- rolling-shutter estimation between `match` and `map`:
// keypoints moved to where an instant exposure would have put them, and
// optionally the photos resampled the same way. docs/notes/rolling-shutter.md.
#include <cstdio>
#include <cstring>
#include <string>

#include "sfm/rs/Apply.h"

using namespace sfm;

static const char* kRsHelp =
    "spirula-sfm rs -- estimate and remove rolling-shutter distortion\n\n"
    "Usage:\n  spirula-sfm rs <MATCHES.BIN> <FEATURE_DIR> -o OUT_FEATURE_DIR [options]\n\n"
    "For photo sequences taken roughly straight down from a moving aircraft with a\n"
    "rolling-shutter sensor. Reads unverified matches (`match --no-verify`), works out\n"
    "how each camera rotated while its sensor was read, and writes the feature files\n"
    "again with the keypoints moved to where an instant exposure would have put them.\n"
    "Run `match` and `map` on OUT_FEATURE_DIR afterwards.\n\n"
    "Options:\n"
    "  -o, --output DIR        corrected feature files (required)\n"
    "  --focal X               focal length in pixels; 0 scans for it          [0]\n"
    "  --images DIR            the photos, for --write-images\n"
    "  --write-images DIR      also write the photos resampled without the distortion\n"
    "  --jpeg-quality N                                                        [95]\n"
    "  --window N              photos solved together                          [60]\n"
    "  --overlap N             photos shared by consecutive windows            [12]\n"
    "  --max-gap N             furthest sequence neighbour used                [5]\n"
    "  --threads N             0 uses every core                               [0]\n"
    "  -h, --help\n";

void printRsHelp(FILE* out) { std::fputs(kRsHelp, out); }

int cmdRs(int argc, char** argv) {
    rs::Options opt;
    std::string matches, featdir, outdir, imagedir, write_images;
    int quality = 95;
    auto fail = [](const std::string& msg) {
        std::fprintf(stderr, "spirula-sfm rs: error: %s\n", msg.c_str());
        return 1;
    };
    for (int i = 0; i < argc; i++) {
        const std::string a = argv[i];
        auto value = [&](std::string& v) {
            if (i + 1 >= argc) return false;
            v = argv[++i];
            return true;
        };
        std::string v;
        if (a == "--help" || a == "-h") { printRsHelp(stdout); return 0; }
        if (a == "-o" || a == "--output") { if (!value(outdir)) return fail(a + ": missing value"); }
        else if (a == "--images") { if (!value(imagedir)) return fail(a + ": missing value"); }
        else if (a == "--write-images") { if (!value(write_images)) return fail(a + ": missing value"); }
        else if (a == "--focal") { if (!value(v)) return fail(a + ": missing value"); opt.focal = std::atof(v.c_str()); }
        else if (a == "--jpeg-quality") { if (!value(v)) return fail(a + ": missing value"); quality = std::atoi(v.c_str()); }
        else if (a == "--window") { if (!value(v)) return fail(a + ": missing value"); opt.window = std::atoi(v.c_str()); }
        else if (a == "--overlap") { if (!value(v)) return fail(a + ": missing value"); opt.overlap = std::atoi(v.c_str()); }
        else if (a == "--max-gap") { if (!value(v)) return fail(a + ": missing value"); opt.max_gap = std::atoi(v.c_str()); }
        else if (a == "--threads") { if (!value(v)) return fail(a + ": missing value"); opt.threads = std::atoi(v.c_str()); }
        else if (a[0] == '-') return fail("unknown option " + a);
        else if (matches.empty()) matches = a;
        else if (featdir.empty()) featdir = a;
        else return fail("unexpected argument '" + a + "'");
    }
    if (matches.empty() || featdir.empty()) return fail("a matches file and a feature directory are required");
    if (outdir.empty()) return fail("-o OUT_FEATURE_DIR is required");
    if (!write_images.empty() && imagedir.empty()) return fail("--write-images needs --images");
    if (opt.window < 2 * opt.overlap + opt.min_run) return fail("--window is too small for --overlap");

    std::string err;
    if (!rs::applyToFeatureDir(readMatches(matches), featdir, outdir, opt, imagedir, write_images,
                               quality, stdout, err))
        return fail(err);
    return 0;
}
