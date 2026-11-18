// rollingstone3d_full.cpp —— 三维滚石模拟（独立单文件，作业要求完整三维版）
//
// 与二维版的关系：核心模拟为真三维高度场物理（球可侧向偏移、沿沟谷汇聚），
// 作业要求的功能全部保留并以三维方式实现：
//   1) 批量随机试验 + 坡脚 0–10 m 路线概率（另加三维落点横向分布）
//   2) 统计结果 CSV 导出（高度/横向概率、逐次指标、轨迹、配置、摘要、口径说明）
//   3) 参数影响：初速与地表四系数正态抽样 CV=0.1；--study 做 5 因子×0.8/1.0/1.2 配对研究
//   4) 地形影响：随机地形生成器（坡段/起伏/沟谷/谷宽，种子可复现）；--terrain-study 多场景配对
//   *）刚性挡板拦截（有限高度，可越顶）与弹性挡板变形阻拦（弹簧-阻尼-等效质量，超位移失效）
// 形状说明：三维接触采用球体；凸多边形要求由二维版满足（三维多面体接触超出教学模型范围）。
//
// 编译：cl /std:c++17 /EHsc /utf-8 rollingstone3d_full.cpp /link EasyXw.lib shell32.lib
// 运行：无参数=交互界面；--batch/--study/--terrain-study=无界面统计；--shot out.png [seed]=自动演示截图
// 界面按键：Space 暂停 R 重置 T 尾迹 N 换地形 [ ] 坡段 - = 起伏 V 沟谷 W 谷宽
//           1/2/3 挡板 无/刚/弹  4/5 挡板位置  6/7 挡板高度  , . 倍速
//           B 批量  E 导出  A 参数研究  G 地形研究  Esc 退出
// 单位 SI（m,s,kg,N,J）。教学模型，未经工程标定。

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef UNICODE
#define UNICODE
#endif
#ifndef _UNICODE
#define _UNICODE
#endif
#include <easyx.h>
#include <windows.h>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>
#include <chrono>
#include <future>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <map>

#pragma comment(lib, "user32.lib")
#pragma comment(lib, "gdi32.lib")
#pragma comment(lib, "shell32.lib")

namespace r3f {

constexpr double G = 9.81;
constexpr double PI = 3.14159265358979323846;
constexpr int W = 1280, H = 800;

// ---------------- 随机三维地形 ----------------
struct Terrain {
    int segments = 10;
    double roughness = 0.6;
    int valleys = 2;
    double valleyWidth = 1.0;
    double length = 120.0, width = 60.0;
    double lengthScale = 1.0, dropScale = 1.0;

    uint32_t seedUsed = 0;
    double drop = 52.0, profileExp = 0.82;
    std::vector<double> nodeZ, nodeJitter;
    struct Valley { double y0, amp, freq, phase, depth, width; };
    std::vector<Valley> vlist;
    double bumpPhase = 0, bumpAmp = 1.0;

    static double lerp(double a, double b, double t) { return a + (b - a) * t; }

    void generate(uint32_t seed) {
        seedUsed = seed;
        uint32_t s = seed * 2654435761u + 12345u;
        auto rnd = [&]() { s = s * 1664525u + 1013904223u; return (s >> 8) / 16777216.0; };
        segments = std::clamp(segments, 3, 24);
        roughness = std::clamp(roughness, 0.0, 1.0);
        valleys = std::clamp(valleys, 0, 3);
        valleyWidth = std::clamp(valleyWidth, 0.4, 2.0);
        drop = (40.0 + rnd() * 24.0) * dropScale;
        profileExp = 0.62 + rnd() * 0.55;
        bumpPhase = rnd() * 6.2831853;
        bumpAmp = 0.6 + rnd() * 1.1;
        nodeZ.assign(segments + 1, 0.0);
        nodeJitter.assign(segments + 1, 0.0);
        for (int i = 0; i <= segments; ++i) {
            const double t = (double)i / segments;
            double z = drop * (1.0 - std::pow(t, profileExp));
            if (i > 0 && i < segments) {
                z += (rnd() - 0.5) * drop * 0.16 * roughness;
                nodeJitter[i] = (rnd() - 0.5) * 6.0 * roughness;
            }
            nodeZ[i] = z;
        }
        nodeZ[0] = drop; nodeZ[segments] = 0.0;
        for (int i = 1; i <= segments; ++i)
            nodeZ[i] = std::min(nodeZ[i], nodeZ[i - 1] - drop * 0.012);
        vlist.clear();
        for (int k = 0; k < valleys; ++k) {
            Valley v{};
            v.y0 = (rnd() - 0.5) * width * 0.55;
            v.amp = 3.0 + rnd() * 8.0;
            v.freq = 0.03 + rnd() * 0.06;
            v.phase = rnd() * 6.2831853;
            v.depth = (0.006 + rnd() * 0.010) * (k == 0 ? 1.25 : 0.8);
            v.width = (9.0 + rnd() * 9.0) * valleyWidth;
            vlist.push_back(v);
        }
        length = 120.0 * lengthScale;
    }

    double baseZ(double x, double& tilt) const {
        const double t = std::clamp(x / length, 0.0, 1.0) * segments;
        int i = (int)t; if (i >= segments) i = segments - 1;
        const double f = t - i;
        tilt = lerp(nodeJitter[i], nodeJitter[i + 1], f);
        return lerp(nodeZ[i], nodeZ[i + 1], f);
    }
    double height(double x, double y) const {
        double tilt = 0;
        double z = baseZ(x, tilt);
        const double t = std::clamp(x / length, 0.0, 1.0);
        for (const auto& v : vlist) {
            const double c = v.y0 + v.amp * std::sin(x * v.freq + v.phase);
            const double d = (y - c) / v.width;
            z += v.depth * v.width * v.width * d * d * 0.5;
        }
        z += tilt * (y / (width * 0.5)) * 0.35;
        const double bumps = bumpAmp * roughness *
            (1.7 * std::sin(x * 0.13 + bumpPhase) * std::cos(y * 0.16)
           + 0.9 * std::sin(x * 0.31 + 1.7) * std::sin(y * 0.27 + bumpPhase));
        const double flat = (t > 0.9) ? (t - 0.9) * 10.0 : 0.0;
        return z * (1.0 - 0.45 * flat) + bumps * (1.0 - 0.6 * flat);
    }
    void grad(double x, double y, double& gx, double& gy) const {
        const double e = 0.05;
        gx = (height(x + e, y) - height(x - e, y)) / (2 * e);
        gy = (height(x, y + e) - height(x, y - e)) / (2 * e);
    }
    double startY() const {
        if (vlist.empty()) return width * 0.18;
        const auto& v = vlist.front();
        const double c = v.y0 + v.amp * std::sin(5.0 * v.freq + v.phase);
        return std::clamp(c + v.width * 0.85, -width * 0.42, width * 0.42);
    }
};

// ---------------- 挡板 ----------------
enum class BarrierKind { None, Rigid, Elastic };
struct Barrier {
    BarrierKind kind = BarrierKind::None;
    double x = 100.0, height = 6.0, restitution = 0.15;
    double stiffness = 2e6, damping = 2e5, mass = 1000.0, maxDeflection = 3.0;
};
struct BarrierState {
    double d = 0, vd = 0, diss = 0, peakForce = 0, maxDef = 0;
    bool failed = false, hit = false, passed = false, intercepted = false;
};

// ---------------- 配置 ----------------
struct Config {
    Terrain terrain;
    Barrier barrier;
    uint32_t terrainSeed = 7;
    double initialSpeed = 2.0, radius = 1.6, density = 2600.0;
    double rollRes = 0.02, friction = 0.45, restitution = 0.42, drag = 1.0;
    double cv = 0.1, dt = 0.002, maxTime = 90.0, sampleInterval = 0.05;
    double detectorX = -1.0, binHeight = 0.5, yBin = 2.0;
    int runs = 200;
    uint32_t seed = 42;
    void finalize() { terrain.generate(terrainSeed); if (detectorX < 0) detectorX = terrain.length; }
};

// ---------------- 球与物理 ----------------
struct Sample3 { double t, x, y, z, vx, vy, vz, speed, bd; };
struct Ball {
    double x = 5, y = 0, z = 0, vx = 0, vy = 0, vz = 0;
    double radius = 1.6, mass = 1.0;
    double rollRes = 0.02, friction = 0.45, rest = 0.42, drag = 1.0;
    double spin = 0, t = 0, restTime = 0, peakSpeed = 0, distance = 0, maxBounce = 0;
    bool resting = false;
    BarrierState bs;

    void place(const Terrain& T, const Config& c, double y0) {
        radius = c.radius;
        mass = c.density * 4.0 / 3.0 * PI * radius * radius * radius;
        rollRes = c.rollRes; friction = c.friction; rest = c.restitution; drag = c.drag;
        x = 5.0; y = y0; z = T.height(x, y) + radius;
        vx = c.initialSpeed; vy = 0; vz = 0;
        spin = t = restTime = peakSpeed = distance = maxBounce = 0;
        resting = false; bs = BarrierState{};
    }
    double speed() const { return std::sqrt(vx * vx + vy * vy + vz * vz); }
    double kinetic() const { return 0.5 * mass * speed() * speed() * 1.4; }

    void step(const Terrain& T, const Barrier& B, double dt) {
        if (resting) { t += dt; return; }
        const double px = x, py = y;
        vz -= G * dt;
        const double sp = speed();
        if (sp > 1e-9 && drag > 0) {
            const double kd = 0.5 * 1.225 * drag * PI * radius * radius / mass;
            const double k = 1.0 / (1.0 + kd * sp * dt);
            vx *= k; vy *= k; vz *= k;
        }
        x += vx * dt; y += vy * dt; z += vz * dt;
        if (y < -T.width * 0.5 + radius) { y = -T.width * 0.5 + radius; vy = std::abs(vy) * 0.3; }
        if (y >  T.width * 0.5 - radius) { y =  T.width * 0.5 - radius; vy = -std::abs(vy) * 0.3; }

        // ---- 地面接触 ----
        const double h = T.height(x, y);
        if (z - radius - h <= 0.0) {
            double gx, gy; T.grad(x, y, gx, gy);
            double nx = -gx, ny = -gy, nz = 1.0;
            const double nl = std::sqrt(nx * nx + ny * ny + nz * nz);
            nx /= nl; ny /= nl; nz /= nl;
            z = h + radius;
            const double vn = vx * nx + vy * ny + vz * nz;
            const bool impact = (vn < -0.6);
            if (vn < 0) {
                const double e = impact ? rest : 0.0;
                const double jn = -(1.0 + e) * vn;
                vx += jn * nx; vy += jn * ny; vz += jn * nz;
            }
            if (impact) {
                const double vn2 = vx * nx + vy * ny + vz * nz;
                double tvx = vx - vn2 * nx, tvy = vy - vn2 * ny, tvz = vz - vn2 * nz;
                const double tv = std::sqrt(tvx * tvx + tvy * tvy + tvz * tvz);
                if (tv > 1e-9) {
                    const double jt = std::min(tv * (2.0 / 7.0), friction * std::abs(vn) * (1.0 + rest));
                    vx -= jt * tvx / tv; vy -= jt * tvy / tv; vz -= jt * tvz / tv;
                }
            } else {
                // 纯滚动约束只在有切向滑动时施加；静置时不锁速度，让重力自然驱动
                const double vn0 = vx * nx + vy * ny + vz * nz;
                const double tvx = vx - vn0 * nx, tvy = vy - vn0 * ny, tvz = vz - vn0 * nz;
                const double tv = std::sqrt(tvx * tvx + tvy * tvy + tvz * tvz);
                if (tv > 0.05) {
                    const double gtn = (-G) * nz;
                    const double gtx = -gtn * nx, gty = -gtn * ny, gtz = -G - gtn * nz;
                    vx -= (2.0 / 7.0) * gtx * dt; vy -= (2.0 / 7.0) * gty * dt; vz -= (2.0 / 7.0) * gtz * dt;
                    if (std::sqrt(gx * gx + gy * gy) > friction) {
                        vx += (2.0 / 7.0) * gtx * dt; vy += (2.0 / 7.0) * gty * dt; vz += (2.0 / 7.0) * gtz * dt;
                    }
                }
            }
            const double sp2 = speed();
            if (sp2 > 1e-9) {
                const double k = std::max(0.0, 1.0 - rollRes * G * dt / sp2);
                vx *= k; vy *= k; vz *= k;
            }
            spin += sp2 * dt / radius;
            const double slope = std::sqrt(gx * gx + gy * gy);
            // 接触且几乎静止即累计停稳（含贴墙/平台情形），避免数值平衡悬停
            if (sp2 < 0.15 && (slope <= rollRes + 0.02 || sp2 < 0.05)) {
                restTime += dt;
                if (restTime > 0.6) { resting = true; vx = vy = vz = 0; }
            } else restTime = 0;
        }
        maxBounce = std::max(maxBounce, z - radius - h);

        // ---- 挡板接触（竖直墙，横跨整个坡宽，有限高度） ----
        if (B.kind != BarrierKind::None && !bs.failed) {
            const bool elastic = (B.kind == BarrierKind::Elastic);
            // 弹性板自身动力学（隐式欧拉）
            if (elastic) {
                const double den = 1 + dt * B.damping / B.mass + dt * dt * B.stiffness / B.mass;
                bs.vd = (bs.vd - dt * B.stiffness * bs.d / B.mass) / den;
                bs.d += bs.vd * dt;
                bs.diss += dt * B.damping * bs.vd * bs.vd;
            }
            const double wallX = B.x + (elastic ? bs.d : 0.0);
            const double top = T.height(B.x, y) + B.height;
            if (x > wallX + radius) bs.passed = true;
            const double pen = radius - (wallX - x);
            if (!bs.passed && pen > 0 && (z - radius) < top) {
                bs.hit = true;
                if (elastic) {
                    const double vrel = vx - bs.vd;
                    if (vrel > 0) {
                        // 接触阻尼（临界阻尼的 0.7 倍）耗散法向微振，保证贴板收敛停稳
                        const double meff = 1.0 / (1.0 / mass + 1.0 / B.mass);
                        const double cd = 2.0 * std::sqrt(B.stiffness * meff) * 0.7;
                        const double jj = meff * vrel * std::max(0.0, 1.0 - cd * dt / meff);
                        vx -= jj / mass;
                        bs.vd += jj / B.mass;
                    }
                    x = wallX - radius;
                    const double F = B.stiffness * bs.d + B.damping * bs.vd;
                    bs.peakForce = std::max(bs.peakForce, std::abs(F));
                    bs.maxDef = std::max(bs.maxDef, std::abs(bs.d));
                    if (std::abs(bs.d) > B.maxDeflection) { bs.failed = true; bs.diss += 0.5 * B.stiffness * bs.d * bs.d + 0.5 * B.mass * bs.vd * bs.vd; bs.vd = 0; }
                } else {
                    if (vx > 0) {
                        const double e = (vx > 0.6) ? B.restitution : 0.0;
                        const double jn = (1.0 + e) * vx;
                        vx -= jn;
                        const double tv = std::sqrt(vy * vy + vz * vz);
                        if (tv > 1e-9) {
                            const double jt = std::min(tv * (2.0 / 7.0), 0.6 * jn);
                            vy -= jt * vy / tv; vz -= jt * vz / tv;
                        }
                    }
                    x = wallX - radius;
                }
            }
        }
        // 倚靠挡板的停稳：低速且持续接触 0.6 s 即判截停（不要求坡面摩擦角条件）
        if (B.kind != BarrierKind::None && bs.hit && !bs.failed && !bs.passed) {
            const double wallX = B.x + (B.kind == BarrierKind::Elastic ? bs.d : 0.0);
            if (speed() < 0.3 && (wallX - x) < radius + 0.2) {
                restTime += dt;
                if (restTime > 0.6) { resting = true; vx = vy = vz = 0; }
            }
        }
        t += dt;
        peakSpeed = std::max(peakSpeed, speed());
        distance += std::hypot(x - px, y - py);
        if (x > T.length + 6 || t > 90.0) { resting = true; vx = vy = vz = 0; }
        if (resting && bs.hit && !bs.failed && !bs.passed) bs.intercepted = true;
    }
};

// ---------------- 单次结果与批量 ----------------
struct RunResult {
    bool reached = false, inBand = false;
    double detHeight = 0, detY = 0, detSpeed = 0, detEnergy = 0, detTime = 0;
    double runoutX = 0, peakSpeed = 0, maxBounce = 0, duration = 0;
    double initSpeed = 0, rollRes = 0, friction = 0, rest = 0, drag = 0;
    std::string status, outcome;
    BarrierState bs;
    std::vector<Sample3> traj;
};
struct Batch {
    Config cfg;
    std::vector<RunResult> runs;
    std::vector<int> heightCounts, yCounts;
    size_t reached = 0, inBand = 0, above = 0, below = 0, notReached = 0;
    size_t intercepted = 0, overtopped = 0, failedPass = 0, hitUnresolved = 0;
};

inline double normalBounded(double mean, double cv, double lo, double hi, std::mt19937& rng) {
    if (cv == 0 || mean == 0) return mean;
    std::normal_distribution<double> nd(mean, std::abs(mean) * cv);
    for (int i = 0; i < 10000; ++i) {
        const double v = nd(rng);
        if (std::isfinite(v) && v >= lo && v <= hi) return v;
    }
    return mean;
}

inline std::string outcomeOf(const RunResult& r, BarrierKind k) {
    if (k == BarrierKind::None) return "disabled";
    if (r.bs.intercepted) return "intercepted";
    if (r.bs.failed && r.bs.passed) return "failed_passage";
    if (r.bs.passed) return "overtopped";
    if (r.bs.failed) return "failed_unresolved";
    if (r.bs.hit) return "hit_unresolved";
    return "not_encountered";
}

inline RunResult simulate(const Config& c, uint32_t trialSeed, bool keepTraj) {
    std::mt19937 rng(trialSeed);
    Config cc = c;
    cc.initialSpeed = normalBounded(c.initialSpeed, c.cv, 0, 60, rng);
    cc.rollRes = normalBounded(c.rollRes, c.cv, 0, 1, rng);
    cc.friction = normalBounded(c.friction, c.cv, 0, 2, rng);
    cc.restitution = normalBounded(c.restitution, c.cv, 0, 1, rng);
    cc.drag = normalBounded(c.drag, c.cv, 0, 10, rng);
    Ball b; b.place(c.terrain, cc, c.terrain.startY());
    RunResult r;
    r.initSpeed = cc.initialSpeed; r.rollRes = cc.rollRes; r.friction = cc.friction;
    r.rest = cc.restitution; r.drag = cc.drag;
    const double detX = c.detectorX;
    double nextSample = 0;
    const int sub = (c.barrier.kind == BarrierKind::Elastic) ? 4 : 1;
    const double h = c.dt / sub;
    const int steps = (int)std::ceil(c.maxTime / h);
    for (int i = 0; i < steps && !b.resting; ++i) {
        const double px = b.x, pz = b.z, py = b.y, pt = b.t;
        const double pvx = b.vx, pvy = b.vy, pvz = b.vz;
        b.step(c.terrain, c.barrier, h);
        r.runoutX = std::max(r.runoutX, b.x);
        if (!r.reached && px < detX && b.x >= detX) {
            const double f = (detX - px) / std::max(1e-12, b.x - px);
            r.reached = true;
            const double zi = pz + (b.z - pz) * f, yi = py + (b.y - py) * f;
            r.detHeight = zi - c.terrain.height(detX, yi);
            r.detY = yi;
            const double sx = pvx + (b.vx - pvx) * f, sy = pvy + (b.vy - pvy) * f, sz = pvz + (b.vz - pvz) * f;
            r.detSpeed = std::sqrt(sx * sx + sy * sy + sz * sz);
            r.detEnergy = 0.5 * b.mass * r.detSpeed * r.detSpeed * 1.4;
            r.detTime = pt + h * f;
            r.inBand = (r.detHeight >= 0 && r.detHeight <= 10);
        }
        if (keepTraj && b.t >= nextSample) {
            r.traj.push_back({ b.t, b.x, b.y, b.z, b.vx, b.vy, b.vz, b.speed(), b.bs.d });
            nextSample += c.sampleInterval;
        }
    }
    r.peakSpeed = b.peakSpeed; r.maxBounce = b.maxBounce; r.duration = b.t;
    r.bs = b.bs;
    r.status = b.resting ? (b.bs.intercepted ? "intercepted" : "stopped") : "time_limit";
    if (b.x > c.terrain.length + 5) r.status = "domain_exit";
    r.outcome = outcomeOf(r, c.barrier.kind);
    return r;
}

inline std::pair<double, double> wilson(size_t k, size_t n) {
    if (!n) return { 0, 0 };
    const double z = 1.95996398454, p = (double)k / n, d = 1 + z * z / n;
    const double c = (p + z * z / (2 * n)) / d;
    const double h = z * std::sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d;
    return { std::max(0.0, c - h), std::min(1.0, c + h) };
}

inline Batch runBatch(const Config& c, int count, uint32_t seed, bool keepTraj) {
    Batch b; b.cfg = c;
    const int nb = (int)std::ceil(10.0 / c.binHeight);
    b.heightCounts.assign(nb, 0);
    const int ny = (int)std::ceil(c.terrain.width / c.yBin);
    b.yCounts.assign(ny, 0);
    for (int i = 0; i < count; ++i) {
        std::seed_seq ss{ seed, (uint32_t)i };
        std::mt19937 rng(ss);
        RunResult r = simulate(c, (uint32_t)rng(), keepTraj);
        if (r.reached) {
            ++b.reached;
            if (r.inBand) { ++b.inBand; int k = (int)(r.detHeight / c.binHeight); k = std::clamp(k, 0, nb - 1); ++b.heightCounts[k]; }
            else if (r.detHeight > 10) ++b.above; else ++b.below;
            int ky = (int)((r.detY + c.terrain.width * 0.5) / c.yBin); ky = std::clamp(ky, 0, ny - 1); ++b.yCounts[ky];
        } else ++b.notReached;
        if (r.outcome == "intercepted") ++b.intercepted;
        else if (r.outcome == "overtopped") ++b.overtopped;
        else if (r.outcome == "failed_passage") ++b.failedPass;
        else if (r.outcome == "hit_unresolved" || r.outcome == "failed_unresolved") ++b.hitUnresolved;
        b.runs.push_back(std::move(r));
    }
    return b;
}

// ---------------- CSV 导出 ----------------
namespace fs = std::filesystem;
inline std::ofstream csv(const fs::path& p) {
    std::ofstream f(p, std::ios::binary);
    f.imbue(std::locale::classic());
    f << std::setprecision(10);
    return f;
}
inline void exportBatch(const Batch& b, const std::string& dir) {
    fs::path out(dir);
    if (fs::exists(out) && !fs::is_empty(out)) throw std::runtime_error("Output folder must be new or empty: " + dir);
    fs::create_directories(out);
    const Config& c = b.cfg;
    const size_t n = b.runs.size();
    auto hp = csv(out / "height_probability.csv");
    hp << "height_lower_m,height_upper_m,count,total,probability,wilson95_low,wilson95_high\n";
    for (size_t i = 0; i < b.heightCounts.size(); ++i) {
        auto ci = wilson(b.heightCounts[i], n);
        hp << i * c.binHeight << ',' << std::min(10.0, (i + 1) * c.binHeight) << ',' << b.heightCounts[i]
           << ',' << n << ',' << (double)b.heightCounts[i] / n << ',' << ci.first << ',' << ci.second << '\n';
    }
    auto yp = csv(out / "lateral_probability.csv");
    yp << "y_lower_m,y_upper_m,count,total,probability\n";
    for (size_t i = 0; i < b.yCounts.size(); ++i)
        yp << -c.terrain.width * 0.5 + i * c.yBin << ',' << -c.terrain.width * 0.5 + (i + 1) * c.yBin
           << ',' << b.yCounts[i] << ',' << n << ',' << (double)b.yCounts[i] / n << '\n';
    auto rr = csv(out / "runs.csv");
    rr << "trial,init_speed,roll_res,friction,restitution,drag,reached,in_band,det_height_m,det_y_m,det_speed,det_energy_J,det_time,runout_x,peak_speed,max_bounce,duration,status,barrier_outcome,max_deflection,peak_force_N,dissipated_J\n";
    for (size_t i = 0; i < n; ++i) {
        const auto& r = b.runs[i];
        rr << i + 1 << ',' << r.initSpeed << ',' << r.rollRes << ',' << r.friction << ',' << r.rest << ',' << r.drag
           << ',' << r.reached << ',' << r.inBand << ',';
        if (r.reached) rr << r.detHeight << ',' << r.detY << ',' << r.detSpeed << ',' << r.detEnergy << ',' << r.detTime;
        else rr << ",,,,";
        rr << ',' << r.runoutX << ',' << r.peakSpeed << ',' << r.maxBounce << ',' << r.duration << ',' << r.status
           << ',' << r.outcome << ',' << r.bs.maxDef << ',' << r.bs.peakForce << ',' << r.bs.diss << '\n';
    }
    auto tj = csv(out / "trajectories.csv");
    tj << "trial,t,x,y,z,vx,vy,vz,speed,barrier_deflection\n";
    for (size_t i = 0; i < n; ++i)
        for (const auto& s : b.runs[i].traj)
            tj << i + 1 << ',' << s.t << ',' << s.x << ',' << s.y << ',' << s.z << ',' << s.vx << ',' << s.vy
               << ',' << s.vz << ',' << s.speed << ',' << s.bd << '\n';
    auto cf = csv(out / "config.csv");
    cf << "parameter,value\nruns," << n << "\nseed," << c.seed << "\ncv," << c.cv
       << "\nterrain_seed," << c.terrainSeed << "\nsegments," << c.terrain.segments
       << "\nroughness," << c.terrain.roughness << "\nvalleys," << c.terrain.valleys
       << "\nvalley_width," << c.terrain.valleyWidth << "\nlength," << c.terrain.length
       << "\ndrop," << c.terrain.drop << "\ninitial_speed," << c.initialSpeed << "\nradius," << c.radius
       << "\nroll_res," << c.rollRes << "\nfriction," << c.friction << "\nrestitution," << c.restitution
       << "\ndrag," << c.drag << "\ndetector_x," << c.detectorX << "\nbin_height," << c.binHeight
       << "\nbarrier," << (int)c.barrier.kind << "\nbarrier_x," << c.barrier.x << "\nbarrier_height," << c.barrier.height
       << "\nbarrier_restitution," << c.barrier.restitution << "\nstiffness," << c.barrier.stiffness
       << "\ndamping," << c.barrier.damping << "\nbarrier_mass," << c.barrier.mass
       << "\nmax_deflection," << c.barrier.maxDeflection << '\n';
    auto bo = csv(out / "barrier_outcomes.csv");
    bo << "outcome,count,total,probability\n";
    std::map<std::string, int> oc;
    for (const auto& r : b.runs) ++oc[r.outcome];
    for (const auto& kv : oc) bo << kv.first << ',' << kv.second << ',' << n << ',' << (double)kv.second / n << '\n';
    auto sm = csv(out / "summary.txt");
    auto ci = wilson(b.inBand, n);
    sm << "runs=" << n << " seed=" << c.seed << "\nreached=" << b.reached << " in_0_10m=" << b.inBand
       << " above=" << b.above << " below=" << b.below << " not_reached=" << b.notReached
       << "\nP(in_0_10m)=" << (n ? (double)b.inBand / n : 0) << " Wilson95%=[" << ci.first << "," << ci.second << "]\n";
    if (c.barrier.kind != BarrierKind::None)
        sm << "intercepted=" << b.intercepted << " overtopped=" << b.overtopped << " failed_passage=" << b.failedPass
           << " unresolved=" << b.hitUnresolved << "\nP(intercepted)=" << (n ? (double)b.intercepted / n : 0) << '\n';
    double my = 0; size_t nr = 0;
    for (const auto& r : b.runs) if (r.reached) { my += r.detY; ++nr; }
    sm << "mean_lateral_y=" << (nr ? my / nr : 0) << " (3D landing spread)\n";
    auto rd = csv(out / "READ_ME.txt");
    rd << "UTF-8 CSV; SI units (m,s,kg,N,J). 3D heightfield model.\n"
       << "Detector: first forward crossing of x=detector_x; height = centre z minus LOCAL terrain height at crossing.\n"
       << "Lateral probability: bin of crossing y across slope width; denominators = reached trials.\n"
       << "Height-band probabilities use ALL trials as denominator. Bins [low,high), last includes 10 m.\n"
       << "Barrier: rigid wall (finite height, overtop possible) / elastic spring-damper-mass panel (fails beyond max deflection).\n"
       << "Truncated normal sampling: sigma=|mean|*CV, out-of-range redrawn. Same seed+config reproduces batch.\n"
       << "Sphere contact in 3D; convex-polygon requirement is fulfilled by the 2D version.\n"
       << "Educational model, not calibrated for engineering design.\n";
}

// ---------------- 研究 ----------------
struct StudyRow { std::string factor, level; double reachP, bandP, speed, energy, latY, peak, runout, dur; size_t n; };
inline StudyRow summarize(const Batch& b, const std::string& f, const std::string& l) {
    StudyRow s; s.factor = f; s.level = l; s.n = b.runs.size();
    double sp = 0, en = 0, ly = 0, pk = 0, ro = 0, du = 0; size_t nr = 0;
    for (const auto& r : b.runs) {
        pk += r.peakSpeed; ro += r.runoutX; du += r.duration;
        if (r.reached) { sp += r.detSpeed; en += r.detEnergy; ly += std::abs(r.detY); ++nr; }
    }
    s.reachP = s.n ? (double)b.reached / s.n : 0;
    s.bandP = s.n ? (double)b.inBand / s.n : 0;
    s.speed = nr ? sp / nr : 0; s.energy = nr ? en / nr : 0; s.latY = nr ? ly / nr : 0;
    s.peak = s.n ? pk / s.n : 0; s.runout = s.n ? ro / s.n : 0; s.dur = s.n ? du / s.n : 0;
    return s;
}
inline void writeStudy(const fs::path& out, const std::vector<StudyRow>& rows, const char* fname) {
    fs::create_directories(out);
    auto f = csv(out / fname);
    f << "factor,level,runs,reach_prob,band_prob,det_speed,det_energy_J,abs_lateral_y,peak_speed,runout_x,duration\n";
    for (const auto& s : rows)
        f << s.factor << ',' << s.level << ',' << s.n << ',' << s.reachP << ',' << s.bandP << ',' << s.speed
          << ',' << s.energy << ',' << s.latY << ',' << s.peak << ',' << s.runout << ',' << s.dur << '\n';
}
inline int runParameterStudy(Config c, int count, uint32_t seed, const std::string& dir) {
    c.finalize();
    fs::path out(dir);
    if (fs::exists(out) && !fs::is_empty(out)) throw std::runtime_error("Output folder must be new or empty");
    Batch base = runBatch(c, count, seed, false);
    std::vector<StudyRow> rows{ summarize(base, "baseline", "1.0") };
    struct Factor { const char* name; double Config::*field; };
    const Factor factors[] = { {"initial_speed",&Config::initialSpeed},{"roll_res",&Config::rollRes},
        {"friction",&Config::friction},{"restitution",&Config::restitution},{"drag",&Config::drag} };
    for (const auto& f : factors)
        for (double m : {0.8, 1.2}) {
            Config cc = c; cc.*f.field *= m;
            Batch b = runBatch(cc, count, seed, false);
            char lv[16]; std::snprintf(lv, 16, "%.1f", m);
            rows.push_back(summarize(b, f.name, lv));
        }
    writeStudy(out, rows, "parameter_study.csv");
    auto f = csv(out / "parameter_findings.txt");
    f << "3D parameter study: 5 factors x 0.8/1.2 mean levels, paired same seed.\n";
    for (const auto& s : rows)
        f << s.factor << " " << s.level << ": reach=" << s.reachP << " band=" << s.bandP
          << " speed=" << s.speed << " |y|=" << s.latY << "\n";
    return (int)rows.size();
}
inline int runTerrainStudy(Config c, int count, uint32_t seed, const std::string& dir) {
    c.finalize();
    fs::path out(dir);
    if (fs::exists(out) && !fs::is_empty(out)) throw std::runtime_error("Output folder must be new or empty");
    Batch base = runBatch(c, count, seed, false);
    std::vector<StudyRow> rows{ summarize(base, "baseline", "1.0") };
    auto scen = [&](const char* name, const char* lv, Config cc) {
        cc.finalize();
        Batch b = runBatch(cc, count, seed, false);
        rows.push_back(summarize(b, name, lv));
    };
    for (double s : {0.8, 1.2}) { Config cc = c; cc.terrain.lengthScale = s; char lv[16]; std::snprintf(lv, 16, "%.1f", s); scen("length", lv, cc); }
    for (double s : {0.8, 1.2}) { Config cc = c; cc.terrain.dropScale = s; char lv[16]; std::snprintf(lv, 16, "%.1f", s); scen("drop", lv, cc); }
    for (double r : {0.2, 0.9}) { Config cc = c; cc.terrain.roughness = r; char lv[16]; std::snprintf(lv, 16, "%.1f", r); scen("roughness", lv, cc); }
    for (int v : {0, 3}) { Config cc = c; cc.terrain.valleys = v; scen("valleys", std::to_string(v).c_str(), cc); }
    writeStudy(out, rows, "terrain_study.csv");
    auto f = csv(out / "terrain_findings.txt");
    f << "3D terrain study: length/drop scale, roughness, valley count; paired same seed; detector follows toe.\n";
    for (const auto& s : rows)
        f << s.factor << " " << s.level << ": reach=" << s.reachP << " band=" << s.bandP
          << " speed=" << s.speed << " |y|=" << s.latY << "\n";
    return (int)rows.size();
}

// ---------------- 渲染 ----------------
struct Camera {
    double yaw = -0.90, pitch = 0.46, dist = 190.0, cx = 60, cy = 0, cz = 22, fov = 880.0;
    bool project(double X, double Y, double Z, double& sx, double& sy, double& depth) const {
        const double dx = X - cx, dy = Y - cy, dz = Z - cz;
        const double cy_ = std::cos(yaw), sy_ = std::sin(yaw), cp = std::cos(pitch), sp = std::sin(pitch);
        const double rx = dx * cy_ + dy * sy_, ry = -dx * sy_ + dy * cy_;
        const double rz = dz * cp - ry * sp, ry2 = dz * sp + ry * cp;
        const double zc = ry2 + dist;
        if (zc < 6.0) return false;
        depth = zc; sx = W * 0.5 + fov * rx / zc; sy = H * 0.52 - fov * rz / zc;
        return true;
    }
};
inline COLORREF shade(double k, int r, int g, int b) {
    k = std::clamp(k, 0.0, 1.0);
    return RGB((int)(r * k), (int)(g * k), (int)(b * k));
}
// ---- 主题色板（与二维版一致的深蓝＋青绿） ----
namespace theme {
    constexpr COLORREF CBG      = RGB(13, 23, 36);
    constexpr COLORREF CPANEL   = RGB(18, 30, 46);
    constexpr COLORREF CCARD    = RGB(24, 38, 56);
    constexpr COLORREF CRULE    = RGB(44, 62, 86);
    constexpr COLORREF CINK     = RGB(236, 243, 250);
    constexpr COLORREF CMUTED   = RGB(148, 170, 194);
    constexpr COLORREF CDIM     = RGB(104, 128, 154);
    constexpr COLORREF CACCENT  = RGB(69, 211, 175);
    constexpr COLORREF CBLUE    = RGB(80, 169, 244);
    constexpr COLORREF CORANGE  = RGB(255, 174, 83);
    constexpr COLORREF CRED     = RGB(240, 98, 88);
}
inline COLORREF mixc(COLORREF a, COLORREF b, double t) {
    t = std::clamp(t, 0.0, 1.0);
    return RGB((int)(GetRValue(a) + (GetRValue(b) - GetRValue(a)) * t),
               (int)(GetGValue(a) + (GetGValue(b) - GetGValue(a)) * t),
               (int)(GetBValue(a) + (GetBValue(b) - GetBValue(a)) * t));
}
inline void skyGradient() {
    for (int y = 0; y < H; y += 4) {
        const double t = (double)y / H;
        const COLORREF c = t < 0.55 ? mixc(RGB(9, 16, 29), RGB(17, 32, 51), t / 0.55)
                                    : mixc(RGB(17, 32, 51), RGB(26, 44, 66), (t - 0.55) / 0.45);
        setfillcolor(c); solidrectangle(0, y, W, y + 4);
    }
}
inline void rrect(int x, int y, int w, int h, int r, COLORREF fill) {
    setfillcolor(fill); setlinecolor(fill);
    fillroundrect(x, y, x + w, y + h, r, r);
}
inline void centerTxt(int x, int y, int w, const wchar_t* s, int size, COLORREF c) {
    settextcolor(c); settextstyle(size, 0, L"微软雅黑"); setbkmode(TRANSPARENT);
    const int tw = textwidth(s);
    outtextxy(x + (w - tw) / 2, y, s);
}
struct TrailPt { double x, y, z; };

inline void renderTerrain(const Terrain& T, const Camera& cam) {
    const int NX = 64, NY = 36;
    struct Cell { double depth; POINT pts[4]; COLORREF col; };
    std::vector<Cell> cells; cells.reserve(NX * NY);
    const double lx = -0.45, ly = -0.55, lz = 0.70, ll = std::sqrt(lx * lx + ly * ly + lz * lz);
    const COLORREF low = RGB(40, 72, 66), mid = RGB(66, 100, 90), high = RGB(116, 136, 126), rock = RGB(102, 116, 126);
    const COLORREF fog = RGB(19, 33, 50);
    for (int i = 0; i < NX; ++i) for (int j = 0; j < NY; ++j) {
        const double x0 = T.length * i / NX, x1 = T.length * (i + 1) / NX;
        const double y0 = -T.width * 0.5 + T.width * j / NY, y1 = -T.width * 0.5 + T.width * (j + 1) / NY;
        const double xs[4] = { x0, x1, x1, x0 }, ys[4] = { y0, y0, y1, y1 };
        Cell c{}; double ds = 0; bool ok = true;
        for (int k = 0; k < 4; ++k) {
            double sx, sy, dp;
            if (!cam.project(xs[k], ys[k], T.height(xs[k], ys[k]), sx, sy, dp)) { ok = false; break; }
            c.pts[k] = { (LONG)sx, (LONG)sy }; ds += dp;
        }
        if (!ok) continue;
        c.depth = ds / 4;
        const double xm = (x0 + x1) * 0.5, ym = (y0 + y1) * 0.5;
        double gx, gy; T.grad(xm, ym, gx, gy);
        double nx = -gx, ny = -gy, nz = 1.0;
        const double nl = std::sqrt(nx * nx + ny * ny + nz * nz);
        const double lam = std::max(0.0, (nx * lx + ny * ly + nz * lz) / (nl * ll));
        const double hgt = std::clamp(T.height(xm, ym) / std::max(T.drop, 1.0), 0.0, 1.0);
        const double steep = std::clamp((std::sqrt(gx * gx + gy * gy) - 0.55) / 0.9, 0.0, 1.0);
        COLORREF base = hgt < 0.5 ? mixc(low, mid, hgt * 2) : mixc(mid, high, (hgt - 0.5) * 2);
        base = mixc(base, rock, steep * 0.45);
        base = shade(0.58 + 0.50 * lam, GetRValue(base), GetGValue(base), GetBValue(base));
        const double fogT = std::clamp((c.depth - 150.0) / 260.0, 0.0, 0.66);
        c.col = mixc(base, fog, fogT);
        cells.push_back(c);
    }
    std::sort(cells.begin(), cells.end(), [](const Cell& a, const Cell& b) { return a.depth > b.depth; });
    setlinestyle(PS_SOLID, 1);
    for (const auto& c : cells) { setfillcolor(c.col); setlinecolor(c.col); fillpolygon(c.pts, 4); }
}
inline void renderBarrier(const Terrain& T, const Barrier& B, const BarrierState& bs, const Camera& cam) {
    if (B.kind == BarrierKind::None) return;
    const bool el = B.kind == BarrierKind::Elastic;
    const double wx = B.x + (el && !bs.failed ? bs.d : 0.0);
    const int NS = 32;
    const COLORREF col = bs.failed ? theme::CRED : (el ? theme::CBLUE : theme::CORANGE);
    setlinestyle(PS_SOLID, 1);
    double topAx = 0, topAy = 0, topBx = 0, topBy = 0; bool haveTop = false;
    for (int i = 0; i < NS; ++i) {
        const double y0 = -T.width * 0.5 + T.width * i / NS, y1 = -T.width * 0.5 + T.width * (i + 1) / NS;
        const double g0 = T.height(B.x, y0), g1 = T.height(B.x, y1);
        double sx[4], sy[4], dp;
        bool ok = true;
        const double X[4] = { wx, wx, wx, wx }, Y[4] = { y0, y1, y1, y0 }, Z[4] = { g0, g1, g1 + B.height, g0 + B.height };
        POINT pts[4];
        for (int k = 0; k < 4; ++k) { if (!cam.project(X[k], Y[k], Z[k], sx[k], sy[k], dp)) { ok = false; break; } pts[k] = { (LONG)sx[k], (LONG)sy[k] }; }
        if (!ok) continue;
        if (bs.failed && i % 2 == 0) continue;   // 失效：画成破口
        const COLORREF seg = shade(0.82 + 0.18 * std::sin(i * 0.7), GetRValue(col), GetGValue(col), GetBValue(col));
        setfillcolor(seg); setlinecolor(seg);
        fillpolygon(pts, 4);
        if (!haveTop) { topAx = sx[3]; topAy = sy[3]; haveTop = true; }
        topBx = sx[2]; topBy = sy[2];
    }
    if (haveTop) {
        setlinecolor(mixc(col, RGB(255, 255, 255), 0.55)); setlinestyle(PS_SOLID, 2);
        line((int)topAx, (int)topAy, (int)topBx, (int)topBy);
    }
}
inline void renderTrail(const std::vector<TrailPt>& tr, const Camera& cam, COLORREF base) {
    if (tr.size() < 2) return;
    setlinestyle(PS_SOLID, 2);
    for (size_t i = 1; i < tr.size(); ++i) {
        double ax, ay, ad, bx, by, bd;
        if (!cam.project(tr[i - 1].x, tr[i - 1].y, tr[i - 1].z, ax, ay, ad)) continue;
        if (!cam.project(tr[i].x, tr[i].y, tr[i].z, bx, by, bd)) continue;
        const double f = (double)i / tr.size();
        setlinecolor(shade(0.5 + 0.5 * f, GetRValue(base), GetGValue(base), GetBValue(base)));
        line((int)ax, (int)ay, (int)bx, (int)by);
    }
}
inline void renderBall(const Ball& b, const Terrain& T, const Camera& cam) {
    double sx, sy, sd;
    if (cam.project(b.x, b.y, T.height(b.x, b.y) + 0.05, sx, sy, sd)) {
        const double r = cam.fov * b.radius / sd;
        setfillcolor(RGB(10, 18, 28));
        solidellipse((int)(sx - r * 1.15), (int)(sy - r * 0.40), (int)(sx + r * 1.15), (int)(sy + r * 0.40));
        setfillcolor(RGB(6, 12, 20));
        solidellipse((int)(sx - r * 0.8), (int)(sy - r * 0.28), (int)(sx + r * 0.8), (int)(sy + r * 0.28));
    }
    double bx, by, bd;
    if (!cam.project(b.x, b.y, b.z, bx, by, bd)) return;
    const double r = std::max(3.0, cam.fov * b.radius / bd);
    for (int i = 16; i >= 1; --i) {
        const double f = (double)i / 16, rr = r * f;
        const double ox = -r * 0.32 * (1 - f), oy = -r * 0.36 * (1 - f);
        setfillcolor(shade(0.34 + 0.72 * (1 - f), 255, 176, 96));
        solidellipse((int)(bx + ox - rr), (int)(by + oy - rr), (int)(bx + ox + rr), (int)(by + oy + rr));
    }
    setfillcolor(mixc(RGB(255, 220, 170), RGB(255, 255, 255), 0.35));
    solidellipse((int)(bx - r * 0.52), (int)(by - r * 0.60), (int)(bx - r * 0.16), (int)(by - r * 0.26));
    setlinecolor(RGB(146, 84, 40)); setlinestyle(PS_SOLID, 2);
    line((int)(bx + std::cos(b.spin) * r * 0.72), (int)(by + std::sin(b.spin) * r * 0.72),
         (int)(bx - std::cos(b.spin) * r * 0.72), (int)(by - std::sin(b.spin) * r * 0.72));
}
inline void txt(int x, int y, const wchar_t* s, int size, COLORREF c) {
    settextcolor(c); settextstyle(size, 0, L"微软雅黑"); setbkmode(TRANSPARENT); outtextxy(x, y, s);
}

// ---------------- 参数面板 ----------------
constexpr int PW = 300;                 // 左侧面板宽
struct PanelRow { const wchar_t* label; wchar_t value[48]; };
inline void fmtVal(wchar_t* dst, double v, int prec = 2) { swprintf(dst, 48, L"%.*f", prec, v); }

inline int buildRows(Config& c, int page, PanelRow* rows) {
    int n = 0;
    auto add = [&](const wchar_t* l, const wchar_t* v) { rows[n].label = l; wcscpy_s(rows[n].value, v); ++n; };
    wchar_t tmp[48];
    if (page == 0) {
        fmtVal(tmp, c.radius); add(L"球半径 m", tmp);
        fmtVal(tmp, c.density, 0); add(L"密度 kg/m³", tmp);
        fmtVal(tmp, 4.0 / 3.0 * PI * c.radius * c.radius * c.radius); add(L"体积 m³", tmp);
        fmtVal(tmp, c.density * 4.0 / 3.0 * PI * c.radius * c.radius * c.radius, 0); add(L"质量 kg", tmp);
        fmtVal(tmp, c.initialSpeed); add(L"初速 m/s", tmp);
        fmtVal(tmp, c.rollRes, 3); add(L"滚动阻力", tmp);
        fmtVal(tmp, c.friction); add(L"摩擦系数", tmp);
        fmtVal(tmp, c.restitution); add(L"恢复系数", tmp);
        fmtVal(tmp, c.drag); add(L"空气阻力 Cd", tmp);
        fmtVal(tmp, c.cv); add(L"变异系数 CV", tmp);
        swprintf(tmp, 48, L"%d", c.runs); add(L"批量次数 N", tmp);
        swprintf(tmp, 48, L"%u", c.seed); add(L"随机种子", tmp);
    } else if (page == 1) {
        add(L"挡板类型", c.barrier.kind == BarrierKind::None ? L"无" : (c.barrier.kind == BarrierKind::Rigid ? L"刚性" : L"弹性"));
        fmtVal(tmp, c.barrier.x, 0); add(L"位置 x m", tmp);
        fmtVal(tmp, c.barrier.height, 1); add(L"高度 m", tmp);
        fmtVal(tmp, c.barrier.restitution); add(L"恢复系数", tmp);
        fmtVal(tmp, c.barrier.stiffness / 1000, 0); add(L"刚度 kN/m", tmp);
        fmtVal(tmp, c.barrier.damping / 1000, 0); add(L"阻尼 kN·s/m", tmp);
        fmtVal(tmp, c.barrier.mass, 0); add(L"等效质量 kg", tmp);
        fmtVal(tmp, c.barrier.maxDeflection, 1); add(L"失效位移 m", tmp);
    } else {
        swprintf(tmp, 48, L"%u", c.terrainSeed); add(L"地形种子", tmp);
        swprintf(tmp, 48, L"%d", c.terrain.segments); add(L"坡段数", tmp);
        fmtVal(tmp, c.terrain.roughness); add(L"起伏强度", tmp);
        swprintf(tmp, 48, L"%d", c.terrain.valleys); add(L"沟谷数", tmp);
        fmtVal(tmp, c.terrain.valleyWidth, 2); add(L"谷宽倍率", tmp);
        fmtVal(tmp, c.terrain.lengthScale, 2); add(L"长度倍率", tmp);
        fmtVal(tmp, c.terrain.dropScale, 2); add(L"落差倍率", tmp);
        fmtVal(tmp, c.terrain.drop, 1); add(L"当前落差 m", tmp);
        fmtVal(tmp, c.detectorX, 0); add(L"检测断面 x", tmp);
        fmtVal(tmp, c.binHeight, 2); add(L"高度分箱 m", tmp);
    }
    return n;
}

// 返回点击命中的行号与方向（-1 左 / +1 右），未命中返回 -1
inline int panelHit(int x, int y, int nrows, int& dir) {
    if (x < 8 || x > PW - 8) return -1;
    // 页签
    if (y >= 40 && y <= 66) { dir = 0; return 100 + (x - 8) / ((PW - 16) / 3); }
    // 2D/3D 切换
    if (y >= 72 && y <= 98) { dir = 0; return 200; }
    const int y0 = 108, rh = 27;
    if (y < y0) return -1;
    const int row = (y - y0) / rh;
    if (row < 0 || row >= nrows) return -1;
    if (x >= PW - 62 && x <= PW - 36) { dir = -1; return row; }
    if (x >= PW - 32 && x <= PW - 8) { dir = +1; return row; }
    return -1;
}

inline void renderPanel(Config& c, int page, int nrows, const PanelRow* rows, bool mode3d, bool busy) {
    setfillcolor(theme::CPANEL); solidrectangle(0, 0, PW, H);
    setfillcolor(theme::CACCENT); solidrectangle(0, 0, 3, H);
    setlinestyle(PS_SOLID, 1); setlinecolor(theme::CRULE); line(PW, 0, PW, H);
    txt(14, 10, L"参数面板", 17, theme::CINK);
    settextstyle(11, 0, L"微软雅黑"); settextcolor(theme::CDIM); setbkmode(TRANSPARENT);
    outtextxy(PW - 78, 16, L"ROLLINGSTONE");
    const wchar_t* tabs[3] = { L"球", L"挡板", L"地形" };
    const int tw = (PW - 16) / 3;
    for (int i = 0; i < 3; ++i) {
        const int x = 8 + i * tw;
        rrect(x, 40, tw - 4, 26, 10, i == page ? theme::CACCENT : theme::CCARD);
        centerTxt(x, 45, tw - 4, tabs[i], 14, i == page ? RGB(10, 26, 22) : theme::CMUTED);
    }
    rrect(8, 72, PW / 2 - 10, 26, 10, mode3d ? theme::CBLUE : theme::CCARD);
    centerTxt(8, 77, PW / 2 - 10, L"3D 视图", 13, mode3d ? RGB(8, 20, 34) : theme::CMUTED);
    rrect(PW / 2 + 2, 72, PW / 2 - 10, 26, 10, !mode3d ? theme::CBLUE : theme::CCARD);
    centerTxt(PW / 2 + 2, 77, PW / 2 - 10, L"2D 剖面", 13, !mode3d ? RGB(8, 20, 34) : theme::CMUTED);
    const int y0 = 108, rh = 27;
    for (int i = 0; i < nrows; ++i) {
        const int y = y0 + i * rh;
        if (i % 2 == 0) { setfillcolor(RGB(21, 34, 51)); solidrectangle(6, y, PW - 6, y + rh); }
        txt(12, y + 5, rows[i].label, 12, theme::CMUTED);
        settextstyle(13, 0, L"微软雅黑"); settextcolor(theme::CACCENT); setbkmode(TRANSPARENT);
        const int vw = textwidth(rows[i].value);
        outtextxy(PW - 70 - vw, y + 5, rows[i].value);
        rrect(PW - 62, y + 3, 26, 20, 8, theme::CCARD);
        rrect(PW - 32, y + 3, 24, 20, 8, theme::CCARD);
        centerTxt(PW - 62, y + 4, 26, L"−", 13, theme::CINK);
        centerTxt(PW - 32, y + 4, 24, L"+", 13, theme::CINK);
    }
    const int by = y0 + nrows * rh + 10;
    const wchar_t* acts[4] = { L"B 批量", L"E 导出", L"A 参研", L"G 地研" };
    const int aw = (PW - 16) / 4;
    for (int i = 0; i < 4; ++i) {
        const int x = 8 + i * aw;
        rrect(x, by, aw - 4, 26, 10, busy ? RGB(35, 42, 52) : (i == 0 ? theme::CACCENT : theme::CCARD));
        centerTxt(x, by + 6, aw - 4, acts[i], 11, busy ? theme::CDIM : (i == 0 ? RGB(10, 26, 22) : theme::CINK));
    }
    setlinestyle(PS_SOLID, 1); setlinecolor(theme::CRULE); line(10, by + 36, PW - 10, by + 36);
    txt(12, by + 44, L"Space 暂停 · R 重置 · T 尾迹", 11, theme::CDIM);
    txt(12, by + 62, L"D 切换2D/3D · 拖动旋转 · 滚轮缩放", 11, theme::CDIM);
    txt(12, by + 80, L", . 倍速 · Esc 退出", 11, theme::CDIM);
}

// 面板调参：page/row/dir -> 修改 cfg；返回是否需要重生地形/重置球
inline bool applyParam(Config& c, int page, int row, int dir, bool& needRegen, bool& needReball) {
    needRegen = false; needReball = false;
    auto cl = [](double v, double lo, double hi) { return std::clamp(v, lo, hi); };
    if (page == 0) {
        switch (row) {
        case 0: c.radius = cl(c.radius + dir * 0.1, 0.3, 4.0); needReball = true; break;
        case 1: c.density = cl(c.density + dir * 100, 500, 8000); needReball = true; break;
        case 2: case 3: return false;   // 体积/质量为派生量
        case 4: c.initialSpeed = cl(c.initialSpeed + dir * 0.5, 0, 30); needReball = true; break;
        case 5: c.rollRes = cl(c.rollRes + dir * 0.005, 0, 0.3); needReball = true; break;
        case 6: c.friction = cl(c.friction + dir * 0.05, 0.05, 1.5); needReball = true; break;
        case 7: c.restitution = cl(c.restitution + dir * 0.05, 0, 0.95); needReball = true; break;
        case 8: c.drag = cl(c.drag + dir * 0.1, 0, 3); needReball = true; break;
        case 9: c.cv = cl(c.cv + dir * 0.02, 0, 0.5); break;
        case 10: c.runs = (int)cl(c.runs + dir * 50, 10, 2000); break;
        case 11: c.seed += (uint32_t)(dir > 0 ? 1 : -1); break;
        }
    } else if (page == 1) {
        switch (row) {
        case 0: c.barrier.kind = (BarrierKind)(((int)c.barrier.kind + dir + 3) % 3); needReball = true; break;
        case 1: c.barrier.x = cl(c.barrier.x + dir * 4, 20, c.terrain.length - 4); needReball = true; break;
        case 2: c.barrier.height = cl(c.barrier.height + dir * 0.5, 1, 20); needReball = true; break;
        case 3: c.barrier.restitution = cl(c.barrier.restitution + dir * 0.05, 0, 0.95); needReball = true; break;
        case 4: c.barrier.stiffness = cl(c.barrier.stiffness + dir * 2.5e5, 2.5e5, 1e7); needReball = true; break;
        case 5: c.barrier.damping = cl(c.barrier.damping + dir * 2.5e4, 0, 1e6); needReball = true; break;
        case 6: c.barrier.mass = cl(c.barrier.mass + dir * 250, 250, 10000); needReball = true; break;
        case 7: c.barrier.maxDeflection = cl(c.barrier.maxDeflection + dir * 0.5, 0.5, 10); needReball = true; break;
        }
    } else {
        switch (row) {
        case 0: c.terrainSeed += (uint32_t)(dir > 0 ? 1 : -1); needRegen = true; break;
        case 1: c.terrain.segments = (int)cl(c.terrain.segments + dir, 3, 24); needRegen = true; break;
        case 2: c.terrain.roughness = cl(c.terrain.roughness + dir * 0.1, 0, 1); needRegen = true; break;
        case 3: c.terrain.valleys = (int)cl(c.terrain.valleys + dir, 0, 3); needRegen = true; break;
        case 4: c.terrain.valleyWidth = cl(c.terrain.valleyWidth + dir * 0.25, 0.4, 2); needRegen = true; break;
        case 5: c.terrain.lengthScale = cl(c.terrain.lengthScale + dir * 0.1, 0.5, 2); needRegen = true; break;
        case 6: c.terrain.dropScale = cl(c.terrain.dropScale + dir * 0.1, 0.5, 2); needRegen = true; break;
        case 7: return false;
        case 8: c.detectorX = cl(c.detectorX + dir * 5, 10, c.terrain.length); needReball = true; break;
        case 9: c.binHeight = cl(c.binHeight + dir * 0.25, 0.25, 2); break;
        }
    }
    return true;
}

// ---------------- 2D 剖面渲染（中心线剖面，物理仍为三维） ----------------
struct Cam2 {
    double scale = 7.0, ox = 0, oz = 0;   // 像素/米，世界原点偏移
    bool project(double X, double Z, double& sx, double& sy) const {
        sx = PW + 30 + (X - ox) * scale;
        sy = H - 170 - (Z - oz) * scale;
        return sx > PW - 40 && sx < W + 60 && sy > -60 && sy < H + 60;
    }
};
inline void render2D(const Terrain& T, const Barrier& B, const BarrierState& bs, const Ball& b,
                     const std::vector<TrailPt>& trail, const std::vector<std::vector<TrailPt>>& batchTrails,
                     const Cam2& cam, double detectorX) {
    // 坡面剖面填充
    const int NS = 240;
    std::vector<POINT> poly;
    double sx, sy;
    for (int i = 0; i <= NS; ++i) {
        const double x = T.length * i / NS;
        if (cam.project(x, T.height(x, 0), sx, sy)) poly.push_back({ (LONG)sx, (LONG)sy });
    }
    if (poly.size() > 2) {
        poly.push_back({ (LONG)(PW + 30 + (T.length - cam.ox) * cam.scale), H });
        poly.push_back({ (LONG)(PW + 30 - cam.ox * cam.scale), H });
        setfillcolor(RGB(38, 56, 62)); setlinecolor(RGB(38, 56, 62));
        fillpolygon(poly.data(), (int)poly.size());
    }
    setlinestyle(PS_SOLID, 2); setlinecolor(theme::CACCENT);
    for (int i = 1; i <= NS; ++i) {
        double ax, ay, bx2, by2;
        const double x0 = T.length * (i - 1) / NS, x1 = T.length * i / NS;
        if (cam.project(x0, T.height(x0, 0), ax, ay) && cam.project(x1, T.height(x1, 0), bx2, by2))
            line((int)ax, (int)ay, (int)bx2, (int)by2);
    }
    // 检测断面 0–10 m
    double dx0, dy0, dx1, dy1;
    if (cam.project(detectorX, T.height(detectorX, 0), dx0, dy0) && cam.project(detectorX, T.height(detectorX, 0) + 10, dx1, dy1)) {
        setlinecolor(theme::CBLUE); setlinestyle(PS_SOLID, 2);
        line((int)dx0, (int)dy0, (int)dx1, (int)dy1);
        txt((int)dx0 - 40, (int)dy1 - 20, L"0–10 m", 12, theme::CBLUE);
    }
    // 挡板剖面
    if (B.kind != BarrierKind::None) {
        const double wx = B.x + (B.kind == BarrierKind::Elastic && !bs.failed ? bs.d : 0);
        double bx0, by0, bx1, by1;
        if (cam.project(wx, T.height(B.x, 0), bx0, by0) && cam.project(wx, T.height(B.x, 0) + B.height, bx1, by1)) {
            const COLORREF col = bs.failed ? theme::CRED : (B.kind == BarrierKind::Elastic ? theme::CBLUE : theme::CORANGE);
            setlinecolor(col); setlinestyle(PS_SOLID, bs.failed ? 2 : 5);
            line((int)bx0, (int)by0, (int)bx1, (int)by1);
        }
    }
    // 批量轨迹（剖面投影）
    setlinestyle(PS_SOLID, 1);
    for (const auto& t : batchTrails) {
        setlinecolor(RGB(72, 116, 150));
        for (size_t i = 1; i < t.size(); ++i) {
            double ax, ay, bx2, by2;
            if (cam.project(t[i - 1].x, t[i - 1].z, ax, ay) && cam.project(t[i].x, t[i].z, bx2, by2))
                line((int)ax, (int)ay, (int)bx2, (int)by2);
        }
    }
    // 当前轨迹
    if (trail.size() > 1) {
        setlinestyle(PS_SOLID, 2); setlinecolor(theme::CACCENT);
        for (size_t i = 1; i < trail.size(); ++i) {
            double ax, ay, bx2, by2;
            if (cam.project(trail[i - 1].x, trail[i - 1].z, ax, ay) && cam.project(trail[i].x, trail[i].z, bx2, by2))
                line((int)ax, (int)ay, (int)bx2, (int)by2);
        }
    }
    // 球
    double cx2, cy2;
    if (cam.project(b.x, b.z, cx2, cy2)) {
        const int r = (int)std::max(3.0, b.radius * cam.scale);
        setfillcolor(theme::CORANGE); setlinecolor(RGB(146, 84, 40));
        solidcircle((int)cx2, (int)cy2, r);
        circle((int)cx2, (int)cy2, r);
        setfillcolor(RGB(255, 220, 170));
        solidcircle((int)(cx2 - r * 0.3), (int)(cy2 - r * 0.35), (int)std::max(1.0, r * 0.22));
        line((int)(cx2 + std::cos(b.spin) * r * 0.7), (int)(cy2 + std::sin(b.spin) * r * 0.7),
             (int)(cx2 - std::cos(b.spin) * r * 0.7), (int)(cy2 - std::sin(b.spin) * r * 0.7));
    }
}

// ---------------- 交互界面 ----------------
inline int runViewer(const Config& startCfg, const wchar_t* shotPath) {
    Config cfg = startCfg; cfg.finalize();
    Ball b; b.place(cfg.terrain, cfg, cfg.terrain.startY());
    Camera cam;
    std::vector<TrailPt> trail;
    std::vector<std::vector<TrailPt>> batchTrails;
    bool paused = false, trailOn = true, hasBatch = false, mode3d = true;
    double speedScale = 1.0;
    int page = 0;
    Cam2 cam2;
    Batch batch;
    std::future<Batch> fut;
    std::wstring note = L"左侧面板调参（球/挡板/地形）；D 切换 2D/3D；B 批量 E 导出 A/G 研究。";
    const bool shot = shotPath != nullptr;

    HWND window = initgraph(W, H);
    SetWindowTextW(window, L"山区滚石 · 三维完整版（统计/挡板/研究）");
    setbkcolor(RGB(9, 17, 30));
    BeginBatchDraw();
    auto last = std::chrono::steady_clock::now();
    double fps = 60; int frames = 0;
    const double dt = 0.002;

    auto reball = [&]() { b.place(cfg.terrain, cfg, cfg.terrain.startY()); trail.clear(); paused = false; };
    auto regen = [&]() { cfg.terrain.generate(cfg.terrainSeed); cfg.detectorX = cfg.terrain.length; reball(); };

    while (IsWindow(window)) {
        auto now = std::chrono::steady_clock::now();
        const double el = std::min(0.05, std::chrono::duration<double>(now - last).count());
        last = now;
        fps = fps * 0.9 + (el > 1e-6 ? 1.0 / el : 60.0) * 0.1;

        if (fut.valid() && fut.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
            batch = fut.get(); hasBatch = true;
            batchTrails.clear();
            for (size_t i = 0; i < batch.runs.size() && i < 30; ++i) {
                std::vector<TrailPt> t;
                for (size_t k = 0; k < batch.runs[i].traj.size(); k += 4)
                    t.push_back({ batch.runs[i].traj[k].x, batch.runs[i].traj[k].y, batch.runs[i].traj[k].z });
                batchTrails.push_back(std::move(t));
            }
            wchar_t buf[256];
            swprintf(buf, 256, L"批量完成：%zu 次，到达 %zu，0–10m %zu（%.1f%%）", batch.runs.size(), batch.reached, batch.inBand,
                     batch.runs.size() ? 100.0 * batch.inBand / batch.runs.size() : 0);
            note = buf;
        }

        if (!shot) {
            ExMessage m;
            static bool drag = false; static int lx = 0, ly = 0;
            PanelRow rows[16];
            const int nrows = buildRows(cfg, page, rows);
            while (peekmessage(&m, EX_MOUSE | EX_KEY)) {
                if (m.message == WM_LBUTTONDOWN) {
                    if (m.x < PW) {
                        int dir = 0;
                        const int hit = panelHit(m.x, m.y, nrows, dir);
                        if (hit >= 200) { mode3d = !mode3d; }
                        else if (hit >= 100) { page = std::clamp(hit - 100, 0, 2); }
                        else if (hit >= 0 && !fut.valid()) {
                            bool rg = false, rb = false;
                            if (applyParam(cfg, page, hit, dir, rg, rb)) {
                                if (rg) { cfg.terrain.generate(cfg.terrainSeed); cfg.detectorX = cfg.terrain.length; }
                                if (rg || rb) reball();
                            }
                        }
                        // 底部动作按钮
                        const int by = 108 + nrows * 27 + 10;
                        if (m.y >= by && m.y <= by + 26 && !fut.valid()) {
                            const int ai = (m.x - 8) / ((PW - 16) / 4);
                            if (ai == 0) { Config cc = cfg; int n = cc.runs; uint32_t sd = cc.seed;
                                fut = std::async(std::launch::async, [cc, n, sd]() { return runBatch(cc, n, sd, true); }); note = L"批量计算中…"; }
                            else if (ai == 1 && hasBatch) { try { exportBatch(batch, "results3d/batch_" + std::to_string(std::chrono::system_clock::now().time_since_epoch().count() / 1000000000)); note = L"已导出到 results3d/"; } catch (...) { note = L"导出失败（目录须新建或为空）"; } }
                            else if (ai == 2) { Config cc = cfg; int n = std::max(30, cc.runs / 4); uint32_t sd = cc.seed;
                                fut = std::async(std::launch::async, [cc, n, sd]() { runParameterStudy(cc, n, sd, "results3d/param_" + std::to_string(std::chrono::system_clock::now().time_since_epoch().count() / 1000000000)); return runBatch(cc, 1, sd, true); }); note = L"参数研究计算中…"; }
                            else if (ai == 3) { Config cc = cfg; int n = std::max(30, cc.runs / 4); uint32_t sd = cc.seed;
                                fut = std::async(std::launch::async, [cc, n, sd]() { runTerrainStudy(cc, n, sd, "results3d/terrain_" + std::to_string(std::chrono::system_clock::now().time_since_epoch().count() / 1000000000)); return runBatch(cc, 1, sd, true); }); note = L"地形研究计算中…"; }
                        }
                    } else { drag = true; lx = m.x; ly = m.y; }
                }
                else if (m.message == WM_LBUTTONUP) drag = false;
                else if (m.message == WM_MOUSEMOVE && drag) {
                    if (mode3d) {
                        cam.yaw += (m.x - lx) * 0.006;
                        cam.pitch = std::clamp(cam.pitch + (m.y - ly) * 0.004, 0.08, 1.35);
                    } else {
                        cam2.ox -= (m.x - lx) / cam2.scale;
                        cam2.oz += (m.y - ly) / cam2.scale;
                    }
                    lx = m.x; ly = m.y;
                }
                else if (m.message == WM_MOUSEWHEEL) {
                    if (mode3d) cam.dist = std::clamp(cam.dist - m.wheel * 0.06, 60.0, 420.0);
                    else cam2.scale = std::clamp(cam2.scale * (m.wheel > 0 ? 1.1 : 0.9), 2.0, 30.0);
                }
                else if (m.message == WM_KEYDOWN) {
                    switch (m.vkcode) {
                    case VK_ESCAPE: EndBatchDraw(); closegraph(); return 0;
                    case VK_SPACE: paused = !paused; break;
                    case 'D': mode3d = !mode3d; break;
                    case 'R': reball(); break;
                    case 'T': trailOn = !trailOn; if (!trailOn) trail.clear(); break;
                    case 'N': cfg.terrainSeed += 1; regen(); break;
                    case VK_OEM_4: cfg.terrain.segments = std::max(3, cfg.terrain.segments - 1); regen(); break;
                    case VK_OEM_6: cfg.terrain.segments = std::min(24, cfg.terrain.segments + 1); regen(); break;
                    case VK_OEM_MINUS: cfg.terrain.roughness = std::max(0.0, cfg.terrain.roughness - 0.1); regen(); break;
                    case VK_OEM_PLUS: cfg.terrain.roughness = std::min(1.0, cfg.terrain.roughness + 0.1); regen(); break;
                    case 'V': cfg.terrain.valleys = (cfg.terrain.valleys + 1) % 4; regen(); break;
                    case 'W': cfg.terrain.valleyWidth = cfg.terrain.valleyWidth >= 2.0 ? 0.5 : cfg.terrain.valleyWidth + 0.25; regen(); break;
                    case '1': cfg.barrier.kind = BarrierKind::None; reball(); break;
                    case '2': cfg.barrier.kind = BarrierKind::Rigid; reball(); break;
                    case '3': cfg.barrier.kind = BarrierKind::Elastic; reball(); break;
                    case '4': cfg.barrier.x = std::max(20.0, cfg.barrier.x - 4); reball(); break;
                    case '5': cfg.barrier.x = std::min(cfg.terrain.length - 4, cfg.barrier.x + 4); reball(); break;
                    case '6': cfg.barrier.height = std::max(1.0, cfg.barrier.height - 0.5); reball(); break;
                    case '7': cfg.barrier.height = std::min(20.0, cfg.barrier.height + 0.5); reball(); break;
                    case VK_OEM_COMMA: speedScale = std::max(0.25, speedScale / 2); break;
                    case VK_OEM_PERIOD: speedScale = std::min(4.0, speedScale * 2); break;
                    case 'B': if (!fut.valid()) { Config cc = cfg; int n = cc.runs; uint32_t sd = cc.seed;
                        fut = std::async(std::launch::async, [cc, n, sd]() { return runBatch(cc, n, sd, true); });
                        note = L"批量计算中…"; } break;
                    case 'E': if (hasBatch) { try { exportBatch(batch, "results3d/batch_" + std::to_string(std::chrono::system_clock::now().time_since_epoch().count() / 1000000000)); note = L"已导出到 results3d/"; } catch (const std::exception& e) { note = L"导出失败（目录须新建或为空）"; } } break;
                    case 'A': if (!fut.valid()) { Config cc = cfg; int n = std::max(30, cc.runs / 4); uint32_t sd = cc.seed;
                        fut = std::async(std::launch::async, [cc, n, sd]() { runParameterStudy(cc, n, sd, "results3d/param_" + std::to_string(std::chrono::system_clock::now().time_since_epoch().count() / 1000000000)); return runBatch(cc, 1, sd, true); });
                        note = L"参数研究计算中…（结果在 results3d/）"; } break;
                    case 'G': if (!fut.valid()) { Config cc = cfg; int n = std::max(30, cc.runs / 4); uint32_t sd = cc.seed;
                        fut = std::async(std::launch::async, [cc, n, sd]() { runTerrainStudy(cc, n, sd, "results3d/terrain_" + std::to_string(std::chrono::system_clock::now().time_since_epoch().count() / 1000000000)); return runBatch(cc, 1, sd, true); });
                        note = L"地形研究计算中…（结果在 results3d/）"; } break;
                    }
                }
            }
        }

        if (!paused) {
            const int sub = (cfg.barrier.kind == BarrierKind::Elastic) ? 4 : 1;
            const double h = dt / sub;
            int steps = std::clamp((int)std::lround(el * speedScale / h), 1, 1600);
            for (int i = 0; i < steps; ++i) {
                b.step(cfg.terrain, cfg.barrier, h);
                if (trailOn && i % (6 * sub) == 0) { trail.push_back({ b.x, b.y, b.z }); if (trail.size() > 1400) trail.erase(trail.begin(), trail.begin() + 200); }
            }
        }
        const double tx = 0.45 * b.x + 0.55 * (cfg.terrain.length * 0.5);
        cam.cx += (tx - cam.cx) * 0.03;
        cam.cy += (b.y * 0.5 - cam.cy) * 0.03;
        cam.cz += (cfg.terrain.drop * 0.34 - cam.cz) * 0.03;
        // 2D 相机自动取景
        cam2.scale = std::clamp((W - PW - 80) / cfg.terrain.length, 2.0, 30.0);
        cam2.ox = -4; cam2.oz = -2;

        skyGradient();
        if (mode3d) {
            renderTerrain(cfg.terrain, cam);
            if (hasBatch) for (const auto& t : batchTrails) renderTrail(t, cam, RGB(96, 148, 182));
            renderBarrier(cfg.terrain, cfg.barrier, b.bs, cam);
            if (trailOn) renderTrail(trail, cam, theme::CACCENT);
            renderBall(b, cfg.terrain, cam);
        } else {
            render2D(cfg.terrain, cfg.barrier, b.bs, b, trail, batchTrails, cam2, cfg.detectorX);
        }
        PanelRow rows[16];
        const int nrows = buildRows(cfg, page, rows);
        renderPanel(cfg, page, nrows, rows, mode3d, fut.valid());

        wchar_t buf[512];
        // 标题栏
        setfillcolor(RGB(12, 21, 34)); solidrectangle(PW, 0, W, 52);
        setfillcolor(theme::CACCENT); solidrectangle(PW, 50, W, 52);
        txt(PW + 18, 12, mode3d ? L"山区滚石 · 三维模拟" : L"山区滚石 · 二维剖面", 21, theme::CINK);
        swprintf(buf, 512, L"seed %u · %s · %.2fx", cfg.terrainSeed, paused ? L"暂停" : L"运行", speedScale);
        settextstyle(13, 0, L"微软雅黑"); settextcolor(theme::CMUTED); setbkmode(TRANSPARENT);
        outtextxy(W - 300, 18, buf);
        // 挡板状态徽标
        const wchar_t* bk = cfg.barrier.kind == BarrierKind::None ? L"无挡板" : (cfg.barrier.kind == BarrierKind::Rigid ? L"刚性挡板" : L"弹性挡板");
        const COLORREF bc = cfg.barrier.kind == BarrierKind::None ? theme::CDIM : (cfg.barrier.kind == BarrierKind::Rigid ? theme::CORANGE : theme::CBLUE);
        rrect(W - 130, 12, 112, 26, 12, theme::CCARD);
        centerTxt(W - 130, 17, 112, bk, 12, bc);
        // 底部统计卡片条
        setfillcolor(RGB(12, 21, 34)); solidrectangle(PW, H - 132, W, H);
        setlinestyle(PS_SOLID, 1); setlinecolor(theme::CRULE); line(PW, H - 132, W, H - 132);
        const int cardY = H - 122, cardH = 62;
        const int areaW = W - PW - 24, cw = (areaW - 3 * 10) / 4;
        auto card = [&](int i, const wchar_t* label, const wchar_t* value, COLORREF vc) {
            const int x = PW + 12 + i * (cw + 10);
            rrect(x, cardY, cw, cardH, 12, theme::CCARD);
            txt(x + 12, cardY + 8, label, 11, theme::CDIM);
            settextstyle(19, 0, L"微软雅黑"); settextcolor(vc); setbkmode(TRANSPARENT);
            outtextxy(x + 12, cardY + 28, value);
        };
        swprintf(buf, 512, L"%.2f s", b.t);
        card(0, L"时间 / 状态", buf, b.resting ? theme::CACCENT : theme::CINK);
        swprintf(buf, 512, L"%.2f m/s", b.speed());
        card(1, L"当前速度", buf, theme::CBLUE);
        swprintf(buf, 512, L"%.2f m/s", b.peakSpeed);
        card(2, L"峰值速度", buf, theme::CORANGE);
        swprintf(buf, 512, L"%.3f MJ", b.kinetic() / 1e6);
        card(3, L"动能", buf, theme::CACCENT);
        // 次级信息行
        const int ly = H - 52;
        swprintf(buf, 512, L"位置 x=%.1f y=%.1f z=%.1f", b.x, b.y, b.z);
        txt(PW + 14, ly, buf, 12, theme::CMUTED);
        swprintf(buf, 512, L"地形 坡段%d 起伏%.2f 沟谷%d 落差%.0fm", cfg.terrain.segments, cfg.terrain.roughness, cfg.terrain.valleys, cfg.terrain.drop);
        txt(PW + 250, ly, buf, 12, theme::CDIM);
        if (cfg.barrier.kind != BarrierKind::None) {
            swprintf(buf, 512, L"挡板 %s 挠度%.2fm 峰值力%.0fkN", b.bs.failed ? L"失效" : (b.bs.intercepted ? L"截停" : (b.bs.passed ? L"越顶" : (b.bs.hit ? L"接触中" : L"未接触"))), b.bs.d, b.bs.peakForce / 1000);
            txt(PW + 560, ly, buf, 12, theme::CORANGE);
        }
        if (hasBatch) {
            swprintf(buf, 512, L"批量%zu次 · 0–10m %zu(%.1f%%) · 截停%zu 越顶%zu", batch.runs.size(), batch.inBand,
                     batch.runs.size() ? 100.0 * batch.inBand / batch.runs.size() : 0, batch.intercepted, batch.overtopped);
            txt(PW + 14, H - 30, buf, 12, theme::CACCENT);
        } else {
            txt(PW + 14, H - 30, note.c_str(), 12, theme::CDIM);
        }
        swprintf(buf, 512, L"%.0f fps", fps);
        settextstyle(12, 0, L"微软雅黑"); settextcolor(theme::CDIM); setbkmode(TRANSPARENT);
        outtextxy(W - 70, H - 30, buf);
        FlushBatchDraw();

        if (shot && ++frames >= 480) {
            saveimage(shotPath);
            std::printf("t=%.2f x=%.1f y=%.1f speed=%.2f peak=%.2f rest=%d barrier=%d hit=%d passed=%d failed=%d\n",
                        b.t, b.x, b.y, b.speed(), b.peakSpeed, (int)b.resting, (int)cfg.barrier.kind,
                        (int)b.bs.hit, (int)b.bs.passed, (int)b.bs.failed);
            EndBatchDraw(); closegraph(); return 0;
        }
        Sleep(1);
    }
    if (fut.valid()) fut.wait();
    EndBatchDraw(); closegraph(); return 0;
}

inline void help() {
    std::printf("rollingstone3d_full | 3D rockfall with statistics/barriers/studies | SI units\n"
        "No args: interactive viewer. --shot out.png [seed]: auto demo + screenshot.\n"
        "--batch --runs 200 --seed 42 --cv 0.1 --out DIR\n"
        "--study / --terrain-study --runs N --out DIR\n"
        "--barrier none|rigid|elastic --barrier-x 100 --barrier-height 6 --barrier-restitution 0.15\n"
        "--stiffness 2e6 --damping 2e5 --barrier-mass 1000 --max-deflection 3\n"
        "--terrain-seed 7 --segments 10 --roughness 0.6 --valleys 2 --valley-width 1.0\n"
        "--speed 2 --radius 1.6 --roll 0.02 --friction 0.45 --restitution 0.42 --drag 1\n"
        "--detector-x -1 (auto toe) --bin-height 0.5\n");
}

} // namespace r3f

int main(int argc, char** argv) {
    using namespace r3f;
    Config c;
    std::string mode = "viewer", out, shotPath;
    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto val = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : ""; };
        auto num = [&]() { return std::atof(val().c_str()); };
        if (k == "--help" || k == "-h") { help(); return 0; }
        else if (k == "--batch") mode = "batch";
        else if (k == "--study") mode = "study";
        else if (k == "--terrain-study") mode = "terrain";
        else if (k == "--shot") { mode = "shot"; shotPath = val(); }
        else if (k == "--runs") c.runs = (int)num();
        else if (k == "--seed") c.seed = (uint32_t)num();
        else if (k == "--cv") c.cv = num();
        else if (k == "--out") out = val();
        else if (k == "--barrier") { std::string v = val(); c.barrier.kind = v == "rigid" ? BarrierKind::Rigid : (v == "elastic" ? BarrierKind::Elastic : BarrierKind::None); }
        else if (k == "--barrier-x") c.barrier.x = num();
        else if (k == "--barrier-height") c.barrier.height = num();
        else if (k == "--barrier-restitution") c.barrier.restitution = num();
        else if (k == "--stiffness") c.barrier.stiffness = num();
        else if (k == "--damping") c.barrier.damping = num();
        else if (k == "--barrier-mass") c.barrier.mass = num();
        else if (k == "--max-deflection") c.barrier.maxDeflection = num();
        else if (k == "--terrain-seed") c.terrainSeed = (uint32_t)num();
        else if (k == "--segments") c.terrain.segments = (int)num();
        else if (k == "--roughness") c.terrain.roughness = num();
        else if (k == "--valleys") c.terrain.valleys = (int)num();
        else if (k == "--valley-width") c.terrain.valleyWidth = num();
        else if (k == "--speed") c.initialSpeed = num();
        else if (k == "--radius") c.radius = num();
        else if (k == "--roll") c.rollRes = num();
        else if (k == "--friction") c.friction = num();
        else if (k == "--restitution") c.restitution = num();
        else if (k == "--drag") c.drag = num();
        else if (k == "--detector-x") c.detectorX = num();
        else if (k == "--bin-height") c.binHeight = num();
        else if (!k.empty() && k[0] != '-' && mode == "shot") c.terrainSeed = (uint32_t)std::strtoul(k.c_str(), nullptr, 10);
        else { std::fprintf(stderr, "Unknown option: %s\n", k.c_str()); return 2; }
    }
    try {
        if (mode == "viewer" || mode == "shot") {
            std::wstring wp(shotPath.begin(), shotPath.end());
            return runViewer(c, mode == "shot" ? wp.c_str() : nullptr);
        }
        c.finalize();
        if (out.empty()) out = "results3d/run_" + std::to_string(std::chrono::system_clock::now().time_since_epoch().count() / 1000000000);
        if (mode == "batch") {
            Batch b = runBatch(c, c.runs, c.seed, true);
            exportBatch(b, out);
            auto ci = wilson(b.inBand, b.runs.size());
            std::printf("runs=%zu reached=%zu in_0_10m=%zu P=%.4f Wilson95%%=[%.4f,%.4f]\n",
                        b.runs.size(), b.reached, b.inBand, (double)b.inBand / b.runs.size(), ci.first, ci.second);
            if (c.barrier.kind != BarrierKind::None)
                std::printf("intercepted=%zu overtopped=%zu failed_passage=%zu unresolved=%zu\n",
                            b.intercepted, b.overtopped, b.failedPass, b.hitUnresolved);
            std::printf("Saved: %s\n", fs::absolute(out).string().c_str());
        }
        else if (mode == "study") { int n = runParameterStudy(c, std::max(30, c.runs), c.seed, out); std::printf("parameter scenarios=%d Saved: %s\n", n, fs::absolute(out).string().c_str()); }
        else if (mode == "terrain") { int n = runTerrainStudy(c, std::max(30, c.runs), c.seed, out); std::printf("terrain scenarios=%d Saved: %s\n", n, fs::absolute(out).string().c_str()); }
        return 0;
    } catch (const std::exception& e) { std::fprintf(stderr, "Error: %s\n", e.what()); return 2; }
}
