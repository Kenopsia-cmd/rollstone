// rollingstone3d.cpp —— 三维滚石物理演示（独立文件，不影响二维版）
//
// 物理：球在三维高度场 z=h(x,y) 上运动。接触时重力沿坡面切向分解驱动，
//       法向速度按恢复系数弹跳，切向按库仑摩擦限幅，滚动阻力沿速度反向；
//       球可沿任意方向运动（侧向偏移、沿沟谷汇聚），是真三维而非剖面动画。
// 渲染：EasyX 纯软件渲染（画家算法地形网格 + Lambert 光照 + 球体径向渐变 + 轨迹尾迹）。
// 交互：鼠标拖动旋转视角，滚轮缩放；Space 暂停，R 重置，N 换地形，
//       T 尾迹开关，+/- 播放倍速，Esc 退出。
// 编译：cl /std:c++17 /EHsc /utf-8 rollingstone3d.cpp /link EasyXw.lib user32.lib gdi32.lib
//       （VS 开发者命令行；EasyX 头/库随 VS 或单独安装）
// 命令行：rollingstone3d.exe --shot out.png [seed]  无窗口自动演示并截图（用于验证）
//
// 单位：SI（m, s, kg）。教学演示，未经工程标定。

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
#include <algorithm>
#include <vector>
#include <string>
#include <chrono>

#pragma comment(lib, "user32.lib")
#pragma comment(lib, "gdi32.lib")

namespace r3 {

constexpr double G = 9.81;
constexpr int W = 1280, H = 800;

// ---------------- 地形高度场 ----------------
// 随机三维地形生成器：主坡为随机折线坡段（上陡下缓可调），横向叠加 1–3 条
// 随机沟谷（位置/宽度/深度/摆动均由种子决定）与随机起伏。同种子同参数完全可复现。
struct Terrain {
    // —— 可调生成参数 ——
    int segments = 10;          // 主坡折线段数 3..24
    double roughness = 0.6;     // 起伏强度 0..1
    int valleys = 2;            // 沟谷数量 0..3
    double valleyWidth = 1.0;   // 沟谷宽度倍率 0.4..2.0
    double length = 120.0, width = 60.0;

    // —— 由种子生成的形态 ——
    uint32_t seedUsed = 0;
    double drop = 52.0;
    double profileExp = 0.82;                 // 主坡凹凸指数
    std::vector<double> nodeZ;                // 坡段节点高程（含首尾）
    std::vector<double> nodeJitter;           // 横向倾斜扰动
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

        drop = 40.0 + rnd() * 24.0;
        profileExp = 0.62 + rnd() * 0.55;
        bumpPhase = rnd() * 6.2831853;
        bumpAmp = 0.6 + rnd() * 1.1;

        // 主坡折线：按幂律分配落差，再对每段施加随机扰动（段数越多越碎）
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
        // 保证整体下降趋势（避免局部反坡把球困住）
        for (int i = 1; i <= segments; ++i)
            nodeZ[i] = std::min(nodeZ[i], nodeZ[i - 1] - drop * 0.012);

        // 随机沟谷
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
    }

    void reseed(uint32_t seed) { generate(seed); }

    // 主坡高程（折线插值）与横向倾斜
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
        // 沟谷：抛物线截面，中心线随 x 正弦摆动
        for (const auto& v : vlist) {
            const double c = v.y0 + v.amp * std::sin(x * v.freq + v.phase);
            const double d = (y - c) / v.width;
            z += v.depth * v.width * v.width * d * d * 0.5;
        }
        // 横向整体倾斜（让球有侧向初驱）
        z += tilt * (y / (width * 0.5)) * 0.35;
        // 随机起伏
        const double bumps = bumpAmp * roughness *
            (1.7 * std::sin(x * 0.13 + bumpPhase) * std::cos(y * 0.16)
           + 0.9 * std::sin(x * 0.31 + 1.7) * std::sin(y * 0.27 + bumpPhase));
        // 末段平台，便于停稳
        const double flat = (t > 0.9) ? (t - 0.9) * 10.0 : 0.0;
        return z * (1.0 - 0.45 * flat) + bumps * (1.0 - 0.6 * flat);
    }

    // 数值梯度（中心差分），用于法向与切向分解
    void grad(double x, double y, double& gx, double& gy) const {
        const double e = 0.05;
        gx = (height(x + e, y) - height(x - e, y)) / (2 * e);
        gy = (height(x, y + e) - height(x, y - e)) / (2 * e);
    }

    // 起点横向位置：放在最主沟谷的一侧，便于观察侧向汇聚
    double startY() const {
        if (vlist.empty()) return width * 0.18;
        const auto& v = vlist.front();
        const double c = v.y0 + v.amp * std::sin(5.0 * v.freq + v.phase);
        return std::clamp(c + v.width * 0.85, -width * 0.42, width * 0.42);
    }
};

// ---------------- 球体状态与物理 ----------------
struct Ball {
    double x = 6, y = 0, z = 0;
    double vx = 0, vy = 0, vz = 0;
    double radius = 1.6;
    double mass = 2600.0 * 4.0 / 3.0 * 3.14159265358979 * 1.6 * 1.6 * 1.6;
    double rollRes = 0.02, friction = 0.45, rest = 0.42, drag = 1.0;
    double spin = 0;           // 视觉滚动角
    bool resting = false;
    double restTime = 0;
    double t = 0;
    double peakSpeed = 0, distance = 0;

    void place(const Terrain& T, double y0) {
        x = 5.0; y = y0;
        z = T.height(x, y) + radius;
        vx = 2.0; vy = 0; vz = 0;
        resting = false; restTime = 0; t = 0;
        peakSpeed = 0; distance = 0; spin = 0;
    }

    // 单步：半隐式欧拉 + 接触冲量。dt 建议 <= 0.004
    void step(const Terrain& T, double dt) {
        if (resting) { t += dt; return; }
        const double px = x, py = y;
        // 重力
        vz -= G * dt;
        // 空气阻力（二次，隐式化处理避免大步长发散）
        const double sp = std::sqrt(vx * vx + vy * vy + vz * vz);
        if (sp > 1e-9 && drag > 0) {
            const double area = 3.14159265358979 * radius * radius;
            const double kd = 0.5 * 1.225 * drag * area / mass;
            const double k = 1.0 / (1.0 + kd * sp * dt);
            vx *= k; vy *= k; vz *= k;
        }
        x += vx * dt; y += vy * dt; z += vz * dt;

        // 边界：横向软反弹（山谷两侧），纵向出界即停
        if (y < -T.width * 0.5 + radius) { y = -T.width * 0.5 + radius; vy = std::abs(vy) * 0.3; }
        if (y >  T.width * 0.5 - radius) { y =  T.width * 0.5 - radius; vy = -std::abs(vy) * 0.3; }

        const double h = T.height(x, y);
        const double gap = z - radius - h;
        if (gap <= 0.0) {
            // 法向（指向坡面外）
            double gx, gy; T.grad(x, y, gx, gy);
            double nx = -gx, ny = -gy, nz = 1.0;
            const double nl = std::sqrt(nx * nx + ny * ny + nz * nz);
            nx /= nl; ny /= nl; nz /= nl;
            // 位置修正：贴回坡面
            z = h + radius;
            const double vn = vx * nx + vy * ny + vz * nz;
            const bool impact = (vn < -0.6);
            if (vn < 0) {
                // 法向：足够快才弹跳，否则视为持续接触
                const double e = impact ? rest : 0.0;
                const double jn = -(1.0 + e) * vn;
                vx += jn * nx; vy += jn * ny; vz += jn * nz;
            }
            if (impact) {
                // 真实碰撞：切向按库仑冲量限幅（实心球滚动等效 2/7）
                const double vn2 = vx * nx + vy * ny + vz * nz;
                double tvx = vx - vn2 * nx, tvy = vy - vn2 * ny, tvz = vz - vn2 * nz;
                const double tv = std::sqrt(tvx * tvx + tvy * tvy + tvz * tvz);
                if (tv > 1e-9) {
                    const double jt = std::min(tv * (2.0 / 7.0), friction * std::abs(vn) * (1.0 + rest));
                    vx -= jt * tvx / tv; vy -= jt * tvy / tv; vz -= jt * tvz / tv;
                }
            } else {
                // 持续接触（纯滚动）：重力切向分量按 5/7 折减，
                // 等价于 a = (5/7) g sinθ；此处扣除已全额施加的 2/7 部分。
                const double gtn = (-G) * nz;              // g_vec·n
                const double gtx = -gtn * nx;
                const double gty = -gtn * ny;
                const double gtz = -G - gtn * nz;
                vx -= (2.0 / 7.0) * gtx * dt;
                vy -= (2.0 / 7.0) * gty * dt;
                vz -= (2.0 / 7.0) * gtz * dt;
                // 陡坡打滑：坡度超过摩擦角时不再维持纯滚动，全额切向加速
                const double slopeNow = std::sqrt(gx * gx + gy * gy);
                if (slopeNow > friction) {
                    vx += (2.0 / 7.0) * gtx * dt;
                    vy += (2.0 / 7.0) * gty * dt;
                    vz += (2.0 / 7.0) * gtz * dt;
                }
            }
            // 滚动阻力：沿速度反向的减速（接触时）
            const double sp2 = std::sqrt(vx * vx + vy * vy + vz * vz);
            if (sp2 > 1e-9) {
                const double dec = rollRes * G * dt;
                const double k = std::max(0.0, 1.0 - dec / sp2);
                vx *= k; vy *= k; vz *= k;
            }
            spin += sp2 * dt / radius;
            // 停稳判据：速度足够小且坡度小于滚动阻力角
            const double slope = std::sqrt(gx * gx + gy * gy);
            if (sp2 < 0.12 && slope <= rollRes + 0.02) {
                restTime += dt;
                if (restTime > 0.6) { resting = true; vx = vy = vz = 0; }
            } else restTime = 0;
        }
        t += dt;
        const double sp3 = std::sqrt(vx * vx + vy * vy + vz * vz);
        peakSpeed = std::max(peakSpeed, sp3);
        distance += std::hypot(x - px, y - py);
        if (x > T.length + 4 || t > 90.0) { resting = true; vx = vy = vz = 0; }
    }

    double speed() const { return std::sqrt(vx * vx + vy * vy + vz * vz); }
    double kinetic() const { return 0.5 * mass * speed() * speed() * (1.0 + 2.0 / 5.0); }
};

// ---------------- 相机与投影 ----------------
struct Camera {
    double yaw = -0.90, pitch = 0.46, dist = 190.0;
    double cx = 60, cy = 0, cz = 22;     // 注视点
    double fov = 880.0;

    // 世界 -> 屏幕（右手系：x 前进，y 横向，z 向上）
    bool project(double X, double Y, double Z, double& sx, double& sy, double& depth) const {
        const double dx = X - cx, dy = Y - cy, dz = Z - cz;
        const double cy_ = std::cos(yaw), sy_ = std::sin(yaw);
        const double cp = std::cos(pitch), sp = std::sin(pitch);
        // 绕 z 旋转 yaw
        const double rx = dx * cy_ + dy * sy_;
        const double ry = -dx * sy_ + dy * cy_;
        // 绕 x 旋转 pitch
        const double rz = dz * cp - ry * sp;
        const double ry2 = dz * sp + ry * cp;
        const double zc = ry2 + dist;      // 相机在 -y 方向 dist 处
        if (zc < 6.0) return false;
        depth = zc;
        sx = W * 0.5 + fov * rx / zc;
        sy = H * 0.52 - fov * rz / zc;
        return true;
    }
};

inline COLORREF shade(double k, int r, int g, int b) {
    k = std::clamp(k, 0.0, 1.0);
    return RGB((int)(r * k), (int)(g * k), (int)(b * k));
}

// ---------------- 渲染 ----------------
struct TrailPt { double x, y, z; };

void renderTerrain(const Terrain& T, const Camera& cam) {
    const int NX = 46, NY = 26;
    struct Cell { double depth; POINT pts[4]; COLORREF col; };
    std::vector<Cell> cells;
    cells.reserve(NX * NY);
    // 光源方向（世界系）
    const double lx = -0.45, ly = -0.55, lz = 0.70;
    const double ll = std::sqrt(lx * lx + ly * ly + lz * lz);
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            const double x0 = T.length * i / NX, x1 = T.length * (i + 1) / NX;
            const double y0 = -T.width * 0.5 + T.width * j / NY;
            const double y1 = -T.width * 0.5 + T.width * (j + 1) / NY;
            const double xs[4] = { x0, x1, x1, x0 }, ys[4] = { y0, y0, y1, y1 };
            Cell c{}; double dsum = 0; bool ok = true;
            for (int k = 0; k < 4; ++k) {
                double sx, sy, dp;
                if (!cam.project(xs[k], ys[k], T.height(xs[k], ys[k]), sx, sy, dp)) { ok = false; break; }
                c.pts[k] = { (LONG)sx, (LONG)sy };
                dsum += dp;
            }
            if (!ok) continue;
            c.depth = dsum / 4;
            const double xm = (x0 + x1) * 0.5, ym = (y0 + y1) * 0.5;
            double gx, gy; T.grad(xm, ym, gx, gy);
            double nx = -gx, ny = -gy, nz = 1.0;
            const double nl = std::sqrt(nx * nx + ny * ny + nz * nz);
            const double lam = std::max(0.0, (nx * lx + ny * ly + nz * lz) / (nl * ll));
            const double hgt = T.height(xm, ym) / std::max(T.drop, 1.0);
            // 高处偏岩灰、低处偏土黄，叠加光照
            const int br = (int)(104 + 46 * (1 - hgt)), bg = (int)(112 + 30 * (1 - hgt)), bb = (int)(120 - 18 * (1 - hgt));
            c.col = shade(0.42 + 0.58 * lam, br, bg, bb);
            cells.push_back(c);
        }
    }
    std::sort(cells.begin(), cells.end(), [](const Cell& a, const Cell& b) { return a.depth > b.depth; });
    setlinestyle(PS_SOLID, 1);
    for (const auto& c : cells) {
        setfillcolor(c.col);
        setlinecolor(c.col);
        fillpolygon(c.pts, 4);
    }
}

void renderTrail(const std::vector<TrailPt>& trail, const Camera& cam) {
    if (trail.size() < 2) return;
    setlinestyle(PS_SOLID, 2);
    for (size_t i = 1; i < trail.size(); ++i) {
        double ax, ay, ad, bx, by, bd;
        if (!cam.project(trail[i - 1].x, trail[i - 1].y, trail[i - 1].z, ax, ay, ad)) continue;
        if (!cam.project(trail[i].x, trail[i].y, trail[i].z, bx, by, bd)) continue;
        const double f = (double)i / trail.size();
        setlinecolor(shade(0.55 + 0.45 * f, 90, 200, 235));
        line((int)ax, (int)ay, (int)bx, (int)by);
    }
}

void renderBall(const Ball& B, const Terrain& T, const Camera& cam) {
    // 阴影：投到坡面
    double sx, sy, sd;
    if (cam.project(B.x, B.y, T.height(B.x, B.y) + 0.05, sx, sy, sd)) {
        const double r = cam.fov * B.radius / sd;
        setfillcolor(RGB(24, 32, 44));
        solidellipse((int)(sx - r), (int)(sy - r * 0.42), (int)(sx + r), (int)(sy + r * 0.42));
    }
    double bx, by, bd;
    if (!cam.project(B.x, B.y, B.z, bx, by, bd)) return;
    const double r = std::max(3.0, cam.fov * B.radius / bd);
    // 径向渐变球体（同心圆近似）
    const int steps = 14;
    for (int i = steps; i >= 1; --i) {
        const double f = (double)i / steps;
        const double rr = r * f;
        const double ox = -r * 0.30 * (1 - f), oy = -r * 0.34 * (1 - f);
        const double k = 0.40 + 0.60 * (1 - f);
        setfillcolor(shade(k, 226, 138, 74));
        solidellipse((int)(bx + ox - rr), (int)(by + oy - rr), (int)(bx + ox + rr), (int)(by + oy + rr));
    }
    // 滚动条纹（体现自转）
    setlinecolor(RGB(120, 66, 30));
    setlinestyle(PS_SOLID, 2);
    const double a = B.spin;
    line((int)(bx + std::cos(a) * r * 0.72), (int)(by + std::sin(a) * r * 0.72),
         (int)(bx - std::cos(a) * r * 0.72), (int)(by - std::sin(a) * r * 0.72));
}

void text(int x, int y, const wchar_t* s, int size, COLORREF c) {
    settextcolor(c);
    settextstyle(size, 0, L"微软雅黑");
    setbkmode(TRANSPARENT);
    outtextxy(x, y, s);
}

void renderHUD(const Ball& B, const Terrain& T, uint32_t seed, bool paused,
               double speedScale, bool trailOn, double fps) {
    wchar_t buf[256];
    setfillcolor(RGB(14, 26, 44));
    solidrectangle(0, 0, W, 54);
    text(20, 14, L"三维滚石 · 高度场物理演示", 22, RGB(238, 242, 247));
    swprintf(buf, 256, L"seed %u   %s   %.2fx", seed, paused ? L"暂停" : L"运行", speedScale);
    text(W - 330, 20, buf, 15, RGB(150, 190, 225));

    setfillcolor(RGB(16, 30, 50));
    solidrectangle(0, H - 122, W, H);
    swprintf(buf, 256, L"随机地形  坡段 %d   起伏 %.2f   沟谷 %d   谷宽 %.1fx   落差 %.1f m   [ ] 段数   - = 起伏   V 沟谷   B 谷宽",
             T.segments, T.roughness, T.valleys, T.valleyWidth, T.drop);
    text(20, H - 114, buf, 14, RGB(126, 196, 178));
    swprintf(buf, 256, L"t = %.2f s    速度 %.2f m/s    峰值 %.2f m/s    动能 %.3f MJ",
             B.t, B.speed(), B.peakSpeed, B.kinetic() / 1e6);
    text(20, H - 88, buf, 16, RGB(226, 234, 242));
    swprintf(buf, 256, L"位置 x=%.1f  y=%.1f  z=%.1f    水平行程 %.1f m    %s",
             B.x, B.y, B.z, B.distance, B.resting ? L"已停稳" : L"运动中");
    text(20, H - 62, buf, 15, RGB(150, 190, 225));
    swprintf(buf, 256, L"拖动旋转 · 滚轮缩放 · Space 暂停 · R 重置 · N 换地形 · T 尾迹(%s) · , . 倍速 · Esc 退出   %.0f fps",
             trailOn ? L"开" : L"关", fps);
    text(20, H - 34, buf, 13, RGB(112, 146, 180));
}

} // namespace r3

// ---------------- 主程序 ----------------
int main(int argc, char** argv) {
    using namespace r3;
    bool shot = false; std::string shotPath = "shot3d.png"; uint32_t seed = 7;
    int seg = 10, val = 2; double rough = 0.6, vw = 1.0;
    for (int i = 1; i < argc; ++i) {
        auto next = [&](double def) { return (i + 1 < argc) ? std::atof(argv[++i]) : def; };
        if (!std::strcmp(argv[i], "--shot") && i + 1 < argc) { shot = true; shotPath = argv[++i]; }
        else if (!std::strcmp(argv[i], "--segments")) seg = (int)next(10);
        else if (!std::strcmp(argv[i], "--roughness")) rough = next(0.6);
        else if (!std::strcmp(argv[i], "--valleys")) val = (int)next(2);
        else if (!std::strcmp(argv[i], "--valley-width")) vw = next(1.0);
        else seed = (uint32_t)std::strtoul(argv[i], nullptr, 10);
    }

    Terrain T;
    T.segments = seg; T.roughness = rough; T.valleys = val; T.valleyWidth = vw;
    T.generate(seed);
    Ball B; B.place(T, T.startY());
    Camera cam;
    std::vector<TrailPt> trail;
    bool paused = false, trailOn = true;
    double speedScale = 1.0;

    initgraph(W, H);
    setbkcolor(RGB(9, 17, 30));
    BeginBatchDraw();

    auto last = std::chrono::steady_clock::now();
    double fps = 60.0;
    int frames = 0;
    const double dt = 0.002;

    while (true) {
        auto now = std::chrono::steady_clock::now();
        const double el = std::min(0.05, std::chrono::duration<double>(now - last).count());
        last = now;
        fps = fps * 0.9 + (el > 1e-6 ? 1.0 / el : 60.0) * 0.1;

        if (!shot) {
            ExMessage m;
            static bool dragging = false; static int lx = 0, ly = 0;
            while (peekmessage(&m, EX_MOUSE | EX_KEY)) {
                if (m.message == WM_LBUTTONDOWN) { dragging = true; lx = m.x; ly = m.y; }
                else if (m.message == WM_LBUTTONUP) dragging = false;
                else if (m.message == WM_MOUSEMOVE && dragging) {
                    cam.yaw += (m.x - lx) * 0.006;
                    cam.pitch = std::clamp(cam.pitch + (m.y - ly) * 0.004, 0.08, 1.35);
                    lx = m.x; ly = m.y;
                }
                else if (m.message == WM_MOUSEWHEEL) {
                    cam.dist = std::clamp(cam.dist - m.wheel * 0.06, 60.0, 420.0);
                }
                else if (m.message == WM_KEYDOWN) {
                    auto regen = [&]() {
                        T.generate(seed);
                        B.place(T, T.startY());
                        trail.clear(); paused = false;
                    };
                    switch (m.vkcode) {
                    case VK_ESCAPE: EndBatchDraw(); closegraph(); return 0;
                    case VK_SPACE: paused = !paused; break;
                    case 'R': B.place(T, T.startY()); trail.clear(); paused = false; break;
                    case 'N': seed += 1; regen(); break;
                    case 'T': trailOn = !trailOn; if (!trailOn) trail.clear(); break;
                    case VK_OEM_4: T.segments = std::max(3, T.segments - 1); regen(); break;   // [
                    case VK_OEM_6: T.segments = std::min(24, T.segments + 1); regen(); break;  // ]
                    case VK_OEM_MINUS: T.roughness = std::max(0.0, T.roughness - 0.1); regen(); break;
                    case VK_OEM_PLUS: T.roughness = std::min(1.0, T.roughness + 0.1); regen(); break;
                    case 'V': T.valleys = (T.valleys + 1) % 4; regen(); break;
                    case 'B': T.valleyWidth = (T.valleyWidth >= 2.0) ? 0.5 : T.valleyWidth + 0.25; regen(); break;
                    case VK_OEM_COMMA: speedScale = std::max(0.25, speedScale / 2); break;
                    case VK_OEM_PERIOD: speedScale = std::min(4.0, speedScale * 2); break;
                    }
                }
            }
        }

        if (!paused) {
            int steps = (int)std::lround(el * speedScale / dt);
            steps = std::clamp(steps, 1, 400);
            for (int i = 0; i < steps; ++i) {
                B.step(T, dt);
                if (trailOn && (i % 6 == 0)) {
                    trail.push_back({ B.x, B.y, B.z });
                    if (trail.size() > 1400) trail.erase(trail.begin(), trail.begin() + 200);
                }
            }
        }
        // 相机跟随（缓动）
        // 相机缓动跟随：注视点在球与坡体中部之间取平衡，保证整坡可见
        const double tx = 0.45 * B.x + 0.55 * (T.length * 0.5);
        cam.cx += (tx - cam.cx) * 0.03;
        cam.cy += (B.y * 0.5 - cam.cy) * 0.03;
        cam.cz += (T.drop * 0.34 - cam.cz) * 0.03;

        cleardevice();
        renderTerrain(T, cam);
        if (trailOn) renderTrail(trail, cam);
        renderBall(B, T, cam);
        renderHUD(B, T, seed, paused, speedScale, trailOn, fps);
        FlushBatchDraw();

        if (shot) {
            if (++frames >= 480) {   // 约 8 秒演示后截图
                saveimage(std::wstring(shotPath.begin(), shotPath.end()).c_str());
                std::printf("seed=%u seg=%d rough=%.2f val=%d vw=%.2f drop=%.1f | t=%.2f x=%.1f y=%.1f speed=%.2f peak=%.2f rest=%d\n",
                            seed, T.segments, T.roughness, T.valleys, T.valleyWidth, T.drop,
                            B.t, B.x, B.y, B.speed(), B.peakSpeed, (int)B.resting);
                EndBatchDraw(); closegraph(); return 0;
            }
        }
        Sleep(1);
    }
}
