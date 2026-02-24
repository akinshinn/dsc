#include <vector>
#include <array>
#include <map>
#include <functional>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <locale>

// ================================================================
// Выбор задачи: задайте номер здесь или передайте -DTASK=2 компилятору
// ================================================================
#ifndef TASK
#  define TASK 2
#endif
// ================================================================

// ---- структуры сетки ----

struct Node { double x, y; };
struct Element { std::array<int, 3> n; };
struct BC { int node, type; double x, y; };

struct Mesh {
    std::vector<Node>    nodes;
    std::vector<Element> elems;
};

Mesh make_mesh(double a, double b, double c, double d, int nx, int ny)
{
    Mesh m;
    m.nodes.resize((nx + 1) * (ny + 1));
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i)
            m.nodes[j * (nx + 1) + i] = { a + i * (b - a) / nx, c + j * (d - c) / ny };

    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            int n0 = j * (nx + 1) + i, n1 = n0 + 1, n2 = n0 + (nx + 1), n3 = n2 + 1;
            m.elems.push_back({ {n0,n1,n3} });
            m.elems.push_back({ {n0,n3,n2} });
        }
    return m;
}

// ================================================================
// ЗАДАЧА 1
//   lambda = 1, q = 0
//   u_exact = x^2 - y^2 + 2
//
//   ГУ по граням:
//     x=a=0 (левая):  тип 1,  u = -y^2 + 2
//     x=b=1 (правая): тип 1,  u = 1 - y^2 + 2
//     y=c=0 (нижняя): тип 1,  u = x^2 + 2
//     y=d=1 (верхняя):тип 2,  -(du/dy) = 2y|_{y=1} = 2
// ================================================================

#if TASK == 1

double lambda_(double, double) { return 1.0; }
double q_src(double, double) { return 0.0; }
double u_exact(double x, double y) { return x * x - y * y + 2.0; }
double u_dir(double x, double y) { return u_exact(x, y); }

// -(lambda * du/dn):
//   верхняя грань y=d=1: нормаль (0,+1), du/dn = du/dy = -2y => q = -(-2y) = 2y
double q_neu(double x, double y) {
    (void)x;
    return 2.0 * y;   // используется только на верхней грани y=1 => 2*1 = 2
}

std::vector<BC> boundary_table(const Mesh& m,
    double a, double b, double c, double d, double tol = 1e-10)
{
    std::vector<BC> bcs;
    for (int i = 0; i < (int)m.nodes.size(); ++i) {
        double x = m.nodes[i].x, y = m.nodes[i].y;
        bool left = std::abs(x - a) < tol;
        bool right = std::abs(x - b) < tol;
        bool bottom = std::abs(y - c) < tol;
        bool top = std::abs(y - d) < tol;
        if (left)   bcs.push_back({ i, 1, x, y });
        else if (right)  bcs.push_back({ i, 1, x, y });
        else if (bottom) bcs.push_back({ i, 1, x, y });
        else if (top)    bcs.push_back({ i, 2, x, y });  // Нейман
    }
    return bcs;
}

const char* TASK_NAME = "Task 1: u = x^2 - y^2 + 2, lambda=1, q=0";
const char* BC_DESC_LEFT = "x=0: Dirichlet  u = -y^2 + 2";
const char* BC_DESC_RIGHT = "x=1: Dirichlet  u = 1 - y^2 + 2";
const char* BC_DESC_BOT = "y=0: Dirichlet  u = x^2 + 2";
const char* BC_DESC_TOP = "y=1: Neumann   -(du/dn) = 2";
constexpr double TEST_TOL = 5e-2;

#endif // TASK == 1

// ================================================================
// ЗАДАЧА 2
//   lambda = x + y,
//   q = -20*cos(5(x+y)) + 50*(x+y)*sin(5(x+y))
//   u_exact = sin(5(x+y)) + 2
//
//   ГУ по граням:
//     x=a=0 (левая):  тип 1,  u = sin(5y) + 2
//     x=b=1 (правая): тип 2,  -(lambda*du/dn) = -(x+y)*5*cos(5(x+y))|_{x=1}
//     y=c=0 (нижняя): тип 1,  u = sin(5x) + 2
//     y=d=1 (верхняя):тип 2,  -(lambda*du/dn) = -(x+y)*5*cos(5(x+y))|_{y=1}
// ================================================================

#if TASK == 2

double lambda_(double x, double y) { return x + y; }

double q_src(double x, double y) {
    double s = x + y;
    return -20.0 * std::cos(5.0 * s) + 50.0 * s * std::sin(5.0 * s);
}

double u_exact(double x, double y) { return std::sin(5.0 * (x + y)) + 2.0; }
double u_dir(double x, double y) { return u_exact(x, y); }

// -(lambda * du/dn):
//   du/dx = du/dy = 5*cos(5(x+y))
//   правая грань  x=b=1: нормаль (+1,0), q = -(x+y)*5*cos(5(x+y))
//   верхняя грань y=d=1: нормаль (0,+1), q = -(x+y)*5*cos(5(x+y))
double q_neu(double x, double y) {
    return -(x + y) * 5.0 * std::cos(5.0 * (x + y));
}

std::vector<BC> boundary_table(const Mesh& m,
    double a, double b, double c, double d, double tol = 1e-10)
{
    std::vector<BC> bcs;
    for (int i = 0; i < (int)m.nodes.size(); ++i) {
        double x = m.nodes[i].x, y = m.nodes[i].y;
        bool left = std::abs(x - a) < tol;
        bool right = std::abs(x - b) < tol;
        bool bottom = std::abs(y - c) < tol;
        bool top = std::abs(y - d) < tol;
        if (left)   bcs.push_back({ i, 1, x, y });  // Дирихле
        else if (right)  bcs.push_back({ i, 2, x, y });  // Нейман
        else if (bottom) bcs.push_back({ i, 1, x, y });  // Дирихле
        else if (top)    bcs.push_back({ i, 2, x, y });  // Нейман
    }
    return bcs;
}

const char* TASK_NAME = "Task 2: u = sin(5(x+y)) + 2, lambda=x+y";
const char* BC_DESC_LEFT = "x=0: Dirichlet  u = sin(5y) + 2";
const char* BC_DESC_RIGHT = "x=1: Neumann   -(lambda*du/dn) = -(x+y)*5*cos(5(x+y))";
const char* BC_DESC_BOT = "y=0: Dirichlet  u = sin(5x) + 2";
const char* BC_DESC_TOP = "y=1: Neumann   -(lambda*du/dn) = -(x+y)*5*cos(5(x+y))";
constexpr double TEST_TOL = 1e-1;

#endif // TASK == 2

// ---- разреженная матрица ----

struct DokMatrix {
    int n;
    std::map<std::pair<int, int>, double> data;
    explicit DokMatrix(int sz) : n(sz) {}
    void add(int i, int j, double v) { data[{i, j}] += v; }
    void set(int i, int j, double v) { data[{i, j}] = v; }
    double get(int i, int j) const {
        auto it = data.find({ i,j });
        return it != data.end() ? it->second : 0.0;
    }
};

struct CsrMatrix {
    int n = 0;
    std::vector<double> val;
    std::vector<int> col, row;
    std::vector<double> mv(const std::vector<double>& x) const {
        std::vector<double> y(n, 0.0);
        for (int i = 0; i < n; ++i)
            for (int k = row[i]; k < row[i + 1]; ++k)
                y[i] += val[k] * x[col[k]];
        return y;
    }
};

CsrMatrix to_csr(const DokMatrix& d)
{
    CsrMatrix c; c.n = d.n;
    c.row.resize(d.n + 1, 0);
    for (auto& [k, v] : d.data) c.row[k.first + 1]++;
    for (int i = 0; i < d.n; ++i) c.row[i + 1] += c.row[i];
    int nnz = c.row[d.n];
    c.val.resize(nnz); c.col.resize(nnz);
    std::vector<int> pos(c.row.begin(), c.row.end());
    for (auto& [k, v] : d.data) {
        int p = pos[k.first]++;
        c.val[p] = v; c.col[p] = k.second;
    }
    return c;
}

static double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0; for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i]; return s;
}

std::vector<double> pcg(const CsrMatrix& A,
    const std::vector<double>& b, double tol = 1e-10, int maxit = 200000)
{
    int n = A.n;
    std::vector<double> Minv(n, 1.0);
    for (int i = 0; i < n; ++i)
        for (int k = A.row[i]; k < A.row[i + 1]; ++k)
            if (A.col[k] == i && std::abs(A.val[k]) > 1e-300)
                Minv[i] = 1.0 / A.val[k];

    std::vector<double> x(n, 0.0), r = b, z(n), p(n);
    for (int i = 0; i < n; ++i) z[i] = Minv[i] * r[i];
    p = z;
    double rz = dot(r, z), bn = std::sqrt(dot(b, b));
    if (bn < 1e-300) return x;

    for (int it = 0; it < maxit; ++it) {
        auto Ap = A.mv(p);
        double pAp = dot(p, Ap);
        if (std::abs(pAp) < 1e-300) break;
        double alpha = rz / pAp;
        for (int i = 0; i < n; ++i) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }
        if (std::sqrt(dot(r, r)) < tol * bn) break;
        for (int i = 0; i < n; ++i) z[i] = Minv[i] * r[i];
        double rz2 = dot(r, z), beta = rz2 / rz;
        for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
        rz = rz2;
    }
    return x;
}

// ---- МКЭ ----

struct Ke { double K[3][3]; double F[3]; };

Ke local_matrix(const std::array<Node, 3>& p,
    const std::array<double, 3>& lam,
    const std::array<double, 3>& q)
{
    Ke ke{};
    double x1 = p[0].x, y1 = p[0].y, x2 = p[1].x, y2 = p[1].y, x3 = p[2].x, y3 = p[2].y;
    double D = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    double S = std::abs(D) / 2.0;
    if (S < 1e-15) return ke;

    double b[3] = { (y2 - y3) / D,(y3 - y1) / D,(y1 - y2) / D };
    double c[3] = { (x3 - x2) / D,(x1 - x3) / D,(x2 - x1) / D };
    double la = (lam[0] + lam[1] + lam[2]) / 3.0;
    double qa = (q[0] + q[1] + q[2]) / 3.0;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            ke.K[i][j] = la * S * (b[i] * b[j] + c[i] * c[j]);
    for (int i = 0; i < 3; ++i) ke.F[i] = qa * S / 3.0;
    return ke;
}

void assemble(const Mesh& mesh,
    const std::vector<double>& lam, const std::vector<double>& qv,
    DokMatrix& K, std::vector<double>& F)
{
    for (auto& el : mesh.elems) {
        std::array<Node, 3>   xy{ mesh.nodes[el.n[0]], mesh.nodes[el.n[1]], mesh.nodes[el.n[2]] };
        std::array<double, 3> la{ lam[el.n[0]], lam[el.n[1]], lam[el.n[2]] };
        std::array<double, 3> qe{ qv[el.n[0]], qv[el.n[1]], qv[el.n[2]] };
        auto ke = local_matrix(xy, la, qe);
        for (int a = 0; a < 3; ++a) {
            F[el.n[a]] += ke.F[a];
            for (int b = 0; b < 3; ++b) K.add(el.n[a], el.n[b], ke.K[a][b]);
        }
    }
}

void apply_dirichlet(DokMatrix& K, std::vector<double>& F,
    const std::vector<std::pair<int, double>>& dbc)
{
    std::map<int, double> m(dbc.begin(), dbc.end());
    for (auto& [key, v] : K.data) {
        int i = key.first, j = key.second;
        if (m.count(j) && !m.count(i)) F[i] -= v * m[j];
    }
    std::vector<std::pair<int, int>> rm;
    for (auto& [key, v] : K.data)
        if (m.count(key.first) || m.count(key.second)) rm.push_back(key);
    for (auto& k : rm) K.data.erase(k);
    for (auto& [idx, u] : m) { K.set(idx, idx, 1.0); F[idx] = u; }
}

void apply_neumann(const Mesh& mesh, const std::vector<BC>& bcs,
    std::vector<double>& F, std::function<double(double, double)> q0)
{
    std::map<int, bool> neu;
    for (auto& bc : bcs) if (bc.type == 2) neu[bc.node] = true;

    std::map<std::pair<int, int>, int> ecnt;
    for (auto& el : mesh.elems)
        for (int k = 0; k < 3; ++k) {
            int a = el.n[k], b = el.n[(k + 1) % 3];
            if (a > b) std::swap(a, b);
            ecnt[{a, b}]++;
        }

    for (auto& [e, cnt] : ecnt) {
        if (cnt != 1) continue;
        int i = e.first, j = e.second;
        if (!neu.count(i) || !neu.count(j)) continue;
        double xi = mesh.nodes[i].x, yi = mesh.nodes[i].y;
        double xj = mesh.nodes[j].x, yj = mesh.nodes[j].y;
        double L = std::sqrt((xj - xi) * (xj - xi) + (yj - yi) * (yj - yi));
        double qi = q0(xi, yi), qj = q0(xj, yj);
        F[i] -= L * (2 * qi + qj) / 6.0;
        F[j] -= L * (qi + 2 * qj) / 6.0;
    }
}

// ---- вывод ----

void write_csv(const std::string& fname, const Mesh& mesh, const std::vector<double>& u)
{
    std::ofstream f(fname);
    if (!f) throw std::runtime_error("failed to open " + fname);
    f << "x,y,u\n" << std::fixed << std::setprecision(10);
    for (int i = 0; i < (int)mesh.nodes.size(); ++i)
        f << mesh.nodes[i].x << ',' << mesh.nodes[i].y << ',' << u[i] << '\n';
}

void write_vtk(const std::string& fname, const Mesh& mesh, const std::vector<double>& u)
{
    std::ofstream f(fname);
    if (!f) throw std::runtime_error("failed to open " + fname);
    int nn = (int)mesh.nodes.size(), ne = (int)mesh.elems.size();
    f << "# vtk DataFile Version 3.0\nFEM\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    f << "POINTS " << nn << " double\n" << std::fixed << std::setprecision(8);
    for (auto& nd : mesh.nodes) f << nd.x << ' ' << nd.y << " 0\n";
    f << "\nCELLS " << ne << ' ' << 4 * ne << '\n';
    for (auto& el : mesh.elems)
        f << "3 " << el.n[0] << ' ' << el.n[1] << ' ' << el.n[2] << '\n';
    f << "\nCELL_TYPES " << ne << '\n';
    for (int i = 0; i < ne; ++i) f << "5\n";
    f << "\nPOINT_DATA " << nn << "\nSCALARS u double 1\nLOOKUP_TABLE default\n";
    for (int i = 0; i < nn; ++i) f << u[i] << '\n';
}

void write_txt(const std::string& fname, const Mesh& mesh,
    const std::vector<double>& u,
    double a, double b, double c, double d,
    int nx, int ny)
{
    std::ofstream f(fname);
    if (!f) throw std::runtime_error("failed to open " + fname);

    f << "================================================\n";
    f << "  " << TASK_NAME << "\n";
    f << "  Equation: div(lambda * grad(u)) + q = 0\n";
    f << "================================================\n\n";

    f << "Domain:   [" << a << ", " << b << "] x [" << c << ", " << d << "]\n";
    f << "Mesh:     " << nx << " x " << ny << "\n";
    f << "Nodes:    " << mesh.nodes.size() << "\n";
    f << "Elements: " << mesh.elems.size() << "\n\n";

    double umin = *std::min_element(u.begin(), u.end());
    double umax = *std::max_element(u.begin(), u.end());
    double uavg = 0; for (auto v : u) uavg += v; uavg /= u.size();

    f << std::fixed << std::setprecision(6);
    f << "u_min = " << umin << "\n";
    f << "u_max = " << umax << "\n";
    f << "u_avg = " << uavg << "\n\n";

    f << "Boundary conditions:\n";
    f << "  " << BC_DESC_LEFT << "\n";
    f << "  " << BC_DESC_RIGHT << "\n";
    f << "  " << BC_DESC_BOT << "\n";
    f << "  " << BC_DESC_TOP << "\n\n";

    f << "------------------------------------------------\n";
    f << "  Node       x           y           u\n";
    f << "------------------------------------------------\n";
    f << std::fixed << std::setprecision(8);
    for (int i = 0; i < (int)mesh.nodes.size(); ++i)
        f << "  " << std::setw(5) << i
        << "  " << std::setw(10) << mesh.nodes[i].x
        << "  " << std::setw(10) << mesh.nodes[i].y
        << "  " << std::setw(12) << u[i] << "\n";
    f << "------------------------------------------------\n";
}

// ---- решение ----

std::vector<double> solve(const Mesh& mesh, double a, double b, double c, double d)
{
    int N = (int)mesh.nodes.size();
    std::vector<double> lam(N), qv(N);
    for (int i = 0; i < N; ++i) {
        lam[i] = lambda_(mesh.nodes[i].x, mesh.nodes[i].y);
        qv[i] = q_src(mesh.nodes[i].x, mesh.nodes[i].y);
    }
    DokMatrix K(N);
    std::vector<double> F(N, 0.0);
    assemble(mesh, lam, qv, K, F);

    auto bcs = boundary_table(mesh, a, b, c, d);
    apply_neumann(mesh, bcs, F, q_neu);

    std::vector<std::pair<int, double>> dbc;
    for (auto& bc : bcs)
        if (bc.type == 1) dbc.push_back({ bc.node, u_dir(bc.x,bc.y) });
    apply_dirichlet(K, F, dbc);

    return pcg(to_csr(K), F);
}

bool run_test()
{
    double a = 0, b = 1, c = 0, d = 1;
#if TASK == 1
    Mesh mesh = make_mesh(a, b, c, d, 20, 20);
#else
    Mesh mesh = make_mesh(a, b, c, d, 60, 60); // задача 2: мельче сетка из-за осцилляций
#endif
    int N = (int)mesh.nodes.size();

    std::vector<double> lam(N), qv(N);
    for (int i = 0; i < N; ++i) {
        lam[i] = lambda_(mesh.nodes[i].x, mesh.nodes[i].y);
        qv[i] = q_src(mesh.nodes[i].x, mesh.nodes[i].y);
    }
    DokMatrix K(N); std::vector<double> F(N, 0.0);
    assemble(mesh, lam, qv, K, F);

    auto bcs = boundary_table(mesh, a, b, c, d);
    apply_neumann(mesh, bcs, F, q_neu);

    std::vector<std::pair<int, double>> dbc;
    const double tol = 1e-10;
    (void)tol;
    for (auto& bc : bcs)
        if (bc.type == 1) dbc.push_back({ bc.node, u_exact(bc.x,bc.y) });
    apply_dirichlet(K, F, dbc);

    auto u = pcg(to_csr(K), F);

    double err = 0;
    for (int i = 0; i < N; ++i)
        err = std::max(err, std::abs(u[i] - u_exact(mesh.nodes[i].x, mesh.nodes[i].y)));

    std::cout << std::scientific << std::setprecision(2);
    std::cout << "test [" << TASK_NAME << "]:  max error = " << err
        << (err < TEST_TOL ? "  ok\n" : "  FAIL\n");
    return err < TEST_TOL;
}

int main()
{

    std::cout << "Running " << TASK_NAME << "\n";
    if (!run_test()) return 1;

    double a = 0, b = 1, c = 0, d = 1;
    int nx = 40, ny = 40;
    std::cout << "mesh " << nx << "x" << ny << "\n";

    Mesh mesh = make_mesh(a, b, c, d, nx, ny);
    auto u = solve(mesh, a, b, c, d);

    double umin = *std::min_element(u.begin(), u.end());
    double umax = *std::max_element(u.begin(), u.end());
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "u: [" << umin << ", " << umax << "]\n";

    write_csv("lab1_result.csv", mesh, u);
    write_vtk("lab1_result.vtk", mesh, u);
    write_txt("lab1_result.txt", mesh, u, a, b, c, d, nx, ny);
    std::cout << "written: lab1_result.csv, lab1_result.vtk, lab1_result.txt\n";
    return 0;
}