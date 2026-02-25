#include <vector>
#include <array>
#include <map>
#include <functional>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

struct Node { double x, y; };
struct Element { std::array<int, 3> n; };
struct BC { int node, type; double x, y; };
struct Mesh { std::vector<Node> nodes; std::vector<Element> elems; };

Mesh make_mesh(double a, double b, double c, double d, int nx, int ny) {
    Mesh m; m.nodes.resize((nx + 1) * (ny + 1));
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i)
            m.nodes[j * (nx + 1) + i] = { a + i * (b - a) / nx, c + j * (d - c) / ny };
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            int n0 = j * (nx + 1) + i, n1 = n0 + 1, n2 = n0 + (nx + 1), n3 = n2 + 1;
            m.elems.push_back({ {n0,n1,n3} }); m.elems.push_back({ {n0,n3,n2} });
        }
    return m;
}

struct DokMatrix {
    int n; std::map<std::pair<int, int>, double> data;
    explicit DokMatrix(int sz) : n(sz) {}
    void add(int i, int j, double v) { data[{i, j}] += v; }
    void set(int i, int j, double v) { data[{i, j}] = v; }
};

struct CsrMatrix {
    int n = 0; std::vector<double> val; std::vector<int> col, row;
    std::vector<double> mv(const std::vector<double>& x) const {
        std::vector<double> y(n, 0.0);
        for (int i = 0; i < n; ++i) for (int k = row[i]; k < row[i + 1]; ++k) y[i] += val[k] * x[col[k]];
        return y;
    }
};

CsrMatrix to_csr(const DokMatrix& d) {
    CsrMatrix c; c.n = d.n; c.row.resize(d.n + 1, 0);
    for (auto& [k, v] : d.data) c.row[k.first + 1]++;
    for (int i = 0; i < d.n; ++i) c.row[i + 1] += c.row[i];
    int nnz = c.row[d.n]; c.val.resize(nnz); c.col.resize(nnz);
    std::vector<int> pos(c.row.begin(), c.row.end());
    for (auto& [k, v] : d.data) { int p = pos[k.first]++; c.val[p] = v; c.col[p] = k.second; }
    return c;
}

static double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0; for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i]; return s;
}

std::vector<double> pcg(const CsrMatrix& A, const std::vector<double>& b, double tol = 1e-10, int maxit = 50000) {
    int n = A.n; std::vector<double> Minv(n, 1.0);
    for (int i = 0; i < n; ++i) for (int k = A.row[i]; k < A.row[i + 1]; ++k)
        if (A.col[k] == i && std::abs(A.val[k]) > 1e-300) Minv[i] = 1.0 / A.val[k];
    std::vector<double> x(n, 0.0), r = b, z(n), p(n);
    for (int i = 0; i < n; ++i) z[i] = Minv[i] * r[i]; p = z;
    double rz = dot(r, z), bn = std::sqrt(dot(b, b)); if (bn < 1e-300) return x;
    for (int it = 0; it < maxit; ++it) {
        auto Ap = A.mv(p); double pAp = dot(p, Ap); if (std::abs(pAp) < 1e-300) break;
        double alpha = rz / pAp;
        for (int i = 0; i < n; ++i) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }
        if (std::sqrt(dot(r, r)) < tol * bn) break;
        for (int i = 0; i < n; ++i) z[i] = Minv[i] * r[i];
        double rz2 = dot(r, z), beta = rz2 / rz;
        for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i]; rz = rz2;
    }
    return x;
}

struct Ke { double K[3][3]; double F[3]; };

Ke local_matrix(const std::array<Node, 3>& p, const std::array<double, 3>& lam, const std::array<double, 3>& q) {
    Ke ke{}; double x1=p[0].x, y1=p[0].y, x2=p[1].x, y2=p[1].y, x3=p[2].x, y3=p[2].y;
    double D = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1), S = std::abs(D) / 2.0;
    if (S < 1e-15) return ke;
    double b[3] = { (y2-y3)/D, (y3-y1)/D, (y1-y2)/D }, c[3] = { (x3-x2)/D, (x1-x3)/D, (x2-x1)/D };
    double la = (lam[0]+lam[1]+lam[2])/3.0, qa = (q[0]+q[1]+q[2])/3.0;
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) ke.K[i][j] = la * S * (b[i]*b[j] + c[i]*c[j]);
    for (int i = 0; i < 3; ++i) ke.F[i] = qa * S / 3.0;
    return ke;
}

void assemble(const Mesh& mesh, const std::vector<double>& lam, const std::vector<double>& qv, DokMatrix& K, std::vector<double>& F) {
    for (auto& el : mesh.elems) {
        std::array<Node, 3> xy{ mesh.nodes[el.n[0]], mesh.nodes[el.n[1]], mesh.nodes[el.n[2]] };
        std::array<double, 3> la{ lam[el.n[0]], lam[el.n[1]], lam[el.n[2]] };
        std::array<double, 3> qe{ qv[el.n[0]], qv[el.n[1]], qv[el.n[2]] };
        auto ke = local_matrix(xy, la, qe);
        for (int a = 0; a < 3; ++a) { F[el.n[a]] += ke.F[a]; for (int b = 0; b < 3; ++b) K.add(el.n[a], el.n[b], ke.K[a][b]); }
    }
}

void apply_dirichlet(DokMatrix& K, std::vector<double>& F, const std::vector<std::pair<int, double>>& dbc) {
    std::map<int, double> m(dbc.begin(), dbc.end());
    for (auto& [key, v] : K.data) { int i = key.first, j = key.second; if (m.count(j) && !m.count(i)) F[i] -= v * m[j]; }
    std::vector<std::pair<int, int>> rm;
    for (auto& [key, v] : K.data) if (m.count(key.first) || m.count(key.second)) rm.push_back(key);
    for (auto& k : rm) K.data.erase(k);
    for (auto& [idx, u] : m) { K.set(idx, idx, 1.0); F[idx] = u; }
}

void apply_neumann(const Mesh& mesh, const std::vector<BC>& bcs, std::vector<double>& F, std::function<double(double,double)> q0) {
    std::map<int, bool> neu; for (auto& bc : bcs) if (bc.type == 2) neu[bc.node] = true;
    std::map<std::pair<int, int>, int> ecnt;
    for (auto& el : mesh.elems) for (int k = 0; k < 3; ++k) {
        int a = el.n[k], b = el.n[(k + 1) % 3]; if (a > b) std::swap(a, b); ecnt[{a, b}]++;
    }
    for (auto& [e, cnt] : ecnt) {
        if (cnt != 1) continue; int i = e.first, j = e.second;
        if (!neu.count(i) || !neu.count(j)) continue;
        double xi=mesh.nodes[i].x, yi=mesh.nodes[i].y, xj=mesh.nodes[j].x, yj=mesh.nodes[j].y;
        double L = std::sqrt((xj-xi)*(xj-xi) + (yj-yi)*(yj-yi)), qi = q0(xi, yi), qj = q0(xj, yj);
        F[i] -= L * (2*qi + qj) / 6.0; F[j] -= L * (qi + 2*qj) / 6.0;
    }
}

double lambda1(double, double) { return 1.0; }
double q_src1(double, double) { return 0.0; }
double u_exact1(double x, double y) { return x*x - y*y + 2.0; }
double q_neu1(double x, double y) { (void)x; return 2.0*y; }

std::vector<BC> bt1(const Mesh& m, double a, double b, double c, double d) {
    std::vector<BC> bcs; double tol = 1e-10;
    for (int i = 0; i < (int)m.nodes.size(); ++i) {
        double x = m.nodes[i].x, y = m.nodes[i].y;
        if (std::abs(x-a) < tol) bcs.push_back({i,1,x,y});
        else if (std::abs(x-b) < tol) bcs.push_back({i,1,x,y});
        else if (std::abs(y-c) < tol) bcs.push_back({i,1,x,y});
        else if (std::abs(y-d) < tol) bcs.push_back({i,2,x,y});
    }
    return bcs;
}

double lambda2(double x, double y) { return x + y; }
double q_src2(double x, double y) { double s = x+y; return -10.0*std::cos(5.0*s) + 50.0*s*std::sin(5.0*s); }
double u_exact2(double x, double y) { return std::sin(5.0*(x+y)) + 2.0; }
double q_neu2(double x, double y) { return -(x+y)*5.0*std::cos(5.0*(x+y)); }

std::vector<BC> bt2(const Mesh& m, double a, double b, double c, double d) {
    std::vector<BC> bcs; double tol = 1e-10;
    for (int i = 0; i < (int)m.nodes.size(); ++i) {
        double x = m.nodes[i].x, y = m.nodes[i].y;
        if (std::abs(x-a) < tol) bcs.push_back({i,1,x,y});
        else if (std::abs(x-b) < tol) bcs.push_back({i,2,x,y});
        else if (std::abs(y-c) < tol) bcs.push_back({i,1,x,y});
        else if (std::abs(y-d) < tol) bcs.push_back({i,2,x,y});
    }
    return bcs;
}

double test(int task, int nx) {
    double a=0, b=1, c=0, d=1; Mesh mesh = make_mesh(a, b, c, d, nx, nx); int N = (int)mesh.nodes.size();
    std::vector<double> lam(N), qv(N);
    auto lf = (task==1) ? lambda1 : lambda2; auto qf = (task==1) ? q_src1 : q_src2;
    auto uf = (task==1) ? u_exact1 : u_exact2; auto nf = (task==1) ? q_neu1 : q_neu2;
    for (int i = 0; i < N; ++i) { lam[i] = lf(mesh.nodes[i].x, mesh.nodes[i].y); qv[i] = qf(mesh.nodes[i].x, mesh.nodes[i].y); }
    DokMatrix K(N); std::vector<double> F(N, 0.0); assemble(mesh, lam, qv, K, F);
    auto bcs = (task==1) ? bt1(mesh,a,b,c,d) : bt2(mesh,a,b,c,d); apply_neumann(mesh, bcs, F, nf);
    std::vector<std::pair<int, double>> dbc;
    for (auto& bc : bcs) if (bc.type == 1) dbc.push_back({bc.node, uf(bc.x, bc.y)});
    apply_dirichlet(K, F, dbc); auto u = pcg(to_csr(K), F);
    double err = 0; for (int i = 0; i < N; ++i) err = std::max(err, std::abs(u[i] - uf(mesh.nodes[i].x, mesh.nodes[i].y)));
    return err;
}

int main() {
    std::vector<int> sizes = {5, 10, 20, 40};
    std::cout << "=== Task 1: u = x^2 - y^2 + 2 ===\n";
    std::cout << std::setw(6) << "N" << std::setw(12) << "h" << std::setw(14) << "Error" << std::setw(10) << "Rate\n";
    double pe1=0, ph1=0;
    for (int n : sizes) {
        double h = 1.0/n, err = test(1, n), rate = (pe1>0) ? std::log(pe1/err)/std::log(ph1/h) : 0;
        std::cout << std::setw(6) << n << std::setw(12) << std::fixed << std::setprecision(4) << h 
                  << std::setw(14) << std::scientific << std::setprecision(2) << err;
        if (pe1>0) std::cout << std::setw(10) << std::fixed << std::setprecision(2) << rate;
        std::cout << "\n"; pe1=err; ph1=h;
    }
    std::cout << "\n=== Task 2: u = sin(5(x+y)) + 2 ===\n";
    std::cout << std::setw(6) << "N" << std::setw(12) << "h" << std::setw(14) << "Error" << std::setw(10) << "Rate\n";
    double pe2=0, ph2=0;
    for (int n : sizes) {
        double h = 1.0/n, err = test(2, n), rate = (pe2>0) ? std::log(pe2/err)/std::log(ph2/h) : 0;
        std::cout << std::setw(6) << n << std::setw(12) << std::fixed << std::setprecision(4) << h 
                  << std::setw(14) << std::scientific << std::setprecision(2) << err;
        if (pe2>0) std::cout << std::setw(10) << std::fixed << std::setprecision(2) << rate;
        std::cout << "\n"; pe2=err; ph2=h;
    }
    return 0;
}
