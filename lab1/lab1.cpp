#include <vector>
#include <array>
#include <map>
#include <functional>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <string>

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
            m.elems.push_back({ {n0, n1, n3} });
            m.elems.push_back({ {n0, n3, n2} });
        }
    return m;
}

struct TaskDef {
    std::string name;
    std::string bc_left, bc_right, bc_bot, bc_top;
    std::function<double(double, double)> lambda;
    std::function<double(double, double)> q_src;
    std::function<double(double, double)> u_exact;
    std::function<double(double, double)> q_neu;
    std::function<std::vector<BC>(const Mesh&, double, double, double, double)> boundary_table;
    double test_tol;
    int test_nx, test_ny;
};

// Задача 1: -div(grad u) = 0,  u = x^2 - y^2 + 2
TaskDef make_task1()
{
    TaskDef t;
    t.name = "Задача 1: u = x^2 - y^2 + 2, lambda=1, q=0";
    t.bc_left = "x=0: Дирихле  u = 2 - y^2";
    t.bc_right = "x=1: Дирихле  u = 3 - y^2";
    t.bc_bot = "y=0: Дирихле  u = x^2 + 2";
    t.bc_top = "y=1: Нейман   -(du/dn) = 2";

    t.lambda = [](double, double) { return 1.0; };
    t.q_src = [](double, double) { return 0.0; };
    t.u_exact = [](double x, double y) { return x * x - y * y + 2.0; };
    // верхняя грань: нормаль (0,+1), du/dy = -2y => -(du/dn) = 2y, при y=1 => 2
    t.q_neu = [](double, double y) { return 2.0 * y; };

    t.boundary_table = [](const Mesh& m, double a, double b, double c, double d) {
        std::vector<BC> bcs;
        const double tol = 1e-10;
        for (int i = 0; i < (int)m.nodes.size(); ++i) {
            double x = m.nodes[i].x, y = m.nodes[i].y;
            if (std::abs(x - a) < tol) bcs.push_back({ i, 1, x, y });
            else if (std::abs(x - b) < tol) bcs.push_back({ i, 1, x, y });
            else if (std::abs(y - c) < tol) bcs.push_back({ i, 1, x, y });
            else if (std::abs(y - d) < tol) bcs.push_back({ i, 2, x, y });
        }
        return bcs;
        };

    t.test_tol = 5e-2; t.test_nx = 20; t.test_ny = 20;
    return t;
}

// Задача 2: -div((x+y) grad u) = q,  u = sin(5(x+y)) + 2
TaskDef make_task2()
{
    TaskDef t;
    t.name = "Задача 2: u = sin(5(x+y)) + 2, lambda=x+y";
    t.bc_left = "x=0: Дирихле  u = sin(5y) + 2";
    t.bc_right = "x=1: Нейман   -(lambda*du/dn) = -(x+y)*5*cos(5(x+y))";
    t.bc_bot = "y=0: Дирихле  u = sin(5x) + 2";
    t.bc_top = "y=1: Нейман   -(lambda*du/dn) = -(x+y)*5*cos(5(x+y))";

    t.lambda = [](double x, double y) { return x + y; };
    t.q_src = [](double x, double y) {
        double s = x + y;
        return -10.0 * std::cos(5.0 * s) + 50.0 * s * std::sin(5.0 * s);
        };
    t.u_exact = [](double x, double y) { return std::sin(5.0 * (x + y)) + 2.0; };
    t.q_neu = [](double x, double y) { return -(x + y) * 5.0 * std::cos(5.0 * (x + y)); };

    t.boundary_table = [](const Mesh& m, double a, double b, double c, double d) {
        std::vector<BC> bcs;
        const double tol = 1e-10;
        for (int i = 0; i < (int)m.nodes.size(); ++i) {
            double x = m.nodes[i].x, y = m.nodes[i].y;
            if (std::abs(x - a) < tol) bcs.push_back({ i, 1, x, y });
            else if (std::abs(x - b) < tol) bcs.push_back({ i, 2, x, y });
            else if (std::abs(y - c) < tol) bcs.push_back({ i, 1, x, y });
            else if (std::abs(y - d) < tol) bcs.push_back({ i, 2, x, y });
        }
        return bcs;
        };

    t.test_tol = 1e-1; t.test_nx = 20; t.test_ny = 20;
    return t;
}

// ---- матрица жёсткости (DOK) ----

struct DokMatrix {
    int n;
    std::map<std::pair<int, int>, double> data;
    explicit DokMatrix(int sz) : n(sz) {}
    void   add(int i, int j, double v) { data[{i, j}] += v; }
    double get(int i, int j) const {
        auto it = data.find({ i,j });
        return it != data.end() ? it->second : 0.0;
    }
};

// ---- метод Гаусса с выбором главного элемента ----

std::vector<double> gauss(std::vector<std::vector<double>> a)
{
    int n = (int)a.size();
    for (int col = 0; col < n; ++col) {
        int pivot = col;
        for (int row = col + 1; row < n; ++row)
            if (std::abs(a[row][col]) > std::abs(a[pivot][col])) pivot = row;
        if (std::abs(a[pivot][col]) < 1e-15)
            throw std::runtime_error("gauss: вырожденная матрица");
        if (pivot != col) std::swap(a[col], a[pivot]);

        double inv = 1.0 / a[col][col];
        for (int row = col + 1; row < n; ++row) {
            double f = a[row][col] * inv;
            if (std::abs(f) < 1e-15) { a[row][col] = 0.0; continue; }
            for (int k = col; k <= n; ++k)
                a[row][k] -= f * a[col][k];
        }
    }
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        double s = a[i][n];
        for (int j = i + 1; j < n; ++j) s -= a[i][j] * x[j];
        x[i] = s / a[i][i];
    }
    return x;
}

std::vector<double> solve_gauss(const DokMatrix& K, const std::vector<double>& F)
{
    int n = K.n;
    std::vector<std::vector<double>> a(n, std::vector<double>(n + 1, 0.0));
    for (auto& [key, v] : K.data) a[key.first][key.second] = v;
    for (int i = 0; i < n; ++i) a[i][n] = F[i];
    return gauss(std::move(a));
}

// ---- сборка СЛАУ ----

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

    double b[3] = { (y2 - y3) / D, (y3 - y1) / D, (y1 - y2) / D };
    double c[3] = { (x3 - x2) / D, (x1 - x3) / D, (x2 - x1) / D };
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
        std::array<double, 3> qe{ qv[el.n[0]],  qv[el.n[1]],  qv[el.n[2]] };
        auto ke = local_matrix(xy, la, qe);
        for (int a = 0; a < 3; ++a) {
            F[el.n[a]] += ke.F[a];
            for (int b = 0; b < 3; ++b) K.add(el.n[a], el.n[b], ke.K[a][b]);
        }
    }
}

// ---- граничные условия ----

// ГУ 1-го рода (Дирихле): метод штрафа
//   K[node][node] += alpha,  F[node] += alpha * u_bc
// Штраф 1e6: компромисс между точностью ГУ и обусловленностью матрицы.
// При больших штрафах (1e10+) метод Гаусса теряет точность на сетках от 10x10.
const double PENALTY = 1e6;

void BC1(int node, double u_bc, DokMatrix& K, std::vector<double>& F)
{
    K.add(node, node, PENALTY);
    F[node] += PENALTY * u_bc;
}

// ГУ 2-го рода (Нейман): интегрирование потока по граничным рёбрам
//   F[i] -= L*(2*qi + qj)/6,  F[j] -= L*(qi + 2*qj)/6
void BC2(const Mesh& mesh, const std::vector<BC>& bcs,
    std::vector<double>& F,
    const std::function<double(double, double)>& q_neu)
{
    std::map<int, bool> neu;
    for (auto& bc : bcs)
        if (bc.type == 2) neu[bc.node] = true;

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
        double qi = q_neu(xi, yi), qj = q_neu(xj, yj);
        F[i] -= L * (2.0 * qi + qj) / 6.0;
        F[j] -= L * (qi + 2.0 * qj) / 6.0;
    }
}

// ---- решение задачи ----

std::vector<double> solve_task(const TaskDef& task, const Mesh& mesh,
    double a, double b, double c, double d)
{
    int N = (int)mesh.nodes.size();
    std::vector<double> lam(N), qv(N);
    for (int i = 0; i < N; ++i) {
        lam[i] = task.lambda(mesh.nodes[i].x, mesh.nodes[i].y);
        qv[i] = task.q_src(mesh.nodes[i].x, mesh.nodes[i].y);
    }

    DokMatrix K(N);
    std::vector<double> F(N, 0.0);
    assemble(mesh, lam, qv, K, F);

    auto bcs = task.boundary_table(mesh, a, b, c, d);
    BC2(mesh, bcs, F, task.q_neu);
    for (auto& bc : bcs)
        if (bc.type == 1) BC1(bc.node, task.u_exact(bc.x, bc.y), K, F);

    return solve_gauss(K, F);
}

// ---- тест ----

bool run_test(const TaskDef& task, double a, double b, double c, double d)
{
    Mesh mesh = make_mesh(a, b, c, d, task.test_nx, task.test_ny);
    auto u = solve_task(task, mesh, a, b, c, d);

    double err = 0.0;
    for (int i = 0; i < (int)mesh.nodes.size(); ++i)
        err = std::max(err, std::abs(u[i] - task.u_exact(mesh.nodes[i].x, mesh.nodes[i].y)));

    std::cout << std::scientific << std::setprecision(2)
        << task.name << ": погрешность = " << err
        << (err < task.test_tol ? "  ok\n" : "  FAIL\n");
    return err < task.test_tol;
}

// ---- вывод ----

void write_txt(const std::string& fname, const TaskDef& task, const Mesh& mesh,
    const std::vector<double>& u,
    double a, double b, double c, double d, int nx, int ny)
{
    std::ofstream f(fname);
    if (!f) throw std::runtime_error("не удалось открыть " + fname);

    double umin = *std::min_element(u.begin(), u.end());
    double umax = *std::max_element(u.begin(), u.end());
    double uavg = 0.0;
    for (double v : u) uavg += v;
    uavg /= u.size();

    f << task.name << "\n"
        << "Область: [" << a << "," << b << "] x [" << c << "," << d << "]\n"
        << "Сетка: " << nx << "x" << ny
        << "  узлов: " << mesh.nodes.size()
        << "  элементов: " << mesh.elems.size() << "\n"
        << "ГУ:\n"
        << "  " << task.bc_left << "\n"
        << "  " << task.bc_right << "\n"
        << "  " << task.bc_bot << "\n"
        << "  " << task.bc_top << "\n"
        << std::fixed << std::setprecision(6)
        << "u_min=" << umin << "  u_max=" << umax << "  u_avg=" << uavg << "\n\n";

    f << std::setw(6) << "узел"
        << std::setw(14) << "x"
        << std::setw(14) << "y"
        << std::setw(16) << "u" << "\n";
    f << std::fixed << std::setprecision(8);
    for (int i = 0; i < (int)mesh.nodes.size(); ++i)
        f << std::setw(6) << i
        << std::setw(14) << mesh.nodes[i].x
        << std::setw(14) << mesh.nodes[i].y
        << std::setw(16) << u[i] << "\n";
}

// ---- main ----

int main()
{
    const double a = 0.0, b = 1.0, c = 0.0, d = 1.0;
    const int nx = 15, ny = 15;

    TaskDef tasks[2] = { make_task1(), make_task2() };

    for (auto& task : tasks)
        if (!run_test(task, a, b, c, d)) return 1;

    for (int t = 0; t < 2; ++t) {
        auto& task = tasks[t];
        Mesh mesh = make_mesh(a, b, c, d, nx, ny);
        auto u = solve_task(task, mesh, a, b, c, d);

        double umin = *std::min_element(u.begin(), u.end());
        double umax = *std::max_element(u.begin(), u.end());
        std::cout << std::fixed << std::setprecision(4)
            << task.name << "\n  u: [" << umin << ", " << umax << "]\n";

        std::string fname = "result_task" + std::to_string(t + 1) + ".txt";
        write_txt(fname, task, mesh, u, a, b, c, d, nx, ny);
        std::cout << "  -> " << fname << "\n";
    }
    return 0;
}