// Метод конечных элементов для уравнения теплопроводности
// div(lambda * grad(u)) + f = 0

#include <vector>
#include <array>
#include <map>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <functional>

using namespace std;

// ============== Структуры данных ==============

struct Node { double x, y; };
struct Triangle { array<int, 3> v; };
struct BoundaryCondition { int node, type; double x, y; };

struct Mesh {
    vector<Node> nodes;
    vector<Triangle> elements;
};

// ============== Генерация сетки ==============

Mesh createMesh(double x0, double x1, double y0, double y1, int nx, int ny) {
    Mesh mesh;
    
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++)
            mesh.nodes.push_back({
                x0 + i * (x1 - x0) / nx,
                y0 + j * (y1 - y0) / ny
            });
    
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            int n0 = j * (nx + 1) + i;
            int n1 = n0 + 1;
            int n2 = n0 + (nx + 1);
            int n3 = n2 + 1;
            mesh.elements.push_back({{n0, n1, n3}});
            mesh.elements.push_back({{n0, n3, n2}});
        }
    
    return mesh;
}

// ============== Разреженная матрица ==============

class SparseMatrix {
public:
    int n;
    map<pair<int,int>, double> data;
    
    SparseMatrix(int size) : n(size) {}
    void add(int i, int j, double val) { data[{i, j}] += val; }
    void set(int i, int j, double val) { data[{i, j}] = val; }
};

// ============== Решатель СЛАУ ==============

double dot(const vector<double>& a, const vector<double>& b) {
    double s = 0;
    for (size_t i = 0; i < a.size(); i++) s += a[i] * b[i];
    return s;
}

vector<double> multiply(const SparseMatrix& A, const vector<double>& x) {
    vector<double> y(A.n, 0.0);
    for (auto& [idx, val] : A.data)
        y[idx.first] += val * x[idx.second];
    return y;
}

vector<double> solveCG(const SparseMatrix& A, const vector<double>& b) {
    int n = A.n;
    vector<double> x(n, 0.0), r = b, p = r;
    double rsold = dot(r, r);
    
    for (int iter = 0; iter < 100000; iter++) {
        vector<double> Ap = multiply(A, p);
        double alpha = rsold / dot(p, Ap);
        
        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        
        double rsnew = dot(r, r);
        if (sqrt(rsnew) < 1e-10) break;
        
        for (int i = 0; i < n; i++)
            p[i] = r[i] + (rsnew / rsold) * p[i];
        
        rsold = rsnew;
    }
    return x;
}

// ============== Локальная матрица жесткости ==============

void localStiffness(const array<Node, 3>& tri, double lambda, double f,
                    double K[3][3], double F[3]) {
    double x1 = tri[0].x, y1 = tri[0].y;
    double x2 = tri[1].x, y2 = tri[1].y;
    double x3 = tri[2].x, y3 = tri[2].y;
    
    double D = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    double area = abs(D) / 2.0;
    
    double b[3] = {(y2 - y3) / D, (y3 - y1) / D, (y1 - y2) / D};
    double c[3] = {(x3 - x2) / D, (x1 - x3) / D, (x2 - x1) / D};
    
    for (int i = 0; i < 3; i++) {
        F[i] = f * area / 3.0;
        for (int j = 0; j < 3; j++)
            K[i][j] = lambda * area * (b[i] * b[j] + c[i] * c[j]);
    }
}

// ============== Сборка глобальной системы ==============

void assemble(const Mesh& mesh, 
              function<double(double,double)> lambda,
              function<double(double,double)> f,
              SparseMatrix& K, vector<double>& F) {
    
    for (auto& elem : mesh.elements) {
        array<Node, 3> tri = {mesh.nodes[elem.v[0]], mesh.nodes[elem.v[1]], mesh.nodes[elem.v[2]]};
        
        double xc = (tri[0].x + tri[1].x + tri[2].x) / 3.0;
        double yc = (tri[0].y + tri[1].y + tri[2].y) / 3.0;
        
        double Ke[3][3], Fe[3];
        localStiffness(tri, lambda(xc, yc), f(xc, yc), Ke, Fe);
        
        for (int i = 0; i < 3; i++) {
            F[elem.v[i]] += Fe[i];
            for (int j = 0; j < 3; j++)
                K.add(elem.v[i], elem.v[j], Ke[i][j]);
        }
    }
}

// ============== Граничные условия ==============

void applyNeumann(const Mesh& mesh, const vector<BoundaryCondition>& bcs,
                  vector<double>& F, function<double(double,double)> q) {
    
    map<int, bool> neumannNodes;
    for (auto& bc : bcs) 
        if (bc.type == 2) neumannNodes[bc.node] = true;
    
    map<pair<int,int>, int> edgeCount;
    for (auto& elem : mesh.elements)
        for (int k = 0; k < 3; k++) {
            int a = elem.v[k], b = elem.v[(k + 1) % 3];
            if (a > b) swap(a, b);
            edgeCount[{a, b}]++;
        }
    
    for (auto& [edge, cnt] : edgeCount) {
        if (cnt != 1) continue;
        int i = edge.first, j = edge.second;
        if (!neumannNodes.count(i) || !neumannNodes.count(j)) continue;
        
        double xi = mesh.nodes[i].x, yi = mesh.nodes[i].y;
        double xj = mesh.nodes[j].x, yj = mesh.nodes[j].y;
        double L = sqrt((xj - xi) * (xj - xi) + (yj - yi) * (yj - yi));
        
        double qi = q(xi, yi), qj = q(xj, yj);
        F[i] -= L * (2 * qi + qj) / 6.0;
        F[j] -= L * (qi + 2 * qj) / 6.0;
    }
}

void applyDirichlet(SparseMatrix& K, vector<double>& F,
                    const vector<pair<int, double>>& dirichlet) {
    map<int, double> fixed(dirichlet.begin(), dirichlet.end());
    
    for (auto& [idx, val] : K.data) {
        int i = idx.first, j = idx.second;
        if (fixed.count(j) && !fixed.count(i))
            F[i] -= val * fixed[j];
    }
    
    vector<pair<int,int>> toRemove;
    for (auto& [idx, val] : K.data)
        if (fixed.count(idx.first) || fixed.count(idx.second))
            toRemove.push_back(idx);
    for (auto& idx : toRemove)
        K.data.erase(idx);
    
    for (auto& [node, val] : fixed) {
        K.set(node, node, 1.0);
        F[node] = val;
    }
}

// ============== Тестовые задачи ==============

struct Task1 {
    static double lambda(double, double) { return 1.0; }
    static double f(double, double) { return 0.0; }
    static double exact(double x, double y) { return x * x - y * y + 2.0; }
    static double neumann(double, double y) { return 2.0 * y; }
    
    static vector<BoundaryCondition> getBCs(const Mesh& m) {
        vector<BoundaryCondition> bcs;
        for (int i = 0; i < (int)m.nodes.size(); i++) {
            double x = m.nodes[i].x, y = m.nodes[i].y;
            if (abs(x) < 1e-10 || abs(x - 1) < 1e-10 || abs(y) < 1e-10)
                bcs.push_back({i, 1, x, y});
            else if (abs(y - 1) < 1e-10)
                bcs.push_back({i, 2, x, y});
        }
        return bcs;
    }
};

struct Task2 {
    static double lambda(double x, double y) { return x + y; }
    static double f(double x, double y) { 
        double s = x + y;
        return -10.0 * cos(5.0 * s) + 50.0 * s * sin(5.0 * s);
    }
    static double exact(double x, double y) { return sin(5.0 * (x + y)) + 2.0; }
    static double neumann(double x, double y) { return -(x + y) * 5.0 * cos(5.0 * (x + y)); }
    
    static vector<BoundaryCondition> getBCs(const Mesh& m) {
        vector<BoundaryCondition> bcs;
        for (int i = 0; i < (int)m.nodes.size(); i++) {
            double x = m.nodes[i].x, y = m.nodes[i].y;
            if (abs(x) < 1e-10 || abs(y) < 1e-10)
                bcs.push_back({i, 1, x, y});
            else if (abs(x - 1) < 1e-10 || abs(y - 1) < 1e-10)
                bcs.push_back({i, 2, x, y});
        }
        return bcs;
    }
};

// ============== Решение задачи ==============

template<typename Task>
double solve(int nx, vector<double>* solution = nullptr, Mesh* outMesh = nullptr) {
    Mesh mesh = createMesh(0, 1, 0, 1, nx, nx);
    int N = mesh.nodes.size();
    
    SparseMatrix K(N);
    vector<double> F(N, 0.0);
    
    assemble(mesh, Task::lambda, Task::f, K, F);
    
    auto bcs = Task::getBCs(mesh);
    applyNeumann(mesh, bcs, F, Task::neumann);
    
    vector<pair<int, double>> dirichlet;
    for (auto& bc : bcs)
        if (bc.type == 1)
            dirichlet.push_back({bc.node, Task::exact(bc.x, bc.y)});
    applyDirichlet(K, F, dirichlet);
    
    vector<double> u = solveCG(K, F);
    
    if (solution) *solution = u;
    if (outMesh) *outMesh = mesh;
    
    double maxErr = 0;
    for (int i = 0; i < N; i++)
        maxErr = max(maxErr, abs(u[i] - Task::exact(mesh.nodes[i].x, mesh.nodes[i].y)));
    
    return maxErr;
}

// ============== Экспорт для Mathematica ==============

template<typename Task>
void exportMathematica(const string& filename, int nx) {
    vector<double> u;
    Mesh mesh;
    solve<Task>(nx, &u, &mesh);
    
    ofstream file(filename);
    file << fixed << setprecision(8);
    
    file << "(* FEM Solution Data *)\n\n";
    
    file << "numSolution = {";
    for (size_t i = 0; i < mesh.nodes.size(); i++) {
        if (i > 0) file << ", ";
        file << "{" << mesh.nodes[i].x << ", " << mesh.nodes[i].y << ", " << u[i] << "}";
    }
    file << "};\n\n";
    
    file << "exactSolution = {";
    for (size_t i = 0; i < mesh.nodes.size(); i++) {
        if (i > 0) file << ", ";
        file << "{" << mesh.nodes[i].x << ", " << mesh.nodes[i].y << ", " 
             << Task::exact(mesh.nodes[i].x, mesh.nodes[i].y) << "}";
    }
    file << "};\n\n";
    
    file << "errorData = {";
    for (size_t i = 0; i < mesh.nodes.size(); i++) {
        if (i > 0) file << ", ";
        double err = abs(u[i] - Task::exact(mesh.nodes[i].x, mesh.nodes[i].y));
        file << "{" << mesh.nodes[i].x << ", " << mesh.nodes[i].y << ", " << err << "}";
    }
    file << "};\n\n";
    
    // Код визуализации
    file << "(* Visualization *)\n";
    file << "pNum = ListPlot3D[numSolution, PlotLabel -> \"Numerical\", ColorFunction -> \"Rainbow\", Mesh -> None];\n";
    file << "pExact = ListPlot3D[exactSolution, PlotLabel -> \"Exact\", ColorFunction -> \"Rainbow\", Mesh -> None];\n";
    file << "pErr = ListPlot3D[errorData, PlotLabel -> \"Error\", ColorFunction -> \"TemperatureMap\", Mesh -> None];\n";
    file << "GraphicsGrid[{{pNum, pExact}, {pErr, SpanFromLeft}}, ImageSize -> 800]\n";
}

// ============== Main ==============

int main() {
    cout << "FEM 2D - Convergence Test\n\n";
    
    vector<int> sizes = {5, 10, 20, 40};
    
    cout << "Task 1: u = x^2 - y^2 + 2, lambda = 1\n";
    cout << setw(8) << "N" << setw(12) << "h" << setw(14) << "Error" << setw(10) << "Order\n";
    cout << string(44, '-') << "\n";
    
    double prevErr = 0, prevH = 0;
    for (int n : sizes) {
        double h = 1.0 / n;
        double err = solve<Task1>(n);
        double order = (prevErr > 0) ? log(prevErr / err) / log(prevH / h) : 0;
        
        cout << setw(8) << n 
             << setw(12) << fixed << setprecision(4) << h
             << setw(14) << scientific << setprecision(2) << err;
        if (prevErr > 0) cout << setw(10) << fixed << setprecision(2) << order;
        cout << "\n";
        
        prevErr = err; prevH = h;
    }
    
    cout << "\nTask 2: u = sin(5(x+y)) + 2, lambda = x + y\n";
    cout << setw(8) << "N" << setw(12) << "h" << setw(14) << "Error" << setw(10) << "Order\n";
    cout << string(44, '-') << "\n";
    
    prevErr = 0; prevH = 0;
    for (int n : sizes) {
        double h = 1.0 / n;
        double err = solve<Task2>(n);
        double order = (prevErr > 0) ? log(prevErr / err) / log(prevH / h) : 0;
        
        cout << setw(8) << n 
             << setw(12) << fixed << setprecision(4) << h
             << setw(14) << scientific << setprecision(2) << err;
        if (prevErr > 0) cout << setw(10) << fixed << setprecision(2) << order;
        cout << "\n";
        
        prevErr = err; prevH = h;
    }
    
    exportMathematica<Task1>("task1_data.m", 20);
    exportMathematica<Task2>("task2_data.m", 30);
    cout << "\nExported: task1_data.m, task2_data.m\n";
    
    return 0;
}
