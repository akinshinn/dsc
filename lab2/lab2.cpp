/*
 * ЛР2 — Построение структуры смежности
 *
 * Программа генерирует треугольную сетку на четырёхугольнике
 * с вершинами (0,0), (10,0), (8,10), (0,7), разбиение 15x10.
 * Сетка сохраняется в файл, считывается обратно,
 * строится динамическая и статическая структуры смежности,
 * результат выводится в файл.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <iomanip>
#include <cmath>

// ---- Структуры данных ----

struct Node { double x, y; };
struct Triangle { int n[3]; };

struct Mesh {
    std::vector<Node> nodes;
    std::vector<Triangle> elems;
};

// ---- Генерация сетки на четырёхугольнике ----
// Билинейное отображение (s,t) in [0,1]^2 -> (x,y)
// P0=(0,0), P1=(10,0), P2=(8,10), P3=(0,7)

Mesh generate_quad_mesh(int nx, int ny) {
    Mesh m;
    double x0 = 0, y0 = 0;   // P0
    double x1 = 10, y1 = 0;  // P1
    double x2 = 8, y2 = 10;  // P2
    double x3 = 0, y3 = 7;   // P3

    // Узлы через билинейное отображение
    m.nodes.resize((nx + 1) * (ny + 1));
    for (int j = 0; j <= ny; ++j) {
        double t = (double)j / ny;
        for (int i = 0; i <= nx; ++i) {
            double s = (double)i / nx;
            double x = (1 - s) * (1 - t) * x0 + s * (1 - t) * x1
                     + s * t * x2 + (1 - s) * t * x3;
            double y = (1 - s) * (1 - t) * y0 + s * (1 - t) * y1
                     + s * t * y2 + (1 - s) * t * y3;
            m.nodes[j * (nx + 1) + i] = {x, y};
        }
    }

    // Каждый прямоугольник делится на 2 треугольника
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            int n0 = j * (nx + 1) + i;
            int n1 = n0 + 1;
            int n2 = n0 + (nx + 1);
            int n3 = n2 + 1;
            m.elems.push_back({{n0, n1, n3}});
            m.elems.push_back({{n0, n3, n2}});
        }

    return m;
}

// ---- Запись сетки в файл ----

void write_mesh(const std::string& fname, const Mesh& m) {
    std::ofstream f(fname);
    f << m.nodes.size() << "\n";
    f << std::fixed << std::setprecision(6);
    for (auto& nd : m.nodes)
        f << nd.x << " " << nd.y << "\n";
    f << m.elems.size() << "\n";
    for (auto& el : m.elems)
        f << el.n[0] << " " << el.n[1] << " " << el.n[2] << "\n";
    std::cout << "Сетка записана в " << fname << "\n";
}

// ---- Чтение сетки из файла ----

Mesh read_mesh(const std::string& fname) {
    std::ifstream f(fname);
    if (!f) {
        std::cerr << "Ошибка: не удалось открыть " << fname << "\n";
        exit(1);
    }
    Mesh m;
    int nn;
    f >> nn;
    m.nodes.resize(nn);
    for (int i = 0; i < nn; ++i)
        f >> m.nodes[i].x >> m.nodes[i].y;
    int ne;
    f >> ne;
    m.elems.resize(ne);
    for (int i = 0; i < ne; ++i)
        f >> m.elems[i].n[0] >> m.elems[i].n[1] >> m.elems[i].n[2];
    std::cout << "Сетка считана из " << fname
              << ": " << nn << " узлов, " << ne << " элементов\n";
    return m;
}

// ---- Динамическая структура смежности ----
// Для каждого узла — множество (set) смежных узлов

std::vector<std::set<int>> build_dynamic_adjacency(const Mesh& m) {
    int n = (int)m.nodes.size();
    std::vector<std::set<int>> adj(n);

    for (auto& el : m.elems) {
        // Каждая пара вершин треугольника — смежные узлы
        for (int i = 0; i < 3; ++i)
            for (int j = i + 1; j < 3; ++j) {
                adj[el.n[i]].insert(el.n[j]);
                adj[el.n[j]].insert(el.n[i]);
            }
    }
    return adj;
}

// ---- Статическая структура смежности (формат CSR) ----

struct StaticAdjacency {
    std::vector<int> xadj;   // массив указателей (размер n+1)
    std::vector<int> adjncy; // массив смежных вершин
};

StaticAdjacency build_static_adjacency(const std::vector<std::set<int>>& dyn) {
    StaticAdjacency sa;
    int n = (int)dyn.size();
    sa.xadj.resize(n + 1);

    sa.xadj[0] = 0;
    for (int i = 0; i < n; ++i)
        sa.xadj[i + 1] = sa.xadj[i] + (int)dyn[i].size();

    sa.adjncy.resize(sa.xadj[n]);
    for (int i = 0; i < n; ++i) {
        int pos = sa.xadj[i];
        for (int nb : dyn[i])
            sa.adjncy[pos++] = nb;
    }
    return sa;
}

// ---- Вывод структуры смежности в файл ----

void write_adjacency(const std::string& fname, const StaticAdjacency& sa) {
    int n = (int)sa.xadj.size() - 1;
    std::ofstream f(fname);

    f << "Количество узлов: " << n << "\n";
    f << "Общее число связей: " << sa.adjncy.size() << "\n\n";
    f << "Структура смежности (узел: список смежных):\n";
    f << "---------------------------------------------\n";

    for (int i = 0; i < n; ++i) {
        f << std::setw(4) << i << ": ";
        for (int k = sa.xadj[i]; k < sa.xadj[i + 1]; ++k)
            f << sa.adjncy[k] << " ";
        f << "\n";
    }
    std::cout << "Структура смежности записана в " << fname << "\n";
}

// ---- Главная функция ----

int main() {
    std::cout << "=== ЛР2: Построение структуры смежности ===\n\n";

    // Генерация сетки
    int nx = 15, ny = 10;
    std::cout << "Генерация сетки " << nx << "x" << ny
              << " на четырёхугольнике (0,0)-(10,0)-(8,10)-(0,7)\n";
    Mesh mesh = generate_quad_mesh(nx, ny);

    // Запись сетки в файл
    write_mesh("mesh.txt", mesh);

    // Чтение сетки из файла (демонстрация считывания)
    Mesh mesh2 = read_mesh("mesh.txt");

    // Построение динамической структуры смежности
    std::cout << "\nПостроение динамической структуры смежности...\n";
    auto dyn_adj = build_dynamic_adjacency(mesh2);

    // Вывод статистики
    int max_deg = 0, min_deg = (int)mesh2.nodes.size();
    for (auto& s : dyn_adj) {
        int d = (int)s.size();
        max_deg = std::max(max_deg, d);
        min_deg = std::min(min_deg, d);
    }
    std::cout << "Степени узлов: min = " << min_deg
              << ", max = " << max_deg << "\n";

    // Построение статической структуры смежности
    std::cout << "Построение статической структуры смежности...\n";
    StaticAdjacency sa = build_static_adjacency(dyn_adj);

    // Вывод в файл
    write_adjacency("adjacency.txt", sa);

    // Вывод первых нескольких узлов для проверки
    std::cout << "\nПервые 5 узлов:\n";
    int show = std::min(5, (int)mesh2.nodes.size());
    for (int i = 0; i < show; ++i) {
        std::cout << "  Узел " << i << " ("
                  << std::fixed << std::setprecision(2)
                  << mesh2.nodes[i].x << ", " << mesh2.nodes[i].y
                  << "): смежные = {";
        for (int k = sa.xadj[i]; k < sa.xadj[i + 1]; ++k) {
            if (k > sa.xadj[i]) std::cout << ", ";
            std::cout << sa.adjncy[k];
        }
        std::cout << "}\n";
    }

    std::cout << "\nГотово.\n";
    return 0;
}
