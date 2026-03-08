/*
 * ЛР3 — Построение корневой структуры уровней
 *
 * Программа считывает сетку из файла, строит структуру смежности,
 * затем выполняет обход в ширину (BFS) от заданного узла
 * и формирует корневую структуру уровней.
 * Результат выводится в файл.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <queue>
#include <string>
#include <iomanip>

// ---- Структуры данных (из ЛР2) ----

struct Node { double x, y; };
struct Triangle { int n[3]; };

struct Mesh {
    std::vector<Node> nodes;
    std::vector<Triangle> elems;
};

struct StaticAdjacency {
    std::vector<int> xadj;
    std::vector<int> adjncy;
};

// ---- Генерация сетки (из ЛР2) ----

Mesh generate_quad_mesh(int nx, int ny) {
    Mesh m;
    double x0 = 0, y0 = 0;
    double x1 = 10, y1 = 0;
    double x2 = 8, y2 = 10;
    double x3 = 0, y3 = 7;

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

// ---- Запись/чтение сетки (из ЛР2) ----

void write_mesh(const std::string& fname, const Mesh& m) {
    std::ofstream f(fname);
    f << m.nodes.size() << "\n";
    f << std::fixed << std::setprecision(6);
    for (auto& nd : m.nodes)
        f << nd.x << " " << nd.y << "\n";
    f << m.elems.size() << "\n";
    for (auto& el : m.elems)
        f << el.n[0] << " " << el.n[1] << " " << el.n[2] << "\n";
}

Mesh read_mesh(const std::string& fname) {
    std::ifstream f(fname);
    if (!f) {
        std::cerr << "Ошибка: не удалось открыть " << fname << "\n";
        exit(1);
    }
    Mesh m;
    int nn; f >> nn;
    m.nodes.resize(nn);
    for (int i = 0; i < nn; ++i)
        f >> m.nodes[i].x >> m.nodes[i].y;
    int ne; f >> ne;
    m.elems.resize(ne);
    for (int i = 0; i < ne; ++i)
        f >> m.elems[i].n[0] >> m.elems[i].n[1] >> m.elems[i].n[2];
    return m;
}

// ---- Построение структуры смежности (из ЛР2) ----

std::vector<std::set<int>> build_dynamic_adjacency(const Mesh& m) {
    int n = (int)m.nodes.size();
    std::vector<std::set<int>> adj(n);
    for (auto& el : m.elems)
        for (int i = 0; i < 3; ++i)
            for (int j = i + 1; j < 3; ++j) {
                adj[el.n[i]].insert(el.n[j]);
                adj[el.n[j]].insert(el.n[i]);
            }
    return adj;
}

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

// ---- Корневая структура уровней (BFS) ----

struct LevelStructure {
    int root;
    std::vector<std::vector<int>> levels; // levels[k] — узлы уровня k
    std::vector<int> level_of;            // level_of[i] — номер уровня узла i
};

LevelStructure build_level_structure(int root, const StaticAdjacency& sa) {
    int n = (int)sa.xadj.size() - 1;
    LevelStructure ls;
    ls.root = root;
    ls.level_of.assign(n, -1);

    std::queue<int> q;
    q.push(root);
    ls.level_of[root] = 0;

    while (!q.empty()) {
        int v = q.front(); q.pop();
        int lv = ls.level_of[v];

        // Добавляем новый уровень при необходимости
        if (lv >= (int)ls.levels.size())
            ls.levels.resize(lv + 1);
        ls.levels[lv].push_back(v);

        // Обход соседей
        for (int k = sa.xadj[v]; k < sa.xadj[v + 1]; ++k) {
            int nb = sa.adjncy[k];
            if (ls.level_of[nb] == -1) {
                ls.level_of[nb] = lv + 1;
                q.push(nb);
            }
        }
    }
    return ls;
}

// ---- Вывод корневой структуры уровней в файл ----

void write_level_structure(const std::string& fname, const LevelStructure& ls) {
    std::ofstream f(fname);
    f << "Корневая структура уровней\n";
    f << "Корневой узел: " << ls.root << "\n";
    f << "Количество уровней: " << ls.levels.size() << "\n\n";

    for (int lv = 0; lv < (int)ls.levels.size(); ++lv) {
        f << "Уровень " << std::setw(3) << lv
          << " (" << std::setw(4) << ls.levels[lv].size() << " узлов): ";
        for (int nd : ls.levels[lv])
            f << nd << " ";
        f << "\n";
    }
    std::cout << "Структура уровней записана в " << fname << "\n";
}

// ---- Главная функция ----

int main() {
    std::cout << "=== ЛР3: Построение корневой структуры уровней ===\n\n";

    // Генерация и сохранение сетки
    int nx = 15, ny = 10;
    Mesh mesh = generate_quad_mesh(nx, ny);
    write_mesh("mesh.txt", mesh);
    std::cout << "Сетка: " << mesh.nodes.size() << " узлов, "
              << mesh.elems.size() << " элементов\n";

    // Построение структуры смежности
    auto dyn_adj = build_dynamic_adjacency(mesh);
    StaticAdjacency sa = build_static_adjacency(dyn_adj);
    std::cout << "Структура смежности построена\n";

    // Построение корневой структуры уровней от узла 0 (нумерация с 0)
    int root = 0;
    std::cout << "\nПостроение структуры уровней от узла " << root << "...\n";
    LevelStructure ls = build_level_structure(root, sa);

    std::cout << "Количество уровней: " << ls.levels.size() << "\n";
    std::cout << "Узлов на первом уровне: " << ls.levels[0].size() << "\n";
    std::cout << "Узлов на последнем уровне: "
              << ls.levels.back().size() << "\n";

    // Вывод в файл
    write_level_structure("levels.txt", ls);

    // Вывод первых нескольких уровней на экран
    std::cout << "\nПервые уровни:\n";
    int show_levels = std::min(5, (int)ls.levels.size());
    for (int lv = 0; lv < show_levels; ++lv) {
        std::cout << "  Уровень " << lv << " (" << ls.levels[lv].size()
                  << " узлов): ";
        int show_nodes = std::min(10, (int)ls.levels[lv].size());
        for (int k = 0; k < show_nodes; ++k)
            std::cout << ls.levels[lv][k] << " ";
        if ((int)ls.levels[lv].size() > show_nodes) std::cout << "...";
        std::cout << "\n";
    }

    std::cout << "\nГотово.\n";
    return 0;
}
