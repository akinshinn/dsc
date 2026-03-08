/*
 * ЛР3 — Построение корневой структуры уровней
 *
 * Программа считывает сетку из файла, строит структуру смежности,
 * затем выполняет обход в ширину (BFS) от заданного узла
 * и формирует корневую структуру уровней.
 * Поддерживаются типы сеток из ЛР3 (Семестр 1):
 *   тип 1 — четырёхугольные элементы (4 узла)
 *   тип 2 — треугольные, диагональ /  (3 узла)
 *   тип 3 — треугольные, диагональ \  (3 узла)
 *   тип 4 — треугольные, 4 треугольника через центр (3 узла)
 * Результат выводится в файл.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <queue>
#include <string>
#include <iomanip>
#include <algorithm>

// ---- Структуры данных (из ЛР2) ----

struct Node { double x, y; };

struct Element {
    std::vector<int> nodes;
};

struct Mesh {
    std::vector<Node>    nodes;
    std::vector<Element> elems;
};

struct StaticAdjacency {
    std::vector<int> xadj;
    std::vector<int> adjncy;
};

// ---- Генерация сетки (из ЛР2) ----

Mesh generate_quad_mesh(int nx, int ny, int type = 2) {
    Mesh m;
    double x0 = 0, y0 = 0;
    double x1 = 10, y1 = 0;
    double x2 = 8, y2 = 10;
    double x3 = 0, y3 = 7;

    int base_nodes = (nx + 1) * (ny + 1);
    m.nodes.resize(base_nodes);
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

    if (type == 4) {
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                int n0 = j * (nx + 1) + i, n1 = n0 + 1;
                int n2 = n0 + (nx + 1), n3 = n2 + 1;
                double cx = (m.nodes[n0].x + m.nodes[n1].x +
                             m.nodes[n2].x + m.nodes[n3].x) / 4.0;
                double cy = (m.nodes[n0].y + m.nodes[n1].y +
                             m.nodes[n2].y + m.nodes[n3].y) / 4.0;
                m.nodes.push_back({cx, cy});
            }
    }

    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            int n0 = j * (nx + 1) + i, n1 = n0 + 1;
            int n2 = n0 + (nx + 1), n3 = n2 + 1;
            switch (type) {
            case 1: m.elems.push_back({{n0, n1, n3, n2}}); break;
            case 2:
                m.elems.push_back({{n0, n1, n3}});
                m.elems.push_back({{n0, n3, n2}});
                break;
            case 3:
                m.elems.push_back({{n0, n1, n2}});
                m.elems.push_back({{n1, n3, n2}});
                break;
            case 4: {
                int nc = base_nodes + j * nx + i;
                m.elems.push_back({{n0, n1, nc}});
                m.elems.push_back({{n1, n3, nc}});
                m.elems.push_back({{n3, n2, nc}});
                m.elems.push_back({{n2, n0, nc}});
                break;
            }
            }
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
    for (auto& el : m.elems) {
        f << el.nodes.size();
        for (int nd : el.nodes) f << " " << nd;
        f << "\n";
    }
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
    for (int i = 0; i < ne; ++i) {
        int np; f >> np;
        m.elems[i].nodes.resize(np);
        for (int j = 0; j < np; ++j)
            f >> m.elems[i].nodes[j];
    }
    return m;
}

// ---- Построение структуры смежности (из ЛР2) ----

std::vector<std::set<int>> build_dynamic_adjacency(const Mesh& m) {
    int n = (int)m.nodes.size();
    std::vector<std::set<int>> adj(n);
    for (auto& el : m.elems) {
        int np = (int)el.nodes.size();
        for (int i = 0; i < np; ++i)
            for (int j = i + 1; j < np; ++j) {
                adj[el.nodes[i]].insert(el.nodes[j]);
                adj[el.nodes[j]].insert(el.nodes[i]);
            }
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
    std::vector<std::vector<int>> levels;
    std::vector<int> level_of;
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

        if (lv >= (int)ls.levels.size())
            ls.levels.resize(lv + 1);
        ls.levels[lv].push_back(v);

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

// ---- Название типа сетки ----

const char* mesh_type_name(int type) {
    switch (type) {
    case 1: return "четырёхугольники";
    case 2: return "треугольники (диагональ /)";
    case 3: return "треугольники (диагональ \\)";
    case 4: return "треугольники (4 через центр)";
    default: return "неизвестный";
    }
}

// ---- Главная функция ----

int main() {
    std::cout << "=== ЛР3: Построение корневой структуры уровней ===\n\n";

    int nx = 15, ny = 10;
    int type = 2;

    Mesh mesh = generate_quad_mesh(nx, ny, type);
    write_mesh("mesh.txt", mesh);
    std::cout << "Сетка: " << mesh.nodes.size() << " узлов, "
              << mesh.elems.size() << " элементов"
              << ", тип " << type << " (" << mesh_type_name(type) << ")\n";

    // Построение структуры смежности
    auto dyn_adj = build_dynamic_adjacency(mesh);
    StaticAdjacency sa = build_static_adjacency(dyn_adj);
    std::cout << "Структура смежности построена\n";

    // Построение корневой структуры уровней от узла 0
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
