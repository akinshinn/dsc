/*
 * ЛР4 — Поиск начального узла
 *
 * Программа находит псевдопериферийный узел графа сетки
 * (алгоритм Джорджа—Лю). Псевдопериферийный узел используется
 * как начальный для алгоритма Катхилла-МакКи.
 * Для найденного узла строится и выводится корневая
 * структура уровней.
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

// ---- Корневая структура уровней (из ЛР3) ----

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
}

// ---- Поиск начального (псевдопериферийного) узла ----
// Алгоритм Джорджа—Лю:
// 1. Выбрать произвольный узел r, построить структуру уровней
// 2. На последнем уровне найти узел v с минимальной степенью
// 3. Построить структуру уровней от v
// 4. Если глубина увеличилась — повторить с шага 2, иначе v — результат

int find_starting_node(const StaticAdjacency& sa) {
    int n = (int)sa.xadj.size() - 1;

    // Степень узла
    auto degree = [&](int v) {
        return sa.xadj[v + 1] - sa.xadj[v];
    };

    // Начинаем с узла 0
    int r = 0;
    LevelStructure ls = build_level_structure(r, sa);
    int depth = (int)ls.levels.size();

    std::cout << "  Итерация: корень=" << r
              << ", глубина=" << depth << "\n";

    while (true) {
        // На последнем уровне ищем узел с минимальной степенью
        auto& last_level = ls.levels.back();
        int best = last_level[0];
        int best_deg = degree(best);
        for (int nd : last_level) {
            int d = degree(nd);
            if (d < best_deg) {
                best = nd;
                best_deg = d;
            }
        }

        // Строим структуру уровней от найденного узла
        LevelStructure ls2 = build_level_structure(best, sa);
        int new_depth = (int)ls2.levels.size();

        std::cout << "  Итерация: корень=" << best
                  << ", глубина=" << new_depth << "\n";

        if (new_depth <= depth) {
            // Глубина не увеличилась — нашли псевдопериферийный узел
            return best;
        }

        // Глубина увеличилась — продолжаем
        r = best;
        ls = ls2;
        depth = new_depth;
    }
}

// ---- Главная функция ----

int main() {
    std::cout << "=== ЛР4: Поиск начального узла ===\n\n";

    // Генерация сетки
    int nx = 15, ny = 10;
    Mesh mesh = generate_quad_mesh(nx, ny);
    write_mesh("mesh.txt", mesh);
    std::cout << "Сетка: " << mesh.nodes.size() << " узлов, "
              << mesh.elems.size() << " элементов\n\n";

    // Структура смежности
    auto dyn_adj = build_dynamic_adjacency(mesh);
    StaticAdjacency sa = build_static_adjacency(dyn_adj);

    // Поиск начального узла
    std::cout << "Поиск псевдопериферийного узла (алгоритм Джорджа—Лю):\n";
    int start = find_starting_node(sa);

    std::cout << "\nНайденный начальный узел: " << start
              << " (координаты: " << std::fixed << std::setprecision(2)
              << mesh.nodes[start].x << ", "
              << mesh.nodes[start].y << ")\n";

    // Построение и вывод структуры уровней для найденного узла
    LevelStructure ls = build_level_structure(start, sa);
    std::cout << "Глубина структуры уровней: " << ls.levels.size() << "\n";

    write_level_structure("levels_start.txt", ls);
    std::cout << "Структура уровней записана в levels_start.txt\n";

    // Для сравнения — структура уровней от узла 0
    LevelStructure ls0 = build_level_structure(0, sa);
    write_level_structure("levels_node0.txt", ls0);
    std::cout << "Для сравнения: структура от узла 0 — глубина "
              << ls0.levels.size() << " (записана в levels_node0.txt)\n";

    std::cout << "\nГотово.\n";
    return 0;
}
