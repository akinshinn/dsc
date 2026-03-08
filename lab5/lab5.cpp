/*
 * ЛР5 — Оптимизация нумерации узлов обратным методом Катхилла-МакКи (RCM)
 *
 * Программа выполняет перенумерацию узлов сетки методом RCM
 * для уменьшения полуширины ленты и оболочки матрицы.
 * Поддерживаются типы сеток из ЛР3 (Семестр 1):
 *   тип 1 — четырёхугольные элементы (4 узла)
 *   тип 2 — треугольные, диагональ /  (3 узла)
 *   тип 3 — треугольные, диагональ \  (3 узла)
 *   тип 4 — треугольные, 4 треугольника через центр (3 узла)
 * Сравниваются параметры до и после перенумерации.
 * Перенумерованная сетка сохраняется в файл.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <queue>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cmath>

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

// ---- Поиск начального узла (из ЛР4) ----

int find_starting_node(const StaticAdjacency& sa) {
    auto degree = [&](int v) {
        return sa.xadj[v + 1] - sa.xadj[v];
    };

    int r = 0;
    LevelStructure ls = build_level_structure(r, sa);
    int depth = (int)ls.levels.size();

    while (true) {
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

        LevelStructure ls2 = build_level_structure(best, sa);
        int new_depth = (int)ls2.levels.size();

        if (new_depth <= depth)
            return best;

        r = best;
        ls = ls2;
        depth = new_depth;
    }
}

// ---- Алгоритм Катхилла-МакКи (CM) ----

std::vector<int> cuthill_mckee(int start, const StaticAdjacency& sa) {
    int n = (int)sa.xadj.size() - 1;
    std::vector<int> order;
    order.reserve(n);
    std::vector<bool> visited(n, false);

    auto degree = [&](int v) {
        return sa.xadj[v + 1] - sa.xadj[v];
    };

    order.push_back(start);
    visited[start] = true;

    for (int head = 0; head < (int)order.size(); ++head) {
        int v = order[head];

        std::vector<int> neighbors;
        for (int k = sa.xadj[v]; k < sa.xadj[v + 1]; ++k) {
            int nb = sa.adjncy[k];
            if (!visited[nb])
                neighbors.push_back(nb);
        }

        std::sort(neighbors.begin(), neighbors.end(),
            [&](int a, int b) { return degree(a) < degree(b); });

        for (int nb : neighbors) {
            visited[nb] = true;
            order.push_back(nb);
        }
    }
    return order;
}

// ---- Обратный метод Катхилла-МакКи (RCM) ----

std::vector<int> reverse_cuthill_mckee(int start, const StaticAdjacency& sa) {
    std::vector<int> order = cuthill_mckee(start, sa);
    std::reverse(order.begin(), order.end());
    return order;
}

// ---- Применение перенумерации к сетке ----

Mesh renumber_mesh(const Mesh& mesh, const std::vector<int>& order) {
    int n = (int)mesh.nodes.size();

    std::vector<int> old_to_new(n);
    for (int new_idx = 0; new_idx < n; ++new_idx)
        old_to_new[order[new_idx]] = new_idx;

    Mesh m2;
    m2.nodes.resize(n);
    for (int i = 0; i < n; ++i)
        m2.nodes[old_to_new[i]] = mesh.nodes[i];

    m2.elems.resize(mesh.elems.size());
    for (int e = 0; e < (int)mesh.elems.size(); ++e) {
        int np = (int)mesh.elems[e].nodes.size();
        m2.elems[e].nodes.resize(np);
        for (int k = 0; k < np; ++k)
            m2.elems[e].nodes[k] = old_to_new[mesh.elems[e].nodes[k]];
    }

    return m2;
}

// ---- Подсчёт параметров матрицы ----

struct MatrixParams {
    int half_bandwidth;
    long long envelope;
};

MatrixParams compute_params(const StaticAdjacency& sa) {
    int n = (int)sa.xadj.size() - 1;
    MatrixParams p;
    p.half_bandwidth = 0;
    p.envelope = 0;

    for (int i = 0; i < n; ++i) {
        int min_j = i;
        for (int k = sa.xadj[i]; k < sa.xadj[i + 1]; ++k) {
            int j = sa.adjncy[k];
            int diff = std::abs(i - j);
            if (diff > p.half_bandwidth)
                p.half_bandwidth = diff;
            if (j < min_j)
                min_j = j;
        }
        p.envelope += (i - min_j);
    }
    return p;
}

// ---- Вывод параметров в файл ----

void write_params(const std::string& fname,
                  const MatrixParams& before, const MatrixParams& after,
                  int start_node) {
    std::ofstream f(fname);
    f << "=== Параметры матрицы до и после перенумерации RCM ===\n\n";
    f << "Начальный узел для RCM: " << start_node << "\n\n";

    f << "До перенумерации:\n";
    f << "  Полуширина ленты: " << before.half_bandwidth << "\n";
    f << "  Оболочка (профиль): " << before.envelope << "\n\n";

    f << "После перенумерации RCM:\n";
    f << "  Полуширина ленты: " << after.half_bandwidth << "\n";
    f << "  Оболочка (профиль): " << after.envelope << "\n\n";

    double bw_ratio = (before.half_bandwidth > 0)
        ? 100.0 * (before.half_bandwidth - after.half_bandwidth) / before.half_bandwidth
        : 0.0;
    double env_ratio = (before.envelope > 0)
        ? 100.0 * (before.envelope - after.envelope) / before.envelope
        : 0.0;

    f << "Уменьшение полуширины: " << std::fixed << std::setprecision(1)
      << bw_ratio << "%\n";
    f << "Уменьшение оболочки:   " << std::fixed << std::setprecision(1)
      << env_ratio << "%\n";
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
    std::cout << "=== ЛР5: Оптимизация нумерации узлов методом RCM ===\n\n";

    int nx = 15, ny = 10;
    int type = 2;

    Mesh mesh = generate_quad_mesh(nx, ny, type);
    write_mesh("mesh_original.txt", mesh);
    std::cout << "Сетка: " << mesh.nodes.size() << " узлов, "
              << mesh.elems.size() << " элементов"
              << ", тип " << type << " (" << mesh_type_name(type) << ")\n";

    // Структура смежности исходной сетки
    auto dyn_adj = build_dynamic_adjacency(mesh);
    StaticAdjacency sa = build_static_adjacency(dyn_adj);

    // Параметры матрицы до перенумерации
    MatrixParams params_before = compute_params(sa);
    std::cout << "\nДо перенумерации:\n";
    std::cout << "  Полуширина ленты: " << params_before.half_bandwidth << "\n";
    std::cout << "  Оболочка (профиль): " << params_before.envelope << "\n";

    // Поиск начального узла
    std::cout << "\nПоиск начального узла...\n";
    int start = find_starting_node(sa);
    std::cout << "Начальный узел: " << start << "\n";

    // Перенумерация RCM
    std::cout << "\nВыполнение RCM...\n";
    std::vector<int> rcm_order = reverse_cuthill_mckee(start, sa);

    // Применение перенумерации к сетке
    Mesh mesh_rcm = renumber_mesh(mesh, rcm_order);

    // Структура смежности перенумерованной сетки
    auto dyn_adj2 = build_dynamic_adjacency(mesh_rcm);
    StaticAdjacency sa2 = build_static_adjacency(dyn_adj2);

    // Параметры матрицы после перенумерации
    MatrixParams params_after = compute_params(sa2);
    std::cout << "\nПосле перенумерации RCM:\n";
    std::cout << "  Полуширина ленты: " << params_after.half_bandwidth << "\n";
    std::cout << "  Оболочка (профиль): " << params_after.envelope << "\n";

    // Вывод улучшения
    double bw_ratio = 100.0 * (params_before.half_bandwidth - params_after.half_bandwidth)
                     / params_before.half_bandwidth;
    double env_ratio = 100.0 * (params_before.envelope - params_after.envelope)
                      / params_before.envelope;
    std::cout << "\nУменьшение полуширины: " << std::fixed << std::setprecision(1)
              << bw_ratio << "%\n";
    std::cout << "Уменьшение оболочки:   " << env_ratio << "%\n";

    // Сохранение результатов
    write_mesh("mesh_rcm.txt", mesh_rcm);
    std::cout << "\nПеренумерованная сетка записана в mesh_rcm.txt\n";

    write_params("params.txt", params_before, params_after, start);
    std::cout << "Параметры записаны в params.txt\n";

    std::cout << "\nГотово.\n";
    return 0;
}
