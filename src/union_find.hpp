#pragma once

#include <vector>
#include <numeric>
using namespace std;

struct UnionFind {
	vector<int> e;
	UnionFind(int n) {
        e.assign(n, -1);
    }
    void reset() {
        e.assign(e.size(), -1);
    }
	int find(int x) {
        if (e[x] < 0)
            return x;
        return e[x] = find(e[x]);
    }
	void join(int a, int b) {
        a = find(a);
        b = find(b);
        // if (a > b) swap(a, b);
        if (a == b) return;
        e[b] = a;
    }
};


struct UnionFindWithUndo {
	vector<int> e;
    vector<pair<int,int>> mods;
	UnionFindWithUndo(int n) {
        e.resize(n);
        iota(e.begin(), e.end(), 0);
    }
	int find(int x) {
        if (e[x] == x)
            return x;
        return find(e[x]);
    }
	void joinBySize(int a, int b) {
        if (a > b) swap(a, b);
        mods.emplace_back(b, e[b]);
        e[b] = a;
    }
	void join(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return;
        if (a > b) swap(a, b);
        mods.emplace_back(b, e[b]);
        e[b] = a;
    }
    void undo(int num = 1) {
        while(num--) {
            auto [a, b] = mods.back();
            e[a] = b;
            mods.pop_back();
        }
    }
};

