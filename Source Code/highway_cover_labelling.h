#ifndef HGHWAY_LABELING_H_
#define HGHWAY_LABELING_H_

#include <sys/time.h>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <fstream>

#include "two_layer_queue.h"

//
// NOTE: Currently only unweighted and undirected graphs are supported.
//

class HighwayLabelling {
 public:
  // Constructs an index from a graph, given as a list of edges.
  HighwayLabelling(std::string filename, int k);
  HighwayLabelling();
  ~HighwayLabelling();
  void deallocate(bool q);

  void ConstructHighwayLabelling(int i, int topk[]);
  void BuildIndex(int topk[]);

  // For reflecting graph changes into highway labelling
  void UpdateLabelling(std::string filename, int prune);

  void insert(int a, int b, uint8_t *P, int *visited_vertices, int prune);
  void insert_prune(int a, int b, uint8_t *P, int *visited_vertices);
  void remove(int a, int b, uint8_t *P, uint8_t *D, int *visited_vertices);

  void SelectLandmarks_HD(int topk[]);
  void RemoveLandmarks(int topk[]);
  long LabellingSize();

  uint8_t query(int r, int v);
  uint8_t min(uint8_t a, uint8_t b);

  uint8_t HL_UB(int s, int t);
  void QueryDistance(std::string pairs, std::string output);

  void LoadLabelling_Pruned(std::string filename);
  void LoadLabelling_Full(std::string filename, int topk[], int k);
  void StoreLabelling(std::string filename);

 private:
  int V;  // total number of vertices
  long E; // total number of edges
  int K; // total number of landmarks

  uint8_t **vertices;
  uint8_t **distances;
  uint8_t **highway;
  uint8_t *C;
  std::vector<std::vector<int> > adj;
  std::map<int, uint8_t> landmarks;

  double GetCurrentTimeSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
  }

  long GetCurrentTimeMicroSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * (uint64_t) 1e6 + tv.tv_usec;
  }

  long GetCurrentTimeMilliSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000LL + tv.tv_usec / 1000;
  }

  // Statistics
  double time_, time_querying_sec_;
  long time_querying_microsec_, time_querying_millisec_;
};

HighwayLabelling::HighwayLabelling() { }

HighwayLabelling::~HighwayLabelling() { }

void HighwayLabelling::deallocate(bool q) {

  for(int i = 0; i < V; i++)
    delete [] distances[i];
  delete [] distances;

  for(int i = 0; i < K; i++)
    delete [] highway[i];
  delete [] highway;

  if(q) {
    for(int i = 0; i < V; i++)
      delete [] vertices[i];
    delete [] vertices;
    delete [] C;
  }
}

HighwayLabelling::HighwayLabelling(std::string filename, int k) {
  K = k; V = 0; E = 0;

  std::ifstream ifs(filename);
  if (ifs.is_open()) {
    ifs >> V >> E;

    adj.reserve(V);

    int v, w, deg;
    while (ifs >> v >> deg) {
      adj[v].reserve(deg);
      for (int i = 0; i < deg; i++) {
        ifs >> w;
        adj[v].push_back(w);
      }
    }
    ifs.close();

    std::cout << "V : " << V << " E : " << E << std::endl << std::endl;
  } else
      std::cout << "Unable to open file" << std::endl;
}

void HighwayLabelling::RemoveLandmarks(int topk[]) {

  for(int i = 0; i < K; i++) {
    for (int v : adj[topk[i]]) {
      adj[v].erase(std::remove(adj[v].begin(), adj[v].end(), topk[i]), adj[v].end());
      adj[v].shrink_to_fit();
    }
    adj[topk[i]].clear();
    adj[topk[i]].shrink_to_fit();
  }
}

long HighwayLabelling::LabellingSize() {
  long size = 0;
  for (int i = 0; i < V; i++) {
    for (int j = 0; j < K; j++) {
      if(distances[i][j] != 111)
        size++;
    }
  }

  return size; //(V + 2 * size) / (1024 * 1024);
}

void HighwayLabelling::ConstructHighwayLabelling(int i, int topk[]) {

  uint8_t *P = new uint8_t[V];
  for(int j = 0; j < V; j++)
    P[j] = 111;

  std::queue<int> que[2];
  que[0].push(topk[i]); que[0].push(-1);
  distances[topk[i]][i] = 0; P[topk[i]] = 0; int use = 0;
  while (!que[0].empty()) {
    int u = que[use].front();
    que[use].pop();

    if(u == -1) {
      use = 1 - use;
      que[use].push(-1);
      continue;
    }

    for (int w : adj[u]) {
      if (P[w] == 111) {
        P[w] = P[u] + 1;
        if(use == 1 || landmarks.count(w) > 0)
          que[1].push(w);
        else {
          que[0].push(w);
          distances[w][i] = P[w];
        }
      }
    }
  }

  for(int j = 0; j < K; j++) {
    highway[i][j] = P[topk[j]];
  }

  delete [] P;
}

void HighwayLabelling::BuildIndex(int topk[]) {

  for(int i = 0; i < K; i++)
    landmarks[topk[i]] = i;

  // Initialization
  distances = new uint8_t*[V];
  for(int i = 0; i < V; i++) {
    distances[i] = new uint8_t[K];
    for(int j = 0; j < K; j++)
      distances[i][j] = 111;
  }

  highway = new uint8_t*[K];
  for(int i = 0; i < K; i++) {
    highway[i] = new uint8_t[K];
    for(int j = 0; j < K; j++)
      highway[i][j] = 111;
  }

  // Start computing Highway Labelling (HL)
  time_ = -GetCurrentTimeSec();
  for (int i = 0; i < K; i++)
    ConstructHighwayLabelling(i, topk);
  time_ += GetCurrentTimeSec();
  std::cout << "Construction Time (sec.): " << time_ << " Labelling Size (MB): " << LabellingSize() << std::endl;
}

void HighwayLabelling::UpdateLabelling(std::string filename, int prune) {

  uint8_t *P = new uint8_t[V];
  uint8_t *D = new uint8_t[V];
  int *visited_vertices = new int[V];
  for(int j = 0; j < V; j++) {
    P[j] = 111; D[j] = 111;
  }

  std::ifstream ifs(filename); int a, b; std::string op;
  time_ = -GetCurrentTimeSec();
  while (ifs >> op >> a >> b) {
    if(op == "EI") {
      adj[a].push_back(b);
      adj[b].push_back(a);

      insert(a, b, P, visited_vertices, prune);
    } else if(op == "ED") {
      adj[a].erase(std::remove(adj[a].begin(), adj[a].end(), b), adj[a].end());
      adj[b].erase(std::remove(adj[b].begin(), adj[b].end(), a), adj[b].end());

      remove(a, b, P, D, visited_vertices);
    }
  }
  time_ += GetCurrentTimeSec();
  std::cout << "Update Time (sec.): " << time_ << " Updated Labelling Size (MB): " <<  LabellingSize() << std::endl;

  ifs.close();
  delete [] P;
  delete [] D;
  delete [] visited_vertices;
}

void HighwayLabelling::insert(int a, int b, uint8_t *P, int *visited_vertices, int prune) {

  int c;
  for(int i = 0; i < K; i++) {
    c = 0;
    std::queue<int> que;

    uint8_t da = query(i, a);
    uint8_t db = query(i, b);

    if(da > db) {
      P[a] = db + 1; que.push(a); visited_vertices[c] = a; c++;
    } else if(da < db) {
      P[b] = da + 1; que.push(b); visited_vertices[c] = b; c++;
    } else
      continue;

    while (!que.empty()) {
      int u = que.front();
      que.pop();

      if(landmarks.count(u) > 0) {
        highway[i][landmarks[u]] = P[u];
	if(prune == 1) continue;
        for (int w : adj[u]) {
          if(P[w] == 111) {
            P[w] = query(i, w);
	    visited_vertices[c] = w; c++;
            if(P[w] >= P[u] + 1) {
              P[w] = P[u] + 1;
              que.push(w);
            }
          }
        }
      } else {
        bool flag = false;  std::vector<int> aff_neighs;
        for (int w : adj[u]) {
          if(P[w] == 111) {
            P[w] = query(i, w);
	    visited_vertices[c] = w; c++;
            if(P[w] == P[u] - 1) {
              if(!flag && distances[w][i] == 111)
                flag = true;
            } else if(P[w] >= P[u] + 1) {
	      aff_neighs.push_back(w);
            }
          } else {
            if(P[w] == P[u] - 1) {
              if(!flag && distances[w][i] == 111)
                flag = true;
            }
          }
        }

        if(flag) {
          distances[u][i] = 111;
	  if(prune == 1) continue;
	} else {
          distances[u][i] = P[u];
	}

	for(int w : aff_neighs) {
          P[w] = P[u] + 1;
          que.push(w);
        }
      }
    }

    for(int j = 0; j < c; j++)
      P[visited_vertices[j]] = 111;
  }
}

void HighwayLabelling::remove(int a, int b, uint8_t *P, uint8_t *D, int *visited_vertices) {

  int dLsToEdge[K][2];
  for(int j = 0; j < K; j++) {
    dLsToEdge[j][0] = query(j, a);
    dLsToEdge[j][1] = query(j, b);
  }

  int c;
  for (int i = 0; i < K; i++) {
    c = 0;
    std::queue<int> que[2];
    std::unordered_set<int> anchors;

    if(dLsToEdge[i][0] > dLsToEdge[i][1]) {
      P[a] = dLsToEdge[i][0];
      que[0].push(a); visited_vertices[c] = a; c++;
    } else if(dLsToEdge[i][0] < dLsToEdge[i][1]) {
      P[b] = dLsToEdge[i][1];
      que[0].push(b); visited_vertices[c] = b; c++;
    } else
      continue;

    que[0].push(-1);
    int use = 0, flag = 0;
    while (!que[flag].empty() || !anchors.empty()) {
      int u = que[use].front();
      que[use].pop();

      if(u == -1) {
        if(!anchors.empty() && use == 0) {
          std::vector<std::pair<int, int> > anchorVs;

          std::unordered_set<int>::iterator itr;
          for (itr = anchors.begin(); itr != anchors.end(); ++itr)
            anchorVs.push_back(std::make_pair(D[*itr], *itr));
          std::sort(anchorVs.begin(), anchorVs.end());

          for(int j = 0; j < anchorVs.size(); j++)
            que[1].push(anchorVs[j].second);

 	        flag = 1;
          anchors.clear();
        }
        use = 1 - use;
        que[use].push(-1);
        continue;
      }

      if(use == 0) {
        // Check for case
        bool case1 = false, case2 = false, check = false; std::vector<int> vcase1, vcase2, vcase3;
        for (int w : adj[u]) {
          if(P[w] == 111) {

            int neigh_d = query(i, w);
            if(neigh_d == P[u] - 1) {
              case1 = true; case2 = false; vcase1.push_back(w);
            } else if(neigh_d == P[u] && !case1) {
              case2 = true; vcase2.push_back(w);
            } else if(neigh_d >= P[u] + 1) {
	            vcase3.push_back(w);
            }
          }
        }

        if(case1) {
          if(D[u] >= P[u]) {
            if(landmarks.count(u) > 0) {
              highway[i][landmarks[u]] = P[u];  highway[landmarks[u]][i] = P[u];
            } else {
	            for(int j = 0; j < vcase1.size(); j++) {
                if(distances[vcase1[j]][i] == 111) {
                  check = true; break;
                }
              }

              if(D[u] > P[u]) {
                if(!check)
                  distances[u][i] = P[u];
                else
                  distances[u][i] = 111;
              } else {
                if(distances[u][i] != 111) {
                  if(check)
                    distances[u][i] = 111;
                }
              }
            }

            D[u] = P[u];
            anchors.insert(u);
          }

          if(landmarks.count(u) == 0 && !check) {
            for(int j = 0; j < vcase3.size(); j++) {
              P[vcase3[j]] = P[u] + 1;
              que[use].push(vcase3[j]);
	            visited_vertices[c] = vcase3[j]; c++;
            }
          }
        } else if(case2) {
          if(D[u] >= P[u] + 1) {
            if(landmarks.count(u) > 0) {
              highway[i][landmarks[u]] = P[u] + 1;  highway[landmarks[u]][i] = P[u] + 1;
            } else {
	            for(int j = 0; j < vcase2.size(); j++) {
                if(distances[vcase2[j]][i] == 111) {
                  check = true; break;
                }
              }

              if(D[u] > P[u] + 1) {
                if(!check)
                  distances[u][i] = P[u] + 1;
                else
                  distances[u][i] = 111;
              } else {
                if(distances[u][i] != 111) {
                  if(check)
                    distances[u][i] = 111;
                }
              }
            }

            D[u] = P[u] + 1;
            anchors.insert(u);
          }

          if(landmarks.count(u) == 0) {
            for(int j = 0; j < vcase3.size(); j++) {
              P[vcase3[j]] = P[u] + 1;
              que[use].push(vcase3[j]);
	            visited_vertices[c] = vcase3[j]; c++;
            }
          }
        } else {
          for(int j = 0; j < vcase3.size(); j++) {
            P[vcase3[j]] = P[u] + 1;
            que[use].push(vcase3[j]);
	          visited_vertices[c] = vcase3[j]; c++;
          }
	      }
      } else if(use == 1) {
        for (int w : adj[u]) {
          if(P[w] != 111) {
            if(D[w] >= D[u] + 1) {
              if(landmarks.count(w) > 0) {
                highway[i][landmarks[w]] = D[u] + 1; highway[landmarks[w]][i] = D[u] + 1;
              } else {
		            if(D[w] > D[u] + 1) {
                  if(distances[u][i] != 111)
                    distances[w][i] = D[u] + 1;
                  else
                    distances[w][i] = 111;
                } else {
                  if(distances[w][i] != 111) {
                    if(distances[u][i] == 111)
                      distances[w][i] = 111;
                  }
                }
              }

              D[w] = D[u] + 1;
              anchors.insert(w);
            }
          }
        }
      }
    }

    if(flag == 0) {
      for(int j = 0; j < c; j++)
        distances[visited_vertices[j]][i] = 111;
    }

    for(int j = 0; j < c; j++) {
      P[visited_vertices[j]] = 111; D[visited_vertices[j]] = 111;
    }
  }
}

void HighwayLabelling::SelectLandmarks_HD(int topk[]) {
  std::vector<std::pair<int, int> > deg(V);
  long sum = 0;
  for (int v = 0; v < V; v++) {
    deg[v] = std::make_pair(adj[v].size(), v);
    sum = sum + adj[v].size();
  }
  std::sort(deg.rbegin(), deg.rend());

  for (int v = 0; v < K; v++)
    topk[v] = deg[v].second;
}

uint8_t HighwayLabelling::min(uint8_t a, uint8_t b) {
  return (a < b) ? a : b;
}

uint8_t HighwayLabelling::query(int r, int v) {

  uint8_t m = 111;
  for(int i = 0; i < K; i++)
    m = min(m, distances[v][i] + highway[r][i]);
  return m;
}

uint8_t HighwayLabelling::HL_UB(int s, int t) {

  uint8_t m = 111; int i, j;
  for(i = 0; i < C[s]; i++) {
    for (j = 0; j < C[t]; j++)
      m = min(m, distances[s][i] + highway[vertices[s][i]][vertices[t][j]] + distances[t][j]);
  }
  return m;
}

void HighwayLabelling::QueryDistance(std::string pairs, std::string output) {
  std::vector<TwoLayerQueue > qque;
  std::vector<uint8_t> qdist[2];

  qdist[0].resize(V, 111);
  qdist[1].resize(V, 111);
  qque.push_back(TwoLayerQueue(V));
  qque.push_back(TwoLayerQueue(V));

  time_querying_millisec_ = 0;
  std::ifstream ifs(pairs); int s = 0, t = 0, total = 0;
  std::ofstream ofs(std::string(output) + std::to_string(K));
  while(ifs >> s >> t) { total++;

    double temp = -GetCurrentTimeMilliSec();

    uint8_t dist_upper = HL_UB(s, t);
    uint8_t res = dist_upper, dis[2] = {0, 0};
    for (int dir = 0; dir < 2; dir++){
      int v = dir == 0 ? s : t;
      qque[dir].clear();
      qque[dir].push(v);
      qque[dir].next();
      qdist[dir][v] = 0;
    }

    while (!qque[0].empty() && !qque[1].empty()) {
      int use = (qque[0].size() <= qque[1].size()) ? 0 : 1;
      dis[use]++;

      if (dis[0] + dis[1] == dist_upper) {
        res = dis[0] + dis[1];
        goto LOOP_END;
      }

      while (!qque[use].empty()) {
        int v = qque[use].front();
        qque[use].pop();

        for (int w : adj[v]) {
          uint8_t &src_d = qdist[    use][w];
          uint8_t &dst_d = qdist[1 - use][w];
          if (src_d != 111) continue;
          if (dst_d != 111) {
            res = qdist[use][v] + 1 + dst_d;
            goto LOOP_END;
          } else {
            qque[use].push(w);
            qdist[use][w] = qdist[use][v] + 1;
          }
        }
      }
      qque[use].next();
    }
    LOOP_END:

    temp += GetCurrentTimeMilliSec();
    time_querying_millisec_ += temp;

    for (int dir = 0; dir < 2; dir++) {
      for (int v : qque[dir])
        qdist[dir][v] = 111;
      qque[dir].clear();
    }

    ofs << s << " " << t << " " << (int) min(res, dist_upper) << "\n";
  }
  std::cout << "Average Query Time (ms) : " << (double) time_querying_millisec_ / total << std::endl;
  ifs.close();
  ofs.close();
}

void HighwayLabelling::LoadLabelling_Pruned(std::string filename) {
  std::ifstream ifs(std::string(filename) + std::to_string(K) + std::string("_index"));

  C = new uint8_t[V];
  vertices = new uint8_t*[V];
  distances = new uint8_t*[V];

  for(int i = 0; i < V; i++) {
    ifs.read((char*)&C[i], sizeof(C[i]));
    vertices[i] = new uint8_t[C[i]];
    distances[i] = new uint8_t[C[i]];
    for(uint8_t j = 0; j < C[i]; j++) {
      ifs.read((char*)&vertices[i][j], sizeof(vertices[i][j]));
      ifs.read((char*)&distances[i][j], sizeof(distances[i][j]));
    }
  }
  ifs.close();

  ifs.open(std::string(filename) + std::to_string(K) + std::string("_highway"));
  highway = new uint8_t*[K];
  for(uint8_t i = 0; i < K; i++) {
    highway[i] = new uint8_t[K];
    for(uint8_t j = 0; j < K; j++)
      ifs.read((char*)&highway[i][j], sizeof(highway[i][j]));
  }
  ifs.close();
}

void HighwayLabelling::LoadLabelling_Full(std::string filename, int topk[], int k) {
  K = k;
  for(int i = 0; i < K; i++)
    landmarks[topk[i]] = i;

  std::ifstream ifs(std::string(filename) + std::to_string(K) + std::string("_index"));
  distances = new uint8_t*[V];
  for(int i = 0; i < V; i++) {
    distances[i] = new uint8_t[K];
    for(int j = 0; j < K; j++)
      distances[i][j] = 111;
  }

  uint8_t C, idx;
  for(int i = 0; i < V; i++) {
    ifs.read((char*)&C, sizeof(C));
    for(uint8_t j = 0; j < C; j++) {
      ifs.read((char*)&idx, sizeof(idx));
      ifs.read((char*)&distances[i][idx], sizeof(distances[i][idx]));
    }
  }
  ifs.close();

  ifs.open(std::string(filename) + std::to_string(K) + std::string("_highway"));
  highway = new uint8_t*[K];
  for(uint8_t i = 0; i < K; i++) {
    highway[i] = new uint8_t[K];
    for(uint8_t j = 0; j < K; j++)
      ifs.read((char*)&highway[i][j], sizeof(highway[i][j]));
  }
  ifs.close();
}

void HighwayLabelling::StoreLabelling(std::string filename) {
  std::ofstream ofs(std::string(filename) + std::to_string(K) + std::string("_index"));

  for(int i = 0; i < V; i++) {
    uint8_t C = 0;
    for(int j = 0; j < K; j++) {
      if(distances[i][j] != 111)
        C++;
    }

    ofs.write((char*)&C, sizeof(C));
    for(uint8_t j = 0; j < K; j++) {
      if(distances[i][j] != 111) {
        ofs.write((char*)&j, sizeof(j));
        ofs.write((char*)&distances[i][j], sizeof(distances[i][j]));
      }
    }
  }
  ofs.close();

  ofs.open(std::string(filename) + std::to_string(K) + std::string("_highway"));
  for(int i = 0; i < K; i++) {
    for(int j = 0; j < K; j++) {
      ofs.write((char*)&highway[i][j], sizeof(highway[i][j]));
    }
  }
  ofs.close();
}

#endif  // PRUNED_LANDMARK_LABELING_H_
