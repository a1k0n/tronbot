#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <limits.h>
#include <map>
#include <vector>

#define TIMEOUT_USEC 850000
#define INITIAL_DEPTH 1
#define DRAW_PENALTY -100
#define VERBOSE 0

// {{{ position
struct position {
  static const char dx[5], dy[5];

  int x,y;
  position() {}
  position(int _x, int _y): x(_x), y(_y) {}
  // direction encoding:
  //  1
  // 4 2
  //  3
  position next(int move) { return position(x+dx[move], y+dy[move]); }
  position prev(int move) { return position(x-dx[move], y-dy[move]); }
};

bool operator==(const position &a, const position &b) { return a.x == b.x && a.y == b.y; }

const char position::dx[5]={0,0,1,0,-1};
const char position::dy[5]={0,-1,0,1,0};
// }}}

// {{{ Map
template <class T> struct Map {
  T *map;
  int width, height;
  Map() { map = NULL; }
  Map(int w, int h) { resize(w,h); }
  void resize(int w, int h) {
    width = w; height = h;
    map = new T[w*h];
    memset(map, 0, w*h*sizeof(T));
  }
  Map(const Map &m) { abort(); } // this shouldn't happen
  ~Map() { if(map) delete[] map; }
  T& operator()(position p) { return map[p.x + p.y*width]; }
  T& operator()(int x, int y) { return map[x + y*width]; }
  T& M(position p) { return map[p.x + p.y*width]; }
  T& M(int x, int y) { return map[x + y*width]; }

  position argmin(void) {
    position p(0,0);
    for(int j=0;j<height;j++) {
      for(int i=0;i<width;i++) {
        if(M(i,j) < M(p))
          p = position(i,j);
      }
    }
    return p;
  }

  void dump(void) {
    for(int j=0;j<height;j++) {
      for(int i=0;i<width;i++) {
        int n = map[i+j*width];
        if(n == 0) fprintf(stderr, "  ");
        else fprintf(stderr, "%2d", n);
      }
      fprintf(stderr, "\n");
    }
  }
};
// }}}

// {{{ gamestate
struct gamestate {
  position p[2]; // position in current state
  int m[2]; // last move made

  // affects map using m[0], m[1]
  gamestate move(Map<char> M) {
    gamestate s = *this;
    M(p[0]) = 1;
    M(p[1]) = 1;
    s.p[0] = p[0].next(m[0]);
    s.p[1] = p[1].next(m[1]);
    return s;
  }
  // undoes effect on map
  void unmove(Map<char> M) {
    M(p[0]) = 0;
    M(p[1]) = 0;
  }
};

// }}}

static Map<char> M;
static Map<int> dp0, dp1;
static gamestate curstate;

// {{{ imported map update garbage from original code
bool map_update()
{
  int x, y, c;
  int map_width, map_height;
  int num_items = fscanf(stdin, "%d %d\n", &map_width, &map_height);
  if (feof(stdin) || num_items < 2) {
    return false;
  }
  if(!M.map) {
    M.resize(map_width, map_height);
    dp0.resize(map_width, map_height);
    dp1.resize(map_width, map_height);
  }
  x = 0;
  y = 0;
  while (y < M.height && (c = fgetc(stdin)) != EOF) {
    switch (c) {
    case '\r':
      break;
    case '\n':
      if (x != M.width) {
	fprintf(stderr, "x != width in Board_ReadFromStream\n");
	return false;
      }
      ++y;
      x = 0;
      break;
    case '#':
      if (x >= M.width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return false;
      }
      M(x,y) = 1;
      ++x;
      break;
    case ' ':
      if (x >= M.width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return false;
      }
      M(x,y) = 0;
      ++x;
      break;
    case '1':
    case '2':
      if (x >= M.width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return false;
      }
      {
        position p(x,y);
        M(p) = 0;
        curstate.p[c - '1'] = p;
        curstate.m[c - '1'] = 0;
        ++x;
      }
      break;
    default:
      fprintf(stderr, "unexpected character %d in Board_ReadFromStream", c);
      return false;
    }
  }
  return true;
}
// }}}

// {{{ connected components algorithm

struct Components {
  Map<int> c;
  std::map<int,int> csize;

  Components(Map<char> &M, gamestate s): c(M.width, M.height)
  {
//    std::map<int,int> equiv;
    int nextclass = 1;
    M(s.p[0]) = 0;
    M(s.p[1]) = 0;
    for(int j=1;j<M.height-1;j++) {
      for(int i=1;i<M.width-1;i++) {
        if(M(i,j)) continue; // wall
        int cup   = c(i, j-1),
            cleft = c(i-1, j);
        if(cup == 0 && cleft == 0) { // new component
          c(i,j) = nextclass++;
        } else if(cup == cleft) { // existing component
          c(i,j) = cup;
        } else { // join components
          // deprecate the higher-numbered component in favor of the lower
          if(cleft == 0 || (cup != 0 && cup < cleft)) {
            c(i,j) = cup;
            if(cleft != 0) _merge(cleft, cup);
          } else {
            c(i,j) = cleft;
            if(cup != 0) _merge(cup, cleft);
          }
        }
      }
    }
    M(s.p[0]) = 1;
    M(s.p[1]) = 1;
#if 0
    dump();
    fprintf(stderr, "equivalences: ");
    for(std::map<int,int>::iterator k=equiv.begin();k!=equiv.end();k++)
      fprintf(stderr, "%d->%d ", k->first, k->second);
    fprintf(stderr, "\n");
#endif
    // now make another pass to compute connected area
    for(int j=1;j<M.height-1;j++) {
      for(int i=1;i<M.width-1;i++) {
        csize[c(i,j)] ++;
      }
    }
  }
  void dump() {
    std::map<int,int>::iterator i;
    for(i=csize.begin();i!=csize.end();i++) {
      fprintf(stderr, "area %d: %d vertices\n", i->first, i->second);
    }
    c.dump();
  }
  int component(const position &p) { return c(p); }
  int connectedarea(int component) { return csize[component]; }
  int connectedarea(const position &p) { return csize[c(p)]; }
private:
#if 0
  int _find_equiv(std::map<int,int> &equiv, int c) {
    while(true) {
      std::map<int,int>::iterator e = equiv.find(c);
      if(e == equiv.end()) break;
      if(c < e->second)
        c = e->second;
      else
        break;
    }
    return c;
  }
#endif
  void _merge(int o, int n) {
    for(int j=1;j<M.height-1;j++)
      for(int i=1;i<M.width-1;i++)
        if(c(i,j) == o)
          c(i,j) = n;
  }
};

// }}}

// {{{ basic geometric stuff
int runout(position p, int dir) {
  int r = 0;
  while(!M(p)) { r++; p = p.next(dir); }
  return r;
}

int degree(position x) {
  return 4 - M(x.next(1)) - M(x.next(2)) - M(x.next(3)) - M(x.next(4));
}

int turn(int d, int n) { return 1+(d-1+n)&3; }
// }}}

// {{{ run timing
long _get_time()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_usec + tv.tv_sec*1000000;
}

static long _timer;
static bool _timed_out = false;
void reset_timer(void) { _timer = _get_time(); _timed_out = false; }
long elapsed_time() { return _get_time() - _timer; }
bool timeout() { _timed_out = elapsed_time() > TIMEOUT_USEC; return _timed_out; }
// }}}

// {{{ Dijkstra's
void dijkstra(Map<int> &d, position s, Components &cp, int component)
{
  static std::vector<position> Q;
  int i,j;
  for(j=0;j<M.height;j++)
    for(i=0;i<M.width;i++) {
      if(cp.c(i,j) != component) continue;
      Q.push_back(position(i,j));
      d(i,j) = INT_MAX;
    }
  d(s) = 0;
  while(!Q.empty()) {
    position u(0,0);
    int min_d = INT_MAX, min_i=0;
    for(i=0;i<(int)Q.size();i++)
      if(d(Q[i]) < min_d) { u = Q[i]; min_d = d(u); min_i = i; }
    Q[min_i] = Q[Q.size()-1]; Q.pop_back();
    if(min_d == INT_MAX) continue; // ??
    for(int m=1;m<=4;m++) {
      position v = u.next(m);
      if(M(v)) continue;
      int alt = 1 + d(u);
      if(alt < d(v))
        d(v) = alt;
    }
  }
#if 0
  // REMOVEME: cleanup for printing really
  for(j=0;j<M.height;j++)
    for(i=0;i<M.width;i++) {
      if(d(i,j) == INT_MAX)
        d(i,j) = 0;
    }
#endif
}

// }}}

// {{{ alpha-beta iterative deepening search

static int evaluations=0;
int _evaluate_board(gamestate s, int player)
{
  Components cp(M, s);
  evaluations++;
#if VERBOSE >= 3
  fprintf(stderr, "evaluating board: \n");
  M(s.p[player]) = 2; M(s.p[player^1]) = 3; M.dump();
  M(s.p[0]) = 1; M(s.p[1]) = 1;
#endif
  int comp;
  if((comp = cp.component(s.p[0])) == cp.component(s.p[1])) {
    dijkstra(dp0, s.p[player], cp, comp);
    dijkstra(dp1, s.p[player^1], cp, comp);
#if VERBOSE >= 3
    Map<int> vor(M.width, M.height);
#endif
    int nodecount = 0;
    for(int j=0;j<M.height;j++)
      for(int i=0;i<M.width;i++) {
        int diff = dp0(i,j) - dp1(i,j);
        // if the opponent's distance is shorter than ours, then this is "their" node
        if(diff>0) { nodecount--; }
        // otherwise it's ours
        if(diff<0) { nodecount++; }
#if VERBOSE >= 3
        vor(i,j) = diff > 0 ? 1 : diff < 0 ? 2 : 0;
#endif
      }
#if VERBOSE >= 3
    dp0.dump();
    dp1.dump();
    vor.dump();
    fprintf(stderr, "player=%d nodecount: %d\n", player, nodecount);
#endif
    return nodecount;
  } else {
    int v = 100*(cp.connectedarea(s.p[player]) -
                 cp.connectedarea(s.p[player^1]));
#if VERBOSE >= 3
    fprintf(stderr, "player=%d connectedarea value: %d\n", player, v);
#endif
    return v;
  }
}

// do an iterative-deepening search on all moves and see if we can find a move
// sequence that cuts off our opponent
int _ab_runs=0;
int _alphabeta(int &move, gamestate s, int player, int a, int b, int itr)
{
  if(s.p[0] == s.p[1]) { return (player == 1 ? -1 : 1) * DRAW_PENALTY; } // crash!  draw!

  if(_timed_out || ((_ab_runs++)&127) == 0 && timeout()) {
#if VERBOSE >= 1
    fprintf(stderr, "timeout; a=%d b=%d itr=%d\n", a,b,itr);
#endif
    return a;
  }

  if(itr == 0) {
    int v = _evaluate_board(s, player);
#if VERBOSE >= 2
    fprintf(stderr, "_alphabeta(itr=%d [%d,%d,%d]|[%d,%d,%d] p=%d a=%d b=%d) -> %d\n", 
            itr, s.p[0].x, s.p[0].y, s.m[0], 
            s.p[1].x, s.p[1].y, s.m[1], player, a,b,v);
#endif
    return v;
  }
#if VERBOSE >= 2
  fprintf(stderr, "_alphabeta(itr=%d [%d,%d,%d]|[%d,%d,%d] p=%d a=%d b=%d)\n", 
          itr, s.p[0].x, s.p[0].y, s.m[0], 
          s.p[1].x, s.p[1].y, s.m[1], player, a,b);
#endif

  // periodically check timeout.  if we do time out, give up, we can't do any
  // more work; whatever we found so far will have to do
  for(int m=1;m<=4 && !_timed_out;m++) {
    if(M(s.p[player].next(m))) // impossible move?
      continue;
    gamestate r = s;
    r.m[player] = m;
    // after both players 0 and 1 make their moves, the game state updates
    if(player == 1) {
      r.p[0] = s.p[0].next(r.m[0]);
      r.p[1] = s.p[1].next(r.m[1]);
      M(r.p[0]) = 1;
      M(r.p[1]) = 1;
    }
    int m_; // next move; discard
    int a_ = -_alphabeta(m_, r, player^1, -b, -a, itr-1);
    if(a_ > a) {
      a = a_;
      move = m;
    }
    // undo game state update
    if(player == 1) {
      M(r.p[0]) = 0;
      M(r.p[1]) = 0;
      r.p[0] = s.p[0];
      r.p[1] = s.p[1];
    }

    if(a >= b) // beta cut-off
      break;
  }
  return a;
}

int next_move_alphabeta()
{
  int itr;
  int bestv = -1000000, bestm=1;
  M(curstate.p[0]) = 1;
  M(curstate.p[1]) = 1;
  reset_timer();
  evaluations=0;
  for(itr=INITIAL_DEPTH;itr<100 && !timeout();itr++) {
    int m;
    int bv=-1000000, bm=1;
    int v = _alphabeta(m, curstate, 0, -10000000, 10000000, itr*2);
#if 0
    if(v >= 5000) {
#if VERBOSE >= 1
      struct timeval tv;
      gettimeofday(&tv, NULL);
      fprintf(stderr, "%d.%06d: v=%d best=%d (m=%d) -> found compelling move\n", (int) tv.tv_sec, (int) tv.tv_usec, v, bestv, m);
#endif
      return m;
    }
#endif
    if(v > bv) { bv = v; bm = m;}
#if VERBOSE >= 1
    struct timeval tv;
    gettimeofday(&tv, NULL);
    //M.dump();
    fprintf(stderr, "%d.%06d: v=%d best=%d (m=%d) @depth %d _ab_runs=%d\n", (int) tv.tv_sec, (int) tv.tv_usec, v, bv, bm, itr*2, _ab_runs);
#endif
    if(v == -10000000) {
      // deeper searching is apparently impossible
      break;
    }
    // update bestv to the results of either the best-so-far this iteration or
    // the best complete iteration (in other words, better estimates in
    // previous iterations are irrelevant if a full exploration of a deeper
    // level doesn't reach them)
    if(!_timed_out || v > bestv) {
      bestv = bv; bestm = bm;
    }
  }
#if VERBOSE >= 1
  long e = elapsed_time();
  float rate = (float)evaluations*1000000.0/(float)e;
  fprintf(stderr, "%d evals in %ld us; %0.1f evals/sec\n", evaluations, e, rate);
  if(e > TIMEOUT_USEC*11/10) {
    fprintf(stderr, "10%% timeout violation: %ld us\n", e);
  }
#endif
  return bestm;
}
// }}}

// {{{ space-filling

// {{{ crappy greedy space filler that doesn't work as well as just wall-following
struct SpaceFiller
{
  int area, move;

  static SpaceFiller fill(position p, int max_area) {
    SpaceFiller f;
    f.area = _fill_area(f.move, p, max_area, 0);
    return f;
  }

#if 0
  static SpaceFiller greedy(position p) {
    SpaceFiller f;
    f.area = 0;
    for(int m=1;m<=4;m++) {
      int nm;
      int a = _fill_greedy(nm, p, 1);
      if(a > f.area) { f.area = a; f.move = nm; }
    return f;
  }
#endif

private:
#if 0
  int _fill_greedy(int &move, position p0, int majordir) {
    // try to fill in row-major or column-major order, looking ahead one row
    // (or column) for the largest opening

  }
  int _find_next_greedy(position p, int dir, int majordir) {
    //  0123456789ab
    //  >----------- ->    ##########...
    //---   ------       ..... 1#####
    // travel in dir, finding all 0->1 transitions in the next column over (by majordir)
    // for each 0->1 transition, find the run length going -dir
    // optimize for this_runlen + next_runlen
    int ow = M(p.next(majordir));
    int n = 0;
    int besta=0, bestn=0;
    while(!M(p.next(dir))) {
      p = p.next(dir);
      int nw = M(p.next(majordir));
      if(nw && nw != ow) {
        int a = n+runout(p.prev(dir).next(majordir), turn(dir, 2));
        turnarounds.push_back(n);
      }
      n++;
    }
  }
#endif

  // this is exponential and sucks
  static int _fill_area(int &move, position p0, int max_area, int lastmove) {
//    printf("_fill_area([%d,%d], area=%d)\n", p0.x, p0.y, max_area);
    if(max_area <= 1) { move = 0; return 1; }
    int besta = 0, bestm = 0;
    for(int m=1;m<=4;m++) {
      int r = 0;
      position p = p0;
      if(m == lastmove) continue;
      if(M(p.next(m))) continue;
      while(!M(p)) { r++; M(p) = 1; p = p.next(m); } // do!
      while(r>0) {
        // undo! (we just ran into a wall, so undo that first...)
        r--;
        p = p.prev(m);
        M(p) = 0;
        if(r>0) {
          // now try longest to shortest...
          int _m;
//          printf("%d: ", m);
          int a = r + _fill_area(_m, p, max_area - r, m);
          if(a >= max_area) { move = m; return a; }
          if(a > besta) { besta = a; bestm = m; }
        }
      }
    }
    move = bestm;
    return besta;
  }
};
// }}}

const int degreescore[] = {0, 4, 3, 2, 1};
int next_move_spacefill()
{
  M(curstate.p[0]) = 1;
  M(curstate.p[1]) = 1;

  gamestate nullstate;
  nullstate.p[0] = nullstate.p[1] = position(0,0);
  Components ca(M, nullstate);

  int bestm=1, bestv=0;
  for(int m=1;m<=4;m++) {
    position p = curstate.p[0].next(m);
    if(M(p)) continue;
    int v = ca.connectedarea(p) + degreescore[degree(p)];
    if(v > bestv) { bestv = v; bestm = m; }
#if VERBOSE >= 1
    fprintf(stderr, "move %d: ca=%d, degree=%d, v=%d\n", m, ca.connectedarea(p), degree(p), v);
#endif
  }


//  SpaceFiller s(curstate.p[0], c.connectedarea(curstate.p[0])/2);
//  return SpaceFiller::greedy(curstate.p[0]).move;
  return bestm;
}
// }}}

int next_move() {
  gamestate nullstate;
  nullstate.p[0] = nullstate.p[1] = position(0,0);
  Components cp(M, nullstate);
#if VERBOSE >= 3
  cp.dump();
#endif
  if(cp.component(curstate.p[0]) == cp.component(curstate.p[1])) {
    // start-midgame: try to cut off our opponent
    return next_move_alphabeta();
  } else {
    // endgame: use up space as efficiently as we can, and hope we have more
    // left than they do.
    return next_move_spacefill();
  }
}

int main() {
  setlinebuf(stdout);
  while (map_update()) {
    printf("%d\n", next_move());
  }
//#if VERBOSE >= 1
//  fprintf(stderr, "%d evaluations\n", evaluations);
//#endif
  return 0;
}

// vim: sw=2:ts=8:et:foldmethod=marker
