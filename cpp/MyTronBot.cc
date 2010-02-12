#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int map_width, map_height;
char *map = NULL;
#define M(p) map[(p).x+(p).y*map_width]

struct position {
  static const int dx[5], dy[5];

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

const int position::dx[5]={0,0,1,0,-1};
const int position::dy[5]={0,-1,0,1,0};


struct gamestate {
  position p[2]; // position in current state
  int m[2]; // last move made

  // affects map using m[0], m[1]
  gamestate move() {
    gamestate s = *this;
    M(p[0]) = 1;
    M(p[1]) = 1;
    s.p[0] = p[0].next(m[0]);
    s.p[1] = p[1].next(m[1]);
    return s;
  }
  // undoes effect on map
  void unmove() {
    M(p[0]) = 0;
    M(p[1]) = 0;
  }
};

gamestate curstate;

gamestate try_move(gamestate r, int player, int move) {
  // only when player 2 moves does the state actually progress...  which means you need to enter both moves.
  gamestate s = r;
  s.m[player] = move;
  if(player == 1) {
    // etc etc
  }
  return s;
}

// {{{ imported map update garbage from original code
void map_update()
{
  int x, y, c;
  int num_items = fscanf(stdin, "%d %d\n", &map_width, &map_height);
  if (feof(stdin) || num_items < 2) {
    exit(0); // End of stream means end of game. Just exit.
  }
  if(!map) map = new char[map_width * map_height];
  memset(map, 0, map_width*map_height);
  x = 0;
  y = 0;
  while (y < map_height && (c = fgetc(stdin)) != EOF) {
    switch (c) {
    case '\r':
      break;
    case '\n':
      if (x != map_width) {
	fprintf(stderr, "x != width in Board_ReadFromStream\n");
	return;
      }
      ++y;
      x = 0;
      break;
    case '#':
      if (x >= map_width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return;
      }
      M(position(x,y)) = 1;
      ++x;
      break;
    case ' ':
      if (x >= map_width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return;
      }
      M(position(x,y)) = 0;
      ++x;
      break;
    case '1':
    case '2':
      if (x >= map_width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return;
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
      return;
    }
  }
}
// }}}

int degree(position x) {
  return 4 - M(x.next(1)) - M(x.next(2)) - M(x.next(3)) - M(x.next(4));
}

int next_move() {
  int bestm = 0;
  int best = 50;
  M(curstate.p[0]) = 1;
  M(curstate.p[1]) = 1;
  for(int m=1;m<=4;m++) {
    position x = curstate.p[0].next(m);
    if(M(x)) continue;
    int d = degree(x);
    if(bestm == 0 || abs(d-2) < best) {
      bestm = m;
      best = d;
    }
  }
  //fprintf(stderr, "moving %d; degree=%d\n", bestm, best);
  return bestm;
}

int main() {
  setlinebuf(stdout);
  while (true) {
    map_update();
    printf("%d\n", next_move());
  }
  return 0;
}
